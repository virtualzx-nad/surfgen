! load coordinate sets

SUBROUTINE readCoordSets()
  use hddata,  only:ncoord,order,getFLUnit
  use progdata,only:nrij,nout,printlvl,natoms,atomCount,atomList,&
       nCoordSets,CoordSet,coordmap,nCoordCond,CoordCond,condRHS
  implicit none
  character(72)            :: comment
  character(4)             :: str
  integer                  :: ios,i,j,k,l,m,n,&
                              CSETFL,       &  !Unit ID for coordinate definition file
                              nAddCond         !number of additional conditions
! temporary variables used to generate definitions
  integer                  :: tmpCoord(4,natoms**4),& ! temporary coordinate definition holder
                              rawCoord(4),ordOOP(4),& ! atoms referenced, before and after ordering
                              prty                    ! parity of OOP angle, just a dummy variable
  integer,dimension(:),allocatable   ::  lhs          ! Left hand side of coordinate order restrictions
  CSETFL=getFLUnit()
! total coordinates count
  ncoord=0
! total coordinate condition count
  nCoordCond=0

! Load coordinate definition file
  if(printlvl>0)print *,"   Reading coordinate set definitions."
  open(unit=CSETFL,file='coord.in',access='sequential',form='formatted',&
    STATUS='OLD',ACTION='READ',POSITION='REWIND',IOSTAT=ios)
  read(CSETFL,1000,IOSTAT=ios) comment
! Get number of coordinate sets and additional conditions
  read(CSETFL,*,IOSTAT=ios) nCoordSets,nAddCond
  if(ios/=0.or.nCoordSets<1)stop"Error reading coord set definitions."

  if(allocated(CoordSet))deallocate(CoordSet)
  allocate(CoordSet(nCoordSets))

  if(printlvl>0)Print *,"      generating",nCoordSets," sets"  
  allocate(lhs(nCoordSets))

  do i=1,nCoordSets
! read in definition for one set of coordinates
    read(CSETFL,1000,IOSTAT=ios) comment
    read(CSETFL,*,IOSTAT=ios)CoordSet(i)%Type,CoordSet(i)%Scaling,CoordSet(i)%Order
    read(CSETFL,*,IOSTAT=ios)CoordSet(i)%atom
    if(ios/=0)stop "Error reading coord set definitions."
! check order settings
    if(CoordSet(i)%Order<0)stop "Error:  Wrong maximum order value"
    if(CoordSet(i)%Order==0.or.CoordSet(i)%Order>=order)then
       CoordSet(i)%Order=order
    else
       nCoordCond=nCoordCond+1
    end if!(CoordSet(i)%Order==0.or.CoordSet(i)%Order>order)
    if(printlvl>0)Print '("     set ",I4," type=",I4 ," scaling=",I4," max order=",I4)',&
                                i,CoordSet(i)%Type,CoordSet(i)%Scaling,CoordSet(i)%Order

    select case(CoordSet(i)%Type)

!----Plain or Scaled Internuclear distance coordinate

      case(0)  !rij or scaled rij.   %atom = A1 , A2 ,  X , X

! field 1 and 2 of atom record the atom group the end point atoms belong to
! field 3 and 4 are not used
          if(CoordSet(i)%atom(1).eq.CoordSet(i)%atom(2))theN
! bond between same type of atoms
! determine number of coordinates contained in the set
            CoordSet(i)%ncoord  =  atomCount(CoordSet(i)%atom(1))  *     &
                        (atomCount(CoordSet(i)%atom(1))-1)/2
! allocate memory for field Coord, definition of all the coordinates in the set
            if(allocated(CoordSet(i)%coord))deallocate(CoordSet(i)%coord)
            allocate(CoordSet(i)%coord(2,CoordSet(i)%ncoord))
! make coord list
            l = 1
            do j=1,atomCount(CoordSet(i)%atom(1))
              do k=j+1,atomCount(CoordSet(i)%atom(1))
                CoordSet(i)%coord(1,l) = atomList(CoordSet(i)%atom(1),j)
                CoordSet(i)%coord(2,l) = atomList(CoordSet(i)%atom(1),k)
                l=l+1
              end do !k
            end do !j
          else   !(CoordSet(i)%atom(1).eq.CoordSet(i)%atom(2))
! bond between different types of atoms
            CoordSet(i)%ncoord  =  atomCount(CoordSet(i)%atom(1))  *     &
                                        atomCount(CoordSet(i)%atom(2))
            allocate(CoordSet(i)%coord(2,CoordSet(i)%ncoord))
! make coord list
            l = 1
            do j=1,atomCount(CoordSet(i)%atom(1))
              do k=1,atomCount(CoordSet(i)%atom(2))
                CoordSet(i)%coord(1,l) = atomList(CoordSet(i)%atom(1),j)
                CoordSet(i)%coord(2,l) = atomList(CoordSet(i)%atom(2),k)
                l=l+1
              end do !k
            end do !j
          end if !(CoordSet(i)%atom(1).eq.CoordSet(i)%atom(2))
! Read in scaling parameters if is scaled Rij
        allocate(CoordSet(i)%Coef(2))
        if(CoordSet(i)%Scaling.ne.0)then
          read(CSETFL,'(2F10.6)',IOSTAT=ios) CoordSet(i)%Coef(:)
          if(ios/=0)stop "Error reading coord set definitions."
        end if !(CoordSet(i)%Scaling.ne.0)


!-- Out of plane angle coordinate


      case(-1) ! OOP.   %atom = A1, A2, A3, A4

        CoordSet(i)%ncoord=0
! Generate all possible 4 combinations.  Only the ones that are cannonically
! ordered will be inserted into the list.
        do j=1,atomCount(CoordSet(i)%atom(1))
          rawCoord(1) =  atomList(CoordSet(i)%atom(1),j) 
          do k=1,atomCount(CoordSet(i)%atom(2))
            rawCoord(2) =  atomList(CoordSet(i)%atom(2),k) 
            if(rawCoord(2)==rawCoord(1))cycle
            do m=1,atomCount(CoordSet(i)%atom(3))
              rawCoord(3) =  atomList(CoordSet(i)%atom(3),m) 
              if(count(rawCoord(3)==rawCoord(1:2))>0)cycle
              do n=1,atomCount(CoordSet(i)%atom(4))
                rawCoord(4) =  atomList(CoordSet(i)%atom(4),n) 
                if(count(rawCoord(4).eq.rawCoord(1:3))>0)cycle  !all indices have to be different
                call reorderOOP(rawCoord,ordOOP,prty)
                if(count(rawCoord.eq.ordOOP).eq.4)then ! the atoms are cannonically ordered
                  CoordSet(i)%ncoord = CoordSet(i)%ncoord+1
                  tmpCoord( :, CoordSet(i)%ncoord) = rawCoord
                end if ! count(...).eq.4
              end do !n
            end do !m
          end do !k
        end do ! j
! allocate field %coord and transfer definition from temporary array to global structure
        if(allocated(CoordSet(i)%coord))deallocate(CoordSet(i)%coord)
        allocate(CoordSet(i)%coord(4,CoordSet(i)%ncoord))
        CoordSet(i)%coord = tmpCoord(:,1:CoordSet(i)%ncoord)

! load scaling factor
        allocate(CoordSet(i)%Coef(2))
        read(CSETFL,'(2F10.6)',IOSTAT=ios) CoordSet(i)%Coef
        if(ios/=0)stop "Error reading coord set definitions."

      case(-2) ! OOP.   %atom = A1, A2, A3, A4

        CoordSet(i)%ncoord=0
! Generate all possible 3 atom combinations.  Only the ones that are cannonically
! ordered will be inserted into the list.
        do j=1,atomCount(CoordSet(i)%atom(1))
          rawCoord(1) =  atomList(CoordSet(i)%atom(1),j) 
          do k=1,atomCount(CoordSet(i)%atom(2))
            rawCoord(2) =  atomList(CoordSet(i)%atom(2),k) 
            if(rawCoord(2)==rawCoord(1))cycle
            do m=1,atomCount(CoordSet(i)%atom(3))
              rawCoord(3) =  atomList(CoordSet(i)%atom(3),m) 
              if(count(rawCoord(3)==rawCoord(1:2))>0)cycle
              do n=1,atomCount(CoordSet(i)%atom(4))
                rawCoord(4) =  atomList(CoordSet(i)%atom(4),n) 
                if(count(rawCoord(4).eq.rawCoord(1:3))>0)cycle  !all indices have to be different
                call reorderOOP2(rawCoord,ordOOP,prty)
                if(count(rawCoord.eq.ordOOP).eq.4)then ! the atoms are cannonically ordered
                  CoordSet(i)%ncoord = CoordSet(i)%ncoord+1
                  tmpCoord( :, CoordSet(i)%ncoord ) = rawCoord
                end if ! count(...).eq.4
              end do !n
            end do !m
          end do !k
        end do ! j
! allocate field %coord and transfer definition from temporary array to global structure
        if(allocated(CoordSet(i)%coord))deallocate(CoordSet(i)%coord)
        allocate(CoordSet(i)%coord(4,CoordSet(i)%ncoord))
        CoordSet(i)%coord = tmpCoord(:,1:CoordSet(i)%ncoord)

! load scaling factor
        allocate(CoordSet(i)%Coef(2))
        read(CSETFL,'(2F10.6)',IOSTAT=ios) CoordSet(i)%Coef
        if(ios/=0)stop "Error reading coord set definitions."


! no coefficient to read for OOP

!---Bond angle coordinates and their periodic scalings

      case(1)  ! angle A1-A2-A3.   %atom = A1, A2, A3, X. 

        CoordSet(i)%ncoord=0
! Generate all possible 4 combinations.  Only the ones that are cannonically
! ordered will be inserted into the list. (ie, A3>A1 when same type, and that
! A1,A2,A3 are all different )
        do j=1,atomCount(CoordSet(i)%atom(1))
          rawCoord(1) =  atomList(CoordSet(i)%atom(1),j)
          do k=1,atomCount(CoordSet(i)%atom(2))
            rawCoord(2) =  atomList(CoordSet(i)%atom(2),k)
            if(rawCoord(2).eq.rawCoord(1))cycle
            do l=1,atomCount(CoordSet(i)%atom(3))
              rawCoord(3) =  atomList(CoordSet(i)%atom(3),l)
              if(CoordSet(i)%atom(1).eq.CoordSet(i)%atom(3)  .and. &
                          rawCoord(3)<=rawCoord(1) .or. rawCoord(3).eq.rawCoord(2))  cycle
              CoordSet(i)%ncoord = CoordSet(i)%ncoord+1
              tmpCoord( 1:3, CoordSet(i)%ncoord ) = rawCoord(1:3)
            end do !l
          end do !k
        end do ! j

! allocate field %coord and transfer definition from temporary array to global structure
        if(allocated(CoordSet(i)%coord))deallocate(CoordSet(i)%coord)
        allocate(CoordSet(i)%coord(3,CoordSet(i)%ncoord))
        CoordSet(i)%coord = tmpCoord(1:3,1:CoordSet(i)%ncoord)

! no scaling parameters for bond angles

!---Torsion angle coordinates and their periodic scalings

      case (2) ! torsion A1-A2-A3-A4.   %atom = A1, A2, A3, A4

        CoordSet(i)%ncoord=0
! cannonical order = A3>A2 and all atom indices are different
        do j=1,atomCount(CoordSet(i)%atom(1))
          rawCoord(1) =  atomList(CoordSet(i)%atom(1),j)
          do k=1,atomCount(CoordSet(i)%atom(2))
            rawCoord(2) =  atomList(CoordSet(i)%atom(2),k)
            if(rawCoord(2).eq.rawCoord(1))cycle
            do m=1,atomCount(CoordSet(i)%atom(3))
              rawCoord(3) =  atomList(CoordSet(i)%atom(3),m)
              if(rawCoord(3).eq.rawCoord(1) .or. rawCoord(3)<=rawCoord(2))cycle
              do n=1,atomCount(CoordSet(i)%atom(4))
                rawCoord(4) =  atomList(CoordSet(i)%atom(4),n)
                if(count(rawCoord(4).eq.rawCoord(1:3))>0)cycle  !all indices have to be different
                CoordSet(i)%ncoord = CoordSet(i)%ncoord+1
                tmpCoord( :, CoordSet(i)%ncoord ) = rawCoord
              end do !n
            end do !m
          end do !k
        end do ! j

! allocate field %coord and transfer definition from temporary array to global structure
        if(allocated(CoordSet(i)%coord))deallocate(CoordSet(i)%coord)
        allocate(CoordSet(i)%coord(4,CoordSet(i)%ncoord))
        CoordSet(i)%coord = tmpCoord(:,1:CoordSet(i)%ncoord)

! no coefficient to read for torsion

!-- Other types

      case default
        print *,"TYPE=",CoordSet(i)%Type
        stop "Error : COORDINATE TYPE NOT SUPPORTED. "
    end select!case(CoordSet(i)%Type)

! %icoord  maps coordinate index inside set to full coordinate list index

    allocate(CoordSet(i)%icoord(CoordSet(i)%ncoord))
    do j=1,CoordSet(i)%ncoord
      CoordSet(i)%icoord(j)=j+ncoord
    end do!j=1,CoordSet(i)%ncoord
    ncoord=ncoord+CoordSet(i)%ncoord

  end do!i=1,nCoordSets

! generate coordinate map from full coordinate list to coordinate set list
  if(allocated(coordmap))deallocate(coordmap)
  allocate(coordmap(ncoord,2))
  k=0
  do i=1,nCoordSets
    do j=1,CoordSet(i)%ncoord
      k=k+1
      coordmap(k,1)=i   !(:,1) = index of set
      coordmap(k,2)=j   !(:,2) = index within a set
    end do!j=1,CoordSet(i)%ncoord
  end do!i=1,nCoordSets

!Generate coordinate inequality conditions for hddata
  allocate(CoordCond(nCoordCond+nAddCond,ncoord))
  allocate(condRHS(nCoordCond+nAddCond))
  nCoordCond=0
  CoordCond=0

! Generate individual set order limits
  do i=1,nCoordSets
    if(CoordSet(i)%Order>=order)cycle
    nCoordCond=nCoordCond+1
    CoordCond(nCoordCond,CoordSet(i)%icoord)=1
    condRHS(nCoordCond)=CoordSet(i)%Order
  end do!i=1,nCoordSets

! Additional set order restrictions
  write(str,'(I4)') nCoordSets+1
  if(printlvl>0)print *,"      Reading ",naddcond," additional order restrictions."
  do i=1,nAddCond
    nCoordCond=nCoordCond+1
    read(CSETFL,*,IOSTAT=ios) lhs,condRHS(nCoordCond)
    if(ios/=0)stop "Error reading addintional coordinate order restrictions."
    if(printlvl>1)print *,"      lhs: ",lhs," rhs:",condRHS(nCoordCond)
    do j=1,nCoordSets
      CoordCond(nCoordCond,CoordSet(j)%icoord)=lhs(j)
    end do
  end do
  close(unit=CSETFL)
  deallocate(lhs)
1000 format(72a)
END SUBROUTINE

!-------------------------------------------------------------------------------
! read coordinate definitions.   coordinates are defined by atom index instead
! of atom group index.  the permutated atoms are generated using the list of 
! feasible permutations.  this subroutine is only
! load coordinate sets
SUBROUTINE readCoords()
  use hddata,  only:ncoord,order,getFLUnit
  use progdata,only:printlvl,nCoordSets,CoordSet,coordmap,nCoordCond,CoordCond,condRHS
  use CNPI, only: nPmt, pmtList
  implicit none
  character(72) :: comment
  character(4)  :: str
  logical       :: found
  integer       :: ios,i,j,k,l,newcrd(4), &
                   CSETFL,         &  !Unit ID for coordinate definition file
                   nAddCond           !number of additional conditions
! temporary variables used to generate definitions
  integer       :: rawCoord(4),ordOOP(4),& ! atoms referenced, before and after ordering
                   prty                    ! parity of OOP angle, just a dummy variable
  integer,dimension(:),allocatable   ::  lhs   ! Left hand side of order restrictions
  integer,dimension(:,:),allocatable ::  tmpCoord   ! temporary coord group list

  CSETFL=getFLUnit()
! total coordinates count
  ncoord=0
! total coordinate condition count
  nCoordCond=0

! Load coordinate definition file
  if(printlvl>0)print *,"   Reading coordinate set definitions."
  open(unit=CSETFL,file='coord.in',access='sequential',form='formatted',&
            STATUS='OLD',ACTION='READ',POSITION='REWIND',IOSTAT=ios)
  read(CSETFL,1000,IOSTAT=ios) comment
! Get number of coordinate sets and additional conditions
  read(CSETFL,*,IOSTAT=ios) nCoordSets,nAddCond
  if(ios/=0.or.nCoordSets<1)stop"Error reading coord set definitions."

  if(allocated(CoordSet))deallocate(CoordSet)
  allocate(CoordSet(nCoordSets))

  if(printlvl>0)Print *,"      generating",nCoordSets," sets"
  allocate(lhs(nCoordSets))

  do i=1,nCoordSets
! read in definition for one set of coordinates
    read(CSETFL,1000,IOSTAT=ios) comment
    read(CSETFL,*,IOSTAT=ios)CoordSet(i)%Type,CoordSet(i)%Scaling,&
                                    CoordSet(i)%Order
! here, atom stores atom index instead of atom group index
    read(CSETFL,*,IOSTAT=ios)CoordSet(i)%atom
    if(ios/=0)stop "Error reading coord set definitions."
! check order settings
    if(CoordSet(i)%Order<0)stop "Error:  Wrong maximum order value"
    if(CoordSet(i)%Order==0.or.CoordSet(i)%Order>=order)then
        CoordSet(i)%Order=order
    else
        nCoordCond=nCoordCond+1
    end if!(CoordSet(i)%Order==0.or.CoordSet(i)%Order>order)
    if(printlvl>0)Print '(5x,"set ",I4," type=",I4 ," scaling=",I4," max order=",I4)',&
            i,CoordSet(i)%Type,CoordSet(i)%Scaling,CoordSet(i)%Order

    select case(CoordSet(i)%Type)

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!               Plain or Scaled Internuclear distance coordinate
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      case(0)  !rij or scaled rij.   %atom = A1 , A2 ,  X , X

      ! field 1 and 2 of atom record the atom index of the end point atoms
      ! field 3 and 4 are not used
        if(CoordSet(i)%atom(1).eq.CoordSet(i)%atom(2))then
        ! bond between same atom is not allowed
            print *,"Bond between the same atom is not allowed!!!"
            stop "Error: Bond between the same atom is not allowed!!!"
        end if  !(CoordSet(i)%atom(1).eq.CoordSet(i)%atom(2))
        ! bond between different atoms
        ! allocate temporary coord set list
        allocate(tmpCoord(2,nPmt))
        CoordSet(i)%ncoord    = 0
        do l=1,nPmt
            newcrd(1)=pmtList(l,CoordSet(i)%atom(1))
            newcrd(2)=pmtList(l,CoordSet(i)%atom(2))
            if(newcrd(1)>newcrd(2))then
                newcrd(3) = newcrd(1)
                newcrd(1) = newcrd(2)
                newcrd(2) = newcrd(3)
            end if
            ! look up if this coordinate is already defined in the list
            found = .false.
            do j=1,CoordSet(i)%ncoord
                if(newcrd(1)==tmpCoord(1,j).and. &
                   newcrd(2)==tmpCoord(2,j) ) then
                    found = .true.
                    exit
                end if
            end do!j
            if(.not.found)then
                CoordSet(i)%ncoord = CoordSet(i)%ncoord+1
                tmpCoord(:,CoordSet(i)%ncoord) = newcrd(1:2)
            end if!.not.found
        end do!l=1,nPmt

        ! allocate coord list
        allocate(CoordSet(i)%coord(2,CoordSet(i)%ncoord))
        CoordSet(i)%coord = tmpCoord(:,1:CoordSet(i)%ncoord)
        deallocate(tmpCoord)

        ! read scaling coefficients
        allocate(CoordSet(i)%Coef(2))
        if(CoordSet(i)%Scaling.ne.0)then
            read(CSETFL,*,IOSTAT=ios) CoordSet(i)%Coef(:)
            if(ios/=0)stop "Error reading coord set definitions."
        end if !(CoordSet(i)%Scaling.ne.0)

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!                   Out of plane angle coordinate
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      case(-1,-2) ! OOP.   %atom = A1, A2, A3, A4

        CoordSet(i)%ncoord=0
        allocate(tmpCoord(4,nPmt))
    ! Generate all possible permutations of the 4 atoms in question.
    ! Coordinates are inserted in cannonical order
        do l=1,nPmt
            rawCoord = pmtList(l,CoordSet(i)%atom)
            if(CoordSet(i)%Type==-1)then
                call reorderOOP(rawCoord,ordOOP,prty)
            else
                call reorderOOP2(rawCoord,ordOOP,prty)
            end if
            ! look up the reordered coordinate in the temporary coordinate list
            found = .false.
            do j=1,CoordSet(i)%ncoord
                if(count(tmpCoord(:,j).eq.ordOOP).eq.4)then
                    found =.true.
                    exit
                end if
            end do!j
            if(.not.found)then
                CoordSet(i)%ncoord = CoordSet(i)%ncoord + 1
                tmpCoord(:,CoordSet(i)%ncoord) = ordOOP
            end if
        end do!l

        ! allocate field %coord and transfer definition from temporary array to global structure
        if(allocated(CoordSet(i)%coord))deallocate(CoordSet(i)%coord)
        allocate(CoordSet(i)%coord(4,CoordSet(i)%ncoord))
        CoordSet(i)%coord = tmpCoord(:,1:CoordSet(i)%ncoord)
        deallocate(tmpCoord)

        ! load scaling factor
        allocate(CoordSet(i)%Coef(2))
        read(CSETFL,*,IOSTAT=ios) CoordSet(i)%Coef
        if(ios/=0)stop "Error reading coord set definitions."

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!           Bond angle coordinates and their periodic scalings
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      case(1)  ! angle A1-A2-A3.   %atom = A1, A2, A3, X.

! atom 1 and 2 are end points and atom 3 is the vertex.  atom 4 is not used

        if(CoordSet(i)%atom(1).eq.CoordSet(i)%atom(2).or.  & ! check validness
            CoordSet(i)%atom(1).eq.CoordSet(i)%atom(3).or.  & ! of definition
            CoordSet(i)%atom(2).eq.CoordSet(i)%atom(3) )then
            ! angle definition cannot involve identical atoms
            print *,"Angle between the same atom is not allowed!!!"
            stop "Error: Angle between the same atom is not allowed!!!"
        end if  !(CoordSet(i)%atom(1).eq.CoordSet(i)%atom(2))

        ! allocate temporary coord set list
        allocate(tmpCoord(3,nPmt))
        CoordSet(i)%ncoord    = 0
        do l=1,nPmt
            newcrd(1) = pmtList(l,CoordSet(i)%atom(1))
            newcrd(2) = pmtList(l,CoordSet(i)%atom(2))
            if(newcrd(1)>newcrd(2))then
                newcrd(3) =  newcrd(1)
                newcrd(1) =  newcrd(2)
                newcrd(2) =  newcrd(3)
            end if
            newcrd(3) = pmtList(l,CoordSet(i)%atom(3))
            ! look up if this coordinate is already defined in the list
            found = .false.
            do j=1,CoordSet(i)%ncoord
                if( newcrd(1)==tmpCoord(1,j).and. &
                    newcrd(2)==tmpCoord(2,j).and. &
                    newcrd(3)==tmpCoord(3,j) ) then
                    found = .true.
                    exit
                end if
            end do!j
            if(.not.found)then
                CoordSet(i)%ncoord = CoordSet(i)%ncoord+1
                tmpCoord(:,CoordSet(i)%ncoord) = newcrd(1:3)
            end if!.not.found
        end do!l=1,nPmt

        ! allocate coord list
        allocate(CoordSet(i)%coord(3,CoordSet(i)%ncoord))
        CoordSet(i)%coord = tmpCoord(:,1:CoordSet(i)%ncoord)
        deallocate(tmpCoord)

        ! load scaling factor
        allocate(CoordSet(i)%Coef(2))
        if(CoordSet(i)%Scaling==2)then
            read(CSETFL,*,IOSTAT=ios) CoordSet(i)%Coef
            if(ios/=0)stop "Error reading coord set definitions."
        else
            CoordSet(i)%Coef = 0d0
        end if

! no scaling parameters for bond angles

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!           Torsion angle coordinates and their periodic scalings
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      case (2) ! torsion A1-A2-A3-A4.   %atom = A1, A2, A3, A4

        CoordSet(i)%ncoord=0
        ! TORSION NOT YET IMPLEMENTD

! no coefficient to read for torsion

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!                   Other types.  Should not get here
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      case default
        print *,"TYPE=",CoordSet(i)%Type
        stop "Error : COORDINATE TYPE NOT SUPPORTED. "
    end select!case(CoordSet(i)%Type)

! %icoord  maps coordinate index inside set to full coordinate list index

    allocate(CoordSet(i)%icoord(CoordSet(i)%ncoord))
    do j=1,CoordSet(i)%ncoord
        CoordSet(i)%icoord(j)=j+ncoord
    end do!j=1,CoordSet(i)%ncoord
    ncoord=ncoord+CoordSet(i)%ncoord

  end do!i=1,nCoordSets

  ! generate coordinate map from full coordinate list to coordinate set list
  if(allocated(coordmap))deallocate(coordmap)
  allocate(coordmap(ncoord,2))
  k=0
  do i=1,nCoordSets
    do j=1,CoordSet(i)%ncoord
        k=k+1
        coordmap(k,1)=i   !(:,1) = index of set
        coordmap(k,2)=j   !(:,2) = index within a set
    end do!j=1,CoordSet(i)%ncoord
  end do!i=1,nCoordSets

  !Generate coordinate inequality conditions for hddata
  allocate(CoordCond(nCoordCond+nAddCond,ncoord))
  allocate(condRHS(nCoordCond+nAddCond))
  nCoordCond=0
  CoordCond=0

  ! Generate individual set order limits
  do i=1,nCoordSets
    if(CoordSet(i)%Order>=order)cycle
    nCoordCond=nCoordCond+1
    CoordCond(nCoordCond,CoordSet(i)%icoord)=1
    condRHS(nCoordCond)=CoordSet(i)%Order
  end do!i=1,nCoordSets

  ! Additional set order restrictions
  write(str,'(I4)') nCoordSets+1
  if(printlvl>0)print *,"      Reading ",naddcond," additional order restrictions."
  do i=1,nAddCond
    nCoordCond=nCoordCond+1
    read(CSETFL,*,IOSTAT=ios) lhs,condRHS(nCoordCond)
    if(ios/=0)stop "Error reading addintional coordinate order restrictions."
    if(printlvl>1)then
        write(*,"(2x,A)",advance='no') "lhs:"
        write(*,"(15I4)") lhs
        print "(2x,A,I5)","rhs:",condRHS(nCoordCond)
    end if
    do j=1,nCoordSets
        CoordCond(nCoordCond,CoordSet(j)%icoord)=lhs(j)
    end do
  end do
  close(unit=CSETFL)
  deallocate(lhs)
1000 format(72a)
END SUBROUTINE readCoords

!-------------------------------------------------------------------------------
! load irrep matrices from irrep.in
SUBROUTINE readIrreps()
  use progdata, only:IRREPFL,printlvl
  use CNPI, only:  nirrep, irrep, deallocIrreps
  IMPLICIT NONE
  character(72)            ::  comment
  character(3)             ::  tpStr
  integer                  ::  i,j,k,d,ordr,ios

  call deallocIrreps()
  open(unit=IRREPFL,file='irrep.in',access='sequential',form='formatted',&
      STATUS='OLD',ACTION='READ',POSITION='REWIND',IOSTAT=ios)
  if(ios/=0)stop "readIrreps:  cannot open file irrep.in"
  read(IRREPFL,1000) comment
  read(IRREPFL,*) nirrep
  if(printlvl>0)print *,"    Loading irreducible representations...."
  if(printlvl>1)print *,"    ",trim(comment)
  if(printlvl>1)print *,"      Number of irreps: ",nirrep
  allocate(irrep(nirrep))
  do i=1,nirrep
    read(IRREPFL,1000) comment
    if(printlvl>2)print '(7x,"Irrep #",I2," : ",a)',i,trim(adjustl(comment))
    read(IRREPFL,*)  d,ordr
    irrep(i)%Dim   = d
    irrep(i)%Order = ordr
    write(tpStr,'(I3)') d
    allocate(irrep(i)%RepMat(ordr,d,d))
    do j=1,ordr
      do k=1,d
        read(IRREPFL,*) irrep(i)%repmat(j,k,1:d)
        if(printlvl>2)print '(11X,'//trim(tpStr)//'F14.10)',irrep(i)%repmat(j,k,1:d)
      end do!k=1,d
      if(printlvl>2)print *,''
    end do!j=1,ordr
  end do!i=1,nirrep
  close(unit=IRREPFL)
1000 format(72a)
END SUBROUTINE

! read reference geometry/energy/internal coordiante/gradients from input file
! also generates global coordiantes and symmetrized polynomial basis
! b matrices are calcualted and reference geometry is symmetrized
SUBROUTINE initialize(jobtype)
  use hddata
  use progdata
  use CNPI
  IMPLICIT NONE
  INTEGER,INTENT(IN)                          :: jobtype
  integer :: i
  double precision  :: eguess(nGroups)

  LOGICAL removed

  if(printlvl>0)print *,"Entering Initialize()"
  if(printlvl>0)print *,"  Generating atom permutations"
  Call genAtomPerm()
  if(jobtype<0)then
    print *,"  All atom permutations:"
    do i=1,nPmt
        write(unit=*,fmt="(8x,A,I5,A)",Advance="NO") "Pmt #",i,":"
        PRINT "(20I3)",pmtList(i,:)
    end do!i=1,nPmt
  end if
  if(printlvl>0)print *,  "  Selecting feasible permutations"
  Call selectAtomPerm(removed)
  if(.not.removed)print *,"  No permutations removed."
  if(jobtype<0)then
    print *,"Finished printing permutations."
    stop
  end if

!  if(removed)then
!    if(printlvl>0)print *,"  Using atom indexes to define coordinates"
    call readCoords()
!  else!removed
!    if(printlvl>0)print *,"  Using atom group index to define coordinates"
!    call readCoordSets()
!  end if!removed
  if(printlvl>0)print *,"  Generating permutations for scaled coordinates"
  call genCoordPerm()
  CALL initHd()
  if(printlvl>0)print *,"  Generating term list"
  Call genTermList(nCoordCond,CoordCond,condRHS)
  if(printlvl>0)print *,"  Partitioning term list"
  Call PartitionTerms()

  if(printlvl>0)print *,"  Generating Maptab"
  call genMaptab()

  !allocate space for Hd  
  do i=1,nGroups
    eguess(i)=(i-1)*10000/AU2CM1
  end do
  call allocateHd(eguess)
  if(inputfl/='')then 
    if(printlvl>0)print *,"  Reading Hd Coefficients from ",trim(inputfl)
    call readHd(inputfl)
  end if!(inputfl/='')
  call printTitle(jobtype)
  
  if(printlvl>0)print *,"Exiting Initialize()"
end SUBROUTINE initialize



!
!
!
!
SUBROUTINE printTitle(jobtype)
  use progdata,only: OUTFILE
  IMPLICIT NONE
  INTEGER,INTENT(IN)          :: jobtype
  CHARACTER(9),dimension(-1:2) :: types = (/'GENSYM   ','POTLIB   ','MAKESURF ',&
                                            'EXTREMA  ' /)

  open(unit=OUTFILE,file='surfgen.out',access='sequential',form='formatted')
  write(OUTFILE,1000)'-----------------------------------------------------------'
  write(OUTFILE,1000)'                                                           '
  write(OUTFILE,1000)'                     surfgen.global                        '
  write(OUTFILE,1000)'                                                           '
  write(OUTFILE,1000)'   * Creates quasi-diabatic Hamiltonian by reproducing     '
  write(OUTFILE,1000)'   ab initio energy gradients and derivative couplings     '
  write(OUTFILE,1000)'   with selected data fitted exactly and the rest in a     '
  write(OUTFILE,1000)'   least-squares sense                                     '
  write(OUTFILE,1000)'                                                           '
  write(OUTFILE,1000)'   * Blocks of Hd are expanded as polynomials of scaled    '
  write(OUTFILE,1000)'   internuclear distance and out-of-plane angle coords     '
  write(OUTFILE,1000)'   that are well-defined globally.                         ' 
  write(OUTFILE,1000)'                                                           '
  write(OUTFILE,1000)'   * Hd is constructed, using projection operators, so     '
  write(OUTFILE,1000)'   that the diabats carry certain irreps of CNPI group     '
  write(OUTFILE,1000)'                                                           '
  write(OUTFILE,1000)'   Xiaolei Zhu,  2010                                      '
  write(OUTFILE,1000)'   Department of Chemistry, Johns Hopkins University       '
  write(OUTFILE,1000)'   based on SURFGEN.X, Michael Schuurman,  2008            '
  write(OUTFILE,1000)'-----------------------------------------------------------'
  write(OUTFILE,1000)''
  write(OUTFILE,1000)'  jobtype = '//types(jobtype)
  write(OUTFILE,1000)''
  write(OUTFILE,1000)''

1000 format(72a)
end SUBROUTINE printTitle

! output surface to file
!
!
!
!

!
!
!
SUBROUTINE cleanup()
  use hddata, only: cleanHdData
  use progdata
  use CNPI
  IMPLICIT NONE
  INTEGER                 :: i

  if(allocated(pmtList))deallocate(pmtList)
  if(allocated(subPerm))deallocate(subPerm)
  if(allocated(nSubPerm))deallocate(nSubPerm)
  if(allocated(coordPerm))deallocate(coordPerm)
  if(allocated(sgnCPerm))deallocate(sgnCPerm)

  CALL deallocIrreps()
  call deallocPCycle()
  CALL cleanHdData
  
  if(allocated(CoordSet))then
    do i=1,nCoordSets
      if(allocated(CoordSet(i)%iCoord))deallocate(CoordSet(i)%iCoord)
      if(allocated(CoordSet(i)%iCoord))deallocate(CoordSet(i)%Coef)
    end do
    deallocate(CoordSet)
  end if!(allocated(CoordSet))

  write(OUTFILE,1000)
  close(OUTFILE)
  return
1000 format(/,/,'  ---------------- EXECUTION COMPLETE -----------------',/)
end SUBROUTINE cleanup

! read COLUMBUS geom file and obtain geometry and atom info
SUBROUTINE readColGeom(gfile,ngeoms,na,atoms,anums,cgeom,masses)
  use hddata, only:  getFLUnit
  IMPLICIT NONE
  CHARACTER(72),INTENT(IN)                              :: gfile
  INTEGER,INTENT(IN)                                    :: na,ngeoms
  CHARACTER(3),dimension(na),INTENT(INOUT)              :: atoms
  DOUBLE PRECISION,dimension(na),INTENT(INOUT)          :: anums,masses
  DOUBLE PRECISION,dimension(3*na,ngeoms),INTENT(INOUT) :: cgeom
  INTEGER                                               :: i,j,k,GUNIT,ios

  GUNIT=getFLUnit()
  open(unit=GUNIT,file=trim(adjustl(gfile)),access='sequential',form='formatted',&
      status='old',action='read',position='rewind',iostat=ios)
  if(ios/=0)then
    print *,"gfile = [", trim(adjustl(gfile)),"]"
    stop"readColGeom: cannot open file for read"
  end if
  do i = 1,ngeoms
   do j = 1,na
    read(GUNIT,*)atoms(j),anums(j),(cgeom(3*(j-1)+k,i),k=1,3),masses(j)
   enddo
  enddo
  close(GUNIT)

  return
END SUBROUTINE readColGeom

! output COLUMBUS geom file
SUBROUTINE writeColGeom(gfile,na,atoms,anums,cgeom,masses)
  IMPLICIT NONE
  CHARACTER(72),INTENT(IN)                       :: gfile
  INTEGER,INTENT(IN)                             :: na
  CHARACTER(3),dimension(na),INTENT(INOUT)       :: atoms
  DOUBLE PRECISION,dimension(na),INTENT(INOUT)   :: anums,masses
  DOUBLE PRECISION,dimension(3*na),INTENT(INOUT) :: cgeom
  INTEGER                                        :: i,j,GUNIT

  GUNIT=11

  open(unit=GUNIT,file=trim(adjustl(gfile)),access='sequential',form='formatted')
  rewind(GUNIT)
  do i = 1,na
   write(GUNIT,1000)adjustl(atoms(i)),anums(i),(cgeom((i-1)*3+j),j=1,3),masses(i)
  enddo

  close(GUNIT)

  return
1000 format(1x,a3,2x,f4.0,1x,3(f13.8,1x),2x,f12.8)
END SUBROUTINE writeColGeom

!
!
!
!
SUBROUTINE readEner(efile,ngeoms,ne,eners)
  use hddata, only: getFLUnit
  IMPLICIT NONE
  CHARACTER(72),INTENT(IN)                              :: efile
  INTEGER,INTENT(IN)                                    :: ngeoms,ne
  DOUBLE PRECISION,dimension(ne,ngeoms),INTENT(INOUT)   :: eners
  INTEGER                                               :: i,j,EUNIT,ios

  EUNIT=getFLUnit()
  open(unit=EUNIT,file=trim(adjustl(efile)),access='sequential',form='formatted',&
    position='rewind',action='read',status='old',iostat=ios)
  if(ios/=0)stop'readEner: cannot open file for read'
  do i = 1,ngeoms
   read(EUNIT,*)(eners(j,i),j=1,ne)
  enddo
  close(EUNIT)
END SUBROUTINE readEner

!
!
!
SUBROUTINE readGrads(gfile,ngrads,na,cgrads)
  use hddata, only: getFLUnit
  IMPLICIT NONE
  CHARACTER(72),INTENT(IN)                              :: gfile
  INTEGER,INTENT(IN)                                    :: ngrads,na
  DOUBLE PRECISION,dimension(3*na,ngrads),INTENT(INOUT) :: cgrads
  INTEGER                                               :: i,j,k,GUNIT,ios
  GUNIT=getFLUnit()
  open(unit=GUNIT,file=trim(adjustl(gfile)),access='sequential',form='formatted',&
   action='read',position='rewind',status='old',iostat=ios)
  if(ios/=0)stop'readGrads:  cannot open file for read'
  do i = 1,ngrads
!   read(GUNIT,*)scr
   do j = 1,na
    read(GUNIT,*)(cgrads(3*(j-1)+k,i),k=1,3)
   enddo
  enddo
  close(GUNIT)
END SUBROUTINE readGrads

!
!
!
!
SUBROUTINE readHessian(gfile,nrc,hess)
  IMPLICIT NONE
  CHARACTER(72),INTENT(IN)                                :: gfile
  INTEGER,INTENT(IN)                                      :: nrc
  DOUBLE PRECISION,dimension(nrc*(nrc+1)/2),INTENT(INOUT) :: hess
  INTEGER                                                 :: i,GUNIT,ileft,nread,ndone

  GUNIT=11
  ndone = 0
  ileft = nrc*(nrc+1)/2
  open(unit=GUNIT,file=trim(adjustl(gfile)),access='sequential',form='formatted')
  do
   if(ileft.eq.0)EXIT
   nread = min(ileft,10)
   read(GUNIT,1000)(hess(ndone+i),i=1,nread)
   ndone = ndone + nread
   ileft = ileft - nread
  enddo
  close(GUNIT)

  return
1000 format(10(F10.6))
end SUBROUTINE readHessian
!
!
!
!
FUNCTION filename(s1,s2,suffix,usefij)
  IMPLICIT NONE
  INTEGER,INTENT(IN)          :: s1,s2
  LOGICAL ,INTENT(IN)         :: usefij
  CHARACTER(10),INTENT(IN)    :: suffix
  CHARACTER(72)               :: filename
  CHARACTER(1)                :: st1,st2

  write(st1,'(i1)')s1
  write(st2,'(i1)')s2

  if(s1.eq.s2)then
   filename = 'cartgrd.drt1.state'//st1//trim(adjustl(suffix))
  else
   if(usefij)then
     filEname = 'cartgrd_total.drt1.state'//st1//'.drt1.state'//st2//trim(adjustl(suffix))
   else
     filename = 'cartgrd.nad.drt1.state'//st1//'.drt1.state'//st2//trim(adjustl(suffix))
   end if
  endif

  filename = trim(adjustl(filename))
END FUNCTION filename


!printMatrix prints a matrix of double precision numbers to file
!ofile   :    output file
!rlabs   :    array(nr) of row label strings
!clabs   :    array(nc) of col label strings
!pcols   :    max number of columns displayed in a line
!nr      :    number of rows
!nc      :    number of columns
!mat     :    matrix to print
!fld,dcml:    printing format (F($fld).($dcml))
SUBROUTINE printMatrix(ofile,rlabs,clabs,pcols,nr,nc,mat,fld,dcml)
  IMPLICIT NONE
  INTEGER,INTENT(IN)                     :: ofile,pcols,nr,nc,fld,dcml
  CHARACTER(16),dimension(nr),INTENT(IN) :: rlabs
  CHARACTER(16),dimension(nc),INTENT(IN) :: clabs
  DOUBLE PRECISION,dimension(nr,nc),INTENT(IN) :: mat
  INTEGER                                :: i,j,k,ilo,ihi,nbatch
  INTEGER                                :: rlen,clen,flen,lspace,rspace
  CHARACTER(4)                           :: rlstr,clstr,lstr,rstr
  CHARACTER(72)                          :: FMT1,FMT2,FMT3,FMT4

  flen = fld
  rlen = 0
  clen = 0
  do i = 1,nr
   j = len_trim(rlabs(i))
   if(j.gt.rlen)rlen = j
  enddo
  do i = 1,nc
   j = len_trim(clabs(i))
   if(j.gt.clen)clen = j
  enddo

  write(rlstr,'(i4)')rlen
  write(clstr,'(i4)')clen
  FMT1 = '(/,'//trim(adjustl(rlstr))//'x)'
  FMT2 = '(a'//trim(adjustl(rlstr))//')'
  rspace = int((flen-clen)/2.)
  lspace = flen - clen - rspace
  if(rspace.lt.0)rspace=0
  if(lspace.lt.1)lspace=1
  write(lstr,'(i4)')lspace
  write(rstr,'(i4)')rspace
  FMT3 = '('//trim(adjustl(lstr))//'x,a'//trim(adjustl(clstr))//','//trim(adjustl(rstr))//'x)'
  write(lstr,'(i4)')flen-1
  write(rstr,'(i4)')dcml
  FMT4 = '(x,F'//trim(adjustl(lstr))//'.'//trim(adjustl(rstr))//')'

  nbatch = Ceiling(1.*nc/pcols)
  do i = 1,nbatch
   ilo = pcols*(i-1)+1
   ihi = min(pcols*i,nc)
   write(ofile,trim(FMT1),advance='no')
   do j = ilo,ihi
    write(ofile,trim(FMT3),advance='no')clabs(j)
   enddo
   write(ofile,1001,advance='no')
   do j = 1,nr
    write(ofile,trim(FMT2),advance='no')rlabs(j)
    do k=ilo,ihi
     write(ofile,trim(FMT4),advance='no')mat(j,k)
    enddo
    write(ofile,1001,advance='no')
   enddo
  enddo

  return
1001 format(/)
end SUBROUTINE printMatrix


