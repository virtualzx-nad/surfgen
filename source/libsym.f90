!********************************************************************
! Module contains information regarding the CNPI symmetry
! of the molecule and polynomial
!********************************************************************
MODULE CNPI

 use hddata, only: order,nblks,pTermDef,TTermDef,ncoord
 use progdata, only: printlvl,T2DList
 use combinatorial

 !----CONSTANTS------------------------------------------
 DOUBLE PRECISION,PARAMETER    :: IrrepMatCutOff=1D-7
 DOUBLE PRECISION,PARAMETER    :: MatProjCutOff=1D-7
 INTEGER,PARAMETER             :: MaxParents = 100

 !----DERIVED TYPES------------------------------------------
 ! Derived type for a permutation cycle
 ! A permutation cycle contains all the terms that are symmetry related to a
 ! specific term, which is homomorphic to a subgroup of the molecular symmetry
 ! group. The structure also contains the group structure of this subgroup and 
 ! its relation with the molecular symmetry group.
 !
 ! Symbols :
 ! G        Molecular symmetry group
 ! R        One symmetry operation of the molecular symmetry group
 ! P_R      The transformation operator corresponding to the symmetry operation
 !          R.  Also known as the Wigner operator of R.  
 ! fA       The root term from which the whole cycle is generated.
 ! W        The set of symmetry operations under which fA is invariant.  In
 !          other words this is the set of operations which take fA back to
 !          itself.  W is a subgroup if G.
 ! U        The set of symmetry operation which takes fA to -fA.  U is a left
 !          coset of W in G.
 ! V        V contains all the unique terms that are symmetry related to fA.
 !          Each element of V correspond to a left coset of (W+U) in G.  V is
 !          therefore also homomorphic to a subgroup of G.  In fact, G is a
 !          direct product group of (U+W) and V.
 ! U,V and W are used to compute the projection of the term to a specific
 ! representation.
 type TPermCycle
   integer                                :: nterms=1
   !term contains all P_R(fA) of l-cosets of U+W
   type(pTermDef),dimension(:),allocatable:: term
   integer,dimension(:),allocatable       :: sgnTerm !sign of each term
   integer,dimension(:),allocatable       :: W       !W={R|P_R(fA)=fA}
   integer                                :: nW    =0!|W|
   integer,dimension(:),allocatable       :: U       !U={R|P_R(fA)=-fA}.|U|=|W|
   integer                                :: nU    =0!|U|.  |U|=0 or |U|=|W|
   integer,dimension(:),allocatable       :: V       !one elem from each l-cosets
   integer                                :: parity=1!+1:even OOPs,-1:odd OOPs
   type(TPermCycle),pointer               :: pNext=>null()!pointer to the next node
 end type

 ! derived type for list of cycles
 type TPCycleList
   integer                                ::  nCycle=0
   type(TPermCycle),pointer               ::  handle
   type(TPermCycle),pointer               ::  last
 end type

 ! derived type for irredicible representations of state groups.
 ! Dim is dimensionality of irrep and Order is order of the group
 ! RepMat contains all representation matrices
 type TIrrep
   INTEGER                                      ::  Dim
   DOUBLE PRECISION,DIMENSION(:,:,:),allocatable::  RepMat
 end type

 ! irrep matrix type for symmetry blocks of Hd
 ! NParents: Number of parent state pairs that can generate this irrep mat
 ! parents:  List of parent state pairs that can generate this irrep mat
 type TBlkRepMat
   integer                                      :: NParents
   integer,dimension(MaxParents,2)              :: parents      
   integer                                      :: ldim,rdim     ! dim of state irrep
   double precision,dimension(:,:,:),allocatable:: RepMat       !irrep matrix
 end type 

!----GLOBAL VARIABLES-----------------------------------
 ! * SYMMETRY PROPERTIES *
 ! array that contains all the relevant irreducibles
 type(TIrrep),DIMENSION(:),allocatable  ::  irrep
 integer                                ::  nirrep=0

 INTEGER,dimension(:,:),allocatable     ::  GrpSym
 INTEGER,dimension(:,:),allocatable     ::  GrpPrty
 INTEGER                                ::  nSymLineUps

! * PERMUTATION PROPERTIES
! nPmt      number of feasible permutations
! nPmtAll   total number of permutations, including infeasible
! pmtList   list of all feasible permutations
! subPerm   all permutations of each atom group
! nSubPerm  number of permutations for each atom group
 INTEGER                                :: nPmt
 INTEGER                                :: nPmtAll
 INTEGER, DIMENSION(:,:), allocatable   :: pmtList !all permutation of atoms
 INTEGER, DIMENSION(:,:,:),allocatable  :: subPerm
 INTEGER, DIMENSION(:),allocatable      :: nSubPerm

! coordPerm  The the effect on the coordinates by permutations.  Permutation 
!            i would permute coordinate j into coordPerm(i,j), up to a change
!            of sign.  Wheter the sign changes or not is recorded in sgnCPerm.
!            Segment coordPerm(0,:) saves the parity of each coordinate.  In
!            other word, the character of inversion on coordinate k is
!            inversion(f_k)=coordPerm(0,k)*f_k
! sgnCPerm   Sign change associated with each permutation.
 INTEGER, DIMENSION(:,:), allocatable   :: coordPerm
 INTEGER, DIMENSION(:,:), allocatable   :: sgnCPerm

! term definitions
 type(TPCycleList),dimension(:),allocatable   ::  permCycle

! identical block symmetry list
! The program will try to detect blocks that carry the same group structure.
! NBlockSym is the number of blocks that have unique symmetry.
! Two lists will be made, one consists of the list of block symmetries each point
! to a block that have unique symmetry properties.
! The second list maps each of the blocks to an entry in the first list.
 INTEGER                                :: NBlockSym
 INTEGER, DIMENSION(:),allocatable      :: blockSymLs,blockSymId

! This is the threshold for determining the symmetry of individual geometries.
! If all coordinate changes are smaller than this value for a certain symmetry
! operation, we 
 DOUBLE PRECISION                       :: GeomSymT
!----SUBROUTINES--------------------------------------------
CONTAINS
 !***********************************************************************
 ! Generate all possible permutations of identical atoms
 ! Permutations of groups of permutationally identical atoms are first generated.
 ! They are then combined to generate the full list of permutations of all atoms.
 SUBROUTINE genAtomPerm()
  use progdata, only: natoms,atomList,atomCount,typeCount
  IMPLICIT NONE
  INTEGER                           ::  i, j, k
  INTEGER                           ::  maxSubPerm,maxCount
  INTEGER, DIMENSION(typeCount)     ::  indices
  CHARACTER(3)                      ::  tpStr

  if(allocated(nSubPerm))deallocate(nSubPerm)
  if(allocated(pmtList))deallocate(pmtList)
  if(allocated(subPerm))deallocate(subPerm)
  allocate(nSubPerm(typeCount))
  do i=1, typeCount
    nSubPerm(i)=factl(atomCount(i))
  end do
  ! determine the number of total permutations and allocate array space
  maxSubPerm=maxval(nSubPerm)
  maxCount=maxval(atomCount)
  nPmt= product(nSubPerm)
  if(printlvl>0)print *,"           nPmt=",nPmt
  allocate(pmtList(nPmt,natoms))
  allocate(subPerm(typeCount,maxSubPerm,maxCount))

  ! generate permutation list for each type of atoms
  do i=1,typeCount
     subPerm(i,:nSubPerm(i),:atomCount(i))=Permutation(atomCount(i),nSubPerm(i),atomList(i,1:atomCount(i)))
  end do

 ! construct total permutations list using subgroup lists
  indices(1:typeCount)=1
  write (tpStr,"(I3)") natoms
  do i=1,nPmt
   do j=1,typeCount
     do k=1, atomCount(j)
       pmtList(i,atomList(j,k))=subPerm(j,indices(j),k)
     end do !k
   end do !j
   indices(1)=indices(1)+1
   do j=2,typeCount
     if(indices(j-1)>nSubPerm(j-1)) then
       indices(j-1)=1
       indices(j)=indices(j)+1
     end if
   end do !j
  end do !i
 ! deallocate arrays
 END SUBROUTINE genAtomPerm

 !***********************************************************************
 ! This subroutine removes unfeasible permutations of atoms from the atom
 ! permutation list.
 !
 ! METHOD
 !   The selection procedure proceeds by defining connectivity relations between
 ! atoms. Permutations that preserves all connectivity conditions are considered
 ! feasible.
 !
 ! USES
 !   cntfl    [progdata] CHARACTER*30
 !          Input file that contains connectivity information
 ! VARIABLES
 !   removed  [out] LOGICAL
 !            True when some permutations are removed.
 ! FUTURE CHANGES
 !   METHOD NOT IMPLEMENTED YET
 !   ADD SUPPORT FOR PRESERVATION OF CHIRALITY OF OPTICAL CENTERS
  SUBROUTINE selectAtomPerm(removed)
    use progdata, only: cntfl,natoms
    IMPLICIT NONE
    LOGICAL,intent(out):: removed
    integer,parameter  :: CNTUNIT=2077
    ! connectivity restrictions used as permutation selection rules
    INTEGER                                :: nBonds
    INTEGER,dimension(:,:), allocatable    :: BondList
    ! masks for selection
    logical ::  feasible, found
    integer ::  ios,i,j,k, newBond(2), nPmtNew
    integer,allocatable :: newList(:,:)

    allocate(newList(nPmt,natoms))
    nPmtAll = nPmt
    removed =.false.
    ! skip selection if the file is not specified
    if(len_trim(cntfl)==0) return

    open(unit=CNTUNIT,file=trim(adjustl(cntfl)),access='sequential',&
        form='formatted',STATUS='OLD',ACTION='READ',POSITION='REWIND',IOSTAT=ios)
    if(ios/=0)then
        print *,"  Failed to open connectivity file."
        print *,"  Skipping feasible permutation selections."
        return
    end if

    ! read connectivity rules
    read(unit=CNTUNIT,fmt=*,IOSTAT=ios) nBonds
    if(ios/=0)then
        print *,"  Failed to read the number of connectivity rules."
        print *,"  Skipping feasible permutation selections."
        return
    end if
    if(allocated(BondList))deallocate(BondList)
    allocate(BondList(2,nBonds))
    if(nBonds>0)read(unit=CNTUNIT,fmt=*)BondList
    close(CNTUNIT)

    ! make sure all bonding relations are in cannonical order
    do i=1,nBonds
        if(BondList(1,i)>BondList(2,i))BondList(:,i)=BondList([2,1],i)
    end do!i=1,nbonds

    ! process rules to find out feasible permutations.  
    ! mask labels the feasible permutations in permutation list.
    nPmtNew = 0
    do i=1,nPmt
        ! Test if all the bonding rules form a group under the current
        ! permutation relation.   
        ! To-Do:  Also verify the preservation of chirality
        feasible =  .true.
        do j=1,nBonds
            newBond(1) = minval(pmtList(i,BondList(:,j)))
            newBond(2) = maxval(pmtList(i,BondList(:,j)))
            !search the bond list for the permutated bond
            found = .false.
            do k=1,nBonds
                if( newBond(1)==Bondlist(1,k) .and. &
                    newBond(2)==Bondlist(2,k))found=.true.
            end do
            if(.not.found)then
                feasible=.false.
                exit
            end if
        end do!j=1,nBonds
        if(feasible)then
            nPmtNew=nPmtNew+1
            newList(nPmtNew,:)=pmtList(i,:)
        end if!feasible
    end do!i=1,nPmt
    deallocate(BondList)

    ! reform the permutation list
    deallocate(pmtlist)
    nPmtAll = nPmt
    nPmt = nPmtNew
    allocate(pmtList(nPmt,natoms))
    pmtList = newList(:nPmt,:)
    deallocate(newList)

    ! log the selected permutations
    print *,"     Number of feasible permutations:",nPmt
    print *,"     Feasible Permutations List:"
    do i=1,nPmt
        write(unit=*,fmt="(8x,A,I5,A)",Advance="NO") "Pmt #",i,":"
        print "(20I3)",pmtList(i,:)
    end do!i=1,nPmt
    print *,""
    if(nPmtAll>nPmt)removed=.true.
  END SUBROUTINE selectAtomPerm

 !***********************************************************************
 !genCoordPerm() generates permutation list for all scaled coordinates
 !-----------------------
 !* Parities of coordinates are stored in coordPerm(0,:), which is 
 !1 for all scaled/simple Rij and -1 for any OOP coordinates.
 !* Permutations are stored in coordPerm(1:,:) with the corresponding 
 !signs of permuted term located in sgnCPerm(1:,:)
 !-----------------------
 !Note that parities are signs may be different, although they are the same
 !in the case of 4 atoms.  Parities are regarding inversion and signs are 
 !regarding permutations. A coordinate can have different signs with respect
 !to different permutations. 

 SUBROUTINE genCoordPerm()

   use progdata, only: CoordSet,coordmap

   IMPLICIT NONE

   integer :: i,j,m,n,k,sgn
   integer :: newCoord(4)  ! the coordinate generated by a permutation

! initialize structures
   if(allocated(coordPerm))deallocate(coordPerm)
   if(allocated(sgnCPerm))deallocate(sgnCPerm)
   allocate(coordPerm(0:nPmt,ncoord))
   allocate(sgnCPerm(nPmt,ncoord))

! determine parities with respect to inversion, and store it in unit 0
! of the permutation list
   do j=1,ncoord
     m=coordmap(j,1) !index of the coordinate set of the current coord
     select case(CoordSet(m)%Type)  !it is an OOP coord
! Rij, bond-bending and linear bendings are symmetric
       case (0,1,-3)
         coordPerm(0,j)=1
! Out-of-plane angle and torsion angles are antisymmetric
       case (-1,2,-2)
         coordPerm(0,j)=-1
     end select ! case(CoordSet(m)%Type)
   end do!j=1,ncoord

!determine effect of permutations on each of the scaled coordinates
   do i=1,nPmt
     do j=1,ncoord
       m=coordmap(j,1) !index of the coordinate set of the current coord
       n=coordmap(j,2) !index of the current coord in the coordinate set

! permutations are generated by different method for each type of coordinates
       select case(CoordSet(m)%Type)

       ! Rij
         case (0)
        ! generate permuted distances.  it will be in the same coordinate set.
           newCoord(1) = pmtList(i,CoordSet(m)%coord(1,n))
           newCoord(2) = pmtList(i,CoordSet(m)%coord(2,n))
           if(newCoord(1)>newCoord(2))then
                k           =  newCoord(2)
                newCoord(2) =  newCoord(1)
                newCoord(1) =  k
           end if
           do k=1,CoordSet(m)%ncoord
           ! find the permuted pair in the coordinate list
             if(count(newCoord(1:2).eq.CoordSet(m)%coord(1:2,k)).eq.2)then
               coordPerm(i,j) = CoordSet(m)%icoord(k)
               sgnCPerm(i,j)  = 1
               exit
             end if ! count(...)==2
           end do !k

       ! OOP
         case (-1,-2)
         ! generate the permuted OOP reference atom numbers in cannonical order
            if(CoordSet(m)%Type==-1)then
                CALL reorderOOP(pmtList(i,CoordSet(m)%coord(:,n)),newCoord, sgn)
            else
                CALL reorderOOP2(pmtList(i,CoordSet(m)%coord(:,n)),newCoord, sgn)
            end if
            do k=1,CoordSet(m)%ncoord
             if(all(newCoord.eq.CoordSet(m)%coord(:,k)))then
               coordPerm(i,j) = CoordSet(m)%icoord(k)
               sgnCPerm(i,j)  = sgn
               exit
             end if ! count(...)==4
            end do ! k

       ! 4 center dot product
         case (-3)
           call reorderDot4(pmtList(i,CoordSet(m)%coord(:,n)),newCoord,sgn)
           do k=1,CoordSet(m)%ncoord
             if(all(newCoord.eq.CoordSet(m)%coord(:,k)))then
               coordPerm(i,j) = CoordSet(m)%icoord(k)
               sgnCPerm(i,j)  = sgn
               exit
             end if ! count(...)==4
           end do ! k

       ! bond angle
         case (1)
        ! generate permuted angles
            newCoord(1) = minval(pmtList(i,CoordSet(m)%coord(1:2,n)))
            newCoord(2) = maxval(pmtList(i,CoordSet(m)%coord(1:2,n)))
            newCoord(3) = pmtList(i,CoordSet(m)%coord(3,n))
            do k=1,CoordSet(m)%ncoord
            ! find the permuted pair in the coordinate list
                if(count(newCoord(1:3).eq.CoordSet(m)%coord(:,k)).eq.3)then
                    coordPerm(i,j) = CoordSet(m)%icoord(k)
                    sgnCPerm(i,j)  = 1
                    exit
                end if ! count(...)==2
            end do !k

        !antisymmetric bend
         case(2)
           call reorderABend(pmtList(i,CoordSet(m)%coord(:,n)),newCoord,sgn)
           do k=1,CoordSet(m)%ncoord
             if(all(newCoord.eq.CoordSet(m)%coord(:,k)))then
               coordPerm(i,j) = CoordSet(m)%icoord(k)
               sgnCPerm(i,j)  = sgn
               exit
             end if ! count(...)==4
           end do ! k

       ! Others are not implemented yet
         case default
           print *," Type = ",CoordSet(m)%Type
           stop "Unsupported coordinate type in genCoordPerm"
       end select !case(CoordSet(m)%Type)
     end do!j=1,ncoord
     if(ncoord>30)then
       print "(3x,A,I2,A,30I4)","PMT[",i,"]",coordPerm(i,1:30)*sgnCPerm(i,1:30)
       print "(10x,30I4)",coordPerm(i,31:)*sgnCPerm(i,31:)
     else
       print "(3x,A,I2,A,30I4)","PMT[",i,"]",coordPerm(i,:)*sgnCPerm(i,:)
     end if
   end do!i=1,nPmt
   
 END SUBROUTINE genCoordPerm

 !***********************************************************************
 !PartTerms partition all terms into permutation cycles using genPermCycle
 ! result is stored in termCycle
 !The field val on termList is used temporarily to identify terms that are
 !already classified into subgroups (Cycles)
 SUBROUTINE PartitionTerms()
  use hddata, only: termList
  implicit none

  INTEGER                       :: i,j,count1,count2,count_rate
  type(TTermDef),pointer        :: pT
  integer                       :: nCycles,nTerms
  
  call deallocPCycle()
  allocate(permCycle(0:order)) !initialize the structure
  call system_clock(COUNT=count1,COUNT_RATE=count_rate)
  do i=0,order
   allocate(permCycle(i)%handle)
   permCycle(i)%last=>permCycle(i)%handle
   nTerms=termList(i)%nTerms
   nCycles=0
   pT=>TermList(i)%handle
   do j=1,nTerms
     pT=>pT%pNext
     pT%val=1
   end do!j=1,TermList(i)%nTerms
   pT=>TermList(i)%handle
   do j=1,nTerms
     pT=>pT%pNext
     if(pT%val<0)cycle
     CALL genPermCycle(pT)
     nCycles=nCycles+1
   end do!j=1,TermList(i)%nTerms
   if(printlvl>1)print  1000, nTerms, nCycles
  end do!i=1,order
  call system_clock(COUNT=count2)
  print 1001,dble(count2-count1)/count_rate
 1000 format("  ",7X,I9," terms partitioned into ",I8," cycles.")
 1001 format("  Terms classified into subgroups (Cycles) in ",f6.2," seconds")
 CONTAINS
   !genPermCycle will generate permutation Cycle for a specific term
   !The subroutine will located pointed by pT, then locate each of its
   !permutations by searching up from the 0th order term in termList
   !using icTerm links.
   !The val field if termList entries are used for marks. Unclassified
   !terms will have val=1.  -1 is the mark for used term, and -2 is the
   !mark for first term of current circle.
   !Permutations leading to the starting term will be added into U and W
   !depending on their sgns.
   !Permutations leading to an unclassified term will be added to array V
   !and the new term will be added to pC%term.
   !Temporary arrays of size nPmt are used before the sizes are determined.
   SUBROUTINE genPermCycle(pT)
    use hddata, only: termList
    implicit none
    TYPE(TTermDef),POINTER :: pT
    type(TTermDef),pointer         :: pTPmt
    type(TPermCycle),pointer       :: pC    !pointer to the new cycle
    integer                        :: i,j,sgn
    type(pTermDef),dimension(nPmt) :: tmpTerm
    integer,dimension(nPmt)        :: tmpU,tmpV,tmpW,tmpSgn
    if(pT%val<0)return
    pT%val=-2       !mark the starting term
    allocate(pC)    !create the new term
    tmpV(1)=1
    tmpTerm(1)%p=>pT
    tmpSgn(1)=1
    do i=1,nPmt
      pTPmt=>termList(0)%handle%pNext
      sgn=1
      do j=1,pT%ord
        if(.not. associated(pTPmt%icTerm(coordPerm(i,pT%coord(j)))%p))then
          print *,"Pmt=",coordPerm(i,:),", term=",pT%coord 
          stop'genPermCycle: Permuted term does not exist'
        end if
        pTPmt=>pTPmt%icTerm(coordPerm(i,pT%coord(j)))%p
        sgn=sgn*sgnCPerm(i,pT%coord(j))
      end do!j=1,pT%ord
      if(pTPmt%val>0)then        !new term encountered
        pC%nTerms=pC%nTerms+1
        tmpV(pC%nTerms)=i
        tmpTerm(pC%nTerms)%p=>pTPmt
        pTPmt%val=-1             !mark the term as used
        tmpSgn(pC%nTerms)=sgn
      elseif(pTPmt%val<-1.5)then !the first term encounted
        if(sgn>0)then
          pC%nW=pC%nW+1
          tmpW(pC%nW)=i
        else!if(sgn>0)
          pC%nU=pC%nU+1
          tmpU(pC%nU)=i
        end if!(sgn>0)
      end if!if(pT%val>0
    end do!i=1,nPmt
    allocate(pC%W(pC%nW))
    allocate(pC%U(pC%nU))
    allocate(pC%V(pC%nTerms))
    allocate(pC%term(pC%nTerms))
    allocate(pC%sgnTerm(pC%nTerms))
    pC%W=tmpW(1:pC%nW)
    pC%U=tmpU(1:pC%nU)
    pC%V=tmpV(1:pC%nTerms)
    pC%term=tmpTerm(1:pC%nTerms)
    pC%sgnTerm=tmpSgn(1:pC%nTerms)
    do j=1,pT%ord
      pC%parity=pC%parity*coordPerm(0,pT%coord(j))
    end do!j=1,pT%ord
    nullify(pC%pNext)
    pT%val=-1
    !link the new cycle to the list
    permCycle(pT%ord)%last%pNext=>pC
    permCycle(pT%ord)%last=>pC
    permCycle(pT%ord)%nCycle=permCycle(pT%ord)%nCycle+1
   END SUBROUTINE genPermCycle
 END SUBROUTINE PartitionTerms

 !***********************************************************************
 !Generating the list of symmetrized polynomial basis for each block
 subroutine genMaptab()
  use hddata, only:nl, nr, RowGrp, ColGrp, nBasis, lnBlock,CpOrder
  use progdata
  implicit none
  INTEGER                          ::  ord,iCyc,ntotal,iblk,broot,i
  type(TPermCycle),pointer         ::  pC
  if(printlvl>0)print *,"    Entering genMaptab..."
  ntotal=0
  do ord=0,order
    if(printlvl>1)print "(A,10X,A,I3)","  ","Order ",ord
 !construct blocks with unique symmetry 
    do i=1,NBlockSym 
      iblk = blockSymLs(i)
      if((RowGrp(iblk).ne.ColGrp(iblk)).and.ord>CpOrder)cycle
      pC=>permCycle(ord)%handle
      do iCyc=1,permCycle(ord)%nCycle
        pC=>pC%pNext
        CALL MatProj(ord,iblk,pC,nl(iblk)*nr(iblk))
      end do !iCyc=1,permCycle(ord)%nCycle
      if(printlvl>1)print 1000,iblk, nBasis(ord,iblk)
    end do !i=1,NBlockSym
  end do !do ord=1,order
 !copy basis to symmetry blocks with repeating symmetry
  if(printlvl>1)print *,"     Cloning basis for blocks with nonunique symmetry"
  do iblk=1,nblks
    broot = blockSymLs(blockSymId(iblk))
    if(broot.ne.iblk) call lnBlock(broot,iblk)
    do ord=0,order
      ntotal = ntotal + nBasis(ord,iblk)
    end do
  end do!iblk
  if(printlvl>0)print *,"    Maptab generated. nBasis=",ntotal
 1000 format(("  ",12X,"blk",I3," has ",I9," matrices"))
 END SUBROUTINE !genMapTab

 !Generate Matrix Projection element g^ij,lr_n
 SUBROUTINE MatProj(n,m,pC,LR)
    USE hddata
    implicit none
    INTEGER, INTENT(IN)                                 ::  n,m,LR
    TYPE(TPermCycle),POINTER                            ::  pC
    double precision,dimension(LR,LR)    ::  XT    ! Transpose of matrix X defined as
                                                   !X=Sum[D^l_ik(T)*D^r_jl(T),{T/in W}]
    double precision                     ::  x ! temp var for element of X matrix
    INTEGER                              ::  i,j,k,l,s,t,nProj,lirr,rirr,ll,rr
    DOUBLE PRECISION,DIMENSION(LR,LR)    ::  TAU     ! Dummy variabls for QR decomposition
    INTEGER,DIMENSION(LR)                ::  JPVT
    INTEGER                              ::  LWORK,info
    DOUBLE PRECISION,DIMENSION(LR*5+1)   ::  WORK
    double precision,dimension(LR,pC%nterms)  ::  coef
    double precision,dimension(LR,pC%nterms,LR*nSymLineUps)  ::  coefSave  ! all added coefficients
    double precision,external :: dnrm2
    double precision   :: ovlp, xproj(LR)
    integer   ::  iAdd
    integer   ::  iSym  ! index of symmetry setup

    ll=nl(m)
    rr=nr(m)
    if(LR/=ll*rr)stop"MatProj: inconsistent LR value"
    iAdd = 0
    do iSym=1,nSymLineUps
      ! check the parity of the current symmetry line up
      ! skip to the next line-up if parity is wrong
      if(pC%parity.ne.GrpPrty(RowGrp(m),iSym)*GrpPrty(ColGrp(m),iSym))cycle
      lirr=GrpSym(RowGrp(m),iSym)
      rirr=GrpSym(ColGrp(m),iSym)
      !initializations
      !Construct X matrix X=Sum[D^l_ik(T)*D^r_jl(T),{T/in W}]
      do i=1,ll
       do j=1,rr
        do k=1,ll
          do l=1,rr
            x=sum(irrep(lirr)%RepMat(pC%W,i,k)*irrep(rirr)%RepMat(pC%W,j,l))&
             -sum(irrep(lirr)%RepMat(pC%U,i,k)*irrep(rirr)%RepMat(pC%U,j,l))
            x=x/nPmt
            CALL chop(x,IrrepMatCutOff)
            XT(k*rr-rr+l,i*rr-rr+j)=x
          end do !l=1,rr
        end do !k=1,ll
       end do !j=1,rr
      end do !i=1,ll
      !if X matrix is 0 then the block is 0 by symmetry
      if(sum(XT*XT)<MatProjCutOff**2)return
      ! Make linearly independent components of matrix X
      nProj=1  !norm(X)>0, at least 1 independent component exists
      if(LR>1)then !if dim>1, linearly independent rows of X are constructed
       LWORK=LR*5+1
       JPVT=0
       CALL DGEQP3(LR,LR,XT,LR,JPVT, TAU, WORK, LWORK, INFO)
       if(INFO/=0)stop"MatProj: Failed to decomposite XT"
       ! Filter out rows whose diagonal element is higher than the threshold.
       do i=2,LR
         if(XT(i,i)>MatProjCutOff)then
           nProj=nProj+1  ! Remove the transformation info output by DGEQP3
           XT(i,1:i-1)=dble(0)  !   and only keep the upper triangle matrix R.
         else   !!if(XT(i,i)>MatProjCutOff
           exit
         end if !if(XT(i,i)>MatProjCutOff
       end do !i=2,LR
      end if !if(LR>1)
      ! fill in the output array
      do t=1,nProj
        if(LR>1)then
          xproj=0d0
          do i=1,LR
            if(JPVT(i)>0)  xproj(JPVT(i))=XT(t,i)
          end do
        else
          xproj(1)=xt(1,1)
        end if
        coef=dble(0)
        do i=1,ll
         do j=1,rr
           do s=1,pC%nTerms
            do k=1,ll
               coef(i+(j-1)*ll,s)=coef(i+(j-1)*ll,s)+ irrep(lirr)%RepMat(pC%V(s),i,k)* &
                    dot_product(xproj(k*rr-rr+1:k*rr),irrep(rirr)%RepMat(pC%V(s),j,:))*pC%sgnTerm(s)
            end do ! k
           end do !do s=1,pC%nterms
         end do !do j=1,rr
        end do !do i=1,ll
        coef = coef / dnrm2(LR*pC%nTerms,coef,1)
        ! verifies if this new basis overlaps with an older basis
        ovlp = 0d0
        do i=1,iAdd
           ovlp = ovlp+sum(coef*coefSave(:,:,i))**2
        end do
        if(ovlp>9.5d-1)cycle !this function is already present. skipping
        CALL add2maptab(n,m,pC%nterms,pC%term,coef)
        iAdd = iAdd+1
        coefSave(:,:,iAdd) = coef
      end do !do t=1,nProj
    end do !iSym
 END SUBROUTINE MatProj
 
 ! Release the memory occupied by entries of permCycle
 SUBROUTINE deallocPCycle()
   IMPLICIT NONE
   integer  :: i,j,status
   type(TPermCycle),pointer ::  pC,pNext
   if(allocated(permCycle))then
     do i=1,order
       pC=>permCycle(i)%handle%pNext
       if(associated(pC))then
         do j=1,permCycle(i)%nCycle
           if(associated(pC))then
             if(allocated(pC%U))deallocate(pC%U)
             if(allocated(pC%W))deallocate(pC%W)
             if(allocated(pC%term))deallocate(pC%term)
             if(allocated(pC%V))deallocate(pC%V)
             pNext=>pC%pNext
             deallocate(pC)
             pC=>pNext
           else 
             print *,"deallocPCycle: unexpected chain end at order:",i,",cycle:",j
           end if!(associated(pC))
         end do !j=1,permCycle(i)%nCycle
         deallocate(permCycle(i)%handle,STAT=status)
       end if
       if(status/=0)print *,"ERROR: failed to deallocate permCycle(",i,")%pCycle. STAT=",status
     end do!i=1,order
     deallocate(permCycle,STAT=status)
     if(status/=0)print *,"ERROR: failed to deallocate permCycle. STAT=",status
   end if!(allocated(permCycle))
 END SUBROUTINE deallocPCycle

 ! release memory occupied by irreducible representations
 SUBROUTINE deallocIrreps()
   IMPLICIT NONE
   integer :: i
   if(allocated(irrep))then
     do i=1,nirrep
        if(allocated(irrep(i)%RepMat))deallocate(irrep(i)%RepMat)
     end do!i=1,nirrep
     deallocate(irrep)
   end if!(allocated(irrep))
 END SUBROUTINE deallocIrreps

! this subroutine identifies the blocks that have unique symmetry and the ones
! that have identical symmetry properties.   
! This is complicated by the existence of multiple symmetry line-ups. The
! subroutine maintains a list of representation matrices of blocks in terms of
! permutation 
 SUBROUTINE getBlockSym()
   USE hddata, only: nblks,nstates,RowGrp,ColGrp
   IMPLICIT NONE
   ! variable definitions
   integer                                      :: NBlkIrrep
   type(TBlkRepMat),dimension(nirrep*nirrep)    :: BlkIrrep
   integer :: i,j,k,id
   integer,dimension(nblks,NSymLineUps)   :: BPerm,BPrty
   integer,dimension(nstates*nstates)     ::tmpBlkLs

   if(allocated(blockSymLs))deallocate(blockSymLs)
   if(allocated(blockSymId))deallocate(blockSymId)
   allocate(blockSymId(nblks))
   ! get all the permutational irrep mats for blocks
   NBlkIrrep = 0
   do i=1,nSymLineUps
     do j=1,nblks
       call addBlkIrrep(GrpSym(RowGrp(j),i),GrpSym(ColGrp(j),i),NBlkIrrep,BlkIrrep)
     end do
   end do!i=1,nSymLineUps

   if(printlvl>1)then
     print "(A,I4)","   Number of Possible Representations for Blocks:",NBlkIrrep
     do i=1,NBlkIrrep
       print "(A,I4)","   Matrices for Representation #",i
       do j=1,nPmt
         print "(A,I4)","     Permutation operation #",j 
         do k=1,BlkIrrep(i)%ldim*BlkIrrep(i)%rdim
           print "(6x,50F10.7)",BlkIrrep(i)%RepMat(j,k,:)
         end do!k
       end do!j
       print *,"     Parent pairs:"
       do j=1,BlkIrrep(i)%NParents
         print "(6x,'L:',I3,'   R:',I3)",BlkIrrep(i)%parents(j,:)
       end do!j
    end do!i
  end if!printlvl>1

  ! get the irrep of each block under each symmetry line-up
  do i=1,nblks
    do j=1,nSymLineUps
      call seekBlockIrrep(GrpSym(RowGrp(i),j),GrpSym(ColGrp(i),j),BPerm(i,j),&
                NBlkIrrep,BlkIrrep)
      BPrty(i,j) = GrpPrty(RowGrp(i),j)*GrpPrty(ColGrp(i),j)
    end do!j
  end do!i

  call sortBlkSymLists(nblks,nSymLineUps,BPerm,BPrty)

   NBlockSym = 0
   do i=1,nblks
     id=0
     do j=1,NBlockSym
       if(any(BPerm(tmpBlkLs(j),:).ne.BPerm(i,:)))cycle
       if(any(BPrty(tmpBlkLs(j),:).ne.BPrty(i,:)))cycle
       id=j
       exit    
     end do
     if(id==0)then
       NBlockSym=NBlockSym+1
       tmpBlkLs(NBlockSym) = i
       print *,"  Block symmetry # ",j
       print "(A,10I4)","     Permutation irrep  : ",BPerm(i,:)
       print "(A,10I4)","     Inversion symmetry : ",BPrty(i,:)
       id = j
       id = NBlockSym
     end if!(id.ne.0)then
     blockSymId(i) = id    
   end do

   allocate(blockSymLs(NBlockSym))
   blockSymLs = tmpBlkLs(1:NBlockSym)

   if(printlvl>1)then
     print "(A,I4)","    Number of existing block symmetries :",NBlockSym
     print *,"    Symmetry to block mappings:" 
     print "(6x,20I3)",blockSymLs
     print *,"    Block to symmetry mappings:" 
     print "(6x,20I3)",blockSymId
   end if

   do i=1,NBlkIrrep
     deallocate(BlkIrrep(i)%RepMat)
   end do
 END SUBROUTINE getBlockSym
!--------------------------------------------------------------------------
! Sort irrep matrix list for each block to canonical order to enable comparison
 SUBROUTINE sortBlkSymLists(nblk,nsym,BPerm,BPrty)
   IMPLICIT NONE
   INTEGER,intent(IN)             :: nblk,nsym
   INTEGER,dimension(nblks,nsym)  :: BPerm,BPrty
   integer ::  i,j,k, tmpPerm,tmpPrty
   do i=1,nblk
   ! repeat any repeating items to 0
     do j=1,nsym-1
       if(BPerm(i,j)==0)cycle
       do k=j+1,nsym
         if(BPerm(i,j)==BPerm(i,k).and.BPrty(i,j)==BPrty(i,k))then
           BPerm(i,k) = 0
           BPrty(i,k) = 0
         end if
       end do!k
     end do!j
   ! sort the list
     do j=1,nsym-1
       do k=nsym,j+1,-1
         if(BPerm(i,k)*BPrty(i,k)>BPerm(i,k-1)*BPrty(i,k-1))then
           tmpPerm = BPerm(i,k) 
           tmpPrty = BPrty(i,k) 
           BPerm(i,k) = BPerm(i,k-1)
           BPrty(i,k) = BPrty(i,k-1)
           BPerm(i,k-1) = tmpPerm 
           BPrty(i,k-1) = tmpPrty
         end if!
       end do!k
     end do!j
   end do
 END SUBROUTINE sortBlkSymLists
!--------------------------------------------------------------------------
! Seek the block irrep matrix table for a pair of parent states
 SUBROUTINE seekBlockIrrep(lstate,rstate,ind,NBlkIrrep,BlkIrrep)
   IMPLICIT NONE
   INTEGER, intent(IN)                          :: lstate,rstate
   INTEGER,intent(IN)                           :: NBlkIrrep
   type(TBlkRepMat),dimension(*),intent(IN)     :: BlkIrrep
   INTEGER, intent(OUT)                         :: ind
   integer  ::  i,j
   ind = 0
   do i=1,NBlkIrrep
     do j=1,BlkIrrep(i)%NParents
       if(BlkIrrep(i)%parents(j,1).ne.lstate)cycle
       if(BlkIrrep(i)%parents(j,2).ne.rstate)cycle
       ind = i
       return
     end do !j=1,BlkIrrep(i)%NParents
   end do !i=1,NBlkIrrep
 END SUBROUTINE seekBlockIrrep
!--------------------------------------------------------------------------
! generate block irrep matrix using state irrep matrix.  Seek the group irreps
! list and if it already exists, add the state pair to its parent list.   If the
! irrep matrix is new, create a new entry.
!
! Arguments
! ---------
! lstate,rstate [in] INTEGER
!               Index of irrep for the left and right parent state.
! NBlkIrrep     [in/out] INTEGER
!               Number of unique irrep matrices for symmetry blocks.
! BlkIrrep      [in/out] TBlkRepMat,dimension(*)
!               List of irrep matrices for the symmetry blocks of Hd
 SUBROUTINE addBlkIrrep(lstate,rstate,NBlkIrrep,BlkIrrep)
   IMPLICIT NONE
   INTEGER,intent(IN)           :: lstate,rstate
   INTEGER,intent(INOUT)        :: NBlkIrrep
   type(TBlkRepMat),dimension(*):: BlkIrrep
   
   integer  :: i,j,k,s1,s2,l,r,ll,rr,id
   double precision,dimension(:,:,:),allocatable  ::  RepMat 
  
   call seekBlockIrrep(lstate,rstate,id,NBlkIrrep,BlkIrrep)
   if(id.ne.0)return  !l,r pair already in the list. exiting
 
   ll = irrep(lstate)%dim
   rr = irrep(rstate)%dim
   allocate(RepMat(nPmt,ll*rr,ll*rr))
  
   ! construct the irrep matrices
   do i=1,nPmt
    l=0
    do j=1,ll
     do k=1,rr 
      l=l+1
      r=0
      do s1=1,ll
       do s2=1,rr
        r=r+1
        RepMat(i,l,r) = irrep(lstate)%RepMat(i,j,s1)*irrep(rstate)%RepMat(i,k,s2)
       end do !s2=1,rr
      end do!s1=1,ll
     end do! k=1,rr 
    end do!j=1,ll
   end do
   ! seek existing irrep table for this new sets of matrices
   id = 0
   do i=1,NBlkIrrep 
     if(BlkIrrep(i)%ldim.ne.ll .or. BlkIrrep(i)%rdim.ne.rr)cycle
     if(maxval((BlkIrrep(i)%RepMat-RepMat)**2)>MatProjCutoff**2)cycle
     id = i
     exit
   end do
   ! check if repmat is new
   if(id==0)then  ! it is a new representation
     NBlkIrrep = NBlkIrrep+1
     BlkIrrep(NBlkIrrep)%ldim = ll
     BlkIrrep(NBlkIrrep)%rdim = rr
     BlkIrrep(NBlkIrrep)%NParents = 1
     BlkIrrep(NBlkIrrep)%parents(1,1) = lstate
     BlkIrrep(NBlkIrrep)%parents(1,2) = rstate
     allocate(BlkIrrep(NBlkIrrep)%RepMat(nPmt,ll*rr,ll*rr))
     BlkIrrep(NBlkIrrep)%RepMat = RepMat
   else!if(id==0) ! located as old
     BlkIrrep(id)%NParents = BlkIrrep(id)%NParents+1
     if(BlkIrrep(id)%NParents>MaxParents)  &
        stop "addBlkIrrep: number of parent state pairs exceeding the limit."
     BlkIrrep(id)%parents(BlkIrrep(id)%NParents,1)=lstate
     BlkIrrep(id)%parents(BlkIrrep(id)%NParents,2)=rstate
   end if!(id==0)then
   
   deallocate(RepMat)
 END SUBROUTINE addBlkIrrep

! This subroutine checks a geometry (in internal coordinates) and returns all the permutation
!-inversions that keep the geometry invariant. 
 SUBROUTINE CheckGeomSym(igeom,nsym,sympmt,syminv)
   IMPLICIT NONE
   DOUBLE PRECISION,dimension(ncoord),intent(IN)    :: igeom
   INTEGER,intent(OUT)                              :: nsym
   INTEGER,dimension(:),allocatable,intent(OUT)     :: sympmt
   INTEGER,dimension(:),allocatable,intent(OUT)     :: syminv
   integer :: i
   double precision :: newcoord(ncoord)
   integer :: sym_pmt(npmt*2),sym_inv(npmt*2)
   nsym=0
   do i=1,npmt
     newcoord(coordPerm(i,:))=igeom*sgnCPerm(i,:)
     if(maxval(abs(newcoord-igeom))<GeomSymT)then
       nsym=nsym+1
       sym_pmt(nsym)=i
       sym_inv(nsym)=1
     end if
     newcoord=newcoord*coordPerm(0,:)
     if(maxval(abs(newcoord-igeom))<GeomSymT)then
       nsym=nsym+1
       sym_pmt(nsym)=i
       sym_inv(nsym)=-1
     end if
   end do 
   allocate(sympmt(nsym))
   sympmt=sym_pmt(1:nsym)
   allocate(syminv(nsym))
   syminv=sym_inv(1:nsym)
 END SUBROUTINE

END MODULE CNPI 
