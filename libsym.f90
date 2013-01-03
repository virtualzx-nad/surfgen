!********************************************************************
! Module contains information regarding the CNPI symmetry
! of the molecule and polynomial
!********************************************************************
MODULE CNPI

 use hddata, only: order,nblks,pTermDef,TTermDef,ncoord
 use progdata, only: printlvl,T2DList
 use combinatorial

 !----DERIVED TYPES------------------------------------------
 ! derived type for permutation cycle
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
   integer                                  ::  nCycle=0
   type(TPermCycle),pointer                 ::  handle
   type(TPermCycle),pointer                 ::  last
 end type

 ! derived type for irredicible representations.
 ! Dim is dimensionality of irrep and Order is order of the group
 ! RepMat contains all representation matrices
 type TIrrep
   INTEGER                                        ::  Dim
   INTEGER                                        ::  Order
   DOUBLE PRECISION,DIMENSION(:,:,:),allocatable  ::  RepMat
 end type

!----CONSTANTS------------------------------------------
 DOUBLE PRECISION,PARAMETER    :: IrrepMatCutOff=1D-7
 DOUBLE PRECISION,PARAMETER    :: MatProjCutOff=1D-5

!----GLOBAL VARIABLES-----------------------------------
 ! * SYMMETRY PROPERTIES *
 ! array that contains all the relevant irreducibles
 type(TIrrep),DIMENSION(:),allocatable        ::  irrep
 integer                                      ::  nirrep=0

 INTEGER,dimension(:),allocatable             ::  grpSym
 INTEGER,dimension(:),allocatable             ::  grpPrty

! * PERMUTATION PROPERTIES
 INTEGER                                :: nPmt
 INTEGER, DIMENSION(:,:), allocatable   :: pmtList !all permutation of atoms
 INTEGER, DIMENSION(:,:,:),allocatable  :: subPerm
 INTEGER, DIMENSION(:),allocatable      :: nSubPerm

 INTEGER, DIMENSION(:,:), allocatable   :: coordPerm
 INTEGER, DIMENSION(:,:), allocatable   :: sgnCPerm

! term definitions
 type(TPCycleList),dimension(:),allocatable   ::  permCycle

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
  if(printlvl>1) print *,"           Atom permutation list"
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
   if(printlvl>1)write(*,"(A,10X,A,"//trim(tpStr)//"I3,A)")  &
                    "  ","[",pmtList(i,:),"  ]"
  end do !i
 ! deallocate arrays
 END SUBROUTINE genAtomPerm

 !***********************************************************************
 ! This subroutine removes unfeasible permutations of atoms from the atom
 ! permutation list.
 !
 ! METHOD
 !   The filter proceeds by defining connectivity relations between atoms.
 ! Permutations that preserves all connectivities are considered feasible.
 !
 ! ARGUMENT
 !   cntfl  [in] CHARACTER*30
 !          Input file that contains connectivity information
 ! FUTURE CHANGES
 !   METHOD NOT IMPLEMENTED YET
 !   ADD SUPPORT FOR PRESERVATION OF CHIRALITY OF OPTICAL CENTERS
  SUBROUTINE filterAtomPerm()
    IMPLICIT NONE

  END SUBROUTINE filterAtomPerm

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

   integer :: i,j,m,n,k,l
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
       case (0,1)
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
           newCoord(1) = minval(pmtList(i,CoordSet(m)%coord(n,1:2)))
           newCoord(2) = maxval(pmtList(i,CoordSet(m)%coord(n,1:2)))
           do k=1,CoordSet(m)%ncoord
           ! find the permuted pair in the coordinate list
             if(count(newCoord(1:2).eq.CoordSet(m)%coord(k,:)).eq.2)then
               coordPerm(i,j) = CoordSet(m)%icoord(k)
               sgnCPerm(i,j)  = 1
               exit
             end if ! count(...)==2
           end do !k

       ! OOP
         case (-1,-2)
         ! generate the permuted OOP reference atom numbers in cannonical order
           CALL reorderOOP(pmtList(i,CoordSet(m)%coord(n,:)),newCoord, sgnCPerm(i,j))
           do k=1,CoordSet(m)%ncoord
             if(all(newCoord.eq.CoordSet(m)%coord(k,:)))then
               coordPerm(i,j) = CoordSet(m)%icoord(k)
               exit
             end if ! count(...)==4
           end do ! k

       ! Others are not implemented yet
         case default
           print *," Type = ",CoordSet(m)%Type
           stop "Unsupported coordinate type in genCoordPerm"
       end select !case(CoordSet(m)%Type)
     end do!j=1,ncoord
     PRINT "(3x,A,I2,A,50I3)","PMT[",i,"]",coordPerm(i,:)
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
    use progdata, only:  natoms
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
      end do!j=1,pT%order
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
  use hddata, only:nl, nr, RowGrp, ColGrp, nBasis
  use progdata
  implicit none
  INTEGER                          ::  ord,iCyc,i,ntotal
  type(TPermCycle),pointer         ::  pC
  if(printlvl>0)print *,"    starting genMaptab..."
  ntotal=0
  do ord=0,order
    if(printlvl>1)print "(A,10X,A,I3)","  ","Order ",ord
    do i=1,nblks
      pC=>permCycle(ord)%handle
      do iCyc=1,permCycle(ord)%nCycle
        pC=>pC%pNext
        if(pC%parity==grpPrty(RowGrp(i))*grpPrty(ColGrp(i))) &
            CALL MatProj(ord,i,pC,nl(i)*nr(i))
      end do !iCyc=1,permCycle(ord)%nCycle
      if(printlvl>1)print 1000,i, nBasis(ord,i)
      ntotal=ntotal+nBasis(ord,i)
    end do !i=0,nblks
  end do !do ord=1,order
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
    ll=nl(m)
    rr=nr(m)
    if(LR/=ll*rr)stop"MatProj: inconsistent LR value"
    lirr=grpSym(RowGrp(m))
    rirr=grpSym(ColGrp(m))
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
      coef=dble(0)
      do i=1,ll
       do j=1,rr
         do s=1,pC%nTerms
          do k=1,ll
             coef(i+(j-1)*ll,s)=coef(i+(j-1)*ll,s)+ irrep(lirr)%RepMat(pC%V(s),i,k)* &
                  sum(XT(t,k*rr-rr+1:k*rr)*irrep(rirr)%RepMat(pC%V(s),j,:))*pC%sgnTerm(s)
          end do ! k
         end do !do s=1,pC%nterms
       end do !do j=1,rr
      end do !do i=1,ll
      CALL add2maptab(n,m,pC%nterms,pC%term,coef)
    end do !do t=1,nProj
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
END MODULE CNPI 
