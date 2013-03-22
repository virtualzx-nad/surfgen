MODULE makesurfdata
  use hddata, only: T3DDList,T2DDList,TDList
  use progdata, only: abpoint
  IMPLICIT NONE

  ! derived type to identify a certain block of a point
  TYPE tblk
    INTEGER :: point
    INTEGER :: i,j
  END TYPE

 ! derived type for eq lists
  TYPE TEqList
   type(tblk),dimension(:),allocatable         :: List
   INTEGER                                     :: length
  END TYPE

  DOUBLE PRECISION                             :: eshift    ! uniform shift on ab initio energies
  DOUBLE PRECISION                             :: deg_cap   ! threshold below which intersection adapted coordinates will be used
  INTEGER                                      :: maxiter
  INTEGER                                      :: mmiter    ! maximum number of micro iterations
  INTEGER                                      :: npoints
  DOUBLE PRECISION                             :: toler     !criteria for convergence
  DOUBLE PRECISION                             :: gcutoff   !gradient cutoff
  DOUBLE PRECISION                             :: exactTol  !Exact equations rank cutoff
  DOUBLE PRECISION                             :: LSETol    !LSE rank cutoff
  DOUBLE PRECISION                             :: flattening!flattening parameter for differential convergence
!gorder: energy difference below which ckl will be ordered according to gradients rather than energies
  DOUBLE PRECISION                             :: gorder
  CHARACTER(72)                                :: outputfl  !name of the output file generated after fit
  CHARACTER(72)                                :: flheader  !comment line in the outputfile
  CHARACTER(72)                                :: ckl_input   !if nonempty, read wave functions from this file
  CHARACTER(72)                                :: ckl_output  !if nonempty, output wave functions to this file

  LOGICAL                                      :: usefij   ! use coupling instead of coupling*dE

! weights for equations of reproduction of energy, energy gradients, and non-adiabatic coupling
  DOUBLE PRECISION                             :: w_energy,w_grad,w_fij 
  DOUBLE PRECISION                             :: nrmediff, ediffcutoff, nrmediff2, ediffcutoff2
! ndiis   : size of diis space
! ndstart : number of iterations before starting DIIS acceleration for wavefunctions
  INTEGER                                      :: ndiis,ndstart

! rmsexcl : specify the threshold for point weight below which the point will be excluded from RMS analysis
!           same format as points.in input
  INTEGER                                      :: rmsexcl

  TYPE(abpoint),dimension(:),allocatable       :: dispgeoms
  type(TEqList)                                :: exclEner,exclGrad,exactEner,exactGrad,exactDiff,enfGO
  INTEGER                                      :: enfDiab ! index of point where diabatic and adiabatic matches
  DOUBLE PRECISION,dimension(:),allocatable    :: ptWeights

  type(TDList),dimension(:),allocatable        :: WVals
  type(TDList),dimension(:),allocatable        :: dWVals    
  type(T2DDList),dimension(:),allocatable      :: WMat
  DOUBLE PRECISION,dimension(:,:,:),allocatable:: ckl

! when followPrev is true, eigenvectors will be phased and ordered to match previous iteration instead of ab initio data
  LOGICAL                                      :: followPrev 

! The following parameters control the construction of local coordinates
! Translational and rotational coordinates of the whole system will be removed, as well as those of the
! dissociated pieces if threshold and decaying parameter for the scaling functions are set properly.

! useInternal  :   transform ab initio gradients to an orthogonal local internal coordinates defined as
! eigenvectors of B^T.B, where B is Wilson's B matrix.
! intGradT     :   Threshold for eigenvalue of B^T.B higher than which the coordinate is considered internal
! intGradS     :   Threshold for eigenvalue of B^T.B lower than which the coordinate will be scaled with a
!                  factor    eval/intGrad,   where eval is the eigenvalue
! gScaleMode:     =0   Do not scale
!                 >0   Scale all coordinates
!                 <0   Scale couplings only
  LOGICAL                                      :: useIntGrad
  DOUBLE PRECISION                             :: intGradT,intGradS 
  INTEGER                                      :: gScaleMode

! nvibs:  number of coordinates along which gradients will be taken.
!  3*natoms when useIntGrad=.FALSE.,  3*natoms-5 when useIntGrad=.TRUE. 
! Individual data point will have a number of internal degree of freedom
! equal to or lower than nvibs, stored in the field abpoint%nvibs
  INTEGER                                      :: nvibs 

! These parameters control the automatic downscaling of high energy points
! EnergyT      Threshold energy above which gradients/energies will be scaled down
! HighEScale   The scaling factor for these high energy data 
  DOUBLE PRECISION,dimension(10)               ::  energyT, highEScale

! These variables store the coefficient matrix and right hand side of
! linear equality constrained least squares equations.
  DOUBLE PRECISION,dimension(:,:),allocatable  ::  NEL
  DOUBLE PRECISION,dimension(:)  ,allocatable  ::  rhs

! maximum allowed change in coefficients in each iteration
  DOUBLE PRECISION  :: maxd, maxED

! scaling factors to exact equations
  DOUBLE PRECISION  :: scaleEx

! scaling factor for the gradients for gradient based methods
  DOUBLE PRECISION  :: gScaler

! scaling factors to exact equations
  logical,dimension(:,:,:),allocatable       :: incgrad,g_exact
  logical,dimension(:,:,:),allocatable       :: incener,e_exact

! map of equations and coefficients
  integer,dimension(:,:),allocatable         :: EqMap
  integer,dimension(:,:),allocatable         :: coefMap

! basis preconditioning and contraction 
  integer,dimension(:),allocatable           :: npb  ! number of primitive basis for each block
  integer,dimension(:),allocatable           :: nbas ! number of reconstructed basis for each block
  double precision                           :: TBas ! theshold for eigenvalue cutoff of the primitive basis overlap matrix (linear dependency)
  double precision                           :: MaxRot  ! maximum rotation of eigenvector of Schrodingers equations at a point between iterations    
  double precision                           :: ecutoff  ! energy threshold above which data will be excluded from basis construction
  double precision                           :: egcutoff ! energy threshold above which gradients will no longer be used in fit
  type(T2DDList),dimension(:),allocatable    :: ZBas ! transformation from primitive to the reconstructed basis for each block
  type(T2DDList),dimension(:),allocatable    :: ZBasI! corresponding backwards transformation

! scaling factor of dij term.  0 gives old result and 1 gives exact gradients
  double precision            :: DijScale
  double precision            :: DijScale2

! index of the iteration to start differential convergence
  integer    :: dfstart

! method used to predict steps
  integer    :: stepMethod 

! tolerance for micro iterations convergence of exact equation error
  double precision   :: ExConv

! number of test point for gradient verifications
  integer    :: ntest

! 
! number of unknown coefficients
  integer    :: ncons
! number of least squares fitting equations
  integer    :: neqs
! number of exact equations
  integer    :: nex

! number of linear steps to search along positive and negative direction
  integer    :: linSteps, linNegSteps

!Hd predictions for all data points
  DOUBLE PRECISION,dimension(:,:,:,:),allocatable    :: fitG
  DOUBLE PRECISION,dimension(:,:),allocatable        :: fitE
 CONTAINS

  ! determine if each equation will be include / excluded / fitted exactly
  SUBROUTINE getPtList()
    use hddata, only: nstates
    use progdata, only: printlvl
    IMPLICIT NONE
    INTEGER i,s1,s2

    if(printlvl>0)print *,"Processing data inclusion/exclusion/exact fit data"
    if(allocated(incgrad))deallocate(incgrad)
    if(allocated(g_exact))deallocate(g_exact)
    if(allocated(incener))deallocate(incener)
    if(allocated(e_exact))deallocate(e_exact)

    allocate(incgrad(npoints,nstates,nstates))
    allocate(g_exact(npoints,nstates,nstates))
    allocate(incener(npoints,nstates,nstates))
    allocate(e_exact(npoints,nstates,nstates))

! Generate destiny table for each piece of ab initio data
    incgrad=.false.
    incener=.false.
    e_exact=.false.
    g_exact=.false.
    do i=1,npoints
      if(abs(ptWeights(i)-1)>1d-5.and.printlvl>1) &
         print "(A,I6,A,F9.4,A,I6)","point ",i," weight=",ptWeights(i), " nvibs=",dispgeoms(i)%nvibs
      if(ptWeights(i)>1D-8)then
        do s1=1,nstates
          do s2=1,s1
            incgrad(i,s1,s2)=.true.
          end do
          incener(i,s1,s1)=.true.
          if(i==enfDiab)incener(i,:,:)=.true.
        end do!s1=1,nstates
      end if
    end do
    do i=1,exclEner%length
      do s1=exclEner%List(i)%i,exclEner%List(i)%j
        do s2=exclEner%List(i)%i,exclEner%List(i)%j
          incener(exclEner%List(i)%point,s1,s2)=.false.
        end do
      end do
    end do
    do i=1,exclGrad%length
      incgrad(exclGrad%List(i)%point,exclGrad%List(i)%i,exclGrad%List(i)%j)=.false.
      incgrad(exclGrad%List(i)%point,exclGrad%List(i)%j,exclGrad%List(i)%i)=.false.
    end do
    do i=1,exactEner%length
     do s1=exactEner%List(i)%i,exactEner%List(i)%j
       do s2=exactEner%List(i)%i,exactEner%List(i)%j
          e_exact(exactEner%List(i)%point,s1,s2)=.true.
       end do
     end do
    end do
    do i=1,exactGrad%length
      g_exact(exactGrad%List(i)%point,exactGrad%List(i)%i,exactGrad%List(i)%j)=.true.
      g_exact(exactGrad%List(i)%point,exactGrad%List(i)%j,exactGrad%List(i)%i)=.true.
    end do

  END SUBROUTINE getPtList
!-----------------------------------------------------------------------------
!This subroutine reads in an array of points, then output a list of equations.
!In the list are the point ID/block ID/gradient direction of each equation.
!This can be used to generate either exact or least square equation list
!maxEqs          :  Largest possible number of equations
!nEqs            :  Number of equations actually generated
!lseMap(maxEqs,4):  The first nEqs rows record the (/point ID,state1,state2,grad ID/)
!                     of the equations generated. gradID=0 for energies
!nex             :  Number of exact equations actually generated
!exactEqs(maxEqs,4) :  Same as lseMap but those equations will be fitted exactly
  SUBROUTINE makeEqMap(maxEqs,lseMap,exactEqs,gradNorm,wvec)
    use hddata, only: nblks,blkmap,nstates
    use progdata, only: AU2CM1
    implicit none
    integer,intent(in)                      :: maxEqs
    integer,dimension(MaxEqs,4),intent(out) :: lseMap,exactEqs
    double precision,dimension(MaxEqs),intent(out)   :: wvec
    DOUBLE PRECISION,dimension(npoints,nvibs,nblks),intent(in)     :: gradNorm

    integer  :: i,j,k,s1,s2,nr,ng
    double precision :: residue,resDir,ediff
    logical :: needHeader

    character(4),dimension(0:2)                     ::  fatesymb 
    character(4),dimension(nstates)                 ::  enerfate
    character(7),dimension(nstates*(nstates+1)/2)   ::  gradfate
    character(3)                                    ::  str1, str2  ! # of ener and grad data


! Output point destiny table
    write(str1,"(I3)")  nstates
    write(str2,"(I3)")  nstates*(nstates+1)/2

    print *, "Destiny of ab initio data at each data point"
    print *, "  LS = Least Squares Fit"
    print *, "  EX = Excluded from Fit"
    print *, "  LM = Enforced using Lagrange multipliers"
    ng = 1
    do s1 = 1,nstates
      write(enerfate(s1),"(X,I2,X)") s1
      do s2=s1,nstates
        write(gradfate(ng),"(x,I2,',',I2,x)"),s1,s2
        ng= ng + 1
      end do !s2
    end do !s1
    print *," POINT ",enerfate,gradfate

    fatesymb(0) = " LS "
    fatesymb(1) = " EX "
    fatesymb(2) = " LM "
    do i=1, npoints
      ng = 1
      do s1 = 1,nstates
        enerfate(s1) = fatesymb(0)
        if(     e_exact(i,s1,s1))  enerfate(s1) = fatesymb(2)
        if(.not.incener(i,s1,s1))  enerfate(s1) = fatesymb(1)
        do s2=1,s1
          gradfate(ng)  =  " "//fatesymb(0)//"  "
          if(     g_exact(i,s1,s2))  gradfate(ng) = " "//fatesymb(2)//"  "
          if(.not.incgrad(i,s1,s2))  gradfate(ng) = " "//fatesymb(1)//"  "
          ng= ng + 1
        end do !s2
      end do !s1
      print "(2X,I5,2X,"//trim(adjustl(str1))//"A4,"//trim(adjustl(str2))//"A7)",& 
                  i,                  enerfate,                        gradfate
    end do !i

! make equations
    needHeader=.true.
    residue=0.
    nr=0
    nEqs=0
    nex=0
    do i=1,npoints
      if(ptWeights(i)<1D-8)cycle             ! Create equations for gradients and derivative couplings
      do s1=1,nstates
        do s2=1,nstates
          if(incgrad(i,s1,s2))then
            do j=1,dispgeoms(i)%nvibs
              if(blkMap(s2,s1)>0)then
               if(gradNorm(i,j,blkMap(s2,s1))>1D-8)then
                if(g_exact(i,s2,s1))then                                    ! Exact equations
                  nex=nex+1
                  if(nex>maxEqs)stop 'makeEqMap: nex>maxeqs'
                  exactEqs(nex,:)=(/i,s1,s2,j/)
                else !(g_exact(i,s2,s1))                                    ! Least squares equations
                  nEqs=nEqs+1
                  if(nEqs>maxEqs)stop 'makeEqMap: neqs>maxeqs'
                  lseMap(nEqs,:)=(/i,s1,s2,j/)
                  if(s1==s2)then                                            ! Energy gradient
                    if(gScaleMode==1)then
                      wvec(nEqs)=ptWeights(i)*w_grad*dispgeoms(i)%scale(j)
                    else
                      wvec(nEqs)=ptWeights(i)*w_grad
                    end if
                    k = 0 
                    do while(dispgeoms(i)%energy(s1)>energyT(k+1))  !  determine the bracket of current energy
                      k = k+1
                      if(k==10)exit
                    end do
                    if(k>0) wvec(nEqs) =  wvec(nEqs)*highEScale(k)
                  else  !(s1==s2)                                           ! Derivative couplings
                    if(gScaleMode==3)then
                      wvec(nEqs)=ptWeights(i)*w_fij*dispgeoms(i)%scale(j)**2
                    else if(gScaleMode>0.or.gScaleMode==-2)then
                      wvec(nEqs)=ptWeights(i)*w_fij*dispgeoms(i)%scale(j)
                    else
                      wvec(nEqs)=ptWeights(i)*w_fij
                    end if
                    k = 0
                    do while(dispgeoms(i)%energy(s1)+dispgeoms(i)%energy(s2)> 2*energyT(k+1))
                                 !  determine the bracket of current energy
                      k = k+1
                      if(k==10)exit
                    end do
                    if(k>0) wvec(nEqs) =  wvec(nEqs)*highEScale(k)
                    if(nrmediff>0)then
                      ediff=abs(dispgeoms(i)%energy(s1)-dispgeoms(i)%energy(s2))*AU2CM1
                      ediff=(ediff+ediffcutoff)/nrmediff
                      wvec(nEqs)=wvec(nEqs)/ediff
                    end if!(nrmediff>0)

                  end if!(s1==s2)
                end if !(g_exact(i,s2,s1))
               else !(gradNorm(i,j,blkMap(s2,s1))>1D-8)                     ! This gradient othorgonal to the basis.
                resDir=dispgeoms(i)%grads(j,s1,s2)                          ! This implies a symmetry zero or an inadequate
                residue=residue+resDir**2                                   ! basis expansion.
                nr=nr+1
                if(abs(resDir)>10*gcutoff)then                              ! Calculate the total residual ab initio gradients 
                  if(needHeader)then                                        ! along these othorgonal coordinates.
                    print *,"    Significant gradients zeroed by symmetry." ! If this residual gradient is large, either the 
                    needHeader=.false.                                      ! ab initio is symmetry broken, or the symmetry 
                  end if                                                    ! setting is incorrect, or basis is too small.
                  print 1000,i,s1,s2,j,resDir
                  print *,"  projection of bmat at point"
                  print 1001,dispgeoms(i)%bmat(:,j)
                end if!(abs(resDir)>10*gcutoff)
               end if!(gradNorm(i,j,blkMap(s2,s1))>0.)
              end if!(blkMap(s2,s1)>0)
            end do!j=1,dispgeom()%nvibs
          end if!(incgrad(i,s1,s2))
        end do !s2=1,s1
      end do !s1=1,nstates
    end do !i=1,npoints
    residue=sqrt(residue/nr)
    if(residue>gcutoff)print *,"   Norm of residue right-hand-side gradients: ",residue
    do i=1,exactDiff%length
        if(incener(exactDiff%List(i)%point,exactDiff%List(i)%i,exactDiff%List(i)%i).and.&
           incener(exactDiff%List(i)%point,exactDiff%List(i)%j,exactDiff%List(i)%j)) then
          nex=nex+1
          if(nex>maxEqs)stop 'makeEqMap: nex>maxeqs'
          exactEqs(nex,:)=(/exactDiff%List(i)%point,exactDiff%List(i)%i,&
                              -exactDiff%List(i)%j,0/)
          incener(exactDiff%List(i)%point,exactDiff%List(i)%i,exactDiff%List(i)%j)=.true.
        end if
    end do!i=1,exactDiff%length
    do i=1,npoints
        if(ptWeights(i)<1D-8)cycle
        do s1=1,nstates
          do s2=1,nstates
            if(incener(i,s1,s2))then
              if(e_exact(i,s1,s2))then                                      ! Exact equations
                nex=nex+1
                if(nex>maxEqs)stop 'makeEqMap: nex>maxeqs'
                exactEqs(nex,:)=(/i,s1,s2,0/)
              else                                                          ! Least squares equations
                nEqs=nEqs+1
                if(nEqs>maxEqs)stop 'makeEqMap: neqs>maxeqs'
                lseMap(nEqs,:)=(/i,s1,s2,0/)
                wvec(nEqs)=ptWeights(i)*w_energy
                k=0
                do while(dispgeoms(i)%energy(s1)>energyT(k+1))
                  k=k+1
                  if(k==10)exit
                end do 
                if(k>0)   wvec(nEqs) =  wvec(nEqs)*highEScale(k)
                if(s1==s2 .and. s2>0 .and. nrmediff2>0)then
                  if(s1>1)then
                    ediff=abs(dispgeoms(i)%energy(s1)-dispgeoms(i)%energy(s1-1))*AU2CM1
                    ediff=(ediff+ediffcutoff2)/nrmediff2
                    if(ediff<1D0)wvec(nEqs)=wvec(nEqs)/ediff
                  end if
                  if(s1<nstates)then
                    ediff=abs(dispgeoms(i)%energy(s1+1)-dispgeoms(i)%energy(s1))*AU2CM1
                    ediff=(ediff+ediffcutoff2)/nrmediff2
                    if(ediff<1D0)wvec(nEqs)=wvec(nEqs)/ediff
                  end if!(s1<nstates)
                end if!(s1==s2 .and. s2>0 .and. nrmediff2>0)
              end if
            end if
          end do !s2=1,nstates
        end do !s1=1,nstates
    end do!i=1,npoints
  1000 format(8X,"Point",I5,", St",2I3,", dir",I3,", grad=",E10.3)
  1001 format(8X,3F14.10)
  END SUBROUTINE
!-----------------------------------------------------------------------------------
!Make H and DH from Hd for a point and solve schrodinger equation to get eigenvectors.
!Diabatrize the states if degenerate.
!Reorder states by best fit of energy gradients if energy difference lower than gorder.
!Then determine phase of eigenvectors to best fit couplings.
  SUBROUTINE updateEigenVec(hvec,follow)
    USE hddata, only:nstates,EvaluateHd3
    USE progdata, only: abpoint,printlvl,AU2CM1
    USE combinatorial
    IMPLICIT NONE

    INTERFACE
      SUBROUTINE fixphase(nvibs,scale,fitgrad,abgrad,ckl,phaseList)
        use hddata, only: nstates
        IMPLICIT NONE
        INTEGER,intent(IN)                                        :: nvibs
        DOUBLE PRECISION,dimension(nvibs,nstates,nstates),intent(inout)           :: fitgrad
        double precision,dimension(nvibs,nstates,nstates),intent(in)              :: abgrad
        double precision,dimension(nvibs),intent(in)              :: scale
        DOUBLE PRECISION,dimension(nstates,nstates),INTENT(INOUT) :: ckl
        integer,dimension(2**nstates,nstates),INTENT(IN)          :: phaseList
      END SUBROUTINE
    END INTERFACE

    DOUBLE PRECISION,DIMENSION(ncons),intent(IN)     :: hvec   ! vector that defines Hd coefficients
    LOGICAL,INTENT(IN),OPTIONAL                      :: follow ! whether state signs and ordering will be
                                                               ! determined by following current ckl

    DOUBLE PRECISION,dimension(nstates,nstates)      :: hmatPT
    DOUBLE PRECISION,dimension(nvibs,nstates,nstates):: dhmatPT
    INTEGER                                          :: i,j,k,INFO
    DOUBLE PRECISION,dimension(5*nstates*nstates)    :: scr
    type(abpoint)                                    :: fitpt
    integer,dimension(2**(nstates-1),nstates)        :: phaseList
    logical                                          :: folPrev
    double precision,dimension(nstates,nstates)      :: cklPrev
    integer,dimension(:,:),allocatable               :: pmtList

    allocate(pmtList(factl(nstates),nstates))
    pmtList=Permutation(nstates,factl(nstates))
    folPrev = .false.
    if(present(follow))folPrev=follow
    if(printlvl>0) print *,"     Updating wave functions"
    allocate(fitpt%energy(nstates))
    allocate(fitpt%grads(nvibs,nstates,nstates))
    if(printlvl>1) print *,"      Making phase list"
    !making permutation list for grad ordering
    ! Generate list of all possible phase factors for fixphase
    phaseList(1,1:nstates)=1
    do i=2,2**(nstates-1)
      phaseList(i,:)=phaseList(i-1,:)
      j=2
      do while(phaseList(i,j)==-1)
        phaseList(i,j)=1
        j=j+1
      end do
      phaseList(i,j)=-1
    end do
    !generating eigenvectors and conform to intersection adapted coordinations/reorder
    do i=1,npoints
      cklPrev = ckl(i,:,:)
      CALL EvaluateHd3(hvec,nBas,npoints,i,nvibs,hmatPt,dhmatPt,WMat)
      ckl(i,:,:)=hmatPt

      call DSYEV('V','U',nstates,ckl(i,:,:),nstates,fitE(i,:),scr,5*nstates*nstates,INFO)
      IF(INFO/=0)then
        print *,"Failed to solve eigenvectors at point",i
        print *,"INFO=",info

        stop 'error solving eigenvalue equations'
      end if
      fitpt%energy=fitE(i,:)      
      if(i==enfDiab)then
        ckl(i,:,:)=dble(0)
        do k=1,nstates
          ckl(i,k,k)=dble(1)
        end do
        fitpt%ndeggrp = 0
      else
        call genEnerGroups(fitpt,deg_cap,nstates)
      end if!(i==enfDiab)
      call OrthGH_Hd(dispgeoms(i),dhmatPt(:dispgeoms(i)%nvibs,:,:),ckl(i,:,:),100,sqrt(gcutoff)/10)
      do k=1,dispgeoms(i)%nvibs
        fitpt%grads(k,:,:)=matmul(transpose(ckl(i,:,:)),matmul(dhmatPt(k,:,:),ckl(i,:,:)))
      end do!k=1,dispgeoms(i)%nvibs
      fitpt%energy=fitE(i,:)
      do k=1,enfGO%Length
        if(enfGO%List(k)%point==i)then
          do j=enfGO%List(k)%i+1,enfGO%List(k)%j
            fitpt%energy(j)=fitpt%energy(enfGO%List(k)%i)
          end do!j=enfGO%List(k)%i+1,enfGO%List(k)%j
        end if!(enfGO%List(k)%point==i)
      end do!k=1,enfGO%Length
      if(folPrev)then
! determine ordering and phase by consistency with previous iteration
        call folPrevCkl(ckl(i,:,:),CklPrev,pmtList,factl(nstates),nstates,i)
        do k=1,dispgeoms(i)%nvibs
          fitG(i,k,:,:)=matmul(transpose(ckl(i,:,:)),matmul(dhmatPt(k,:,:),ckl(i,:,:)))
        end do!k=1,dispgeoms(i)%nvibs
        fitG(i,dispgeoms(i)%nvibs+1:nvibs,:,:) = 0d0
      else
! determine ordering and phase by comparing with ab initio data
        call genEnerGroups(fitpt,gorder/AU2CM1,nstates)
        call gradOrder(i,fitpt,dispgeoms(i),ckl(i,:,:),pmtList,factl(nstates),w_energy,w_grad)
        call fixphase(dispgeoms(i)%nvibs,dispgeoms(i)%scale(:dispgeoms(i)%nvibs),fitpt%grads(1:dispgeoms(i)%nvibs,:,:),&
                 dispgeoms(i)%grads(1:dispgeoms(i)%nvibs,:,:),ckl(i,:,:),phaseList)
        fitpt%grads(dispgeoms(i)%nvibs+1:,:,:)=dble(0)
        fitG(i,:,:,:)=fitpt%grads
      end if !folPrev
      do k=1,nstates
        fitE(i,k)= dble(0)
        do j=1,nstates
          fitE(i,k)=fitE(i,k)+ckl(i,j,k)*dot_product(ckl(i,:,k),hmatPt(j,:))
        end do
      end do
    end do!i=1,npoints
    deallocate(fitpt%energy)
    deallocate(fitpt%grads)
    deallocate(fitpt%deg_groups)
  END SUBROUTINE
!-----------------------------------------------------------------------------------
! determine ordering and phase of Ckl by comparing with previous ckl
! cklnew   [in/out] DOUBLE PRECISION,dimension(nstates,nstates)
!          On input:  new eigenvectors of Hd before reordering
!          On output:  new eigenvectors of Hd with correct ordering and phases
! cklold   [in] DOUBLE PRECISION,dimension(nstates,nstates)
!          Eigenvector of a previous Hd.  New eigenvectors are first reordered to gain
!          maximum overlap with these old vectors, then adjusted to the same phase
! pmts     [in] INTEGER,dimension(nPmts,nstates)
!          List of all permutations that will be attempted by the reordering procedure
! nPmts    [in] INTEGER
!          Number of permutations
! nstates  [in] INTEGER
!          Number of electronic states
! PtId     [in] INTEGER
!          Point index for output purpose
  SUBROUTINE folPrevCkl(cklnew,cklold,pmts,npmts,nstates,ptId)
    USE progdata, ONLY:  printlvl
    IMPLICIT NONE
    INTEGER,INTENT(IN)                                        :: nPmts,nstates,PtId
    INTEGER,DIMENSION(nPmts,nstates),INTENT(IN)               :: pmts
    DOUBLE PRECISION,DIMENSION(nstates,nstates),INTENT(INOUT) :: cklnew
    DOUBLE PRECISION,DIMENSION(nstates,nstates),INTENT(IN)    :: cklold

    integer  ::  i,j
    double precision   ::  ototal,&    ! total overlap of a specific ordering scheme
                           omax,  &    ! maximum overlap so far
                           ovlp        ! overlap between individual vector
    double precision,dimension(nstates,nstates)::  cklmax  

!  first iterate through all possible permutations to find the one with maximum overlap
    omax = -1d0
    do i=1,nPmts
      ototal= sum([(abs(dot_product(cklnew(:,Pmts(i,j)),cklold(:,j))),j=1,nstates)])
      if(ototal>omax)then
        omax = ototal
        do j=1,nstates
          cklmax(:,j) = cklnew(:,Pmts(i,j))
        end do
      end if
    end do!i=1,nPmts
    cklnew = cklmax
! reset phase of wave functions to make them consistent with previous iteration
    do i=1,nstates
      ovlp = dot_product(cklnew(:,i),cklold(:,i))
      if(printlvl>1.and.ovlp**2<0.6d0)print "(A,I5,A,I2,A,F6.2)",&
                 "  State matched with small overlap at point ",ptid," state ",i,", ovlp=",ovlp
      if(ovlp<0) cklnew(:,i) =  -cklnew(:,i)
    end do !i=1,nstates
  END SUBROUTINE folPrevCkl
!-----------------------------------------------------------------------------------
!Generate the matrix A of the linear coefficients between Hd basis and RHS values
  SUBROUTINE makecoefmat(AMat)
    use hddata, only: nblks,nstates,RowGrp,ColGrp,offs,grpLen,nBasis
    use progdata, only: printlvl
    IMPLICIT NONE
    DOUBLE PRECISION,dimension(nex+neqs,ncons),INTENT(OUT)               :: AMat
    INTEGER      :: i,j,l,r,l0,r0,pt,s1,s2,g,as1,as2,iBlk,iBasis,ll,rr
    double precision  :: ediff,error
    double precision  :: wij(npoints,nstates,nstates),Dij(npoints,nstates,nstates)
    DOUBLE PRECISION,dimension(:,:,:,:),allocatable  :: VIJ
    double precision  :: cc(npoints)
    integer  :: count1,count2,count_rate
    INTEGER    ::  Id(npoints),Jd(npoints),ndeg,Pd(npoints),ind1
    DOUBLE PRECISION :: gvec(nvibs),hvec(nvibs)

! construct Wij, Dij and Vij matrices for all data points
    allocate(VIJ(npoints,nvibs,nstates,nstates))    

    if(printlvl>0)print *,"   Making equation coefficients."
    call system_clock(COUNT=count1,COUNT_RATE=count_rate)

    error = dble(0)
    j = 0
    do iBlk = 1,nblks
     l0 = offs(RowGrp(iblk))
     ll = GrpLen(RowGrp(iblk))
     r0 = offs(ColGrp(iblk))
     rr = GrpLen(ColGrp(iblk))
     do iBasis=1,nBas(iBlk)
      j=j+1
      WIJ = dble(0)
      VIJ = dble(0)
      DIJ = dble(0)
      ndeg = 0
      do s1=1,nstates
        do s2=s1,nstates
          do r= 1, rr
            do l= 1, ll
   ! calculate the multiplication of ckl, with Hermitianization for off diagonal blocks
              if(RowGrp(iblk).eq.ColGrp(iblk))then
                cc = ckl(:,l+l0,s1)*ckl(:,r+r0,s2)
              else !(RowGrp(iblk).eq.ColGrp(iblk))
                cc = ckl(:,l+l0,s1)*ckl(:,r+r0,s2) + ckl(:,r+r0,s1)*ckl(:,l+l0,s2)
              end if !(RowGrp(iblk).eq.ColGrp(iblk))

   ! get the values of W and V
              ind1 = ((r-1)*ll+l-1)*npoints*(nvibs+1)+1
              WIJ(:,s1,s2) = WIJ(:,s1,s2)+cc*WMat(iBlk)%List(ind1:ind1+(npoints-1)*(nvibs+1):nvibs+1,iBasis)
   ! construct the Hamiltonian contribution of VIJ 
              do i=1,nvibs
                ind1 = ind1+1
                VIJ(:,i,s1,s2) = VIJ(:,i,s1,s2) + cc * WMat(iBlk)%List(ind1:ind1+(npoints-1)*(nvibs+1):nvibs+1,iBasis)
              end do
            end do !  l= 1, ll
          end do !r= 1, rr

          if(s1.ne.s2) then
            if(abs(DijScale*DijScale2)>1D-30)then
  ! construct DIJ from WIJ.  It is non zero only on the off diagonal
              do i=1,npoints
                if(i==enfDiab)cycle !no rotation for reference point
                ediff = fitE(i,s2)-fitE(i,s1)
                if(abs(ediff)<deg_cap)then
                  DIJ(i,s1,s2)=dble(0)
                  ndeg = ndeg + 1
                  if(ndeg>npoints)STOP"TOO MANY DEGENERATE POINTS"
                  Id(ndeg)=s1
                  Jd(ndeg)=s2
                  Pd(ndeg)=i
                else!abs(ediff)<deg_cap
                  DIJ(i,s1,s2)=WIJ(i,s1,s2)/ediff*DijScale*DijScale2
                end if!abs(ediff)<deg_cap
              end do !i=1,npoints
            end if!
  ! fill the lower triangle of D and W
            WIJ(:,s2,s1)     =  WIJ(:,s1,s2)
            DIJ(:,s2,s1)     = -DIJ(:,s1,s2)
          end if !(J.ne.I)
  ! the Dij Contribution of VIJ
          do g=1,nvibs
            do l=1,nstates
              VIJ(:,g,s1,s2) = VIJ(:,g,s1,s2) + DIJ(:,l,s1)*fitG(:,g,l,s2) + DIJ(:,l,s2)*fitG(:,g,l,s1)     
            end do
          end do
          if(s1.ne.s2)VIJ(:,:,s2,s1)=VIJ(:,:,s1,s2)
        end do !s2
      end do ! s1
  ! construct DIJ contribution to degenerate points
      if(abs(DijScale*DijScale2)>1D-30)then
        do i=1,ndeg
          pt = Pd(i)
          s1 = Id(i)
          s2 = Jd(i) 
          if(abs(DIJ(pt,s1,s2)).gt.1D-10)STOP "DIJ CONSTRUCTION FOR DEG POINT INCONSISTENT"
          gvec = fitG(pt,:,s1,s1)-fitG(pt,:,s2,s2)
          hvec = fitG(pt,:,s1,s2)
          DIJ(pt,s1,s2)=(dot_product(VIJ(pt,:,s1,s1)-VIJ(pt,:,s2,s2),hvec)+dot_product(gvec,VIJ(pt,:,s1,s2)))&
                   /(4*dot_product(hvec,hvec)-dot_product(gvec,gvec))
          DIJ(pt,s2,s1)=-DIJ(pt,s1,s2)
          VIJ(pt,:,s1,s2) = VIJ(pt,:,s1,s2) + DIJ(pt,s1,s2)*gvec
          VIJ(pt,:,s2,s1) = VIJ(pt,:,s1,s2) 
          VIJ(pt,:,s1,s1) = VIJ(pt,:,s1,s1) - DIJ(pt,s1,s2)*hvec*2
          VIJ(pt,:,s2,s2) = VIJ(pt,:,s2,s2) + DIJ(pt,s1,s2)*hvec*2
        end do! I=1,ndeg
      end if !DijScale /= 0
! Make AMAT from W and V matrices
      do i=1,nex+neqs
        AMat(i,j)=dble(0)
        pt=eqMap(i,1)    ! point index
        s1=eqMap(i,2)    ! state index 1
        s2=eqMap(i,3)    ! state index 2
        g =eqMap(i,4)    ! gradient component index
        if(g==0)then !iGrad==0  => its an energy fit
          do l=1,GrpLen(RowGrp(iblk))
            do r=1,GrpLen(ColGrp(iblk))
              as1=l+offs(RowGrp(iblk))
              as2=r+offs(ColGrp(iblk))
              if(s2<0) then! equation for energy difference
                AMat(i,j)=WIJ(pt,s1,s1)-WIJ(pt,-s2,-s2)
              else         ! equation for energy 
                if(s1==s2)then
                  AMat(i,j)=WIJ(pt,s1,s2)
                else
                  AMat(i,j)=WIJ(pt,s1,s2)*(1-DijScale*DijScale2)
                end if 
              end if !s2<0
            end do!r
          end do!l
        else                                      ! equation for energy gradients or derivative couplings
          AMat(i,j)=VIJ(pt,g,s1,s2)
        end if!if(eqMap(i,4)==0)
      end do!i=1,neqs+nex
     end do! ibasis
    end do! iblk
    deallocate(VIJ)
    call system_clock(COUNT=count2)
    if(printlvl>0)print "(7X,I12,A,F7.2,A)",&
         ncons*(nex+neqs)," coefficients constructed after ",(count2-count1)/dble(count_rate),"s"
  end SUBROUTINE makecoefmat
!-----------------------------------------------------------------------------------
!generate rhs vectors for Hd fitting
  SUBROUTINE makebvec(bvec,diff)
    use progdata, only: printlvl
    IMPLICIT NONE
    DOUBLE PRECISION,dimension(nex+neqs),INTENT(INOUT)             :: bvec
    LOGICAL, INTENT(IN)                                            :: diff
    INTEGER                                                :: i,pt,s1,s2,g
    if(printlvl>0)print *,"   Making right hand side vectors."
    do i = 1,nex+neqs
      pt=eqmap(i,1)
      s1=eqmap(i,2)
      s2=eqmap(i,3)
      g =eqmap(i,4)
      if(g==0)then !is an energy
        if(s2>0)then !its adiabatic energy or adiabatic off diagonal terms(0)
          if(s1==s2)then
            bvec(i)=dispgeoms(pt)%energy(s1)
            if(diff)bvec(i)=bvec(i)-fitE(pt,s1)
          else 
            bvec(i)=dble(0)
          end if
        else  !it is adiabatic energy difference
          bvec(i)=dispgeoms(pt)%energy(s1)-dispgeoms(pt)%energy(-s2)
          if(diff)bvec(i)=bvec(i)-fitE(pt,s1)+fitE(pt,-s2)
        end if
      else !is gradient
        bvec(i)=dispgeoms(pt)%grads(g,s1,s2)
        if(diff)bvec(i)=bvec(i)-fitG(pt,g,s1,s2)
      end if!g==0
    enddo!i=1,nex+neqs
  end SUBROUTINE makebvec

 !calculate error of fit
 SUBROUTINE getError(nrmener,avgener,nrmgrad,avggrad)
   USE progdata, ONLY : printlvl
   USE hddata, ONLY: nstates
   IMPLICIT NONE
   DOUBLE PRECISION,intent(OUT)                   :: nrmener,avgener,nrmgrad,avggrad

   INTEGER  :: incdata,j,k,l,NEx_E,NEx_grd,nvpt,inc_e,inc_cp
   DOUBLE PRECISION  :: gnrm,dgrd,rmsexclthreshold,dE_exact,dG_exact,dcp,nrmcp,nrmdcp,ncp,de1,de2

   if(rmsexcl>=0)then
     rmsexclthreshold = 0D0
   else
     rmsexclthreshold = -1D0/rmsexcl
   end if
   if(printlvl>1)print *,"rmsexclthreshold = ",rmsexclthreshold
    nrmener = 0D0
    nrmgrad = 0D0
    avgener = 0D0
    avggrad = 0D0
    nrmcp  = 0d0
    nrmdcp = 0d0
    incdata = 0
    inc_e   = 0
    inc_cp  = 0
    NEx_E   = 0
    NEx_grd = 0
    dE_exact  = 0D0
    dG_exact  = 0D0
    
    do j=1,npoints
       nvpt = dispgeoms(j)%nvibs
       if(ptWeights(j)>rmsexclthreshold)then
       do k = 1,nstates
         do l=k,nstates
           if(l>k.and.(incgrad(j,k,l).or.incgrad(j,l,k)))then
             de1 = max(abs(dispgeoms(j)%energy(k)-dispgeoms(j)%energy(l)),1D-7 )
             de2 = max(abs(fitE(j,k)-fitE(j,l)),1D-7 )
             inc_cp=inc_cp+1
             dcp    =  dot_product(  &
                   dispgeoms(j)%grads(:nvpt,k,l)/de1-fitG(j,:nvpt,k,l)/de2, &
                 ( dispgeoms(j)%grads(:nvpt,k,l)/de1-fitG(j,:nvpt,k,l)/de2)*dispgeoms(j)%scale(1:nvpt)  )
             ncp = dot_product(dispgeoms(j)%grads(:nvpt,k,l)                          , &
                               dispgeoms(j)%grads(:nvpt,k,l)*dispgeoms(j)%scale(:nvpt) ) /de1**2
             if(dcp>1d-1*ncp.and.printlvl>2.and.ncp>1d0 ) print "(5x,A,I5,A,F10.3,A,E12.4)", &
                       "Large coupling error at point ",j," : ",sqrt(dcp/ncp)*100,"% out of ", sqrt(ncp)
             nrmdcp = nrmdcp+dcp
             nrmcp  = nrmcp +ncp
           end if! coupling included
           if(g_exact(j,k,l))then
             dgrd    =  dot_product(  &
                   dispgeoms(j)%grads(:nvpt,k,l)-fitG(j,:nvpt,k,l), &
                 ( dispgeoms(j)%grads(:nvpt,k,l)-fitG(j,:nvpt,k,l))*dispgeoms(j)%scale(1:nvpt)  )
             dG_exact=dG_exact+dgrd
             NEx_grd=NEx_grd+1
           end if
         end do!l
         if(e_exact(j,k,k))then
           dE_exact= dE_exact + (dispgeoms(j)%energy(k)-fitE(j,k))**2
           NEx_e = NEx_e+1
         end if!e_exact(j,k,k)
         if(dispgeoms(j)%energy(k)<energyT(1))then
           nrmener = nrmener +    (dispgeoms(j)%energy(k)-fitE(j,k))**2
           avgener = avgener + abs(dispgeoms(j)%energy(k)-fitE(j,k))
           inc_e=inc_e+1
           dgrd    =  dot_product(  &
                dispgeoms(j)%grads(:nvpt,k,k)-fitG(j,:nvpt,k,k), &
              ( dispgeoms(j)%grads(:nvpt,k,k)-fitG(j,:nvpt,k,k) )*dispgeoms(j)%scale(1:nvpt)  )
           gnrm = dot_product(dispgeoms(j)%grads(:nvpt,k,k),dispgeoms(j)%grads(:nvpt,k,k)*dispgeoms(j)%scale(1:nvpt))
           dgrd = dgrd / gnrm 
           if((dgrd>5D-2.or.gnrm*dgrd>1D-5.and.gnrm<1D-4).and.printlvl>2) print "(5x,A,I5,A,F10.3,A,E12.4)", &
                       "Large gradient error at point ",j," : ",dgrd*100,"% out of ", gnrm
           if(gnrm>1D-4) then
             nrmgrad = nrmgrad + dgrd
             avggrad = avggrad + sqrt(dgrd)
             incdata = incdata + 1
           else
              if(printlvl>4)print "(5x,A,I5,A,E14.4,A,E14.4)",&
                       "Small gradients excluded at point ",j,".  Norm of grad =",gnrm,",Norm of error = ",dgrd*gnrm
           end if ! (gnrm>1D-4)
         end if
       enddo
     end if
    end do
    if(printlvl>1)print *,"    ",incdata," point/states included in gradient RMS analysis"
    if(printlvl>1)print *,"    ",inc_e  ," point/states included in enery RMS analysis"
    nrmener = Sqrt(nrmener/inc_e)
    nrmgrad = Sqrt(nrmgrad/incdata)
    avgener = avgener/inc_e
    avggrad = avggrad/incdata
    if(printlvl>1)print *,"    RMS error to exact energy equations:    ",Sqrt(dE_exact/NEx_e)
    if(printlvl>1)print *,"    RMS error to exact gradient equations:  ",Sqrt(dG_exact/NEx_grd)
    if(printlvl>1)print "(2(A,E12.4))","    RMS error of derivative couplings:   ",&
                        Sqrt(nrmdcp/inc_cp),"  RMS norm of couplings:",Sqrt(nrmcp/inc_cp)
 END SUBROUTINE

  !----------------------------------------------------------------------------------------
  !Calculate numerical gradients with respect to the fitting coefficients
  SUBROUTINE getCGrad(nex,neqs,nc,cilm,dCi,dLambda,lag,weight,jaco)
    USE progdata, ONLY : printlvl
    USE hddata,   ONLY : RowGrp, ColGrp, offs, GrpLen, nstates,nBasis,nBlks
    IMPLICIT NONE
    INTEGER, INTENT(IN)                                                  :: nex, neqs, nc
    DOUBLE PRECISION,dimension(nex+nc), INTENT(IN)                       :: cilm
    DOUBLE PRECISION,dimension(neqs),INTENT(IN)                          :: weight
    DOUBLE PRECISION,dimension(nc ), INTENT(OUT)                         :: dCi
    DOUBLE PRECISION,dimension(nex), INTENT(OUT)                         :: dLambda
    DOUBLE PRECISION,                INTENT(OUT)                         :: Lag    
    DOUBLE PRECISION,dimension(nex,nc), INTENT(OUT)                      :: jaco
    
!  These intermediate quantities are constructed for each coefficient i as vectors over
!  all the points, because they are used by multiple equations that involve the same point.
!  W^IJ_i =<f^I|W_i|f^J>
!  D^IJ_i = W^IJ_i / (E^J-E^I)
!  V^IJ_i =<f^I|del W_i|f^J>
!  U^IJ_i = V^IJ_i + sum_K[D^KI_i*h^KJ+D^KJ_i*h^KI]
    DOUBLE PRECISION,dimension(:,:,:),allocatable   ::  WIJ,  DIJ
    DOUBLE PRECISION,dimension(:,:,:,:),allocatable ::  VIJ,  UIJ

! PmQ is a vector that contain differences between desired values and fit predictions,
! for exact and least squares equations in the first and second section respectively, 
! PmQ = (P[ls] - Q[ls])          index nex+1 to nex+neqs
! PmQ = P[ex] - Q[ex]              index 1 to nex
    DOUBLE PRECISION,dimension(neqs+nex)     ::  PmQ

! dQ is the derivatives of the Hd predictions of the quantities being fitted with respect
! to the current coefficient
    DOUBLE PRECISION,dimension(neqs+nex)     ::  dQ

!  ico : coefficient index.    ieq : equation index.   ipt : point index.   igd : grad index
    INTEGER    ::  ico, ieq, ipt, igd

!  temporary vector for basis values and gradients and for multiplication of two C^kl 
    DOUBLE PRECISION  ::  cc(npoints)
!  timing variables
    DOUBLE PRECISION  ::  t1,t2,t3,t4,t5,t6
    INTEGER           ::  count0,count1,count2,count_rate
!  variables that specify coefficient
    INTEGER    ::  iBlk,iBasis
!  variables that specify equation
    INTEGER    ::  pt, s1, s2, g,ind1
    INTEGER    ::  I,J,K,L,K0,nK, L0,nL,Id(npoints),Jd(npoints),ndeg,Pd(npoints)
!  other local variables
    DOUBLE PRECISION :: EIJ,diffsq
!  used to store g and h vectors for degenerate case
    DOUBLE PRECISION :: gvec(nvibs),hvec(nvibs)
!  external BLAS subroutine
    double precision,external :: ddot
    
    if(printlvl>0) print *,"   Entering getCGrad.  Exact coefficient derivatives will be evaluated."
    call system_clock(COUNT=count0,COUNT_RATE=count_rate)

! Initialization of parameters
    allocate(WIJ(npoints,nstates,nstates))
    allocate(VIJ(npoints,nvibs,nstates,nstates))
    allocate(UIJ(npoints,nvibs,nstates,nstates))
    allocate(DIJ(npoints,nstates,nstates))

    dCi     = dble(0)
    dLambda = dble(0)
    PmQ     = dble(0)
    diffsq  = dble(0)

! calculate difference vectors PmQls and PmQex
    if(printlvl>2) print *,"     Generating difference vectors."
    do ieq = 1,nex + neqs
      pt=eqmap(ieq,1)
      s1=eqmap(ieq,2)
      s2=eqmap(ieq,3)
      g =eqmap(ieq,4)
      if(g==0)then !is an energy
       if(s2>0)then !its adiabatic energy or adiabatic off diagonal terms(0)
         if(s1==s2)then
           PmQ(ieq) = dispgeoms(pt)%energy(s1) - fitE(pt,s1)
         else
           PmQ(ieq) = dble(0)
         end if
       else  !(s2>0)    it is adiabatic energy difference
         PmQ(ieq) = dispgeoms(pt)%energy(s1)-dispgeoms(pt)%energy(-s2)-fitE(pt,s1)+fitE(pt,-s2)
       end if !(s2>0)
     else !(g==0)   is gradient
       PmQ(ieq) = dispgeoms(pt)%grads(g,s1,s2) - fitG(pt,g,s1,s2)
     end if!(g==0)
    end do !ieq
! <<<--- done calculating PmQ
    call system_clock(COUNT=count2,COUNT_RATE=count_rate)
    t1 = dble(count2-count0)/count_rate
    t2 = 0
    t3 = 0
    t4 = 0
    t5 = 0
    t6 = 0

    if(printlvl>2) print *,"     Calculating derivatives of Lagrangian."
    ico = 0
    do iBlk=1,nblks
     do iBasis=1,nBas(iBlk)
      ico = ico+1
      K0 = offs(RowGrp(iblk))
      nK = GrpLen(RowGrp(iblk))
      L0 = offs(ColGrp(iblk))
      nL = GrpLen(ColGrp(iblk))
      UIJ = dble(0)
      WIJ = dble(0)
      VIJ = dble(0)
      DIJ = dble(0)
      ndeg = 0
   ! construct vectors W,V and D for the current coefficient
      do I=1,nstates
        do J=I,nstates
          do L= 1, nL
            do K= 1, nK
              call system_clock(COUNT=count1,COUNT_RATE=count_rate)      
   ! calculate the multiplication of ckl, with Hermitianization for off diagonal blocks
              if(RowGrp(iblk).eq.ColGrp(iblk))then 
                cc = ckl(:,K+K0,I)*ckl(:,L+L0,J)
              else !(RowGrp(iblk).eq.ColGrp(iblk))
                cc = ckl(:,K+K0,I)*ckl(:,L+L0,J) + ckl(:,L+L0,I)*ckl(:,K+K0,J)
              end if !(RowGrp(iblk).eq.ColGrp(iblk))
              call system_clock(COUNT=count2,COUNT_RATE=count_rate)      
              t2 = t2 + dble(count2-count1)/count_rate

   ! get the values of W and V
              ind1 = ((L-1)*nK+K-1)*npoints*(nvibs+1)+1
              WIJ(:,I,J) = WIJ(:,I,J)+cc*WMat(iBlk)%List(ind1:ind1+(npoints-1)*(nvibs+1):nvibs+1,iBasis)
   ! construct the Hamiltonian contribution of VIJ 
              do igd=1,nvibs
                ind1 = ind1+1
                VIJ(:,igd,I,J) = VIJ(:,igd,I,J) + cc * WMat(iBlk)%List(ind1:ind1+(npoints-1)*(nvibs+1):nvibs+1,iBasis)
              end do!igd
              call system_clock(COUNT=count1,COUNT_RATE=count_rate)      
              t3 = t3 + dble(count1-count2)/count_rate
            end do !  K= 1, nK
          end do !L= 1, nL

          call system_clock(COUNT=count1,COUNT_RATE=count_rate)
          if(J.ne.I) then 
   ! construct DIJ from WIJ.  It is non zero only on the off diagonal
            do ipt=1,npoints
              if (ipt==enfDiab) cycle
              EIJ = fitE(ipt,J)-fitE(ipt,I)
              if(abs(EIJ).le.deg_cap)then 
                DIJ(ipt,I,J)=0
                ndeg = ndeg + 1
                if(ndeg>npoints)STOP"TOO MANY DEGENERATE POINTS"
                Id(ndeg)=I
                Jd(ndeg)=J
                Pd(ndeg)=ipt
              else
                DIJ(ipt,I,J)=WIJ(ipt,I,J)/EIJ*DijScale
              end if
            end do !ipt=1,npoints


   ! fill the lower triangle of D and W
            WIJ(:,J,I)     = WIJ(:,I,J)
            DIJ(:,J,I)     = -DIJ(:,I,J)
            VIJ(:,:,J,I)   = VIJ(:,:,I,J)
          end if !(J.ne.I)
   ! calculate the second piece of VIJ with DIJ and hIJ(=fitG)
          do igd = 1,nvibs
            UIJ(:,igd,I,J) = VIJ(:,igd,I,J)
            do K = 1,nstates
              UIJ(:,igd,I,J) = UIJ(:,igd,I,J) + DIJ(:,K,I)*fitG(:,igd,K,J) + DIJ(:,K,J)*fitG(:,igd,K,I)
            end do ! K = 1,nstates
            if(I.ne.J)  UIJ(:,igd,J,I) = UIJ(:,igd,I,J)  !complete lower triangle
          end do!igd=1,nvibs 
          call system_clock(COUNT=count2,COUNT_RATE=count_rate)
          t3 = t3 + dble(count2-count1)/count_rate
        end do !J=I,nstates
      end do !I=1,nstates
      do ipt=1,ndeg
        pt = Pd(ipt)
        I  = Id(ipt)
        J  = Jd(ipt) 
        if(abs(DIJ(pt,I,J)).gt.1D-10)STOP "DIJ CONSTRUCTION FOR DEG POINT INCONSISTENT"
        gvec = fitG(pt,:,I,I)-fitG(pt,:,J,J)
        hvec = fitG(pt,:,I,J)
        DIJ(pt,I,J)=(dot_product(UIJ(pt,:,I,I)-UIJ(pt,:,J,J),hvec)+dot_product(gvec,UIJ(pt,:,I,J)))&
                   /(4*dot_product(hvec,hvec)-dot_product(gvec,gvec))
        DIJ(pt,J,I)=-DIJ(pt,I,J)
        UIJ(pt,:,I,J) = UIJ(pt,:,I,J) + DIJ(pt,I,J)*gvec
        UIJ(pt,:,J,I) = UIJ(pt,:,I,J) 
        UIJ(pt,:,I,I) = UIJ(pt,:,I,I) - DIJ(pt,I,J)*hvec*2
        UIJ(pt,:,J,J) = UIJ(pt,:,J,J) + DIJ(pt,I,J)*hvec*2
      end do! I=1,ndeg
      call system_clock(COUNT=count1,COUNT_RATE=count_rate)
 ! calculate derivative vector with respect to the current variable
      do ieq = 1,nex + neqs
        pt=eqmap(ieq,1)
        s1=eqmap(ieq,2)
        s2=eqmap(ieq,3)
        g =eqmap(ieq,4)
        if(g==0)then !is an energy
         if(s2>0)then !its adiabatic energy or adiabatic off diagonal terms(0)
           if(s1==s2)then
             dQ(ieq) = WIJ(pt,s1,s1)
           else 
             dQ(ieq) = WIJ(pt,s1,s2)*(1-DijScale)
           end if
         else  !(s2>0)    it is adiabatic energy difference
           dQ(ieq) = WIJ(pt,s1,s1) - WIJ(pt,-s2,-s2)
         end if !(s2>0)
       else !(g==0)   is gradient
         dQ(ieq) = UIJ(pt,g,s1,s2)
       end if!(g==0)
      end do !ieq
      call system_clock(COUNT=count2,COUNT_RATE=count_rate)
      t4 = t4 + dble(count2-count1)/count_rate
      
!  construct derivative of the Lagrangian with respect to coefficient
      dCi(ico) = - dot_product(PmQ(nex+1:)*weight,dQ(nex+1:)*weight)+dot_product(cilm(nc+1:),dQ(1:nex))+flattening*cilm(ico)
      jaco(1:nex,ico) =  dQ(1:nex)
      call system_clock(COUNT=count1,COUNT_RATE=count_rate)
      t5 = t5 + dble(count1-count2)/count_rate
     end do! iBasis
    end do ! iBlk 

! construct derivative with respect to Lagrange multipliers
    dLambda = - PmQ(1:nex)
    lag =    sum(weight*weight*PmQ(nex+1:)*PmQ(nex+1:))/2-sum(cilm(nc+1:)*PmQ(1:nex))+&
             sum(cilm(1:nc)*cilm(1:nc))*flattening/2
    deallocate(WIJ)
    deallocate(VIJ)
    deallocate(DIJ)
    deallocate(UIJ)

    call system_clock(COUNT=count2)
    if(printlvl>0)then
      print *,"   Gradient/Hessian calculation finished."
      print *,"   Execution time decomposition for getCGrad"
      print *,"   Total   PmQ   c1*c2    W,V,D    dQ   PmQ.dQ  "
      print "(6F8.2)",dble(count2-count0)/count_rate,t1,t2,t3,t4,t5
    end if
  END SUBROUTINE getCGrad
! 

!------------------------------------------------------------------
! Try to obtain the best values of Lagrange multipliers
! Since only Lagrange multipliers are changed, there is no need to 
! update Hd or reconstruct error vector.  
! The optimal values is a saddle point rather than minimum of the 
! Lagrangian, (or in the case of unsatisfied constraints, minimum
! of the norm of the gradient of Lagrangian, the subroutine attempts
! to minimize norm of gradient of L.
! On output, the rows of B are orthongalized by the transformation matrix Q
  SUBROUTINE optLag(nc,nex,B,ldb,dCi,asol,Borth)
    USE progdata, ONlY: printlvl
    IMPLICIT NONE
    INTEGER, INTENT(IN)                                         :: nex, nc, ldb
    DOUBLE PRECISION,dimension(ldb,nc),intent(IN)               :: B
    DOUBLE PRECISION,dimension(ldb,nc),intent(OUT)              :: Borth
    DOUBLE PRECISION, INTENT(INOUT),dimension(nc)               :: dCi
    DOUBLE PRECISION,dimension(nex+nc),INTENT(INOUT)            :: asol

    double precision,dimension(nex,nex)   ::  BBt
    double precision,dimension(nex)       ::  Bg, tmp, dLag,w
    double precision   ::  WORK(2+7*nex+3*nex**2)
    integer :: LWORK,LIWORK,i,INFO,IWORK(4+6*nex),c1,c2,crate
    LWORK  = 2+7*nex+3*nex**2
    LIWORK = 4+6*nex
    ! Bg = B.grdv
    if(nex<=0)return
    if(printlvl>1)then
      print *,"  Calculating optimal values of Lagrange multipliers..."
      call system_clock(COUNT=c1,COUNT_RATE=crate)
    end if
    CALL DGEMV('N',nex,nc,dble(-1),B,ldb,dCi,int(1),dble(0),Bg,int(1))

    ! BBt = B.B^T
    CALL DSYRK('U','N',nex,nc,dble(1),B,ldb,dble(0),BBt,nex)
    CALL DSYEVD('V','U',nex,BBt,nex,W,WORK,LWORK,IWORK,LIWORK,INFO)
    CALL DGEMV('T',nex,nex,dble(1),BBt,nex,Bg,int(1),dble(0),tmp,int(1))
    do i=1,nex
      if(abs(w(i))>exactTol)then
        tmp(i)=tmp(i)/w(i)
      else
        tmp(i)=dble(0)
      end if
    end do!i=1,nex
    CALL DGEMV('N',nex,nex,dble(1),BBt,nex,tmp,int(1),dble(0),dLag,int(1))
    asol(nc+1:)=asol(nc+1:)+dLag
    do i=1,nex
      if(abs(asol(i+nc))>1D0.and.printlvl>1)then
       print "(4x,A,I4,A,F8.2,A,4I5)",&
                 "Large Lagrange multiplier. eq",i,",val=",asol(i+nc),",eqmap=",eqmap(i,:)
      end if
    end do!i=1,nex
    ! update new gradient
    CALL DGEMV('T',nex,nc,dble(1),B,ldb,dLag,int(1),dble(1),dCi,int(1))
    if(printlvl>1)then
      call system_clock(COUNT=c2)
      print "(A,F6.2,A)","  Lagrange multipliers and Lagrangian gradients updated in ",(c2-c1)/dble(crate),"s"
    end if
    ! B'=Q^T.B
    CALL DGEMM('T','N',nex,nc,nex,dble(1),BBt,nex,B,ldb,dble(0),Borth,ldb)
    do i=1,nex
      if(abs(w(i))>exactTol) Borth(i,:)=Borth(i,:)/sqrt(w(i))
    end do!i=1,nex
  END SUBROUTINE optLag
!------------------------------------------------------------------------------------------------------
! subroutine for testing gradients for Hd and/or gradients for Lagrangian
! An eight point 8th order finite difference scheme is used to evaluate the gradients
! 
  SUBROUTINE testSurfgen(disp,hvec,wtvec)
    USE progdata, only : natoms, printlvl
    USE hddata, only : nstates, updateHd,linearizeHd
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN)                  :: disp
    DOUBLE PRECISION, INTENT(IN),dimension(ncons) :: hvec
    DOUBLE PRECISION, DIMENSION(ncons),INTENT(IN) :: wtvec
    integer           ::  i,  n, j, k, prtl
    double precision  ::  Ea(nstates), dHa(3*natoms,nstates,nstates)   ,&
                          hmat(nstates,nstates), dcgrads(3*natoms,nstates,nstates)
    double precision,dimension(3*natoms) ::  cgeom,dgeom,cgeom0
    double precision                     :: hnorm, ind, maxratio
    double precision, dimension(ncons+nex)   ::  asol0, asol, dsol
    double precision,external            :: dnrm2
    double precision  :: DNumHd(nstates,nstates), DAnaHd(nstates,nstates),&   ! Numerical and analytical gradients
                         DNumHa(nstates),DAnaHa(nstates), &                               ! for Hd in both representations
                         DNumL, DAnaL     ! for Hd and Lagrangian
    double precision  :: lag, dCi(ncons), dL(nex)
    double precision, dimension(nex,ncons)                      :: jaco
    integer, parameter:: ndisp=4
    double precision  :: coef(ndisp)=(/ 4/5D0 , -1/5D0, 4/105D0, -1/280D0 /)
 
    if(ntest<=0)return
    prtl = printlvl
    printlvl = 0
    hnorm = dnrm2(ncons,hvec,int(1))
    print *,"Testing gradients of Hd along 3 random pathes"
    print *,"  Cartesian displacement size :", disp
    call init_random_seed()
    maxratio = 0d0
    do i=1,ntest
      print *,"  displacement",i," :  evaluating analytical gradients"
    ! initialize random Hd, geometry and displacements
    ! Hd is initialized as a random offset from existing Hd
      call random_number(asol0)
      asol0 = asol0-5d-1
      asol0 = asol0/dnrm2(ncons,asol0,int(1))*hnorm*5D-2+hvec
      call updateHd(asol0,coefmap,ncons)
      call linearizeHd()
    ! geometry is taken as a random offset from a random data point
      call random_number(cgeom0)
      cgeom0=(cgeom0-5d-1)*1d-1
      call random_number(ind)
      cgeom0=cgeom0+dispgeoms(int(ind*npoints+1))%cgeom
    ! generate a random direction for geometry displacement
      call random_number(dgeom)
      dgeom = (dgeom-5d-1)
      dgeom = dgeom/dnrm2(3*natoms,dgeom,int(1))*disp

    ! generate a random displacement in coefficient space
      call random_number(dsol)
      dsol = dsol - 5D-1
      dsol = dsol/dnrm2(ncons,dsol,int(1))*disp

    ! analytical evaluation of Hd gradients
      call getCartHd(cgeom0,Ea,dHa,hmat,dcgrads)
      DAnaHd = 0d0
      DAnaHa = 0d0
      do j=1,nstates
        DAnaHa(j) = dot_product(dHa(:,j,j),dgeom)
        DAnaHd(j,j) = dot_product(dcgrads(:,j,j),dgeom)
        do k=j+1,nstates
          DAnaHd(j,k) = dot_product(dcgrads(:,j,k),dgeom)
          DAnaHd(k,j) = DAnaHd(j,k)
        end do
      end do
      DAnaHd = DAnaHd/disp
      DAnaHa = DAnaHa/disp

    ! analytical evaluation of Lag gradients
      call updateEigenVec(asol0)
      CALL getCGrad(nex,neqs,ncons,asol0,dCi,dL,lag,wtvec,jaco)
      DAnaL = dot_product(dsol(1:ncons),dCi)
      DAnaL = DAnaL/disp

      print *,""
      print *,"Generating Hd Numerical Gradients (Ana-Num=Diff)"
    ! numerical evaluation of Hd gradients
      DNumHa = 0d0
      DNumHd = 0d0
      do n = 1,ndisp
        cgeom = cgeom0 + dgeom*n
        call getCartHd(cgeom,Ea,dHa,hmat,dcgrads)
        DNumHa = DNumHa+coef(n)*Ea
        DNumHd = DNumHd+coef(n)*hmat
        cgeom = cgeom0 - dgeom*n
        call getCartHd(cgeom,Ea,dHa,hmat,dcgrads)
        DNumHa = DNumHa-coef(n)*Ea
        DNumHd = DNumHd-coef(n)*hmat
      end do
      DNumHa = DNumHa/disp
      DNumHd = DNumHd/disp
      !print out comparisons  
      print *,"Diabatic Representation  : "
      do j=1,nstates
        do k=j,nstates
          print "(2(A,I4),3(A,E14.7))","  Block(",j,",",k,"):  ",DAnaHd(j,k)," - ",DNumHd(j,k)," = ",DAnaHd(j,k)-DNumHd(j,k)
          maxratio = max(maxratio,abs((DAnaHd(j,k)-DNumHd(j,k))/DNumHd(j,k)))
        end do
      end do
      print *,"Adiabatic Representation  : "
      do j=1,nstates
        print "(A,I4,3(A,E14.7))","  State ",j,":   ",DAnaHa(j)," - ",DNumHa(j)," = ",DAnaHa(j)-DNumHa(j)
        maxratio = max(maxratio,abs((DAnaHa(j)-DNumHa(j))/DNumHa(j)))
      end do

    ! numerical evaluation of Lag gradients
      print *,""
      Print *,"Generating Lagrangian Gradients (Ana-Num=Diff)"
      DNumL = 0d0
      do n = 1,ndisp
        asol = asol0 + dsol*n
        call updateEigenVec(asol)
        CALL getCGrad(nex,neqs,ncons,asol,dCi,dL,lag,wtvec,jaco)
        DNumL = DNumL + coef(n)*lag
        asol = asol0 - dsol*n
        call updateEigenVec(asol)
        CALL getCGrad(nex,neqs,ncons,asol,dCi,dL,lag,wtvec,jaco)
        DNumL = DNumL - coef(n)*lag
      end do
      DNumL = DNumL/disp 
   
    ! Comparisons
      print "(3(A,E14.7))","   ",DAnaL," - ",DNumL," = ",DNumL-DAnaL
      maxratio = max(maxratio,abs((DNumL-DAnaL)/DNumL))
    end do!i=1,3
    call updateHd(hvec,coefmap,ncons)
    print "(A,F10.5,A)"," Maximum percentage error : ",maxratio*100,"%"
    print *,""
    printlvl = prtl
  END SUBROUTINE

! Construction of gradients of cost function and project out the component that overlaps
! with exact equation gradients
  SUBROUTINE GradProj(AMat,bvec,dcex,pgrad) 
    USE progdata, only: printlvl
    IMPLICIT NONE
    DOUBLE PRECISION,INTENT(IN)           ::  AMat(neqs+nex,ncons),bvec(neqs+nex)
    DOUBLE PRECISION,INTENT(OUT)          ::  dcex(ncons),pgrad(ncons)
    
    double precision,dimension(nex,nex)   ::  BBt
    double precision,dimension(nex,ncons) ::  BOrth  
    double precision,dimension(nex)       ::  bex,w ,bexw
    double precision                      ::  WORK(2+7*nex+3*nex**2), nrmG,nrmEx,nrmPG,dp(ncons)
    integer                               ::  LWORK,LIWORK,i,INFO,IWORK(4+6*nex), NNull
    double precision,external             ::  dnrm2

    LWORK  = 2+7*nex+3*nex**2
    LIWORK = 4+6*nex
    if(printlvl>1) print *,"Performing Gradient Projections. "

    ! BBt = B.B^T
    CALL DSYRK('U','N',nex,ncons,dble(1),AMat,nex+neqs,dble(0),BBt,nex)
    CALL DSYEVD('V','U',nex,BBt,nex,W,WORK,LWORK,IWORK,LIWORK,INFO)

    ! transform rhs of exact equations 
    CALL DGEMV('T',nex,nex,dble(1),BBt,nex,bvec,int(1),dble(0),bex,int(1))

    ! transform exact equation gradients
    CALL DGEMM('T','N',nex,ncons,nex,dble(1),BBt,nex,AMat,nex+neqs,dble(0),BOrth,nex)
    
    NNull = 0
    do i=1,nex
      if(abs(w(i))>exactTol)then
        w(i)   = 1/w(i)
      else
        w(i)   = dble(0)
        NNull = i
      end if
    end do
    if(printlvl>2) print *,"  Dimensionaliy of Null Space: ",NNull
    !make exact displacement vector
    bexw=bex*w
    CALL DGEMV('T',nex,ncons,dble(1),BOrth,nex,bexw,int(1),dble(0),dcex,int(1))

    !calculation the gradient of cost function
    CALL DGEMV('T',neqs,ncons,dble(1),AMat(nex+1,1),neqs+nex,bvec(nex+1),int(1),dble(0),pgrad,int(1))
    nrmG = dnrm2(ncons,pgrad,int(1))
    dp   = pgrad
    !orthorgonalize it with respect to exact displacements
    do I=NNull+1,nex
      pgrad = pgrad - BOrth(i,:)*dot_product(pgrad,BOrth(i,:))*w(i)
    end do !i=NNull+1,nex
    dp   = dp-pgrad
    nrmPG = dnrm2(ncons,pgrad,int(1))
    nrmEx = dnrm2(ncons,dp,int(1))
    if(printlvl>1)print "(3(A,E12.5))","  Norm of gradients:",nrmG,", Overlapping component:",nrmEx,", Orthogonal component:",nrmPG
  END SUBROUTINE GradProj

!----------------------------------------------------------------------------------------
! updates Hd then evaluate the RMS error from LSE and exact equations
  SUBROUTINE evaluateError(asol,weight,LSErr,ExErr)
    USE progdata, only: printlvl
    USE hddata, only: updateHd
    IMPLICIT NONE
    DOUBLE PRECISION,dimension(ncons),INTENT(IN)  ::  asol
    DOUBLE PRECISION,dimension(neqs),INTENT(IN)   ::  weight
    DOUBLE PRECISION,INTENT(OUT)                  ::  LSErr, ExErr
    double precision,dimension(nex+neqs)                       :: bvec
    integer  ::  plvl
    
    plvl =  printlvl
    printlvl = 0
    CALL updateEigenVec(asol,followPrev)
    CALL makebvec(bvec,.true.)
    bvec(nex+1:)=bvec(nex+1:)*weight
    bvec = bvec*bvec
    if(plvl>1) &
    print "(2(A,E13.6,A,4I4))","Max Error:  Ex=",sqrt(maxval(bvec(1:nex))), "@",eqmap(maxloc(bvec(1:nex)),:), &
                                        ",  LS=",sqrt(maxval(bvec(nex+1:))),"@",eqmap(maxloc(bvec(nex+1:))+nex,:)
    LSErr = sqrt(sum(bvec(nex+1:))/sum(weight*weight))
!    print "(A,I6,A,4I4,A,E12.5)","eq#",maxloc(bvec(1+nex:))+nex,":",eqmap(nex+maxloc(bvec(1+nex:)),:),", err=",sqrt(maxval(bvec(1+nex:)))
!PRINT *,"WEIGHT:",weight(maxloc(bvec(1+nex:))+nex)
    if(nex==0)then
      ExErr = 0d0
    else
      ExErr = sqrt(sum(bvec(:nex))/nex)
    end if
    printlvl = plvl
  END SUBROUTINE evaluateError
END MODULE
!---------------------------------------------
! store ckl to a file
SUBROUTINE writeCkl(cklfl)
  use makesurfdata, only: npoints,ckl
  IMPLICIT NONE
  CHARACTER(72),INTENT(IN) :: cklfl

  integer  ios, i
  open(unit=7421,file=trim(adjustl(cklfl)),access='sequential',form='formatted',&
     status='replace',action='write',position='rewind',iostat=ios)
  if(ios/=0)then
    print *,"writeCkl: failed to open ckl file ",trim(adjustl(cklfl))
    return
  end if
  
  do i=1,npoints
    write(7421,"(4E32.24)")ckl(i,:,:)
  end do!i=1,npoints
  print *," wave functions exported to file ",trim(adjustl(cklfl))
  close(unit=7421)
END SUBROUTINE writeCkl
!---------------------------------------------
SUBROUTINE readCkl(cklfl)
  use makesurfdata, only: npoints, ckl
  IMPLICIT NONE
  CHARACTER(72),INTENT(IN) :: cklfl

  integer  ios, i
  open(unit=7421,file=trim(adjustl(cklfl)),access='sequential',form='formatted',&
     status='old',action='read',position='rewind',iostat=ios)
  if(ios/=0)then
    print *,"readCkl: failed to open ckl file ",trim(adjustl(cklfl))
    return
  end if

  do i=1,npoints
    read(7421,*)ckl(i,:,:)
  end do
  print *," wave functions loaded from file ",trim(adjustl(cklfl))
  close(unit=7421)
END SUBROUTINE readCkl
!---------------------------------------------
! Generate the basis for the fit using linear combinations of primitive basis functions
!   * Values and gradients of primitive functions are evaluated at each data point and put into inter-
!     mediate structures basisVal and basisGrad.  They will later be released.
!   * Basis are generated separately for each block.
!   * Linear combinations that define fitting basis is chosen to be the eigenvectors of the dot-product
!     matrix of primitive functions in the space spanned by all values and gradient data.   
!   * In a one state fitting procedure with equal weights assigned, this generates a orthonormal basis 
!     that allows determination of the fit directly from gradients.   In a multistate scenario with non-
!     equal weights, such basis is much more well conditioned than the primitive functions themselves
!   * Eigenvalues that correspond to low occupation numbers are removed.  These components do not contri-
!     bute to any data points at any diabatrization schemes.  There will also be null space that arise 
!     from each individual d values (wave functions), which will be removed from fitting procedure.
!   
SUBROUTINE genBasis(gradNorm)
  use hddata, only: T3DDList,nl,nr,nblks,nBasis,EvalRawTerms,EvaluateBasis2,deallocDVal,EvaluateVal
  use progdata, only: printlvl
  use makesurfdata, only:  npoints,dispgeoms,nvibs,npb,nbas,dWVals,WVals,ncons,coefMap,ptWeights,energyT,highEScale,&
                           TBas,ZBas,ZBasI,WMat,w_grad,w_energy,w_fij,incener,incgrad,e_exact,g_exact
  IMPLICIT NONE
  DOUBLE PRECISION,dimension(npoints,nvibs,nblks),intent(OUT):: gradNorm
  integer  :: i,j,k,count1,count2,count_rate,pv1,n1,n2
  integer  :: ll,rr,p,nb
  double precision,allocatable,dimension(:)   :: WORK, eval
  double precision,allocatable,dimension(:,:) :: evec, pbas, dmat, pbasw 
  double precision   ::  gn0, wt
  integer,allocatable,dimension(:)  :: IWORK,ISUP
  integer  :: INFO, ng0

  double precision,external :: dnrm2

  call system_clock(COUNT_RATE=count_rate)
  gradNorm=dble(0)

  ! allocate global info arrays

  if(allocated(npb))deallocate(npb)
  allocate(npb(nblks))
  if(allocated(nbas))deallocate(nbas)
  allocate(nbas(nblks))
  if(allocated(ZBas))then
    do i=1,ubound(ZBas,1)
      if(associated(ZBas(i)%List))deallocate(ZBas(i)%List)
    end do
    deallocate(ZBas)
  end if!allocated(ZBas
  allocate(ZBas(nblks))
  if(allocated(ZBasI))then
    do i=1,ubound(ZBasI,1)
      if(associated(ZBasI(i)%List))deallocate(ZBasI(i)%List)
    end do
    deallocate(ZBasI)
  end if!allocated(ZBasI)
  allocate(ZBasI(nblks))
  if(allocated(WMat))then
    do i=1,ubound(WMat,1)
      if(associated(WMat(i)%List))deallocate(WMat(i)%List)
    end do
    deallocate(WMat)
  end if
  allocate(WMat(nblks))

  ! counting number of functions per block and make map for each block
  do k=1,nblks
    npb(k)= 0 
    do i=1,ncons
      if(coefMap(i,2)==k)npb(k)=npb(k)+1
    end do
  end do

! allocate memory for WVals and dWVals global storage for basis values and gradients in local basis
  if(allocated(WVals))then
    do k=1,ubound(WVals,1)
      if(allocated(WVals(k)%List))deallocate(WVals(k)%List)
    end do
    deallocate(WVals)
  end if
  allocate(WVals(1:nblks))
  if(allocated(dWVals))then
    do k=1,ubound(dWVals,1)
      if(allocated(dWVals(k)%List))deallocate(dWVals(k)%List)
    end do
    deallocate(dWVals)
  end if
  allocate(dWVals(1:nblks))

  if(printlvl>0) print "(/,A,/,A,E15.7)"," Generating preconditioned basis free of d-independent linear dependencies.",&
                     "     Eigenvalue cutoff=",TBas
! construct new basis and store them for each block 
  do k=1,nblks
    if(printlvl>0) print "(/,A,I2,A,I5,A)"," Constructing intermediate basis for block ",K," with ",npb(k)," matrices"
    if(npb(k)==0)then
        print *,"WARNING: No basis function is present for block ",K,".  Skipping basis construction."
        nBas(k)=0
        allocate(ZBas(k)%List(1,1))     ! just a place holder
        allocate(ZBasI(k)%List(1,1))
        allocate(WMat(k)%List(pv1,nb))
        gradNorm(:,:,k) = 0d0
        cycle
    end if
    call system_clock(COUNT=count2)
    if(printlvl>1)write(*,"(A)",advance='no'),"   Evaluating primitive basis... "
    ll = nl(k)
    rr = nr(k)
! pbas contains the values and gradients of all primitive basis matrices.
!   row number is the index for surface quantities, looping through blocks>point>energy,each gradient component
! pbasw are equations weighed by 
    allocate(pbas(npoints*(nvibs+1)*ll*rr,npb(k)))
    allocate(pbasw(npoints*(nvibs+1)*ll*rr,npb(k)))

! evaluate and tabulate value and gradients of basis at all data points
    pbas = 0D0
    n1 = 1
    do i=1,npoints
     call EvalRawTerms(dispgeoms(i)%igeom)
! In EvaluateVal, the pbas matrix is reshaped so that it is cognized as a [npoints*(nvibs+1),ll*rr*npb(k)] dimensional
! matrix.   This is equivalent because loops through ll and rr are treated as inner-most in the subroutine
     CALL EvaluateVal(pbas(n1,1),npoints*(nvibs+1),k,nvibs,dispgeoms(i)%bmat)
! create weighed data matrix
     j=0
     do while(dispgeoms(i)%energy(1)>energyT(j+1))
       j=j+1
       if(j==10)exit
     end do
     wt = ptWeights(i)
     if(j>0) wt = wt*highEScale(j)
     if(all(.not.incener(i,:,:)).and.all(.not.e_exact(i,:,:)).or.abs(w_energy*wt)<1d-3)then
       pbasw(n1,:) = 0d0
     else
       pbasw(n1,:) = pbas(n1,:)*wt!(w_energy*wt)
     end if
     if(all(.not.incgrad(i,:,:)).and.all(.not.g_exact(i,:,:)).or.abs(max(w_grad,w_fij)*wt)<1d-3)then
       pbasw(n1+1:n1+nvibs,:) = 0d0
     else
       pbasw(n1+1:n1+dispgeoms(i)%nvibs,:) = pbas(n1+1:n1+dispgeoms(i)%nvibs,:)*wt!(max(w_grad,w_fij)*wt)
       pbasw(n1+dispgeoms(i)%nvibs+1:n1+nvibs,:) =0d0
     end if
     n1 = n1 + (nvibs+1)*ll*rr
    end do!i=1,npoints

    call system_clock(COUNT=count1)
    if(printlvl>1)print 1001,"finished in ",dble(count1-count2)/count_rate," s"

    if(TBas>0)then  ! construct primitive basis
    ! construct overlap between primitive basis
    ! dmat is the dot product (overlap) matrix of primitive basis matrices on the space spanned by the aforementioned
    !   surface quantities.   dmat=pbasw^T.pbasw
        allocate(dmat(npb(k),npb(k)))
        if(printlvl>1)write(*,"(A)",advance='no'),"   Constructing and diagonalizing intermediate basis overlap matrix..."
        CALL DSYRK('U','T',npb(k),npoints*(nvibs+1)*ll*rr,dble(1),pbasw,&
             npoints*(nvibs+1)*ll*rr,dble(0),dmat,npb(k))
        deallocate(pbasw)

        ! diagonalize overlap matrix to obtain orthogonal basis
        allocate(WORK(1))
        allocate(IWORK(1))
        allocate(ISUP(2*npb(k)))
        allocate(evec(npb(k),npb(k)))
        allocate(eval(npb(k)))
        CALL DSYEVR('V','A','U',npb(k),dmat,npb(k),0d0,0d0,0,0,TBas/100,n1,eval,&
              evec,npb(k),ISUP,WORK,int(-1),IWORK,int(-1),INFO)
        IF(INFO/=0)stop "genBasis:  DSYEVR workspace query failed."
        n1 = int(WORK(1))
        n2 = int(IWORK(1))
        deallocate(WORK)
        deallocate(IWORK)
        allocate(WORK(n1))
        allocate(IWORK(n2))
        CALL DSYEVR('V','A','U',npb(k),dmat,npb(k),0.,0.,0,0,0d0,p,eval,&
              evec,npb(k),ISUP,WORK,n1,IWORK,n2,INFO)
        if(INFO/=0)stop "genBasis:  DSYEVR FAILED."
        
        deallocate(WORK)
        deallocate(IWORK)
        deallocate(ISUP)
        deallocate(dmat)
        call system_clock(COUNT=count2)
        if(printlvl>1)print 1001, " finished in",dble(count2-count1)/count_rate," s"
    ! construct transformations to the new basis 
    ! count number of eigenvalues larger than threshold and invert the square roots of the valid ones
    ! ZBas is forward transformation and ZBasI is backwards transformation
        nb= 0 
        do i=npb(k),1,-1
          if(eval(i)<TBas)exit
          nb=nb+1
          eval(i)=sqrt(eval(i))  
        end do
        nBas(k) = nb
        if(printlvl>0)print "(2(A,I6))"," Size of intermediate basis for block ",k," : ",nb
        if(printlvl>0)print "(/,A)"," Constructing fitting basis from ab initio data." 
    ! form the transformation matrix.Z=U.L^-1/2
        allocate(ZBas(k)%List(npb(k),nb))
        allocate(ZBasI(k)%List(nb,npb(k)))
        do i=npb(k),npb(k)-nb+1,-1
          ZBas(k)%List(:,npb(k)-i+1) =evec(:,i)/eval(i)
          ZBasI(k)%List(npb(k)-i+1,:)=evec(:,i)*eval(i)
        end do
        deallocate(evec)
        deallocate(eval)

    ! get the values and gradients in the transformed basis
        if(printlvl>1)write(*,"(A)",advance='no'),"   Generating values and gradients in transformed basis..." 
        pv1 = npoints*(1+nvibs)*ll*rr
    ! generate new coefficient matrix in the basis.  WMat = pbas.ZBas
        allocate(WMat(k)%List(pv1,nb))
        call DGEMM('N','N',pv1,nb,npb(k),1d0,pbas,pv1,ZBas(k)%List,npb(k),0d0,WMat(k)%List,pv1)  
        deallocate(pbas)
        call system_clock(COUNT=count1)
        if(printlvl>1)print 1001, " finished in",dble(count1-count2)/count_rate," s"
!-----------------------------------------------------------------
    else  ! TBas<=0, use the raw basis itself
        if(printlvl>1)print *,"  TBas<=0, skipping null space removal. Using intermediate basis in fit."
        nBas(k) = npb(k)
        allocate(ZBas(k)%List(npb(k),npb(k)))
        allocate(ZBasI(k)%List(npb(k),npb(k)))
        Zbas(k)%List=0d0
        ZBasI(k)%List=0D0
        do i=1,npb(k)
            ZBas(k)%List(i,i) =  1D0
            ZBasI(k)%List(i,i)=  1D0
        end do
        pv1 = npoints*(1+nvibs)*ll*rr
        ! generate new coefficient matrix in the basis.  WMat = pbas
        allocate(WMat(k)%List(pv1,npb(k)))
        WMat(k)%List=pbas
        deallocate(pbasw)
        deallocate(pbas)
    end if!(TBas>0)

! verify if the resulting basis yield a orthogonal coefficient matrix
!    if(printlvl>2)then
!      count2 = count
!      print *,"    Verifying the overlap matrix in new basis."
!      allocate(dmat(nb,nb))
!      CALL DSYRK('U','T',nb,pv1,dble(1),WMat(k)%List,pv1,dble(0),dmat,nb)
!      dmin = dmat(1,1)
!      dmax = dmat(1,1)
!      omin = minval(dmat(1,2:))
!      i    = 1
!      j    = minloc(dmat(1,2:),1)
!      omax = maxval(dmat(1,2:))
!      p    = 1
!      q    = maxloc(dmat(1,2:),1)
!      do l=2,nb
!        DMIN = MIN(DMIN,DMAT(L,L))
!        dmax = max(dmax,dmat(l,l))
!        pval=minval(dmat(l,l+1:))
!        if(pval<omin)then
!          omin =  pval
!          i = l
!          j = minloc(dmat(l,l+1:),1)+l
!        end if
!        pval=maxval(dmat(l,l+1:))
!        if(pval>omax)then
!          omax =  pval
!          p = l
!          q = maxloc(dmat(l,l+1:),1)+l
!        end if
!      end do
!      print "(2(A,F9.4))","       Diagonal elements     : min=",dmin," max=",dmax
!      print "(2(A,F9.4,A,I5,A,I5,A))","       Off-diagonal elements : min=",omin,"@(",I,',',J,"),",&
!                                                              " max=",omax,"@(",P,',',Q,")"
!      deallocate(dmat)
!      call system_clock(COUNT=count1)
!      if(printlvl>1)print 1001, "     verification finished in",dble(count1-count2)/count_rate," s"
!    end if

! generate gradNorm and zero out components that are excluded
    gn0 = 0d0
    ng0 = 0
    if(nBas(k)>0)then
        do i=1,npoints
          do j=1,nvibs
            gradNorm(i,j,k) = dnrm2(ll*rr*nBas(k),WMat(k)%List((i-1)*(nvibs+1)+j+1,1),npoints*(nvibs+1))
          end do
          if(dispgeoms(i)%nvibs<nvibs)then
            gn0= gn0 + sum(gradNorm(i,dispgeoms(i)%nvibs+1:nvibs,k)**2)
            gradNorm(i,dispgeoms(i)%nvibs+1:nvibs,k) = 0d0
            ng0 = ng0 + nvibs-dispgeoms(i)%nvibs
          end if
        end do
    else  !nBas(k)>0
        gradNorm(:,:,k) = 0d0
    end if!nBas(k)>0
    if(printlvl>1) print "(3x,A,E12.5,A,I6,A)","Total norm of contribution zeroed out:",sqrt(gn0)," over ",ng0," equations"
    call system_clock(COUNT=count2)
    if(printlvl>1)print 1001, "     Total contribution to each equation generated in ",dble(count2-count1)/count_rate," s"
  end do! k=1,nblks
1001 format(a,f7.2,a)
END SUBROUTINE genBasis!
!---------------------------------------------
! convert angle and determinant to 2x2 matrix
SUBROUTINE PutAngle(ckl,theta,sg)
  IMPLICIT NONE
  DOUBLE PRECISION,intent(INOUT) :: ckl   (2,2)
  DOUBLE PRECISION,intent(IN)    :: theta 
  INTEGER         ,intent(IN)    :: sg   

  ckl(1,1) = cos(theta)
  ckl(2,2) = ckl(1,1)
  ckl(2,1) = sin(theta)
  ckl(1,2) = -ckl(2,1)
  if(sg<0) ckl(:,2)=-ckl(:,2)
END SUBROUTINE PutAngle
!---------------------------------------------
! convert an 2x2 orthogonal matrix into angles
SUBROUTINE GetAngles(ckl,npoints,theta,sg)
  IMPLICIT NONE
  INTEGER,intent(IN)          :: npoints
  DOUBLE PRECISION,intent(IN) :: ckl   (npoints,2,2)
  DOUBLE PRECISION,intent(OUT):: theta (npoints)
  INTEGER         ,intent(OUT):: sg    (npoints)
  integer i
  do i=1,npoints
    if(ckl(i,1,1)*ckl(i,2,2)-ckl(i,1,2)*ckl(i,2,1)<0)then
 ! eigenvector matrix has a determinant of -1
      sg(i) = -1
      if(ckl(i,2,1)>0)then
        theta(i)  = acos(ckl(i,1,1))
      else
        theta(i)  = -acos(ckl(i,1,1))
      end if
    else ! special orthogonal
      sg(i) = 1
      if(ckl(i,2,1)>0)then
        theta(i)  = acos(ckl(i,1,1))
      else
        theta(i)  = -acos(ckl(i,1,1))
      end if
    end if
  end do
END SUBROUTINE GetAngles
!---------------------------------------------
SUBROUTINE testCoord(geom,step)
  use hddata, only : ncoord
  use progdata, only:natoms,coordmap,CoordSet
  IMPLICIT NONE
  DOUBLE PRECISION,dimension(3*natoms),intent(in)::geom
  double precision,intent(in) :: step

  INTEGER :: I,j,k,m,n
  double precision ::   dgeom(3*natoms),cgrad,    &
                        bmat(ncoord,3*natoms),igeom(ncoord), &
                        bm2(ncoord,3*natoms)
  integer,parameter :: maxstep=2
  double precision,parameter :: fdcoef(maxstep)=[2d0/3,-1d0/12]
  print *,"Testing Coordinate Gradients."
  print *,"initial Cartesian geometry: "
  j=1
  do i=1,natoms
    print "(3F12.5)",geom(j:j+2)
    j=j+3
  end do
  call buildWBmat(geom,igeom,bmat)
  do i=1,ncoord
    m=coordmap(i,1) !index of set
    n=coordmap(i,2) !index in set
    print *,""
    print *,"Testing Coordinate ",i," (set",m,",ind",n,")"
    print *,"  Coord Type=",CoordSet(m)%Type,"  , Scaling Mode=",CoordSet(m)%Scaling
    print *,"  Atoms: ",CoordSet(m)%coord(:,n)
    print "(A)","   Cart# Analytical Numerical  Difference "
    do j=1,3*natoms
    ! calculating numerical gradient of coord i with respect to cartesian j
        cgrad = 0d0
        do k=1,maxstep
            dgeom = geom
            dgeom(j)=dgeom(j)+step*k
            call buildWBmat(dgeom,igeom,bm2)
            cgrad = cgrad+fdcoef(k)*igeom(i)
            dgeom = geom
            dgeom(j)=dgeom(j)-step*k
            call buildWBmat(dgeom,igeom,bm2)
            cgrad = cgrad-fdcoef(k)*igeom(i)
        end do!k
        cgrad=cgrad/step
        if(abs(cgrad-bmat(i,j))>1d-8)  &
            print "(3x,I4,2x,3E11.4)",j,bmat(i,j),cgrad,cgrad-bmat(i,j)
    end do!j=1,3*natoms
  end do!i=1,ncoord
END SUBROUTINE testCoord
!---------------------------------------------
! Fit Hd according to Ab initio Data
SUBROUTINE makesurf()
  use hddata, only: nstates,ncoord,nblks,order,getHdvec,nl,nr,writeHd,nBasis,updateHd,&
                    makecoefmap,ExtractHd,getFLUnit,EvaluateHd3
  use progdata, only: printlvl,OUTFILE,AU2CM1,natoms,PI
  use makesurfdata
  use rdleclse
  use DIIS
  IMPLICIT NONE
  INTEGER                                        :: i,j,k,l,count1,count2,count_rate
  INTEGER                                        :: m,maxeqs, plvl
  INTEGER                                        :: iter,uerrfl,miter
  DOUBLE PRECISION,dimension(:),allocatable      :: bvec,asol,asol1,asol2,dsol,dCi,dLambda,dCi2,hvec,bvecP
  DOUBLE PRECISION,dimension(:,:),allocatable    :: AMat,jaco,jaco2
  DOUBLE PRECISION,dimension(npoints,nvibs,nblks):: gradNorm
  DOUBLE PRECISION,dimension(npoints,nvibs*2)    :: gradtable
  DOUBLE PRECISION,dimension(npoints,2*nstates)  :: enertable
  INTEGER,DIMENSION(:,:),ALLOCATABLE             :: tmpEqMap,tmpExactEq
  DOUBLE PRECISION,dimension(:),allocatable      :: tmpW,weight, grdv1,grdv2,dis0
  DOUBLE PRECISION                               :: adif,nrmener,avgener,lag,gmin,nrmG, dmin, dinc, disp
  double precision                               :: LSErr,ExErr,LSE0,LSE1,LSE2,beta, denom
  DOUBLE PRECISION                               :: nrmgrad,avggrad,nrmD,stepl
  CHARACTER(4)                                   :: c1,c2
  CHARACTER(16),dimension(npoints)               :: rlabs
  CHARACTER(16),dimension(nvibs*2)               :: clabs
  double precision,dimension(npoints,nstates,nstates)   :: errGrad,errGradh
  double precision,dimension(ncoord)             :: errgrd
  double precision,dimension(ncoord,3*natoms)    :: binv
  CHARACTER(50)                                  :: fmt
  double precision            ::  dener(nstates,nstates)
  integer                     ::  ios, status
  double precision, external  ::  dnrm2
  logical                     ::  diff  !whether differential convergence will be used
  logical                     ::  loadC !whether initial CKL will be loaded from input file
  DOUBLE PRECISION,dimension(nstates,nstates)            :: hmatPT
  DOUBLE PRECISION,dimension(nvibs,nstates,nstates)      :: dhmatPT
  
  integer  :: ncon_total
  double precision,dimension(0:maxiter,npoints)   ::   theta  ! rotation angle
  integer,dimension(0:maxiter,npoints)            ::   sg     ! determinant of eigenvectors at each point
  double precision,dimension(npoints)   ::   theta2

!   call testCoord(dispgeoms(1)%cgeom,1d-5)  !perform testings for coordinate definitions
!   stop                                     !  only use these for new coordinates

  if(printlvl>0)print *,"Entering makesurf"
  call getPtList()

  if(printlvl>0)print*,"   Making list of unknown coefficients."
  ncons = 0
  do i = 1,nblks
   do j= 0,order
     ncons=ncons+nBasis(j,i)
   end do
  enddo
  allocate(coefMap(ncons,3))
  call makeCoefMap(0,order,coefMap,ncons)

  if(printlvl>0)print*,"        ncons=",ncons

  if(allocated(ckl))deallocate(ckl)
  allocate(ckl(npoints,nstates,nstates))
  if(allocated(fitE))deallocate(fitE)
  if(allocated(fitG))deallocate(fitG)
  allocate(fitE(npoints,nstates))
  allocate(fitG(npoints,nvibs,nstates,nstates))

  call genBasis(gradNorm)

  ncon_total = ncons
  ncons = sum(nBas)
  if(printlvl>0)then
    print *,""
    print *,"  Total number of basis matrices:",ncons
    print *,""
    print *,"Making list of fitting equations."
  end if!printlvl>0
  maxeqs=npoints*nstates*(nstates+1)/2*(nvibs+1)
  allocate(tmpEqMap(maxeqs,4))
  allocate(tmpExactEq(maxeqs,4))
  allocate(tmpW(maxeqs))
  call makeEqMap(maxeqs,tmpEqMap,tmpExactEq,gradNorm,tmpW)
  allocate(EqMap(nex+neqs,4))
  allocate(weight(neqs))
  weight=tmpW(1:neqs)
  eqMap(1:nex,:)=tmpExactEq(1:nex,:)
  eqmap(nex+1:nex+neqs,:)=tmpEqMap(1:neqs,:)
  deallocate(tmpW)
  deallocate(tmpEqMap)
  deallocate(tmpExactEq)
  if(printlvl>0)print*,"    Number of Least Squares Equations:",neqs
  if(printlvl>0)print*,"    Number of Exact Equations:        ",nex
  call printSurfHeader(ncons,neqs,nex)
! allocate spaces for arrays and matrices 
  if((neqs+nex)<ncons.and.printlvl>0)print *,"    Warning: neqs<ncons, data set might be inadequate."
  print 1001,"    Memory required to store coefficient matrix:",(neqs+nex)*ncons*7.62939453125D-6," MB"
  allocate(hvec(ncon_total))               ! Hd coefficients vector in original form
  allocate(bvec(neqs+nex))                 !rhs for exact and LSE equations
  allocate(bvecP(neqs+nex))                !rhs for exact and LSE equations from previous iteration
  allocate(asol(ncons+nex))                ! solution vector
  allocate(asol1(ncons+nex))               ! old solution vector
  allocate(asol2(ncons+nex))               ! old solution vector
  allocate(dsol(ncons))                    ! change of solution vector
  allocate(grdv1(ncons))                   ! gradients of Lagrangian of the last iteration 
  allocate(grdv2(ncons))                   ! gradients of Lagrangian(with respect to coefficients)
  allocate(dis0(ncons))                    ! displacement of the previous iteration 
  allocate(jaco(nex,ncons))                ! second derivatives of Lagrangian for exact-LSE jacobian block
  allocate(jaco2(nex,ncons))
  allocate(dCi(ncons))
  allocate(dCi2(ncons))
  allocate(dLambda(nex))
  allocate(NEL(ncons+nex,ncons+nex),STAT=status)
  if(status/=0) stop 'makesurf  :  failed to allocate memory for NEL'
  allocate(rhs(ncons+nex))

  iter = 0
  asol2 = dble(0)
  call getHdvec(hvec,coefmap,ncon_total)
  ! convert primitive h expansion to block orthogonal basis expansions
  count1 = 0
  do k=1,nblks
    count2 = 0
    do i=1,ncon_total
      if(coefMap(i,2)==k)then
        count2 = count2+1
        asol2(count1+1:count1+nBas(k))=asol2(count1+1:count1+nBas(k))+ZBasI(k)%List(:,count2)*hvec(i)
      end if
    end do
    count1=count1+nBas(k)
  end do!k

  !----------------------------------
  ! Begin self-consistent iterations
  !----------------------------------
!generate initial eigen vectors
  if(ckl_input/='')then
    loadC = .false.
    if(printlvl>0)print 1004
    call readCkl(ckl_input) 
    if(printlvl>1)print *,"    Energies from initial Hd and eigenvectors:"
    do i=1,npoints
      CALL EvaluateHd3(asol2,nBas,npoints,i,nvibs,hmatPt,dhmatPt,WMat)
      call OrthGH_Hd(dispgeoms(i),dhmatPt(:dispgeoms(i)%nvibs,:,:),ckl(i,:,:),100,sqrt(gcutoff)/10)
      do k=1,dispgeoms(i)%nvibs
        fitG(i,k,:,:)=matmul(transpose(ckl(i,:,:)),matmul(dhmatPt(k,:,:),ckl(i,:,:)))
      end do!k=1,dispgeoms(i)%nvibs
      do k=1,nstates
        fitE(i,k)= dble(0)
        do j=1,nstates
          fitE(i,k)=fitE(i,k)+ckl(i,j,k)*dot_product(ckl(i,:,k),hmatPt(j,:))
        end do
      end do
      if(printlvl>1)PRINT "(I6,10(x,F15.4))",I,fitE(i,:)*au2cm1
    end do
  else
    loadC = .true.
    if(printlvl>0)print 1003
    printlvl=printlvl-1
    call updateEigenVec(asol2)
    printlvl=printlvl+1
  endif!(ckl_input/='')then
  call testSurfgen(1D-5,asol2,weight)
  asol = asol2
  CALL getCGrad(nex,neqs,ncons,asol,dCi,dLambda,lag,weight,jaco)
  CALL optLag(ncons,nex,jaco,nex,dCi,asol,jaco2)
  nrmG = sqrt(dot_product(dCi,dCi)+dot_product(dLambda,dLambda))
  grdv1 = 0d0
  grdv2 =-dCi
  print "(3(A,E15.7))","Gradients for coef block: ",dnrm2(ncons,dCi,int(1)),",lag block:",dnrm2(nex,dLambda,int(1)),&
            ", total:",sqrt(dot_product(dCi,dCi)+dot_product(dLambda,dLambda))

  CALL initDIISg(ndiis,ndstart,ncons,asol,dCi)
  CALL getError(nrmener,avgener,nrmgrad,avggrad)
  adif = dnrm2(ncons,asol,int(1))
  write(OUTFILE,1005)iter,adif,nrmgrad*100,avggrad*100,nrmener*AU2CM1,avgener*AU2CM1
  adif=toler+1
  diff = .false.
  ! ----------  MAIN LOOP ---------------
  if(stepmethod==1)then
    diff=.true.
    dfstart=-1
  end if
  call GetAngles(ckl,npoints,theta(0,:),sg(0,:))
  do while(iter<maxiter.and.adif>toler)
   iter = iter + 1
   if(printlvl>0)print 1000,iter
   
   if(iter == dfstart)then
     print *,"  Starting differential convergence..."
     diff = .true.
     write(OUTFILE,*)"  Starting differential convergence..."
   end if


   select case(stepMethod)
    ! method 1:  gradient projection.   
    ! Derivatives of cost function is first evaluated.  
    ! Gradients of exact equations are orthogonalized and are used to project out
    ! the overlapping components of cost function gradients.
    ! In the exact dimensions linear step is taken to yield exact results.  
    ! In the othorgonal space, line search is performed to minimize fitting error
    case(1) 
     !Make coefficients vector for least squares and exact equations
     allocate(AMat(neqs+nex,ncons))    !Coefficients for least squares equations
                                        ! and exact equations. Exact block is on the top
     call makecoefmat(AMat)
     call makebvec(bvec,diff)
     ! get projected gradient and displacement in exact block.
     call evaluateError(asol,weight,LSErr,ExErr)
     if(printlvl>1)PRINT "(A,2E13.6)","  RMS Errors before taking step(LSE,Ex): ",LSErr,ExErr
     miter = 0
     plvl = printlvl
     if(printlvl>1)PRINT *,"   Taking steps for exact equations."
     do while(ExErr>ExConv.and.miter<mmiter)
       miter = miter+1
       dCi2 = dsol
       call GradProj(AMat,bvec,dsol,dCi)  !dsol: exact; dCi: LS
       nrmD = dnrm2(ncons,dsol,int(1))
       if(nrmD>maxED)then
         dsol = dsol/nrmD*maxED
         nrmD=maxED
       end if
       if(miter==1)then
         dCi2=dsol
         printlvl = 0
       end if
       PRINT "(A,I4,3(A,E13.6),A,F8.2)","   miter=",miter,", disp=",nrmD,&
           ", RMSE(LS)=",LSErr,", RMSE(Ex)=",ExErr,&
           ", angle= ",ACOS(dot_product(dsol,dCi2)/sqrt(dot_product(dsol,dsol)*dot_product(dCi2,dCi2)))*5.729578D1
       asol = asol+dsol*scaleEx
       call evaluateError(asol,weight,LSErr,ExErr)
       call makecoefmat(AMat)
       call makebvec(bvec,.true.)
     end do! while(ExErr>ExConv .and. miter<mmiter)
     if(plvl>1)PRINT "(A,2E13.6)","  Error after initial steps(LSE,Ex): ",LSErr,ExErr
     LSE0 = LSErr
     LSE1 = LSErr
     do i=1,neqs
       AMat(nex+i,:)=AMat(nex+i,:)*weight(i)
       bvec(nex+i)  =bvec(nex+i)  *weight(i)
     end do
     printlvl = plvl
     call GradProj(AMat,bvec,dsol,dCi)  !dsol: exact; dCi: LS
     PRINT *,"  Taking step along reduced gradient direction."
     printlvl = 0
     asol1=asol
     ! determine proper bond
     dCi=dCi/dnrm2(ncons,dCi,int(1))
     stepl=toler+1D-30
     asol=asol1+stepl*dCi
     call evaluateError(asol,weight,LSErr,ExErr)
     print "(A,2E13.6)","   initial step error:",LSErr,ExErr
     do while(LSErr<LSE1.and.stepl<maxD)
       stepl=stepl*linSteps
       asol=asol1+stepl*dCi
       LSE1=LSErr
       call evaluateError(asol,weight,LSErr,ExErr)
       print "(A,2E13.6)","   trial step error:  ",LSErr,ExErr
     end do!while(stepl>toler)
     print "(A,E13.6)","   crude step size:  ",stepl
     asol=asol1+stepl*dCi
     call evaluateError(asol,weight,LSErr,ExErr)
     call makecoefmat(AMat)
     call makebvec(bvec,.true.)
     do i=1,neqs
       AMat(nex+i,:)=AMat(nex+i,:)*weight(i)
       bvec(nex+i)  =bvec(nex+i)  *weight(i)
     end do
     call GradProj(AMat,bvec,dsol,dCi2)  !dsol: exact; dCi: LS
     dCi=(dCi+dCi2+dsol)/2
     print *,"   movement direction adjusted."
     if(stepl>maxD)then
       if(LSErr>=LSE1)stepl=stepl/linSteps
       if(stepl>maxD)stepl=maxD
       asol=asol1+stepl*dCi
       print *,"Step too large, resized to maxD=",maxD
     else !if stepl>maxD
       stepl=stepl/linSteps
       if(stepl<toler)then
         asol=asol1+stepl*dCi
         PRINT "(A,E13.6)","Step size smaller than convergence tolerance. "
       else!stepl<toler
         do i=2,linSteps
           asol=asol1+stepl*dCi*i
           call evaluateError(asol,weight,LSE2,ExErr)
           print "(A,2E13.6)","   linear step error:",LSErr,ExErr
           if(LSE2>LSE1)then
             exit
           else
             LSE0=LSE1
             LSE1=LSE2
           end if
         end do!i=2,linSteps
         stepl=stepl/2/(LSE0-2*LSE1+LSE2)* &
          ((2*i-1)*LSE0+(4-4*i)*LSE1+(2*i-3)*LSE2)
         asol = asol1+stepl*dCi
         print "(A,3E13.6)","  Interpolating between points with errors:",LSE0,LSE1,LSE2
         call evaluateError(asol,weight,LSErr,ExErr)
         print "(A,2E13.6)","  Final errors : ",LSErr,ExErr
       end if!stepl<toler
     end if!(stepl>maxD)
     printlvl = plvl
     deallocate(AMat)

    ! displacement determined by multidimensional conjugate gradient approach
    case (2)
      if(printlvl>0)print *," Predicting step using generalized conjugate gradients..."
      if(iter==1)then
        if(printlvl>1)print *,"   Initial step taken as gradients..." 
        dsol = grdv2 
        dis0 = dsol 
      else 
        denom = dot_product(grdv1,grdv1)
        if(denom<1D-30)then
          beta = 0d0
        else
          beta = dot_product(grdv2,grdv2-grdv1)/denom
        end if 
        beta = max(0d0,beta)
        if(printlvl>1)print *,"   Factor beta=",beta 
        dsol = grdv2 + beta*dis0
        dis0 = dsol
      end if
      dsol = dsol*gscaler

    ! Default direction method is solving pseudo-normal equations for the linearized
    !   least squares problem with exact equations achieved by lagrange multipliers    
    case default
     !Make coefficients vector for least squares and exact equations
     allocate(AMat(neqs+nex,ncons))    !Coefficients for least squares equations
                                        ! and exact equations. Exact block is on the top
     call makecoefmat(AMat)
     call makebvec(bvec,diff)
     if(printlvl>2.and.diff)then
       if(iter>1.and.dnrm2(neqs+nex,bvecP,1)<dnrm2(neqs+nex,bvec,1))then
         print *,"  Overall error increasing.  Tabulating change in equations..."
         ! look for equations where error got larger
         do i=1,nex+neqs
          if((abs(bvec(i))-abs(bvecP(i)))/abs(bvecP(i))>1d-1.and.abs((abs(bvec(i))-abs(bvecP(i))))>1d-3)then 
           if(eqmap(i,4)==0)then
             if(eqmap(i,2)==eqmap(i,3))then
               print "(A,I6,I2,A,F10.1,A,E12.5,A,SP,F9.2,A)",&
                    "energy(pt,s):",eqmap(i,1:2),",val=",&
                    (dispgeoms(eqmap(i,1))%energy(eqmap(i,2)))*au2cm1,&
                    ",err: ",bvec(i),"(",(abs(bvec(i))-abs(bvecP(i)))/abs(bvecP(i))*100,"%)"
             else
               if(eqmap(i,3)>0)then !off diagonal block values
                 print "(A,I6,2I2,A,2F10.1,A,E12.5,A,SP,F9.2,A)",&
                    "E-diff(pt,s1,s2):",eqmap(i,1:3),",energies=",&
                    au2cm1*(dispgeoms(eqmap(i,1))%energy(abs(eqmap(i,2:3)))),&
                    ",err: ",bvec(i),"(",(abs(bvec(i))-abs(bvecP(i)))/abs(bvecP(i))*100,"%)"
               else   ! energy differences
                 print "(A,I6,2I2,A,2F10.1,A,E12.5,A,SP,F9.2,A)",&
                    "off-diag(pt,s1,s2):",eqmap(i,1:3),",energies=",&
                    au2cm1*(dispgeoms(eqmap(i,1))%energy(eqmap(i,2:3))),&
                    ",err: ",bvec(i),"(",(abs(bvec(i))-abs(bvecP(i)))/abs(bvecP(i))*100,"%)"
               end if
             end if
           else
             if(eqmap(i,2)==eqmap(i,3))then
               print "(A,I6,2I2,A,F10.1,A,E12.5,A,SP,F9.2,A)",&
                    "grad(pt,s,g):",eqmap(i,[1,2,4]),",state energie=",&
                    au2cm1*(dispgeoms(eqmap(i,1))%energy(eqmap(i,2))),&
                    ",err: ",bvec(i),"(",(abs(bvec(i))-abs(bvecP(i)))/abs(bvecP(i))*100,"%)"
             else
               print "(A,I6,3I2,A,2F10.1,A,E12.5,A,SP,F9.2,A)",&
                    "coupling(pt,s1,s2,g):",eqmap(i,:),",state energies=",&
                    au2cm1*(dispgeoms(eqmap(i,1))%energy(eqmap(i,2:3))),&
                    ",err: ",bvec(i),"(",(abs(bvec(i))-abs(bvecP(i)))/abs(bvecP(i))*100,"%)"
             end if
           end if!eqmap(i,4)==0
          end if! if change is large enough
         end do!i
       end if!iter>1 .and. 
     end if!printlvl>2 .and. diff
     bvecP = bvec
     ! create linear equality constrained normal equations
     do i=1,neqs
       AMat(nex+i,:)=AMat(nex+i,:)*weight(i)
       bvec(nex+i)  =bvec(nex+i)  *weight(i)
     end do
      
     Print *,"   Making linear equality constrained normal equations."
     do i=1,nex
       AMat(i,:) = AMat(i,:)*scaleEx
       bvec(i)   = bvec(i)*scaleEx
     end do

     CALL DSYRK('U','T',ncons,neqs,dble(1),AMat(nex+1,1),&
         neqs+nex,dble(0),NEL(nex+1,nex+1),ncons+nex)
     do i=nex+1,nex+ncons-1
       CALL DCOPY(ncons+nex-i,NEL(i,i+1),ncons+nex,NEL(i+1,i),int(1))
     end do
  
     if(diff)then
       print "(2(A,E12.5))"," SumSq error for exact block:",dnrm2(nex,bvec,1),", LSE block:",dnrm2(neqs,bvec(nex+1),1)
       rhs(1:nex)=-dLambda*scaleEx
       rhs(nex+1:)=-dCi
     else
       !construct rhs for normal eqs block y2=A**T.y
       rhs(1:nex)=bvec(1:nex)
       CALL DGEMV('T',neqs,ncons,dble(1),AMat(nex+1,1),neqs+nex,&
          bvec(nex+1),int(1),dble(0),rhs(nex+1),INT(1))
     end if!diff
  
     !construct Lagrange Multipliers block
       do i=1,nex
       NEL(i,nex+1:nex+ncons)=AMat(i,:)
       NEL(nex+1:nex+ncons,i)=AMat(i,:)
     end do
     !zero out L-L block
     NEL(1:nex,1:nex)=dble(0)
 ! 
     deallocate(AMat)
  
     if(diff)then
       ! solve normal equations for change in coefficients
       call system_clock(COUNT=count1)
       CALL solve(ncons,neqs,nex,NEL,rhs,&
              exacttol,lsetol,dsol,printlvl)
     else
       dsol = asol(1:ncons)
       ! solve normal equations
       call system_clock(COUNT=count1)
       CALL solve(ncons,neqs,nex,NEL,rhs,&
              exacttol,lsetol,asol,printlvl)
       dsol = asol(1:ncons)-dsol
     end if

     CALL cleanArrays()

     call system_clock(COUNT=count2,COUNT_RATE=count_rate)
     if(printlvl>0)print 1002,dble(count2-count1)/count_rate  
   end select !step method

   ! Perform line search
   if(stepMethod.ne.1)then 
     nrmD = dnrm2(ncons,dsol,int(1))
     if(printlvl>0)print "(A,E15.7)","   Original size of displacement : ",nrmD
     ! scaling the whole displacement vector is larger than the maximum size
     if(nrmD>maxd)then
       if(printlvl>0)print "(2(A,E12.5))","  Displacement vector scaled from",nrmD," to ",maxD
       dsol=dsol*maxD/nrmD
       nrmD = maxD
     end if 
     if(linSteps+linNegSteps<=0)then
       asol(1:ncons)=asol2(1:ncons)+dsol
     else !linsteps<=0
       dmin = 0d0
       disp = 0d0
       dinc = 1d0
       gmin = nrmG 
       plvl=printlvl
       printlvl=0
       print *,"  Performing ",linSteps," positive displacements."
       do i=1,linSteps
         disp = disp+dinc
         asol(1:ncons)=asol2(1:ncons)+dsol*disp
         call evaluateError(asol,weight,LSErr,ExErr)
         CALL getError(nrmener,avgener,nrmgrad,avggrad)
         CALL getCGrad(nex,neqs,ncons,asol,dCi,dLambda,lag,weight,jaco)
         CALL optLag(ncons,nex,jaco,nex,dCi,asol,jaco2)
         nrmG = sqrt(dot_product(dCi,dCi)+dot_product(dLambda,dLambda))
         if(plvl>1)print "(A,I3,A,E14.7,A,2E12.5,A,2E12.5,A,E16.9)","    #",i," d=",disp*nrmD,": evalEr=",LSErr,ExErr,&
               ",getEr=",nrmener,nrmgrad,",GLag=",nrmG
         if(nrmG<gmin)then
           gmin = nrmG
           dmin = disp
         else if(i<linSteps)then
           !roll back and shrink steps
           disp = disp - dinc 
           dinc = dinc / (1+linSteps-i)
         end if
       end do
!  Negative displacements
       disp = 0d0
       dinc =-1d0
       print *,"  Performing ",linNegSteps," negative displacements."   
       do i=1,linNegSteps
         disp = disp+dinc
         asol(1:ncons)=asol2(1:ncons)+dsol*disp
         call evaluateError(asol,weight,LSErr,ExErr)
         CALL getError(nrmener,avgener,nrmgrad,avggrad)
         CALL getCGrad(nex,neqs,ncons,asol,dCi,dLambda,lag,weight,jaco)
         CALL optLag(ncons,nex,jaco,nex,dCi,asol,jaco2)
         nrmG = sqrt(dot_product(dCi,dCi)+dot_product(dLambda,dLambda))
         if(plvl>1)print "(A,I3,A,E14.7,A,2E12.5,A,2E12.5,A,E16.9)","    #",i," d=",disp*nrmD,": evalEr=",LSErr,ExErr,&
               ",getEr=",nrmener,nrmgrad,",GLag=",nrmG
         if(nrmG<gmin)then
           gmin = nrmG
           dmin = disp
         else if(i<linSteps)then
           !roll back and shrink steps
           disp = disp - dinc
           dinc = dinc / (1+linSteps-i)
         end if
       end do

       if(plvl>1)print "(2(A,E12.5))","   optimal step length:",dmin,", norm of grad=",gmin
       asol(1:ncons)=asol2(1:ncons)+dsol*dmin
       printlvl=plvl
     end if!linSteps<=0

     ! make new eigenvectors
     call updateEigenVec(asol,followPrev)
     CALL getCGrad(nex,neqs,ncons,asol,dCi,dLambda,lag,weight,jaco)
     CALL optLag(ncons,nex,jaco,nex,dCi,asol,jaco2)
     grdv1 = grdv2
     grdv2 = -dCi
     nrmG = sqrt(dot_product(dCi,dCi)+dot_product(dLambda,dLambda))
     print "(3(A,E15.7))","Gradients for coef block: ",dnrm2(ncons,dCi,int(1)),",lag block:",dnrm2(nex,dLambda,int(1)),&
            ", total:",sqrt(dot_product(dCi,dCi)+dot_product(dLambda,dLambda))
   end if
   call GetAngles(ckl,npoints,theta(iter,:),sg(iter,:))
   theta2 = theta(iter,:) - theta(iter-1,:)
   PRINT *,"Maximum allowed change in diabatrization angle is ",maxRot
   do i=1,npoints
     if(theta2(i)>PI)theta2(i)=theta2(i)-2*PI
     if(theta2(i)<-PI)theta2(i)=theta2(i)+2*PI
     if(abs(theta2(i))>maxRot.or.sg(iter,i).ne.sg(iter-1,i))then
       if(abs(theta2(i))>1d-1) &
            print "(A,I5,A,2F10.2,A,SP,F7.3,SS,A,L1)"," point ",I,",E1,2=", (dispgeoms(i)%energy(1:2))*au2cm1,  &
                   ", rotation =",theta2(i),", change=",sg(iter,i).ne.sg(iter-1,i)
       if(maxRot>1D-30) then
         theta(iter,i) = theta(iter-1,i)+sign(maxRot,theta2(i))
         if(maxRot>0) sg(iter,i)=sg(iter-1,i)
         print *,"BEFORE:",ckl(i,:,:)
         call PutAngle(ckl(i,:,:),theta(iter,i),sg(iter,i))
         print *,"AFTER :",ckl(i,:,:)
       end if
     end if
   end do!i=1,npoints (printing diabatrization angles)
   if(printlvl>1)print *,"  Pushing coefficients into DIIS data set."
   CALL pushDIISg(asol,dCi,asol1)
   if(printlvl>1)print *,"  Size of change through DIIS procedure :" ,dnrm2(ncons,asol1-asol,1)
   asol = asol1
   CALL getError(nrmener,avgener,nrmgrad,avggrad)
   adif=DNRM2(ncons,asol-asol2,int(1))
   asol2=asol
   ! write iteration information to output file
   write(OUTFILE,1005)iter,adif,nrmgrad*100,avggrad*100,nrmener*AU2CM1,avgener*AU2CM1
   print *,"   Norm of coefficients:  ",dnrm2(ncons,asol,int(1))
  enddo !while(iter<maxiter.and.adif>toler)
  !----------------------------------
  ! End of self-consistent iterations
  !----------------------------------
  if(adif.le.toler)then
   write(OUTFILE,1006)iter
  else
   write(OUTFILE,1007)
  endif
  if(printlvl>0) print *,"Points where diabatrization angles are changing significantly"
  do i=1,npoints
    if(maxval(theta(0:iter,I))-minval(theta(0:iter,I))>1d-1) &
         print "(I5,15F7.3)",I,theta(iter-min(iter,10):iter,I)
  end do

  !------------------------------------------------------------
  ! Write coefficients of final surface to file
  !------------------------------------------------------------
  if(outputfl/='')then
  ! convert hd to primitive form 
    count1 = 0
    do k=1,nblks
      count2 = 0
      do i=1,ncon_total
        if(coefMap(i,2)==k)then
          count2 = count2+1
          hvec(i) = dot_product(ZBas(k)%List(count2,:),asol2(count1+1:count1+nBas(k)))
        end if
      end do
      count1=count1+nBas(k)
    end do!k
  !----------------------------------
  ! save h vector to file
    call updateHd(hvec,coefmap,ncon_total)
    if(printlvl>1)print *,"  Exporting Hd coefficients to ",trim(outputfl)
    call writehd(outputfl,flheader,.false.)
  end if
  !------------------------------------------------------------
  ! Write final eigenvectors to file 
  !------------------------------------------------------------
  if(ckl_output/='')then
    if(printlvl>1)print *,"  Exporting eigenvectors to ",trim(ckl_output)
    call writeCkl(ckl_output)
  end if
  !------------------------------------------------------------
  ! Compare ab initio and computed derivative couplings
  !------------------------------------------------------------
  errGrad=dble(0)
  errGradh=dble(0)
  do j = 1,npoints
   write(c2,'(i4)')j
   rlabs(j) = ' GM '//trim(adjustl(c2))
  enddo!j=1,npoints
  if(useIntGrad)then
    do i = 1,nvibs
     write(c1,'(i2)')i
     clabs(2*i-1) = ' F[x'//trim(adjustl(c1))//']AB'
     clabs(2*i)   = ' F[x'//trim(adjustl(c1))//']Err'
    end do
  else
    do i = 1,natoms
     write(c1,'(i2)')i
     clabs(6*i-5) = ' F[x'//trim(adjustl(c1))//']AB'
     clabs(6*i-4) = ' F[x'//trim(adjustl(c1))//']Err'
     clabs(6*i-3) = ' F[y'//trim(adjustl(c1))//']AB'
     clabs(6*i-2) = ' F[y'//trim(adjustl(c1))//']Err'
     clabs(6*i-1) = ' F[z'//trim(adjustl(c1))//']AB'
     clabs(6*i)   = ' F[z'//trim(adjustl(c1))//']Err'
    enddo!i=1,natoms
  end if!useIntGrad
  write(OUTFILE,1014)
  !Write ab initio and fit gradients at all input geommetries for all blocks
  do i = 1,nstates
    write(OUTFILE,1016)i
    do l = 1,npoints
     do m = 1,nvibs
      gradtable(l,2*m-1) = dispgeoms(l)%grads(m,i,i)
      gradtable(l,2*m)   = fitG(l,m,i,i)-dispgeoms(l)%grads(m,i,i)
      if(m<=dispgeoms(l)%nvibs)errGrad(l,i,i)=errGrad(l,i,i)+gradtable(l,2*m)**2!*dispgeoms(l)%scale(m)
     enddo !m=1,dispgeoms(i)%nvibs
     errGrad(l,i,i)=sqrt(errGrad(l,i,i))
    enddo!l=1,npoints
    call printMatrix(OUTFILE,rlabs,clabs,int(8 ),npoints,2*nvibs,gradtable,int(13),int(8))
  enddo!i=1,nstates

  do i = 1,nstates
   do j = 1,i-1
    write(OUTFILE,1015)j,i
    do l = 1,npoints
     do m = 1,nvibs
      gradtable(l,2*m-1) = dispgeoms(l)%grads(m,i,j)
      gradtable(l,2*m)   = fitG(l,m,i,j)-dispgeoms(l)%grads(m,i,j)
      if(m<=dispgeoms(l)%nvibs)then
          errGrad(l,i,j)=errGrad(l,i,j)+((&
                      fitG(l,m,i,j)/abs(fitE(l,j)-fitE(l,i))- &
                      dispgeoms(l)%grads(m,i,j)/abs(dispgeoms(l)%energy(j)-dispgeoms(l)%energy(i))   &
                      )  )**2!*dispgeoms(l)%scale(m)
          errGradh(l,i,j)=errGradh(l,i,j)+&
                 (fitG(l,m,i,j)-dispgeoms(l)%grads(m,i,j))**2!*dispgeoms(l)%scale(m)
          errGradh(l,j,i)=errGradh(l,j,i)+&
                 (fitG(l,m,i,i)-fitG(l,m,j,j)-&
                 dispgeoms(l)%grads(m,i,i)+dispgeoms(l)%grads(m,j,j))**2!*dispgeoms(l)%scale(m)
      end if
     enddo !m=1,dispgeoms(i)%nvibs
     errGrad(l,i,j)=sqrt(errGrad(l,i,j))
     errGrad(l,j,i)=errGrad(l,i,j)
     errGradh(l,i,j)=sqrt(errGradh(l,i,j))
     errGradh(l,j,i)=sqrt(errGradh(l,j,i))/2
    enddo!l=1,npoints
    call printMatrix(OUTFILE,rlabs,clabs,int(8 ),npoints,2*nvibs,gradtable,int(13),int(8))
   enddo!j=1,i-1
  enddo!i=1,nstates
  !------------------------------------------------------------------
  ! Compare ab initio and computed energies
  !------------------------------------------------------------------
  write(OUTFILE,1017)
  do i = 1,nstates
   write(c1,'(i2)')i
   clabs(2*i-1) = '  E['//trim(adjustl(c1))//'](AB)'
   clabs(2*i)   = 'DE['//trim(adjustl(c1))//'](FIT)'
  enddo!i=1,nstates
  do j = 1,npoints
   do k = 1,nstates
    enertable(j,2*k-1) = (dispgeoms(j)%energy(k))*AU2CM1
    enertable(j,2*k)   = fitE(j,k)*AU2CM1 - enertable(j,2*k-1)
   enddo!k=1,nstates
  enddo!j=1,npoints
  call printMatrix(OUTFILE,rlabs,clabs,2*nstates,npoints,2*nstates,enertable,int(19),int(9))
  print *," Weighed Error Norms for Energy Gradients and Couplings"
  print *," Gradients are given in the following order:"
  print "(10(I4,'-',I4,',',4X),I4,'-',I4)",(((/j,k/),k=1,j),j=1,nstates)
  write(fmt,'("(A,2X,A,2X,",I2,"A)")'),(nstates+1)*nstates
  print trim(fmt)," PT ","  WT  ",(("   ErrG  ",k=1,j),j=1,nstates),&
                           ((" Ab Grd  ",k=1,j),j=1,nstates)
  fmt=""
  write(fmt,'("(I5,2X,F6.3,2X,",I2,"(X,E12.5))")'),(nstates+1)*nstates
  do i=1,npoints
    do k=1,nstates
      dener(k,k)=1D0
      do l=k+1,nstates
        dener(k,l)= abs(dispgeoms(i)%energy(k)-dispgeoms(i)%energy(l))
        dener(l,k)= dener(k,l)
      end do
    end do
    print trim(fmt),i,ptWeights(i),((errGrad(i,j,k),k=1,j),j=1,nstates),&
           ((dnrm2(dispgeoms(i)%nvibs,dispgeoms(i)%grads(:,j,k),int(1))/dener(k,j),k=1,j),j=1,nstates)
  end do

  print *,""
  print *,"Gradients and their errors of hij instead of fij(g:err,ab,fit;h:err,ab,fit)"
  fmt=""
  write(fmt,'("(I5,2X,",I2,"(X,E12.5))")'),(nstates-1)*nstates*3
  do i=1,npoints
    print trim(fmt),i,( (errGradh(i,j,k),k=1,j-1) ,j=1,nstates),&
       ((dnrm2(dispgeoms(i)%nvibs,dispgeoms(i)%grads(:,j,k),int(1)), &
          k=1,j-1),j=1,nstates), &
       ((dnrm2(dispgeoms(i)%nvibs,fitG(i,:,j,k),int(1)), &
          k=1,j-1),j=1,nstates), &
                      ( (errGradh(i,k,j),k=1,j-1) ,j=1,nstates),&
       ((dnrm2(dispgeoms(i)%nvibs,dispgeoms(i)%grads(:,j,j)-dispgeoms(i)%grads(:,k,k),int(1))/2, &
          k=1,j-1),j=1,nstates),  &
       ((dnrm2(dispgeoms(i)%nvibs,fitG(i,:,j,j)-fitG(i,:,k,k),int(1))/2, &
          k=1,j-1),j=1,nstates)
  end do

  !------------------------------------------------------------------
  ! out fitting error to error.log
  !------------------------------------------------------------------
  uerrfl = getFLUnit()
  open(unit=uerrfl,file='error.log',access='sequential',form='formatted',&
     status='replace',action='write',position='rewind',iostat=ios)
  if(ios/=0)print *,"FAILED TO CREATE FILE error.log"
  write(c1,'(i2)')nstates
  write(c2,'(i4)')ncoord
  do i=1,npoints
    write(uerrfl,"("//trim(c1)//"E15.7)")    &
            fitE(i,1:nstates)-(dispgeoms(i)%energy(1:nstates))
    do j=1,nstates
      errgrd = dble(0)
      call ginv(ncoord,dispgeoms(i)%nvibs,dispgeoms(i)%bmat,ncoord,binv,1D-2)
      do k=1,dispgeoms(i)%nvibs
        if(dispgeoms(i)%scale(k)>1D-1)&
            errgrd=errgrd+binv(:,k)*(fitG(i,k,j,j)-dispgeoms(i)%grads(k,j,j))
      end do    
      write(uerrfl,"("//trim(c2)//"E15.7)") errgrd
    end do
  end do 
  close(uerrfl)

  if(printlvl>0)print *,"    deallocating arrays"
  !------------------------------------------------------------------
  ! deallocate arrays
  !------------------------------------------------------------------
  CALL cleanDIIS()
  deallocate(asol)
  deallocate(asol2)
  deallocate(NEL)
  deallocate(rhs)
  deallocate(bvec)
  deallocate(ckl)
  if(printlvl>0)print *,"Exiting makesurf()"
 return
1000 format(/,2X,"  ITERATION ",I3)
1001 format(a,f7.2,a)
1002 format(4X,"Hd coefficients solved after ",F8.2," seconds.")
1003 format(3X,"Generating initial eigenvectors.")
1004 format(3X,"Loading initial eigenvectors.")
1005 format(1x,'Iteration',i4,": Delta=",E9.2,", d[g]%=",f7.2,   &
                ", <d[g]%>=",f7.2,"%, d[E]=",f8.2,", <d[E]>=",f8.2)
1006 format(/,2x,'Computation Converged after ',i5,' Iterations')
1007 format(/,2x,'Computation did not Converge')
1014 format(/,2x,'REPRODUCTION OF LEAST-SQUARES DATA ---------')
1015 format(/,2x,'Derivative Coupling, States: [',i3,',',i3,']')
1016 format(/,2x,'Energy Gradients, State [',i3,']')
1017 format(/,2x,'Energies')
END SUBROUTINE makesurf

!
!
!
!
!
!
!
!
SUBROUTINE printSurfHeader(cons,eqs,nexact)
  use progdata, only: OUTFILE
  use makesurfdata,only: maxiter,toler,gcutoff
  IMPLICIT NONE
  INTEGER,INTENT(IN)                  :: cons
  INTEGER,INTENT(IN)                  :: eqs,nexact

  write(OUTFILE,1009)
  write(OUTFILE,1002)cons
  write(OUTFILE,1003)eqs
  write(OUTFILE,1007)nexact

  write(OUTFILE,1004)maxiter
  write(OUTFILE,1005)toler
  write(OUTFILE,1006)gcutoff
  write(OUTFILE,1000)''
  write(OUTFILE,1000)''

  return
1009 format(2x,'-----------------------------------------------------',/, &
            2x,'--------  Fitting Surface to Ab Initio Data  --------',/, &
            2x,'-----------------------------------------------------',/)
1000 format(72a)
1002 format(2x,'Number of Coefficients: ',i5)
1003 format(2x,'Number of Least Squares Equations:    ',i5)
1004 format(2x,'Maximum Number of iterations:       ',i10)
1005 format(2x,'Coefficient convergence tolerance:  ',f10.8)
1006 format(2x,'Threshold for gradient values:      ',f10.8)
1007 format(2x,'Number of Exact Equations:    ',i5,/)
end SUBROUTINE printSurfHeader

!-----------------------------------------------------------------------------------
! determined the phases of Wavefunctions that best reproduce ab initio couplings
  SUBROUTINE fixphase(nvibs,scale,fitgrad,abgrad,ckl,phaseList)
    use hddata, only: nstates
    IMPLICIT NONE
    INTEGER,intent(IN)                                        :: nvibs
    DOUBLE PRECISION,dimension(nvibs,nstates,nstates),intent(inout)           :: fitgrad
    double precision,dimension(nvibs,nstates,nstates),intent(in)              :: abgrad
    double precision,dimension(nvibs),intent(in)              :: scale
    DOUBLE PRECISION,dimension(nstates,nstates),INTENT(INOUT) :: ckl
    integer,dimension(2**(nstates-1),nstates),INTENT(IN)      :: phaseList

    INTEGER                                                   :: j,k,l,minID
    integer,dimension(nstates,nstates)                        :: phaseMat
    double precision,dimension(nstates,nstates)               :: diff
    double precision                                          :: errNorm,minerr
    minID=0
    do j=1,2**(nstates-1)
      do k=1,nstates
        phaseMat(k,1:nstates)=phaseList(j,k)*phaseList(j,1:nstates)
      end do
      errNorm=dble(0)
      do k=1,nvibs
        diff=phaseMat*fitgrad(k,:,:)-abgrad(k,:,:)
        do l=2,nstates
          errNorm=errNorm+dot_product(diff(l,1:l-1),diff(l,1:l-1))*scale(k)
        end do
      end do !k=1,nvibs
      if(minID==0.or.errNorm<minerr)then
        minID=j
        minerr=errNorm
      end if!(minID==0.or.errNorm<minerr)
    end do !j=1,2**nstates
    do k=1,nstates
      phaseMat(k,1:nstates)=phaseList(minID,k)*phaseList(minID,1:nstates)
      ckl(:,k)=phaseList(minID,k)*ckl(:,k)
    end do
    do k=1,nvibs
      fitgrad(k,:,:)=phaseMat*fitgrad(k,:,:)
    end do!k=1,nvibs
  end SUBROUTINE fixphase
!---------------------------------------------------------------------
!reorder degenerate eigenvectors to best fit ab initio gradients
!---------------------------------------------------------------------
SUBROUTINE gradOrder(ptid,fitpt,abpt,ckl,pmtList,LDP,w_en,w_grd)
  use combinatorial
  use hddata,only: nstates
  use progdata, only: abpoint,printlvl
  use makesurfdata, only: incgrad,incener,e_exact,g_exact
  IMPLICIT NONE
  type(abpoint),INTENT(IN)                                   :: abpt
  type(abpoint),INTENT(INOUT)                                :: fitpt
  DOUBLE PRECISION,DIMENSION(nstates,nstates),intent(INOUT)  :: ckl
  DOUBLE PRECISION,INTENT(IN)                                :: w_en,w_grd
  INTEGER,INTENT(IN)                                         :: ptid,LDP
  INTEGER,DIMENSION(LDP,nstates),INTENT(IN)                  :: pmtList

  integer :: i,j,ldeg,udeg,ndeg, pmt, best_p, m
  double precision :: min_err, err, err1
  DOUBLE PRECISION,DIMENSION(nstates,nstates)  :: ckl_new
  DOUBLE PRECISION,DIMENSION(abpt%nvibs)            :: diff
  DOUBLE PRECISION,DIMENSION(abpt%nvibs,nstates,nstates) :: gradnew
  double precision,dimension(nstates)          :: en_new
  double precision   ::  shift

  do i=1,fitpt%ndeggrp
    ldeg=fitpt%deg_groups(i,1)
    ndeg=fitpt%deg_groups(i,2)
    udeg=ldeg+ndeg-1
    shift = -(fitpt%energy(ldeg)+fitpt%energy(udeg)-abpt%energy(ldeg)-abpt%energy(udeg))/2
    do pmt=1,factl(ndeg)
      err=dble(0)
      do j=ldeg,udeg
        m=pmtList(pmt,j-ldeg+1)+ldeg-1
        if(incgrad(ptid,j,j).or.g_exact(ptid,j,j))then
          diff=fitpt%grads(:abpt%nvibs,m,m)-abpt%grads(:abpt%nvibs,j,j)
          err=err+dot_product(diff,diff*abpt%scale(:abpt%nvibs))*w_grd
        end if
        if(incener(ptid,j,j).or.e_exact(ptid,j,j))  err=err+ (fitpt%energy(m)-abpt%energy(j)+shift)**2 *w_en
      end do
      if(pmt==1)then
        min_err=err
        best_p=pmt
        err1=err
      else if (err<min_err)then
        min_err=err
        best_p=pmt
      end if
    end do!pmt=1,factl(ndeg)
    ckl_new=ckl
    en_new =fitpt%energy
    gradnew=fitpt%grads(:abpt%nvibs,:,:)
    do j=ldeg,udeg
      m=pmtList(best_p,j-ldeg+1)+ldeg-1
      ckl_new(:,j)=ckl(:,m)
      en_new(j)   = fitpt%energy(m)
      gradnew(:,:,j)=fitpt%grads(:abpt%nvibs,:,m)
    end do
    if(printlvl>0.and.best_p/=1)print 1000,ptid,ldeg,&
                    ldeg+ndeg-1,sqrt(err1),sqrt(min_err)
    fitpt%grads(:abpt%nvibs,:,:)=gradnew
    fitpt%energy=en_new
    do j=ldeg,udeg
      m=pmtList(best_p,j-ldeg+1)+ldeg-1
      gradnew(:,j,:)=fitpt%grads(:abpt%nvibs,m,:)
    end do
    ckl=ckl_new
    fitpt%grads(:abpt%nvibs,:,:)=gradnew
  end do !i=1,fitpt%ndeggrp
1000 format(6X,"Point ",I4," States ",I2," to ",I2,&
          " ordered by gradients. Error:",E11.5,"->",E11.5)
END SUBROUTINE

! read input parameters for makesurf
SUBROUTINE readMakesurf(INPUTFL)
  USE hddata,only: nstates
  USE progdata,only: natoms, AU2CM1
  USE makesurfdata
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: INPUTFL
  NAMELIST /MAKESURF/       npoints,maxiter,toler,gcutoff,gorder,exactTol,LSETol,outputfl,       &
                            flheader,ndiis,ndstart,enfDiab,followPrev,  &
                            w_energy,w_grad,w_fij,ediffcutoff,nrmediff,ediffcutoff2,nrmediff2,rmsexcl,     &
                            useIntGrad,intGradT,intGradS,gScaleMode,energyT,highEScale,maxd,scaleEx,    &
                            dijscale,dfstart,stepMethod,ntest,ExConv,linSteps,maxED,mmiter,flattening,     &
                            linNegSteps, gscaler,ckl_output,ckl_input,    &
                            TBas,ecutoff,egcutoff,maxRot,dijscale2, usefij, deg_cap, eshift
  npoints   = 0
  gscaler   = 1d-5
  maxiter   = 3
  followprev= .false.
  usefij    = .true.
  deg_cap   = 1d-5
  mmiter    = 10
  ecutoff   = 1d0
  egcutoff  = 6d-1
  TBas      = 1D-6
  maxRot    = 0D0
  eshift    = 0d0
  ckl_input = ''
  ckl_output= 'ckl.out'
  ntest     = 0
  linSteps  = 0
  linNegSteps = 0
  ExConv    = 1D-5
  stepMethod= 0
  dfstart   = 0
  dijscale  = 1d0
  dijscale2 = 0d0
  useIntGrad= .true.
  intGradT  = 1D-3
  scaleEx   = 1D0
  intGradS  = 1D-1
  maxd      = 1D0
  maxed     = 1D-2
  energyT   = 1D30
  highEScale = 1D0
  gScaleMode = 2
  toler     = 1D-3
  gorder    = 1D-3
  gcutoff   = 1D-14
  exactTol  = 1D-12
  LSETol    = 1D-7
  flattening= 1D-8
  outputfl  = ''
  flheader  = '----'
  ndiis     = 10
  ndstart   = 10
  enfDiab   = 0
  w_energy  = dble(1)
  w_grad    = dble(1)
  w_fij     = dble(1)
  nrmediff    = 2.D4
  ediffcutoff = 20.
  nrmediff2   = 100.
  ediffcutoff2= 1.
  rmsexcl     = 0
  read(unit=INPUTFL,NML=MAKESURF)
  if(useIntGrad)then
    nvibs=3*natoms-6
  else
    nvibs=3*natoms
  end if!useIntGrad
  if(linNegSteps<0) linNegSteps=0
  if(linSteps<0)  linSteps=0
  energyT = energyT / AU2CM1
  if(nstates>nstates.or.nstates<1)nstates=nstates
END SUBROUTINE

! load pointwise options
SUBROUTINE readDispOptions(exclEner,exclGrad,exactEner,exactGrad,exactDiff,enfGO,ptWt)
  use hddata, only:nstates
  use progdata, only: PTFL,printlvl
  use makesurfdata,only:TEqList,npoints
  IMPLICIT NONE
  TYPE(TEqList),INTENT(OUT):: exclEner,exclGrad,exactEner,exactGrad,exactDiff,enfGO
  double precision,dimension(npoints),intent(out)  :: ptWt
  character(72)            :: comment
  integer                  :: ptid,i,j,ios,tmp,lsID
  character(2)             :: job
  INTEGER,DIMENSION(6,npoints*nstates*(nstates+1)/2,3)::tpLists !EE,EG,LG,LE,LD,GO
  INTEGER,DIMENSION(6)                                ::tpCount 
  double precision   ::  wt

  ptWt=dble(1)
  if(printlvl>0)print *,"   Reading point specific options."
  open(unit=PTFL,file='points.in',access='sequential',form='formatted',&
    status='old',action='read',position='rewind',iostat=ios)
  if(ios/=0)then
    print *,"Cannot open points.in.  Skipping..."
    return
  end if!ios/=0
  read(PTFL,1000,IOSTAT=ios) comment
  read(PTFL,1000,IOSTAT=ios) comment
  tpCount=0
  do while(ios==0)
    read(PTFL,1001,IOSTAT=ios) job,ptid,i,j
    if(ios/=0)exit
    if(i>j)then
      Print *,"Warning: i>j in points.in.  i and j switched"
      tmp=i
      i=j
      j=tmp
    end if
    if(job=='WT'.or.job=='wt'.or.job=='Wt'.or.job=='wT')then
      if(ptid<0)then
        wt=-1/dble(ptid)
      else
        wt=ptid
      end if
      if(printlvl>0)print *,"    Job:",job,", Point Range",i,j,", wt=",wt
      if(i<1.or.j>npoints)then
        print *,"Point index out of range. "
        i=max(i,1)
        j=min(j,npoints)
      end if
      ptWt(i:j)=wt
      cycle
    end if
    if(printlvl>0)print *,"    Job:",job,", Point",ptid,",states",i,j
    if(ptid<1.or.ptid>npoints)then
      print *,"Point index out of range."
      cycle
    end if
    if(i>nstates.or.j<1)then
      Print *,"Warning: state id out of range in points.in. input removed"
      cycle
    end if
    if(ptid<1.or.ptid>npoints)then
      Print *,"Warning: point id out of range in points.in. input removed"
      cycle
    end if
    select case (job)
      case ('EG','eg','eG','Eg')
        lsID=1
      case ('EE','ee','eE','Ee')
        lsID=2
      case ('LG','lg','lG','Lg')
        lsID=3
      case ('LE','le','lE','Le')
        lsID=4
      case ('LD','ld','lD','Ld')
        lsID=5
      case ('GO','go','gO','Go')
        lsID=6
      case default
        print *,"Warning: control characters not recognized in points.in"
        cycle
    end select
    tpCount(lsID)=tpCount(lsID)+1
    tpLists(lsID,tpCount(lsID),1)=ptid
    tpLists(lsID,tpCount(lsID),2)=i
    tpLists(lsID,tpCount(lsID),3)=j
  end do
  call fillout(exclGrad,1)
  call fillout(exclEner,2)
  call fillout(exactGrad,3)
  call fillout(exactEner,4)
  call fillout(exactDiff,5)
  call fillout(enfGO,6)
  close(unit=PTFL)
1000 format(a72)
1001 format(a2,x,i5,i5,i5)
CONTAINS
  SUBROUTINE fillout(list,job)
    USE progdata, only:  OUTFILE
    IMPLICIT NONE
    TYPE(TEqList),INTENT(OUT):: list
    INTEGER,INTENT(IN)       :: job

    CHARACTER(40),DIMENSION(7) :: titles
    integer   :: i

    titles(1)="Gradients Excluded"
    titles(2)="Energies Excluded"
    titles(3)="Gradients Fitted Exactly"
    titles(4)="Energies Fitted Exactly"
    titles(5)="Energy Differences Fitted Exactly"
    titles(6)="Gradient Ordering Enforced"
    titles(7)="Point Excluded"

    list%length=tpCount(job)
    if(allocated(list%List))deallocate(list%List)
    allocate(list%List(tpCount(job)))
    if(tpCount(job)>0) then
      write(OUTFILE,*) trim(titles(job))
      write (OUTFILE,1004)
      write (OUTFILE,1002)
      write (OUTFILE,1004)
      do i=1,tpCount(job)
        list%List(i)%point=tpLists(job,i,1)
        list%List(i)%i    =tpLists(job,i,2)
        list%List(i)%j    =tpLists(job,i,3)
        write(OUTFILE,1003) tpLists(job,i,:)
      end do
      write (OUTFILE,1004)
      write (OUTFILE,'(/)')
    end if
  1002 format(4X,"POINT  STATE1  STATE2")
  1003 format(4X,I4,3X,I4,4X,I4)
  1004 format(4X,"---------------------")
  END SUBROUTINE
END SUBROUTINE
!
!
!
SUBROUTINE printDisps(type,npts)
  use hddata, only:  ncoord,nstates
  use progdata
  use makesurfdata
  IMPLICIT NONE
  INTEGER,INTENT(IN)              :: type
  INTEGER,INTENT(IN)              :: npts

  INTEGER                         :: i,j,k,npr
  DOUBLE PRECISION,dimension(ncoord,npts) :: disps
  CHARACTER(4)                    :: str
  CHARACTER(16),dimension(ncoord) :: rlabs
  CHARACTER(16),dimension(npts)   :: clabs

  if(type/=0)then
   write(OUTFILE,1008)
  else
   write(OUTFILE,1000)
  endif

  do i = 1,ncoord
   write(str,'(i4)')i
   rlabs(i) = ' w[ '//trim(adjustl(str))//' ]'
  enddo

  do i = 1,npts
   write(str,'(i4)')i
   clabs(i) = 'POINT'//trim(adjustl(str))
  enddo

  npr = 10

  do j = 1,npts
   do k = 1,ncoord
    disps(k,j) = dispgeoms(j)%igeom(k)
   enddo
  enddo
  call printMatrix(OUTFILE,rlabs,clabs,npr,ncoord,npts,disps,int(11),int(6))

  write(OUTFILE,'(/,2x,"Ab initio Energies")')
  write(str,'(i4)') nstates
  do j=1,npts
    write(OUTFILE,'(6x,i5,'//trim(str)//'F14.2)')j,dispgeoms(j)%energy*AU2CM1
  end do  
  return

1000 format(2x,'-----------------------------------------------------',/, &
            2x,'--------    Writing Displacements to File    --------',/, &
            2x,'-----------------------------------------------------',/)
1008 format(2x,'-----------------------------------------------------',/, &
            2x,'--------   Reading Displacements from File   --------',/, &
            2x,'-----------------------------------------------------',/)
end SUBROUTINE printDisps
!------------------------------------------------------------------------------
!Read in geometry and ab initio data for a set of points which will be used
!in the Hd fitting procedure.
!
SUBROUTINE readdisps()
  use hddata
  use progdata,only:natoms,printlvl
  use makesurfdata
  IMPLICIT NONE
  INTEGER                                      :: j,k,l
  DOUBLE PRECISION,dimension(3*natoms,npoints) :: cgrads,cgeoms
  DOUBLE PRECISION,dimension(nstates,npoints)  :: eners
  character(72)                                :: filename,infile
  character(10)                                :: suffix
  character(3),dimension(natoms)               :: atoms
  double precision,dimension(natoms)           :: anums,masses

  if(printlvl>0)print '(7X,A,I5,A)','Reading',npoints,' displacements'
  if(allocated(dispgeoms))deallocate(dispgeoms)
  allocate(dispgeoms(npoints))
  do j = 1,npoints
     allocate(dispgeoms(j)%cgeom(3*natoms))
     allocate(dispgeoms(j)%igeom(ncoord))
     allocate(dispgeoms(j)%energy(nstates))
     allocate(dispgeoms(j)%grads(3*natoms,nstates,nstates))
     allocate(dispgeoms(j)%bmat(ncoord,3*natoms))
     allocate(dispgeoms(j)%lmat(3*natoms,3*natoms))
     allocate(dispgeoms(j)%scale(3*natoms))
     allocate(dispgeoms(j)%eval(3*natoms))
  enddo

  if(printlvl>0)print 1000,'Reading geom.all'
  infile = 'geom.all'
  call readColGeom(infile,npoints,natoms,atoms,anums,cgeoms,masses)
  if(printlvl>0.and.useIntGrad)print 1000,'local internal coordinates will be constructed.'
  do j = 1,npoints
   dispgeoms(j)%id = j
   dispgeoms(j)%cgeom(1:3*natoms) = cgeoms(1:3*natoms,j)
  enddo

  if(printlvl>0)print 1000,'Reading energy.all'
  infile = 'energy.all'
  call readEner(infile,npoints,nstates,eners)
  do j = 1,npoints
   dispgeoms(j)%energy = eners(:,j)+eshift
   call genEnerGroups(dispgeoms(j),deg_cap,nstates)
  enddo

  if(printlvl>0)print 1000,'Reading gradients and couplings'
  do j = 1,nstates
   do k = 1,j
    suffix='.all'
    infile = filename(k,j,suffix,usefij)
    if(printlvl>0) print 1000,"loading COLUMBUS gradients from "//trim(adjustl(infile))
    call readGrads(infile,npoints,natoms,cgrads)
    do l = 1,npoints
      if(j/=k.and.usefij)then
        dispgeoms(l)%grads(:,j,k)=cgrads(:,l)*(eners(j,l)-eners(k,l))
      else
        dispgeoms(l)%grads(:,j,k)=cgrads(:,l)
      end if
  !    call removeTransRot(dispgeoms(l)%grads(:,j,k),dispgeoms(l)%cgeom)
      if(j/=k)dispgeoms(l)%grads(:,k,j)=dispgeoms(l)%grads(:,j,k)
    enddo
   enddo
  enddo

  if(printlvl>0)print 1000,"Generating displacement Wilson's B-Matrices"
  if(printlvl>0.and.useIntGrad)print 1000,'Local internal coordinates will be constructed from B-matrices'
  do l = 1,npoints
    if(printlvl>0)print *,""
    if(printlvl>0)print *,"      Constructing Wilson B matrix at point ",l," :"
    call buildWBMat(dispgeoms(l)%cgeom,dispgeoms(l)%igeom,dispgeoms(l)%bmat)
    if(printlvl>0)then
        print *,"  Internal Coordinates"
        print "(15F8.3)",dispgeoms(l)%igeom
    end if
    if(printlvl>3)then
        print *,"  Wilson B Matrix in nascent coordinate system"
        do j=1,ncoord
            print "(15F8.3)",dispgeoms(l)%bmat(j,:)
        end do
    end if
    if(printlvl>0)print *,"      Constructing local coordinate system for point ",l
    CALL makeLocalIntCoord(dispgeoms(l),nstates,useIntGrad,intGradT,intGradS,nvibs,gScaleMode)
    if(printlvl>1)then
       print *,"  Wilson B Matrix in fitting coordinate system"
       do j=1,ncoord
         print "(15F8.3)",dispgeoms(l)%bmat(j,:)
       end do
       do j=1,nstates
         print *,"  Gradients for state",j," in fitting coordinate system"
         print "(15F8.3)",dispgeoms(l)%grads(:,j,j)
       end do
       do j=1,nstates-1
         do k=j+1,nstates
           print *,"  Couplings for block",j,k," in fitting coordinate system"
           print "(15F8.3)",dispgeoms(l)%grads(:,j,k)
         end do
       end do
    end if
  end do

 ! generate intersection adapated coordinates
  do l=1,npoints
   call OrthGH_ab(dispgeoms(l),100,sqrt(gcutoff)/10)
  end do

  call printDisps(int(1),npoints)

  if(allocated(ptWeights))deallocate(ptWeights)
  allocate(ptWeights(npoints))
  call readDispOptions(exclEner,exclGrad,exactEner,exactGrad,exactDiff,enfGO,ptWeights)

1000 format(9X,A)
end SUBROUTINE readdisps
