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
!
!Input File Options
!
!    The user can provide a series of search paths, and the program will look for
!  energy, gradient and coupling input in each of these search paths, if it exists.
!    If a file is missing in that directory, these data will be marked as nonexist
!  and will not be included in the fit set.
!    This feature is added so that you don't have to make up data then specify the
!  piece of data to be nonexist in points.in.   You may also opt to change the file
!  name patterns for each input file but normally you don't have to.   
!    Note that these file name patterns have special characters such as $ denoting
!  index of an involved state.
!    Also note that there is an optional note file that will be read-in and its 
!  content will be included in output files as comments.
  INTEGER,parameter                            :: MaxSearchPaths=999 ! Maximum number of search pathes that can be provided by user
  INTEGER                                      :: NSearchPaths
  CHARACTER(1000),dimension(MaxSearchPaths)    :: SPNotes
  CHARACTER(255),dimension(MaxSearchPaths)     :: SearchPath    ! The program will look for input files in these directories
  LOGICAL,dimension(:,:),allocatable           :: hasEner       ! specifies the presence of data at each point
  LOGICAL,dimension(:,:,:),allocatable         :: hasGrad       ! specifies the presence of data at each point
  
  CHARACTER(255),dimension(MaxSearchPaths)      :: notefptn,gmfptn,enfptn,grdfptn,cpfptn  ! naming pattern of input files.

! General options
  DOUBLE PRECISION                             :: eshift    ! uniform shift on ab initio energies
  DOUBLE PRECISION                             :: deg_cap   ! threshold below which intersection adapted coordinates will be used
  INTEGER                                      :: maxiter
  INTEGER                                      :: npoints
  DOUBLE PRECISION                             :: toler     !criteria for convergence
  DOUBLE PRECISION                             :: gcutoff   !gradient cutoff
  DOUBLE PRECISION                             :: exactTol  !Exact equations rank cutoff
  DOUBLE PRECISION                             :: LSETol    !LSE rank cutoff
  DOUBLE PRECISION                             :: flattening!flattening parameter for differential convergence
!gorder: energy difference below which ckl will be ordered according to gradients rather than energies
  DOUBLE PRECISION                             :: gorder
  CHARACTER(255)                               :: outputfl  !name of the output file generated after fit
  CHARACTER(255)                               :: flheader  !comment line in the outputfile
  CHARACTER(255)                               :: ckl_input   !if nonempty, read wave functions from this file
  CHARACTER(255)                               :: ckl_output  !if nonempty, output wave functions to this file
  CHARACTER(255)                               :: guide       !if nonempty, ordering guide wave functions to this file
  LOGICAL                                      :: orderall    !tell the ordering to use data that are excluded in points.in

! restartdir specifies a directory where hd coefficients for each iteration will
! be saved. particularly useful if the error turns up at one point and you want
! to restart from a specific iteration, or when the job was interrupted. 
  CHARACTER(255)                               :: restartdir  !if nonempty, output wave functions to this file

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

! whether to generate error and geometry information for error analysis in
! potlib libraries. geometries will be stored in `refgeom`, and error info will
! be stored in `error.log`
  LOGICAL                                      :: printError

  TYPE(abpoint),dimension(:),allocatable       :: dispgeoms
  type(TEqList)                                :: exclEner,exclGrad,exactEner,exactGrad,exactDiff,enfGO
  INTEGER                                      :: enfDiab ! index of point where diabatic and adiabatic matches
  DOUBLE PRECISION,dimension(:),allocatable    :: ptWeights

  type(TDList),dimension(:),allocatable        :: WVals
  type(TDList),dimension(:),allocatable        :: dWVals    
  type(T2DDList),dimension(:),allocatable      :: WMat
  DOUBLE PRECISION,dimension(:,:,:),allocatable:: ckl,   cklguide

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

! diagonal hessian approximation.  when >0, coef-coef block of the Lagrangian
! hessian is approximated by an identity matrix times this number
  DOUBLE PRECISION                             ::  diagHess

! maximum allowed change in coefficients in each iteration
  DOUBLE PRECISION  :: maxd

! scaling factors to exact equations
  DOUBLE PRECISION  :: scaleEx

! energy above which gradients will be automatically removed from fitting
! equations
  DOUBLE PRECISION  :: gradcutoff

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
  double precision                           :: ecutoff  ! energy threshold above which data will be excluded from basis construction
  double precision                           :: egcutoff ! energy threshold above which gradients will no longer be used in fit
  type(T2DDList),dimension(:),allocatable    :: ZBas ! transformation from primitive to the reconstructed basis for each block

! scaling factor of dij term.  0 gives old result and 1 gives exact gradients
  double precision            :: DijScale

! index of the iteration to start differential convergence
  integer    :: dfstart

! number of unknown coefficients after null-space reduction (actually used in fit)
  integer    :: ncons
! number of unknown coefficients before null-space reduction (hd storage)
  integer    :: ncon_total
! number of least squares fitting equations
  integer    :: neqs
! number of exact equations
  integer    :: nex

! number of linear steps refinements 
  integer    :: linSteps
! relative step length convergence threshold for linear step refinement
  double precision  :: dconv

! Total norm of gradients by all basis function for physical gradient component
! at each specific data point
  DOUBLE PRECISION,dimension(:,:,:),allocatable     :: gradNorm

! Weights for least squars equations
  DOUBLE PRECISION,dimension(:),allocatable         :: weight

!Hd predictions for all data points
  DOUBLE PRECISION,dimension(:,:,:,:),allocatable   :: fitG
  DOUBLE PRECISION,dimension(:,:,:),allocatable     :: fitE
 CONTAINS

  ! determine if each equation will be include / excluded / fitted exactly
  SUBROUTINE getPtList()
    use hddata, only: nstates
    use progdata, only: printlvl
    IMPLICIT NONE
    INTEGER i,j,s1,s2

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
      if(ptWeights(i)>1D-8.or.orderall==.false.)then
        do s1=1,nstates
          do s2=1,s1
            if(dispgeoms(i)%energy(s1,s1)<gradcutoff)incgrad(i,s1,s2)=hasGrad(i,s1,s2)
          end do
          incener(i,s1,s1)=hasEner(i,s1)
          if(i==enfDiab)incener(i,:,:)=.true.
        end do!s1=1,nstates
        ! add off-diagonal elements between states within degeneracy groups
        do j=1,dispgeoms(i)%ndeggrp
          s1 = dispgeoms(i)%deg_groups(j,1)
          s2 = s1-1+dispgeoms(i)%deg_groups(j,2)
 ! If all energies are available in the degeneracy group, add off-diagonal elements to the fitting set
 ! because these states are subject to rotations to form intersection adapted coordinates 
          if(all(hasEner(i,s1:s2)))  incener(i,s1:s2,s1:s2)=.true.
        end do !j
      end if!pWeights(i)>1D-8
    end do!i=1,npoints
    do i=1,exclEner%length
      s1=exclEner%List(i)%i
      s2=exclEner%List(i)%j
      incener(exclEner%List(i)%point,s1:s2,s1:s2)=.false.
    end do
    do i=1,exclGrad%length
      incgrad(exclGrad%List(i)%point,exclGrad%List(i)%i,exclGrad%List(i)%j)=.false.
      incgrad(exclGrad%List(i)%point,exclGrad%List(i)%j,exclGrad%List(i)%i)=.false.
    end do
    do i=1,exactEner%length
     do s1=exactEner%List(i)%i,exactEner%List(i)%j
       do s2=exactEner%List(i)%i,exactEner%List(i)%j
          if(hasEner(exactEner%List(i)%point,s1).and.hasEner(exactEner%List(i)%point,s2))then
            e_exact(exactEner%List(i)%point,s1,s2)=.true.
          else
            print "(A,I6,2(A,I3),A)","WARNING : Cannot enforce exact energy at point ",&
                exactEner%List(i)%point,", block(",s1,",",s2,"). Ab initio data not present."
          end if
       end do
     end do
    end do
    do i=1,exactGrad%length
      if(hasGrad(exactGrad%List(i)%point,exactGrad%List(i)%i,exactGrad%List(i)%j))then
        g_exact(exactGrad%List(i)%point,exactGrad%List(i)%i,exactGrad%List(i)%j)=.true.
        g_exact(exactGrad%List(i)%point,exactGrad%List(i)%j,exactGrad%List(i)%i)=.true.
      else
        print "(A,I6,2(A,I3),A)","WARNING : Cannot enforce exact gradients at point ",&
            exactGrad%List(i)%point,", block(",exactGrad%List(i)%i,",",&
            exactGrad%List(i)%j,"). Ab initio data not present."
      end if
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
!Exact equations will be sorted by point index
  SUBROUTINE makeEqMap(maxEqs,lseMap,exactEqs,wvec)
    use hddata, only: blkmap,nstates
    use progdata, only: AU2CM1
    implicit none
    integer,intent(in)                      :: maxEqs
    integer,dimension(MaxEqs,4),intent(out) :: lseMap,exactEqs
    double precision,dimension(MaxEqs),intent(out)   :: wvec
    integer   ::  tmp(4)

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
      do s2= 1,s1     
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
                    do while(dispgeoms(i)%energy(s1,s1)>energyT(k+1))  !  determine the bracket of current energy
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
                    do while(dispgeoms(i)%energy(s1,s1)+dispgeoms(i)%energy(s2,s2)> 2*energyT(k+1))
                                 !  determine the bracket of current energy
                      k = k+1
                      if(k==10)exit
                    end do
                    if(k>0) wvec(nEqs) =  wvec(nEqs)*highEScale(k)
                    if(nrmediff>0)then
                      ediff=abs(dispgeoms(i)%energy(s1,s1)-dispgeoms(i)%energy(s2,s2))*AU2CM1
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
                do while(dispgeoms(i)%energy(s1,s1)>energyT(k+1))
                  k=k+1
                  if(k==10)exit
                end do 
                if(k>0)   wvec(nEqs) =  wvec(nEqs)*highEScale(k)
                if(s1==s2 .and. s2>0 .and. nrmediff2>0)then
                  if(s1>1)then
                    ediff=abs(dispgeoms(i)%energy(s1,s1)-dispgeoms(i)%energy(s1-1,s1-1))*AU2CM1
                    ediff=(ediff+ediffcutoff2)/nrmediff2
                    if(ediff<1D0)wvec(nEqs)=wvec(nEqs)/ediff
                  end if
                  if(s1<nstates)then
                    ediff=abs(dispgeoms(i)%energy(s1+1,s1+1)-dispgeoms(i)%energy(s1,s1))*AU2CM1
                    ediff=(ediff+ediffcutoff2)/nrmediff2
                    if(ediff<1D0)wvec(nEqs)=wvec(nEqs)/ediff
                  end if!(s1<nstates)
                end if!(s1==s2 .and. s2>0 .and. nrmediff2>0)
              end if
            end if
          end do !s2=1,nstates
        end do !s1=1,nstates
    end do!i=1,npoints
    ! sort exact equations by point index
    do i=1,nex-1
      do j=nex,i+1,-1
        if(exactEqs(i,1)>exactEqs(j,1))then
          tmp           = exactEqs(i,:)
          exactEqs(i,:) = exactEqs(j,:)
          exactEqs(j,:) = tmp
        end if
      end do!j
    end do!i
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
      SUBROUTINE fixphase(nvibs,scale,fitgrad,abgrad,ckl,phaseList,hasGrad)
        use hddata, only: nstates
        IMPLICIT NONE
        INTEGER,intent(IN)                                        :: nvibs
        DOUBLE PRECISION,dimension(nvibs,nstates,nstates),intent(inout)           :: fitgrad
        double precision,dimension(nvibs,nstates,nstates),intent(in)              :: abgrad
        double precision,dimension(nvibs),intent(in)              :: scale
        DOUBLE PRECISION,dimension(nstates,nstates),INTENT(INOUT) :: ckl
        LOGICAL,dimension(nstates,nstates),intent(in)             :: hasGrad
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
    double precision,dimension(nstates)              :: eval,ovlp

    allocate(pmtList(factl(nstates),nstates))
    pmtList=Permutation(nstates,factl(nstates))
    folPrev = .false.
    if(present(follow))folPrev=follow
    if(printlvl>0) print *,"     Updating wave functions"
    allocate(fitpt%energy(nstates,nstates))
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
      CALL EvaluateHd3(hvec,nBas,npoints,i,nvibs,hmatPt,dhmatPt,WMat)
      if(i==enfDiab)then
        ckl(i,:,:)=0d0
        do k=1,nstates
          ckl(i,k,k) = 1d0
        end do!k
        fitG(i,:,:,:)=dhmatPt
        fitE(i,:,:) = hmatPt
        fitG(i,dispgeoms(i)%nvibs+1:nvibs,:,:) = 0d0
        call fixphase(dispgeoms(i)%nvibs,dispgeoms(i)%scale(:dispgeoms(i)%nvibs),fitG(i,1:dispgeoms(i)%nvibs,:,:),&
                 dispgeoms(i)%grads(1:dispgeoms(i)%nvibs,:,:),ckl(i,:,:),phaseList,hasGrad(i,:,:))
        cycle
      end if!i==enfDiab
      fitpt%lb = dispgeoms(i)%lb
      fitpt%ub = dispgeoms(i)%ub
      fitpt%id = i
      cklPrev = ckl(i,:,:)
      ckl(i,:,:)=hmatPt
      call DSYEV('V','U',nstates,ckl(i,:,:),nstates,eval,scr,5*nstates*nstates,INFO)
      IF(INFO/=0)then
        print *,"Failed to solve eigenvectors at point",i
        print *,"INFO=",info

        stop 'error solving eigenvalue equations'
      end if
      call OrthGH_Hd(dispgeoms(i),dhmatPt(:dispgeoms(i)%nvibs,:,:),ckl(i,:,:),100,gcutoff,hasGrad(i,:,:))
      do k=1,dispgeoms(i)%nvibs
        fitpt%grads(k,:,:)=matmul(transpose(ckl(i,:,:)),matmul(dhmatPt(k,:,:),ckl(i,:,:)))
      end do!k=1,dispgeoms(i)%nvibs
      fitE(i,:,:)=matmul(transpose(ckl(i,:,:)),matmul(hmatPt,ckl(i,:,:)))
      fitpt%energy=fitE(i,:,:)
      do k=1,enfGO%Length
        if(enfGO%List(k)%point==i)then
          do j=enfGO%List(k)%i+1,enfGO%List(k)%j
            fitpt%energy(j,j)=fitpt%energy(enfGO%List(k)%i,enfGO%List(k)%i)
          end do!j=enfGO%List(k)%i+1,enfGO%List(k)%j
        end if!(enfGO%List(k)%point==i)
      end do!k=1,enfGO%Length
      if(folPrev)then
! determine ordering and phase by consistency with previous iteration
        call folPrevCkl(ckl(i,:,:),cklPrev,pmtList,factl(nstates),nstates,i)
        do k=1,dispgeoms(i)%nvibs
          fitG(i,k,:,:)=matmul(transpose(ckl(i,:,:)),matmul(dhmatPt(k,:,:),ckl(i,:,:)))
        end do!k=1,dispgeoms(i)%nvibs
        fitG(i,dispgeoms(i)%nvibs+1:nvibs,:,:) = 0d0
      else
        if(cklguide(i,1,1)<-1d2.or.dispgeoms(i)%ndeggrp>0)then
! determine ordering and phase by comparing with ab initio data
          call genEnerGroups(fitpt,gorder/AU2CM1)
          call gradOrder(i,fitpt,dispgeoms(i),ckl(i,:,:),pmtList,factl(nstates),w_energy,w_grad,w_fij,nrmediff,ediffcutoff)
        else
! determine ordering by matching guide points
          call guideOrder(i,ckl(i,:,:),cklguide(i,:,:),pmtlist,factl(nstates))
        end if
        call fixphase(dispgeoms(i)%nvibs,dispgeoms(i)%scale(:dispgeoms(i)%nvibs),fitpt%grads(1:dispgeoms(i)%nvibs,:,:),&
                 dispgeoms(i)%grads(1:dispgeoms(i)%nvibs,:,:),ckl(i,:,:),phaseList,hasGrad(i,:,:))
        fitpt%grads(dispgeoms(i)%nvibs+1:,:,:)=dble(0)
        fitG(i,:,:,:)=fitpt%grads
      end if !folPrev
      fitE(i,:,:)=matmul(transpose(ckl(i,:,:)),matmul(hmatPt,ckl(i,:,:)))
      if(ptWeights(i)<1d-8.and.(.not.(printlvl>3)))cycle
      do k=1,nstates
        ovlp(k) = abs(dot_product(cklPrev(:,k),ckl(i,:,k)))
      end do!k
      if(sum(cklPrev**2)>8d-1.and.minval(ovlp)<8d-1)then
        print "(A,I5)","small overlap between iterations at pt ",i
        print *,"before:"
        do k=1,nstates
          print "(10F10.5)",cklPrev(:,k)
        end do
        print *,"after:"
        do k=1,nstates
          print "(10F10.5)",ckl(i,:,k)
        end do
      end if!
    end do!i=1,npoints
    deallocate(fitpt%energy)
    deallocate(fitpt%grads)
    if(allocated(fitpt%deg_groups)) deallocate(fitpt%deg_groups)
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
! rhs is constructed for diff=.false.
! ndep is the number of linear dependencies among exact equations
  SUBROUTINE makeNormalEquations(NE,rhs,wt,ndep)
    use hddata, only: nblks,nstates,nBasis
    use progdata, only: printlvl
    IMPLICIT NONE
    DOUBLE PRECISION,dimension(neqs),INTENT(IN)                 :: wt
    DOUBLE PRECISION,dimension(nex+ncons,nex+ncons),INTENT(OUT) :: NE
    DOUBLE PRECISION,dimension(nex+ncons),INTENT(OUT)           :: rhs
    INTEGER,INTENT(OUT)                                         :: ndep
   
    integer,parameter ::  bars = 80
    integer           ::  i,j,  pc, pc_last, nEqPt,nExPt,  neqTot, nexTot
    double precision,dimension(:,:),allocatable  ::  AMat
    double precision,dimension(:),allocatable    ::  bvec
    integer           ::  LDA

    ! initialization.   Only a strip of A matrix is stored
    LDA = nstates**2*(nvibs+2)
    allocate(AMat(LDA,ncons))
    allocate(bvec(LDA))
    NE  = 0d0
    rhs = 0d0
    neqTot = 0
    nexTot = 0

    !print the progress bar
    if(printlvl>0)then
      print *,"   Constructing normal equations"
      write (*,'(4x)',advance='no')
      do i=1,bars
        write(*,'(a)',advance='no'),"_"
      end do
      print *,""
      write (*,'(a,$)')"    "
    end if

    !main look, add up contributions from each data point
    pc_last=0
    do i=1,npoints
      ! construct the strip of a matrix for the current point
      call makeCoefMatPt(i,AMat,LDA,nEqPt,nExPt,wt,bvec)
      ! fill in exact block of rhs and normal equations
      do j=1,nExPt
        NE(j+nexTot,nex+1:nex+ncons)= AMat(j,:)*scaleEx
        rhs(j+nexTot)               = bvec(j)  *scaleEx
      end do!j
      ! accumulate rhs vector
       CALL DGEMV('T',nEqPt,ncons,1d0,AMat(nExPt+1,1),LDA,&
          bvec(nExPt+1),int(1),dble(1),rhs(nex+1),INT(1))

      ! fill in the LSE block
      CALL DSYRK('U','T',ncons,nEqPt,1d0,AMat(nExPt+1,1),&
         LDA,dble(1),NE(nex+1,nex+1),ncons+nex)
      
      neqTot = neqTot+nEqPt
      nexTot = nexTot+nExPt  
      ! print the progression bar
      if(printlvl>0)then
        pc = i*bars/npoints  ! check percentage finished and print out to screen
        if(pc>pc_last)then
           write(*,"(A,$)"),'>'
           pc_last = pc
        end if!pc>pc_last
      end if!printlvl>0
    end do!i
    if(printlvl>0)print *,"      [ DONE ]"
    deallocate(AMat)
    deallocate(bvec)
! reform exact block to remove linear dpendencies
    ndep = 0
    if(nex>0)then
      if(printlvl>0)print *,"  Removing linear dependencies from exact equations..."
      call removeDependency(NE,rhs,ndep)
      if(printlvl>0)print *,"  Number of linear dependencies in exact equations: ",ndep 
    end if
  END SUBROUTINE makeNormalEquations
!-----------------------------------------------------------------------------------
! remove linear dependencies from the exact equations block and return the
! number of such linear dependencies.
  subroutine removeDependency(Mat, rhs,ndep)
    implicit none
    double precision,dimension(nex+ncons,nex+ncons),intent(INOUT):: Mat
    double precision,dimension(nex+ncons),intent(INOUT)          :: rhs
    integer,intent(OUT)                                          :: ndep

    double precision,dimension(ncons,nex) :: newEqs
    double precision,dimension(nex)       :: newRhs
    double precision,dimension(nex,nex) :: ovlp !overlap matrix between exact equations
    double precision,dimension(2*nex*(nex+2)) ::work
    double precision,dimension(nex) :: eval
    integer :: lwork,info,i,ngood 
    
    lwork =2*nex*(nex+2)
    call dsyrk('U','N',nex,ncons,1d0,Mat(1,nex+1),ncons+nex, 0d0, ovlp,nex) 
    call dsyev('V','U',nex,ovlp,nex,eval,work,lwork,info)
    if(info/=0)then
      print *,"failed to do eigenvalue decomposition.  info=",info
      stop "fatal error in makesurf::removeDependency()"
    end if
    ngood = 0
    do i=1,nex
      if(eval(i)>exactTol)then
        ngood=ngood+1
        call dgemv('T',nex,ncons,1d0,Mat(1,nex+1),nex+ncons,ovlp(1,i),1,0d0,newEqs(1,ngood),1)
        newrhs(ngood) = dot_product(rhs(1:nex),ovlp(:,i))
      end if
    end do
    !reform the rhs and equations
    do i=1,ngood
      Mat(nex+1-i,nex+1:)=newEqs(:,i)
      rhs(nex+1-i)       =newRhs(i)
    end do
    ndep = nex-ngood
  end subroutine removeDependency
!-----------------------------------------------------------------------------------
!Generate the matrix A of the linear coefficients between Hd basis and RHS values
!for a specific data point
!Also returns the number of exact and LSE equations for that point
  SUBROUTINE makeCoefMatPt(pt,AMat,LDA,nEqPt,nExPt,wt,bvec)
    use hddata, only: nblks,nstates,RowGrp,ColGrp,offs,grpLen,nBasis
    IMPLICIT NONE
    INTEGER,INTENT(IN)                                  :: pt, LDA
    DOUBLE PRECISION,dimension(LDA,ncons),INTENT(OUT)   :: AMat
    INTEGER,INTENT(OUT)                                 :: nEqPt,nExPt
    DOUBLE PRECISION,dimension(LDA),INTENT(OUT)         :: bvec
    DOUBLE PRECISION,dimension(neqs),INTENT(IN)         :: wt
    INTEGER      :: l,r,l0,r0,s1,s2,g,iBlk,ll,rr,offset,i
    double precision  :: ediff(nstates,nstates),error
    double precision,dimension(:,:,:),allocatable    :: WIJ,DIJ 
    double precision,dimension(:,:,:,:),allocatable  :: VIJ
    double precision,dimension(:),allocatable        :: dv
    double precision  ::  cc
    INTEGER           ::  Id(nstates**2),Jd(nstates**2),ndeg,ind1,pv1
    DOUBLE PRECISION  ::  gvec(nvibs),hvec(nvibs)
    double precision,external :: dnrm2
    INTEGER           ::  ptEqMap(LDA,3)
    double precision  ::  wtv(LDA)

! construct degeneracy list and energy difference table
    ndeg = 0
    if(enfDiab.ne.pt)then
      do s1=1,nstates
        do s2=s1+1,nstates
          ediff(s1,s2) = fitE(pt,s2,s2)-fitE(pt,s1,s1)
        end do
      end do
    end if

    nEqPt = 0 
    nExPt = 0 
! construct equations map for the point
! exact equations
    do i=1,nex
      if(eqMap(i,1).ne.pt)cycle
      nExPt = nExPt+1
      if(nExPt>LDA) stop "Error! makeCoefMatPt: Lower dimension of matrix too small"
      ptEqMap(nExPt,1:3)=eqMap(i,2:4) 
    end do
! least squares equations
    do i=1,neqs
      if(eqMap(i+nex,1).ne.pt)cycle
      nEqPt = nEqPt+1
      if(nExPt+nEqPt>LDA) stop "Error! makeCoefMatPt: Lower dimension of matrix too small"
      ptEqMap(nExPt+nEqPt,1:3)=eqMap(i+nex,2:4) 
      wtv(nEqPt)=wt(i)
    end do
    if(nExPt+nEqPt==0)return

! construct Wij, Dij and Vij matrices for all blocks
    error = dble(0)
    offset = 0
    AMat = 0d0
    do iBlk = 1,nblks
     l0 = offs(RowGrp(iblk))
     ll = GrpLen(RowGrp(iblk))
     r0 = offs(ColGrp(iblk))
     rr = GrpLen(ColGrp(iblk))
     pv1 =  npoints*(1+nvibs) *ll*rr
     allocate(VIJ(nBas(iBlk),nvibs,nstates,nstates))    
     allocate(WIJ(nBas(iBlk),nstates,nstates))    
     allocate(DIJ(nBas(iBlk),nstates,nstates))    
     allocate( dv(nBas(iBlk)))
     WIJ = dble(0)
     VIJ = dble(0)
     DIJ = dble(0)
     do s1=1,nstates
       do s2=s1,nstates
         do r= 1, rr
           do l= 1, ll
   ! calculate the multiplication of ckl, with Hermitianization for off diagonal blocks
              if(RowGrp(iblk).eq.ColGrp(iblk))then
                cc = ckl(pt,l+l0,s1)*ckl(pt,r+r0,s2)
              else !(RowGrp(iblk).eq.ColGrp(iblk))
                cc = ckl(pt,l+l0,s1)*ckl(pt,r+r0,s2) + ckl(pt,r+r0,s1)*ckl(pt,l+l0,s2)
              end if !(RowGrp(iblk).eq.ColGrp(iblk))

   ! get the values of W and V
              ind1 = (((r-1)*ll+l-1)*npoints+pt-1)*(nvibs+1)+1
              CALL DAXPY(nBas(iBlk),cc,WMat(iBlk)%List(ind1,1),pv1,WIJ(1,s1,s2),1)
   ! construct the Hamiltonian contribution of VIJ 
              do i=1,nvibs
                ind1 = ind1+1
                CALL DAXPY(nBas(iBlk),cc,WMat(iBlk)%List(ind1,1),pv1,VIJ(1,i,s1,s2),1)
              end do!i=1,nvibs
           end do !  l= 1, ll
         end do !r= 1, rr

         if(s1.ne.s2) then
           if(abs(DijScale)>1D-10.and.pt/=enfDiab)then
  ! construct DIJ from WIJ.  It is non zero only on the off diagonal
  !no rotation for reference point
             if(dispgeoms(pt)%ndeggrp==0) DIJ(:,s1,s2)=WIJ(:,s1,s2)/ediff(s1,s2)*DijScale
           end if!Dij/=0 and not enforcing diabat
  ! fill the lower triangle of D and W
           WIJ(:,s2,s1)     =  WIJ(:,s1,s2)
           DIJ(:,s2,s1)     = -DIJ(:,s1,s2)
         end if !(s1.ne.s2)
       end do!s2
     end do!s1

  ! construct DIJ contribution to degenerate points
     if(abs(DijScale)>1D-10.and.dispgeoms(pt)%ndeggrp>0.and.pt/=enfDiab)then
        do l=1,nBas(iBlk)
          call getDegDij(pt,WIJ(l,:,:),VIJ(l,:,:,:),DIJ(l,:,:))
        end do
     end if !DijScale /= 0
  ! the Dij Contribution of VIJ
     do s1=1,nstates
       do s2=s1,nstates
         do l=1,nstates
            if(abs(fitE(pt,l,s2))>1d-20)WIJ(:,s1,s2)=WIJ(:,s1,s2)+DIJ(:,l,s1)*fitE(pt,l,s2)
            if(abs(fitE(pt,l,s1))>1d-20)WIJ(:,s1,s2)=WIJ(:,s1,s2)+DIJ(:,l,s2)*fitE(pt,l,s1)
         end do
         if(s1.ne.s2)WIJ(:,s2,s1)=WIJ(:,s1,s2)
         do g=1,nvibs
            do l=1,nstates
              if(l.ne.s1) call DAXPY(nBas(iBlk),fitG(pt,g,l,s2),DIJ(1,l,s1),1,VIJ(1,g,s1,s2),1)
              if(l.ne.s2) call DAXPY(nBas(iBlk),fitG(pt,g,l,s1),DIJ(1,l,s2),1,VIJ(1,g,s1,s2),1)
      !        VIJ(:,g,s1,s2) = VIJ(:,g,s1,s2) + DIJ(:,l,s1)*fitG(pt,g,l,s2) + DIJ(:,l,s2)*fitG(pt,g,l,s1)     
            end do
         end do!g
         if(s1.ne.s2)VIJ(:,:,s2,s1)=VIJ(:,:,s1,s2)
       end do !s2
     end do ! s1
! Make AMAT from W and V matrices
     do i=1,nEqPt+nExPt
        s1=ptEqMap(i,1)    ! state index 1
        s2=ptEqMap(i,2)    ! state index 2
        g =ptEqMap(i,3)    ! gradient component index
        if(g==0)then !iGrad==0  => its an energy fit
          do l=1,GrpLen(RowGrp(iblk))
            do r=1,GrpLen(ColGrp(iblk))
              if(s2<0) then! equation for energy difference
                 AMat(i,offset+1:offset+nBas(iBlk))=WIJ(:,s1,s1)-WIJ(:,-s2,-s2)
              else         ! equation for energy 
                 AMat(i,offset+1:offset+nBas(iBlk))=WIJ(:,s1,s2)
              end if !s2<0
            end do!r
          end do!l
        else                                      ! equation for energy gradients or derivative couplings
          AMat(i,offset+1:offset+nBas(iBlk))=VIJ(:,g,s1,s2)
        end if!if(eqMap(i,4)==0)
     end do!i=1,nEqPt+nExPt
     offset=offset+nBas(iBlk)
     deallocate(dv)
     deallocate(VIJ)
     deallocate(WIJ)
     deallocate(DIJ)
    end do! iblk
    ! construct bvec
    do i=1,nEqPt+nExPt
        s1=ptEqMap(i,1)    ! state index 1
        s2=ptEqMap(i,2)    ! state index 2
        g =ptEqMap(i,3)    ! gradient component index
        if(g==0)then !iGrad==0  => its an energy fit
          if(s2<0) then! equation for energy difference
            bvec(i)=dispgeoms(pt)%energy(s1,s1)-dispgeoms(pt)%energy(-s2,-s2)
          else         ! equation for energy 
            bvec(i)=dispgeoms(pt)%energy(s1,s2)
          end if !s2<0
        else            ! equation for energy gradients or derivative couplings
          bvec(i)=dispgeoms(pt)%grads(g,s1,s2)
        end if!g==0
    end do!i
    ! scale with weight factors
    do i=1,nEqPt
      CALL DSCAL(ncons,wtv(i),AMat(i+nExPt,1),LDA)
      bvec(nExPt+i)=bvec(nExPt+i)*wtv(i)
    end do
  end SUBROUTINE makeCoefMatPt
!-----------------------------------------------------------------------------------
!generate rhs vectors for Hd fitting
  SUBROUTINE makebvec(bvec,diff)
    IMPLICIT NONE
    DOUBLE PRECISION,dimension(nex+neqs),INTENT(INOUT)             :: bvec
    LOGICAL, INTENT(IN)                                            :: diff
    INTEGER                                                :: i,pt,s1,s2,g
    do i = 1,nex+neqs
      pt=eqmap(i,1)
      s1=eqmap(i,2)
      s2=eqmap(i,3)
      g =eqmap(i,4)
      if(g==0)then !is an energy
        if(s2>0)then !its adiabatic energy or adiabatic off diagonal terms(0)
            bvec(i)=dispgeoms(pt)%energy(s1,s2)
            if(diff)bvec(i)=bvec(i)-fitE(pt,s1,s2)
        else  !it is adiabatic energy difference
          bvec(i)=dispgeoms(pt)%energy(s1,s1)-dispgeoms(pt)%energy(-s2,-s2)
          if(diff)bvec(i)=bvec(i)-fitE(pt,s1,s1)+fitE(pt,-s2,-s2)
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
   DOUBLE PRECISION  :: gnrm,dgrd,rmsexclthreshold,dE_exact,dG_exact,dcp,nrmcp,nrmdcp,ncp,de1,de2,dcp2
   double precision,external :: dnrm2

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
           if(hasGrad(j,k,l).and.l>k.and.(incgrad(j,k,l).or.incgrad(j,l,k)))then
             de1 = max(abs(dispgeoms(j)%energy(k,k)-dispgeoms(j)%energy(l,l)),1D-7 )
             de2 = max(abs(fitE(j,k,k)-fitE(j,l,l)),1D-7 )
             inc_cp=inc_cp+1
             dcp    =  dot_product(  &
                   dispgeoms(j)%grads(:nvpt,k,l)/de1-fitG(j,:nvpt,k,l)/de2, &
                 ( dispgeoms(j)%grads(:nvpt,k,l)/de1-fitG(j,:nvpt,k,l)/de2)*dispgeoms(j)%scale(1:nvpt)  )
             ncp = dot_product(dispgeoms(j)%grads(:nvpt,k,l)                          , &
                               dispgeoms(j)%grads(:nvpt,k,l)*dispgeoms(j)%scale(:nvpt) ) /de1**2
             dcp2 =dnrm2(nvpt,(dispgeoms(j)%grads(:nvpt,k,l)-fitG(j,:nvpt,k,l))*sqrt(dispgeoms(j)%scale(:nvpt)),1)/&
                        dnrm2(nvpt,dispgeoms(j)%grads(:nvpt,k,l)*sqrt(dispgeoms(j)%scale(:nvpt)),1)
             if(dcp>4d-2*ncp.and.dcp2>3d-2.and.printlvl>2.and.ncp>1d0.and.(dispgeoms(j)%energy(k,k)<energyT(1).or.printlvl>3))&
                        print "(4x,A,I5,A,2I2,A,F9.2,A,E12.4,A,F9.2,A)",    &
                        "Large coupling error at pt",j," bkl",k,l,": ",sqrt(dcp/ncp)*100,"% out of ", sqrt(ncp),", ",dcp2*100,"% cp*dE"
             nrmdcp = nrmdcp+dcp
             nrmcp  = nrmcp +ncp
           end if! coupling included
           if(g_exact(j,k,l).and.hasGrad(j,k,l))then
             dgrd    =  dot_product(  &
                   dispgeoms(j)%grads(:nvpt,k,l)-fitG(j,:nvpt,k,l), &
                 ( dispgeoms(j)%grads(:nvpt,k,l)-fitG(j,:nvpt,k,l))*dispgeoms(j)%scale(1:nvpt)  )
             dG_exact=dG_exact+dgrd
             NEx_grd=NEx_grd+1
           end if
           if(e_exact(j,k,l).and.hasEner(j,k).and.hasEner(j,l))then
             dE_exact= dE_exact + (dispgeoms(j)%energy(k,l)-fitE(j,k,l))**2
             NEx_e = NEx_e+1
           end if!e_exact(j,k,l)
         end do!l
         if(dispgeoms(j)%energy(k,k)<energyT(1))then
           if(hasEner(j,k))then
             nrmener = nrmener +    (dispgeoms(j)%energy(k,k)-fitE(j,k,k))**2
             avgener = avgener + abs(dispgeoms(j)%energy(k,k)-fitE(j,k,k))
             inc_e=inc_e+1
           end if
           if(hasGrad(j,k,k))then
             dgrd    =  dot_product(  &
                        dispgeoms(j)%grads(:nvpt,k,k)-fitG(j,:nvpt,k,k), &
                        ( dispgeoms(j)%grads(:nvpt,k,k)-fitG(j,:nvpt,k,k) )*dispgeoms(j)%scale(1:nvpt)  )
             gnrm = dot_product(dispgeoms(j)%grads(:nvpt,k,k),dispgeoms(j)%grads(:nvpt,k,k)*dispgeoms(j)%scale(1:nvpt))
             dgrd = dgrd / gnrm 
             if(((dgrd>1D-2.and.gnrm>1d-4).or.(gnrm*dgrd>1D-4.and.gnrm<=1D-4)).and.printlvl>2.and.(dispgeoms(j)%energy(k,k)<energyT(1).or.printlvl>3)) &
                        print "(4x,A,I5,A,I2,A,F10.3,A,E12.4)", &
                        "Large gradient error at pt",j," state",k," : ",sqrt(dgrd)*100,"% out of ", sqrt(gnrm)
             if(gnrm>1D-4) then
                nrmgrad = nrmgrad + dgrd
                avggrad = avggrad + sqrt(dgrd)
                incdata = incdata + 1
             else
                if(printlvl>4)print "(5x,A,I5,A,E14.4,A,E14.4)",&
                       "Small gradients excluded at point ",j,".  Norm of grad =",sqrt(gnrm),",Norm of error = ",sqrt(dgrd*gnrm)
             end if ! (gnrm>1D-4)
           end if!hasGrad
         end if!dispgeoms(j)%energy(k)<energyT(1)
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
 END SUBROUTINE getError

  !----------------------------------------------------------------------------------------
  !Calculate numerical gradients with respect to the fitting coefficients
  SUBROUTINE getCGrad(cilm,dCi,dLambda,lag,jaco)
    USE progdata, ONLY : printlvl
    USE hddata,   ONLY : RowGrp, ColGrp, offs, GrpLen, nstates,nBasis,nBlks
    IMPLICIT NONE
    DOUBLE PRECISION,dimension(nex+ncons), INTENT(IN)               :: cilm
    DOUBLE PRECISION,dimension(ncons), INTENT(OUT)                  :: dCi
    DOUBLE PRECISION,dimension(nex), INTENT(OUT)                    :: dLambda
    DOUBLE PRECISION,                INTENT(OUT)                    :: Lag
    DOUBLE PRECISION,dimension(nex,ncons), INTENT(OUT)              :: jaco
    
!  These intermediate quantities are constructed for each coefficient i as vectors over
!  all the points, because they are used by multiple equations that involve the same point.
!  W^IJ_i =<f^I|W_i|f^J>
!  D^IJ_i = W^IJ_i / (E^J-E^I)
!  V^IJ_i =<f^I|del W_i|f^J>
!  U^IJ_i = V^IJ_i + sum_K[D^KI_i*h^KJ+D^KJ_i*h^KI]
    DOUBLE PRECISION,dimension(:,:,:),allocatable   ::  WIJ,  DIJ
    DOUBLE PRECISION,dimension(:,:,:,:),allocatable ::  VIJ

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
    DOUBLE PRECISION :: diffsq
!  used to store g and h vectors for degenerate case
    DOUBLE PRECISION :: gvec(nvibs),hvec(nvibs)
!  external BLAS subroutine
    double precision,external :: ddot
    
    if(printlvl>0) print *,"   Entering getCGrad.  Exact coefficient derivatives will be evaluated."
    call system_clock(COUNT=count0,COUNT_RATE=count_rate)

! Initialization of parameters
    allocate(WIJ(npoints,nstates,nstates))
    allocate(VIJ(npoints,nvibs,nstates,nstates))
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
           PmQ(ieq) = dispgeoms(pt)%energy(s1,s1) - fitE(pt,s1,s1)
         else
           PmQ(ieq) = dispgeoms(pt)%energy(s1,s2) - fitE(pt,s1,s2)
         end if
       else  !(s2>0)    it is adiabatic energy difference
         PmQ(ieq) = dispgeoms(pt)%energy(s1,s1) - dispgeoms(pt)%energy(-s2,-s2) - &
                            ( fitE(pt,s1,s1) - fitE(pt,-s2,-s2) )
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
              !call DAXPY(npoints,cc,WMat(iBlk)%List(ind1,iBasis),nvibs+1,WIJ(1,I,J),1)
               WIJ(:,I,J) = WIJ(:,I,J)+cc*WMat(iBlk)%List(ind1:ind1+(npoints-1)*(nvibs+1):nvibs+1,iBasis)
   ! construct the Hamiltonian contribution of VIJ 
              do igd=1,nvibs
                ind1 = ind1+1
   !             call DAXPY(npoints,cc,WMat(iBlk)%List(ind1,iBasis),nvibs+1,VIJ(1,igd,I,J),1)
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
              if (ipt==enfDiab) cycle           ! point has forced eigenvectors
              if(dispgeoms(ipt)%ndeggrp>0)cycle ! point has degeneracy
              DIJ(ipt,I,J)=WIJ(ipt,I,J)/(fitE(ipt,J,J)-fitE(ipt,I,I))
            end do !ipt=1,npoints

   ! fill the lower triangle of D and W
            WIJ(:,J,I)     = WIJ(:,I,J)
            DIJ(:,J,I)     = -DIJ(:,I,J)
            VIJ(:,:,J,I)   = VIJ(:,:,I,J)
          end if !(J.ne.I)
          call system_clock(COUNT=count2,COUNT_RATE=count_rate)
          t3 = t3 + dble(count2-count1)/count_rate
        end do!J
      end do!I
  
      call system_clock(COUNT=count1,COUNT_RATE=count_rate)
   ! construct DIJ for points that have degeneracies 
      do ipt=1,npoints
        if(ipt==enfDiab.or.dispgeoms(ipt)%ndeggrp==0)cycle
        call getDegDij(ipt,WIJ(ipt,:,:),VIJ(ipt,:,:,:),DIJ(ipt,:,:))
      end do! ipt
   ! update the wave function dependent part of Wij and Vij
      do I=1,nstates
        do J=I,nstates
          do k=1,nstates
            if(k.ne.i) WIJ(:,I,J) = WIJ(:,I,J)+DIJ(:,K,I)*fitE(:,K,J)
            if(k.ne.j) WIJ(:,I,J) = WIJ(:,I,J)+DIJ(:,K,J)*fitE(:,K,I)
          end do
          if(I.ne.J) WIJ(:,J,I)=WIJ(:,I,J)
   ! calculate the second piece of VIJ with DIJ and hIJ(=fitG)
          do igd = 1,nvibs
            do K = 1,nstates
              if(k.ne.i)VIJ(:,igd,I,J) = VIJ(:,igd,I,J) + DIJ(:,K,I)*fitG(:,igd,K,J)
              if(k.ne.j)VIJ(:,igd,I,J) = VIJ(:,igd,I,J) + DIJ(:,K,J)*fitG(:,igd,K,I)
            end do ! K = 1,nstates
            if(I.ne.J)  VIJ(:,igd,J,I) = VIJ(:,igd,I,J)  !complete lower triangle
          end do!igd=1,nvibs 
        end do !J=I,nstates
      end do !I=1,nstates
      call system_clock(COUNT=count2,COUNT_RATE=count_rate)
      t3 = t3 + dble(count2-count1)/count_rate
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
             dQ(ieq) = WIJ(pt,s1,s2) 
           end if
         else  !(s2>0)    it is adiabatic energy difference
           dQ(ieq) = WIJ(pt,s1,s1) - WIJ(pt,-s2,-s2)
         end if !(s2>0)
       else !(g==0)   is gradient
         dQ(ieq) = VIJ(pt,g,s1,s2)
       end if!(g==0)
      end do !ieq
      call system_clock(COUNT=count2,COUNT_RATE=count_rate)
      t4 = t4 + dble(count2-count1)/count_rate
      
!  construct derivative of the Lagrangian with respect to coefficient
      dCi(ico) = - dot_product(PmQ(nex+1:)*weight,dQ(nex+1:)*weight)+dot_product(cilm(ncons+1:),dQ(1:nex))+flattening*cilm(ico)
      jaco(1:nex,ico) =  dQ(1:nex)
      call system_clock(COUNT=count1,COUNT_RATE=count_rate)
      t5 = t5 + dble(count1-count2)/count_rate
     end do! iBasis
    end do ! iBlk 

! construct derivative with respect to Lagrange multipliers
    dLambda = - PmQ(1:nex)
    lag =    sum(weight*weight*PmQ(nex+1:)*PmQ(nex+1:))/2-sum(cilm(ncons+1:)*PmQ(1:nex))+&
             sum(cilm(1:ncons)*cilm(1:ncons))*flattening/2
    deallocate(WIJ)
    deallocate(VIJ)
    deallocate(DIJ)

    call system_clock(COUNT=count2)
    if(printlvl>0)then
      print *,"   Gradient/Hessian calculation finished."
      print *,"   Execution time decomposition for getCGrad"
      print *,"   Total   PmQ   c1*c2    W,V,D    dQ   PmQ.dQ  "
      print "(6F8.2)",dble(count2-count0)/count_rate,t1,t2,t3,t4,t5
    end if
  END SUBROUTINE getCGrad
!------------------------------------------------------------------
! Get the Dij at a point of degeneracy by solving linear equations
! The subroutine suppose that the d^I.H.d^J is already correctly stored in fitE
! and the fit gradients in fitG.  
! WIJ and VIJ for the point for the coefficient also needs to be provided.
  SUBROUTINE getDegDij(pt,WIJ,VIJ,DIJ) 
    use hddata, only:  nstates
    IMPLICIT NONE
    INTEGER,intent(IN)          :: pt
    DOUBLE PRECISION,intent(IN) :: WIJ(nstates,nstates),VIJ(nvibs,nstates,nstates)
    DOUBLE PRECISION,intent(OUT):: DIJ(nstates,nstates)
    integer :: nDij
    double precision :: EMat(nstates*(nstates-1)/2,nstates*(nstates-1)/2)
    double precision,dimension(nstates*(nstates-1)/2) :: RVec, sv
    integer :: i,j,k,l,lb,ub
    integer :: grpind(nstates)  ! specifies which degeneracy group a state belongs to
    integer :: smap(nstates,nstates),dmap(nstates*(nstates-1)/2,2)  !mapping between dij and state indices 
    integer :: sgn(nstates,nstates)
    double precision,dimension (nvibs) :: gvec, hvec
    double precision,dimension(:),allocatable :: WORK
    integer :: INFO,LWORK,rank
    integer,allocatable,dimension(:) :: IWORK

    nDij = nstates*(nstates-1)/2
    LWORK = 1000+nDij*nDij*100
    allocate(WORK(LWORK))
    allocate(IWORK(nDij*nDij*5+20*nDij))
    grpind = 0
    ! construct state groupings index table
    do i=1,dispgeoms(pt)%ndeggrp
      lb = dispgeoms(pt)%deg_groups(i,1)
      ub = dispgeoms(pt)%deg_groups(i,2)+lb-1
      grpind(lb:ub) = i
    end do!i=1,dispgeoms(pt)%ndeggrp
    ! construct maps between dij and state indices
    smap = 0
    sgn  = 0
    k    = 0
    do i=1,nstates
      do j=i+1,nstates
        k=k+1
        dmap(k,1) = i
        dmap(k,2) = j
        sgn(i,j)  = 1
        smap(i,j) = k
        sgn(j,i)  = -1
        smap(j,i) = k
      end do!j
    end do!i

    ! construct linear equations for dij
    EMat = 0d0
    do k=1,nDij
      i=dmap(k,1)
      j=dmap(k,2)
      if(grpind(i)==grpind(j).and.grpind(i).ne.0)then
        hvec = fitG(pt,:,I,J)
        gvec = fitG(pt,:,J,J)-fitG(pt,:,I,I)
        rvec(k) = dot_product(VIJ(:,I,J),gvec) + &
                  dot_product(VIJ(:,J,J)-VIJ(:,I,I),hvec)
!        EMat(k,k)=dot_product(gvec,gvec)-4*dot_product(hvec,hvec)
        do l=1,nstates
          !if(l==i.or.l==j)cycle
          if(l/=i)EMat(k,smap(l,i))=EMat(k,smap(l,i))+sgn(l,i)*(  &
                -dot_product(fitG(pt,:,l,j),gvec)+2*dot_product(fitG(pt,:,l,i),hvec))
          if(l/=j)EMat(k,smap(l,j))=EMat(k,smap(l,j))+sgn(l,j)*(  &
                -dot_product(fitG(pt,:,l,i),gvec)-2*dot_product(fitG(pt,:,l,j),hvec))  
        end do!l
      else!i and j are in same degeneracy group
        rvec(k) = - WIJ(i,j) 
        do l=1,nstates
           if(l.ne.i)EMat(k,smap(l,i))=EMat(k,smap(l,i))+fitE(pt,l,j)*sgn(l,i)
           if(l.ne.j)EMat(k,smap(l,j))=EMat(k,smap(l,j))+fitE(pt,l,i)*sgn(l,j)
        end do
      end if!i and j are in same degeneracy group
    end do!k

    !solve the Dij values
    CALL DGELSD(nDij,nDij,1,EMat,nDij,rvec,nDij,sv,1d-12,rank,WORK,LWORK,IWORK,INFO)
    !if(rank<nDij)print "(A,I4,A,I5)","Reduced rank at point ",pt,", rank=",rank
    if(INFO/=0)then
       print *,"Failed to solve linear equations in getDegDij"
       print *,"Point index:",pt
       print *,"ckl at point :"
       do i=1,nstates
         print "(10F10.5)",ckl(pt,:,i)
       end do
       print *,"EMat  :"
       do i=1,nstates*(nstates-1)/2
         print "(20F10.5)",EMat(:,i)
       end do
       print *,"INFO=",info
       stop
    end if
    
    ! pack the solution into Dij matrix
    DIJ = 0D0
    do k=1,nDij
      i=dmap(k,1)
      j=dmap(k,2)
      DIJ(i,j) = rvec(k)
      DIJ(j,i) = - rvec(k)
    end do!k 
    deallocate(WORK)
    deallocate(IWORK)
  END SUBROUTINE getDegDij

!------------------------------------------------------------------
! Try to obtain the best values of Lagrange multipliers
! Since only Lagrange multipliers are changed, there is no need to 
! update Hd or reconstruct error vector.  
! The optimal values is a saddle point rather than minimum of the 
! Lagrangian, (or in the case of unsatisfied constraints, minimum
! of the norm of the gradient of Lagrangian, the subroutine attempts
! to minimize norm of gradient of L.
! On output, the rows of B are orthongalized by the transformation matrix Q
  SUBROUTINE optLag(B,ldb,dCi,asol,Borth)
    USE progdata, ONlY: printlvl
    IMPLICIT NONE
    INTEGER, INTENT(IN)                                         :: ldb
    DOUBLE PRECISION,dimension(ldb,ncons),intent(IN)            :: B
    DOUBLE PRECISION,dimension(ldb,ncons),intent(OUT)           :: Borth
    DOUBLE PRECISION, INTENT(INOUT),dimension(ncons)            :: dCi
    DOUBLE PRECISION,dimension(nex+ncons),INTENT(INOUT)         :: asol

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
    CALL DGEMV('N',nex,ncons,dble(-1),B,ldb,dCi,int(1),dble(0),Bg,int(1))

    ! BBt = B.B^T
    CALL DSYRK('U','N',nex,ncons,dble(1),B,ldb,dble(0),BBt,nex)
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
    asol(ncons+1:)=asol(ncons+1:)+dLag
    do i=1,nex
      if(abs(asol(i+ncons))>1D0.and.printlvl>1)then
       print "(4x,A,I4,A,F8.2,A,4I5)",&
                 "Large Lagrange multiplier. eq",i,",val=",asol(i+ncons),",eqmap=",eqmap(i,:)
      end if
    end do!i=1,nex
    ! update new gradient
    CALL DGEMV('T',nex,ncons,dble(1),B,ldb,dLag,int(1),dble(1),dCi,int(1))
    if(printlvl>1)then
      call system_clock(COUNT=c2)
      print "(A,F6.2,A)","  Lagrange multipliers and Lagrangian gradients updated in ",(c2-c1)/dble(crate),"s"
    end if
    ! B'=Q^T.B
    CALL DGEMM('T','N',nex,ncons,nex,dble(1),BBt,nex,B,ldb,dble(0),Borth,ldb)
    do i=1,nex
      if(abs(w(i))>exactTol) Borth(i,:)=Borth(i,:)/sqrt(w(i))
    end do!i=1,nex
  END SUBROUTINE optLag

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
!---------------------------------------------
! store ckl to a file
  SUBROUTINE writeCkl(cklfl)
    use hddata, only: nstates
    IMPLICIT NONE
    CHARACTER(255),INTENT(IN) :: cklfl

    integer  ios, i,j
    open(unit=7421,file=trim(adjustl(cklfl)),access='sequential',form='formatted',&
       status='replace',action='write',position='rewind',iostat=ios)
    if(ios/=0)then
        print *,"writeCkl: failed to open ckl file ",trim(adjustl(cklfl))
        return
    end if
  
    do i=1,npoints
        write(7421,"(I7)")I
        do j=1,nstates
           write(7421,"(81E32.24)"),ckl(i,:,J)
        end do!j
    end do!i=1,npoints
    print "(2A)"," wave functions exported to file ",trim(adjustl(cklfl))
    close(unit=7421)
  END SUBROUTINE writeCkl
!---------------------------------------------
  SUBROUTINE readCkl(cklfl)
    use hddata, only: nstates
    use progdata, only: printlvl
    IMPLICIT NONE
    CHARACTER(255),INTENT(IN) :: cklfl

    double precision ::  cklbuffer(nstates,nstates),ovlp(nstates)
    integer  ios, i, k, pt, ccount

    if(cklfl=='')return
    print "(A)","  Loading eigenvectors from file "//trim(adjustl(cklfl))
    open(unit=7421,file=trim(adjustl(cklfl)),access='sequential',form='formatted',&
        status='old',action='read',position='rewind',iostat=ios)
    if(ios/=0)then
        print *,"readCkl: failed to open ckl file ",trim(adjustl(cklfl))
        return
    end if

    ios = 0
    ccount=0
    do while (ios==0)
        read(7421,*,iostat=ios) pt,cklbuffer
        if(ios==0.and.pt<=npoints.and.pt>0)then
          do k=1,nstates
            ovlp(k) = abs(dot_product(cklbuffer(:,k),ckl(pt,:,k)))
          end do!k
          if(minval(ovlp)<8d-1.and.ckl(pt,1,1)>-1d2)then
            print "(A,I4)","Small overlap between loaded and automatic Ckl at pt ",pt
            print *,"loaded:"
            do k=1,nstates
              print "(10F10.5)",cklbuffer(:,k)
            end do
            print *,"automatic:"
            do k=1,nstates
              print "(10F10.5)",ckl(pt,:,k)
            end do
          end if!
          ckl(pt,:,:)=cklbuffer
          ccount=ccount+1
        end if
    end do
    if(printlvl>1)  print "(3x,I7,A)",ccount," eigenvectors loaded."
    close(unit=7421)
  END SUBROUTINE readCkl
!---------------------------------------------
! Generate the basis for the fit using linear combinations of primitive basis functions
!   * Values and gradients of primitive functions are evaluated at each data point and put into inter-
!     mediate structures basisVal and basisGrad.  They will later be released.
!   * Basis are generated separately for each symmetry unique block.  The constructed basis are linked
!     to non-unique blocks with pointers.
!   * Linear combinations that define fitting basis is chosen to be the eigenvectors of the dot-product
!     matrix of primitive functions in the space spanned by all values and gradient data.   
!   * In a one state fitting procedure with equal weights assigned, this generates a orthonormal basis 
!     that allows determination of the fit directly from gradients.   In a multistate scenario with non-
!     equal weights, such basis is much more well conditioned than the primitive functions themselves
!   * Eigenvalues that correspond to low occupation numbers are removed.  These components do not contri-
!     bute to any data points at any diabatrization schemes.  There will also be null space that arise 
!     from each individual d values (wave functions), which will be removed from fitting procedure.
!   
  SUBROUTINE genBasis()
  use hddata, only: T3DDList,nl,nr,nblks,nBasis,EvalRawTerms,EvaluateBasis2,deallocDVal,EvaluateVal
  use progdata, only: printlvl
  use cnpi,only: blockSymLs,blockSymId, NBlockSym
  IMPLICIT NONE
  integer  :: i,j,k,l,count1,count2,count_rate,pv1,n1,n2, broot
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
!**************************************************
!*   First construct symmetry unique blocks       *
!**************************************************
  if(printlvl>0) print *,"Generating basis for symmetry unique blocks."
  do l=1,NBlockSym
    k = blockSymLs(l)
    if(printlvl>0) print "(/,A,I2,A,I5,A)"," Constructing intermediate basis for block ",K," with ",npb(k)," matrices"
    ll = nl(k)
    rr = nr(k)
    if(npb(k)==0)then
        print *,"WARNING: No basis function is present for block ",K,".  Skipping basis construction."
        nBas(k)=0
        allocate(ZBas(k)%List(1,1))     ! just a place holder
        allocate(WMat(k)%List(1,1))
        gradNorm(:,:,k) = 0d0
        cycle
    end if
    call system_clock(COUNT=count2)
    if(printlvl>1)write(*,"(A)",advance='no'),"   Evaluating primitive basis... "
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
     do while(dispgeoms(i)%energy(1,1)>energyT(j+1))
       j=j+1
       if(j==10)exit
     end do
     wt = 1d0
     if(ptWeights(i)<1d-5) wt = 0d0
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
    ! ZBas is transformation from old to new basis 
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
        do i=npb(k),npb(k)-nb+1,-1
          ZBas(k)%List(:,npb(k)-i+1) =evec(:,i)
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
        Zbas(k)%List=0d0
        do i=1,npb(k)
            ZBas(k)%List(i,i) =  1D0
        end do
        pv1 = npoints*(1+nvibs)*ll*rr
        ! generate new coefficient matrix in the basis.  WMat = pbas
        allocate(WMat(k)%List(pv1,npb(k)))
        WMat(k)%List=pbas
        deallocate(pbasw)
        deallocate(pbas)
    end if!(TBas>0)

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
  end do! l=1,NBlockSym
!**************************************************
!*   Clone pointers to basis to nonunique blocks  *
!**************************************************
  print *,"Linking basis storage to non-unique blocks"
  do l=1,nblks
    broot = blockSymLs(blockSymId(l))
    if(broot.eq.l)cycle ! it is a unique block.  skipping linking
    ZBas(l)%List =>  ZBas(broot)%List
    WMat(l)%List =>  WMat(broot)%List
    npb(l)       =   npb(broot)
    nbas(l)      =   nbas(broot)
    gradNorm(:,:,l) = gradNorm(:,:,broot)
  end do!l
  1001 format(a,f7.2,a)
  END SUBROUTINE genBasis!
!---------------------------------------------
! initialize variables for makesurf
! Allocate memory for shared structures and construct coefficient and
! equations maps.
  SUBROUTINE initMakesurf()
    use hddata, only: nblks,nstates,order,nBasis,makeCoefMap
    use progdata, only: printlvl
    IMPLICIT NONE
    INTEGER,DIMENSION(:,:),ALLOCATABLE             :: tmpEqMap,tmpExactEq
    DOUBLE PRECISION,dimension(:),allocatable      :: tmpW
    integer    :: i,j,maxeqs

    call getPtList()
    if(allocated(gradNorm))deallocate(gradNorm)
    allocate(gradNorm(npoints,nvibs,nblks))
    if(allocated(ckl))deallocate(ckl)
    allocate(ckl(npoints,nstates,nstates))
    if(allocated(cklguide))deallocate(cklguide)
    allocate(cklguide(npoints,nstates,nstates))
    if(allocated(fitE))deallocate(fitE)
    if(allocated(fitG))deallocate(fitG)
    allocate(fitE(npoints,nstates,nstates))
    allocate(fitG(npoints,nvibs,nstates,nstates))


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

    print *,"Generating fitting basis"
    call genBasis()

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
    call makeEqMap(maxeqs,tmpEqMap,tmpExactEq,tmpW)
    allocate(EqMap(nex+neqs,4))
    if(allocated(weight))deallocate(weight)
    allocate(weight(neqs))
    weight=tmpW(1:neqs)
    eqMap(1:nex,:)=tmpExactEq(1:nex,:)
    eqmap(nex+1:nex+neqs,:)=tmpEqMap(1:neqs,:)
    deallocate(tmpW)
    deallocate(tmpEqMap)
    deallocate(tmpExactEq)
    if(printlvl>0)print*,"    Number of Least Squares Equations:",neqs
    if(printlvl>0)print*,"    Number of Exact Equations:        ",nex
    ! allocate spaces for arrays and matrices
    if((neqs+nex)<ncons.and.printlvl>0)print *,"    Warning: neqs<ncons, data set might be inadequate."
  END SUBROUTINE initMakesurf
!----------------------------------------------------------------------
! solve the quasi-Newton equations with a diagonal Hessian approximation
! Using such approximation avoids the construction and solution of Normal
! equations, which can be very CPU and memory intensive
!
! ARGUMENTS
! ---------
! g     [in] DOUBLE PRECISION,dimension(nex+ncons)
!       Gradients of the Lagrangian
!       First ncons elements are the coefficient gradients and the rest
!       are gradients with respect to Lagrange multipliers
! B     [in] DOUBLE PRECISION,dimension(nex,ncons)
!       Exact block of the Lagrangian Hessian
! x     [out] DOUBLE PRECISION,dimension(nex+ncons)
!       Calculated step vector
! d     [in] DOUBLE PRECISION
!       Diagonal elements of the approximate hessian for the LSE block
!
! METHOD
! ------
! The coefficient-coefficient block of the Lagrangian Hessian is approximated 
! as a diagonal matrix. The exact-coefficient block is included exactly. 
! The Newton-Raphson equations becomes
! ( d*I   B^T  )   ( x )   = - ( g1 )
! ( B     0    )   ( y )       ( g2 )
! Expanding these yields
!   d*x+B^T.y = - g1
!   B x = - g2
! Therefore
! y = (B.B^T)^(-1)(d*g2-B.g1)
! x = -(1/d)( g1+B^T.y)
  SUBROUTINE solve_diagHess(g1,g2,B,x,d)
    use progdata, only : printlvl
    IMPLICIT NONE
    DOUBLE PRECISION,dimension(ncons),    intent(IN)   :: g1
    DOUBLE PRECISION,dimension(nex  ),    intent(IN)   :: g2
    DOUBLE PRECISION,dimension(nex,ncons),intent(IN)   :: B
    DOUBLE PRECISION,dimension(ncons),intent(OUT)      :: x
    DOUBLE PREcision,intent(IN)                        :: d

    double precision,dimension(nex,nex)   :: BBt, evec
    double precision,dimension(nex)       :: rhs1,eval,rhs2,y
    integer   ::  i,nev,LWORK,LIWORK,INFO
    integer   ::  isuppz(2*nex)
    double precision,dimension(:),allocatable :: WORK
    integer,dimension(:),allocatable          :: IWORK
    double precision  ::  tmp(ncons)

    ! no exact equations, then x=-g1/d
    if(nex==0)then
      x=g1/(-d)
      return
    end if

    ! BBt = B.B^t    
    call DSYRK('U','N',nex,ncons,1d0,B,nex,0d0,BBt,nex)

    ! calculate rhs1 = d*g2-B.g1
    CALL DCOPY(nex,rhs1,1,g2,1)
    CALL DGEMV('N',nex,ncons,-1d0,B,nex,g1,1,d*scaleEx,rhs1,1)

    ! solve (B.B^t)y = d*g2-B.g1 =rhs1
    ! DSYEVR is used to perform eigenvalue decomposition of B.Bt
    ! The decomposition is used to invert the matrix and solve the equations.
    ! Decomposition is performed to deal with cases where B.Bt is rank
    ! deficient.  This can occur when exact equations have linear dependencies.
    ! For exampling, putting two equations that cannot be simultaneously
    ! satisfied, or having energy different and energy fitted exactly at the
    ! same time.
    !-------------------------
    ! perform work space query
    CALL DSYEVR('V','A','U',nex,BBt,nex,0.,0.,0,0,exacttol/10,nev,eval,evec,nex,&
                ISUPPZ,tmp,-1,LIWORK,-1,INFO)
    if(INFO/=0)stop"solve_diagHess: work space query failed"
    !-------------------------
    ! allocate workspace
    LWORK = int(tmp(1))
    allocate(WORK(LWORK),stat=info)
    if(info/=0)stop"solve_diagHess: failed to allocate work space."
    allocate(IWORK(LIWORK),stat=info)
    if(info/=0)stop"solve_diagHess: failed to allocate iwork space."
    !-------------------------
    ! do the eigenvalue decomposition  BBt=Q.eval.Q^T
    CALL DSYEVR('V','A','U',nex,BBt,nex,0.,0.,0,0,exacttol/1d8,nev,eval,evec,nex,&
                ISUPPZ,WORK,LWORK,IWORK,LIWORK,INFO)
    if(info/=0)stop"solve_diagHess: failed to perform eigenvalue decomposition"
    if(nex/=nev)print*,"WARNING: solve_diagHess: number of eigenvalues and equations do not match"
    !-------------------------
    ! release workspace 
    deallocate(WORK)
    deallocate(IWORK)
    ! perform inversion rhs2 = Q^T.rhs1
    CALL DGEMV('T',nex,nex,1d0,evec,nex,rhs1,1,0d0,rhs2,1)
    ! rhs2=eval^-1.rhs2
    nev = 0
    do i=1,nex
      if(abs(eval(i))>exacttol)then
        rhs2(i) = rhs2(i)/eval(i)
        nev = nev+1
      else
        rhs2(i) = 0d0
      end if
    end do!i=1,nex
    if(printlvl>1.and.nev<nex)print "(4x,A,I5)",&
                "Exact equations are rank deficient:nNull=",nex-nev 
    ! y = eval.rhs2
    CALL DGEMV('N',nex,nex,1d0,evec,nex,rhs2,1,0d0,y,1)
    ! x =-(1/d)( g1+B^T.y)
    x = g1
    CALL DGEMV('T',nex,ncons,-1/d,B,nex,y,1,-1/d,x,1)
!PRINT *,"VERIFYING SUBROUTINE"
!   B x = - g2
!rhs1 = g2*scaleEx
!CALL DGEMV('N',nex,ncons,1d0,B,nex,x,1,1d0,rhs1,1)
!PRINT *,"RMSE IN g2:",SQRT(DOT_PRODUCT(RHS1,RHS1)/NEX)
!   d*x+B^T.y = - g1
!tmp = g1
!CALL DAXPY(NCONS,d,x,1,tmp,1)
!CALL DGEMV('T',nex,ncons,1d0,B,nex,y,1,1d0,tmp,1)
!PRINT *,"RMSE IN g1:",SQRT(DOT_PRODUCT(TMP,TMP)/NCONS)
  END SUBROUTINE solve_diagHess
END MODULE
!---------------------------------------------
! Transform Hd with a block-by-block transformation ZBas
! This is used by the null-space to do both forward and backward transformation.
! Arguments:
! TRANS     [in] CHARACTER*1
!           Specifies whether the forward or backward transformation
!           will be taken.  Case insensitive.
!           'F' or 'N'      Forward transformation
!           'B' or 'T'      Backward transformation
! hvec_old  [in] DOUBLE PRECISION(*)
!           Original Hd vector before transformation
! hvec_new  [in] DOUBLE PRECISION(*)
!           New Hd vector after transformation
SUBROUTINE tranHd(TRANS,hvec_old,hvec_new)
    use makesurfdata, only: ncons,nex,ncon_total,ZBas,nBas,coefMap
    use hddata,only: nblks
    IMPLICIT NONE
    CHARACTER,INTENT(IN)        :: TRANS
    DOUBLE PRECISION,INTENT(IN) :: hvec_old(*)
    DOUBLE PRECISION,INTENT(OUT):: hvec_new(*)

    integer ::  i,count1,count2,k

    select case (TRANS)
        case('N','n','F','f')
          hvec_new(1:ncons+nex)  = 0d0
          count1 = 0
          do k=1,nblks
            count2 = 0
            do i=1,ncon_total
              if(coefMap(i,2)==k)then
                count2 = count2+1
                hvec_new(count1+1:count1+nBas(k))=hvec_new(count1+1:count1+nBas(k))+ZBas(k)%List(count2,:)*hvec_old(i)
              end if
            end do
            count1=count1+nBas(k)
          end do!k

        case('T','t','B','b')
          count1 = 0
          do k=1,nblks
            count2 = 0
            do i=1,ncon_total
              if(coefMap(i,2)==k)then
                count2 = count2+1
                hvec_new(i) = dot_product(ZBas(k)%List(count2,:),hvec_old(count1+1:count1+nBas(k)))
              end if
            end do
            count1=count1+nBas(k)
          end do!k
        case default
    end select !case (TRANS)

END SUBROUTINE tranHd
!---------------------------------------------
! Fit Hd according to Ab initio Data
SUBROUTINE makesurf()
  use hddata, only: nstates,ncoord,getHdvec,nl,nr,writeHd,nBasis,updateHd,&
                    ExtractHd,getFLUnit,EvaluateHd3
  use progdata, only: printlvl,OUTFILE,AU2CM1,natoms,PI
  use makesurfdata
  use rdleclse
  use DIIS
  IMPLICIT NONE
  INTEGER                                        :: i,j,k,l,count1,count2,count_rate
  INTEGER                                        :: m
  INTEGER                                        :: iter,uerrfl
  DOUBLE PRECISION,dimension(:),allocatable      :: bvec,asol,asol1,dsol,dCi,dLambda,hvec
  DOUBLE PRECISION,dimension(:,:),allocatable    :: jaco,jaco2
  DOUBLE PRECISION,dimension(npoints,nvibs*2)    :: gradtable
  DOUBLE PRECISION,dimension(npoints,2*nstates)  :: enertable
  DOUBLE PRECISION                               :: adif,nrmener,avgener,lag,gmin, dmin, dinc, disp
  double precision                               :: LSErr,ExErr
  DOUBLE PRECISION                               :: nrmgrad,avggrad,nrmD
  CHARACTER(4)                                   :: c1,c2
  CHARACTER(16),dimension(npoints)               :: rlabs
  CHARACTER(16),dimension(nvibs*2)               :: clabs
  double precision,dimension(npoints,nstates,nstates)   :: errGrad,errGradh
  double precision,dimension(ncoord)             :: errgrd
  double precision,dimension(ncoord,3*natoms)    :: binv
  CHARACTER(255)                                 :: fmt, fn
  double precision            ::  dener(nstates,nstates)
  integer                     ::  ios, status, ndep
  double precision, external  ::  dnrm2
  logical                     ::  diff  !whether differential convergence will be used
  DOUBLE PRECISION,dimension(nstates,nstates)            :: hmatPT
  DOUBLE PRECISION,dimension(nvibs,nstates,nstates)      :: dhmatPT
  
  integer,dimension(0:maxiter,npoints)            ::   sg     ! determinant of eigenvectors at each point
  double precision  :: NaN
  character(3)                                    ::  str1, str2  ! # of ener and grad data




  if(printlvl>0)print *,"Entering makesurf"
  NaN  = 0
  NaN  = NaN/NaN
  call initMakesurf()
  call printSurfHeader(ncon_total,ncons,neqs,nex)
  print 1001,"    Memory required to store coefficient matrix:",(nvibs+1)*nstates*(nstates+1)*ncons*7.62939453125D-6/2," MB"
  allocate(hvec(ncon_total))               ! Hd coefficients vector in original form
  allocate(bvec(neqs+nex))                 !rhs for exact and LSE equations
  allocate(asol(ncons+nex))                ! solution vector
  allocate(asol1(ncons+nex))               ! old solution vector
  allocate(dsol(ncons))                    ! change of solution vector
  allocate(jaco(nex,ncons))                ! second derivatives of Lagrangian for exact-LSE jacobian block
  allocate(jaco2(nex,ncons))
  allocate(dCi(ncons))
  allocate(dLambda(nex))
  if(diagHess<=0)then
    allocate(NEL(ncons+nex,ncons+nex),STAT=status)
    if(status/=0) stop 'makesurf  :  failed to allocate memory for NEL'
  end if
  allocate(rhs(ncons+nex))


  iter = 0
  call getHdvec(hvec,coefmap,ncon_total)
  call tranHd('F',hvec,asol1)   ! convert primitive h expansion to block orthogonal basis expansions

  ! load guide ckls
  ckl(:,1,1)=-1d5
  call readCkl(guide)
  cklguide = ckl
  ckl = 0d0
  !----------------------------------
  ! Begin self-consistent iterations
  !----------------------------------
!generate initial eigen vectors
  if(printlvl>0)print 1003
  printlvl=printlvl-1
  call updateEigenVec(asol1)
  printlvl=printlvl+1
  call readCkl(ckl_input) 
  if(printlvl>1)print *,"    Energies from initial Hd and eigenvectors:"
  do i=1,npoints
    CALL EvaluateHd3(asol1,nBas,npoints,i,nvibs,hmatPt,dhmatPt,WMat)
    if(i.ne.enfDiab)call OrthGH_Hd(dispgeoms(i),dhmatPt(:dispgeoms(i)%nvibs,:,:),ckl(i,:,:),100,gcutoff,hasGrad(i,:,:))
    do k=1,dispgeoms(i)%nvibs
      fitG(i,k,:,:)=matmul(transpose(ckl(i,:,:)),matmul(dhmatPt(k,:,:),ckl(i,:,:)))
    end do!k=1,dispgeoms(i)%nvibs
    fitE(i,:,:) = matmul(transpose(ckl(i,:,:)),matmul(hmatPt,ckl(i,:,:)))
    if(printlvl>1)then
      print "(A,I5)","POINT #",i
      print *,"Hd Matrix"
      do k=1,nstates
        print "(10(x,F15.4))",hmatPt(:,k)
      end do
      print *,"Fit Energy Matrix"
      do k=1,nstates
        print "(10(x,F15.4))",fitE(i,:,k)
      end do
    end if
  end do
  asol = asol1
  CALL getCGrad(asol,dCi,dLambda,lag,jaco)
  CALL optLag(jaco,nex,dCi,asol,jaco2)
  print "(3(A,E15.7))","Gradients for coef block: ",dnrm2(ncons,dCi,int(1)),",lag block:",dnrm2(nex,dLambda,int(1)),&
            ", total:",sqrt(dot_product(dCi,dCi)+dot_product(dLambda,dLambda))

  !CALL initDIISg(ndiis,ndstart,ncons,asol,dCi)
  call evaluateError(asol,weight,LSErr,ExErr)
  CALL getError(nrmener,avgener,nrmgrad,avggrad)
  adif = dnrm2(ncons,asol,int(1))
  write(OUTFILE,1005)iter,adif,nrmgrad*100,avggrad*100,nrmener*AU2CM1,avgener*AU2CM1
  adif=toler+1
  diff = .false.
  ! ----------  MAIN LOOP ---------------
  do while(iter<maxiter.and.adif>toler)
   ! Write coefficients of iteration to restart file
   if(trim(restartdir)/='')then !>>>>>>>>>>>>>>>>>>>
     c1=""
     write(c1,"(I4)"),iter
     fn = trim(restartdir)//'/hd.data.'//trim(adjustl(c1))
     if(printlvl>1)print "(A)","  Exporting Hd coefficients to "//fn
     call tranHd('B',asol1,hvec)                ! convert hd to nascent basis
     call updateHd(hvec,coefmap,ncon_total)     ! convert vector hd to list form
     call writeHd(fn,flheader,.false.)          ! save h vector to file
     ! save ckl
     fn =  trim(restartdir)//'/ckl.'//trim(adjustl(c1))
     call writeCkl(fn)
   end if!restartdir/=''   <<<<<<<<<<<<<<<<<<<<<<<
   iter = iter + 1
   if(printlvl>0)print 1000,iter
   
   ! check if should start differential convergence
   if(iter == dfstart)then
     print *,"  Starting differential convergence..."
     diff = .true.
     write(OUTFILE,*)"  Starting differential convergence..."
   end if

   ! get errors 
   if(printlvl>0)then
     printlvl=printlvl-10
     call makebvec(bvec,.true.)
     print "(5x,A)","RMS Errors of Fitting Equations"
     print "(5x,3(A,E12.5))","Exact equations: ",dnrm2(nex,bvec,1)/sqrt(dble(nex)),&
                             ", LSE (unweighted):",dnrm2(neqs,bvec(nex+1),1)/sqrt(dble(neqs)), &
                             ", LSE (weighted): ",dnrm2(neqs,bvec(nex+1:)*weight,1)/dnrm2(neqs,weight,1)
     printlvl=printlvl+10
   end if

   if(diagHess>0)then
     call solve_diagHess(dCi,dLambda,jaco2,dsol,diagHess)
   else !diagHess>0
   ! construct normal equations matrix and right hand side vector
   !-------------------------------------------------------------
     call system_clock(COUNT=count1)
     if(diagHess<=0)call makeNormalEquations(NEL,rhs,weight,ndep)
     call system_clock(COUNT=count2,COUNT_RATE=count_rate)
     if(printlvl>0)print 1004,dble(count2-count1)/count_rate  

   ! solve the equations to get estimated change in coefficients 
   !-------------------------------------------------------------
     call system_clock(COUNT=count1)
     if(diff)then
       rhs(1:nex)=-dLambda*scaleEx
       rhs(nex+1:)=-dCi
       ! solve normal equations for change in coefficients
       CALL solve(ncons,neqs,nex,ndep,NEL,rhs,exacttol,lsetol,dsol,printlvl)
     else
       dsol = asol(1:ncons)
       ! solve normal equations
       CALL solve(ncons,neqs,nex,ndep,NEL,rhs,exacttol,lsetol,asol,printlvl)
       dsol = asol(1:ncons)-dsol
     end if
     CALL cleanArrays()
     call system_clock(COUNT=count2,COUNT_RATE=count_rate)
     if(printlvl>0)print 1002,dble(count2-count1)/count_rate  
   end if!diagHess>0

   ! Perform coefficient stepping or line search to get new Hd
   !-------------------------------------------------------------
   nrmD = dnrm2(ncons,dsol,int(1))
   if(printlvl>0)print "(A,E15.7)","   Original size of displacement : ",nrmD
   ! scaling the whole displacement vector is larger than the maximum size
   if(nrmD>maxd)then
       if(printlvl>0)print "(2(A,E12.5))","  Displacement vector scaled from",nrmD," to ",maxD
       dsol=dsol*maxD/nrmD
       nrmD = maxD
   end if 
   if(linSteps<=0)then
      asol(1:ncons)=asol1(1:ncons)+dsol
      call evaluateError(asol,weight,LSErr,ExErr)
      CALL getError(nrmener,avgener,nrmgrad,avggrad)
      CALL getCGrad(asol,dCi,dLambda,lag,jaco)
      CALL optLag(jaco,nex,dCi,asol,jaco2)
   else !linsteps<=0
     asol(1:ncons)=asol1(1:ncons)
     call takeLinStep(asol,dsol,linSteps,dCi,dLambda,jaco2,dconv)
     CALL getError(nrmener,avgener,nrmgrad,avggrad)
   end if
   ! Print out error analysis
   !-------------------------------------------------------------
   print "(3(A,E15.7))","Gradients for coef block: ",dnrm2(ncons,dCi,int(1)),",lag block:",dnrm2(nex,dLambda,int(1)),&
            ", total:",sqrt(dot_product(dCi,dCi)+dot_product(dLambda,dLambda))
   adif=DNRM2(ncons,asol-asol1,int(1))
   asol1=asol
   ! write iteration information to output file
   write(OUTFILE,1005)iter,adif,nrmgrad*100,avggrad*100,nrmener*AU2CM1,avgener*AU2CM1
   print *,"   Norm of coefficients:  ",dnrm2(ncons,asol,int(1))
   !-------------------------------------------------------------
   ! Write final eigenvectors to file  
   !-------------------------------------------------------------
  enddo !while(iter<maxiter.and.adif>toler)
  !----------------------------------
  ! End of self-consistent iterations
  !----------------------------------
  !print final errors
  call makebvec(bvec,.true.)
  print "(5x,A)","Final RMS Errors of Fitting Equations"
  print "(5x,3(A,E12.5))","Exact equations: ",dnrm2(nex,bvec,1)/sqrt(dble(nex)),&
                             ", LSE (unweighted):",dnrm2(neqs,bvec(nex+1),1)/sqrt(dble(neqs)), &
                             ", LSE (weighted): ",dnrm2(neqs,bvec(nex+1:)*weight,1)/dnrm2(neqs,weight,1)
  if(adif.le.toler)then
   write(OUTFILE,1006)iter
  else
   write(OUTFILE,1007)
  endif

  !------------------------------------------------------------
  ! Write coefficients of final surface to file
  !------------------------------------------------------------
  if(outputfl/='')then
    call tranHd('B',asol1,hvec)             ! convert hd to nascent basis
    call updateHd(hvec,coefmap,ncon_total)  ! transform vector hd to list form
    if(printlvl>1)print *,"  Exporting Hd coefficients to ",trim(outputfl)
    call writeHd(outputfl,flheader,.false.) ! save h vector to file
  end if
  !------------------------------------------------------------
  ! Write final eigenvectors to file 
  !------------------------------------------------------------
  if(ckl_output/='') call writeCkl(ckl_output)
  !------------------------------------------------------------
  ! Compare ab initio and computed derivative couplings
  !------------------------------------------------------------
  errGrad=dble(0)
  errGradh=dble(0)
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
    if(any(hasGrad(:,i,i)))then
      write(OUTFILE,1016)i
      k=0
      do l = 1,npoints
        if(.not.hasGrad(l,i,i))cycle
        k=k+1
        do m = 1,nvibs
          gradtable(k,2*m-1) = dispgeoms(l)%grads(m,i,i)
          gradtable(k,2*m)   = fitG(l,m,i,i)-dispgeoms(l)%grads(m,i,i)
          if(m<=dispgeoms(l)%nvibs)errGrad(l,i,i)=errGrad(l,i,i)+gradtable(k,2*m)**2!*dispgeoms(l)%scale(m)
        enddo !m=1,dispgeoms(i)%nvibs
        errGrad(l,i,i)=sqrt(errGrad(l,i,i))
        write(c2,'(i4)')l
        rlabs(k) = ' GM '//trim(adjustl(c2))
      enddo!l=1,npoints
      call printMatrix(OUTFILE,rlabs,clabs,int(8 ),npoints,k,2*nvibs,gradtable,int(13),int(8))
    end if!any(hasGrad)
  enddo!i=1,nstates

  do i = 1,nstates
   do j = 1,i-1
    if(any(hasGrad(:,i,j)))then
      write(OUTFILE,1015)j,i
      k=0
      do l = 1,npoints
        if(.not.hasGrad(l,i,j))cycle
        k=k+1
        do m = 1,nvibs
          gradtable(k,2*m-1) = dispgeoms(l)%grads(m,i,j)
          gradtable(k,2*m)   = fitG(l,m,i,j)-dispgeoms(l)%grads(m,i,j)
          if(m<=dispgeoms(l)%nvibs)then
            errGrad(l,i,j)=errGrad(l,i,j)+(fitG(l,m,i,j)-dispgeoms(l)%grads(m,i,j))**2
            errGradh(l,i,j)=errGradh(l,i,j)+&
                 (fitG(l,m,i,j)-dispgeoms(l)%grads(m,i,j))**2
            errGradh(l,j,i)=errGradh(l,j,i)+&
                 (fitG(l,m,i,i)-fitG(l,m,j,j)-&
                 dispgeoms(l)%grads(m,i,i)+dispgeoms(l)%grads(m,j,j))**2
          end if!(m<=dispgeoms(l)%nvibs)
        enddo !m=1,dispgeoms(i)%nvibs
        errGrad(l,i,j)=sqrt(errGrad(l,i,j)) 
        errGrad(l,j,i)=errGrad(l,i,j)
        errGradh(l,i,j)=sqrt(errGradh(l,i,j))
        errGradh(l,j,i)=sqrt(errGradh(l,j,i))/2
        write(c2,'(i4)')l
        rlabs(k) = ' GM '//trim(adjustl(c2))
      enddo!l=1,npoints
      call printMatrix(OUTFILE,rlabs,clabs,int(8 ),npoints,k,2*nvibs,gradtable,int(13),int(8))
    end if!any(hasGrad(i,j))
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
    if(hasEner(j,k))then
        enertable(j,2*k-1) = (dispgeoms(j)%energy(k,k))*AU2CM1
        enertable(j,2*k)   = fitE(j,k,k)*AU2CM1 - enertable(j,2*k-1)
    else
        enertable(j,2*k-1) = NaN
        enertable(j,2*k)   = NaN
    end if
   enddo!k=1,nstates
   write(c2,'(i4)')j
   rlabs(j) = ' GM '//trim(adjustl(c2))
  enddo!j=1,npoints
  call printMatrix(OUTFILE,rlabs,clabs,2*nstates,npoints,npoints,2*nstates,enertable,int(14),int(4))
  print *," Weighed Error Norms for Energy Gradients and Couplings"
  print *," Gradients are given in the following order:"
  print "(10(I4,'-',I4,',',4X),I4,'-',I4)",(((/j,k/),k=1,j),j=1,nstates)
  write(fmt,'("(A,2X,A,2X,",I2,"A)")'),(nstates+1)*nstates
  print trim(fmt)," PT ","  WT  ",(("   ErrG  ",k=1,j),j=1,nstates),&
                           ((" Ab Grd  ",k=1,j),j=1,nstates)
  fmt=""
  write(fmt,'("(I5,2X,F6.3,2X,",I2,"(2(X,E11.4),2x))")'),(nstates+1)*nstates/2
  do i=1,npoints
    do k=1,nstates
      dener(k,k)=1D0
      do l=k+1,nstates
        dener(k,l)= abs(dispgeoms(i)%energy(k,k)-dispgeoms(i)%energy(l,l))
        dener(l,k)= dener(k,l)
      end do
    end do
    print trim(fmt),i,ptWeights(i),(  &
          ([errGrad(i,j,k),dnrm2(dispgeoms(i)%nvibs,dispgeoms(i)%grads(:,j,k),1)]/dener(j,k),k=1,j) &
                               ,j=1,nstates)
  end do

  print *,""
  print "(A)","Fit and and ab initio g and h vectors between each pairs of states and their errors"
  fmt=""
  write(fmt,'("(I5,2X,",I2,"(3(X,E11.4),2X))")'),(nstates-1)*nstates/2
  print *,"h vector: (err,fit,ab)"
  do i=1,npoints
    print trim(fmt),i,( (errGradh(i,j,k),k=1,j-1) ,j=1,nstates),&
       ((dnrm2(dispgeoms(i)%nvibs,dispgeoms(i)%grads(:,j,k),int(1)), &
          k=1,j-1),j=1,nstates), &
       ((dnrm2(dispgeoms(i)%nvibs,fitG(i,:,j,k),int(1)), &
          k=1,j-1),j=1,nstates)
  end do
  print *,"g vector: (err,fit,ab)"
  do i=1,npoints
    print trim(fmt),i,( (errGradh(i,k,j),k=1,j-1) ,j=1,nstates),&
       ((dnrm2(dispgeoms(i)%nvibs,dispgeoms(i)%grads(:,j,j)-dispgeoms(i)%grads(:,k,k),int(1))/2, &
          k=1,j-1),j=1,nstates),  &
       ((dnrm2(dispgeoms(i)%nvibs,fitG(i,:,j,j)-fitG(i,:,k,k),int(1))/2, &
          k=1,j-1),j=1,nstates)
  end do

  if(printError)then
  !------------------------------------------------------------------
  ! output fitting error to error.log and geom info to refgeom
  !------------------------------------------------------------------
    uerrfl = getFLUnit()
    open(unit=uerrfl,file='error.log',access='sequential',form='formatted',&
       status='replace',action='write',position='rewind',iostat=ios)
    if(ios/=0)print *,"FAILED TO CREATE FILE error.log"
    write(c2,'(i4)')ncoord
    do i=1,npoints
      do j=1,nstates
        if(hasEner(i,j))then
          write(uerrfl,"(E15.7)",advance='no')  fitE(i,j,j)-(dispgeoms(i)%energy(j,j))
        else
          write(uerrfl,"(E15.7)",advance='no')  0d0
        end if
      end do!j
      write(uerrfl,*) ""
      do j=1,nstates
        errgrd = dble(0)
        if(hasGrad(i,j,j))then
          call ginv(ncoord,dispgeoms(i)%nvibs,dispgeoms(i)%bmat,ncoord,binv,1D-2)
          do k=1,dispgeoms(i)%nvibs
              if(dispgeoms(i)%scale(k)>1D-1)&
                  errgrd=errgrd+binv(:,k)*(fitG(i,k,j,j)-dispgeoms(i)%grads(k,j,j))
          end do
        end if!hasGrad
        write(uerrfl,"("//trim(c2)//"E15.7)") errgrd
      end do
    end do 
    close(uerrfl)
    open(unit=uerrfl,file='refgeom',access='sequential',form='formatted',&
      status='replace',action='write',position='rewind',iostat=ios)
    if(ios/=0)print *,"FAILED TO CREATE FILE error.log"
    do i=1,npoints
      write(uerrfl,*) ""
      write(uerrfl,"(3F20.15)") dispgeoms(i)%cgeom
    end do
    close(uerrfl)
  !------------------------------------------------------------------
  !  generate fitinfo.csv which contains fitting info for each point
  !  in csv format which can be read in by an analysis program 
  !------------------------------------------------------------------
! Output point destiny table
    open(unit=uerrfl,file='fitinfo.csv',access='sequential',form='formatted',&
       status='replace',action='write',iostat=ios)
    if(ios/=0)then
      print *,"Failed to open file fitinfo.csv for write."
    else
   ! write out header which contains the number of states and energy threshold
      write(unit=uerrfl,fmt='(I8,F15.4)'),nstates,energyT(1)*au2cm1
      ! format :
      ! ID,weight,EnerInc,GrdInc,EnerErr,GrdErr
      do i=1,npoints
        write(unit=uerrfl,fmt='(I8,",",F12.6)',advance='no')i,ptWeights(i)
        write(str1,"(I3)")  nstates+nstates*(nstates+1)/2
        write(unit=uerrfl,fmt='('//str1//'(",",L))',advance='no'),(incener(i,j,j),j=1,nstates),&
                ((incgrad(i,j,k),k=j,nstates),j=1,nstates)
        write(str1,"(I3)")  nstates
        write(unit=uerrfl,fmt='('//str1//'(",",E24.16))',advance='no'),&
                (dispgeoms(i)%energy(j,j)*AU2CM1,j=1,nstates)
        write(unit=uerrfl,fmt='('//str1//'(",",E24.16))',advance='no'),&
                ((fitE(i,j,j)-dispgeoms(i)%energy(j,j))*AU2CM1,j=1,nstates)
        write(str1,"(I3)")  nstates*(nstates+1)/2
        write(unit=uerrfl,fmt='('//str1//'(",",E24.16))',advance='no'),&
                ((dnrm2(nvibs,dispgeoms(i)%grads(:,j,k),int(1)),k=j,nstates),j=1,nstates)
        write(unit=uerrfl,fmt='('//str1//'(",",E24.16))',advance='no'),&
                (( errGrad(i,j,k),k=j,nstates),j=1,nstates)
        write(unit=uerrfl,fmt='(A)') " "
      end do
      close(uerrfl)
    endif
  end if

  if(printlvl>0)print *,"    deallocating arrays"
  !------------------------------------------------------------------
  ! deallocate arrays
  !------------------------------------------------------------------
  CALL cleanDIIS()
  deallocate(asol)
  deallocate(asol1)
  if(diagHess<=0)deallocate(NEL)
  deallocate(rhs)
  deallocate(bvec)
  deallocate(ckl)
  if(printlvl>0)print *,"Exiting makesurf()"
1000 format(/,2X,"  ITERATION ",I3)
1001 format(a,f9.2,a)
1002 format(4X,"Hd coefficients solved after ",F8.2," seconds.")
1003 format(3X,"Generating initial eigenvectors.")
1004 format(4X,"Normal equations constructed after ",F8.2," seconds.")
1005 format(1x,'Iteration',i4,": Delta=",E9.2,", d[g]%=",f7.2,   &
                ", <d[g]%>=",f7.2,"%, d[E]=",f8.2,", <d[E]>=",f8.2)
1006 format(/,2x,'Computation Converged after ',i5,' Iterations')
1007 format(/,2x,'Computation did not Converge')
1014 format(/,2x,'REPRODUCTION OF LEAST-SQUARES DATA ---------')
1015 format(/,2x,'Derivative Coupling, States: [',i3,',',i3,']')
1016 format(/,2x,'Energy Gradients, State [',i3,']')
1017 format(/,2x,'Energies')
END SUBROUTINE makesurf

! Print job description to surfgen.out
SUBROUTINE printSurfHeader(totcons,cons,eqs,nexact)
  use progdata, only: OUTFILE
  use makesurfdata,only: maxiter,toler,gcutoff
  IMPLICIT NONE
  INTEGER,INTENT(IN)                  :: totcons,cons
  INTEGER,INTENT(IN)                  :: eqs,nexact

  write(OUTFILE,1009)
  write(OUTFILE,1001)totcons
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
1001 format(2x,'Number of Nascent Coefficients:                 ',i7)
1002 format(2x,'Number of Coefficients after Nullspace Removal: ',i7)
1003 format(2x,'Number of Least Squares Equations:              ',i7)
1004 format(2x,'Maximum Number of iterations:                   ',i10)
1005 format(2x,'Coefficient convergence tolerance:              ',f12.8)
1006 format(2x,'Threshold for gradient values:                  ',f12.8)
1007 format(2x,'Number of Exact Equations:                      ',i5,/)
end SUBROUTINE printSurfHeader

!-----------------------------------------------------------------------------------
! determined the phases of Wavefunctions that best reproduce ab initio couplings
  SUBROUTINE fixphase(nvibs,scale,fitgrad,abgrad,ckl,phaseList,hasGrad)
    use hddata, only: nstates
    IMPLICIT NONE
    INTEGER,intent(IN)                                        :: nvibs
    DOUBLE PRECISION,dimension(nvibs,nstates,nstates),intent(inout)  :: fitgrad
    double precision,dimension(nvibs,nstates,nstates),intent(in)     :: abgrad
    double precision,dimension(nvibs),intent(in)              :: scale
    LOGICAL,dimension(nstates,nstates),intent(in)             :: hasGrad
    DOUBLE PRECISION,dimension(nstates,nstates),INTENT(INOUT) :: ckl
    integer,dimension(2**(nstates-1),nstates),INTENT(IN)      :: phaseList

    INTEGER                                                   :: j,k,l,m,minID
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
          do m=1,l-1
            if(hasGrad(l,m))  errNorm = errNorm+ diff(l,m)**2 *scale(k)
          end do!m
        end do!l
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
!both gradients and couplings are considered when doing the ordering
!---------------------------------------------------------------------
SUBROUTINE gradOrder(ptid,fitpt,abpt,ckl,pmtList,LDP,w_en,w_grd,w_cp,cpE1,cpE2)
  use combinatorial
  use hddata,only: nstates
  use progdata, only: abpoint,printlvl
  use makesurfdata, only:incgrad,incener,e_exact,g_exact,hasEner,hasGrad,orderall
  IMPLICIT NONE
  type(abpoint),INTENT(IN)                                   :: abpt
  type(abpoint),INTENT(INOUT)                                :: fitpt
  DOUBLE PRECISION,DIMENSION(nstates,nstates),intent(INOUT)  :: ckl
  DOUBLE PRECISION,INTENT(IN)                                :: w_en,w_grd,w_cp,cpE1,cpE2
  INTEGER,INTENT(IN)                                         :: ptid,LDP
  INTEGER,DIMENSION(LDP,nstates),INTENT(IN)                  :: pmtList

  integer :: i,j,ldeg,udeg,ndeg, pmt, best_p, m
  double precision :: min_err, err, err1
  DOUBLE PRECISION,DIMENSION(nstates,nstates)       :: ckl_new
  DOUBLE PRECISION,DIMENSION(abpt%nvibs)            :: diff,diff2
  DOUBLE PRECISION,DIMENSION(abpt%nvibs,nstates,nstates) :: gradnew
  double precision,dimension(nstates,nstates)       :: en_new
  double precision   ::  shift,de,ndiff
  integer :: reord(nstates)
  if(printlvl>3.and.fitpt%ndeggrp>0)print "(A,I4)","Reordering point",ptid
  min_err=0d0
  err1   =0d0
  best_p =0
  do i=1,fitpt%ndeggrp
    ldeg=fitpt%deg_groups(i,1)
    ndeg=fitpt%deg_groups(i,2)
    udeg=ldeg+ndeg-1
    if(printlvl>3)then
      print "(A,I3,A,2I3)","Permutation group",i,", range ",ldeg,udeg
    end if
    if(hasEner(ptid,ldeg).and.hasEner(ptid,udeg))then
        shift = -(fitpt%energy(ldeg,ldeg)+fitpt%energy(udeg,udeg)-abpt%energy(ldeg,ldeg)-abpt%energy(udeg,udeg))/2
    else
        shift = 0d0
    end if
    do pmt=1,factl(ndeg)
      err=dble(0)
      ! establish new ordering and store in array reord
      do j=1,nstates
        reord(j)=j
      end do
      do j=ldeg,udeg
        reord(j)=pmtList(pmt,j-ldeg+1)+ldeg-1
      end do
      if(printlvl>3)print "(A,10I3)","  New order: ",reord
      do j=ldeg,udeg
      ! calculate error contributions from gradients 
        if(hasGrad(ptid,j,j).and.(incgrad(ptid,j,j).or.g_exact(ptid,j,j).or.orderall))then
          diff=fitpt%grads(:abpt%nvibs,reord(j),reord(j))-abpt%grads(:abpt%nvibs,j,j)
          err=err+dot_product(diff,diff)*w_grd**2
        end if
      ! calculate error contributions from energies 
        if(hasEner(ptid,j).and.(incener(ptid,j,j).or.e_exact(ptid,j,j).or.orderall))  &
                    err=err+ (fitpt%energy(reord(j),reord(j))-abpt%energy(j,j))**2*w_en**2
      ! add coupling contrubutions
        do m=1,nstates
           if(m==j)cycle
           if(hasGrad(ptid,m,j).and.(incgrad(ptid,m,j).or.g_exact(ptid,m,j).or.orderall))then
             de = abs(abpt%energy(j,j)-abpt%energy(m,m))
             if(de<cpE2)de=cpE2
             diff=fitpt%grads(:abpt%nvibs,reord(m),reord(j))-abpt%grads(:abpt%nvibs,m,j)
             diff2=fitpt%grads(:abpt%nvibs,reord(m),reord(j))+abpt%grads(:abpt%nvibs,m,j)
             ndiff = min(dot_product(diff,diff),dot_product(diff2,diff2)) 
             if(m>=ldeg.and.m<=udeg)then
               err=err+ndiff*(w_cp*cpE1/de)**2/2
             else
               err=err+ndiff*(w_cp*cpE1/de)**2
             end if
           end if
        end do
      end do!j=ldeg,udeg
      if(printlvl>3)  print "(A,I3,A,E12.5)","pmt",pmt,",err=",sqrt(err)
      if(pmt==1)then
        min_err=err
        best_p=pmt
        err1=err
      else if (err<min_err)then
        min_err=err
        best_p=pmt
      end if
    end do!pmt=1,factl(ndeg)
    if(printlvl>3)  print "(A,I3)","best permutation:",best_p
    ckl_new=ckl
    en_new =fitpt%energy
    gradnew=fitpt%grads(:abpt%nvibs,:,:)
    do j=ldeg,udeg
      m=pmtList(best_p,j-ldeg+1)+ldeg-1
      ckl_new(:,j)   =ckl(:,m)
      en_new(:,j)    = fitpt%energy(:,m)
      gradnew(:,:,j) =fitpt%grads(:abpt%nvibs,:,m)
    end do
    if(printlvl>0.and.best_p/=1)print 1000,ptid,ldeg,&
                    ldeg+ndeg-1,sqrt(err1),sqrt(min_err)
    fitpt%grads(:abpt%nvibs,:,:)=gradnew
    fitpt%energy=en_new
    do j=ldeg,udeg
      m             =pmtList(best_p,j-ldeg+1)+ldeg-1
      gradnew(:,j,:)=fitpt%grads(:abpt%nvibs,m,:)
      en_new(j,:)   =fitpt%energy(m,:)
    end do
    ckl=ckl_new
    fitpt%grads(:abpt%nvibs,:,:)=gradnew
  end do !i=1,fitpt%ndeggrp
1000 format(6X,"Point ",I4," States ",I2," to ",I2,&
          " ordered by gradients. Error:",E11.4,"->",E11.4)
END SUBROUTINE gradOrder

! reorder states by maximizing absolute overlap with guide ckl
SUBROUTINE guideOrder(ptid,realckl,guideckl,pmtList,LDP)
  use combinatorial
  use hddata, only: nstates
  use progdata, only: printlvl
  use makesurfdata, only: ptWeights
  IMPLICIT NONE
  DOUBLE PRECISION,dimension(nstates,nstates),intent(IN)     :: guideckl
  DOUBLE PRECISION,dimension(nstates,nstates),intent(INOUT)  :: realckl
  INTEGER,INTENT(IN)                                         :: LDP
  INTEGER,DIMENSION(LDP,nstates),INTENT(IN)                  :: pmtList
  double precision :: reordered(nstates,nstates), ovlp, maxovlp
  integer  :: ptid,i,j, maxi

  maxovlp = -1d0
  do i=1,factl(nstates)
    ovlp = 0d0
    do j=1,nstates
      ovlp = ovlp+abs(dot_product(realckl(:,pmtList(i,j)),guideckl(:,j)))
    end do
    if(ovlp>maxovlp)then
      maxovlp = ovlp
      maxi = i
    end if
  end do!i
  do j=1,nstates
    reordered(:,j)=realckl(:,pmtList(maxi,j))
  end do !j
  if(printlvl>3.and.maxi>1)then
    PRINT *,"Guide ordering at point ",PTID
    PRINT *,"before:"
    PRINT "(4F10.6)",REALCKL
    PRINT *,"after:"
    PRINT "(4F10.6)",reordered
  else
     if(printlvl>0.and.maxi>1.and.ptWeights(ptid)>1d-8)print 1000,ptid,maxovlp/nstates
  end if
  realckl = reordered
1000 format(6x,"POINT ",I4," reordered by overlap with guide. Overlap=",F8.5)
END SUBROUTINE guideOrder

! read input parameters for makesurf
SUBROUTINE readMakesurf(INPUTFL)
  USE hddata,only: nstates
  USE progdata,only: natoms, AU2CM1
  USE makesurfdata
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: INPUTFL
  integer :: i
  NAMELIST /MAKESURF/ npoints,maxiter,toler,gcutoff,gorder,exactTol,LSETol,outputfl,TBas,ecutoff,egcutoff, guide,&
                      flheader,ndiis,ndstart,enfDiab,followPrev,w_energy,w_grad,w_fij,usefij, deg_cap, eshift, &
                      ediffcutoff,nrmediff,ediffcutoff2,nrmediff2,rmsexcl,useIntGrad,intGradT,intGradS,gScaleMode,  &
                      energyT,highEScale,maxd,scaleEx, ckl_output,ckl_input,dijscale,  diagHess, dconv, printError, &
                      dfstart,linSteps,flattening,searchPath,notefptn,gmfptn,enfptn,grdfptn,cpfptn,restartdir,orderall,&
                      gradcutoff
                      
  npoints   = 0
  gradcutoff=100000
  printError= .false.
  maxiter   = 3
  orderall  = .true.
  diagHess  = -1d0
  dconv     = 1d-4
  followprev= .false.
  usefij    = .true.
  deg_cap   = 1d-5
  ecutoff   = 1d0
  egcutoff  = 6d-1
  TBas      = 1D-6
  eshift    = 0d0
  ckl_input = ''
  ckl_output= 'ckl.out'
  guide = ''
  linSteps  = 0
  dfstart   = 0
  dijscale  = 1d0
  useIntGrad= .true.
  restartdir= ''
  intGradT  = 1D-3
  scaleEx   = 1D0
  intGradS  = 1D-1
  maxd      = 1D0
  energyT   = 1D30
  highEScale = 1D0
  gScaleMode = 0
  toler     = 1D-3
  gorder    = 1D-3
  gcutoff   = 1D-14
  exactTol      = 1D-12
  LSETol        = 1D-7
  flattening    = 1D-8
  outputfl      = ''
  flheader      = '----'
  ndiis         = 10
  ndstart       = 10
  enfDiab       = 0
  w_energy      = dble(1)
  w_grad        = dble(1)
  w_fij         = dble(1)
  nrmediff      = 2.D4
  ediffcutoff   = 20.
  nrmediff2     = 100.
  ediffcutoff2  = 1.
  rmsexcl       = 0
  SearchPath    = ''
  SearchPath(1) ='.'
  notefptn      ='note'
  enfptn        ='energy.all'
  gmfptn        ='geom.all'
  grdfptn       ='cartgrd.drt1.state$.all'
  cpfptn        ='cartgrd.nad.drt1.state$.drt2.state$.all'
  read(unit=INPUTFL,NML=MAKESURF)
  if(useIntGrad)then
    nvibs=3*natoms-6
  else
    nvibs=3*natoms
  end if!useIntGrad

  ! get the number of search paths
  NSearchPaths=MaxSearchPaths
  do i=1,MaxSearchPaths
    SearchPath(i)=adjustl(SearchPath(i))
    if(trim(SearchPath(i))=='')then
        NSearchPaths=i-1
        exit
    end if
  end do
  if(NSearchPaths==0) STOP "Empty search path list.  No data to fit."

  ! round up inappropriate option values and do unit changes
  if(linSteps<0)  linSteps=0
  energyT = energyT / AU2CM1
  gradcutoff = gradcutoff/ AU2CM1
  if(nstates>nstates.or.nstates<1)nstates=nstates
END SUBROUTINE readMakesurf

! load pointwise options
SUBROUTINE readDispOptions(exclEner,exclGrad,exactEner,exactGrad,exactDiff,enfGO,ptWt)
  use hddata, only:nstates
  use progdata, only: PTFL,printlvl
  use makesurfdata,only:TEqList,npoints
  IMPLICIT NONE
  TYPE(TEqList),INTENT(OUT):: exclEner,exclGrad,exactEner,exactGrad,exactDiff,enfGO
  double precision,dimension(npoints),intent(out)  :: ptWt
  character(255)           :: comment
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
    exclEner%Length=0
    exclGrad%Length=0
    exactEner%Length=0
    exactGrad%Length=0
    exactDiff%Length=0
    enfGO%Length=0
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
  END SUBROUTINE fillout
END SUBROUTINE readDispOptions
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
  DOUBLE PRECISION :: NaN, enbuffer(nstates)

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
  call printMatrix(OUTFILE,rlabs,clabs,npr,ncoord,ncoord,npts,disps,int(11),int(6))

  write(OUTFILE,'(/,2x,"Ab Initio Energies")')
  write(str,'(i4)') nstates
  NaN = 0
  NaN = NaN/NaN
  do j=1,npts
    do k=1,nstates
      if(.not. hasEner(j,k))then
        enbuffer(k)=NaN
      else
        enbuffer(k) = dispgeoms(j)%energy(k,k)*AU2CM1
      end if
    end do
    write(OUTFILE,'(6x,i5,'//trim(str)//'F14.2)')j,enbuffer
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
! This function gets the number of state index place holders ('$') in an input
! filename pattern.
INTEGER FUNCTION PatternNInd(pattern)
    IMPLICIT NONE
    CHARACTER(255),intent(IN) :: pattern
    integer :: i, nind

    nind = 0
    do i=1,len_trim(pattern)
        if(pattern(i:i)=='$')nind=nind+1
    end do!i
    PatternNInd=nind
END FUNCTION PatternNInd
!------------------------------------------------------------------------------
!Read in geometry and ab initio data for a set of points which will be used
!in the Hd fitting procedure.
!
SUBROUTINE readdisps()
  use hddata
  use progdata,only:natoms,printlvl
  use makesurfdata
  IMPLICIT NONE
  INTEGER                                      :: i,j,k,l,ios
  DOUBLE PRECISION,dimension(3*natoms,npoints) :: cgrads,cgeoms
  DOUBLE PRECISION,dimension(nstates,npoints)  :: eners
  character(255)                               :: filename,infile
  character(3),dimension(natoms)               :: atoms
  double precision,dimension(natoms)           :: anums,masses
  integer,external  :: PatternNInd
! ptinfile : number of points in a geometry input file.
! npts     : number of points that were already included and have proper data
! nnew     : number of new geometries added
! lb, ub   : energies for states between lb and ub are given in energy file
  integer  :: ptinfile, npts, nnew,lb,ub
  double precision  :: NaN
  double precision,external :: dnrm2

  NaN = 0d0
  NaN = NaN/NaN
  if(printlvl>0)print '(4X,A,I5,A)','Reading',npoints,' displacements'
  if(allocated(dispgeoms))deallocate(dispgeoms)
  allocate(dispgeoms(npoints))
  do j = 1,npoints
     allocate(dispgeoms(j)%cgeom(3*natoms))
     allocate(dispgeoms(j)%igeom(ncoord))
     allocate(dispgeoms(j)%energy(nstates,nstates))
     allocate(dispgeoms(j)%grads(3*natoms,nstates,nstates))
     dispgeoms(j)%grads=NaN
     allocate(dispgeoms(j)%bmat(ncoord,3*natoms))
     allocate(dispgeoms(j)%lmat(3*natoms,3*natoms))
     allocate(dispgeoms(j)%scale(3*natoms))
     allocate(dispgeoms(j)%eval(3*natoms))
  enddo
  allocate(hasEner(npoints,nstates))
  allocate(hasGrad(npoints,nstates,nstates))

  hasEner = .false.
  hasGrad = .false.
  SPNotes = ''
  npts    = 0
  do i=1,NSearchPaths
    if(printlvl>0)print "(A)",'  Searching input file in path <'//trim(SearchPath(i))//'>'
    ! read note for the search path
    if(trim(notefptn(i))/='')then
        infile = trim(SearchPath(i))//'/'//trim(adjustl(notefptn(i)))
        open(unit=7423,file=trim(infile),access='sequential',form='formatted',&
                status='old',action='read',position='rewind',iostat=ios)
        if(ios==0)then
            read(unit=7423,fmt="(A)",IOSTAT=ios)SPNotes(i)
            if(ios==0.and.printlvl>0)print "(4X,A)",trim(SPNotes(i))
        end if
        close(unit=7423)
    end if

    ! read geometries
    if(PatternNInd(gmfptn(i))>0.or.len_trim(gmfptn(i))==0)then
        print *,"invalid geometry file naming pattern.  skipping directory..."
        cycle
    end if
    infile = trim(SearchPath(i))//'/'//trim(adjustl(gmfptn(i)))
    if(printlvl>1)print 1001,'- Reading geometry input <'//trim(infile)//'>'
    ptinfile=npoints
    call readColGeom(infile,ptinfile,natoms,atoms,anums,cgeoms,masses)
    if(printlvl>1)print 1002,"found ",ptinfile, " geometries."
    if(ptinfile+npts>npoints)then
        print *,"WARNING: Amount of data in input file exceeded point count. Data truncated."
        ptinfile=npoints-npts
    end if
    ! we do not yet increase the number of total points <npts> because we still
    ! need the old total as offset for other data inputs
    do j = 1,ptinfile
        dispgeoms(j+npts)%id = j+npts
        dispgeoms(j+npts)%cgeom(1:3*natoms) = cgeoms(1:3*natoms,j)
    enddo!j=1,ptinfile
    nnew = ptinfile

    ! read energy data
    infile = trim(SearchPath(i))//'/'//trim(adjustl(enfptn(i)))
    if(len_trim(enfptn(i))>0.and.PatternNInd(gmfptn(i))==0)then
        if(printlvl>1)print 1001,'- Reading energy input <'//trim(infile)//'>'
        ptinfile=npoints
        call readEner(infile,ptinfile,nstates,eners,lb,ub)  !k and l specify range of data
        if(ptinfile>nnew)then
            print *,"WARNING : energy file contains more entries than geometry file! Data truncated."
            ptinfile=nnew
        elseif (ptinfile<nnew)then
            if(printlvl>1)print *,"WARNING : energy file contains less entries than geometry file! "
        end if
        if(ptinfile==0)then
            print *,"WARNING : Cannot find any energy data in current search path. Skipping..."
            cycle
        end if
        if(printlvl>3.or.(printlvl>1.and.ptinfile<nnew))print 1002,"found ",ptinfile," energy data entries."
        do j = 1,ptinfile
            dispgeoms(j+npts)%energy = NaN
            dispgeoms(j+npts)%energy(lb:ub,lb:ub) = 0d0 
            do k=lb,ub
                dispgeoms(j+npts)%energy(k,k) = eners(k,j)+eshift
            end do
            dispgeoms(j+npts)%lb = lb
            dispgeoms(j+npts)%ub = ub
            call genEnerGroups(dispgeoms(j+npts),deg_cap)
            hasEner(j+npts,lb:ub)=.true.
        enddo
    else !invalid energy file name
        if(printlvl>1)print 1001,"- Skipping energy file input! Filename not present or invalid."
    end if! valid energy file name

    ! read in gradients and couplings data
    if(printlvl>1)print 1001,'- Reading gradients and couplings'
    do j = lb,ub
      do k = lb,ub
        infile = trim(SearchPath(i))//'/'//filename(j,k,grdfptn(i),cpfptn(i))
        if(printlvl>3) print 1000,"Searching for gradients in <"//trim(infile)//">..."
        ptinfile=npoints
        call readGrads(infile,ptinfile,natoms,cgrads)
        if(ptinfile==0)then
            if(printlvl>3)print 1000,"  gradient data not found."
            cycle
        endif!
        if(ptinfile>nnew)then
            print *,"WARNING : gradient file contains more entries than geometry file! Data truncated."
            ptinfile=nnew
        elseif (ptinfile<nnew)then
            print *,"WARNING : gradient file contains less entries than geometry file! "
        end if
        if(printlvl>3.or.(printlvl>1.and.ptinfile<nnew))print 1002,"found ",ptinfile, " gradient data entries." 
        do l = 1,ptinfile
          if(j/=k.and.usefij)then
            dispgeoms(l+npts)%grads(:,j,k)=cgrads(:,l)*(eners(j,l)-eners(k,l))
          else
            dispgeoms(l+npts)%grads(:,j,k)=cgrads(:,l)
          end if
        !    call removeTransRot(dispgeoms(l)%grads(:,j,k),dispgeoms(l)%cgeom)
          if(j/=k)dispgeoms(l+npts)%grads(:,k,j)=dispgeoms(l+npts)%grads(:,j,k)
          hasGrad(l+npts,j,k)=.true.
          hasGrad(l+npts,k,j)=.true.
        enddo!l
      enddo!k
    enddo!j
    npts=npts+nnew
    if(npts>npoints)stop "readdisps BUG: Shouldn't get here"
    if(npts==npoints)exit
  end do !i=1,NSearchPaths
  if(npts<npoints)then
    print *,"WARNING : Found less data points than specified.  Adjusting npoints to ",npts
    npoints=npts
  end if

  ! process the data that have been read in
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
    if(printlvl>1)then
        print *,"  Contribution of each nascent coordinate to B matrix"
        print "(15F8.3)",(dnrm2(3*natoms,dispgeoms(l)%bmat(j,1),ncoord),j=1,ncoord)
    end if
    if(printlvl>0)print *,"      Constructing local coordinate system for point ",l
    CALL makeLocalIntCoord(dispgeoms(l),nstates,useIntGrad,intGradT,intGradS,nvibs,gScaleMode)
    if(printlvl>3)then
       print *,"  Wilson B Matrix in fitting coordinate system"
       do j=1,ncoord
         print "(15F8.3)",dispgeoms(l)%bmat(j,:)
       end do
    end if
    if(printlvl>1)then
       print *,"  Energy gradients in fitting coordinate system:"
       do j=1,nstates
         if(.not.hasGrad(l,j,j))cycle
         write(unit=*,fmt="(A,I3,A)",advance='no')"   state ",j," : "
         print "(15F8.3)",dispgeoms(l)%grads(:,j,j)
       end do
       print *,"  Derivative couplings in fitting coordinate system:"
       do j=1,nstates-1
         do k=j+1,nstates
           if(.not.hasGrad(l,j,k))cycle
           write(unit=*,fmt="(A,I2,A,I2,A)",advance='no')"   block(",j,",",k,") : "
           print "(15F8.3)",dispgeoms(l)%grads(:,j,k)
         end do!k
       end do!j
    end if!printlvl>1
  end do!l

 ! generate intersection adapated coordinates
  if(printlvl>0) print *,"Transforming degenerate points to intersection adapted coordinates:"    
  do l=1,npoints
   call OrthGH_ab(dispgeoms(l),100,gcutoff,hasGrad(l,:,:))
  end do

  call printDisps(int(1),npoints)

  if(allocated(ptWeights))deallocate(ptWeights)
  allocate(ptWeights(npoints))
  call readDispOptions(exclEner,exclGrad,exactEner,exactGrad,exactDiff,enfGO,ptWeights)

1000 format(7X,A)
1001 format(4X,A)
1002 format(10X,A,I6,A)
end SUBROUTINE readdisps
