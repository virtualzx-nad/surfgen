!---------------------------------------------
! Take linear step along vector dvec, then acoording to the new gradients,
! refine the estimation of the optimal step length along this line.  
! Repeat the procedure for a fixed number of times or until the point with
! minimum Lagrangian norm is found
! [ Arguments ]
! hvec  [in/out] DOUBLE PRECISION,dimension(ncons)
!       On input: Current Hd coefficient vector
!       On output: Updated Hd coefficient vector
! dvec  [in/out] DOUBLE PRECISION,dimension(ncons)
!       On input: initial dispacement vector for hvec
!       On output: estimated best displacement vec
! nStep [in] INTEGER
!       Number of linear steps to take
! dCi   [in/out] DOUBLE PRECISION,dimension(ncons)
!       On input: Coefficient part of the Lagrangian gradient at the 
!       initial Hd
!       On output: Coefficient part of the Lagragian gradient after linear
!       steps
! dLam  [in/out] DOUBLE PRECISION,dimension(nex)
!       On input: Lagrange multipliers part of the Lagrangian gradients 
!       of the initial Hd (lambda)
!       On output: Lagrange multipliers part of the Lagrangian gradients
!       of the Hd after linear steps.
! bestB [out] DOUBLE PRECISION,dimension(nex,ncons)
!       Optimized coef-multiplier block of Hessian of the optimal displacement
! dtol  [in] DOUBLE PRECISION
!       Convergence threshold for step length refinement.  When change in
!       relative step length is smaller than this number, convergence considered
!       to have been reached.
  SUBROUTINE takeLinStep(hvec, dvec, nStep, dCi, dLam,bestB,dtol)
    use progdata, only: printlvl
    use makesurfdata
    IMPLICIT NONE
    DOUBLE PRECISION,intent(INOUT)       ::  hvec(ncons),dvec(ncons)
    INTEGER,intent(IN)                   ::  nStep
    DOUBLE PRECISION,intent(INOUT)       ::  dCi(ncons),dLam(nex)
    DOUBLE PRECISION,intent(OUT)         ::  bestB(nex,ncons)
    DOUBLE PRECISION,intent(IN)          ::  dtol
 
! An linear two point extrapolation is used to refine steps. The subroutine stores
! the best two hd vectors and their Lagrangian norms. The two best ones will be
! used to do the interpolation/extrapolation
! Note that the second best is forced to be adjacent to the best one.  That is, 
! the second displacement is chosen as the lower gradient displacement in one of
! the two displacements adjacent to the best.
! The subroutine stores all the displacements that has been performed, as well
! as the gradients at each of the dispacement.
    double precision,dimension(nStep+2)        :: d,g
    double precision,dimension(ncons,nStep+2)  :: dCiLs
    double precision,dimension(nex,nStep+2)    :: dLamLs
    integer :: origin(2,nStep+2)
    integer          :: i, plvl, ndisp,best
    double precision,external     :: dnrm2
! The seek list for displacements.  seek(i) gives the index of the smallest
! displacement that is larger than displacement i.  A value of -1 means the 
! current element is the largest displacement.  seek(0) is the index of the 
! smallest displacement.
    integer,dimension(0:nStep+2)   :: seek
    integer :: nbest,m1,m2
! This list specifies a number of lengths to try.  These trials will all be
! performed before the program use the best two points to perform the
! extrapolation/interpolations.
    integer,parameter                   :: nDivide = 1
    double precision,dimension(nDivide) :: TrialStack
    integer                             :: nTrials

    logical :: found,converged
    double precision  :: dbest
! initialization to have initial point and initial displacement set up to be the
! first two points.
    origin     = 0
    converged  = .false.
    plvl       = printlvl
    printlvl   = 0 
    ndisp      = 1
    nbest      = 1
    seek(0)    = 1
    seek(1)    = -1 
    d(1)       = 0d0
    dCiLs(:,1) = dCi
    dLamLs(:,1)= dLam
    g(1)       = sqrt(dot_product(dCi,dCi)+dot_product(dLam,dLam))
    if(plvl>1)print "(A,E14.7,A,E14.7)",&
             "    Size of displacement=",dnrm2(ncons,dvec,1),", initial grad:",g(1)
    call addDisp(1d0)

    nTrials = nDivide
    do i=1,nDivide
       TrialStack(i) = i/dble(nDivide+1)
    end do
    ! perform nStep refines of the linear displacement
    do i=1,nStep
      if(nTrials>0)then  !do the interpolations specified in the trial stacks
        call addDisp(TrialStack(nTrials))
        nTrials=nTrials-1
      else ! stack is empty , use the two best to do interpolations
        call getBestDisp(dbest,m1,m2,found)
        if(.not.found)exit
        call addDisp(dbest)
        origin(1,ndisp)=m1
        origin(2,ndisp)=m2
      end if 
      if(converged)exit
    end do!i
    if(plvl>1)print "(5x,2(A,E13.5))",&
                "    Final disp scaling factor ",d(nbest)," with gradient",g(nbest)
    dCi  = dCiLs(:,nbest)
    dLam = dLamLs(:,nbest)
    dvec = d(nbest)*dvec
    hvec = hvec+dvec
    printlvl = plvl
  CONTAINS
  !------
    subroutine addDisp(dnew)
      implicit none
      double precision,intent(in) :: dnew
      double precision :: hnew(ncons+nex),lag,jaco(nex,ncons),jaco2(nex,ncons)
      integer  :: iseek,ilast
      double precision :: nrmener, avgener, nrmgrad, avggrad, LSErr,ExErr
      if(ndisp==nStep+2)stop"bug in takeLinStep::addDisp()  :  ndisp overflow"
      ndisp=ndisp+1
      d(ndisp) = dnew     
      hnew(1:ncons) = hvec+dnew*dvec
      hnew(ncons+1:)= 0d0
      call evaluateError(hnew,weight,LSErr,ExErr)
      CALL getError(nrmener,avgener,nrmgrad,avggrad)
      CALL getCGrad(hnew,dCiLs(1,ndisp),dLamLs(1,ndisp),lag,jaco)
      CALL optLag(jaco,nex,dCiLs(1,ndisp),hnew,jaco2)
      g(ndisp) = sqrt(sum(dCiLs(1:ncons,ndisp)**2)+sum(dLamLs(1:nex,ndisp)**2))
      ! fild the insertion point for seek list 
      ilast=0
      iseek=seek(0)
      do while(d(iseek)<dnew)
        ilast = iseek
        iseek = seek(iseek)
        if(iseek<0)exit
      end do
      ! insert the index of current 
      seek(ndisp) = iseek
      seek(ilast)  = ndisp
      if(g(nbest)>g(ndisp))then
        if(abs(d(nbest)-d(ndisp))<dtol)converged=.true.
        nbest=ndisp
        bestB=jaco2
      end if
      if(plvl>1)print "(A,I3,A,F8.5,A,2E12.5,A,2E12.5,A,E16.9)",&
                        "    disp#",ndisp-1,"dnew=",dnew,": evalEr=",LSErr,ExErr,&
               ",getEr=",nrmener,nrmgrad,",GLag=",g(ndisp)
    end subroutine addDisp
  !------
  ! search through the estimated displacement between any adjacent data points
  ! to find the estimated best displacement
    subroutine getBestDisp(dnew,m1,m2,found)
      implicit none
      double precision,intent(out) :: dnew
      logical,intent(out)          :: found
      integer,intent(out)          :: m1,m2
      integer :: iseek,i
      logical :: success
      double precision  ::  gbest, dcurr, gcurr
      gbest = -1
      found = .false.
      iseek = seek(0)
      do while(seek(iseek)>0)  ! seek through all adjacent elements
        success = .true.
        do i=1,ndisp
          if(origin(1,i)==iseek.and.origin(2,i)==seek(iseek))success=.false.
        end do
        if(success)then
          call interpStep(iseek,seek(iseek),dcurr,gcurr,success)
          if(success.and.gbest<0.or.gcurr<gbest)then!new displacement is better
            gbest = gcurr
            dnew  = dcurr
            m1 = iseek
            m2 = seek(iseek)
            found = .true.
          end if!new displacement is better
        end if!new pair
        iseek = seek(iseek)
      end do !while(seek(iseek)>0)
!print "(A,E13.5,A,2I4)","estimated best gradient",gbest,"from pts",m1,m2
    end subroutine getBestDisp
  !------
  ! estimate the best displacement point using data of step m and n.  also gives
  ! the estimated gradients
  ! The estimation is made by assuming linearity of Lagrangian gradients with
  ! respect to displacements.
    subroutine interpStep(m,n,dnew,gnew,success)
      implicit none
      integer, intent(in)           :: m,n
      double precision,intent(out)  :: dnew,gnew
      logical,intent(out)           :: success
      double precision,dimension(nex+ncons) :: dgvec,g0vec 
      double precision  :: d2
      double precision,external  :: dnrm2
      dgvec(ncons+1:ncons+nex) = dLamLs(:,m)-dLamLs(:,n)
      dgvec(1:ncons)           = dCiLs(:,m) -dCiLs(:,n)
      dgvec = dgvec/(d(m)-d(n))
      d2  =  sum(dgvec**2)
      if(d2>1d-12)then
        g0vec(1:ncons)           = d(m)*dCiLs(:,n)-d(n)*dCiLs(:,m)
        g0vec(ncons+1:ncons+nex) = d(m)*dLamLs(:,n)-d(n)*dLamLs(:,m)
        g0vec   = g0vec/(d(m)-d(n))
        dnew = -dot_product(g0vec,dgvec)/d2
        gnew = dnrm2(nex+ncons,g0vec+dnew*dgvec,1)
        success = .true.
      else
      ! its a constant gradient displacement.  can't make any predictions
        success = .false.
      end if
    end subroutine interpStep
  END SUBROUTINE takeLinStep
