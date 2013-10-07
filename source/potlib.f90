! POTLIB interface for surfgen
! Module stores variables related to pot() calls
MODULE potdata
! dcoordls    INTEGER,dimension(*)
!             Specifies the set of curvilinear coordinates that will
!             be used to calculate the distance factor.   Must contain
!             NRij entries and in the same order as the Rij definition
!             generated in libsym
     INTEGER,PARAMETER  ::  MaxNDCoord= 200
     INTEGER            ::  ndcoord
     INTEGER            ::  dcoordls(MaxNDCoord)

! diagonal shift potential.  currently only supports one coordinate functions.
     INTEGER,PARAMETER  ::  MaxNfTerms=100
     INTEGER            ::  nfterms  ! number of shift functions
     INTEGER            ::  fterm(MaxNfTerms)  ! index of each shift function
     DOUBLE PRECISION   ::  fcoef(MaxNfTerms)  ! expansion coefficient
     DOUBLE PRECISION   ::  forig(MaxNfTerms)  ! place where shift is 0

! how much the energy will be shifted
     DOUBLE PRECISION  :: eshift

     LOGICAL       ::  initialized = .false.
     INTEGER       ::  GUNIT       = -1
     INTEGER       ::  NEval       = 0
     INTEGER       ::  MUnit  
!# of points before geometry is recorded in molden output 
     INTEGER       ::  molden_p
!# of points before starting to record in molden file
     INTEGER       ::  m_start
     INTEGER       ::  nrec = 0

! switch on calculation of min distance to a set of know geometries
     LOGICAL       ::  calcmind
! maximum deviation amount points of the current trajectory
     DOUBLE PRECISION :: mdev

! name of file that contains reference geometries
     CHARACTER(255)::  gflname
     DOUBLE PRECISION :: mindcutoff
     INTEGER       ::  nrpts   !# of points in reference geom file
! reference geometry in scaled rij coordinates
     DOUBLE PRECISION,dimension(:,:),allocatable   :: rgeom

! register for eigenvectors to the Schrodinger equations
     DOUBLE PRECISION,dimension(:,:),allocatable   :: evecstore

! energy and gradient fitting error for all the reference geometries
     DOUBLE PRECISION,dimension(:,:),allocatable   :: enererrdata
     DOUBLE PRECISION,dimension(:,:,:),allocatable :: graderrdata

! atom labels of each of the atoms for molden file generation
     CHARACTER(3),dimension(20)                   :: atomlabels

!variables used in the search for minimum distance 
! lastrgeom   Scaled rij geometry of the last point evaluated
! ldbounds    Lower bounds of distances from the last point to all ref points
! udbounds    Upper bounds of distances from the last point to all ref points
     DOUBLE PRECISION,dimension(:),allocatable,PRIVATE    :: lastrgeom
     DOUBLE PRECISION,dimension(:),allocatable,PRIVATE    :: ldbounds,udbounds
     LOGICAL,PRIVATE :: FirstPt

! this parameter enables parsing and distance checks
! can be modified through subroutines enableParsing() and disableParsing()
     LOGICAL                                ::  parsing
! this parameter specifies if parsing is temporarily paused
     LOGICAL                                ::  paused 

! These parameters controls the estimation of energy errors over evaluation points
! calcErr     Whether or not energy error will be estimated
! errflname   Name of the file that contains fitting energy error and gradient error
!             of all the data points
     LOGICAL                                ::  calcErr
     CHARACTER(255)                         ::  errflname

! Surface number and time data used for outputfile generation
     DOUBLE PRECISION                       ::  timetraj
     INTEGER                                ::  isurftraj

! Switch for timing test outputs
! When set to .true., EvaluateSurfgen will time each of the following procedures:
!    buildWBMat EvalRawTerms EvalHdDirect DSYEVR int2cart CpScaling analysis
! and report time spent in each of the sections in ms
     LOGICAL                                ::  timeEval

! Threshold for a coordinate to be considered vanished.   Norm of B-matrix
! elements is used to check if a coordinate is vanished.   The distance output
! in trajdata file is the RMS distance over number of non-vanishing coordinates
! at the current evaluation point.
! Number of non-vanishing coordinates is also appended in the output file. 
     DOUBLE PRECISION                       ::  cvanish
    
! EvalCount stores the evaluation counts.   
! Use ResetEvalCount() and GetEvalCount() to reset and retrieve the counter
     INTEGER                                ::  EvalCount
!-----------------------------------------------------------------------------
CONTAINS
!-----------------------------------------------------------------------------
! initialze the module
     SUBROUTINE init(ncoord)
       IMPLICIT NONE
       INTEGER,INTENT(IN) :: ncoord
       integer :: i
       double precision, external :: dnrm2
       if(allocated(lastrgeom))deallocate(lastrgeom)
       if(allocated(ldbounds))deallocate(ldbounds)
       if(allocated(udbounds))deallocate(udbounds)
       allocate(lastrgeom(ncoord))
       allocate(ldbounds(nrpts))
       allocate(udbounds(nrpts))
   ! set the origin to be initial geometry and set distances accordingly
       lastrgeom=dble(0)
       do i=1,nrpts
         ldbounds(i)=dble(0)
         udbounds(i)=dble(0)
       end do
       firstPt=.true.
       timetraj  = -1D0
       isurftraj = 0
       mdev = -1D0
     END SUBROUTINE init
!-----------------------------------------------------------------------------
! calculate the distance, using a set of internal coordinate of choice,
! between two points
!
! Internal coordinates used to evaluate the distance is specified with 
! `dcoordls`.   
! Distance will be calculated using all possible permutations
! of the atoms and the smallest one will be used.
! coordinate permutation list generated in CNPI module is used 
! if ptid>0, energy error will be estimated and output through estErr
  SUBROUTINE getdist(rgeom1,rgeom2,min_d2,ptid,estErr)
    USE CNPI,     ONLY: coordPerm,nPmt
    USE hddata,   ONLY: ncoord,nstates
    IMPLICIT NONE

    DOUBLE PRECISION,dimension(ncoord),intent(IN)     :: rgeom1,rgeom2
    DOUBLE PRECISION,intent(OUT)                      :: min_d2
    INTEGER, intent(IN)                               :: ptid
    DOUBLE PRECISION,dimension(nstates),intent(INOUT) :: estErr

    double precision,dimension(ndcoord)         :: rgeomp,minrgeomp    
    double precision                            :: d2
    DOUBLE PRECISION,dimension(ncoord)          :: rgeomtmp
    integer  :: i
    rgeomtmp = rgeom1(coordperm(1,dcoordls(1:ndcoord)))
    rgeomp=rgeomtmp(dcoordls(1:ndcoord))-rgeom2(dcoordls(1:ndcoord))
    min_d2=dot_product(rgeomp,rgeomp)
    do i=2,nPmt
      rgeomtmp = rgeom1(coordperm(i,dcoordls(1:ndcoord)))
      rgeomp=rgeomtmp(dcoordls(1:ndcoord))-rgeom2(dcoordls(1:ndcoord))
      d2=dot_product(rgeomp,rgeomp)
      if(d2<min_d2)then
        min_d2=d2
        minrgeomp = rgeomp
      end if
    end do
    min_d2=dsqrt(min_d2)
    if(ptid>0)then   ! estimate energy error
      if(calcerr)then
        do i=1,nstates
          estErr(i)=enererrdata(ptid,i)+&
             dot_product(graderrdata(ptid,i,dcoordls(1:ndcoord)),minrgeomp)
        end do
      else
        estErr = 0d0
      end if
    end if!ptid>0
  END SUBROUTINE getdist
!-----------------------------------------------------------------------------
! Measures the minimum distance of a point defined by cartesian 
! geometry of all the atoms to a set of reference points.
!-----------------------------------------------------------------------------
! ALGORITHM   The subroutine elaborates to minimize distance calculations.
! The module keeps a list of lower bounds and upper bounds of distances 
! to all the reference geometries of the last point, as well as the geometry of
! the last point.    When a new point is supplied, the distance to the last 
! point is first calculated, and the bounds updated using triangular inequality.
! The upper bound of least distance(max(minD)) is the lowest of the upper bounds.
! A list of points where the lower bound is less than max(minD) is generated.
! The subroutine then proceeds to execute a loop.  
! Each iteration, distance to the point with least lower bound is calculated, the
! range of distance to this point updated.  This distance will serve as the new
! max(minD) if it is lower, or the point will be eliminated.  
! The list is updated each iteration, until minD is found.
!     This algorithm works the Best when the newly point is always very close 
! to the last point evaluated.  This is always the case for trajectory calculations.
!-----------------------------------------------------------------------------
! cgeom    [input] DOUBLE PRECISION,dimension(3*natoms)
!          Cartesian geometry of all the atoms for the point to be tested
! min_d    [output] DOUBLE PRECISION
!          Minimum distance to existing points
! ptid     [output] INTEGER
!          Index of the point associated with distance min_d
! eerr     [output]  DOUBLE PRECISION
!          Energy error estimated from nearest reference point fitting error
!-----------------------------------------------------------------------------
  SUBROUTINE getmind(currgeom,min_d,ptid,eerr)
    USE progdata, ONLY: natoms,AU2CM1
    USE hddata,   ONLY: ncoord,nstates
    IMPLICIT NONE
    DOUBLE PRECISION,dimension(ncoord),intent(IN)   :: currgeom
    DOUBLE PRECISION,intent(OUT)                    :: min_d
    INTEGER,intent(OUT)                             :: ptid
    DOUBLE PRECISION,dimension(nstates),intent(OUT) :: eerr

    double precision :: dlpt     !distance to last point
    double precision :: mmd      !max(minD)
    integer  ::  WL(nrpts)  ! waiting list for candidates that worth evaluation
    integer  ::  NWL        ! # of points in WL
    double precision  ::  d ! distance to the point with least lower bound
    integer  :: i,j
    double precision,parameter :: vsmall = 1D-8  ! resolution of distances
    
    if(firstpt)then
      do i=1,nrpts
        CALL getdist(currgeom,rgeom(i,:),d,int(0),eerr)
        ldbounds(i)=d
        udbounds(i)=d
        if(i.eq.1.or.d<min_d)then
           min_d= d
           ptid = i
        end if
      end do
      lastrgeom=currgeom
      firstpt=.false.
      return
    end if!firstpt
! Evaluate distance to last point
    CALL getdist(currgeom,lastrgeom,dlpt,int(0),eerr) 
! Update lastrgeom and distance bounds 
    lastrgeom=currgeom
    do i=1,nrpts
      if(dlpt>udbounds(i))then
        ldbounds(i)=dlpt-udbounds(i)
      else if(dlpt<ldbounds(i)) then
        ldbounds(i)=ldbounds(i)-dlpt
      else
        ldbounds(i)=dble(0)
      end if!(dlpt>udbounds(i))
    end do !i=1,nrpts
    udbounds=udbounds+dlpt
! prepare the candidate list
    mmd = minval(udbounds)
    nwl = 0
    do i=1,nrpts
      if(ldbounds(i)<mmd+vsmall)then
        nwl     = nwl+1
        wl(nwl) = i
      end if !(ldbounds(i)<mmd+vsmall)
    end do!i=1,npts
! Evaluate these points that may be lowest 
    min_d = mmd
    ptid  = 0 
    do i=1,nwl
      if(ldbounds(wl(i))>min_d)cycle
      CALL getdist(currgeom,rgeom(wl(i),:),d,int(0),eerr)
      ldbounds(wl(i))=d
      udbounds(wl(i))=d
      if(d<min_d.or.ptid.eq.0)then
        min_d = d
        ptid  = wl(i) 
      end iF
    end do !while (nwl>1)
    call getdist(currgeom,rgeom(ptid,:),d,ptid,eerr)
    eerr = eerr *AU2CM1
  END SUBROUTINE getmind
!--------------------------------------------------------------------
! Reset the evaluation counter
  SUBROUTINE ResetEvalCount()
     EvalCount = 0
  END SUBROUTINE ResetEvalCount
!--------------------------------------------------------------------
! Retrieve the evaluation counter
  SUBROUTINE GetEvalCount(count)
     INTEGER,INTENT(OUT)  :: count
     count= EvalCount
  END SUBROUTINE GetEvalCount
END MODULE potdata

!--------------------------------------------------------------------
SUBROUTINE finalizeSurfgen()
  USE potdata
  INTEGER   CNT
  if(mdev>0) print "(A,E10.3)","Maximum deviation of last trajectory:",mdev
  call GetEvalCount(cnt)
  print "(A,I8)","Evaluation count for last trajectory:",cnt
END SUBROUTINE
!--------------------------------------------------------------------
! create the files corresponding to a new trajectory
SUBROUTINE openTrajFile(itraj)
  USE hddata, ONLY: nstates
  USE progdata, ONLY: natoms,atoms
  USE potdata
  IMPLICIT NONE
  INTEGER,intent(IN)     ::  itraj
  integer cnt
  character(4)   ::  str

! output analysis of previous trajectories.
  if(itraj>0.and.mdev>0)then
    print "(A,I4,A,E10.3)","Maximum deviation for trajectory ",itraj," : ",mdev 
    mdev = -1D0
    call GetEvalCount(cnt)
    print "(A,I8,A)","Last trajectory finished with ",cnt," evaluations. "
    call ResetEvalCount()
  end if

  write(str,"(I4)")itraj 

! close old files
 
  if(molden_p>0)close(MUNIT)

! prepare output for geometry, energy and distance print out
  if(GUNIT<0)then
    GUNIT = 427
  else
    close(GUNIT)
  end if
  open(GUNIT,file='trajdata'//trim(adjustl(str))//'.csv',status='REPLACE',action='write')

! prepare molden output file
  if(molden_p>0)then
    open(MUnit,file='molden'//trim(adjustl(str))//'.all',status='REPLACE',action='write')
    write(MUnit,'(A)')" [Molden Format]"
    write(MUnit,'(A)')" [GEOMETRIES] XYZ"
  end if !molden_p>0

  write(str,"(I4)")natoms 
  write(GUNIT,&
           "(I5,',',I5,"//trim(adjustl(str))//"(',',A3))")natoms,nstates,atomlabels(atoms(1:natoms))
END SUBROUTINE

!---------------------------------------------------------------------
! initializes Hd on entrance 
SUBROUTINE prepot
     USE hddata, ONLY: getFLUnit, ncoord,nstates
     USE potdata
     IMPLICIT NONE
     integer   ::jobtype

     character(72) :: ver

     print *,"Initializing surfgen"

     call getver(ver)

     print "(2A)","Program version: ",trim(ver)
     call readginput(jobtype)
     call initialize(jobtype)

     if(calcmind)call loadRefGeoms()
     if(allocated(evecstore))deallocate(evecstore)
     allocate(evecstore(nstates,nstates)) 
     paused = .false. 
     if(parsing)then 
       GUNIT = 427
       open(GUNIT,file='trajdata0.csv',status='REPLACE',action='write')
     else
       GUNIT = -1
     end if
! prepare molden output file
     if(molden_p>0)then
        MUNIT=getFLUnit()
        open(MUnit,file='molden0.all',status='REPLACE',action='write')
        write(MUnit,'(A)')" [Molden Format]"
        write(MUnit,'(A)')" [GEOMETRIES] XYZ"
     end if

! initialize important variables
     print *," Initializing module potlib..."
     initialized = .true.
     NEval       = 0
     nrec        = 0
     call init(ncoord)
     print *," Initialization complete."
    
     call ResetEvalCount()
     return
END SUBROUTINE prepot
!---------------------------------------------------------------------
! samething as prepot just a different name
SUBROUTINE initPotential()
  CALL prepot
END SUBROUTINE
!---------------------------------------------------------------------
 ! getData returns the Hd predicted energies and gradients
 ! cartgeom [input] DOUBLE PRECISION, dimension(3*natoms)
 !          Cartesian coordinates of all the atoms.
 ! E        [output] DOUBLE PRECISION, dimension(nstates)
 !          Adiabatic energy predicted by Hd for each of the states
 ! dH       [output] DOUBLE PRECISION, dimension(3*natoms,nstates,nstates)
 !          Matrix elements of derivative of Hd <i|dH/dR|j>
 ! If distance calculation is enabled, each iteration getMinD() will be
 ! called to locate the nearest data point within the point set used to
 ! fit Hd.  The distance is reported.
 ! These numbers are printed out into the geom.all file, along with the
 ! geometry and adiabatic energies, for analysis after simulation.
 ! If molden_p>0, molden output will also be generated every molden_p
 ! evaluations, starting from the (m_start+1)th evaluation.
 subroutine getdata(cartgeom,E,dH)
   USE progdata, ONLY: natoms
   USE hddata, ONLY: nstates
   implicit none
   real*8,intent(in),dimension(3*natoms)                  ::  cartgeom
   real*8,intent(out),dimension(nstates)                  ::  E
   real*8,intent(out),dimension(3*natoms,nstates,nstates) ::  dH
   DOUBLE PRECISION,dimension(nstates,nstates)            ::  hmat
   DOUBLE PRECISION,dimension(3*natoms,nstates,nstates)   ::  dcgrads

   CALL EvaluateSurfgen(cartgeom,E,dH,hmat,dcgrads)
 end subroutine
!---------------------------------------------------------------------
! Wrapper for fotran 77 callers where array size may need to be fixed
! parameters.  The subroutine takes input arrays with the size passed
! as argument, repack them then calls EvaluateSurfgen using repacked
! arrays.  
!
! Arguments
! ---------
! mnatm [in] INTEGER
!       Maximum number of atoms. Used as dimensionality of cgeos,Ga and Gd
! mnst  [in] INTEGER 
!       Maximum number of states, used as dimensionality of arrays
! cgeom [in] DOUBLE PRECISION, dimension(3*mnat)
!       Cartesian geometry of the molecule.
! Ea    [out] DOUBLE PRECISION,dimension(mnst)
!       Adiabatic energies of all electronic states.
! Ga    [out] DOUBLE PRECISION,dimension(3*mnat,mnst,mnst)
!       The diagonal contains adiabatic energy gradients, and off-diagonal 
!       contains derivative coupling vectors.
! Ed    [out] DOUBLE PRECISION,dimension(mnst,mnst)
!       Diabatic energy matrix
! Gd    [out] DOUBLE PRECISION,dimension(3*mnat,mnst,mnst)
!       Gradient of diabatic energies.
 subroutine EvaluateSurfgen77(mnatm,mnst,cgeom,Ea,Ga,Ed,Gd)
   USE progdata, ONLY: natoms
   USE hddata, ONLY: nstates
   integer, intent(in)          ::  mnatm,mnst
   double precision,intent(in)  ::  cgeom(3*mnatm)
   double precision,intent(out) ::  Ea(mnst),Ed(mnst,mnst)
   double precision,intent(out),dimension(3*mnatm,mnst,mnst) :: Ga,Gd
! local packed variables
   double precision :: cgeompck(3*natoms),edpck(nstates,nstates)
   double precision,dimension(3*natoms,nstates,nstates) :: gapck,gdpck

   if(natoms>mnatm.or.nstates>mnst)&
     stop "EvaluateSurgen: Maximum on atom or state count insufficient."
   cgeompck = cgeom(1:3*natoms)
   call EvaluateSurfgen(cgeompck,Ea,gapck,edpck,gdpck)
   Ed(1:nstates,1:nstates) = edpck
   Ga(1:3*natoms,1:nstates,1:nstates) = gapck
   Gd(1:3*natoms,1:nstates,1:nstates) = gdpck
 end subroutine EvaluateSurfgen77
!---------------------------------------------------------------------
! POTLIB interface in compliance with ANT 09's NH3 potential 
! works only for 2 state systems
! Xcart          [input] DOUBLE PRECISION, dimension(3*natoms)
!                Cartesian coordinates of all the atoms.  Coordinates
!                of atom i are stored in units i*3-2 to i*3.
! U11,U22,U12    [output] DOUBLE PRECISION
!                Diabatic energies and off-diagonal terms
! V1,V2          [output] DOUBLE PRECISION
!                Adiabatic energies of both states
! gU11,gU22,gU12 [output] DOUBLE PRECISION, dimension(3*natoms)
!                Gradients of U11,U22 and U12
! gV1,gV2        [output] DOUBLE PRECISION, dimension(3*natoms)
SUBROUTINE pot(Xcart,U11,U22,U12,V1,V2,gU11,gU22,gU12,gV1,gV2,gV12)
  USE progdata, ONLY : natoms,switchdiab
  USE potdata
  USE hddata, ONLY:  nstates
  IMPLICIT NONE
  DOUBLE PRECISION,dimension(3*natoms),intent(IN)   ::  Xcart
  DOUBLE PRECISION,dimension(3*natoms),intent(OUT)  ::  gV1,gV2,gV12,gU11,gU22,gU12
  DOUBLE PRECISION,intent(OUT)                      ::   V1, V2, U11, U22, U12

  double precision,dimension(nstates,nstates)           ::  hmat
  double precision,dimension(3*natoms,nstates,nstates)  :: cgrads,dcgrads
  double precision,dimension(nstates)                   ::  energy
 
  CALL EvaluateSurfgen(Xcart,energy,cgrads,hmat,dcgrads) 
  V1=energy(1)+eshift
  V2=energy(2)+eshift
  gV1=cgrads(:,1,1)
  gV2=cgrads(:,2,2)
  gV12 = cgrads(:,1,2)
  if(switchdiab)then !if the first two diabats should be switched to match adiabats
    U22  = hmat(1,1)+eshift
    U12  = hmat(1,2)
    U11  = hmat(2,2)+eshift
    gU22 = dcgrads(:,1,1)
    gU12 = dcgrads(:,1,2)
    gU11 = dcgrads(:,2,2)
  else
    U11  = hmat(1,1)+eshift
    U12  = hmat(1,2)
    U22  = hmat(2,2)+eshift
    gU11 = dcgrads(:,1,1)
    gU12 = dcgrads(:,1,2)
    gU22 = dcgrads(:,2,2)
  end if
  return
END SUBROUTINE pot
! --------------------------------------------------------------
SUBROUTINE getEvec(evec)
  use hddata, only: nstates
  use potdata, only: evecstore
  DOUBLE PRECISION, dimension(nstates,nstates),intent(out):: evec
  evec = evecstore
END SUBROUTINE getEvec
! --------------------------------------------------------------
![Description]
! EvaluateSurfgen evaluates Hd at a given Cartesian coordinate, returning adiabatic
! energies, energy gradients, derivative couplings, diabatic Hamiltonian matrix
! and gradients of each of the blocks of Hd.   Vectors are given in the same
! Cartesian coordinates.
!
! This is intended to be the interface subroutine used by other programs.
!  *PARAMETER "SWITCHDIAB" NOT IMPLEMENTED FOR THIS SUBROUTINE
!
![Arguments]
! cgeom    [input] DOUBLE PRECISION,dimension(3*natoms)
!          Cartesian geometries in a flattened 1-dimensional array.  Coordinate
!          of atom ($i) is located from items ($i)*3-2 to ($i)*3
! energy   [output] DOUBLE PRECISION,dimension(natoms)
!          Adiabatic energies of each of the states.
! cgrads   [output] DOUBLE PRECISION,dimension(3*natoms,nstates,nstates)
!          The vector cgrads(:,state1,state2) contains (in Cartesian coordinate)
!          Adiabatic energy gradients for state1 when state1.eq.state2
!          Derivative coupling between state1 and state2 when state1.ne.state2
! hmat     [output] DOUBLE PRECISION,dimension(nstates,nstates)
!          Value of diabatic Hamiltonian matrix
! dcgrads  [output] DOUBLE PRECISION,dimension(3*natoms,nstates,nstates)
!          Cartesian gradients each of the blocks of Hd.
!
![Dependencies]
! MODULES         SOURCE
!   progdata      progdata.f90
!   hddata        hddtat.f90
!   potdata       potlib.f90
! SUBROUTINES     SOURCE
!   getmind       potlib.f90
!   prepot        potlib.f90
!   buildWBmat    libinternal.f90
!   EvalRawTerms  hddata.f90
!   EvalHdDirect  hddata.f90
!   DSYEVR        LAPACK
SUBROUTINE EvaluateSurfgen(cgeom,energy,cgrads,hmat,dcgrads)
  USE progdata, ONLY : natoms,atoms
  USE potdata
  USE hddata, ONLY:  nstates,ncoord,EvalRawTermsL,EvalHdDirect
  IMPLICIT NONE
  DOUBLE PRECISION,dimension(3*natoms),intent(IN)                  ::  cgeom
  DOUBLE PRECISION,dimension(nstates,nstates),intent(OUT)          ::  hmat
  DOUBLE PRECISION,dimension(3*natoms,nstates,nstates),intent(OUT) ::  cgrads,dcgrads
  DOUBLE PRECISION,dimension(nstates),intent(OUT)                  ::  energy

  double precision,dimension(nstates,nstates)              ::  evec,hmat2,dhtmp
  double precision,dimension(nstates)                      ::  eerr
  double precision,dimension(ncoord,nstates,nstates)       ::  dhmat
  double precision,dimension(ncoord)                ::  igeom
  double precision,dimension(ncoord,3*natoms)       ::  bmat 
  integer                                           :: m,LWORK,LIWORK,INFO
  integer,dimension(nstates*2)                      :: ISUPPZ
  double precision,dimension(nstates*(nstates+26))  :: WORK
  integer,dimension(nstates*10)                     :: IWORK
  double precision  ::  dvec(3),dX(3),minwlist(natoms),dw
  double precision    :: bohr2ang,  mind, minRij,dcsX,csX,eWall,xlist(natoms*(natoms-1)/2),gfactor,dtR
  double precision   :: comP(3), fnorm(nstates),de
  integer   :: i,j,k,ptid,count1,count2,count_rate,minK,minI,minJ,tmp
  integer   :: counter = 1   ! count the number of evaluations
  character(4)  ::  str,str2
  double precision  ::  teval(7)
  double precision  ::  f,df(3*natoms)  !diagonal shift
  double precision,external  :: dnrm2
  integer   :: nth, nc ! number of nonvanishing coordinates
  integer,external  :: omp_get_max_threads

  LWORK  = nstates*(nstates+26)
  LIWORK = nstates*10
  bohr2ang = 0.529177249

  if(.not.initialized)then
    Print *,"Error :   Hd needs to be initialized before evaluation."
    print *,"Please call initPotential() or prepot() once before calling EvaluateSurfgen()"
    stop "Execution aborted."
  end if

  nth=omp_get_max_threads()
  EvalCount = EvalCount+1
  if(timeeval) call system_clock(COUNT=count1,COUNT_RATE=count_rate)
  call buildWBMat(cgeom,igeom,bmat,.false.)
  if(timeeval)then
    call system_clock(COUNT=count2)
    teval(1) = dble(count2-count1)/count_rate*1000
  end if!timeeval

! calculate raw polynomial terms
  call EvalRawTermsL(igeom)
  if(timeeval)then
    call system_clock(COUNT=count1)
    teval(2) = dble(count1-count2)/count_rate*1000
  end if!timeeval

 ! construct Hd and its derivatives
  call EvalHdDirect(hmat,dhmat)
  if(timeeval)then
    call system_clock(COUNT=count2)
    teval(3) = dble(count2-count1)/count_rate*1000
  end if!timeeval

  ! convert gradients and couplings to cartesian coordinates
  dcgrads=dble(0) 
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(I,J,dhtmp)
  do i=1,3*natoms
    do j=1,ncoord
      call daxpy(nstates*nstates,bmat(j,i),dhmat(j,1,1),ncoord,dcgrads(i,1,1),3*natoms)
    end do !j=1,ncoord
  end do !i=1,3*natoms
!$OMP END PARALLEL DO
  if(timeeval)then
    call system_clock(COUNT=count1)
    teval(4) = dble(count1-count2)/count_rate*1000
  end if!timeeval

 ! calculate diagonal shift factor and its derivatives
 f = eshift
 df = 0
 do i=1,nfterms
   f = f+fcoef(i)*(1-igeom(fterm(i))/forig(i))**2
  df =df-2*fcoef(i)*bmat(fterm(i),:)*(1-igeom(fterm(i))/forig(i))/forig(i)
 end do 
 ! perform shift for diagonals
 do i=1,nstates
   hmat(i,i)      = hmat(i,i)+f
   dcgrads(:,i,i) = dcgrads(:,i,i) +df
 end do
 if(timeeval)then
    call system_clock(COUNT=count2)
    teval(6) = dble(count2-count1)/count_rate*1000
 end if!timeeval

  call mkl_set_num_threads(1)
 ! generate eigenvectors and energies at current geometry
  hmat2 = hmat
  CALL DSYEVR('V','A','U',nstates,hmat2,nstates,dble(0),dble(0),0,0,1D-15,m,&
            energy,evec,nstates,ISUPPZ,WORK,LWORK,IWORK,LIWORK, INFO )
  evecstore=evec
  call mkl_set_num_threads(nth)

  do i=1,3*natoms
    cgrads(i,:,:)=matmul(transpose(evec),matmul(dcgrads(i,:,:),evec))
  end do!i=1,ncoord
  do i=1,nstates-1
    do j=i+1,nstates
      de=energy(j)-energy(i)
      if(abs(de)<1d-30)de=1d-30
      call dscal(3*natoms,1/de,cgrads(1,i,j),1)
      call dcopy(3*natoms,cgrads(1,i,j),1,cgrads(1,j,i),1)
!      cgrads(:,i,j) = cgrads(:,i,j)/de
!      cgrads(:,j,i) = cgrads(:,i,j)
    end do
  end do

  if(timeeval)then
    call system_clock(COUNT=count1)
    teval(5) = dble(count1-count2)/count_rate*1000
  end if!timeeval

  if(parsing)then
    NEval=NEval+1
    write(str,"(I4)") natoms*3
    write(str2,"(I4)") nstates
    if(GUNIT<0)call openTrajFile(0)
    if(calcmind.and.ndcoord>0)then
       !get the number of nonvanishing coordinates
       if(cvanish<=0)then
          nc=ndcoord
       else  
          nc=0
          do j=1,ndcoord
            if(dnrm2(3*natoms,bmat(dcoordls(j),1),ncoord)>cvanish)nc=nc+1
          end do
       end if 
       !calculation minimum distances to all reference points
       call getmind(igeom,mind,ptid,eerr)
       mind=mind/sqrt(dble(nc))
       mdev = max(mdev,mind)
       fnorm = 0d0
       if(isurftraj<=nstates.and.isurftraj>0)then
         do j=1,nstates
           if(j==isurftraj)cycle
           fnorm(j) = dnrm2(3*natoms,cgrads(1,j,isurftraj),1)
         end do
       end if
       if(calcerr)then
         write(GUNIT,&
           "(F16.3,',',I9,',',"//trim(str)//"(F12.7,','),"//trim(str2)//"(E16.8,','),F14.8,',',I6,"//&
                                                        trim(str2)//"(',',F16.4),',',I3,"//trim(str2)//"(',',E16.8),',',I4)")&
            timetraj,NEval,                cgeom,                         energy  ,  mind    , ptid,   eerr, isurftraj, &
                       fnorm, nc
       else
         write(GUNIT,&
           "(F16.3,',',I9,',',"//trim(str)//"(F12.7,','),"//trim(str2)//"(E16.8,','),F14.8,',',I6,',',I3,"&
                 //trim(str2)//"(',',E16.8),',',I4)")&
            timetraj, NEval,      cgeom,               energy  ,  mind   , ptid  , isurftraj,fnorm,nc
       end if
    else!calcmind
       write(GUNIT,"(F16.3,',',I9,"//trim(str)//"(',',F12.7),"//trim(str2)//"(',',E16.8),',',I3,"//trim(str2)//"(',',E16.8))")&
                  timetraj,   NEval                ,  cgeom                    ,  energy   , isurftraj, fnorm
    end if !calcmind
    if(molden_p>0.and.NEval-m_start>nrec*molden_p.and.parsing)then
      ! output molden geometries
      nrec=nrec+1
      write(MUnit,'(I3,A)') natoms,"  /coord "
      write(MUnit,'(I5,A)') nrec,"  /current iter "
      do i=1, natoms
        write(MUnit,"(' "//trim(adjustl(atomlabels(atoms(i))))//" ',3F11.6)") cgeom(i*3-2:i*3)*bohr2ang
      end do
    end if!(molden_p>0.and.NEval-m_start>nrec*molden_p.and.parsing)
  end if!(parsing)
  if(timeeval)then
    call system_clock(COUNT=count2)
    teval(7) = dble(count2-count1)/count_rate*1000
    print *,"Execution time of subroutines(ms)"
    print *," buildWBMat EvalRawTerms EvalHdDirect  int2cart   DSYEVR   DiagShift   analysis"
    print   "(7F11.3)",teval
  end if!timeeval
END SUBROUTINE EvaluateSurfgen
!----------------------------------------------------------------------------------
![Description]
! getEnergy retrieves adiabatic energies at a certain geometry
! This subroutine runs much faster than the full version.
! It does not calculate gradients of Hd, nor does not need to perform
! coordinate transformation from internal to cartesian, nor does it
! perform distance and fitting error checks.
!
![Arguments]
! cgeom    [input] DOUBLE PRECISION,dimension(3*natoms)
!          Cartesian geometries in a flattened 1-dimensional array.  Coordinate
!          of atom ($i) is located from items ($i)*3-2 to ($i)*3
! energy   [output] DOUBLE PRECISION,dimension(natoms)
!          Adiabatic energies of each of the states.
!
![Dependencies]
! MODULES         SOURCE
!   progdata      progdata.f90
!   hddata        hddtat.f90
!   potdata       potlib.f90
! SUBROUTINES     SOURCE
!   prepot        potlib.f90
!   buildWBmat    libinternal.f90
!   makehmat      hddata.f90
!   DSYEVR        LAPACK
SUBROUTINE getEnergy(cgeom,energy)
  USE progdata, ONLY : natoms,atoms
  USE potdata
  USE hddata, ONLY:  nstates,ncoord,makehmat
  IMPLICIT NONE
  DOUBLE PRECISION,dimension(3*natoms),intent(IN)                  ::  cgeom
  DOUBLE PRECISION,dimension(nstates,nstates)                      ::  hmat
  DOUBLE PRECISION,dimension(nstates),intent(OUT)                  ::  energy

  double precision,dimension(nstates,nstates)       ::  evec
  double precision,dimension(ncoord)                ::  igeom
  double precision,dimension(ncoord,3*natoms)       ::  bmat
  integer                                           :: m,LWORK,LIWORK,INFO
  integer,dimension(nstates*2)                      :: ISUPPZ
  double precision,dimension(nstates*(nstates+26))  :: WORK
  integer,dimension(nstates*10)                     :: IWORK
  double precision    :: bohr2ang,  mind
  integer   :: i,j,ptid
  integer   :: counter = 1   ! count the number of evaluations
  character(4)  ::  str,str2

  LWORK  = nstates*(nstates+26)
  LIWORK = nstates*10

  if(.not.initialized)then
       print *,"WARNING:  Initializing Hd for pot().  "
       call prepot()
  end if
  call buildWBMat(cgeom,igeom,bmat,.false.)
 ! construct Hd
  call makehmat(igeom,hmat)
 ! generate eigenvectors and energies at current geometry
  CALL DSYEVR('V','A','U',nstates,hmat,nstates,dble(0),dble(0),0,0,1D-12,m,&
            energy,evec,nstates,ISUPPZ,WORK,LWORK,IWORK,LIWORK, INFO )

END SUBROUTINE getEnergy

!-----------------------------------------------------------------------------------
SUBROUTINE pauseParsing()
  USE potdata,only: parsing,paused
  paused = parsing .or. paused
  parsing = .false.
END SUBROUTINE
SUBROUTINE resumeParsing()
  USE potdata,only: parsing,paused
  parsing = parsing .or. paused
  paused = .false.
END SUBROUTINE
SUBROUTINE enableParsing()
  USE potdata, ONLY: parsing
  parsing = .true.
END SUBROUTINE
SUBROUTINE disableParsing()
  USE potdata, ONLY: parsing
  parsing = .false.
END SUBROUTINE
!-----------------------------------------------------------------------------------
! read general input from input file (surfgen.in) for prepot() 
! specific inputs for potential library mode are located in namelist POTLIB
SUBROUTINE readginput(jtype) 
  use potdata
  use progdata 
  use hddata, only: initGrps,ncoord,nstates,order,getFLUnit,CpOrder 
  use CNPI, only: irrep,GrpPrty,GrpSym,nSymLineUps 
  IMPLICIT NONE 
  INTEGER,INTENT(INOUT)           :: jtype 
  INTEGER                         :: i,j,mkadiabat,ios,k 
  INTEGER                         :: nGrp,jobtype 
  INTEGER,dimension(10)           :: surface,updatehess 
  INTEGER,dimension(50)           :: ci,atmgrp 
  DOUBLE PRECISION,dimension(500) :: minguess,mexguess 
  INTEGER,DIMENSION(20,MAX_ALLOWED_SYM)        :: groupSym,groupPrty 
  DOUBLE PRECISION,dimension(50)  :: e_guess
   
  NAMELIST /GENERAL/        jobtype,natoms,order,nGrp,groupsym,groupprty,&
                            printlvl,inputfl,atmgrp,nSymLineUps,cntfl,CpOrder
  NAMELIST /POTLIB/         molden_p,m_start,switchdiab,calcmind,gflname,nrpts, cvanish,&
                            mindcutoff, atomlabels,dcoordls,errflname, ndcoord,&
                            timeeval,parsing,eshift,calcErr,nfterms,fterm,fcoef,forig

  atomlabels = ''

  do i=1,MaxNDCoord
    dcoordls(i) = i
  end do
  ndcoord = 3
  nSymLineUps = 1
  CpOrder    = -1
  jtype      = 0 
  natoms     = 2 
  printlvl   = 1 
  cvanish    = 0d0
  inputfl    = ''
  eshift     = dble(0) 
  switchdiab = .false.
  
  print *,"Entering readinput()." 
 !----------- GENERAL SECTION ----------------! 
  open(unit=INPUTFILE,file='surfgen.in',access='sequential',form='formatted',& 
       IOSTAT=ios,POSITION='REWIND',ACTION='READ',STATUS='OLD') 
  if(ios/=0)then 
    print *,"readinput: cannot open file surfgen.in.  IOSTAT=",ios 
  end if!ios/=0 
 
  read(unit=INPUTFILE,NML=GENERAL) 
  if(printlvl>0)print *,"    readinput():  Control parameters read from surfgen.in" 
  jtype  = jobtype 
  if(CpOrder<0) CpOrder = order
 
  call genAtomList(atmgrp) 
   
  call readIrreps() 
  if(allocated(GrpSym))deallocate(GrpSym) 
  if(allocated(GrpPrty))deallocate(GrpPrty) 
  allocate(GrpSym(nGrp,nSymLineUps)) 
  allocate(GrpPrty(nGrp,nSymLineUps)) 
  GrpSym=groupsym(1:nGrp,1:nSymLineUps) 
  GrpPrty=groupprty(1:nGrp,1:nSymLineUps) 
  call initGrps(nGrp,irrep(GrpSym(:,1))%Dim) 
 
  if(jobtype.ne.0)print *,"WARNING:  Calling prepot() with jobtype.ne.0" 
  print *,"   Reading POTLIB related parameters"
  molden_p     = 0
  m_start      = 100
  nfterms = 0
  calcmind     = .false.
  parsing      = .false.
  calcErr      = .false.
  gflname      = 'refgeom'
  errflname    = 'error.log'
  nrpts        = 20
  timeeval     = .false.
  mindcutoff   = 1D-5
  nfterms      = 0
  fcoef        = 0d0
  forig        = 1d0
  read(unit=INPUTFILE,NML=POTLIB)
  print 1000,m_start,molden_p
  close(unit=INPUTFILE) 
  if(printlvl>0)print *,"Exiting readinput..." 
  return 
1000 format("Point to start MOLDEN recording ",I4,"   Record every ",I4," Points")
end SUBROUTINE readginput 

! Read in the list of reference geometries used to construct Hd
SUBROUTINE loadRefGeoms()
  USE progdata, ONLY: natoms
  USE potdata
  USE hddata, ONLY: ncoord,getFLUnit,nstates
  IMPLICIT NONE
  integer              :: i, j, uerrfl, ios
  double precision     :: cgeom(3*natoms,nrpts),igeom(ncoord),bmat(ncoord,3*natoms)
  character(3),dimension(natoms)               :: atoms
  double precision,dimension(natoms)           :: anums,masses
  character(3)    ::  ncstr
  integer         :: fid
  PRINT *,"   Minimum distances to a set of reference geometries will be calculated."
  allocate(rgeom(nrpts,ncoord))
  PRINT *,"   Loading ",nrpts," geometries from file ",trim(adjustl(gflname))
  fid = getFLUnit()
  open(unit=fid,file=trim(adjustl(gflname)),access='sequential',form='formatted',&
      status='old',action='read',position='rewind',iostat=ios)
  if(ios/=0)stop "failed to open geometry input file"
  read(unit=fid,fmt=*,iostat=ios)cgeom
  if(ios/=0)stop "failed to read geometries" 
  close(unit=fid)
  print *,"   ",nrpts," geometries loaded."
  PRINT *,"   Converting to internal coordinates"
  do i=1,nrpts
    call buildWBMat(cgeom(:,i),igeom,bmat,.false.)
    rgeom(i,:)=igeom
  end do
 ! Load point error information if availale
  if(calcErr)then
    print *,"Energy error will be estimated at evaluation points"
    uerrfl = getFLUnit()
    open(unit=uerrfl,file=trim(adjustl(errflname)),access='sequential',  &
       form='formatted',IOSTAT=ios,POSITION='REWIND',ACTION='READ',STATUS='OLD')
    if(ios/=0)then
      print *,"loadRefGeoms: cannot open error log file.  IOSTAT=",ios
      stop "program terminated"
    end if!ios/=0
    allocate(enererrdata(nrpts,nstates))
    allocate(graderrdata(nrpts,nstates,ncoord))
    write(ncstr,"(I3)")ncoord
    do i=1,nrpts
      read(uerrfl,*) enererrdata(i,:)
      do j=1,nstates
        read(uerrfl,*) graderrdata(i,j,:)
      end do
    end do
    close(uerrfl)
  end if!calcErr
END SUBROUTINE

! This subroutine is used to set time and surface index data for output
SUBROUTINE setTrajData(ttime, tisurf)
  USE potdata, ONLY:  timetraj, isurftraj
  IMPLICIT NONE
  DOUBLE PRECISION,intent(IN)  :: ttime
  INTEGER, intent(IN)          :: tisurf
  
  timetraj  = ttime
  isurftraj = tisurf
END SUBROUTINE setTrajData

SUBROUTINE getInfo(nat,nst)
  use hddata, only:   nstates
  use progdata, only: natoms
  integer,intent(out)  :: nat,nst

  nat = natoms
  nst = nstates
END SUBROUTINE
