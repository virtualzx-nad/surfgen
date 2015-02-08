! This file contains the parallel wrapper.
! use c preprocessor macro variable PARALLEL to enable compilation with
! ScaLapack and MPI.   When the variable is not defined, normal Lapack
! subroutines are used and message passing between nodes are removed.

! psurfgen contains information regarding parallel execution, as well as
! local portion of global matrices.   In serial runs, the array
! themselves are stored here.
MODULE psurfgen
  IMPLICIT NONE
  ! number of processors, rows and columns of the process grid
  INTEGER            ::  nproc, nprow, npcol  
  ! rank and coordinate of the current process
  INTEGER            ::  myrow, mycol, rank
  ! context index
  INTEGER            ::  context
  ! dimensionality of blocks.  number of columns in block will be set to the
  ! same as number of variables + lagrange multipliers, and therefore should not
  ! be used. 
  INTEGER            ::  rpb, cpb

! local portion of global arrays.  these arrays are distributed onto each
! process are uniquely stored with ScaLapack using 2D-block cyclic distrubution.

  ! local portion of normal equations matrix
  DOUBLE PRECISION,dimension(:,:),allocatable  ::  NEL
  ! local portion of raw equation matrix 
  DOUBLE PRECISION,dimension(:,:),allocatable  ::  AMat
  ! local portion of rhs vector
  DOUBLE PRECISION,dimension(:),allocatable    ::  rhs
  ! local portion of raw rhs vector
  DOUBLE PRECISION,dimension(:),allocatable    ::  bvec

! data regarding point partition
  !  Number of points on this process
  INTEGER         ::  nlocalPt
  !  Max number of points per process
  INTEGER         ::  PpP
  !  Range of local points
  INTEGER         ::  LocalPtLB,LocalPtUB
 
! Point specific information that are partitioned into each process.  Each
! process contain data for a set of data points.   
  

! shared information that are stored redundantly on each processes.

END MODULE psurfgen
!------------------------------------------------------------
! Partition points to processes
SUBROUTINE PartitionPt(npoints)
  implicit none
  include "mpif.h"
  integer,intent(IN)  :: npoints
  integer  ::  ierr
!if multi-process, sync number of points
#ifdef PARALLE
  Call MPI_Bcast(npoints,1, MPI_INTEGER,0,&
           MPI_COMM_WORLD,ierr)
  if(ierr.ne.0)then
    print *,"Failed to broadcast/receive npoints at rank",rank
    call blacs_abort(context,ierr)  
  end if
#endif
  ! calculate number of points per process
  PpP=ceiling(dble(npoints)/nproc)
  ! get point range
  LocalPtLB    =  rank*PpP+1
  LocalPtUB    =  (rank+1)*PpP
  if(LocalPtUB>npoints)LocalPtUB=npoints
  ! calculate number of points allocated to this process
  nLocalPt = LocalPtUB-LocalPtLB+1
END SUBROUTINE PartitionPt 
!------------------------------------------------------------
! wrapper for MPI/BLACS barrier function
SUBROUTINE pbarrier()
  IMPLICIT NONE
#ifdef PARALLEL
  call blacs_barrier(context,'A')
#endif  
END SUBROUTINE pbarrier
!------------------------------------------------------------
! Initialize the program and determine if the program is compiled as
! parallel version.  
SUBROUTINE init_psurfgen()
  implicit none
  integer  ::  nr,nc

! initialize parallel environment
#ifdef PARALLEL
  ! get number of processes and 
  call BLACS_pinfo(rank,nproc)
  ! one dimensional grid is used
  nprow = nproc
  npcol = 1
  rpb   = 4
  ! initialize grid
  call sl_init(context,nprow,npcol)
  ! get my coordinates 
  call blacs_gridinfo(context,nr,nc,myrow,mycol)
#else
  nproc = 1
  nprow = 1 
  npcol = 1
  myrow = 0
  mycol = 0
  rank  = 0
  context = -1
  rpb   = 1
#endif
END SUBROUTINE init_psurfgen
!------------------------------------------------------------
! Initialize local and global arrays
SUBROUTINE initArrays 
END SUBROUTINE 
!------------------------------------------------------------
! Synchronize important data 
SUBROUTINE syncPSurfgen()
END SUBROUTINE syncPSurfgen
