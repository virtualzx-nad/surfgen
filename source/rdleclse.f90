!This module contains multiple RDLECLSE(Rank-Deficient Linear
!Equility Constrained Least Squares Equations) solver.
!The memory management subroutine are also contained within this module
!Required memory for the operation is automatically allocated and managed and
!is invisible from outside unless a fatal problem occurs.
MODULE rdleclse
  implicit none
  private
  !temporary arrays
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE   :: WORK,TAU
  INTEGER,DIMENSION(:),ALLOCATABLE            :: IWORK,JPVT
  INTEGER               :: LWORK=0,LTAU=0,LIWORK=0,LJPVT=0

  public :: solve,cleanArrays

CONTAINS

  !--------------------------------------------------------------
  ! Normal Equations with Lagrange Multipliers Solved by EVD
  !--------------------------------------------------------------
  ! minimize |A.x-y|**2 with restrictions B.x=z
  ! Here the restrictions are enforced by Lagrange multipliers
  ! Total equation becomes
  ! ( 0   B     )( l ) = ( z )        ...       (1)
  ! ( B' A'A+tI )( x )   (A'y)
  ! Undetermined unknowns are driven to zero by term pI
  ! Note that the left hand side matrix is not positive definite.
  ! Eigenvalue decomposition is performed for the lhs
  ! ( 0   B     )= M = U.SIGMA.U^T
  ! ( B' A'A+tI )
  ! Lagrange Multipliers and solutions are obtained by
  ! ( l ) = U.SIGMA^-1.U^T ( z )
  ! ( x )                  (A'y)
  !--------------------------------------------------------------
  ! Arguments
  ! m        (input) INTEGER
  !          The number of unknown coefficients.
  ! nlse    (input) INTEGER
  !          Number of least squares equations
  ! nex      (input) INTEGER
  !          Number of exact equations, which are independent
  ! A        (input) DOUBLE PRECISION array,dimension(r_ex+nlse,m)
  !          Coefficient matrices.  The top r_ex rows contain the
  !          exact equation block B and bottom nlse rows contain
  !          least squares block A
  ! rhs      (input) DOUBLE PRECISION array,dimension(r_ex+nlse)
  !          Right hand side vector of the equation.  First r_ex
  !          elements comes from exact equations.
  ! tol_ex   (input) DOUBLE PRECISION
  !          Tolerance of accuracy for exact equations.  Used as
  !          criteria for singularity when calculating the inverse
  !          of diagonal matrix.
  ! t        (input) DOUBLE PRECISION
  !          Flattening factor added to the diagonal elements of
  !          normal equations block.  This term drives undetermined
  !          unknown coefficients to zero.
  ! sol      (output)DOUBLE PRECISION array,dimension(m)
  !          Solution of the constrained least squares equations.
  ! printlvl (input) INTEGER
  !          Level of output to default I/O.
  SUBROUTINE solve(m,nlse,nex,NEL,rhs,tol_ex,t,sol,printlvl)
    IMPLICIT NONE
    INTEGER,VALUE,INTENT(IN)      :: m,nex,nlse,printlvl
    DOUBLE PRECISION,INTENT(INOUT):: NEL(nex+m,nex+m)
    DOUBLE PRECISION,INTENT(INOUT):: rhs(nex+m)
    DOUBLE PRECISION,INTENT(IN)   :: tol_ex,t
    DOUBLE PRECISION,INTENT(OUT)  :: sol(m)
    !temporary arrays for Langrange Multiplier Normal Equations
    !evec       :  Eigenvectors of NEL
    !eval       :  Eigenvalues of NEL
    double precision,dimension(:,:),allocatable     :: evec
    double precision,dimension(:)  ,allocatable     :: eval
    integer,dimension(:),allocatable                :: ISUPPZ
integer,dimension(:),allocatable                :: ipiv
    integer  :: nsstart,nsend  !starting and ending index of null space
    integer  :: i,n,nx,rank, INFO
    double precision,external :: dnrm2

    if(printlvl>0)print 1000,m,nlse,nex
    n  = nex+nlse
    nx = nex+m

    allocate(ipiv(nx),STAT=INFO)
    IF(INFO/=0) STOP "solve: failed to allocate memory for ipiv"

    do i=nex+1,nx
      NEL(i,i)=NEL(i,i)+t
    end do

    CALL allocArrays(1,0,0,0)
    call dsysv('U',nx,1,NEL,nx,ipiv,rhs,nx,work,-1,info)
    if(info/=0)stop "dsysv: space query failed"
    
    CALL allocArrays(int(WORK(1)),0,0,0)
    call dsysv('U',nx,1,NEL,nx,ipiv,rhs,nx,work,lwork,info)
    if(info/=0)stop "dsysv: failed to solve linear equations"

    sol = rhs(nex+1:nx)
  1000 format(4X,'Solving Linear Equility Constrained Least Squares Equations',/,&
              6X,'Number of Unkown:',I7,', LSE:',I7,', Exact Eq:',I4)
  END SUBROUTINE solve

!----------------------------------------------------------------
! internal subroutines used to handle memory allocations
!----------------------------------------------------------------
  ! allocArrays sets the workspace for the module scratch arrays
  ! WORK, IWORK, JPVT and TAU
  ! It compares the required volume of array with existing length.
  ! When needed, new workspace is allocated and size increased
  ! A nonpositive size will be ignored
  ! SYNTAX
  ! CALL allocArrays(l_work,l_iwork,l_jpvt,l_tau)
  subroutine allocArrays(l_work,l_iwork,l_jpvt,l_tau)
    implicit none
    INTEGER,VALUE,INTENT(IN):: l_work,l_iwork,l_jpvt,l_tau
    integer :: status
    status=0
    !Allocating WORK
    if(l_work>LWORK)then
      if(allocated(WORK))deallocate(WORK,STAT=status)
      if(status/=0)then
        print *,"STAT=",status
        stop "rdleclse.allocArrays: cannot deallocate WORK."
      end if!(status/=0)
      allocate(WORK(l_work),STAT=status)
      if(status/=0)then
        print *,"STAT=",status
        stop "rdleclse.allocArrays: failed to allocate memory for WORK."
      end if!(status/=0)
      LWORK=l_work
    end if!(L_WORK>LWORK.AND.L_WORK>0)
    !Allocating IWORK
    if(l_iwork>LIWORK)then
      if(allocated(IWORK))deallocate(IWORK,STAT=status)
      if(status/=0)then
        print *,"STAT=",status
        stop "rdleclse.allocArrays: cannot deallocate IWORK."
      end if!(status/=0)
      allocate(IWORK(l_iwork),STAT=status)
      if(status/=0)then
        print *,"STAT=",status
        stop "rdleclse.allocArrays: failed to allocate memory for IWORK."
      end if!(status/=0)
      LIWORK=l_iwork
    end if!(L_IWORK>LIWORK.AND.L_IWORK>0)
    !Allocating JPVT
    if(l_jpvt>LJPVT)then
      if(allocated(JPVT))deallocate(JPVT,STAT=status)
      if(status/=0)then
        print *,"STAT=",status
        stop "rdleclse.allocArrays: cannot deallocate JPVT."
      end if!(status/=0)
      allocate(JPVT(l_jpvt),STAT=status)
      if(status/=0)then
        print *,"STAT=",status
        stop "rdleclse.allocArrays: failed to allocate memory for JPVT."
      end if!(status/=0)
      LJPVT=l_jpvt
    end if!(L_JPVT>LJPVT.AND.L_JPVT>0)
    if(allocated(JPVT))JPVT=0
    !Allocating TAU
    if(l_TAU>LTAU)then
      if(allocated(TAU))deallocate(TAU,STAT=status)
      if(status/=0)then
        print *,"STAT=",status
        stop "rdleclse.allocArrays: cannot deallocate TAU."
      end if!(status/=0)
      allocate(TAU(l_TAU),STAT=status)
      if(status/=0)then
        print *,"STAT=",status
        stop "rdleclse.allocArrays: failed to allocate memory for TAU."
      end if!(status/=0)
      LTAU=l_TAU
    end if!(L_TAU>LTAU.AND.L_TAU>0)
  end subroutine allocArrays

  !----------------------------------------------------
  ! cleanArrays release all memories allocated to scratch space
  ! and set all upper bounds to zero
  subroutine cleanArrays()
    implicit none
    integer    ::  status
    status=0
    if(allocated(WORK))deallocate(WORK,STAT=status)
    if(status/=0)print *,"rdleclse.cleanArrays:  cannot to deallocate WORK"
    if(allocated(IWORK))deallocate(IWORK,STAT=status)
    if(status/=0)print *,"rdleclse.cleanArrays:  cannot to deallocate IWORK"
    if(allocated(JPVT))deallocate(JPVT,STAT=status)
    if(status/=0)print *,"rdleclse.cleanArrays:  cannot to deallocate JPVT"
    if(allocated(TAU))deallocate(TAU,STAT=status)
    if(status/=0)print *,"rdleclse.cleanArrays:  cannot to deallocate TAU"
    LWORK=0
    LIWORK=0
    LJPVT=0
    LTAU=0
  end subroutine cleanArrays
end module rdleclse
