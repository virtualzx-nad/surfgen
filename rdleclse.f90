!This module contains multiple RDLECLSE(Rank-Deficient Linear
!Equility Constrained Least Squares Equations) solvers.
!All solvers are interfaced through one
MODULE rdleclse
  implicit none
  private
  !temporary arrays
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE   :: WORK,TAU
  INTEGER,DIMENSION(:),ALLOCATABLE            :: IWORK,JPVT
  INTEGER               :: LWORK=0,LTAU=0,LIWORK=0,LJPVT=0
  DOUBLE PRECISION,DIMENSION(1)               :: TEST
  INTEGER                                     :: INFO

  public :: solve,cleanArrays

CONTAINS

  subroutine solve(method,m,nlse,nex,NEL,rhs,&
          tol_ex,tol_lse,sol,printlvl,INFO)
  !--------------------------------------------------------------
  !This subroutine examines the dimensionalities of LSE and exact
  !blocks and calls the corresponding subroutine depending on the
  !value of method.  Working procedure is as follows
  !--------------------------------------------------------------
  !1. Rank revealing QR decompositions for exact block to remove
  !inconsistencies among exact equations
  !2. If Normal Equations approach is not selected, perform rank
  !revealing QR decomposition for (L' E') with E' block as lead
  !columns to remove underdetermined unknowns
  !3. Calls the subroutines with new L and E blocks
  !4. Transform the solution back to old basis if neccesary
  !--------------------------------------------------------------
  ! Method Keywords:
  ! NEL-SVD     Normal equations with Lagrange Multipliers solved
  !             by Singular Value Decomposition
  ! NEL-EVD     Normal equations with Lagrange Multipliers solved
  !             by Eigenvalue Decomposition with RRR algorithm
    implicit none
    CHARACTER(32),INTENT(IN)       :: method
    DOUBLE PRECISION,INTENT(IN)    :: tol_lse,tol_ex
    INTEGER,INTENT(IN)             :: m,nlse,nex,printlvl
    DOUBLE PRECISION,INTENT(INOUT) :: NEL(nex+m,nex+m)
    DOUBLE PRECISION,INTENT(INOUT) :: rhs(nex+m)
    DOUBLE PRECISION,INTENT(OUT)   :: sol(m)
    INTEGER         ,INTENT(OUT)   :: INFO
    double precision, external     :: dnrm2
    integer                        :: rank
    integer         :: i, lc(nex)
    character(32)   :: meth
    double precision:: sv(nex)
    if(printlvl>0)then
      PRINT 1000
      PRINT 1001,m,nlse,nex
    end if!(printlvl)
    meth=strUpCase(32,method)
    if(meth/='NEL-EVD'.and. nex==0)then
      Print 1100
      meth=trim(adjustl('NEL-EVD'))
    end if!(meth/='NE'.and. nex==0)
    select case(meth)
      case('NEL-SVD')
      !Lagrange Multipliers for E and Normal equations A'Ax=A'y for L
      !solved by SVD
        CALL NELSVD(m,nex,nlse,NEL,rhs,tol_ex,&
                      tol_lse,rank,sol,printlvl,INFO)
        if(INFO/=0)then
          print *,'solve:  NELSVD Failed.'
          return
        end if!(INFO/=0)
      case default !NE, or anything else
      !Lagrange Multipliers for E and Normal equations A'Ax=A'y for L
      !solved by EVD with RRR algorithm
        if(meth/='NEL-EVD')print 1004
        meth='NEL-EVD'
        CALL NELEVD(m,nex,nlse,NEL,rhs,tol_ex,&
                      tol_lse,rank,sol,printlvl,INFO)
        if(INFO/=0)then
          print *,'solve:  NELEVD Failed.'
          return
        end if!(INFO/=0)
    end select!case(method)
    if(printlvl>0)then
      print 1998,nex+m-rank,DNRM2(m,sol,int(1))/sqrt(dble(m))
    end if!(printlvl>0)
  1000 format(4X,'Solving Linear Equility Constrained Least Squares Equations')
  1100 format(6X,'NEX==0, Switched to normal equations approach.')
  1001 format(6X,'Number of Unkown:',I7,', LSE:',I5,', Exact Eq:',I4)
  1003 format(6X,'Exact Equations are Linearly Independant')
  1004 format(6X,'Method not recognized. Using normal equations with EVD.')
  1008 format(8X,"Exact equations are linear dependant")
  1998 format(6X,'Undetermined Unknowns:',I5,4X,&
                 'RMS Sol Vector:',F7.3)
  1999 format(6X,A,'Rank=',I5,', RMS Error=',E8.2)
  2000 format(6X,'No LSE present.  Exact equations will be solved by DGELSD via SVD.')
  2001 format(6X,'Rank of Exact Equations: ',I5,', RMS Error=',E8.2,', Condition Number ',E8.2)
  end subroutine solve

!----------------------------------------------------------------
! these are private subroutines each solve the problem using a
! different algothrism
!----------------------------------------------------------------
  ! Normal Equations with Lagrange Multipliers Solved by SVD
  !--------------------------------------------------------------
  !   minimize |A.x-y|**2 with restrictions B.x=z
  ! Here the restrictions are enforced by Lagrange multipliers
  ! Total equation becomes
  ! ( 0   B     )( l ) = ( z )        ...       (1)
  ! (-B' A'A+tI )( x )   (A'y)
  ! Undetermined unknowns are driven to zero by term pI
  ! Note that the left hand side matrix is not positive definite.
  ! singular value decomposition is performed for the lhs
  ! ( 0   B     )= M = U.SIGMA.V^T
  ! ( B' A'A+tI )
  ! Lagrange Multipliers and solutions are obtained by
  ! ( l ) = V.SIGMA^-1.U^T ( z )
  ! ( x )                  (A'y)
  !--------------------------------------------------------------
  ! THE EVD SUBROUTINE RUNS FASTER AND USE LESS SPACE IN ALL CASES
  ! THIS SUBROUTINE IS KEPT ONLY FOR CONSISTENCY WITH MATHEMATICA 
  ! VERSION OF THE SCRIPT  
  !--------------------------------------------------------------
  ! Arguments
  ! m        (input) INTEGER
  !          The number of unknown coefficients.
  ! nex      (input) INTEGER
  !          Number of exact equations, which are independent
  ! nlse    (input) INTEGER
  !          Number of least squares equations
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
  ! rank     (output)INTEGER
  !          Rank of total matrix (Lagrange multipliers and NE)
  ! sol      (output)DOUBLE PRECISION array,dimension(m)
  !          Solution of the constrained least squares equations.
  ! printlvl (input) INTEGER
  !          Level of output to default I/O.
  ! INFO     (output)INTEGER
  !          Error code.  0 if finished successfully.
  SUBROUTINE NELSVD(m,nex,nlse,NEL,rhs,tol_ex,t,rank,sol,printlvl,INFO)
    IMPLICIT NONE
    INTEGER,VALUE,INTENT(IN)      :: m,nex,nlse,printlvl
    INTEGER,INTENT(OUT)           :: rank,INFO
    DOUBLE PRECISION,INTENT(INOUT):: NEL(nex+m,nex+m)
    DOUBLE PRECISION,INTENT(INOUT):: rhs(nex+m)
    DOUBLE PRECISION,INTENT(IN)   :: tol_ex,t
    DOUBLE PRECISION,INTENT(OUT)  :: sol(m)
    !temporary arrays for Langrange Multiplier Normal Equations
    !NEL        :  Normal equations with least squares
    !rhsL       :  Right hand side which includes Lagrange multipliers
    !U,VT       :  Left and right singular vectors of NEL
    !sval       :  Singular values of NEL
    double precision,dimension(:,:),allocatable     :: U,VT
    double precision,dimension(:)  ,allocatable     :: sval
    integer  :: i,n,nx
    if(printlvl>0)print 1000
    n=nex+nlse
    nx=nex+m

    !perform diagonal shift for LSE block
    do i=nex+1,nx
      NEL(i,i)=NEL(i,i)+t
    end do

    !Allocate space for Eigenvalue Decomposition results
    allocate(sval(nx),STAT=INFO)
    IF(INFO==0)allocate(U(nx,nx),STAT=INFO)
    IF(INFO==0)allocate(VT(nx,nx),STAT=INFO)
    IF(INFO/=0)THEN
      PRINT *,"solve: failed to allocate memory for SVD temp matrices"
      RETURN
    END IF
    ! Perform SV Decomposition
    CALL DGESVD('A','A',nx,nx,NEL,nx,sval,U,nx,VT,nx,TEST,int(-1),INFO)
    IF(INFO/=0)THEN
      PRINT *,"solve:  DGESVD workspace query failed."
      RETURN
    END IF
    CALL allocArrays(int(TEST(1)),0,0,0)
    CALL DGESVD('A','A',nx,nx,NEL,nx,sval,U,nx,VT,nx,WORK,LWORK,INFO)
    IF(INFO/=0)THEN
      PRINT *,"solve:  DGESVD failed for NEL-SVD."
      RETURN
    END IF
    ! WORK=U**T.y'
    CALL allocArrays(nx,0,0,0)
    CALL DGEMV('T',nx,nx,dble(1),U,nx,rhs,int(1),&
              dble(0),WORK,int(1))
    ! WORK=S**-1*WORK
    rank=0
    do i=1,nx
      if(abs(sval(i))<tol_ex)then
         WORK(i)=dble(0)
      else!if(abs(sval(i))<tol_ex)
         WORK(i)=WORK(i)/sval(i)
         rank=rank+1
      end if!(abs(sval(i))<tol_ex)
    end do!i=1,n
    ! x=V.WORK
    CALL DGEMV('T',nx,nx,dble(1),Vt,nx,WORK,int(1),&
              dble(0),rhs,int(1))
    sol=rhs(nex+1:nx) 
    ! Output values of Lagrange Multipliers
    if(printlvl>1.and. nex>0)then
      print 1001
      print 1002,rhs(1:nex)
    end if!(printlvl>1.and. nex>0)
    deallocate(sval)
    deallocate(U)
    deallocate(Vt)
  1000 format(6X,'Solving Normal Equations with Lagrange Multipliers.',/&
              8x,'Equations are solved by SVD')
  1001 format(6X,'Value of Lagrange Multipliers:')
  1002 format(6X,7F9.2)
  1003 format(6X,"Memory required for intermediate matrices:",F8.2,"Mb")
  END SUBROUTINE NELSVD

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
  ! See NELSVD
  SUBROUTINE NELEVD(m,nex,nlse,NEL,rhs,tol_ex,t,rank,sol,printlvl,INFO)
    IMPLICIT NONE
    INTEGER,VALUE,INTENT(IN)      :: m,nex,nlse,printlvl
    INTEGER,INTENT(OUT)           :: rank,INFO
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
    integer  :: nsstart,nsend  !starting and ending index of null space
    integer  :: i,n,nx

    if(printlvl>0)print 1000
    n  = nex+nlse
    nx = nex+m

    allocate(ISUPPZ(nx*2),STAT=INFO)
    IF(INFO/=0)THEN
      PRINT *,"solve: failed to allocate memory for ISUPPZ"
      RETURN
    END IF

    do i=nex+1,nx
      NEL(i,i)=NEL(i,i)+t
    end do

    !Allocate space for Eigenvalue Decomposition results
    allocate(eval(nx),STAT=INFO)
    IF(INFO==0)allocate(evec(nx,nx),STAT=INFO)
    IF(INFO/=0)THEN
     PRINT *,"solve: failed to allocate memory for EVD matrices"
     RETURN
    END IF

    ! Perform symmetric eigenvalue decomposition
    CALL allocArrays(1,1,0,0)
    CALL DSYEVR('V','A','U',nx,NEL,nx,0.,0.,0,0,tol_ex,i,eval,&
          evec,nx,ISUPPZ,WORK,int(-1),IWORK,int(-1),INFO)
    IF(INFO/=0)THEN
      PRINT *,"solve:  DSYEVR workspace query failed."
      RETURN
    END IF
    CALL allocArrays(int(WORK(1)),IWORK(1),0,0)
    PRINT 1004,(LWORK*dble(2)+LIWORK)/262144.
    CALL DSYEVR('V','A','U',nx,NEL,nx,0.,0.,0,0,tol_ex,i,eval,&
          evec,nx,ISUPPZ,WORK,LWORK,IWORK,LIWORK,INFO)
    IF(INFO/=0)THEN
      PRINT *,"solve:  DSYEVR failed for NELEVD."
      RETURN
    END IF
    ! WORK=U**T.y'
    CALL allocArrays(nx,0,0,0)
    CALL DGEMV('T',nx,nx,dble(1),evec,nx,rhs,int(1),&
              dble(0),WORK,int(1))
    ! WORK=S**-1*WORK
    rank=0
    nsstart = nex+1
    nsend = nx
if(printlvl>3)print "(15E10.2)",eval
    do i=1,nx
      if(abs(eval(i))<tol_ex)then
         WORK(i)=dble(0)
      else!if(abs(eval(i))<tol_ex)
         WORK(i)=WORK(i)/eval(i)
         rank=rank+1
      end if!(abs(eval(i))<tol_ex)
      if(i>nex)then
        if(eval(i)>t+tol_ex.and.eval(i-1)<t+tol_ex)nsend=i-1
        if(eval(i)>t-tol_ex.and.eval(i-1)<t-tol_ex)nsstart=i
      end if
    end do!i=1,n
    ! x=U.WORK
    CALL DGEMV('N',nx,nx,dble(1),evec,nx,WORK,int(1),&
              dble(0),rhs,int(1))
    sol=rhs(nex+1:nx) 
    ! Output values of Lagrange Multipliers
    if(printlvl>1.and. nex>0)then
      print 1001
      print 1002,rhs(1:nex)
    end if!(printlvl>1.and. nex>0)
    NEL(1:nsend-nsstart+1,1:m) = evec(nex+1:nx,nsstart:nsend)
    deallocate(evec)
    if(printlvl>2)then
      print *,"Performing analysis of the null space"
      print *," Null space range:",nsstart," to ",nsend
      !print *,"  Eigenvalues of normal equations matrix"
      !print "('  ',10E14.5)",eval
      if(printlvl>4)then
        print *,"  Performing QR-decomposition of null space"
        JPVT = 0
        CALL allocArrays(4*nx,0,nx,nx)
        call DGEQP3(nsend-nsstart+1,m,NEL,nx,jpvt,tau,WORK,LWORK,INFO)
PRINT *,"dependent part: jpvt="
PRINT "(15I5)",JPVT(1:nsend-nsstart+1)
PRINT *,"independent part: jpvt="
PRINT "(15I5)",JPVT(nsend-nsstart+2:m)
PRINT *,"UPPER LEFT CORNER OF OUT MAT:"
PRINT "(10E14.5)",NEL(1:10,1:10)
open(unit=2999,file='coeflist',access='sequential',form='formatted',&
   status='replace',action='write',position='rewind',iostat=INFO)
if(INFO/=0)print *,"PROBLEM OPENING COEFLIST OUTPUT.  IOS=",INFO
WRITE(2999,"(I20)")JPVT(nsend-nsstart+2:m)
close(2999)
      end if
    end if
    deallocate(eval)
    deallocate(ISUPPZ)
  1000 format(6X,'Solving Normal Equations with Lagrange Multipliers.',/&
              8x,'Equations are solved by EVD')
  1001 format(6X,'Value of Lagrange Multipliers:')
  1002 format(6X,7F9.2)
  1003 format(6X,"Memory required for intermediate matrices:",F8.2,"MB")
  1004 format(6X,"Memory required for symmetric EVD:    ",F12.2,"MB")
  END SUBROUTINE NELEVD

subroutine dummytest(a,lda,nnull)
  integer :: lda,nnull
  double precision :: a(lda,nnull)

  print *,"INDUMMY:",A(7,6)
  stop
end subroutine dummytest
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
    JPVT=0
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

  !----------------------------------------------------
  ! Rank Revealing QR
  !----------------------------------------------------
  subroutine revealRank(M,N,A,NLC,LC,piv,tol,Q,rank,printlvl,INFO)
    IMPLICIT NONE
    LOGICAL,VALUE,   INTENT(IN)   :: piv
    INTEGER,VALUE,   INTENT(IN)   :: M,N,NLC,printlvl
    INTEGER,         INTENT(OUT)  :: INFO,rank
    INTEGER,         INTENT(IN)   :: LC(NLC)
    DOUBLE PRECISION,INTENT(IN)   :: tol
    DOUBLE PRECISION,INTENT(INOUT):: A(M,N)
    DOUBLE PRECISION,INTENT(OUT)  :: Q(M,M)
    integer :: i,j
    if(printlvl>1)print 1000
    CALL allocArrays(0,0,N,min(N,M))
    CALL DGEQP3(M,N,A,M,JPVT,TAU,TEST,int(-1),INFO)
    IF(INFO/=0)THEN
      print *,"revealRank: workspace query failed for DGEQP3."
      return
    END IF!INFO/=0
    CALL allocArrays(int(TEST(1)),0,0,0)
    if(NLC>0)JPVT(LC)=1  !Set leading columns
    CALL DGEQP3(M,N,A,M,JPVT,TAU,WORK,LWORK,INFO)
    IF(INFO/=0)THEN
      print *,"revealRank: cannot perform QRE decomposition."
      return
    END IF!INFO/=0
    !determine the effective rank
    rank=NLC
    do i=NLC+1,min(M,N)
      if(Abs(A(i,i))>tol)then
        rank=i
      else!if(exblk(i,i)>tol_ex)
        cycle
      end if!(exblk(i,i)>tol_ex)
    end do!i=1,min(M,N)
    !generate matrix orthogonal Q
    Q(:,1:min(M,N))=A(:,1:min(M,N))
    Q(:,N+1:m)=dble(0)
    !replace entries below the diagonal with 0
    do i=1,rank-1
      A(i+1:rank,i)=dble(0)
    end do
    !pivot back if pivot==.false.
    if(.not. piv)then
       do i=1,n
         do while(JPVT(i)/=i)
           j=JPVT(i)
           CALL DSWAP(rank,A(1:rank,j),int(1),A(1:rank,i),int(1))
           JPVT(i)=JPVT(j)
           JPVT(j)=j
         end do!while(JPVT(j)/=i)
       end do
    end if!.not. pivot
    CALL DORGQR(M,min(M,N),min(M,N),Q,M,TAU,TEST,int(-1),INFO)
    if(INFO/=0)then
      print *,"revealRank: workspace query failed for DORGQR"
      RETURN
    end if
    CALL allocArrays(int(TEST(1)),0,0,0)
    CALL DORGQR(M,min(M,N),min(M,N),Q,M,TAU,WORK,LWORK,INFO)
    if(INFO/=0)then
      print *,"revealRank: DORGQR failure."
      RETURN
    end if

  1000 format(6x,'Performing Ranking Revealing QR')
  end subroutine revealRank

!--------------------------------
!converts lower case characters to upper case
  function strUpCase(l,str1) result(str2)
    implicit none
    integer,intent(in)       :: l
    character(l),intent(in)  :: str1
    character(l)             :: str2
    integer :: i
    str2=trim(adjustl(str1))
    do i=1,len_trim(str2)
      if(iachar(str2(i:i))>96.and. iachar(str2(i:i))<123)&
          str2(i:i)=achar(iachar(str2(i:i))-32)
    end do
    return
  end function strUpCase
end module rdleclse
