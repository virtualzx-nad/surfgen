! Partition states at a certain geometry into groups that are energetically close.
! Group states into minimal sets so that the energy difference between any two states
! from different groups is greater than de. Input energies are supposed to be sorted.
! Energy groups are only generated for the range state1 to state2
SUBROUTINE genEnerGroups(pt,de)
  use progdata, only: abpoint
  IMPLICIT NONE
  type(abpoint),INTENT(INOUT)     :: pt
  DOUBLE PRECISION,INTENT(IN)     :: de

  integer                         :: i,ngrp,cnt
  integer,dimension(pt%ub,2)      :: deg_table

  ngrp=0
  cnt=1
  do i=pt%lb+1,pt%ub
    if(pt%energy(i,i)-pt%energy(i-1,i-1)>de)then
      if(cnt>1)then
        ngrp=ngrp+1
        deg_table(ngrp,2)=cnt
        deg_table(ngrp,1)=i-cnt
      end if
      cnt=0
    end if
    cnt=cnt+1
  end do
  if(cnt>1)then
    ngrp=ngrp+1
    deg_table(ngrp,2)=cnt
    deg_table(ngrp,1)=pt%ub-cnt+1
  end if

  if(allocated(pt%deg_groups))deallocate(pt%deg_groups)
  allocate(pt%deg_groups(ngrp,2))
  pt%ndeggrp=ngrp
  pt%deg_groups=deg_table(1:ngrp,:)
END SUBROUTINE

! Obtain the rotation angle beta that orthogonalize g and h
SUBROUTINE orthgh(n,g,h,beta)
  IMPLICIT NONE
  INTEGER,intent(IN)          :: n
  DOUBLE PRECISION,intent(IN) :: g(n),h(n)
  DOUBLE PRECISION,intent(OUT):: beta
  double precision :: denom,num
  num=dot_product(g,h)*2
  denom=-dot_product(g,g)+dot_product(h,h)
  if(abs(denom)<1D-10)then
PRINT *,"DENORM VERY SMALL.  NO ROTATION IS PERFORMED"
      beta=0D0
  else
      beta=atan(num/denom)/4
  end if
END SUBROUTINE orthgh

!---------------------------------------------------------------------
! Chop() sets the value of the variable to 0 if it is smaller
! than a threshold
SUBROUTINE chop(dpvar,threshold)
 implicit none

 double precision, intent(INOUT)    ::  dpvar
 double precision, intent(in)       ::  threshold

 if(abs(dpvar)<abs(threshold))dpvar=0

END SUBROUTINE chop


!---sort the elements of an array-------------------------
SUBROUTINE sort(nelements,data,ordering)
 implicit none
 integer,intent(IN)                         ::  nelements
 integer,dimension(nelements),intent(INOUT) ::  data
 integer,intent(IN)                         ::  ordering

 integer      :: i,j,temp

 do i=nelements-1,1,-1
  do j=i,nelements-1
    if(data(j)*ordering>data(j+1)*ordering)then
      temp=data(j)
      data(j)=data(j+1)
      data(j+1)=temp
    end if
  end do
 end do

 return
END SUBROUTINE sort

!-------------------------------------------------------------------
! SCHDMIT() Performs Schdmit orthorgonalization for rows of matrix A
!-------------------------------------------------------------------
! Each row in A will be orthogonalized to all row vectors above it.
! After removal, if the norm of residue vector is smaller than $(tol),
! it will be permuted to the bottom of A and considered 0. Otherwise,
! The residue vector will be normalized
!-------------------------------------------------------------------
! m,n       (input) INTEGER
!           Number of rows/columns
! A         (input/output) DOUBLE PRECISION,dimension(m,n)
!           On entry, rows of A contain vectors to be orthogonalized
!           On exit, first $(rank) rows of A are orthonormal and the
!           other rows are zero
! d         (input) INTEGER
!           Dimension of vector.  Only first $(d) columns will be 
!           used when calculating inner products and included in 
!           orthogonalization.  Other columns will simply be carried.  
! tol       (input) DOUBLE PRECISION
!           Tolerance for the norm of residue to be considered zero.
! rank      (output) INTEGER
!           Number of orthonormal vectors contained in A
! nstart    (input) INTEGER
!           The number of row to start orthogonalization.  The rows
!           above nstart will be supposed to be orthogonal and will
!           only be normalized.
SUBROUTINE Schdmit(m,n,A,d,tol,rank,nstart)
  IMPLICIT NONE
  INTEGER,intent(IN)                           :: m,n,d,nstart
  DOUBLE PRECISION,dimension(m,n),intent(INOUT):: A
  DOUBLE PRECISION,intent(IN)                  :: tol
  INTEGER,intent(OUT)                          :: rank
  integer           :: i
  double precision  :: p(m),nrm
  double precision,external :: DNRM2
  rank=m
  do i=1,nstart-1
    A(i,:)=A(i,:)/DNRM2(d,A(i,1:d),int(1))
  end do
  i=nstart
  do while(i<=rank)
    !Compute the inner product with all row vectors above
    if(i>1)then
      CALL DGEMV('N',i-1,d,dble(1),A,m,A(i,1),m,dble(0),p,int(1))
      CALL DGEMV('T',i-1,n,dble(-1),A,m,p,int(1),dble(1),A(i,1),m)
    end if!(i>1)
    nrm=DNRM2(d,A(i,1),m)
    if(nrm<tol)then
      CALL DSWAP(n,A(i,1),m,A(rank,1),m)
      rank=rank-1
    else!if(nrm<tol)
      A(i,:)=A(i,:)/nrm
      i=i+1
    end if!(nrm<tol)
  end do!while(i<m)
END SUBROUTINE Schdmit

!DSYINV inverts a symmetric double precision matrix
SUBROUTINE DSYINV(UPLO,N,A,LDA,INVA,LDI,TOL,NNEG,NZERO,ENFPD)
  IMPLICIT NONE
  CHARACTER(1),intent(IN)     :: UPLO
  INTEGER,intent(IN)          :: N,LDA,LDI
  DOUBLE PRECISION,intent(IN) :: A(LDA,N),TOL
  DOUBLE PRECISION,intent(OUT):: INVA(LDI,N)
  LOGICAL,intent(IN)          :: ENFPD
  INTEGER,INTENT(OUT)         :: NNEG,NZERO
  double precision,dimension(:),allocatable     :: WORK
  integer,dimension(:),allocatable              :: IWORK
  integer          :: LWORK,LIWORK,ISUPPZ(2*N),INFO,itst(1),nev,j
  double precision :: dtst(1),eval(N),evec(N,N)

  CALL DSYEVR('V','A',UPLO,n,A,LDA,dble(0),dble(0),0,0,TOL,nev, &
         eval,evec,n,ISUPPZ,dtst,int(-1),itst,int(-1), INFO )
  if(INFO/=0)STOP"DSYINV: Workspace query failed."
  LWORK  = int(dtst(1))
  LIWORK = itst(1)
  allocate(WORK(LWORK))
  allocate(IWORK(LIWORK))
  CALL DSYEVR('V','A',UPLO,n,A,LDA,dble(0),dble(0),0,0,TOL,nev, &
             eval,evec,n,ISUPPZ,WORK,LWORK,IWORK,LIWORK, INFO )
  if(INFO/=0)STOP"DSYINV: eigenvalue decomposition failed."
  nneg=0
  nzero=n-nev
  invA=dble(0)
  do j=1,nev
    if(abs(eval(j))>1D-4)then
      if(eval(j)<0)nneg=nneg+1
      if(eval(j)<0.and. ENFPD)eval(j)=dble(1)
      CALL DSYR(UPLO,n,1/eval(j),evec(:,j),int(1),invA,LDI)
    else
      nzero=nzero+1
    end if
  end do
  deallocate(WORK)
  deallocate(IWORK)
END SUBROUTINE

SUBROUTINE init_random_seed()
  INTEGER :: i, n, clock
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed
          
  CALL RANDOM_SEED(size = n)
  ALLOCATE(seed(n))
          
  CALL SYSTEM_CLOCK(COUNT=clock)
          
  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  CALL RANDOM_SEED(PUT = seed)
          
  DEALLOCATE(seed)
END SUBROUTINE
! --------------------------------------------------------------
![Description]
! getCartHd evaluates Hd at a given Cartesian coordinate, returning adiabatic
! energies, energy gradients, derivative couplings, diabatic Hamiltonian matrix
! and gradients of each of the blocks of Hd.   Vectors are given in the same
! Cartesian coordinates.
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
!   buildWBmat    libinternal.f90
!   EvalRawTerms  hddata.f90
!   EvalHdDirect  hddata.f90
!   DSYEVR        LAPACK
SUBROUTINE getCartHd(cgeom,energy,cgrads,hmat,dcgrads)
  USE progdata, ONLY : natoms
  USE hddata, ONLY:  nstates,ncoord,EvalRawTermsL,EvalHdDirect
  IMPLICIT NONE
  DOUBLE PRECISION,dimension(3*natoms),intent(IN)                  ::  cgeom
  DOUBLE PRECISION,dimension(nstates,nstates),intent(OUT)          ::  hmat
  DOUBLE PRECISION,dimension(3*natoms,nstates,nstates),intent(OUT) ::  cgrads,dcgrads
  DOUBLE PRECISION,dimension(nstates),intent(OUT)                  ::  energy

  double precision,dimension(nstates,nstates)              ::  evec,hmat2
  double precision,dimension(ncoord,nstates,nstates)       ::  dhmat
  double precision,dimension(ncoord)                ::  igeom
  double precision,dimension(ncoord,3*natoms)       ::  bmat 
  integer                                           :: m,LWORK,LIWORK,INFO
  integer,dimension(nstates*2)                      :: ISUPPZ
  double precision,dimension(nstates*(nstates+26))  :: WORK
  integer,dimension(nstates*10)                     :: IWORK
  integer   :: i,j

  LWORK  = nstates*(nstates+26)
  LIWORK = nstates*10

  call buildWBMat(cgeom,igeom,bmat,.false.)

! calculate raw polynomial terms
  call EvalRawTermsL(igeom)

 ! construct Hd and its derivatives
  call EvalHdDirect(hmat,dhmat)

  ! convert gradients and couplings to cartesian coordinates
  dcgrads=dble(0) 
  do i=1,3*natoms
    do j=1,ncoord
      dcgrads(i,:,:)=dcgrads(i,:,:)+ dhmat(j,:,:)*bmat(j,i)
    end do !j=1,ncoord
  end do !i=1,3*natoms

 ! generate eigenvectors and energies at current geometry
  hmat2 = hmat
  CALL DSYEVR('V','A','U',nstates,hmat2,nstates,dble(0),dble(0),0,0,1D-16,m,&
            energy,evec,nstates,ISUPPZ,WORK,LWORK,IWORK,LIWORK, INFO )
  do i=1,3*natoms
    cgrads(i,:,:)=matmul(transpose(evec),matmul(dcgrads(i,:,:),evec))
  end do!i=1,ncoord

END SUBROUTINE getCartHd

!------------------------------------------------------------------------------
! Rotate degenerate ab initio adiabatic data so that g and h vectors between any
! pairs of two degenerate states are orthogonalized.  The rotation is made so that
! total norm of off diagonal gradients (h vectors) is minimized.
!pt     : ab initio data of the point to perform transformation 
!maxiter: maximum number of jacobi rotations
!toler  : convergence tolerance for rotation angle beta
!hasGrad: specifies if gradient/couple data is available for a block
![method]
! This subroutine tries orthogonalize all g-h vectors by an algorithm similar to 
! the Jacobi's method for eigenvalues.  The Jacobi transformations here extremize 
! the norm of the corresponding coupling block rather than eliminating them.   
SUBROUTINE OrthGH_ab(pt,maxiter,toler,hasGrad)
  USE hddata, only: nstates
  USE progdata, only: abpoint,printlvl
  IMPLICIT NONE
  TYPE(abpoint),INTENT(INOUT)                           :: pt
  INTEGER,INTENT(IN)                                    :: maxiter
  LOGICAL,DIMENSION(nstates,nstates),INTENT(in)         :: hasGrad
  DOUBLE PRECISION,INTENT(IN)                           :: toler

  integer           ::  igrp,i,j,iter,ldeg,udeg
  integer           ::  mi,mj  !location of maximum rotation
  double precision, dimension(pt%nvibs,nstates,nstates) :: newgrad
  double precision, dimension(nstates,nstates)          :: newener 
  double precision, dimension(nstates,nstates)          :: beta   !required rotation 
  double precision           :: max_b,t, c,s
  ! allowedRot stores the infomation whether rotation between two specific states 
  ! will not cause a gradient data to mix into a gradient that has not data
  LOGICAL,dimension(nstates,nstates)  :: allowedRot

  beta=dble(0)  
  if(pt%ndeggrp>0.and.printlvl>1)print *,"     Transforming ab initio data for point ",pt%id
  if(printlvl>2.and. pt%ndeggrp>0)then
    print *,"      Couplings and gradients before transformation" 
    do i=1,nstates
      do j=1,i
        print 1002,i,j
        if(hasGrad(i,j))then
            print 1003,pt%grads(:pt%nvibs,i,j)
        else
            print "(10X,A)","Data not available"
        end if
      end do
    end do!i=1,nstates
  end if!(printlvl>2)
  do igrp=1,pt%ndeggrp  ! in case there are multiple groups of degeneracies, loops over all of them       
    max_b=-1
    ldeg=pt%deg_groups(igrp,1)
    udeg=pt%deg_groups(igrp,2)+pt%deg_groups(igrp,1)-1
    iter=0
    ! check if rotations among these states are allowed.
    if(all(hasGrad))then
      allowedRot(ldeg:udeg,ldeg:udeg)=.true.
    else
      allowedRot(ldeg:udeg,ldeg:udeg)=.false.
      if(printlvl>2)write (*,"(6X,A)",advance='no') "Allowed rotations:"
      do i=ldeg,udeg
          do j=ldeg,i-1
            ! both gradients and coupling between then has to be present for them to be rotatable
            if(hasGrad(i,i).and.hasGrad(i,j).and.hasGrad(j,j))then
              ! the available gradients list has to be identical for them
              if(all(hasGrad(i,:).eqv.hasGrad(j,:)))then
                  if(printlvl>2)write (*,"(' (',I1,',',I1,') ')",advance='no') I,J
                  allowedRot(i,j)=.true.
                  allowedRot(j,i)=.true.
              end if
            end if!hasGrad ii, ij, jj
          end do!j
      end do!i=ldeg,udeg
      if(.not.any(allowedRot(ldeg:udeg,ldeg:udeg)))then
        print *," none"
        print "(4X,A)","Missing data forbit any rotations.  Skipping transformation of degenerate group."
        cycle
      else
        print *,""
      end if
    end if
    !build normalized g and h vectors from Hd and compute rotation angles the first time
    do i=ldeg,udeg 
      do j=ldeg,i-1
        if(.not.allowedRot(i,j))then
          beta(i,j) = 0d0
          cycle
        end if
        beta(i,j) = getBeta(i,j)
        if(abs(beta(i,j))>max_b)then
          max_b=abs(beta(i,j))
          mi=i
          mj=j
        end if!(abs(beta(i,j))>max_b)
      end do!j=ldeg,i-1
    end do!i=ldeg,udeg
    if(printlvl>2)print "(6X,A,E12.5,2(A,I1))","max(|beta|)=",max_b,&
                " ,block:",mi,",",mj
    do while(iter<maxiter.and.max_b>toler)
      iter=iter+1
      t=beta(mi,mj)
      c = cos(t)
      s = sin(t)
! Gnew = J^T.G.J.   Gnew_ij=Jki*Gkl*Jlj
      newgrad = pt%grads(:pt%nvibs,:,:)
      newener = pt%energy
      do i=1,nstates
        newgrad(:,i,mi)=pt%grads(:pt%nvibs,i,mi)*c-pt%grads(:pt%nvibs,i,mj)*s
        newgrad(:,i,mj)=pt%grads(:pt%nvibs,i,mj)*c+pt%grads(:pt%nvibs,i,mi)*s
        newener(i,mi)=pt%energy(i,mi)*c-pt%energy(i,mj)*s
        newener(i,mj)=pt%energy(i,mj)*c+pt%energy(i,mi)*s
      end do
      pt%grads(:pt%nvibs,:,:) = newgrad
      pt%energy               = newener
      do i=1,nstates
        pt%grads(:pt%nvibs,mi,i)=newgrad(:,mi,i)*c-newgrad(:,mj,i)*s
        pt%grads(:pt%nvibs,mj,i)=newgrad(:,mj,i)*c+newgrad(:,mi,i)*s 
        pt%energy(mi,i)=newener(mi,i)*c-newener(mj,i)*s
        pt%energy(mj,i)=newener(mj,i)*c+newener(mi,i)*s 
      end do
      !update rotation angles
      do i=mj+1,udeg
        if(allowedRot(i,mj))beta(i,mj)=getBeta(i,mj)
      end do
      do j=ldeg,mj-1
        if(allowedRot(mi,j))beta(mi,j)=getBeta(mi,j)
      end do
      do j=mj+1,mi-1
        if(allowedRot(mi,j))beta(mi,j)=getBeta(mi,j)
      end do
      max_b=-1
      do i=ldeg,udeg
        do j=ldeg,i-1
          if(.not.allowedRot(i,j))cycle
          if(abs(beta(i,j))>max_b)then
            max_b=abs(beta(i,j))
            mi=i
            mj=j
          end if!(abs(beta(i,j))>max_b)
        end do!j=ldeg,i-1
      end do!i=ldeg,udeg      
      if(printlvl>2)print *,"   iteration ",iter,", max(|beta|)=",max_b
    end do!while(iter<maxiter.and.max_b>toler)do
    if(printlvl>1.and.iter>0)then
        if(max_b<toler)then
          print 1001,iter
        else
          print 1000,iter,max_b
        end if!(max_b<toler)
        print *,"      Energies after transformation" 
        do i=1,nstates
          print "(20F12.7)",pt%energy(i,:)
        end do
        print *,"      Couplings and gradients after transformation" 
        do i=1,nstates
          do j=1,i
            write (*,1002) i,j
            if(hasGrad(i,j))then
                print 1003,pt%grads(:pt%nvibs,i,j)
            else
                print "(10X,A)","Data not available"
            end if
          end do
        end do!i=1,nstates
    end if!(printlvl>0)
  end do!igrp=1,pt%ndeggrp
1000 format(7X,"no convergence after ",I4," iterations, max residue angle=",F8.2)
1001 format(7X,"convergence after ",I4," iterations")
1002 format(14X,'states(',I2,',',I2,'):')
1003 format(12X,3E16.7)
CONTAINS
  FUNCTION getBeta(i,j) RESULT(beta)
    USE progdata, only : abpoint,pi
    IMPLICIT NONE
    INTEGER,INTENT(IN)         :: i,j
    DOUBLE PRECISION           :: beta

    double precision, dimension(pt%nvibs)                :: g,h

    g=(pt%grads(:pt%nvibs,i,i)-pt%grads(:pt%nvibs,j,j))/2
    h=pt%grads(:pt%nvibs,i,j)
    CALL orthgh(pt%nvibs,g,h,beta)
  END FUNCTION getBeta
END SUBROUTINE OrthGH_Ab 

!------------------------------------------------------------------------------
! Rotate Hd generated eigenvectors so that the resulting g and h vectors between any
! pairs of two degenerate states are orthogonalized.  g and h vectors are identified 
! by their overlap with the ab initio g and h vectors
! pt     : ab initio data of the point to perform transformation, already orthongalized 
! dhmat  : gradient of diabatic Hamiltonian matrix
! ckl    : eigenvectors corresponding to rotation that orthogonalize g and h vectors 
! maxiter: maximum number of jacobi rotations
! toler  : convergence tolerance for rotation angle beta
![method]
! This subroutine tries orthogonalize all g-h vectors by an algorithm similar to 
! the Jacobi's method for eigenvalues.  The Jacobi transformations here extremize 
! the norm of the corresponding coupling block rather than eliminating them.   
SUBROUTINE OrthGH_Hd(pt,dhmat,ckl,maxiter,toler,hasGrad)
  USE hddata, only: nstates
  USE progdata, only: abpoint,printlvl
  IMPLICIT NONE
  TYPE(abpoint),INTENT(IN)                                        :: pt
  INTEGER,INTENT(IN)                                              :: maxiter
  DOUBLE PRECISION,INTENT(IN)                                     :: toler
  LOGICAL,DIMENSION(nstates,nstates),INTENT(in)         :: hasGrad
  DOUBLE PRECISION,dimension(pt%nvibs,nstates,nstates),INTENT(IN) :: dhmat
  DOUBLE PRECISION,dimension(nstates,nstates),INTENT(INOUT)       :: ckl

  integer           ::  igrp,i,j,iter,ldeg,udeg
  integer           ::  mi,mj  !location of maximum rotation
  double precision, dimension(nstates,nstates)         :: jacobi    !wave function rotation matrix
  double precision, dimension(nstates,nstates)         :: beta      !desired rotation angles
  double precision, dimension(pt%nvibs,nstates,nstates):: gAb, hAb,&!ab initio g and h vectors
                                                          gradnew   !Hd gradients after rotation
  double precision           :: max_b,t
  ! allowedRot stores the infomation whether rotation between two specific states 
  ! will not cause a gradient data to mix into a gradient that has not data
  LOGICAL,dimension(nstates,nstates)  :: allowedRot
 
  if(pt%ndeggrp<=0)return 
! initialize eigenvectors and rotated gradients
  do i=1,pt%nvibs
    gradnew(i,:,:)=matmul(matmul(transpose(ckl),dhmat(i,:,:)),ckl)
  end do
  beta=dble(0)  
  if(printlvl>2)print *,"     transforming Hd eigenvectors for point ",pt%id
  if(printlvl>2)then
    print *,"      reference ab initio gradients" 
    do i=1,nstates
      do j=1,i
        print 1002,i,j
        if(hasGrad(i,j))then
            print 1003,pt%grads(:pt%nvibs,i,j)
        else
            print "(10X,A)","Data not available"
        end if
      end do
    end do!i=1,nstates
    print *,"      gradients before transformation" 
    do i=1,nstates
      do j=1,i
        print 1002,i,j
        print 1003,dhmat(:,i,j)
      end do
    end do!i=1,nstates
  end if!(printlvl>2)
! tabulate ab initio g and h vectors  
  do i=1,nstates
    do j=1,i-1
      gAb(:,i,j) = (pt%grads(1:pt%nvibs,i,i)-pt%grads(1:pt%nvibs,j,j))/2
      hAb(:,i,j) =  pt%grads(1:pt%nvibs,i,j)
      gAb(:,j,i) =-gAb(:,i,j)
      hAb(:,j,i) = hAb(:,i,j)
    end do
  end do!i
! generate rotations
  do igrp=1,pt%ndeggrp  ! in case there are multiple groups of degeneracies, loops over all of them       
    max_b=-1
    ldeg=pt%deg_groups(igrp,1)
    udeg=pt%deg_groups(igrp,2)+pt%deg_groups(igrp,1)-1
    if(printlvl>2)print "(3(A,I3))","  Degenerate group ",igrp,", states ",ldeg," to ",udeg
    iter=0

    ! check if rotations among these states are allowed.
    if(all(hasGrad))then
      allowedRot(ldeg:udeg,ldeg:udeg)=.true.
    else
      allowedRot(ldeg:udeg,ldeg:udeg)=.false.
      if(printlvl>2)write (*,"(6X,A)",advance='no') "Allowed rotations:"
      do i=ldeg,udeg
          do j=ldeg,i-1
            ! both gradients and coupling between then has to be present for
            ! them to be rotatable
            if(hasGrad(i,i).and.hasGrad(i,j).and.hasGrad(j,j))then
              ! the available gradients list has to be identical for them
              if(all(hasGrad(i,:).eqv.hasGrad(j,:)))then
                  if(printlvl>2)write (*,"(' (',I1,',',I1,') ')",advance='no') I,J
                  allowedRot(i,j)=.true.
                  allowedRot(j,i)=.true.
              end if
            end if!hasGrad ii, ij, jj
          end do!j
      end do!i=ldeg,udeg
      if(.not.any(allowedRot(ldeg:udeg,ldeg:udeg)))then
        if(printlvl>2)then
           print *," none"
           print "(4X,A)","Missing data forbit any rotations.  Skipping transformation of degenerate group."
        end if
        cycle
      else
        if(printlvl>2)print *,""
      end if
    end if

    !build normalized g and h vectors from Hd and compute rotation angles the first time
    do i=ldeg,udeg 
      do j=ldeg,i-1
        beta(i,j) = getBetaHd(i,j)
        if(abs(beta(i,j))>max_b)then
          max_b=abs(beta(i,j))
          mi=i
          mj=j
        end if!(abs(beta(i,j))>max_b)
      end do!j=ldeg,i-1
    end do!i=ldeg,udeg
    if(printlvl>2)print "(6x,A,E15.7,A,I3,A,I3,A)","max(|beta|)=",max_b,&
                        "(",mi,",",mj,")"
    do while(iter<maxiter.and.max_b>toler)
      iter=iter+1
      t=beta(mi,mj)
      jacobi=dble(0)  !construct Givens rotation matrix
      do i=1,nstates
        jacobi(i,i)=dble(1)
      end do
      jacobi(mi,mi)=cos(t)
      jacobi(mj,mj)=jacobi(mi,mi)
      jacobi(mi,mj)=sin(t)
      jacobi(mj,mi)=-jacobi(mi,mj)
      do i=1,pt%nvibs
        gradnew(i,:,:)=matmul(matmul(transpose(jacobi),gradnew(i,:,:)),jacobi)
      end do
      if(printlvl>2)then
        print *,"  Ab initio g vector:"
        print "(20F7.3)",(pt%grads(:pt%nvibs,mi,mi)-pt%grads(:pt%nvibs,mj,mj))/2
        print *,"  Fit g vector:"
        print "(20F7.3)",(gradnew(:pt%nvibs,mi,mi)-gradnew(:pt%nvibs,mj,mj))/2
        print *,"  Ab initio h vector:"
        print "(20F7.3)",pt%grads(:pt%nvibs,mi,mj)
        print *,"  Fit h vector:"
        print "(20F7.3)",gradnew(:pt%nvibs,mi,mj)
      end if
      ckl = matmul(ckl,jacobi)
      !update rotation angles
      do i=mj+1,udeg 
        beta(i,mj)=getBetaHd(i,mj)
      end do
      do j=ldeg,mj-1
        beta(mi,j)=getBetaHd(mi,j)
      end do
      do j=mj+1,mi-1
        beta(mi,j)=getBetaHd(mi,j)
      end do
      max_b=-1
      do i=ldeg,udeg
        do j=ldeg,i-1
          if(abs(beta(i,j))>max_b)then
            max_b=abs(beta(i,j))
            mi=i
            mj=j
          end if!(abs(beta(i,j))>max_b)
        end do!j=ldeg,i-1
      end do!i=ldeg,udeg      
      if(printlvl>2)print "(4x,A,I4,A,E15.7,A,I3,A,I3,A)","iteration ",iter,", max(|beta|)=",max_b,&
                        "(",mi,",",mj,")"
    end do!while(iter<maxiter.and.max_b>toler)do
    if(iter==0.and.printlvl>2)then
        print *,"  Ab initio g vector:"
        print "(30F7.3)",(pt%grads(:pt%nvibs,mi,mi)-pt%grads(:pt%nvibs,mj,mj))/2
        print *,"  Fit g vector:"
        print "(30F7.3)",(gradnew(:pt%nvibs,mi,mi)-gradnew(:pt%nvibs,mj,mj))/2
        print *,"  Ab initio h vector:"
        print "(30F7.3)",pt%grads(:pt%nvibs,mi,mj)
        print *,"  Fit h vector:"
        print "(30F7.3)",gradnew(:pt%nvibs,mi,mj)
    end if
    if(printlvl>2.and.iter>0)then
        if(max_b<toler)then
          print 1001,iter
        else
          print 1000,iter,max_b
        end if!(max_b<toler)
        print *,"      gradients after transformation" 
        do i=1,nstates
          do j=1,i
            write (*,1002) i,j
            write (*,1003) gradnew(:pt%nvibs,i,j)
          end do
        end do!i=1,nstates
    end if!(printlvl>0)
  end do!igrp=1,pt%ndeggrp
1000 format(7X,"no convergence after ",I4," iterations, max residue angle=",F8.2)
1001 format(7X,"convergence after ",I4," iterations")
1002 format(14X,'states(',I2,',',I2,'):')
1003 format(12X,3E20.11)
CONTAINS
  FUNCTION getBetaHd(i,j) RESULT(beta)
    USE progdata, only : abpoint,pi
    IMPLICIT NONE
    INTEGER,INTENT(IN)         :: i,j
    DOUBLE PRECISION           :: beta

    double precision, dimension(pt%nvibs)                :: g,h
    double precision  ::  errg, errh, errcurr, errrot

    if(.not.allowedRot(i,j))then
      beta = 0d0
      return
    end if
    g=(gradnew(:,i,i)-gradnew(:,j,j))/2
    h=gradnew(:,i,j)
    CALL orthgh(pt%nvibs,g,h,beta)
    errg = min(dot_product(g-gAb(:,i,j),g-gAb(:,i,j)),dot_product(g+gAb(:,i,j),g+gAb(:,i,j)))
    errh = min(dot_product(h-hAb(:,i,j),h-hAb(:,i,j)),dot_product(h+hAb(:,i,j),h+hAb(:,i,j)))
    errcurr=0
    if(hasGrad(i,j)) errcurr=errh
    if(hasGrad(i,i).and.hasGrad(j,j)) errcurr=errcurr+errg
    errrot =0
    errg = min(dot_product(g-hAb(:,i,j),g-hAb(:,i,j)),dot_product(g+hAb(:,i,j),g+hAb(:,i,j)))
    errh = min(dot_product(h-gAb(:,i,j),h-gAb(:,i,j)),dot_product(h+gAb(:,i,j),h+gAb(:,i,j)))
    if(hasGrad(i,j)) errrot=errg
    if(hasGrad(i,i).and.hasGrad(j,j)) errrot=errrot+errh
    if(errrot<errcurr)then
      if(printlvl>2)print *,"    G and H vectors are switched."
      if(beta>0)then
          beta=beta-pi/4
      else
          beta=beta+pi/4
      end if
    end if
  END FUNCTION getBetaHd
END SUBROUTINE 
!------------------------------------------------------------------------------
! Rotate degenerate ab initio adiabatic data so that g and h vectors between any
! pairs of two degenerate states are orthogonalized.  The rotation is made so that
! total norm of off diagonal gradients (h vectors) is minimized.
! This subroutine is used as interface for external programs. 
!dh             : gradients in the current representation
!h              : hamiltonian in the current representation
!lstate,ustate  : range of degenerate states
!nvibs          : vibrational degrees of freedom
!
!maxiter: maximum number of jacobi rotations
!toler  : convergence tolerance for rotation angle beta
!hasGrad: specifies if gradient/couple data is available for a block
![method]
! This subroutine tries orthogonalize all g-h vectors by an algorithm similar to 
! the Jacobi's method for eigenvalues.  The Jacobi transformations here extremize 
! the norm of the corresponding coupling block rather than eliminating them.   
SUBROUTINE OrthogonalizeGH(h,dh,lstate,ustate,nvibs,maxiter,toler)
  USE hddata, only: nstates
  USE progdata, only: abpoint,printlvl
  IMPLICIT NONE
  DOUBLE PRECISION,INTENT(INOUT),DIMENSION(nvibs,nstates,nstates) :: dh
  DOUBLE PRECISION,INTENT(INOUT),DIMENSION(nstates,nstates)       :: h
  INTEGER,INTENT(IN)                                    :: maxiter,lstate,ustate,nvibs
  DOUBLE PRECISION,INTENT(IN)                           :: toler

  integer           ::  i,j,iter
  integer           ::  mi,mj  !location of maximum rotation
  double precision, dimension(nvibs,nstates,nstates)    :: newgrad
  double precision, dimension(nstates,nstates)          :: newener 
  double precision, dimension(nstates,nstates)          :: beta   !required rotation 
  double precision           :: max_b,t, c,s

  beta=dble(0)  
  max_b=-1
  iter=0
  !build normalized g and h vectors from Hd and compute rotation angles the first time
  do i=lstate,ustate 
      do j=lstate,i-1
        beta(i,j) = getBeta(i,j)
        if(abs(beta(i,j))>max_b)then
          max_b=abs(beta(i,j))
          mi=i
          mj=j
        end if!(abs(beta(i,j))>max_b)
      end do!j=lstate,i-1
    end do!i=lstate,ustate
    do while(iter<maxiter.and.max_b>toler)
      iter=iter+1
      t=beta(mi,mj)
      c = cos(t)
      s = sin(t)
! Gnew = J^T.G.J.   Gnew_ij=Jki*Gkl*Jlj
      newgrad = dh 
      newener = h
      do i=1,nstates
        newgrad(:,i,mi)=dh(:,i,mi)*c-dh(:,i,mj)*s
        newgrad(:,i,mj)=dh(:,i,mj)*c+dh(:,i,mi)*s
        newener(i,mi)=h(i,mi)*c-h(i,mj)*s
        newener(i,mj)=h(i,mj)*c+h(i,mi)*s
      end do
      dh     = newgrad
      h      = newener
      do i=1,nstates
        dh(:,mi,i)=newgrad(:,mi,i)*c-newgrad(:,mj,i)*s
        dh(:,mj,i)=newgrad(:,mj,i)*c+newgrad(:,mi,i)*s 
        h(mi,i)=newener(mi,i)*c-newener(mj,i)*s
        h(mj,i)=newener(mj,i)*c+newener(mi,i)*s 
      end do
      !update rotation angles
      do i=mj+1,ustate
        beta(i,mj)=getBeta(i,mj)
      end do
      do j=lstate,mj-1
        beta(mi,j)=getBeta(mi,j)
      end do
      do j=mj+1,mi-1
        beta(mi,j)=getBeta(mi,j)
      end do
      max_b=-1
      do i=lstate,ustate
        do j=lstate,i-1
          if(abs(beta(i,j))>max_b)then
            max_b=abs(beta(i,j))
            mi=i
            mj=j
          end if!(abs(beta(i,j))>max_b)
        end do!j=lstate,i-1
      end do!i=lstate,ustate 
    end do!while(iter<maxiter.and.max_b>toler)do
CONTAINS
  FUNCTION getBeta(i,j) RESULT(beta)
    USE progdata, only : abpoint,pi
    IMPLICIT NONE
    INTEGER,INTENT(IN)         :: i,j
    DOUBLE PRECISION           :: beta

    double precision, dimension(nvibs)                :: g,h

    g=(dh(:,i,i)-dh(:,j,j))/2
    h=dh(:,i,j)
    CALL orthgh(nvibs,g,h,beta)
  END FUNCTION getBeta
END SUBROUTINE OrthogonalizeGH 


!---print out geometry info, including bond lenths, angles and torsions---
! Arguments:
! INTEGER			[in] natoms
!				Number of atoms in the molecule
! DOUBLE PRECISION(3*natoms)	[in] geom
!				geometry input
! CHARACTER*3(3*natoms)		[in] aname
!				name of elements
! DOUBLE PRECISION(3*natoms)	[in] anum
!				atomic number of each atom
! DOUBLE PRECISION(3*natoms)	[in] masses
!				Nuclear mass of each atom
! DOUBLE PRECISION		[in] TLen
!				Threshold for bond-length
! LOGICAL			[in] ShortList
!				Whether to generate only a short list of coordinates, 
!				instead those from all connected atoms
subroutine analysegeom(natoms,geom,aname,anum,masses,TLen,ShortList)
  implicit none
  integer, intent(in)          ::  natoms
  character*3,intent(in)       ::  aname(natoms)
  double precision,intent(in)  ::  anum(natoms),masses(natoms)
  double precision,intent(in)  ::  geom(3,natoms)
  double precision,intent(in)  ::  TLen
  logical,intent(in)           ::  ShortList
  double precision, parameter  ::  bohr2ang=0.529177249d0
  integer   ::  i,j,k,l
  double precision  ::  distmat(natoms,natoms),  d,d1(3),d2(3),d3(3),cpd(3)
  double precision,external  ::   dnrm2
  logical :: hasOOP(natoms)
  print *,"Cartesian Geometries in Atomic Units"
  do i=1,natoms
     print "(x,a3,1x,f4.1,3F14.8,F14.8)",aname(i),anum(i),geom(:,i),masses(i)
  end do
  distmat = 0d0
  print "(/,A)","   Atom1   Atom2   Bond Length(Ang)"
  do i=1,natoms
    do j=i+1,natoms
      d = dnrm2(3,geom(:,i)-geom(:,j),1)*bohr2ang
      distmat(i,j) = d
      distmat(j,i) = d
      if(d<TLen)print "(2x,I5,3x,I5,3x,F12.8)",i,j,d
    end do
  end do

  print "(/,A)","   Atom1   Atom2   Atom3   Bond Angle (Degrees)"
  do i=1,natoms
    do j=1,natoms
     if(j==i .or. distmat(i,j)>TLen)cycle
     d1 = geom(:,j)-geom(:,i)
     d1 = d1/dnrm2(3,d1,1)
     do k=j+1,natoms
       if(k==i .or. distmat(i,k)>TLen)cycle
       d2 = geom(:,k)-geom(:,i)
       d2 = d2/dnrm2(3,d2,1)
       print "(2x,3(I5,3x),F12.4)",J,I,K, 90/Acos(0d0)* &
           ACOS(dot_product(d1,d2))
     end do!k
    end do !j
  end do   !i    

  hasOOP = .false.
  print "(/,A)","   Atom1   Atom2   Atom3   Atom4   Torsion Angle (Degrees)"
  do i=1,natoms
    do j=1,natoms
     if(j==i)cycle
     d1 = geom(:,j)-geom(:,i)
     d1 = d1/dnrm2(3,d1,1)
     do k=j+1,natoms
       if(k==i .or. distmat(j,k)>TLen .or. (distmat(i,j)>TLen .and. distmat(i,k)>TLen) )cycle
       d2 = geom(:,k)-geom(:,j)
       d2 = d2/dnrm2(3,d2,1)
       do l=i+1,natoms
         if(l==j .or. l==k  .or. (hasOOP(l).and.hasOOP(i)) .or. &
             (distmat(j,l)>TLen .and. distmat(k,l)>TLen)  )cycle
         d3 = geom(:,l)-geom(:,j)
         d3 = d3/dnrm2(3,d3,1)
         cpd(1) = d1(3)*(d2(2)-d3(2))+d2(3)*d3(2)-d2(2)*d3(3)+d1(2)*(d3(3)-d2(3))
         cpd(2) =-d2(3)*d3(1)+d1(3)*(d3(1)-d2(1))+d1(1)*(d2(3)-d3(3)) +d2(1)*d3(3)
         cpd(3) = d1(2)*(d2(1)-d3(1))+d2(2)*d3(1)-d2(1)*d3(2)+d1(1)*(d3(2)-d2(2))
         print "(2x,4(I5,3x),F12.4)",I,J,K,L, 90/Acos(0d0)* &
           asin((-d1(3)*d2(2)*d3(1)+d1(2)*d2(3)*d3(1)+d1(3)*d2(1)*d3(2)       &
                -d1(1)*d2(3)*d3(2)-d1(2)*d2(1)*d3(3)+d1(1)*d2(2)*d3(3))/      &
               dnrm2(3,cpd,1))
         if(ShortList)then
            hasOOP(l)=.true.
            hasOOP(i)=.true.
            if(distmat(i,j)<TLen.and.distmat(l,j)<TLen)hasOOP(k)=.true.
            if(distmat(i,k)<TLen.and.distmat(l,k)<TLen)hasOOP(j)=.true.
         end if
       end do! l
     end do!k
    end do !j
  end do   !i    
end subroutine analysegeom

