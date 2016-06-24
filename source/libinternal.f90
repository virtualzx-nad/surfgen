!---------------------------------------------------
! Generate the generalized inverse of B matrix:
! ginv(B).B**T=I
! if B=U.D.V**T, where U,V are orthogonal and D is diagonal, then
! ginv(B)=U.D**-1.V**T
SUBROUTINE ginv(m,n,bmat,ldb,invB,tol)
  implicit none

  INTEGER,INTENT(IN)                            :: m,n,ldb
  DOUBLE PRECISION,dimension(ldb,n),INTENT(IN)  :: bmat
  DOUBLE PRECISION,dimension(m,n),INTENT(OUT)   :: invB
  DOUBLE PRECISION,INTENT(IN)                   :: tol
  DOUBLE PRECISION,dimension(m,m)               :: U
  DOUBLE PRECISION,dimension(n,n)               :: VT
  double precision,dimension(min(m,n))          :: w
  double precision,dimension(:),allocatable     :: work
  INTEGER                                       :: INFO,i,lwork,rk
  
  invB=bmat(1:m,:)
  ! B=U.D.V**T
  CALL DGESVD('S','S',m,n,invB,ldb,w,U,m,VT,n,W,-1,INFO)
  if(INFO/=0)then
    print *,"INFO=",info
    stop "ginv: workspace query failed."
  end if
  LWORK=int(W(1))
  allocate(WORK(LWORK))
  CALL DGESVD('S','S',m,n,invB,ldb,w,U,m,VT,n,WORK,LWORK,INFO)
  deallocate(WORK)
  if(INFO/=0)then
    print *,"INFO=",info
    stop "ginv: DGESVD failed."
  end if
  rk=min(m,n)
  do i=1,min(m,n)
    if(abs(w(i))<tol)then
      rk=i-1
      exit
    else
      CALL DSCAL(m,1/w(i),U(:,i),int(1))
    end if
  end do
  ! generating binv using the factorization
  CALL DGEMM('n','n',m,n,rk,dble(1),U,m,VT,n,dble(0),invB,m)
END SUBROUTINE ginv

!***********************************************************************
! swaps the index of two atoms and change the parity
SUBROUTINE swapInd(i1,i2,parity)
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: i1,i2,parity
  integer tmp
  tmp = i1
  i1  = i2
  i2  = tmp
  parity  = - parity
END SUBROUTINE swapInd

!***********************************************************************
!reorder tetrahedron oop atoms to canonical order and returns the parity of perm
!this subroutine is used for tetrahedron coordinates where any permutations of 
!the 4 atoms yield the same coordinate, up to a sign.
SUBROUTINE reorderOOP(oldOOP,newOOP,parity)
  implicit none
  INTEGER,DIMENSION(4),INTENT(IN)   ::  oldOOP
  INTEGER,DIMENSION(4),INTENT(OUT)  ::  newOOP
  INTEGER,INTENT(OUT)               ::  parity
  integer          ::  loc, i
  parity=1
  newOOP=oldOOP
  do i=4,2,-1
    loc   = maxloc(newOOP(1:i),1)
    if(loc/=i) call swapInd(newOOP(loc),newOOP(i),parity)
  end do
END SUBROUTINE reorderOOP

!***********************************************************************
!reorder umbrella oop atoms to canonical order and returns the parity of perm
!this subroutine is used for umbrella coordinate where any permutation of 3
!in the 4 atoms yield the same coordinate, up to a sign.
SUBROUTINE reorderOOP2(oldOOP,newOOP,parity)
  implicit none
  INTEGER,DIMENSION(4),INTENT(IN)   ::  oldOOP
  INTEGER,DIMENSION(4),INTENT(OUT)  ::  newOOP
  INTEGER,INTENT(OUT)               ::  parity
  integer          ::  loc, i
  parity=1
  newOOP=oldOOP
  do i=4,3,-1
    loc   = maxloc(newOOP(2:i),1)+1
    if(loc/=i) call swapInd(newOOP(loc),newOOP(i),parity)
  end do
END SUBROUTINE reorderOOP2

!***********************************************************************
!reorder atom reference of torsion angle to canonical order and returns
! the parity of permutation
SUBROUTINE reorderTORS(oldTOR,newTOR,parity)
  implicit none
  INTEGER,DIMENSION(4),INTENT(IN)   ::  oldTOR
  INTEGER,DIMENSION(4),INTENT(OUT)  ::  newTOR
  INTEGER,INTENT(OUT)               ::  parity
  if(oldTOR(2)>oldTor(3))then
    newTOR(1)=oldTOR(4)
    newTOR(2)=oldTOR(3)
    newTOR(3)=oldTOR(2)
    newTOR(4)=oldTOR(1)
    parity = -1
  else
    parity = 1
    newTOR = oldTOR
  end if
END SUBROUTINE reorderTORS

!***********************************************************************
! reorder atom reference in 4center dot product to dictionary order
! and returns the parity of the permutations
! Dictionary order would be a1<a2, a3<a4, and a1<a3
SUBROUTINE reorderDot4(oldDot4,newDot4,parity)
  IMPLICIT NONE
  INTEGER,DIMENSION(4),INTENT(IN)  :: oldDot4
  INTEGER,DIMENSION(4),INTENT(OUT) :: newDot4
  INTEGER,INTENT(OUT)              :: parity

  parity = 1
  newDot4=oldDot4
  if(newDot4(1)>newDot4(2))call swapInd(newDot4(1),newDot4(2),parity)
  if(newDot4(3)>newDot4(4))call swapInd(newDot4(3),newDot4(4),parity)
  if(newDot4(1)>newDot4(3))then
        call swapInd(newDot4(1),newDot4(3),parity)
        call swapInd(newDot4(2),newDot4(4),parity)
  end if
END SUBROUTINE reorderDot4

!***********************************************************************
! reorder atoms for anti-symmetric bends
! only atoms 3 and 4 needs to be ordered
SUBROUTINE reorderABend(oldAtms,newAtms,parity)
  IMPLICIT NONE
  INTEGER,DIMENSION(4),INTENT(IN)  :: oldAtms
  INTEGER,DIMENSION(4),INTENT(OUT) :: newAtms
  INTEGER,INTENT(OUT)              :: parity

  parity = 1
  newAtms=oldAtms
  if(newAtms(3)>newAtms(4))call swapInd(newAtms(3),newAtms(4),parity)
END SUBROUTINE

!
! Builds a b-matrix given a set of internal coordinates definitions
!
!
SUBROUTINE buildWBmat(cgeom,igeom,bmat)
  use hddata, only: ncoord
  use progdata, only: natoms,coordmap,CoordSet
  IMPLICIT NONE
  DOUBLE PRECISION,dimension(3*natoms),INTENT(IN)           :: cgeom
  DOUBLE PRECISION,dimension(ncoord),INTENT(OUT)            :: igeom
  DOUBLE PRECISION,dimension(ncoord,3*natoms),INTENT(OUT)   :: bmat
  INTEGER                                :: i,j,offs,m,n,atms(4),cmap(4)

! bval     Derivative of a unscaled coordinate with respect to cartesian 
!            coordinates of the four reference atoms.
  DOUBLE PRECISION     ::  bval(12), w(3), dwdR(3,12), dw(12),ss,fs,s(3),igm

! intialize bmat
  bmat = 0d0
  do i=1,ncoord
    m=coordmap(i,1) !index of set
    n=coordmap(i,2) !index in set
    select case(CoordSet(m)%Type)

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     Plain or scaled Rij coordinates
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      case(0)   

        call calcwij(CoordSet(m)%Scaling,CoordSet(m)%coord(1,n),CoordSet(m)%coord(2,n),&
                      CoordSet(m)%coef,cgeom,igeom(i),bval)
        do j=1,2
            offs = 3*(CoordSet(m)%coord(j,n)-1)
            bmat(i,offs+1:offs+3) =  bval(3*j-2:3*j)
        end do!j=1,2

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     OOP ANGLE
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      case(-1)  !oop

        select case(CoordSet(m)%Scaling)

          case(0)  ! Reciprocal scaling   Det/(Product[rij])^a
            call  oop(natoms,CoordSet(m)%coord(1,n),cgeom,igeom(i),bval,coordset(m)%coef(1))
            igeom(i) = igeom(i)*coordset(m)%coef(2)
            do j=1,4
              offs = 3*(CoordSet(m)%coord(j,n)-1)
              bmat(i,offs+1:offs+3) =  bval(3*j-2:3*j)*coordset(m)%coef(2)
            end do!j=1,4
            
          case default  ! Exponential scaling using one of the linear scaling types.
            call oop2(natoms,CoordSet(m)%coord(1,n),cgeom,igeom(i),bval,coordset(m)%coef(1),&
                                                    CoordSet(m)%Scaling,coordset(m)%coef(2))
            do j=1,4
              offs = 3*(CoordSet(m)%coord(j,n)-1)
              bmat(i,offs+1:offs+3) =  bval(3*j-2:3*j)
            end do!j=1,4

        end select !case(CoordSet(m)%Scaling)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     Special OOP angles
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      case(-2)  !special oop angles
        if(CoordSet(m)%Scaling>0)then
            call  oop(natoms,CoordSet(m)%coord(1,n),cgeom,igeom(i),bval,5d-1)
            ! scaling with exponentials
            dwdR = 0d0
            do j=2,4
                call calcwij(coordset(m)%scaling,CoordSet(m)%coord(1,n), &
                     CoordSet(m)%coord(j,n),coordset(m)%coef,cgeom,w(j-1),dw)
                if(abs(w(j-1))<1D-30) w(j-1)=sign(1D-30,w(j-1))
                dwdR(j-1,1:3) = dw(1:3)
                dwdR(j-1,j*3-2:j*3) = dw(4:6)
            end do
            s(1) = w(2)*w(3)
            s(2) = w(3)*w(1)
            s(3) = w(1)*w(2)
            ss = sum(s)
            fs = w(1)*s(1)/ss
            bval=bval*fs+igeom(i)*(s(1)**2*dwdR(1,:)+s(2)**2*dwdR(2,:)+s(3)**2*dwdR(3,:))/ss**2
            igeom(i)=igeom(i)*fs
        else
            call oop3(natoms,CoordSet(m)%coord(1,n),cgeom,igeom(i),bval,&
                coordset(m)%coef(1),CoordSet(m)%Scaling,coordset(m)%coef(2))
        end if
        do j=1,4
          offs = 3*(CoordSet(m)%coord(j,n)-1)
          bmat(i,offs+1:offs+3) =  bval(3*j-2:3*j)
        end do!j=1,4

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     Four-center dot-product coordinate
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     case (-3) ! (a1-a2).(a3-a4)/norm(a1-a2)/norm(a3-a4) *scaling
        call dot4(natoms,CoordSet(m)%coord(1,n),cgeom,igeom(i),bval,&
                Coordset(m)%coef,CoordSet(m)%Scaling)
        do j=1,4
            offs = 3*(CoordSet(m)%coord(j,n)-1)
            bmat(i,offs+1:offs+3) =  bval(3*j-2:3*j)
        end do!j=1,4

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     Scaled bond angle
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      case (1) ! angle
        call bend(natoms,CoordSet(m)%coord(1,n),cgeom,igeom(i),bval,&
                coordset(m)%coef(1),CoordSet(m)%Scaling,coordset(m)%coef(2))
        do j=1,3
            offs = 3*(CoordSet(m)%coord(j,n)-1)
            bmat(i,offs+1:offs+3) =  bval(3*j-2:3*j)
        end do!j=1,3

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     Antisymmetric bendings of the following motion
!                1 ->
!                |
!                2
!              /   \
!             3     4
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      case(2)
        cmap=[1,3,2,4]
        atms=CoordSet(m)%coord(cmap,n) 
        call bend(natoms,atms,cgeom,igeom(i),bval,&
                coordset(m)%coef(1),CoordSet(m)%Scaling,coordset(m)%coef(2))
        do j=1,3
            offs = 3*(CoordSet(m)%coord(cmap(j),n)-1)
            bmat(i,offs+1:offs+3) =  bval(3*j-2:3*j)
        end do!j=1,3
        cmap=[1,4,2,3]
        atms=CoordSet(m)%coord(cmap,n) 
        call bend(natoms,atms,cgeom,igm,bval,&
                coordset(m)%coef(1),CoordSet(m)%Scaling,coordset(m)%coef(2))
        igeom(i)=igeom(i)-igm
        do j=1,3
            offs = 3*(CoordSet(m)%coord(cmap(j),n)-1)
            bmat(i,offs+1:offs+3) = bmat(i,offs+1:offs+3)-bval(3*j-2:3*j)
        end do!j=1,3

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     other types.  shouldn't get here
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      case default
        print *,"CoordSet(",m,")%Type    = ",CoordSet(m)%Type
        stop "Unsupported coordinate type"

    end select!case(CoordSet(m)%Type)

  end do!i=1,ncoord
end SUBROUTINE buildWBmat


!********************************************************************************************
!This subroutine calculates the bond length between two reference atoms and the derivative of
! this length with respect to the cartesian coordinates of the two atoms.
SUBROUTINE stre(na,a1,a2,cgeom,qval,bval)
  IMPLICIT NONE
  INTEGER,INTENT(IN)                           :: na,a1,a2
  DOUBLE PRECISION,INTENT(IN)                  :: cgeom(3*na)
  DOUBLE PRECISION,INTENT(INOUT)               :: qval
  DOUBLE PRECISION,dimension(12),INTENT(INOUT) :: bval
  INTEGER                                      :: i
  DOUBLE PRECISION,dimension(3)                :: vec1

  call dvec(na,a1,a2,cgeom,vec1,qval)
  do i = 1,3
   bval(i)   = vec1(i)
   bval(3+i) = -vec1(i)
   bval(6+i) = 0D0
   bval(9+i) = 0D0
  enddo

  return
end SUBROUTINE stre

!********************************************************************************************
! OOP angle values and derivatives
SUBROUTINE oop(natoms,atms,cgeom,qval,bval,scale)
  IMPLICIT NONE
  INTEGER,INTENT(IN)                              :: natoms
  INTEGER,DIMENSION(4),INTENT(IN)                 :: atms
  DOUBLE PRECISION,DIMENSION(natoms*3),INTENT(IN) :: cgeom
  DOUBLE PRECISION,INTENT(IN)                     :: scale
  DOUBLE PRECISION,INTENT(OUT)                    :: qval
  DOUBLE PRECISION,DIMENSION(12),INTENT(OUT)      :: bval

  INTEGER           :: i, j, sgns(6,4)
  DOUBLE PRECISION  :: denom, dnorm2(6), dvecs(3,6),&
                       geom(4,3), ddvec(3,3)
  DOUBLE PRECISION  :: factor  ! normalization

! The oriented direction vectors between all 6 atom pairs
  sgns(1,:)=(/ 1, 0, 0,-1/)
  sgns(2,:)=(/ 0, 1, 0,-1/)
  sgns(3,:)=(/ 0, 0, 1,-1/)
  sgns(4,:)=(/ 1, 0,-1, 0/)
  sgns(5,:)=(/ 0, 1,-1, 0/)
  sgns(6,:)=(/ 1,-1, 0, 0/)

! extract the geometry of the 4 involved atoms
  do i=1,4
    geom(i,:)=cgeom(atms(i)*3-2:atms(i)*3)
  end do

! calculate displacement vectors between these atoms
  dvecs=transpose(matmul(sgns,geom))
  do i=1,6
    dnorm2(i)=dot_product(dvecs(:,i),dvecs(:,i))
  end do

! calculate value of scaled OOP angle and its derivatives
  denom=exp(log(product(dnorm2))/2*scale)
  if(abs(denom)<1d-30)stop "OOP angle cannot be defined when atoms coincide."
  qval=det3(dvecs)/denom
  do i=1,3
    do j=1,3
      ddvec(i,j)=sum(sgns(:,i)*dvecs(j,:)/dnorm2(:))
    end do
  end do
  bval(1:3)=cross(dvecs(:,2),dvecs(:,3))/denom-qval*scale*ddvec(1,:)
  bval(4:6)=cross(dvecs(:,3),dvecs(:,1))/denom-qval*scale*ddvec(2,:)
  bval(7:9)=cross(dvecs(:,1),dvecs(:,2))/denom-qval*scale*ddvec(3,:)
  bval(10:12)=-bval(1:3)-bval(4:6)-bval(7:9)

! calculate the normalization factor
  factor = 2D0 ** (scale*3-3)*3D0**(scale*3)
  qval = qval * factor
  bval = bval * factor

CONTAINS

  ! this function calculates the value of a 3x3 determinant
  FUNCTION det3(m)
    IMPLICIT NONE
    double precision ::  det3, m(:,:)
    det3=(m(1,2)*m(2,3)-m(1,3)*m(2,2))*m(3,1)+&
         (m(1,3)*m(2,1)-m(1,1)*m(2,3))*m(3,2)+&
         (m(1,1)*m(2,2)-m(1,2)*m(2,1))*m(3,3)
  END FUNCTION det3

  !this function calculate the cross product of two 3D vectors
  FUNCTION cross(v1,v2)RESULT(v3)
    implicit none
    double precision,dimension(3) :: v1,v2,v3
    v3(1)=-v1(3)*v2(2)+v1(2)*v2(3)
    v3(2)=v1(3)*v2(1)-v1(1)*v2(3)
    v3(3)=-v1(2)*v2(1)+v1(1)*v2(2)
  END FUNCTION cross

END SUBROUTINE oop

!
!
!
SUBROUTINE dvec(na,a1,a2,cgeom,vec,vnorm)
  IMPLICIT NONE
  INTEGER, INTENT(IN)                         :: na,a1,a2
  DOUBLE PRECISION,INTENT(IN)                 :: cgeom(3*na)
  DOUBLE PRECISION,INTENT(INOUT)              :: vnorm
  DOUBLE PRECISION,dimension(3),INTENT(INOUT) :: vec

  INTEGER                                     :: i
  DOUBLE PRECISION                            :: denom

  vnorm = dble(0)
  do i = 1,3
   vec(i) = cgeom(3*(a1-1)+i)-cgeom(3*(a2-1)+i)
   vnorm = vnorm + vec(i)**2
  enddo

  vnorm = Sqrt(vnorm)
  if(vnorm<1d-30)then
   denom = dble(1)
  else
   denom = vnorm
  endif

  vec=vec/denom

  return
end SUBROUTINE dvec

!********************************************************************************
! CalcWij  
! Purpose
!         Calculate the value and gradient of one scaled Rij coordinate
! Arguments
! scaling  [input] INTEGER
!          Type of coordinate scaling method
! a1,a2    [input] INTEGER
!          Indices of end point atoms that defines the distance
! coef     [input] DOUBLE PRECISION,dimension(2)
!          Scaling coefficients
! cgeom    [input] 
SUBROUTINE calcwij(scaling,a1,a2,coef,cgeom,w,dwdR)
  use progdata, only: natoms
  IMPLICIT NONE
  INTEGER,INTENT(IN)                                :: scaling,a1,a2
  DOUBLE PRECISION,dimension(2),INTENT(IN)          :: coef
  DOUBLE PRECISION,dimension(3*natoms),INTENT(IN)   :: cgeom
  DOUBLE PRECISION,INTENT(OUT)                      :: w
  DOUBLE PRECISION,dimension(12),INTENT(OUT)        :: dwdR

! fval     Value of current coordinate in its unscaled form
! fbmat    Cartesian gradients of current coordinate in its unscaled form
! bval     Derivative of a unscaled coordinate with respect to cartesian 
!            coordinates of the four reference atoms.

  DOUBLE PRECISION     ::  fval, bval(12), g

!---Plain or scaled Rij coordinates--->>>

  call stre(natoms,a1,a2,cgeom,fval,bval)
  if(scaling<0) stop  "Invalid scaling type in CalcWij"
  select case(scaling)
    !  no scalings
    case(0)
      w    = (fval-coef(2))*coef(1)
      dwdR = bval*coef(1)

    !  Morse functions   Exp(-c1*(r-c2))
    case(1)
      w    =  exp(-coef(1)*(fval-coef(2)))
      dwdR = -bval*w*coef(1)

    !  Gaussian functions Exp(-c1*(r-c2)**2)
    case(2)
      w    =  exp(-coef(1)*(fval-coef(2))**2)
      dwdR = -bval*w*coef(1)*2*(fval-coef(2))

    !  Screened Columb potential Exp(-c1*(r-c2))/r     (Yukawa, leading term)
    case(3)
      g    =  exp(-coef(1)*(fval-coef(2)))
      w    =  g/(fval+1D-40)
      dwdR = -bval*w*(coef(1)+1/(fval+1D-40))

    !  Long range term of screened Columb w=Exp[c1*(c2-x)]*(x/c2)^(c1*c2)
    !  its derivative: w*c1*(c2-x)/x*x'
    case(4)
      w    = exp(coef(1)*(coef(2)-fval))*(fval/coef(2))**(coef(1)*coef(2))
      dwdR = bval*w*coef(1)*(coef(2)-fval)/fval

    !  Lennard Jones functions (c2/r)**c1
    case(5)
      w    =  (coef(2)/fval)**coef(1)
      dwdR = -bval*w*coef(1)/fval

    !  Shifted(chasmless) Yukawa exp(-c1*(r-c2))/(r+c2)
    case(6)
      g    =  exp(-coef(1)*(fval-coef(2)))
      w    =  g/(fval+coef(2))
      dwdR = -bval*w*(coef(1)+1/(fval+coef(2)))

    ! segmentation function tanh
    case(7)
      w    = tanh((fval-coef(2))/coef(1))/2
      dwdR = bval/(cosh((fval-coef(2))/coef(1))**2*coef(1))/2

    case default
      print *,"scaling = ",scaling
      stop "Unsupported bond distance scaling method."

  end select !case(scaling)
  if(abs(w)<1d-30)w=1d-30
END SUBROUTINE calcwij

!************************************************************************************
! SUBROUTINE oop2 
! OOP angle values and derivatives, with the scaling method of distance coordinates
! f= g*Product[ Scale[rij] ]
! where g is the triple product of displacement vectors
!************************************************************************************
! natom s [input] INTEGER
!         Total number of atoms
! atms    [input] INTEGER, dimension(4) 
!         The indices of the four reference atoms
! cgeom   [input] DOUBLE PRECISION, dimension(3*natoms)
!         Cartesian geometry of all the atoms
! qval    [output] DOUBLE PRECISION
!         Value of the scaled out-of-plane coordinate
! bval    [output] DOUBLE PRECISION, dimension(12)
!         Derivatives of the OOP coordinate with respect to the cartesian coordinate
! scale   [input] DOUBLE PRECISION
!         Scaling factor used in the exponetials 
! sctype  [input] INTEGER
!         Type of the scaling function used.
SUBROUTINE oop2(natoms,atms,cgeom,qval,bval,scale,sctype,factor)
  IMPLICIT NONE
  INTEGER,INTENT(IN)                              :: natoms,sctype
  INTEGER,DIMENSION(4),INTENT(IN)                 :: atms
  DOUBLE PRECISION,DIMENSION(natoms*3),INTENT(IN) :: cgeom
  DOUBLE PRECISION,INTENT(OUT)                    :: qval
  DOUBLE PRECISION,INTENT(IN)                     :: scale
  DOUBLE PRECISION,DIMENSION(12),INTENT(OUT)      :: bval
  double precision,intent(IN)                     :: factor 

  INTEGER           :: i, j, k, sgns(6,4), tp
  DOUBLE PRECISION  :: dnorm2(6), dvecs(6,3),   &
                       geom(4,3), ddvec(3,3),   dw(12),       &
                       coef(2), dwdR(6,12), w(6),  denom, ddenrm
! eprod = product of the scaling exponentials
! fval  = value of unscaled oop coordinate (the triple product)
! fbmat = derivatives of unscaled OOP coordinate
  DOUBLE PRECISION  :: eprod, fval, fbmat(12),efactor

! The oriented direction vectors between all 6 atom pairs
  sgns(1,:)=(/ 1, 0, 0,-1/)
  sgns(2,:)=(/ 0, 1, 0,-1/)
  sgns(3,:)=(/ 0, 0, 1,-1/)
  sgns(4,:)=(/ 1, 0,-1, 0/)
  sgns(5,:)=(/ 0, 1,-1, 0/)
  sgns(6,:)=(/ 1,-1, 0, 0/)

! calculate the scaling factors and their derivatives
  coef(1) = scale
  coef(2) = factor
  tp  = sctype
  if(sctype<0) tp = 0
  k = 1
  dwdR=0d0
  do i=1,3
   do j=i+1,4
     call calcwij(tp,atms(i),atms(j),coef,cgeom,w(k),dw)
     dwdR(k,i*3-2:i*3) = dw(1:3)
     dwdR(k,j*3-2:j*3) = dw(4:6)
     k = k+1
   end do !j
  end do


! extract the geometry of the 4 involved atoms
  do i=1,4
    geom(i,:)=cgeom(atms(i)*3-2:atms(i)*3)
  end do

! calculate displacement vectors between these atoms
  dvecs=matmul(sgns,geom)
  do i=1,6
    dnorm2(i)=dot_product(dvecs(i,:),dvecs(i,:))
  end do

! calculate value of unscaled OOP angle and its derivatives
  fval=det3(dvecs(1,1),6)
  do i=1,3
    do j=1,3
      ddvec(i,j)=sum(sgns(:,i)*dvecs(:,j)/dnorm2(:))
    end do
  end do
  CALL cross(dvecs(2,1),6,dvecs(3,1),6,fbmat(1),1)
  CALL cross(dvecs(3,1),6,dvecs(1,1),6,fbmat(4),1)
  CALL cross(dvecs(1,1),6,dvecs(2,1),6,fbmat(7),1)
  fbmat(10:12)=-fbmat(1:3)-fbmat(4:6)-fbmat(7:9)

! calculate the scaled oop
  if(sctype<0)then
    denom = (factor/sum(w))**(scale+3)
    qval = fval * denom
    bval = fbmat* denom - qval*(scale+3)/sum(w)*&
              (dwdR(1,:)+dwdR(2,:)+dwdR(3,:)+dwdR(4,:)+dwdR(5,:)+dwdR(6,:))
  else
    bval = 0d0
    eprod = product(w)
    qval = fval * eprod
    do i=1,6
      efactor = 1D0
      do j=1,6
        if(i.ne.j) then
          efactor=efactor*w(j)
        end if
      end do
      bval = bval + dwdR(i,:)*efactor
    end do
    bval = fbmat*eprod + fval*bval
  end if
CONTAINS

  ! this function calculates the value of a 3x3 determinant
  FUNCTION det3(m,ldm)
    IMPLICIT NONE
    integer :: ldm
    double precision ::  det3, m(ldm,3)
    det3=(m(1,2)*m(2,3)-m(1,3)*m(2,2))*m(3,1)+&
         (m(1,3)*m(2,1)-m(1,1)*m(2,3))*m(3,2)+&
         (m(1,1)*m(2,2)-m(1,2)*m(2,1))*m(3,3)
  END FUNCTION det3

  !this function calculate the cross product of two 3D vectors
  SUBROUTINE cross(v1,iv1,v2,iv2,v3,iv3)
    implicit none
    integer,intent(in)  ::  iv1,iv2,iv3  !increment of v1,v2 and v3
    double precision,dimension(*) :: v1,v2,v3
    integer             ::  i12,i13,i22,i23,i32,i33
    i12=1+iv1
    i13=1+2*iv1
    i22=1+iv2
    i23=1+2*iv2
    i32=1+iv3
    i33=1+2*iv3

    v3(  1)=-v1(i13)*v2(i22)+v1(i12)*v2(i23)
    v3(i32)= v1(i13)*v2(  1)-v1(  1)*v2(i23)
    v3(i33)=-v1(i12)*v2(  1)+v1(  1)*v2(i22)
  END SUBROUTINE cross

END SUBROUTINE oop2

!*******************************************************************************
! SUBROUTINE oop3 
! OOP angle values and derivatives, with the expoential scaling method
! f= C1*g*Product[Scale[rij]]
! where a is the scaling factor, g is the triple product of unit displacement vectors
! here, one of the four atoms is made a special vertex and only the three distances
! launching from this atom is used for scaling.  This the three atom atoms are
! permutationally equivalent while this one is special.
! ******************************************************************************
! natom s [input] INTEGER
!         Total number of atoms
! atms    [input] INTEGER, dimension(4) 
!         The indices of the four reference atoms
! cgeom   [input] DOUBLE PRECISION, dimension(3*natoms)
!         Cartesian geometry of all the atoms
! qval    [output] DOUBLE PRECISION
!         Value of the scaled out-of-plane coordinate
! bval    [output] DOUBLE PRECISION, dimension(12)
!         Derivatives of the OOP coordinate with respect to the cartesian coordinate
! scale   [input] DOUBLE PRECISION
!         Scaling factor used in the exponetials 
! sctype  [input] INTEGER
!         Type of the scaling function used.
SUBROUTINE oop3(natoms,atms,cgeom,qval,bval,scale,sctype,factor)
  IMPLICIT NONE
  INTEGER,INTENT(IN)                              :: natoms,sctype
  INTEGER,DIMENSION(4),INTENT(IN)                 :: atms
  DOUBLE PRECISION,DIMENSION(natoms*3),INTENT(IN) :: cgeom
  DOUBLE PRECISION,INTENT(OUT)                    :: qval
  DOUBLE PRECISION,INTENT(IN)                     :: scale
  DOUBLE PRECISION,DIMENSION(12),INTENT(OUT)      :: bval
  double precision,intent(IN)                     :: factor 

  INTEGER           :: i, j, sgns(6,4)
  DOUBLE PRECISION  :: dnorm2(6), dvecs(6,3), s(3),  pw, tdv(3,6),  &
                       geom(4,3), ddvec(3,3), ss,    dw(12),    &
                       coef(2), dwdR(3,12), w(3), dwdR2(3,12), w2(3)
! eprod = product of the scaling exponentials
! fval  = value of unscaled oop coordinate (the triple product)
! fbmat = derivatives of unscaled OOP coordinate
  DOUBLE PRECISION  :: eprod, fval, fbmat(12),a,p

! The oriented direction vectors between all 6 atom pairs
  sgns(1,:)=(/ 1, 0, 0,-1/)
  sgns(2,:)=(/ 0, 1, 0,-1/)
  sgns(3,:)=(/ 0, 0, 1,-1/)
  sgns(4,:)=(/ 1, 0,-1, 0/)
  sgns(5,:)=(/ 0, 1,-1, 0/)
  sgns(6,:)=(/ 1,-1, 0, 0/)

! calculate the three scaling factors and their derivatives
! scaled distances for the three bonds A1-A2, A1-A3 and A1-A4 are used
  dwdR = 0d0
  coef(1)   = 1d0
  coef(2)   = 0d0
  do i=2,4
! scaling mode=0   (unscaled).  unscaled distances are produced here
     call calcwij(int(0),atms(1),atms(i),coef,cgeom,w(i-1),dw)
     if(abs(w(i-1))<1D-30) w(i-1)=sign(1D-30,w(i-1))
     dwdR(i-1,1:3) = dw(1:3)
     dwdR(i-1,i*3-2:i*3) = dw(4:6)
  end do

! here are the scaled distances
  coef(1) =scale
  coef(2) =factor
  dwdR2 = 0d0
  if(sctype.ne.0)then
    do i=2,4
     call calcwij(abs(sctype),atms(1),atms(i),coef,cgeom,w2(i-1),dw)
     if(abs(w2(i-1))<1D-50) w2(i-1)=sign(1D-50,w2(i-1))
     dwdR2(i-1,1:3)= dw(1:3)
     dwdR2(i-1,i*3-2:i*3)= dw(4:6)
    end do
  end if

! extract the geometry of the 4 involved atoms
  do i=1,4
    geom(i,:)=cgeom(atms(i)*3-2:atms(i)*3)
  end do

! calculate displacement vectors between these atoms
  dvecs=matmul(sgns,geom)
  do i=1,6
    dnorm2(i)=dot_product(dvecs(i,:),dvecs(i,:))
  end do

! calculate value of unscaled OOP coordinate and its derivatives
  fval=det3(dvecs)
  do i=1,3
    do j=1,3
      ddvec(i,j)=sum(sgns(:,i)*dvecs(:,j)/dnorm2(:))
    end do
  end do
  tdv=transpose(dvecs)
  fbmat(1:3)=cross(tdv(:,2),tdv(:,3))
  fbmat(4:6)=cross(tdv(:,3),tdv(:,1))
  fbmat(7:9)=cross(tdv(:,1),tdv(:,2))
  fbmat(10:12)=-fbmat(1:3)-fbmat(4:6)-fbmat(7:9)

  pw = product(w)
  fval=fval/pw
  fbmat=fbmat/pw-fval*(dwdR(1,:)/w(1)+dwdR(2,:)/w(2)+dwdr(3,:)/w(3))

! calculate the scaled oop
  if(sctype.lt.0)then
! use the harmonic mean of scaled distances to scale OOP
    s(1) = w2(2)*w2(3)
    s(2) = w2(3)*w2(1)
    s(3) = w2(1)*w2(2)
    ss = sum(s)
    a = w2(1)*s(1)/ss
    qval= fval*a
    bval=fbmat*a+fval*(s(1)**2*dwdR2(1,:)+s(2)**2*dwdR2(2,:)+s(3)**2*dwdR2(3,:))/ss**2
  else if(sctype.eq.0)then
! reciprocal or order scaling
    do i=1,3
        if(abs(w(i))<1D-20)w(i)=sign(1D-20,w(i))
    end do
    p=coef(2)**3/product(w)
    a=p**coef(1)
    qval = a*fval
    bval = a*fbmat-coef(1)*qval*(dwdR(1,:)/w(1)+dwdR(2,:)/w(2)+dwdR(3,:)/w(3))
  else
! multiply all 3 scaled length scaling functions
    eprod = product(w2)
    qval = fval*eprod
    bval = fbmat*eprod+fval*(w2(2)*w2(3)*dwdR2(1,:)+w2(1)*w2(3)*dwdR2(2,:)+w2(1)*w2(2)*dwdR2(3,:))
  end if

CONTAINS

  ! this function calculates the value of a 3x3 determinant
  FUNCTION det3(m)
    IMPLICIT NONE
    double precision ::  det3, m(:,:)
    det3=(m(1,2)*m(2,3)-m(1,3)*m(2,2))*m(3,1)+&
         (m(1,3)*m(2,1)-m(1,1)*m(2,3))*m(3,2)+&
         (m(1,1)*m(2,2)-m(1,2)*m(2,1))*m(3,3)
  END FUNCTION det3

  !this function calculate the cross product of two 3D vectors
  FUNCTION cross(v1,v2)RESULT(v3)
    implicit none
    double precision,dimension(3) :: v1,v2,v3
    v3(1)=-v1(3)*v2(2)+v1(2)*v2(3)
    v3(2)=v1(3)*v2(1)-v1(1)*v2(3)
    v3(3)=-v1(2)*v2(1)+v1(1)*v2(2)
  END FUNCTION cross

END SUBROUTINE oop3

!*******************************************************************************
! bending angle with specific scaling
!*******************************************************************************
SUBROUTINE BEND(natoms,atms,cgeom,qval,bval,scale,sctype,factor)
  IMPLICIT NONE
  INTEGER,INTENT(IN)                              :: natoms,sctype
  INTEGER,DIMENSION(3),INTENT(IN)                 :: atms
  DOUBLE PRECISION,DIMENSION(3,natoms),INTENT(IN) :: cgeom
  DOUBLE PRECISION,INTENT(OUT)                    :: qval
  DOUBLE PRECISION,INTENT(IN)                     :: scale
  DOUBLE PRECISION,DIMENSION(9),INTENT(OUT)       :: bval
  double precision,intent(IN)                     :: factor

  double precision ::  v, dv(9), d1(3),d2(3), d12, d22, dd12(9), dd22(9)
  double precision ::  p, dp(9), tmp, pw
  integer,dimension(3,9) ::   dd1,dd2    ! derivatives of d1 and d2
  integer  ::  i

  dd1(1,:) =[1,0,0,0,0,0,-1,0,0]
  dd1(2,:) =[0,1,0,0,0,0,0,-1,0]
  dd1(3,:) =[0,0,1,0,0,0,0,0,-1]
  dd2(1,:) =[0,0,0,1,0,0,-1,0,0]
  dd2(2,:) =[0,0,0,0,1,0,0,-1,0]
  dd2(3,:) =[0,0,0,0,0,1,0,0,-1]

  ! calculate the dot product of the two borders of the angle and its derivative
  ! v is the dot product, and dv is the derivative

  d1 = cgeom(:,atms(1))-cgeom(:,atms(3))
  d2 = cgeom(:,atms(2))-cgeom(:,atms(3))
  v = dot_product(d1,d2)
  do i=1,9
    dv(i) = dot_product(d1,dd2(:,i))+dot_product(dd1(:,i),d2)
  end do!i=1,9

  ! d12 and d22 are the square norm of the two borders, and dd12 and dd22 are the
  ! derivatives of these two quantities.
  d12 = dot_product(d1,d1)
  d22 = dot_product(d2,d2)
  do i=1,9
    dd12(i) = 2*dot_product(d1,dd1(:,i))
    dd22(i) = 2*dot_product(d2,dd2(:,i))
  end do!i=1,9

  ! calculate the cosine of angle p = v/sqrt(d12*d22)
  p = v/sqrt(d12*d22)
  dp =  ( dv - v/2*(dd12/d12+dd22/d22)  ) / sqrt(d12*d22)

  ! calculate final value and derivatives according to scaling method
  select case(sctype)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!       Cos theta, p
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    case ( 1 )
        qval = p
        bval = dp

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!       Raw angle. scaled by `scale` and shifted by `factor`
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  obtained as acos(p).    d acos(p)= - dp / Sqrt(1-p^2)
    case ( 0 )
        qval = (acos(p)-factor)*scale
        bval = -dp/sqrt(1-p**2)*scale

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!       Scaled cosine, p/( 1+exp[c1*(d1^2+d2^2-c2^2)] )
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!   q = p/( 1+exp[c1*(d1^2+d2^2-c2^2)] )
    case  (2)
        pw=scale*(d12+d22-factor**2)
        if(pw>70d0)then
          qval = 0d0
          bval = 0d0
        else
          tmp  = exp(pw)
          qval = p/(1+tmp)
          bval = (  dp  -  scale*tmp*p/(1+tmp)*(dd12+dd22)  ) /(1+tmp)
        end if
  end select!case(sctype)
END SUBROUTINE BEND

!*******************************************************************************
! 4 center dot product with scalings 
! f = (a1-a2).(a3-a4) / norm(a1-a2) /norm(a3-a4)
! scalings : =0 no scaling
!            >0 multiply the product of 4 scaled distances d13,d14,d23,d24
!            <0 multiply the harmonic mean of 4 scaled distances d13,d14,d23,d24
!*******************************************************************************
SUBROUTINE dot4(natoms,atms,cgeom,fval,bval,coef,scal)
  IMPLICIT NONE
  INTEGER,INTENT(IN)                              :: natoms,scal
  INTEGER,DIMENSION(4),INTENT(IN)                 :: atms
  DOUBLE PRECISION,DIMENSION(3,natoms),INTENT(IN) :: cgeom
  DOUBLE PRECISION,INTENT(OUT)                    :: fval
  DOUBLE PRECISION,DIMENSION(2),INTENT(IN)        :: coef
  DOUBLE PRECISION,DIMENSION(12),INTENT(OUT)       :: bval

  double precision  :: f, r12,r34  ! unscaled dot prod and the two distances 
  double precision,dimension(12) ::b, w12,w34 ! gradients of f, r12 and r34
  double precision,dimension(3)  ::d12, d34  !d12=a1-a2, d34=d3-d4
  double precision  :: w(4) , dwdR(12,4)  , prod
  double precision  :: coef2(2)
  integer           :: i,j,k

  d12 = cgeom(:,atms(1))-cgeom(:,atms(2))
  d34 = cgeom(:,atms(3))-cgeom(:,atms(4))

  f = dot_product(d12,d34)
  b(1:3)   = d34
  b(4:6)   =-d34
  b(7:9)   = d12
  b(10:12) =-d12

  w12 = 0d0
  coef2(1)  = 1d0
  coef2(2)  = 0d0
  call calcwij(0,atms(1),atms(2),coef2,cgeom,r12,w12)
  if(abs(r12)<1D-30) r12=sign(1D-30,r12)
  f = f/r12
  b = b/r12 -f/r12*w12
  
  w34 = 0d0
  call calcwij(0,atms(3),atms(4),coef2,cgeom,r34,w34)
  if(abs(r34)<1D-30) r34=sign(1D-30,r34)
  f = f/r34
  b = b/r34
  b(7:12) = b(7:12) - f/r34*w34(1:6)

  k=1
  dwdR = 0d0
  if(scal/=0)then
    do i=1,2
      do j=3,4
        call calcwij(abs(scal),atms(i),atms(j),coef,cgeom,w(k),w12)
        if(abs(w(k))<1d-30)w(k)=sign(1D-30,w(k))
        dwdR(i*3-2:i*3,k) = w12(1:3) 
        dwdR(j*3-2:j*3,k) = w12(4:6) 
        k=k+1
      end do!j
    end do!i
  end if

  ! perform scaling
  !  0   : unscaled
  ! >0   : scaling using products of scaled dists 
  ! <0   : scaling using harmonic mean of scaled dists
  select case(scal)
      case(0)
        fval = f
        bval = b
      case(1:)
        prod = product(w)
        fval = f*prod
        bval = b*prod +  fval*  &
                (dwdR(:,1)/w(1)+dwdR(:,2)/w(2)+dwdR(:,3)/w(3)+dwdR(:,4)/w(4))
      case(:-1)
        prod = 1/sum(1/w)
        fval = 4*f*prod
        bval = 4*b*prod + fval*prod*&
                (dwdR(:,1)/w(1)**2+dwdR(:,2)/w(2)**2+dwdR(:,3)/w(3)**2+dwdR(:,4)/w(4)**2)
      case default
  end select! case(scal)
END SUBROUTINE
