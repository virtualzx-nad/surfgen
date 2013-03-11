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
!reorder oop atoms to canonical order and returns the parity of perm
SUBROUTINE reorderOOP(oldOOP,newOOP,parity)
  implicit none
  INTEGER,DIMENSION(4),INTENT(IN)   ::  oldOOP
  INTEGER,DIMENSION(4),INTENT(OUT)  ::  newOOP
  INTEGER,INTENT(OUT)               ::  parity
  integer          ::  loc, tmp, i
  parity=1
  newOOP=oldOOP
  do i=4,3,-1
    loc   = maxloc(newOOP(2:i),1)+1
    if(loc/=i)then
      tmp   = newOOP(loc)
      newOOP(loc)  =newOOP(i)
      newOOP(i)  = tmp
      parity  = -parity
    end if
  end do
END SUBROUTINE reorderOOP

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


!
! Transform a gradient in cartesian coordinates to one in internal
! coordinates
!
SUBROUTINE cart2intgrad(cgrad,bmat,igrad)
  use hddata, only:  ncoord
  use progdata, only: natoms
  IMPLICIT NONE
  DOUBLE PRECISION,dimension(3*natoms),INTENT(IN)         :: cgrad
  DOUBLE PRECISION,dimension(ncoord),INTENT(INOUT)        :: igrad
  DOUBLE PRECISION,dimension(ncoord,3*natoms),INTENT(IN)  :: bmat

  DOUBLE PRECISION,dimension(ncoord,3*natoms)             :: invB
  DOUBLE PRECISION,dimension(3*natoms,1)      :: tmp

  tmp(1:3*natoms,1)=cgrad
  call ginv(ncoord,3*natoms,bmat,ncoord,invB,1D-7)
  ! igrad=invB.cgrad
  igrad=reshape(matmul(invB,tmp),(/ncoord/))
end subroutine cart2intgrad

!builds b matrix then calculate eigenvectors of B^T.B
SUBROUTINE evecBTB(dim,cgeom,evals,evecs,cweight)
  USE progdata, only: natoms
  USE hddata, only: ncoord
  IMPLICIT NONE
  INTEGER,intent(IN)                                   :: dim
  DOUBLE PRECISION,dimension(dim),INTENT(IN)           :: cgeom
  DOUBLE PRECISION,dimension(natoms),INTENT(IN)        :: cweight
  DOUBLE PRECISION,dimension(dim),INTENT(OUT)          :: evals
  DOUBLE PRECISION,dimension(dim,dim),INTENT(OUT)      :: evecs
  
  double precision, dimension(ncoord)                  :: igeom
  double precision, dimension(ncoord,3*natoms)         :: bmat
  double precision, dimension(3*natoms,3*natoms)       :: btb
  double precision,dimension(45*natoms*natoms)  :: scr
  integer                      ::  INFO,i

  call buildWBmat(cgeom(:3*natoms),igeom,bmat,.false.)
  do i=1,natoms
    bmat(:,i*3-2:i*3)=bmat(:,i*3-2:i*3)/sqrt(cweight(i))
  end do
  btb=matmul(transpose(bmat),bmat)
  call DSYEV('V','U',3*natoms,btb,3*natoms,evals(:3*natoms),scr,45*natoms*natoms,INFO)
  evecs(:3*natoms,:3*natoms)=btb
END SUBROUTINE 
!
! Builds a b-matrix given a set of internal coordinates definitions
!
!
SUBROUTINE buildWBmat(cgeom,igeom,bmat)
  use hddata, only: ncoord
  use progdata, only: natoms,nrij,nout,coordmap,CoordSet
  IMPLICIT NONE
  DOUBLE PRECISION,dimension(3*natoms),INTENT(IN)           :: cgeom
  DOUBLE PRECISION,dimension(ncoord),INTENT(OUT)            :: igeom
  DOUBLE PRECISION,dimension(ncoord,3*natoms),INTENT(OUT)   :: bmat
  INTEGER                                       :: i,j,offs,m,n

! bval     Derivative of a unscaled coordinate with respect to cartesian 
!            coordinates of the four reference atoms.
  DOUBLE PRECISION     ::  bval(12)
  DOUBLE PRECISION     ::  w(3), dwdR(3,3*natoms),coef(2),s(3),ss,fs

  do i=1,ncoord
    m=coordmap(i,1) !index of set
    n=coordmap(i,2) !index in set

    select case(CoordSet(m)%Type)

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     Plain or scaled Rij coordinates
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

      case(0)   
        
        call calcwij(CoordSet(m)%Scaling,CoordSet(m)%coord(n,1),CoordSet(m)%coord(n,2),&
                      CoordSet(m)%coef,cgeom,igeom(i),bmat(i,:))

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     OOP ANGLE
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      case(-1)  !oop

        select case(CoordSet(m)%Scaling)

          case(0)  ! Reciprocal scaling   Det/(Product[rij])^a
            bmat(i,:) = dble(0)
            call  oop(natoms,CoordSet(m)%coord(n,:),cgeom,igeom(i),bval,coordset(m)%coef(1))
            do j=1,4
              offs = 3*(CoordSet(m)%coord(n,j)-1)
              bmat(i,offs+1:offs+3) =  bval(3*j-2:3*j)
            end do!j=1,4
            
          case default  ! Exponential scaling using one of the linear scaling types.
            bmat(i,:) = dble(0)
            call oop2(natoms,CoordSet(m)%coord(n,:),cgeom,igeom(i),bval,coordset(m)%coef(1),&
                                                    CoordSet(m)%Scaling,coordset(m)%coef(2))
            do j=1,4
              offs = 3*(CoordSet(m)%coord(n,j)-1)
              bmat(i,offs+1:offs+3) =  bval(3*j-2:3*j)
            end do!j=1,4

        end select !case(CoordSet(m)%Scaling)
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     Special OOP angles
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      case(-2)  !special oop angles
        bmat(i,:) = dble(0)
        call oop3(natoms,CoordSet(m)%coord(n,:),cgeom,igeom(i),bval,&
                coordset(m)%coef(1),CoordSet(m)%Scaling,coordset(m)%coef(2))
        do j=1,4
          offs = 3*(CoordSet(m)%coord(n,j)-1)
          bmat(i,offs+1:offs+3) =  bval(3*j-2:3*j)
        end do!j=1,4

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!     Scaled bond angle
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      case (1) ! angle
        bmat(i,:) = dble(0)
        call bend(natoms,CoordSet(m)%coord(n,:),cgeom,igeom(i),bval,&
        coordset(m)%coef(1),CoordSet(m)%Scaling,coordset(m)%coef(2))
        do j=1,3
            offs = 3*(CoordSet(m)%coord(n,j)-1)
            bmat(i,offs+1:offs+3) =  bval(3*j-2:3*j)
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
  DOUBLE PRECISION  :: denom, dnorm2(6), dvecs(6,3),&
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
  dvecs=matmul(sgns,geom)
  do i=1,6
    dnorm2(i)=dot_product(dvecs(i,:),dvecs(i,:))
  end do

! calculate value of scaled OOP angle and its derivatives
  denom=exp(log(product(dnorm2))/2*scale)
  if(denom==0)stop "OOP angle cannot be defined when atoms coincide."
  qval=det3(dvecs(1:3,:))/denom
  do i=1,3
    do j=1,3
      ddvec(i,j)=sum(sgns(:,i)*dvecs(:,j)/dnorm2(:))
    end do
  end do
  bval(1:3)=cross(dvecs(2,:),dvecs(3,:))/denom-qval*scale*ddvec(1,:)
  bval(4:6)=cross(dvecs(3,:),dvecs(1,:))/denom-qval*scale*ddvec(2,:)
  bval(7:9)=cross(dvecs(1,:),dvecs(2,:))/denom-qval*scale*ddvec(3,:)
  bval(10:12)=-bval(1:3)-bval(4:6)-bval(7:9)

! calculate the normalization factor
  factor = 2D0 ** (scale*3-3)*3D0**(scale*3)
  qval = qval * factor
  bval = bval * factor

CONTAINS

  ! this function calculates the value of a 3x3 determinant
  FUNCTION det3(m)
    IMPLICIT NONE
    double precision ::  det3, m(3,3)
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
  if(vnorm.eq.0)then
   denom = dble(1)
  else
   denom = vnorm
  endif

  vec=vec/denom

  return
end SUBROUTINE dvec


!
!
!
SUBROUTINE cprod(vec1,vec2,cpvec,cpnorm)
  IMPLICIT NONE
  DOUBLE PRECISION,dimension(3),INTENT(IN)    :: vec1,vec2
  DOUBLE PRECISION,INTENT(INOUT)              :: cpnorm
  DOUBLE PRECISION,dimension(3),INTENT(INOUT) :: cpvec

  INTEGER                                     :: i
  DOUBLE PRECISION                            :: denom

  cpvec(1) = vec1(2)*vec2(3) - vec1(3)*vec2(2)
  cpvec(2) = vec1(3)*vec2(1) - vec1(1)*vec2(3)
  cpvec(3) = vec1(1)*vec2(2) - vec1(2)*vec2(1)

  cpnorm = 0
  do i = 1,3
   cpnorm = cpnorm + cpvec(i)**2
  enddo

  cpnorm = Sqrt(cpnorm)
  if(cpnorm.eq.0) then
   denom = 1.
  else
   denom = cpnorm
  endif

  do i = 1,3
   cpvec(i) = cpvec(i)/denom
  enddo

  return
end SUBROUTINE cprod


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
  DOUBLE PRECISION,dimension(3*natoms),INTENT(OUT)  :: dwdR
  INTEGER                                       :: i,j,offs

! fval     Value of current coordinate in its unscaled form
! fbmat    Cartesian gradients of current coordinate in its unscaled form
! bval     Derivative of a unscaled coordinate with respect to cartesian 
!            coordinates of the four reference atoms.

  DOUBLE PRECISION     ::  fval, bval(12), g, fbmat(3*natoms)

!---Plain or scaled Rij coordinates--->>>

  fbmat = dble(0)
  call stre(natoms,a1,a2,cgeom,fval,bval)
  offs = 3*(a1-1)
  fbmat(offs+1:offs+3) =  bval(1:3)
  offs = 3*(a2-1)
  fbmat(offs+1:offs+3) =  bval(4:6)

  select case(scaling) 
    !  no scalings
    case(0)
      w    = fval
      dwdR = fbmat

    !  Morse functions   Exp(-c1*(r-c2))
    case(1)
      w    =  exp(-coef(1)*(fval-coef(2)))
      dwdR = -fbmat*w*coef(1)

    !  Gaussian functions Exp(-c1*(r-c2)**2)
    case(2)
      w    =  exp(-coef(1)*(fval-coef(2))**2) 
      dwdR = -fbmat*w*coef(1)*2*(fval-coef(2))

    !  Screened Columb potential Exp(-c1*(r-c2))/r     (Yukawa, leading term)
    case(3)
      g    =  exp(-coef(1)*(fval-coef(2)))
      w    =  g/(fval+1D-40)
      dwdR = -fbmat*w*(coef(1)+1/(fval+1D-40))

    !  Long range term of screened Columb  Exp(-r/c2)*(r/c2)**c1
    case(4)
      w    = exp(1-fval/coef(2))*(fval/coef(2))**coef(1)
      dwdR = fbmat*w*(coef(1)/fval-1/coef(2))        

    !  Lennard Jones functions (c2/r)**c1
    case(5)
      w    =  (coef(2)/fval)**coef(1)
      dwdR = -fbmat*w*coef(1)/fval

    !  Shifted(chasmless) Yukawa exp(-c1*(r-c2))/(r+c2)
    case(6)
      g    =  exp(-coef(1)*(fval-coef(2)))
      w    =  g/(fval+coef(2))
      dwdR = -fbmat*w*(coef(1)+1/(fval+coef(2)))

   
    case default
      print *,"scaling = ",scaling
      stop "Unsupported bond distance scaling method."

  end select !case(scaling)
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
                       geom(4,3), ddvec(3,3),          &
                       coef(2), dwdR(6,12), w(6),  denom, ddenrm
  DOUBLE PRECISION, parameter :: vsmall = 1D-50
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
  do i=1,3
   do j=i+1,4
     call calcwij(tp,i,j,coef,cgeom,w(k),dwdR(k,:))
     k = k+1
   end do !j
  end do


! extract the geometry of the 4 involved atoms
  do i=1,4
    geom(i,:)=cgeom(atms(i)*3-2:atms(i)*3)
  end do

! calculate displacement vectors between these atoms
  dvecs=matmul(sgns,geom)

! calculate value of unscaled OOP angle and its derivatives
  fval=det3(dvecs(1:3,:))
  do i=1,3
    do j=1,3
      ddvec(i,j)=sum(sgns(:,i)*dvecs(:,j)/dnorm2(:))
    end do
  end do
  fbmat(1:3)=cross(dvecs(2,:),dvecs(3,:))
  fbmat(4:6)=cross(dvecs(3,:),dvecs(1,:))
  fbmat(7:9)=cross(dvecs(1,:),dvecs(2,:))
  fbmat(10:12)=-fbmat(1:3)-fbmat(4:6)-fbmat(7:9)

! calculate the scaled oop
  if(sctype<0)then
    denom = sum(w)**scale
    qval = fval * factor / denom
    bval = fbmat* factor / denom - qval*scale/sum(w)*&
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
  FUNCTION det3(m)
    IMPLICIT NONE
    double precision ::  det3, m(3,3)
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

  INTEGER           :: i, j, k, sgns(6,4)
  DOUBLE PRECISION  :: dnorm2(6), dvecs(6,3), s(3),     &
                       geom(4,3), ddvec(3,3), ss,       &
                       coef(2), dwdR(3,12), w(3), dwdR2(3,12), w2(3)
! eprod = product of the scaling exponentials
! fval  = value of unscaled oop coordinate (the triple product)
! fbmat = derivatives of unscaled OOP coordinate
  DOUBLE PRECISION  :: eprod, fval, fbmat(12),efactor,a,p,p1

! The oriented direction vectors between all 6 atom pairs
  sgns(1,:)=(/ 1, 0, 0,-1/)
  sgns(2,:)=(/ 0, 1, 0,-1/)
  sgns(3,:)=(/ 0, 0, 1,-1/)
  sgns(4,:)=(/ 1, 0,-1, 0/)
  sgns(5,:)=(/ 0, 1,-1, 0/)
  sgns(6,:)=(/ 1,-1, 0, 0/)

! calculate the three scaling factors and their derivatives
! scaled distances for the three bonds A1-A2, A1-A3 and A1-A4 are used
  do i=2,4
! scaling mode=0   (unscaled).  unscaled distances are produced here
     call calcwij(int(0),int(1),i,coef,cgeom,w(i-1),dwdR(i-1,:))
  end do

! here are the scaled distances
  coef(1) =scale
  coef(2) =factor 
  if(sctype.gt.0)then
    do i=2,4
     call calcwij(sctype,int(1),i,coef,cgeom,w2(i-1),dwdR2(i-1,:))
    end do
  end if

! extract the geometry of the 4 involved atoms
  do i=1,4
    geom(i,:)=cgeom(atms(i)*3-2:atms(i)*3)
  end do

! calculate displacement vectors between these atoms
  dvecs=matmul(sgns,geom)

! calculate value of unscaled OOP coordinate and its derivatives
  fval=det3(dvecs(1:3,:))
  do i=1,3
    do j=1,3
      ddvec(i,j)=sum(sgns(:,i)*dvecs(:,j)/dnorm2(:))
    end do
  end do
  fbmat(1:3)=cross(dvecs(2,:),dvecs(3,:))
  fbmat(4:6)=cross(dvecs(3,:),dvecs(1,:))
  fbmat(7:9)=cross(dvecs(1,:),dvecs(2,:))
  fbmat(10:12)=-fbmat(1:3)-fbmat(4:6)-fbmat(7:9)

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
    a = exp(3*coef(1)*log(coef(2)))
    p=product(w)
    if(p<1D-10)p=1D-10
    p1 = exp(coef(1)*log(p))
    qval = fval*a/p1
    bval = fbmat*a/p1-coef(1)*qval/p*(dwdR(1,:)*w(2)*w(3)+dwdR(2,:)*w(1)*w(3)+dwdR(3,:)*w(1)*w(2))
  else
! multiply all 3 scaled length scaling functions
    eprod = product(w2)
    qval = fval*eprod
    bval = fbmat*eprod+qval*(w2(2)*w2(3)*dwdR2(1,:)+w2(1)*w2(3)*dwdR2(2,:)+w2(1)*w2(2)*dwdR2(3,:))
  end if

CONTAINS

  ! this function calculates the value of a 3x3 determinant
  FUNCTION det3(m)
    IMPLICIT NONE
    double precision ::  det3, m(3,3)
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
  double precision ::  p, dp(9), denom, tmp
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
!       Raw angle.
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!  obtained as acos(p).    d acos(p)= - dp / Sqrt(1-p^2)
    case ( 0 )
        qval = acos(p)
        bval = -dp/sqrt(1-p**2)

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!       Scaled cosine, p/( 1+exp[c1*(d1^2+d2^2-c2^2)] )
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!   q = p/( 1+exp[c1*(d1^2+d2^2-c2^2)] )
    case  (2)
        tmp  = exp(scale*(d12+d22-factor**2))
        qval = p/(1+tmp)
        bval = (  dp  +  scale*tmp*p/(1+tmp)*(dd12+dd22)  ) /(1+tmp)
  end select!case(sctype)
END SUBROUTINE BEND
