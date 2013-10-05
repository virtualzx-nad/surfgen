module opttools 
 implicit none
 contains


!---print out geometry info ---
subroutine analysegeom(natoms,geom)
  implicit none
  integer, intent(in)   :: natoms
  double precision,intent(in)  ::  geom(3,natoms)
  double precision, parameter  ::  bohr2ang=0.529177249d0
  integer   ::  i,j,k,l
  double precision  ::  distmat(natoms,natoms), TLen = 2.2D0, d,d1(3),d2(3),d3(3),cpd(3)
  double precision,external  ::   dnrm2
  print *,"Cartesian Geometries in Atomic Units"
  print "(2x,3F15.10)",geom
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

  print "(/,A)","   Atom1   Atom2   Atom3   Atom4   Out-of-Plane Angle (Degrees)"
  do i=1,natoms
    do j=1,natoms
     if(j==i .or. distmat(i,j)>TLen)cycle
     d1 = geom(:,j)-geom(:,i)
     d1 = d1/dnrm2(3,d1,1)
     do k=j+1,natoms
       if(k==i .or. distmat(i,k)>TLen)cycle
       d2 = geom(:,k)-geom(:,i)
       d2 = d2/dnrm2(3,d2,1)
       do l=k+1,natoms
         if(l==i .or. distmat(i,l)>TLen)cycle
         d3 = geom(:,l)-geom(:,i)
         d3 = d3/dnrm2(3,d3,1)
         cpd(1) = d1(3)*(d2(2)-d3(2))+d2(3)*d3(2)-d2(2)*d3(3)+d1(2)*(d3(3)-d2(3))
         cpd(2) =-d2(3)*d3(1)+d1(3)*(d3(1)-d2(1))+d1(1)*(d2(3)-d3(3)) +d2(1)*d3(3)
         cpd(3) = d1(2)*(d2(1)-d3(1))+d2(2)*d3(1)-d2(1)*d3(2)+d1(1)*(d3(2)-d2(2))
         print "(2x,4(I5,3x),F12.4)",I,J,K,L, 90/Acos(0d0)* &
           asin((-d1(3)*d2(2)*d3(1)+d1(2)*d2(3)*d3(1)+d1(3)*d2(1)*d3(2)       &
                -d1(1)*d2(3)*d3(2)-d1(2)*d2(1)*d3(3)+d1(1)*d2(2)*d3(3))/      &
               dnrm2(3,cpd,1))   
       end do! l
     end do!k
    end do !j
  end do   !i    
end subroutine analysegeom


!---calculate harmonic frequencies from hessian matrix
subroutine getFreq(natoms,masses,hess,w)
  implicit none
  integer,intent(in)          :: natoms
  double precision,intent(in) :: masses(natoms),hess(3*natoms,3*natoms)
  double precision,intent(out):: w(3*natoms)
  double precision,  parameter  :: amu2au=1.822888484514D3,au2cm1=219474.6305d0
  integer  :: i,j
  double precision  :: sqrm,  hmw(3*natoms,3*natoms), tmp(1)
  double precision,dimension(:),allocatable :: WORK
  integer,dimension(:),allocatable :: IWORK
  integer           :: LIWORK, LWORK, itmp(1),INFO
  ! convert hessian into mass weighed coordinates
  hmw = hess/amu2au
  do i=1,natoms
    sqrm = sqrt(masses(i))
    do j=i*3-2,i*3
      hmw(j,:)=hmw(j,:)/sqrm
      hmw(:,j)=hmw(:,j)/sqrm
    end do
  end do 
  !calculate eigenvalues of hmw
  call DSYEVD('N','U',3*natoms,hmw,3*natoms,w,tmp,-1,itmp,-1,INFO)
  if(info/=0)print *,"DSYEVD allocation investigation failed.  info=",info
  LWORK = int(tmp(1))
  LIWORK= itmp(1)
  allocate(WORK(LWORK))
  allocate(IWORK(LIWORK))
  
  call DSYEVD('N','U',3*natoms,hmw,3*natoms,w,WORK,LWORK,IWORK,LIWORK,INFO)
  if(info/=0)print *,"DSYEVD failed.  info=",info
  do i=1,3*natoms
    if(w(i)<0)then
      w(i) = -sqrt(-w(i))*au2cm1
    else 
      w(i) = sqrt(w(i))*au2cm1
    end if
  end do
end subroutine getFreq


!---calculate hessian at a certain geometry
subroutine calcHess(natoms,cgeom,nstate,istate,stepsize,hessian,skip)
  implicit none
  integer, intent(in)          :: natoms, nstate,istate
  logical, intent(in),optional :: skip(natoms)
  double precision,intent(in)  :: stepsize
  double precision,intent(in)  :: cgeom(3*natoms)
  double precision,intent(out) :: hessian(3*natoms,3*natoms)
  double precision   ::  mdif
  
  integer   ::   i,  j, skipdisp(natoms*3)
  double precision  ::  dispgeom(3*natoms), dgrd(3*natoms)
  real*8    ::  h(nstate,nstate),cg(3*natoms,nstate,nstate),dcg(3*natoms,nstate,nstate),e(nstate)
  skipdisp=.false.
  if(present(skip))then
    do i=1,natoms
      if(skip(i))skipdisp(i*3-2:i*3)=.true.
    end do
  end if
  hessian = 0d0
  do i=1,3*natoms
    if(skipdisp(i))cycle
    dispgeom=cgeom
    dispgeom(i)=dispgeom(i) - stepsize
    call EvaluateSurfgen(dispgeom,e,cg,h,dcg)
    dgrd  =-cg(:,istate,istate)
    dispgeom=cgeom
    dispgeom(i)=dispgeom(i) + stepsize
    call EvaluateSurfgen(dispgeom,e,cg,h,dcg)
    dgrd = dgrd+cg(:,istate,istate)
    hessian(i,:)= dgrd/2/stepsize
  end do!o=1,3*natoms
  do i=1,3*natoms
    if(skipdisp(i))hessian(:,i)=0d0
  end do
  mdif = maxval(abs(hessian-transpose(hessian)))
  if(mdif>1d-5)print *,"maximum hermiticity breaking : ",mdif
  hessian = (hessian+transpose(hessian))/2
end subroutine calcHess


!----search for minimum on adiabatic surfaces 
subroutine findmin(natoms,nstate,cgeom,isurf,maxiter,shift,Etol,Stol)
  implicit none
  integer, intent(in)                                 ::  natoms,isurf,maxiter,nstate
  double precision,dimension(3*natoms),intent(inout)  ::  cgeom
  double precision,intent(in)                         ::  shift,Etol,Stol

  real*8    ::  h(nstate,nstate),cg(3*natoms,nstate,nstate),dcg(3*natoms,nstate,nstate),e(nstate)
  double precision,dimension(3*natoms)  ::  grad, b1, b2, w
  double precision,dimension(3*natoms,3*natoms)  :: hess,hinv
  double precision,dimension(:),allocatable :: WORK
  integer,dimension(:),allocatable :: IWORK
  integer           :: LIWORK, LWORK, itmp(1),INFO  
  integer  ::  iter  , i
  double precision            :: nrmG, nrmD, tmp(1)
  double precision, external  :: dnrm2
  double precision,  parameter  :: amu2au=1.822888484514D3,au2cm1=219474.6305d0
  double precision, parameter   :: MAXD = 1D-2

  ! initialize work spaces
  call DSYEVD('V','U',3*natoms,hess,3*natoms,w,tmp,-1,itmp,-1,INFO)
  if(info/=0)print *,"DSYEVD allocation failed.  info=",info
  LWORK = int(tmp(1))
  LIWORK= itmp(1)
  allocate(WORK(LWORK))
  allocate(IWORK(LIWORK))

  print "(A,I4,A,I4,A)","Searching for minimum on surface ",isurf," in ",maxiter," iterations."
  do iter=1,maxiter
     call EvaluateSurfgen(cgeom,e,cg,h,dcg)
     grad = cg(:,isurf,isurf)
     nrmG=dnrm2(3*natoms,grad,1)
     call calcHess(natoms,cgeom,nstate,isurf,1D-4,hess)
     hinv = hess
     ! invert the hessian
     call DSYEVD('V','U',3*natoms,hess,3*natoms,w,WORK,LWORK,IWORK,LIWORK,INFO)
     if(info/=0)print *,"DSYEVD failed.  info=",info
     ! hess. w. hess^T = hess_old
     ! so x= hess_old^-1.b = hess. w^-1. hess^T. b
     ! b1= hess^T.g       =>
     call DGEMV('T',3*natoms,3*natoms,1d0,hess,3*natoms,grad,1,0d0,b1,1)
     ! b1' = w^-1*b1
     do i=1,3*natoms
      if(abs(w(i))<shift)then
         b1(i)=dble(0)
      else!
         b1(i)=b1(i)/w(i)
      end if!
     end do!
     ! b2 = hess.b1'
     call DGEMV('N',3*natoms,3*natoms,1d0,hess,3*natoms,b1,1,0d0,b2,1)
     ! cgeom' = cgeom - H^-1.g
     nrmD=dnrm2(3*natoms,b2,1)
     if(nrmD>maxD)then
       b2=b2/nrmD*maxD
       nrmD=maxD
     end if
     print 1000,iter,e(isurf)*au2cm1,nrmG,nrmD
     cgeom = cgeom - b2
     if(nrmG<Etol.and.nrmD<Stol)then
       print *,"Optimization converged"
       return
     end if
  end do 
1000 format("Iteration ",I4,": E=",F20.4,", |Grd|=",E12.5,", |Disp|=",E12.5)
end subroutine findmin
end module opttools 

! problems:  wrong bound for array r() in subroutine 
! rx1 and rx2 have single precision default implicit type
! variable in subroutine set as intent(INOUT) as default and r(6) got switched every run
! inappropriate bound for grid

program findcp

 use opttools 
implicit none

  integer,parameter :: natm = 13, nst=4
! gmin is the global min
  real*8  :: gmin(39)=(/	-1.71199414,	0.00000000,	-0.00354718,	&
				-0.38703243,	0.00000000,	-2.28591809,	&
				2.25530537,	0.00000000,	-2.27653197,	&
				3.58567656,	0.00000000,	-0.00071621,	&
				2.24421854,	0.00000000,	2.27787307,	&
				-0.39002962,	0.00000000,	2.28562351 ,	&
				-4.25258741,	0.00000000,	0.10849745,	&
				-1.40667035,	0.00000000,	-4.06095466,	&
				3.26109586,	0.00000000,	-4.05562411,	&
				5.62831685,	0.00000000,	0.00217672,	&
				3.24992294,	0.00000000,	4.05714529,	&
				-1.44534888,	0.00000000,	4.03356585,	&
				-4.94259551,	0.00000000,	-1.58803693 	/)

  real*8  :: dis(39)=(/	-2.18996537,	0.00000000,	0.00000000,	&
			-0.74939941,	0.00000000,	2.34677175,	&
			1.85640018,	0.00000000,	2.31756577,	&
			3.18392345,	0.00000000,	0.00000000,	&
			1.85640018,	0.00000000,	-2.31756577,	&
			-0.74939941,	0.00000000,	-2.34677175,	&
			-4.51954534,	0.00000000,	0.00000000,	&
			-1.80877667,	0.00000000,	4.09382748,	&
			2.90247794,	0.00000000,	4.07428219,	&
			5.22705676,	0.00000000,	0.00000000,	&
			2.90247794,	0.00000000,	-4.07428219,	&
			-1.80877667,	0.00000000,	-4.09382748,	&
			30.00000000,	0.00000000,	0.00000000	/) 
! the minimum energy intersection
  real*8  :: mex01(39) =(/	-1.74652425,	0.00000000,	0.07042240,	&
				-0.39925857,	0.00000000,	-2.23489349,	&
				2.22489298,	0.00000000,	-2.26280385,	&
				3.59825936,	0.00000000,	-0.00272005,	&
				2.29561690,	0.00000000,	2.30445676,	&
				-0.33931315,	0.00000000,	2.37489408,	&
				-4.21991527,	0.00000000,	0.09775378,	&
				-1.47344284,	0.00000000,	-3.96658552,	&
				3.19890493,	0.00000000,	-4.05833294,	&
				5.64292066,	0.00000000,	-0.04200458,	&
				3.36128528,	0.00000000,	4.05356284,	&
				-1.30715914,	0.00000000,	4.17552016,	&
				-5.73132025,	0.00000000,	-3.00077829	/)
  real*8  :: mex12(39) =(/	-1.71136176,	0.00000000,	-0.00965084,	&
				-0.36782782,	0.00000000,	-2.36235648,	&
				2.22626539,	0.00000000,	-2.32736351,	&
				3.55065740,	0.00000000,	0.00664434,	&
				2.22568282,	0.00000000,	2.33210381,	&
				-0.37264378,	0.00000000,	2.35055752,	&
				-4.11963171,	0.00000000,	0.11212956,	&
				-1.42605725,	0.00000000,	-4.10703334,	&
				3.26775400,	0.00000000,	-4.08237091,	&
				5.59261842,	0.00000000,	0.00434343,	&
				3.26037605,	0.00000000,	4.08939387,	&
				-1.46705199,	0.00000000,	4.07053107,	&
				-4.95464207,	0.00000000,	-1.57906004	/)

  double precision,parameter :: au2cm1=219474.6305d0
  double precision, parameter  ::  bohr2ang=0.529177249d0

  double precision  ::  masses(natm),  hess(natm*3,natm*3), w(3*natm)
  double precision  ::  cgeom(3*natm),hvec(3*natm)
  integer  :: i,j, isurf
  logical :: skip(natm)= .false.
  double precision  ::  di, dj, x,y, r, r0, dr,dh
  real*8    :: h(nst,nst),cg(3*natm,nst,nst),dcg(3*natm,nst,nst),e(nst),evec(nst,nst)
  double precision,dimension(3)::d1,d2,d3,cpd
  double precision,external :: dnrm2
  character*9,dimension(4)  ::  fnames

  call initPotential()

! Plot branching space
  fnames(1)  = "surf1.csv"
  fnames(2)  = "surf2.csv"
  fnames(3)  = "surf3.csv"
  fnames(4)  = "surf4.csv"
  d1 = mex12(37:39)-mex12(19:21)
  r0 = dnrm2(3,d1,1)
  d1 = d1/r0
  dr = 6d-2
  dh = 1d-2
  print *,"r0=",r0*bohr2ang
  do isurf=1,4
   open(unit=123+isurf,file=fnames(isurf),access='sequential',action='write')
   do i=-5 ,40 
      cgeom = mex12
      r = r0+i*dr
      cgeom(37:39) = cgeom(19:21)+r*d1
      call EvaluateSurfgen(cgeom,e,cg,h,dcg)
      call getEvec(evec)
      write(unit=123+isurf,fmt=1000),r*bohr2ang,0d0,e(isurf)*au2cm1,evec(:,isurf)
      hvec = dcg(:,2,3)
      hvec = hvec / dnrm2(3*natm,hvec,1)
      do j=1,35
        cgeom = cgeom+hvec*dh
        call EvaluateSurfgen(cgeom,e,cg,h,dcg)
        call getEvec(evec)
        write(unit=123+isurf,fmt=1000),r*bohr2ang,dh*j,e(isurf)*au2cm1,evec(:,isurf)
      end do
   end do!i
   close(123+isurf)
  end do

stop

  masses = (/	12.,12.,12.,12.,12.,12.,15.99491464,	&
		1.00782504,1.00782504,1.00782504,1.00782504,1.00782504,1.00782504 /) 

1000 format(2(F12.7,', '),F20.5,4(',',F12.8))
  cgeom = gmin
  isurf = 1
!  do i=1,300
!    call EvaluateSurfgen(cgeom,e,cg,h,dcg)
!         d1 = cgeom(4:6)-cgeom(1:3)
!         d2 = cgeom(7:9)-cgeom(1:3)
!         d3 = cgeom(10:12)-cgeom(1:3)
!         d3 = d3/dnrm2(3,d3,1)
!         cpd(1) = d1(3)*(d2(2)-d3(2))+d2(3)*d3(2)-d2(2)*d3(3)+d1(2)*(d3(3)-d2(3))
!         cpd(2) =-d2(3)*d3(1)+d1(3)*(d3(1)-d2(1))+d1(1)*(d2(3)-d3(3)) +d2(1)*d3(3)
!         cpd(3) = d1(2)*(d2(1)-d3(1))+d2(2)*d3(1)-d2(1)*d3(2)+d1(1)*(d3(2)-d2(2))
!    print "(I4,4F14.6)",i,e*au2cm1,(e(2)-e(1))*au2cm1, 90/Acos(0d0)* &
!           asin((-d1(3)*d2(2)*d3(1)+d1(2)*d2(3)*d3(1)+d1(3)*d2(1)*d3(2)       &
!                -d1(1)*d2(3)*d3(2)-d1(2)*d2(1)*d3(3)+d1(1)*d2(2)*d3(3))/      &
!               dnrm2(3,cpd,1))   
!    cgeom = cgeom-cg(:,2,2)*1D-3
!    if((e(2)-e(1))*au2cm1<1)exit
!  end do
  print "(/,A)","-------------- Initial Geometry ------------"
  call analysegeom(natm,cgeom)
  print "(/,A)","----------- Geometry Optimizations ---------"
  call findmin(natm,nst,cgeom,isurf,100,1d-2,1d-10,1d-5)
  print "(/,A)","--------------- Final Geometry -------------"
  call analysegeom(natm,cgeom)
  print "(/,A)","------------ Harmonic Frequencies ----------"
  call calcHess(natm,cgeom,nst,isurf,1D-3,hess,skip)
  call getFreq(natm,masses,hess,w)
  do i=1,3*natm
    print "(I5,F12.2)",i,w(i)
  end do
end
