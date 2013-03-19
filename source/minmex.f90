!minmex Perform intersection and critical point searches
SUBROUTINE minmex()
  use progdata, only: nmin,nmex,minstates,mexstates,minstart,mexstart,&
                      h_update,xscale,sscale
  IMPLICIT NONE
  integer :: i
  do i=1,nmin
    CALL findcp(minstart(:,i),minstates(i),h_update(i))
  end do
  do i=1,nmex
    CALL findmex(mexstart(:,i),mexstates(i,:),xscale,sscale)
  end do
END SUBROUTINE

! FINDCP() try and find an stable point on ab initio surface
! A quasi-Newton method is used.  Hessians are constructed from Hd if
! requested, and inverse Hessian can be updated by BFGS update.
! When Hessian is not calculated, unitary matrix is used instead.
! Extrapolations are done in cartesian coordinates.  Singularities
! corresponding to translations and rotations will be removed when
! hessian is inverted.
! GDIIS will be implemented in the future
! geom      (input/output) DOUBLE PRECISION,dimension(3*natoms)
!           On entry geom contains the cartesian geometry to start the
!           quasi-Newton search.
!           On exit, geom contains the geometry of optimized critical
!           point if the search has converged, or the last geometry
!           if it has not converged yet.
! surf      (input) INTEGER
!           The surface to look for critical point
! update    (input) LOGICAL
!           Whether or not to perform BFGS update to inverse Hessian
SUBROUTINE findcp(geom,surf,update)
  use progdata, only : natoms,opttoler,optiter,maxstep,h_recal,AU2CM1,enforcepd
  use hddata, only: nstates
  IMPLICIT NONE
  DOUBLE PRECISION,intent(INOUT):: geom(3*natoms)
  INTEGER         ,intent(IN)   :: surf
  LOGICAL         ,intent(IN)   :: update
  double precision,dimension(3*natoms,3*natoms) :: hess,ihess
  double precision,dimension(3*natoms)          :: dgeom,grad,yTB,y
  double precision :: step,force,sTy,a,b,e(nstates)
  integer          :: i,h_age, & !age of hessian matrix
                      nneg,nzero  !num of negative/zero eval in hess
  double precision, external   :: DNRM2

  print 1000,surf
  print 1100,geom
  print 1200,update,h_recal,opttoler,optiter
  if(h_recal>0)then
    h_age=h_recal
  else!if(usehs)
    ihess=dble(0)
    do i=1,3*natoms
      ihess(i,i)=dble(1)
    end do
    h_age=1
  end if!(usehs)
  do i=1,optiter
    CALL cgrad(geom,surf,e,grad)
    force=DNRM2(3*natoms,grad,int(1))
    if(force<opttoler)exit
    if(h_age==h_recal)then
    ! recalculate hessian
      print 1007
      CALL chess(geom,surf,hess)
    ! invert hessian
      CALL DSYINV("U",3*natoms,hess,3*natoms,ihess,3*natoms,1D-9,nneg,nzero,enforcepd)
      h_age=0
      print 1001,nneg,nzero
    else if(update.and.i>1)then
      ! update inverse hessian using BFGS method
      !s_k=dgeom, y=grad-gradl
      y=grad-y
      yTB=matmul(y,ihess)
      sTy=dot_product(dgeom,y)
      a=(sTy+dot_product(yTB,y))/sTy**2
      b=-1/sTy
      CALL DSYR('U',3*natoms,a,dgeom,int(1),ihess,3*natoms)
      CALL DSYR2('U',3*natoms,b,dgeom,int(1),yTB,int(1),ihess,3*natoms)
      print 1003
    end if
    y=grad
    !calculate displacement vector using inverse Hessian: d=-H^-1.g
    CALL DSYMV('U',3*natoms,dble(-1),ihess,3*natoms,grad,int(1),&
                                               dble(0),dgeom,int(1))
    h_age=h_age+1
    step=DNRM2(3*natoms,dgeom,int(1))
    if(step>maxstep)then
      print 1002,step
      dgeom=dgeom/step*maxstep
      step=maxstep
    end if
    geom=geom+dgeom
    print 1005,i,e(surf)*AU2CM1,force,step
  end do
  CALL cgrad(geom,surf,e,grad)
  if(DNRM2(3*natoms,grad,int(1))<opttoler)then
    print 1104,i
  else
    print 1004,i
  end if
  print 1204
  print 1100,geom
  print 1304
  print 1101,e*AU2CM1
1000 format(/,2X,"Searching for Critical Point on Adiabatic Surface ",I4,/,&
            4X, "Startng Geometry: ")
1100 format(4X,3F16.8)
1101 format(4X,3F16.2)
1200 format(4X,"BFGS Update to Inverse Hessian:",L,4X,"Re-evaluate Hessian Every ",I4,"Steps",/,&
            4X,"Gradient Convergence Criteria",F8.4,2X,"Max Number of Iternations:",I4)
1001 format(4X,"Negative eigenvalues:",I4,2X," Singularities:",I4)
1002 format(4X,"Step size over cap. |deltaR|=",F10.4)
1003 format(4X,"BFGS update performed on inverse Hessian matrix.")
1004 format(2X,"Critical Point Not Found after ",I4," Iterations")
1104 format(2X,"Critical Point Found after ",I4," Iterations")
1204 format(4X,"Final Geometry:")
1304 format(4X,"Final Energies:")
1005 format(2X,"Iter",I4," : Energy=",F14.4," Force=",F12.8,"  Step=",F12.8)
1007 format(4X,"Reconstructing Hessian Matrix")
END SUBROUTINE

!---------------------------------------------------------------
! search for intersection
SUBROUTINE findmex(geom,dstates,xscale,sscale)
  use progdata, only: natoms,optiter,AU2CM1,printlvl,maxstep,degtoler,opttoler,enforcepd
  use hddata
  IMPLICIT NONE
  DOUBLE PRECISION,dimension(3*natoms),intent(INOUT) :: geom
  INTEGER,dimension(2),intent(IN)                 :: dstates
  DOUBLE PRECISION,intent(IN)                     :: xscale,sscale
  double precision,dimension(3*natoms)            :: g,h,s,dgeom,gT
  double precision    :: e(nstates),beta,sforce,deg,step,gnrm,hnrm
  double precision    :: grad(3*natoms,nstates,nstates),center(3)
  double precision,dimension(3*natoms-7,3*natoms-7)::hess1,hess2
  double precision    :: orthvec(8+3*natoms,3*natoms),disp(2+3*natoms)
  double precision,dimension(nstates,nstates)      :: givens
  integer   :: i,j,m,n,rank,nneg,nzero
  double precision,external :: DNRM2
  m=minval(dstates)
  n=maxval(dstates)
  if(n/=m+1)then
    print *, "optimization of intersection with more than 2 states is not implemented"
    stop
  end if!(dstates(2)/=dstates(1))
  if(printlvl>0)then
    print 1000,m,n,opttoler,optiter,xscale,sscale
    print 1100,geom
  end if
  !initialize the Givens rotation matrix to identity matrix
  givens=dble(0)
  do j=1,nstates
    givens(j,j)=dble(1)
  end do
  ! orthvec contains three translational and three rotational vectors, g and h
  ! vectors, and all 3*natoms unit Cartesian displacement vectors.
  ! It is used to construct a set of Cartesian displacements that is orthogonal to
  ! trans-rot and degeneracy lifting directions.
  orthvec=dble(0)
  ! construct translational vectors
  do i=1,3
    do j=1,natoms
      orthvec(i,j*3+i-3)=dble(1)
    end do
  end do
  disp=dble(0)
  do i=1,optiter
    !obtain diabatic Hd and coupling gradients
    CALL ccoupl(geom,e,grad)
    if(printlvl>0)print 1002,i,e(m)*AU2CM1,e(n)*AU2CM1
    deg=abs(e(m)-e(n))*AU2CM1
    !find 2-state diabatrization that orthogonalizes g and h
    CALL orthgh(3*natoms,(grad(:,m,m)-grad(:,n,n))/2,grad(:,m,n),beta)
    if(printlvl>1.and. abs(beta)>0.01)print 1003,beta
    givens(m,m)=cos(beta)
    givens(n,n)=givens(n,n)
    givens(m,n)=sin(beta)
    givens(n,m)=-givens(m,n)
    do j=1,3*natoms
      grad(j,:,:)=matmul(givens,matmul(grad(j,:,:),transpose(givens)))
    end do
    ! get g and h vectors in diabatic representation
    s=(grad(:,m,m)+grad(:,n,n))/2
    g=(grad(:,m,m)-grad(:,n,n))/2
    h=grad(:,m,n)
    gnrm=DNRM2(3*natoms,g,int(1))
    hnrm=DNRM2(3*natoms,h,int(1))
    if(printlvl>2)then
      print 1103,gnrm,hnrm,dot_product(g,h)
      print *,"      GVEC"
      print 1100,g
      print *,"      HVEC"
      print 1100,h
    end if
    !compute Hessian matrix of degenerate states in seam coordinates if
    !seam coordinates are not frozen
    ! shift the origin to geometric center
    do j=1,3
      center(j)=sum(geom(j:j+3*natoms-3:3))/natoms
    end do
    do j=1,natoms
      geom(j*3-2:j*3)=geom(j*3-2:j*3)-center
    end do
    if(printlvl>2)print 1004,DNRM2(int(3),center,int(1))
    ! construct rotational vectors
    orthvec(4:,:)=dble(0)
    orthvec(4,2:natoms*3-1:3)=-geom(3:natoms*3:3)
    orthvec(4,3:natoms*3:3)  = geom(2:natoms*3-1:3)
    orthvec(5,1:natoms*3-2:3)=-geom(3:natoms*3:3)
    orthvec(5,3:natoms*3:3)  = geom(1:natoms*3-2:3)
    orthvec(6,1:natoms*3-2:3)=-geom(2:natoms*3-1:3)
    orthvec(6,2:natoms*3-1:3)= geom(1:natoms*3-2:3)
    !copy g and h vectors to vector list
    orthvec(7,:)=g
    orthvec(8,:)=h
    do j=1,3*natoms
      orthvec(8+j,j)=dble(1)
    end do
    !obtain internal directions that are orthogonal to g and h
    CALL Schdmit(3*natoms+8,3*natoms,orthvec,3*natoms,1D-5,rank,9)
    if(rank/=3*natoms)&
              print *,"WARNING: rank/=3*natoms in findmex. rank=",rank
    disp(1)=-(e(m)-e(n))/2*xscale/gnrm*cos(beta*2)
    disp(2)=-(e(n)-e(m))/2*xscale/hnrm*sin(beta*2)
    if(abs(sscale)>0)then
      !obtain average Hessian of states m and n
      CALL thess(geom,m,rank-8,orthvec(9:rank,:),hess1,3*natoms-7)
      CALL thess(geom,n,rank-8,orthvec(9:rank,:),hess2,3*natoms-7)
      hess1=(hess1+hess2)/2
      !invert hessian
      CALL DSYINV('U',rank-8,hess1,3*natoms-7,hess2,3*natoms-7,1D-9,nneg,nzero,enforcepd)
      if(printlvl>1.and.(nneg>0.or. nzero>0))PRINT 1005,NNEG,NZERO
      !calculate gradient of average energy in transformed space
      CALL DGEMV('N',rank-8,3*natoms,dble(1),orthvec(9:rank,:),rank-8,&
                                                  s,int(1),dble(0),gT,int(1))
      sforce=DNRM2(rank-8,gT,int(1))
      !calculate displacement vector in transformed space
      CALL DSYMV('U',rank-8,-sscale,hess2,3*natoms-7,gT,int(1),dble(0),&
                                                             disp(3:),int(1))
    else
      disp(3:)=dble(0)
      sforce=dble(0)
    end if! (abs(sscale)>0)
    if(sforce<opttoler*sscale .and. deg<degtoler)exit
    !calculate displacement in original atom centered cartesians
    CALL DGEMV('T',rank-6,3*natoms,dble(1),orthvec(7:rank,:),rank-6,&
                                            disp,int(1),dble(0),dgeom,int(1))
    step=DNRM2(3*natoms,dgeom,int(1))
    if(step>maxstep)then
      if(printlvl>0)print 1007,step
      dgeom=dgeom/step*maxstep
      step=maxstep
    end if
    if(printlvl>0)print 1006,sforce,deg,step
    geom=geom+dgeom
    if(printlvl>2)then
      print 1008
      print 1100,geom
    end if
  end do!i=1,optiter
  if(printlvl>0)then
    if(sforce<opttoler*sscale .and. deg<degtoler)then
      print 1109,i
    else!(sforce<opttoler*sscale .and. deg<degtoler)
      print 1009,i
    end if!(sforce<opttoler*sscale .and. deg<degtoler)
    print 1010
    print 1100,geom
    print 1110
    print 1101,e*AU2CM1
  end if!(printlvl>0)
1000 format(/,2X,"Searching for MEX between States ",I4,I4,/,&
           4X,"Seam Force Convergence Tolerance",F10.4,"Max Number of Iterations",I4,&
           4X,"g-h Scaling Factor",F10.4,4x,"Seam Scaling Factor",F10.4,&
           4X,"Startng Geometry: ")
1100 format(4X,3F16.8)
1101 format(4X,3F16.2)
1002 format(4X,"Iteration ",I4," Adiabatic Energies:",2F16.2," cm-1")
1003 format(6X,"Adiabatic to Diabatic Rotation BETA=",F10.2)
1103 format(6X,"|G|=",F14.8,", |H|=",F14.8,", G.H=",F14.8)
1004 format(6X,"Norm of displacement of molecular center ",F16.8)
1005 format(6X,"Seam Hessian has ",I4," negative evals and ",I4," Singularities")
1006 format(6X,"Seam Force ",F10.4,4X,"Degeneracy ",F10.4,"cm-1",4X,"Step Len ",F9.4)
1007 format(6X,"Step size over cap. |deltaR|=",F10.4)
1008 format(6X,"Displaced geometry: ")
1009 format(2X,"Minimum Energy Intersection Not Found after ",I4," Iterations")
1109 format(2X,"Minimum Energy Intersection Found after ",I4," Iterations")
1010 format(4X,"Final Geometry:")
1110 format(4X,"Final Energies:")
END SUBROUTINE findmex

!---------------------------------------------------------------
! These subroutines calculates energies and/or energy gradients 
! and/or non-adiabatic couplings and/or hessian matrix of adiabatic
! electronic states, defined by the solving Schrodinger equations
! for Hd at the gemetry. Geometry, gradients and hessians are all 
! in cartesian coordinate.
!--------------------------------------------------------------
! CENERGY calculate energies of a state
! cgeom      (input) DOUBLE PRCISION(3*natoms)
!            cartesian geometry
! e          (output) DOUBLE PRECISION(nstates)
!            adiabatic energies
SUBROUTINE cenergy(cgeom,e)
  use progdata ,only: natoms
  use hddata
  IMPLICIT NONE
  DOUBLE PRECISION,dimension(3*natoms),intent(in) :: cgeom
  DOUBLE PRECISION,dimension(nstates),intent(out) :: e
  double precision,dimension(ncoord)              :: igeom
  double precision,dimension(ncoord,3*natoms)     :: bmat
  double precision,dimension(nstates,nstates)     :: hmat,evec
  integer                                         :: m,LWORK,LIWORK,INFO
  integer,dimension(nstates*2)                    :: ISUPPZ
  double precision,dimension(nstates*(nstates+26)):: WORK
  integer,dimension(nstates*10)                   :: IWORK

  LWORK  = nstates*(nstates+26)
  LIWORK = nstates*10
  call buildWbmat(cgeom,igeom,bmat,.false.)
  call makehmat(igeom,hmat)
  CALL DSYEVR('V','A','U',nstates,hmat,nstates,dble(0),dble(0),0,0,1D-12,m,e,&
             evec,nstates,ISUPPZ,WORK,LWORK,IWORK,LIWORK, INFO )
END SUBROUTINE cenergy

!------------------------------------------------------------
! CGRAD calculate gradient of a state
! cgeom      (input) DOUBLE PRCISION(3*natoms)
!            cartesian geometry
! surf       (input) INTEGER
!            The state to evaluate
! grad       (output) DOUBLE PRECISION(3*natoms)
!            Gradients of current state in cartesian coordinates. 
SUBROUTINE cgrad(cgeom,surf,e,grad)
  use progdata ,only: natoms
  use hddata
  IMPLICIT NONE
  DOUBLE PRECISION,dimension(3*natoms),intent(in) :: cgeom
  INTEGER,intent(in)                              :: surf
  DOUBLE PRECISION,dimension(3*natoms),intent(out):: grad
  double precision,dimension(ncoord)              :: igeom
  double precision,dimension(nstates),intent(out) :: e
  double precision,dimension(ncoord,3*natoms)     :: bmat
  double precision,dimension(nstates,nstates)     :: hmat,evec
  double precision,dimension(ncoord,nstates,nstates)     :: dhmat
  integer                                         :: m,k,LWORK,LIWORK,INFO
  integer,dimension(nstates*2)                    :: ISUPPZ
  double precision,dimension(nstates*(nstates+26)):: WORK
  integer,dimension(nstates*10)                   :: IWORK
  double precision                                :: igrad

  LWORK  = nstates*(nstates+26)
  LIWORK = nstates*10
  call buildWbmat(cgeom,igeom,bmat,.false.)
  call makehmat(igeom,hmat)
  CALL DSYEVR('V','A','U',nstates,hmat,nstates,dble(0),dble(0),0,0,1D-12,m,e,&
             evec,nstates,ISUPPZ,WORK,LWORK,IWORK,LIWORK, INFO )
  ! calculate energy gradients
  grad=dble(0)
  call makedhmat(igeom,dhmat)
  do m=1,ncoord
    igrad=dble(0)
    do k=1,nstates
      igrad=igrad+evec(k,surf)*sum(dhmat(m,:,k)*evec(:,surf))
    end do
    do k=1,3*natoms
      grad(k)=grad(k)+igrad*bmat(m,k)
    end do
  end do
END SUBROUTINE cgrad

!------------------------------------------------------------
! CCOUPL calculate coupling gradients between any 2 states in
! adiabatic basis
! cgeom      (input) DOUBLE PRCISION(3*natoms)
!            cartesian geometry
! hmat       (output) DOUBLE PRECISION(nstates,nstates)
!            Value of Hd matrix
! coup       (output) DOUBLE PRECISION(3*natoms,nstates,nstates)
!            Adiabatic coupling gradients in cartesian coordinates.
SUBROUTINE ccoupl(cgeom,e,coupl)
  use progdata ,only: natoms
  use hddata
  IMPLICIT NONE
  DOUBLE PRECISION,dimension(3*natoms),intent(in) :: cgeom
  DOUBLE PRECISION,dimension(3*natoms,nstates,nstates) &
                                      ,intent(out):: coupl
  DOUBLE PRECISION,dimension(nstates),intent(out) :: e
  double precision,dimension(nstates,nstates)     :: hmat
  double precision,dimension(ncoord)              :: igeom
  double precision,dimension(ncoord,3*natoms)     :: bmat
  double precision,dimension(nstates,nstates)     :: evec
  double precision,dimension(ncoord,nstates,nstates)     :: dhmat
  integer                                         :: m,k,LWORK,LIWORK,INFO
  integer,dimension(nstates*2)                    :: ISUPPZ
  double precision,dimension(nstates*(nstates+26)):: WORK
  integer,dimension(nstates*10)                   :: IWORK

  LWORK  = nstates*(nstates+26)
  LIWORK = nstates*10
  call buildWbmat(cgeom,igeom,bmat,.false.)
  call makehmat(igeom,hmat)
  CALL DSYEVR('V','A','U',nstates,hmat,nstates,dble(0),dble(0),0,0,1D-10,m,e,&
             evec,nstates,ISUPPZ,WORK,LWORK,IWORK,LIWORK, INFO )
  coupl=dble(0)
  call makedhmat(igeom,dhmat)
  do m=1,ncoord
    dhmat(m,:,:)=matmul(evec,matmul(dhmat(m,:,:),transpose(evec)))
    do k=1,3*natoms
      coupl(k,:,:)=coupl(k,:,:)+dhmat(m,:,:)*bmat(m,k)
    end do
  end do
END SUBROUTINE ccoupl

!------------------------------------------------------------
! CHESS calculate hessian matrix 
! cgeom      DOUBLE PRCISION(3*natoms),input  cartesian geometry
! surf       INTEGER,input   The state to evaluate
! hess       DOUBLE PRECISION(3*natoms,3*natoms),output  
!            hessian matrix in cartesian coordinates. 
SUBROUTINE chess(cgeom,surf,hess)
  use progdata ,only: natoms
  use hddata
  IMPLICIT NONE
  DOUBLE PRECISION,dimension(3*natoms),intent(in) :: cgeom
  INTEGER,intent(in)                              :: surf
  DOUBLE PRECISION,dimension(3*natoms,3*natoms),intent(out):: hess
  double precision,dimension(3*natoms)            :: geom,pgrad,mgrad
  double precision  ::  disp, s,d,e(nstates)
  integer           ::  m,n
  disp = 1D-6
  ! calculate energy gradients for displacements alone each direction
  do m=1,3*natoms
    geom=cgeom
    geom(m)=geom(m)+disp
    CALL cgrad(geom,surf,e,pgrad)
    geom=cgeom
    geom(m)=geom(m)-disp
    CALL cgrad(geom,surf,e,mgrad)
    hess(m,:)=(pgrad-mgrad)/(2*disp)
  end do
  ! symmetrize hessian matrix
  do m=2,3*natoms
    do n=1,m-1
      d=(hess(m,n)-hess(n,m))/2
      s=(hess(m,n)+hess(n,m))/2
      hess(m,n)=s
      hess(n,m)=s
    end do
  end do
END SUBROUTINE chess

!------------------------------------------------------------
! THESS calculate Hessian matrix in a transformed Cartesian space
! cgeom      (input) DOUBLE PRCISION,dimension(3*natoms)
!            cartesian geometry where Hessian matrix will be generated
! surf       (input) INTEGER
!            The state to be evaluated
! n          (input) INTEGER
!            Dimensionality of the transformed space
! tmat       (input)DOUBLE PRECISION,dimension(n,3*natoms)
!            Transformation matrix from Cartesian coordinates
! hess       (output)DOUBLE PRECISION(LDH,n)
!            Computed hessian matrix in transformed Cartesian space
! LDH        (input) INTEGER
!            Leading dimension of hess
SUBROUTINE thess(cgeom,surf,n,tmat,hess,LDH)
  use progdata ,only: natoms
  use hddata
  IMPLICIT NONE
  DOUBLE PRECISION,dimension(3*natoms),intent(in)  :: cgeom
  INTEGER,intent(in)                               :: surf,n,LDH
  DOUBLE PRECISION,dimension(n,3*natoms),intent(in):: tmat
  DOUBLE PRECISION,dimension(LDH,n),intent(out)    :: hess
  double precision,dimension(3*natoms)             :: pgradc,mgradc,geom
  double precision  ::  disp, s,d,e(nstates)
  integer           ::  i,j
  disp = 1D-5
  ! calculate energy gradients for displacements alone each direction
  do i=1,n
    geom=cgeom
    geom=geom+tmat(i,:)*disp
    CALL cgrad(geom,surf,e,pgradc)
    geom=cgeom
    geom=geom-tmat(i,:)*disp
    CALL cgrad(geom,surf,e,mgradc)
    CALL DGEMV('n',n,3*natoms,1/(2*disp),tmat,n,pgradc-mgradc,int(1),dble(0),hess(i,:),int(1))
  end do
  ! symmetrize hessian matrix
  do i=2,n
    do j=1,i-1
      d=(hess(j,i)-hess(i,j))/2
      s=(hess(j,i)+hess(i,j))/2
      if(abs(d)>1D-2*abs(s).and.abs(s)>1D-2)print 1000,d,i,j
      hess(j,i)=s
      hess(i,j)=s
    end do
  end do
1000 format(4X,"Symmetry breaking ",F10.4," in Hessian block ",I3,",",I3)
END SUBROUTINE thess
