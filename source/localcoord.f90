  !---------------------------------------------------------------------
  ! Construct a unitary transformation L that transforms cartesian
  ! coordinate into an orthongal system defined by the eigenvectors
  ! of matrix B^T.B, where B is Wilson's B matrix.
  ! Eigenvectors will be ordered with descending eigenvalues.  Coordinates
  ! with an eigenvalue larger than intGT will be considered internal.
  ! Number of internal degrees of freedom will be stored in field nvibs 
  !---------------------------------------------------------------------
  ! This subroutine is used to construct a local coordinate system that
  ! is used for the fitting.  All gradients and couplings will be fitted
  ! with this coordinate, which is made form linear combinations of the 
  ! 3*N cartesian coordinates, but has translation/rotation and dissappearing
  ! vibrations from dissociation separated from vibrational coordinates.
  ! The coordinates will also reflect the symmetry of the point, because 
  ! right singular vectors of the B matrix is used, which carries the full
  ! symmetry of the point as long as the internal coordinates are symmetric
  ! (which is always the case).
  ! Note: v2.7 the subroutine will now try to rotate the vectors which have
  ! degenerate eigenvalues so that they carry proper symmetry. 
  SUBROUTINE makeLocalIntCoord(ptdata,nstates,useIntG,intGT,intGS,nvibs)
    USE hddata, ONLY: ncoord
    USE progdata, ONLY: abpoint,natoms,printlvl
    IMPLICIT NONE
    double precision, parameter  :: degcap=1d-3 ! threshold for degeneracy
    TYPE(abpoint),intent(inout)  :: ptdata
    INTEGER,intent(IN)           :: nstates,nvibs
    LOGICAL,intent(IN)           :: useIntG
    DOUBLE PRECISION,intent(IN)  :: intGT,intGS
    integer                      :: i,n1,n2
    double precision,dimension(3*natoms,3*natoms) :: btb      !B^T.B
    double precision,dimension(3*natoms)          :: ev,grad  !eigenvalues of BtB
    double precision,dimension(45*natoms*natoms)  :: scr
    double precision             ::  nrmerr, alpha
    integer                      ::  INFO
    character(3)                 ::  str
   ! Number of degeneracies among coordinates
    integer     :: ndeg
   ! The list of coordinate ranges that are degenerate.  
    integer     :: degcoords(nvibs,2) 
   ! The starting index of the degenerate block. Coordinates from this index to
   ! index i-1 are all degenerate.
    integer     ::  degstart
    double precision,external :: dnrm2

    if(.not. useIntG)then
      ptdata%nvibs=nvibs
      ptdata%lmat=dble(0)
      do i=1,nvibs
        ptdata%lmat(i,i)=dble(1)
      end do
      ptdata%scale = dble(1)
      return
    end if

    ! Construct matrix B^T.B and get its eigenvalues
    call DSYRK('U','T',3*natoms,ncoord,1d0,ptdata%bmat,ncoord,0d0,btb,3*natoms)
    call DSYEV('V','U',3*natoms,btb,3*natoms,ev,scr,45*natoms*natoms,INFO)
    !reorder the eigenvectors in descending order
    do i=1,3*natoms/2
      CALL DSWAP(3*natoms,btb(:,i),int(1),btb(:,3*natoms+1-i),int(1))
    end do
     
    ndeg=0
    degstart=1
    ptdata%bmat=matmul(ptdata%bmat,btb)
    ptdata%lmat=btb
    do i=1,3*natoms
      ptdata%eval(i)=ev(3*natoms+1-i)
    end do
    if(printlvl>1)then
      print *,"  Eigenvalues of B^T.B"
      print "(10F10.4)",ptdata%eval
    end if
    ptdata%scale=dble(1)
    do i=1,3*natoms
      if(ev(3*natoms+1-i)>intGT)then
        !update the list of degenerate modes
        if(i>1)then
          if(abs(ptdata%eval(i)-ptdata%eval(i-1))>degcap)then
          !this mode is not degenerate with the previous  
             if(degstart<i-1)then !there are more than one degenerate coordiante
               ndeg=ndeg+1
               degcoords(ndeg,1)=degstart
               degcoords(ndeg,2)=i-1
             end if
             degstart=i
          end if
        end if
        ptdata%nvibs=i
        if(i>nvibs)then
        ! Degrees of freedom greater than 3N-5.  Something must have gone
        ! wrong if it ever gets here.
          print *,"PT",ptdata%id,"NVIBS(PT) = ",i,", natoms = ",natoms
          print *, "GEOM:"
          print "(3F11.7)",ptdata%CGEOM

          print *, "B Matrix "
          do n1=1,ncoord
            print "(12E11.2)",ptdata%bmat(n1,:)
          end do
          print *,"EV:",ev
          stop "Error: Number of degrees of freedom exceeded theoretical maximum"
        end if
        if(ev(3*natoms+1-i)<intGS)ptdata%scale(i)=ev(3*natoms+1-i)/intGS
      else
        ptdata%scale(i:)=dble(0)
        if(degstart<i-1)then
          ndeg=ndeg+1
          degcoords(ndeg,1)=degstart
          degcoords(ndeg,2)=i-1
        end if
        exit
      end if
    end do !i=1,3*natoms
    if(printlvl>1)then
      print "(A,I7)","     Number of internal degrees of freedom : ",ptdata%nvibs
      print "(A,I7)","     Number of degenerate internal motions : ",ndeg
    end if

    ! construct local coordinates
    do n1=1,nstates
      do n2=1,nstates
        do i=1,3*natoms
          grad(i)=dot_product(btb(:,i),ptdata%grads(:,n1,n2)) 
        end do
        ptdata%grads(:,n1,n2)=grad
        nrmerr=dnrm2(3*natoms-ptdata%nvibs,grad(ptdata%nvibs+1:),1)
        if(nrmerr*0d0 .eq. 0d0 .and. nrmerr.eq.nrmerr)then
            if(printlvl>1.and.nrmerr>1D-10.or.printlvl>2)print 1000,"residual norm of external gradients for block",&
                     n1,n2,dnrm2(3*natoms-ptdata%nvibs,grad(ptdata%nvibs+1:),1)
            if(nrmerr>1D-5)print "(4X,12F9.4)",grad(ptdata%nvibs+1:)
        end if
      end do
    end do ! n1=1,nstates

    ! Make linear combinations of degenerate internal coordinates to symmetrize
    ! them.
    do i=1,ndeg
      if(printlvl>2)print "(3X,2(A,I4))","Symmetrizing coordinate range ",degcoords(i,1)," to ",degcoords(i,2)
      
    end do

    ! Print out information about the final coordinate system
    if(ptdata%nvibs<3*natoms-6.or.printlvl>1)then
      write(str,'(I3)') 3*natoms
      print "(7X,A,"//trim(str)//"F7.3)","scaling factors of local coordinates: ",ptdata%scale
      if(ptdata%nvibs<3*natoms-6) then
        print *,"      Local coordinate system has reduced dimensionalities."
        if(printlvl>2)then
          print *,"       Cartesian vectos of vanishing coordinates:"
          do n1=ptdata%nvibs+1,3*natoms
             print "(A,I3,A,F9.5)","    Fitting coordinate ",n1," ,eigenvalue=",ptdata%eval(n1)
             print "(3F12.7)",btb(:,n1)
          end do
          print *,"       Overlaps of vanishing coordinates with nascent internals:"
          do n1=ptdata%nvibs+1,3*natoms
             print "(A,I3,A,F9.5)","    Fitting coordinate ",n1," ,eigenvalue=",ptdata%eval(n1)
             print "(12F9.4)",(dot_product(btb(:,n1),ptdata%bmat(n2,:)),n2=1,ncoord)
          end do
        end if
      end if
    end if
  1000 format(8x,A,I3,",",I3," : ",E11.4)
  END SUBROUTINE makeLocalIntCoord
