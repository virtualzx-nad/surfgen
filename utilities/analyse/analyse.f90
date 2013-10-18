! analyse 
!-------------------
! read columbus geometries and show internal coordinates 
!-------------------
! This is part of the surfgen program package.  To compile this program,
! first generate surfgen library with `make libs` in surfgen main directory,
! then run the `install.sh` script in the program's directory.
!-------------------
! (C) 2013 Yarkony Group, Johns Hopkins University
!-------------------
! Oct 2013        Xiaolei Zhu        Created 
!-------------------
program analyse 
  implicit none
  integer::   ios
  integer,parameter :: MAXGEOMS=50000
  character*300 geomfile,str 
  integer :: natm , nst, ngeoms
  integer :: i
  double precision,allocatable  :: anum(:),masses(:),  hess(:,:), w(:)
  character*3,allocatable       :: aname(:)
  double precision,allocatable  :: cgeom(:,:)

  print *," ***************************************** "
  print *," *    analyse.x                          * "
  print *," ***************************************** "
  print *," Minimum energy crossing search on Surfgen potential surfaces"
  call initPotential()
  call getInfo(natm,nst)
  call DisableEnergyScaling()

! allocate arrays
  allocate(masses(natm))
  allocate(anum(natm))
  allocate(aname(natm))
  allocate(hess(natm*3,natm*3))
  allocate(w(3*natm))
  allocate(cgeom(3*natm,MAXGEOMS))

  print "(A,I6)","Number of Atoms:  ",natm
  print "(A,I6)","Number of States: ",nst

! process arguments
! synposis:    analyse.x geomfile 
! Default values:
! geomfile        geom
  call get_command_argument(number=1,value=geomfile,status=ios)
  if(ios<-1)  &     ! input argument larger than container
      stop "Filename truncated.  Please use less than 300 characters"
  if(ios>0) then
    print *,"No filename suppied.  Using default."
    write(geomfile,"(A)"), "geom"
  end if

  print *,"Reading input from input file "//trim(geomfile)
  ngeoms=MAXGEOMS
  call readColGeom(geomfile,ngeoms,natm,aname,anum,cgeom,masses)
 
  print "(A,I10,A)","Analysing ",ngeoms," geometries"
  ! print initial geometry information
  do i=1,ngeoms
    print "(A,I10,A)","== Geometry No. ",i," ==="
    call analysegeom(natm,cgeom(1,i),aname,anum,masses,2.2d0,.True.)
  end do

! deallocate arrays
  deallocate(masses)
  deallocate(anum)
  deallocate(aname)
  deallocate(hess)
  deallocate(w)
  deallocate(cgeom)
end program

