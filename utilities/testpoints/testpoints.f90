!testpoints
!------
!Point testing utility for surfgen.
!
! Written by Xiaolei Zhu, May 2013
! Yarkony group, Johns Hopkins Univeristy
!
! For details see README file.
!
!This program uses surfgen evaluation library for construction and evaluation of
!Hd, as well as geometry input. 
program testpoints
  implicit none
  integer,parameter  ::  MaxGeoms = 10000
  character(255)     ::  geomfl 
  integer            ::  npts, i,j,k,natoms,nstates,ios,ptid

  character(3),dimension(:),allocatable        :: atoms
  double precision,dimension(:),allocatable    :: anums,masses,e
  double precision,dimension(:,:),allocatable  :: cgeoms,h
  double precision,dimension(:,:,:),allocatable:: cg,dcg
  double precision,external :: dnrm2
  print *,"-------------------------------------------------"
  print *,"Entering testpoints, a surfgen point testing utility"
  print *,""
  print *,"  This program is part of the surfgen program"
  print *,"  2013 Yarkony group, Johns Hopkins University"
  print *,"-------------------------------------------------"

  print *,""
  print *,"Initializing potential"
  call initPotential()
  call getinfo(natoms,nstates)

! allocate arrays
  allocate(atoms(natoms))
  allocate(anums(natoms))
  allocate(masses(natoms))
  allocate(cgeoms(3*natoms,MaxGeoms))
  allocate(e(nstates))
  allocate(h(nstates,nstates))
  allocate(cg(3*natoms,nstates,nstates))
  allocate(dcg(3*natoms,nstates,nstates))

  print *," Number of states:",nstates
  print *," Number of natoms:",natoms
 
  call get_command_argument(number=1,value=geomfl,status=ios)
  if(ios/=0)then
    print *,"Cannot get file name from command line.  Using default."
    write(geomfl,"(A)")"geom.all"
  end if
  print "(A,I8)","Processing geometry input. Maximum number of geometries allowed:",MaxGeoms
  print "(A)","  Filename:"//trim(geomfl)
  npts = MaxGeoms
  call readColGeom(geomfl,npts,natoms,atoms,anums,cgeoms,masses)

  print "(A,I8)","  Number of geometries found in file: ",npts

  do i=1,npts
    print *,""
    print *,"Cartesian geometry"
    print "(3F18.10)",cgeoms(:,i)
    print "(2x,A,I8,A)","Hd predictions for point #",i," :"
    call EvaluateSurfgen(cgeoms(1,i),e,cg,h,dcg)
    print *,""
    print *,"Quasi-diabatic Hamiltonian"
    do j=1,nstates
      print "(2x,10F24.15)",h(j,:)
    end do
    print *,""
    print *,"Adiabatic energy(a.u.)"
    print "(2x,10F24.15)",e
    print *,"Adiabatic energy(cm-1)"
    print "(2x,10F24.15)",e*219474.6305d0
    print *,""
    print *,"Cartesian Gradients and Couplings in Adiabatic Representation"
    do j=1,nstates
      do k=1,j
         print "(A,I3,A,I3,A,E13.5)"," block (",j,",",k,"),|F|=",dnrm2(3*natoms,cg(1,j,k),1)
         print "(3F18.10)",cg(:,j,k)
      end do!k
    end do!j
    print *,""
    print *,"Cartesian Gradients in Diabatic Representation"
    do j=1,nstates
      do k=1,j
         print "(A,I3,A,I3,A)"," block (",j,",",k,")"
         print "(3F18.10)",dcg(:,j,k)
      end do!k
    end do!j
  end do!i

  ! get neighboring point index
  call getNeighbor(ptid)
  if(ptid/=0) print *,"Index of Closest Data Point : ", ptid
end program testpoints
