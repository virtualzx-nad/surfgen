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
  integer,parameter  ::  MaxGeoms = 10000,MaxAtoms=200
  character(255)     ::  geomfl 
  integer            ::  npts, i,j,k,natoms,nstates
  character(3),dimension(MaxAtoms)               :: atoms
  double precision,dimension(MaxAtoms)           :: anums,masses

  double precision,dimension(:),allocatable    :: e
  DOUBLE PRECISION,dimension(:,:),allocatable  :: cgeoms,h
  double precision,dimension(:,:,:),allocatable:: cg,dcg
  print *,"-------------------------------------------------"
  print *,"Entering testpoints, a surfgen point testing utility"
  print *,""
  print *,"  This program is part of the surfgen program"
  print *,"  2013 Yarkony group, Johns Hopkins University"
  print *,"-------------------------------------------------"

  geomfl = 'geom.all'
  print *,""
  print *,"Initializing potential"
  call initPotential()
  call getinfo(natoms,nstates)
  allocate(cgeoms(3*natoms,MaxGeoms))
  allocate(e(nstates))
  allocate(h(nstates,nstates))
  allocate(cg(3*natoms,nstates,nstates))
  allocate(dcg(3*natoms,nstates,nstates))

  print *," Number of states:",nstates
  print *," Number of natoms:",natoms
 
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
         print "(A,I3,A,I3,A)"," block (",j,",",k,")"
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
end program testpoints
