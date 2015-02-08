!testpoints
!------
!Point testing utility for surfgen.
!
! Written by Xiaolei Zhu, May 2013
! Modified by Christopher L Malbon August 2014
! Yarkony group,
! Johns Hopkins Univeristy
!
! For details see README file.
!
!This program uses surfgen evaluation library for construction and evaluation of

!Hd, as well as geometry input. 
program testpoints
  use hddata, only: nstates
  use progdata, only: abpoint
  implicit none
  integer,parameter  ::  MaxGeoms = 10000
  character(255)     ::  geomfl 
  integer            ::  npts, i,j,k,natoms,ios,ptid

  character(3),dimension(:),allocatable        :: atoms
  double precision,dimension(:),allocatable    :: anums,masses,e
  double precision,dimension(:,:),allocatable  :: cgeoms,h
  double precision,dimension(:,:,:),allocatable:: cg,dcg
  double precision,external   :: dnrm2
  
  real*8,dimension(:),allocatable     :: ab_energy
  real*8,dimension(:,:),allocatable :: ab_grads
  character(len=200)                  :: inputfile
  integer                             :: npoints
  logical,dimension(:,:),allocatable  :: hasgrad
  real*8                              :: degTol
  character(len=200)                  :: perfRot
  integer                             :: st1read,st2read
  ! This is the data structure necessary to use the subroutine OrthGH_ab
  ! It is located in progdata.
  type(abpoint)   :: pt
  ! Only testing for 1 point
  npoints=1
  pt%id=1

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
  allocate(ab_energy(nstates))
  allocate(ab_grads(3*natoms,nstates*(nstates+1)/2))
! Allocate arrays for pt
  allocate(pt%energy(nstates,nstates))
  allocate(pt%grads(3*natoms,nstates,nstates))
  allocate(hasgrad(nstates,nstates))

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
    ! Call 
  end do!i

  call get_command_argument(number=2,value=perfRot,status=ios)
  if(ios.ne.0) then ! Do not perform GH rotation
        stop
  end if
  ! For OrthGH_ab, we need to generate all of the information required for
  ! type(pt), or, at least, the inforamation OrthGH_ab needs. It needs:
  !   pt%...
  !    ndeggrp
  !    deg_groups
  !    grads
  !    energy
  !    nvibs

  ! Read in ab initio energies
  ! The file should be called "energy.file"
  call readEner("energy.file",npoints,nstates,ab_energy,st1read,st2read)
  if(npoints .eq. 0)then
    print "(x,A)", " No energy data found."
    stop "Terminating."
  end if
  
  ! Read in ab initio gradients and couplings
  ! The gradient file should be called "abgrad.s$i.file"
  ! The coupling files should be called "abnadc.s$i.s$j.file"
  do i=1,nstates
      do j=i,nstates
            inputfile=filename(i,j,"abgrad.s$.file","abnadc.s$.s$.file")
            call readGrads(inputfile,npoints,natoms,ab_grads(1,(i*(i-1)/2 + j)))
            if(npoints .eq. 0)then
                  print "(x,A)", " No gradient data found."
                  cycle
            end if
      end do
  end do
  
  ! Move information for readGrads and readEner into pt, which is the type
  ! required for orthogonalization of GH
  do i=1,nstates
      pt%energy(i,i)=ab_energy(i)
      do j=i,nstates
            pt%grads(:,i,j)=ab_grads(:,(i*(i-1)/2+j))
            hasgrad(i,j)=.true.           !We have this gradient
            if (i.ne.j)then
                  pt%grads(:,j,i)=pt%grads(:,i,j)
                  hasgrad(j,i)=.true.     !We have this gradient
            end if
      end do
  end do
  
  ! We now need to generate pt%nvibs
  ! This is accomplished by calling makeLocalIntCoord. 
  !call makeLocalIntCoord(pt,nstates,.true.,1D-3,1D-1,3*natoms-6,0)
  pt%nvibs=natoms*3     ! Cartesians
  ! We now need to generate pt%ndeggrp and pt%deg_groups
  ! The subroutine genDegs will accopmlish this for us
  degTol=1D-1
  pt%lb=st1read
  pt%ub=st2read
  call genEnerGroups(pt,degTol)

  ! Now call the orthogonalization subroutine, OrthGH_ab
  ! OrthGH_ab(pt,maxiter,toler,hasGrad)
  call OrthGH_ab(pt,100,1d-8,hasgrad)

  ! get neighboring point index
  call getNeighbor(ptid)
  if(ptid/=0) print *,"Index of Closest Data Point : ", ptid         
contains
      FUNCTION filename(s1,s2,grdptn,cpptn)
            IMPLICIT NONE
            INTEGER,INTENT(IN)          :: s1,s2
            CHARACTER(LEN=*),INTENT(IN) :: grdptn,cpptn
            CHARACTER(255)              :: filename
            CHARACTER(1)                :: st1,st2
            integer :: i

            write(st1,'(i1)')s1
            write(st2,'(i1)')s2

            if(s1.eq.s2)then
                  filename=grdptn
                  i=index(filename,'$')
                  if(i>0)filename(i:i) = st1
                  if(i==0)filename=''
            else
                  filename=cpptn
                  i=index(filename,'$')
                  if(i>0)filename(i:i) = st1
                  i=index(filename,'$')
                  if(i>0)filename(i:i) = st2
                  if(i==0)filename=''
            endif

            filename = trim(adjustl(filename))
      END FUNCTION filename
end program testpoints
