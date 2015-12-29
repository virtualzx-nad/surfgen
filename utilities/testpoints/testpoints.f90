!testpoints
!------
!Point testing utility for surfgen.
!
! Written by Xiaolei Zhu, May 2013
! Yarkony group, Johns Hopkins Univeristy
!
! For details see README file.
!
!-- Change Log --
! 12-05-2015 - CLM - Added input file routine.
!                    Added printing of s vector.
!                    Added molden output of g/h/s vectors.
! 12-18-2015 - CLM - Molden output is mass-weighted if selected
!                     in input.
!                    Added print level implementation.
!----------------
!This program uses surfgen evaluation library for construction and evaluation of
!Hd, as well as geometry input. 
program testpoints
  implicit none
  integer,parameter  ::  MaxGeoms = 10000
  real*8,  parameter ::  au2ang = 0.529177249
  character(255)     ::  geomfl 
  integer            ::  npts, i,j,k,l,m,natoms,nstates,ios,ptid

  character(3),dimension(:),allocatable        :: atoms
  double precision,dimension(:),allocatable    :: anums,masses,e,dvec,norms
  double precision,dimension(:,:),allocatable  :: cgeoms,h
  double precision,dimension(:,:,:),allocatable:: cg,dcg
  double precision,external :: dnrm2

  logical :: printm, mweight ! print molden output, print molden as 
                             ! mass-weighted vectors
  integer :: printl          ! print level
  real*8,  dimension(:,:,:),allocatable :: hvec, svec, gvec !h/s/g vectors
  integer, dimension(:),    allocatable :: mex

  print *,"-------------------------------------------------"
  print *,"Entering testpoints, a surfgen point testing utility"
  print *,""
  print *,"  This program is part of the surfgen program"
  print *,"  2013 Yarkony group, Johns Hopkins University"
  print *,"-------------------------------------------------"

  print *,""
  print *,"Checking for input file"
  allocate(mex(2))
  call read_input_file(printm, mweight, printl, mex)
  print *,""
  print *,"Initializing potential"
  call initPotential()
  call getinfo(natoms,nstates)

! allocate arrays
  allocate(atoms(natoms))
  allocate(anums(natoms))
  allocate(masses(natoms))
  allocate(cgeoms(3*natoms,MaxGeoms))
  allocate(dvec(3*natoms))
  allocate(e(nstates))
  allocate(norms(nstates))
  allocate(h(nstates,nstates))
  allocate(cg(3*natoms,nstates,nstates))
  allocate(dcg(3*natoms,nstates,nstates))
  allocate(hvec(3*natoms,nstates,nstates))
  allocate(svec(3*natoms,nstates,nstates))
  allocate(gvec(3*natoms,nstates,nstates))

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
    print *,"Diabatic energy(cm-1)"
    print "(2x,10F24.15)",(h(j,j)*219474.6305d0,j=1,nstates)
    print *,""
    print *,"Adiabatic energy(a.u.)"
    print "(2x,10F24.15)",e
    print *,"Adiabatic energy(cm-1)"
    print "(2x,10F24.15)",e*219474.6305d0
    print *,""
    print *,"Norms of g vectors"
    print "(8x,30I17)",(j,j=1,nstates-1)
    do j=2,nstates
      do k=1,j-1
        dvec=cg(:,j,j)-cg(:,k,k)
        gvec(:,j,k)=dvec/2
        norms(k)=dnrm2(3*natoms,dvec,1)/2
        if (mweight) then
          do l=1,natoms
            do m=1,3
              gvec((l-1)*3+m,j,k)=gvec((l-1)*3+m,j,k)/dsqrt(masses(l))
            end do
          end do
        end if
      end do
      print "(I6,30E17.7)",j,norms(1:j-1)
    end do
    print *,""
    print *,"Norms of h vectors"
    print "(8x,30I17)",(j,j=1,nstates-1)
    do j=2,nstates
      do k=1,j-1
        dvec=cg(:,j,k)*(e(j)-e(k))
        hvec(:,j,k)=dvec
        norms(k)=dnrm2(3*natoms,dvec,1)
        if (mweight) then
          do l=1,natoms
            do m=1,3
              hvec((l-1)*3+m,j,k) = hvec((l-1)*3+m,j,k)/dsqrt(masses(l))
              end do
          end do
        end if
      end do
      print "(I6,30E17.7)",j,norms(1:j-1)
    end do
    print *,""
    print *,"Norms of s vectors"
    print "(8x,30I17)",(j,j=1,nstates-1)
    do j=2,nstates
      do k=1,j-1
        dvec=cg(:,j,j)+cg(:,k,k)
        svec(:,j,k)=dvec/2
        norms(k)=dnrm2(3*natoms,dvec,1)/2
        if (mweight) then
          do l=1,natoms
            do m=1,3
              svec((l-1)*3+m,j,k) = svec((l-1)*3+m,j,k)/dsqrt(masses(l))
            end do
          end do
        end if
      end do
      print "(I6,30E17.7)",j,norms(1:j-1)
    end do
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
    if (printm .and. i .eq. npts .and. mex(1) .ne. 0) then
            if (mweight) then
                print *, "ALERT: mass-weighted g, h and s being used!"
            end if
            call print_molden_output(gvec, hvec, svec, nstates, &
                natoms, atoms, cgeoms(:,i), mex, printl)    
            call rot_g_h_vectors(gvec(:,mex(2),mex(1)),hvec(:,mex(2),mex(1)),&
                natoms) 
            call express_s_in_ghplane(gvec(:,mex(2),mex(1)), &
                hvec(:,mex(2),mex(1)), svec(:,mex(2),mex(1)), natoms)
    end if
  end do!i

  ! get neighboring point index
  call getNeighbor(ptid)
  if(ptid/=0) print *,"Index of Closest Data Point : ", ptid

contains
  ! rot_g_h_vectors: rotate g and h vectors
  subroutine rot_g_h_vectors(g, h, na)
          implicit none
          integer,intent(in) :: na
          real*8, dimension(3*na), intent(inout) :: g, h
          real*8, dimension(3*na) :: tmp1, tmp2
          real*8 :: gh, hh, gg, atin, theta
          real*8, external :: ddot

          gh = ddot(na*3,g,1,h,1)
          if (gh .lt. 1d-8) return
          hh = ddot(na*3,h,1,h,1)
          gg = ddot(na*3,g,1,g,1)
          atin = (2 * gh) / (hh - gg)
          theta= datan(atin) / 4.0

          tmp1 = g
          tmp2 = h
          print "(A,f13.5)"," g.h = ", gh
          print "(A,f8.3)"," Rotating by ", theta
          g = dcos(2 * theta) * tmp1 - dsin(2 * theta) * tmp2
          h = dsin(2 * theta) * tmp1 + dcos(2 * theta) * tmp2
          gh = ddot(na*3,g,1,h,1)
          print "(A,f13.5)"," g.h = ", gh
          
          return
  end subroutine rot_g_h_vectors

  ! express_s_in_ghplane: gives components of s in g/h plane.
  subroutine express_s_in_ghplane(g, h, s, na)
          implicit none
          integer, intent(in) :: na
          real*8, dimension(3*na), intent(in) :: g, h, s
          real*8, dimension(3*na) :: gx, hy
          real*8, external :: ddot
          real*8, dimension(2) :: sxy
          
          gx = g / dsqrt(ddot(3*na,g,1,g,1))
          hy = h / dsqrt(ddot(3*na,h,1,h,1))

          sxy(1) = ddot(3*na,s,1,gx,1)
          sxy(2) = ddot(3*na,s,1,hy,1)
          
          print "('||s|| = ',f15.8)", dnrm2(3*na,s,1)
          print "('S vector in g-h plane:')"
          print "('S(x,y) = (',f15.8,',',f15.8,')')", sxy(1), sxy(2)
        
          return
  end subroutine express_s_in_ghplane

  ! read_input_file: reads input file for testpoints.x
  subroutine read_input_file(printm, mweight, printl, mex)
          implicit none
          logical, intent(out) :: printm, mweight
          integer, intent(out) :: printl
          integer, dimension(2), intent(inout) :: mex

          character(25) :: flname = 'testpoints.in'
          integer       :: flunit = 21, ios
          integer, parameter :: fl_not_found = 29

          namelist /testoutput/ printm, mweight, printl, mex
          printm = .false.
          mweight= .false.
          printl = 0
          mex(1) = 0
          mex(2) = 0

          ! Open input file. If file is not found, set printm = .false.
          ! and return. Else, read the file.
          open(file = flname, unit = flunit, status = 'old', &
                  action = 'read', iostat = ios)
          if (ios .eq. 0) then
                  read(unit = flunit, nml = testoutput)
                  return
          else
                  print "(a,i5,a)", "  **Error ", ios, &
                          " occurred opening input file!**"
                  return
          end if
          return
  end subroutine read_input_file

  ! print_molden_output: print g/h/s vectors in molden format
  subroutine print_molden_output(gv, hv, sv, ns, na, a, g, mex, pl)
          implicit none
          integer, intent(in) :: ns, na, pl
          integer, dimension(2),          intent(in) :: mex
          character(3), dimension(na),    intent(in) :: a
          real*8,  dimension(3*na,ns,ns), intent(in) :: gv, hv, sv
          real*8,  dimension(3*na),       intent(in) :: g
          real*8,  dimension(:),         allocatable :: fq
          character(25) :: flnm = 'vectors.molden'
          integer       :: flun = 22, ios, i, j, nfq
          
          ! Allocate dummy frequencies. There will be n(n+1)/2, were 
          ! n = states. Open file. Print vectors in order: g, h, s. 
          !Lower triangle is printed.
          
          nfq = 3 * ns * (ns + 1) / 2
          allocate(fq(nfq))
          fq = 1d0

          open(file = flnm, unit = flun, status = 'unknown', action = 'write',&
                  iostat = ios)
          if (ios .ne. 0) then
                  print *, "**Error ", ios, " occurred opening vectors.molden!"
                  return
          end if        
          write(flun,"(1x,'[Molden Format]')")

          call print_molden_freq(flun, fq, 3)
          call print_molden_geom(flun, g, a, na)
          call print_molden_vecs(flun, gv, hv, sv, ns, na, nfq, mex, pl)

          close(flun)
          return
  end subroutine print_molden_output

  ! print_molden_vecs: print vectors in molden format. Expects g/h/s
  subroutine print_molden_vecs(u, gv, hv, sv, ns, na, nfq, mex, pl)
          implicit none
          integer, intent(in) :: nfq, na, ns, u, pl
          integer, dimension(2), intent(in) :: mex
          real*8,  dimension(3*na,ns,ns), intent(in) :: gv, hv, sv
          real*8 :: ngv, nhv, nsv
          integer :: i, j, p
          
          j = mex(1)
          i = mex(2)
          ngv = dnrm2(3*na,gv(:,i,j),1)
          nhv = dnrm2(3*na,hv(:,i,j),1)
          nsv = dnrm2(3*na,sv(:,i,j),1)

          write(u,"(1x,'[FR-NORM-COORD]')")
          p = 1
          write(u,"(1x,'vibration',i24)") p
          if (pl > 0) then
                  print "(a,2i2,a)", &
                          " g-vector for ", mex(1), mex(2)," intersection:"
                  print "(3f12.5)", gv(:,i,j)
          end if
          write(u,"(3f12.5)") gv(:,i,j)/ngv
          p = 2
          write(u,"(1x,'vibration',i24)") p
          if (pl > 0) then
                  print "(a,2i2,a)", &
                          " h-vector for ", mex(1), mex(2)," intersection:"
                  print "(3f12.5)", hv(:,i,j)
          end if
          write(u,"(3f12.5)") hv(:,i,j)/nhv
          p = 3
          write(u,"(1x,'vibration',i24)") p
          write(u,"(3f12.5)") sv(:,i,j)/nsv
          return
  end subroutine print_molden_vecs


  ! print_molden_geom: print geometry in molden format, file is open.
  subroutine print_molden_geom(u, g, a, na)
          implicit none
          integer, intent(in) :: u, na
          character(3), dimension(na),intent(in) :: a
          real*8,  dimension(3 * na), intent(in) :: g
          integer :: i

          write(u,"(1x,'[FR-COORD]')")
          do i = 1, na
                write(u,"(1x,a3,3f13.6)") a(i), g((i-1)*3+1), g((i-1)*3+2),&
                        g((i-1)*3+3)
          end do
          return
  end subroutine print_molden_geom

  ! print_molden_freq: print frequencies to molden file, file is opened.
  subroutine print_molden_freq(u, f, nf)
          implicit none
          integer, intent(in) :: u, nf
          real*8,  dimension(nf), intent(in) :: f
          integer :: i

          write(u,"(1x,'[FREQ]')")
          do i = 1, nf
                write(u,"(f10.2)") f(i)
          end do

          return
  end subroutine print_molden_freq
end program testpoints
