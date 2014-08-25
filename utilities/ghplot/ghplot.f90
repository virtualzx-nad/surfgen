!-------------------------------------------
! ghplot   :   Utility program for surfgen
!-------------------------------------------
! DESCRIPTION
!
! This program calculates the energies on a two dimensional subspace and
! tabulate the numbers, along with the coordinates, in comma separated
! value(csv) files.   
! There is also option to enable output of the diabatic contributions. 
! Each electronic state will enjoy a separate output file.  This is 
! usually used to make plots of the branching plane.
!
! USAGE
!
! You have to change the constants at the beginning of the file then
! recompile for your system.
!
! PLOTTING
!
! Currently the program takes an initial geometry and then do cartesian
! displacements along two vectors.  The number of atoms, electronic states, as
! well as the initial geometries and two displacement vectors needs to be
! specified.   Below is an example for plotting the branching plane of S1-S2 MEX
! in phenol, where there are 4 electronic states.  
!
! VARIABLES THAT NEED CHANGES
!
! The following variables need to be changed in the program before you can plot
! for your own systems.
!
! natm         INTEGER
!              Number of atoms

! nst          INTEGER
!              Number of electronic states

! cgeom0       DOUBLE PRECISION, dimension(3*natm)
!              Initial geometry of the plot.

! vx, vy       DOUBLE PREICISION, dimension(3*natm)
!              The two displacement vectors.  Both will be normalized.

! dx, dy       DOUBLE PRECISION
!              Size of displacement to take between grid points.

! xstart, xend, ystart, yend  
!              INTEGER
!              Index of initial and final displacement for the two directions. 

! printevec    LOGICAL 
!              Whether or not to print out also the eigenvectors at each point.
program ghplot
implicit none

! Information section  >>>>>>>>>>>>>>>>>>>>>>>>>>>|  change these 
  integer,parameter :: natm = 5, nst=3          !||
! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<|

! definition section  ======shouldn't need to change here
  double precision,parameter :: au2cm1=219474.6305d0
  double precision, parameter  ::  bohr2ang=0.529177249d0
  double precision,dimension(3*natm)  ::  cgeom,vx, vy, cgeom0,cgeom1
  integer  :: i,j, isurf, fid
  double precision  ::  dx, dy, xstart, xend, ystart, yend
  real*8    :: h(nst,nst),cg(3*natm,nst,nst),dcg(3*natm,nst,nst),e(nst),evec(nst,nst)
  double precision,external :: dnrm2
  character*9,dimension(nst)  ::  fnames
  character*50   :: fstr
  logical :: printevec

! Initializations
  call initPotential()

! >>>>>>> Initialize the necessary variables. >>>>>>>>>>>>>>>>>>>   CHANGE THESE
! cgeom0 is the initial geometry of the origin
! the S1-S2 MEX of phenol
  cgeom0    =(/ 1.05357173,    0.41993810,   -0.26545773, &
            -1.31116182,    0.02919745,   -0.00291586, &
            -1.97087721,   -1.87000337,   -0.57933605, &
            -3.02278578,    1.88979365,    0.96496662, &
             2.02264645,   -0.99777181,   -0.94564838 /)

! vx and vy are displacement vectors  
  vx = 0d0
  vx(10:12)=cgeom0(10:12)-cgeom0(4:6)
  print *,"r0=",dnrm2(3*natm,vx,1)*bohr2ang
  call EvaluateSurfgen(cgeom0,e,cg,h,dcg)
  vy = dcg(:,2,3)

! dx and dy are size of the displacement between grid points 
  dx = 6d-2
  dy = 6d-2

! starting and ending indicies for the plotting
  xstart=-6
  xend  = 6
  ystart= -6
  yend  = 6

! printevec dictates whether the eigenvectors will be recorded.
  printevec =.true.

! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
! Plotting program.  Most likely anything below here does not require changing. 

! setup filenames
  do i=1,nst
    write(fnames(i),"(A,I1,A)")"surf",i,".csv"
  end do

! normalizing the two vectors
  vx = vx/dnrm2(3*natm,vx,1)
  vy = vy/dnrm2(3*natm,vy,1)

! construct format string

 if(printevec)then
   write (fstr,"(A,I1,A)"),"(2(F12.7,', '),F20.5,",nst,"(',',F12.8))"
 else
   write (fstr,"(A)"),"(2(F12.7,', '),F20.5)"
 end if

! Plot branching space
  print *,"Plotting surfaces: "
  do isurf=1,nst
   print "(A,I2,2A)"," Generating data for state ",isurf,": ",fnames(isurf)
   fid = 123+isurf
   open(unit=fid,file=fnames(isurf),access='sequential',action='write')
   do i=xstart,xend 
      cgeom1 = cgeom0+i*dx*vx
      do j=ystart,yend
        cgeom = cgeom1+j*dy*vy
        call EvaluateSurfgen(cgeom,e,cg,h,dcg)
        if(printevec)then
          call getEvec(evec)
          write(unit=fid,fmt=fstr),i*dx*bohr2ang,j*dy*bohr2ang,e(isurf)*au2cm1,evec(:,isurf)
        else
          write(unit=fid,fmt=fstr),i*dx*bohr2ang,j*dy*bohr2ang,e(isurf)*au2cm1
        end if
      end do
   end do!i
   close(fid)
  end do
 end 
