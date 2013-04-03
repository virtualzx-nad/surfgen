
!******************************************************************************
! MAIN PROGRAM
program surfgen
    use progdata
    implicit none
    integer   ::jobtype

    print *,"Entering surfgen.global"
    call readinput(jobtype)
    call initialize(jobtype)

    select case(jobtype)
     case(1)
      if(printlvl>0)print *,"Reading displacement data..."
      call readdisps()
      if(printlvl>0)print *,"Making surface..."
      call makesurf()
    end select

  call cleanup()

end program surfgen
! read general input from input file (surfgen.in)
! and read in job specific input according to jtype
SUBROUTINE readinput(jtype)
  use progdata
  use hddata, only: initGrps,nstates,order
  use CNPI, only: irrep,GrpPrty,GrpSym,nSymLineUps,nirrep
  IMPLICIT NONE
  INTEGER,INTENT(INOUT)           :: jtype
  INTEGER                         :: i,j,ios
  INTEGER                         :: nGrp,jobtype
  INTEGER,dimension(10)           :: surface,updatehess
  INTEGER,dimension(50)           :: ci,atmgrp
  DOUBLE PRECISION,dimension(500) :: minguess,mexguess
  INTEGER,DIMENSION(20,MAX_ALLOWED_SYM)        :: groupsym,groupPrty

  NAMELIST /GENERAL/        jobtype,natoms,order,nGrp,groupsym,groupprty,&
                            switchdiab, printlvl,inputfl,atmgrp,cntfl,nSymLineUps

  nSymLineUps = 1
  jtype      = 0
  natoms     = 0
  order      = 2
  printlvl   = 1
  cntfl      = 'connect.in'
  inputfl    = 'hd.data'
  switchdiab = .false.
  print *,"Entering readinput()."
 !----------- GENERAL SECTION ----------------!
  open(unit=INPUTFILE,file='surfgen.in',access='sequential',form='formatted',&
       IOSTAT=ios,POSITION='REWIND',ACTION='READ',STATUS='OLD')
  if(ios/=0)then
    print *,"readinput: cannot open file surfgen.in.  IOSTAT=",ios
  end if!ios/=0

  read(unit=INPUTFILE,NML=GENERAL)
  if(printlvl>0)print *,"    readinput():  Control parameters read from surfgen.in"
  jtype  = jobtype

  call genAtomList(atmgrp)

  if(jtype>=0)then  ! jobtype<0 will only print out symmetry operations
      call readIrreps()
      if(allocated(GrpSym))deallocate(GrpSym)
      if(allocated(GrpPrty))deallocate(GrpPrty)
      allocate(GrpSym(nGrp,nSymLineUps))
      allocate(GrpPrty(nGrp,nSymLineUps))
      GrpSym=groupsym(1:nGrp,1:nSymLineUps)
      GrpPrty=groupprty(1:nGrp,1:nSymLineUps) 
      ! get the number of allowed symmetries for each group
      if(nSymLineUps<1)stop"ERROR: There has to be at least one set of symmetry setups"
      do i=1,nGrp
        do j=1,nSymLineUps
          if(GrpSym(i,j)<1.or.GrpSym(i,j)>nirrep.or.irrep(GrpSym(i,j))%Dim.ne.irrep(GrpSym(i,1))%Dim &
                .or. abs(GrpPrty(i,j)).ne.1 ) stop "Error: Incorrect symmetry input."
        end do!j
      end do
      call initGrps(nGrp,irrep(GrpSym(:,1))%Dim)
  end if!(jtype>=0)

  if(printlvl>0)print *,"    Reading job specific inputs."
  !------------ JOBTYPE SECTION ------------!

  select case(jtype)

   !-----------------------------------------------------
   ! Printing symmetry operations ( permutations ) only
   !----------------------------------------------------
   case(-1)
    if(printlvl>0)print *,"    readinput():  jobtype<0.  Printig permutations."

   !-----------------------------------------------------
   ! Potlib interface.  Do nothing.
   !----------------------------------------------------
   case(0)
    if(printlvl>0)print *,"    readinput():  jobtype=0.  Ready for evaluation."


   !-----------------------------------------------------
   ! Read in displacemnts, construct surface
   !----------------------------------------------------
   case(1)
    if(printlvl>0)print *,"    readinput():  jobtype=1.  Reading displacements"
    CALL readMAKESURF(INPUTFILE)
  end select

  close(unit=INPUTFILE)


  if(printlvl>0)print *,"Exiting readinput..."
  return
end SUBROUTINE readinput
