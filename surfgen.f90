
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
     case(2)
      if(printlvl>0)print *,"Performing critical point searches..."
      call minmex()
     case(3)
      if(printlvl>0)print *,"Loading geometries..."
      call loadgeom()
    end select

  call cleanup()

end program surfgen
! read general input from input file (surfgen.in)
! and read in job specific input according to jtype
SUBROUTINE readinput(jtype)
  use progdata
  use hddata, only: initGrps,ncoord,nstates,order
  use CNPI, only: irrep,grpPrty,grpSym
  IMPLICIT NONE
  INTEGER,INTENT(INOUT)           :: jtype
  INTEGER                         :: i,j,mkadiabat,ios
  INTEGER                         :: nGrp,jobtype
  INTEGER,dimension(10)           :: surface,updatehess
  INTEGER,dimension(50)           :: ci,atmgrp
  DOUBLE PRECISION,dimension(500) :: minguess,mexguess
  INTEGER,DIMENSION(10)           :: groupSym,groupPrty
  DOUBLE PRECISION,dimension(20)  :: e_guess

  NAMELIST /GENERAL/        jobtype,natoms,order,nGrp,groupsym,groupprty,usefij,switchdiab,&
                            printlvl,deg_cap,inputfl,eshift,atmgrp,use_eguess,e_guess
  NAMELIST /MINMEX/         nmin,nmex,minguess,mexguess,optiter,opttoler,updatehess,h_recal,surface,&
                             ci,maxstep,xscale,sscale,degtoler,enforcepd
  NAMELIST /LOADGEOM/       isloop,loopst,ngeoms,geomfl,calchess,outputdir,outputdiab
 
  jtype      = 0
  natoms     = 2
  printlvl   = 1
  inputfl    = ''
  usefij     = .true.
  deg_cap    = 1D-7
  eshift     = dble(0)
  use_eguess = .false.
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
  
  call readIrreps()
  if(allocated(GrpSym))deallocate(GrpSym)
  if(allocated(grpPrty))deallocate(grpPrty)
  allocate(GrpSym(nGrp))
  allocate(grpPrty(nGrp))
  GrpSym=groupsym(1:nGrp)
  grpPrty=groupprty(1:nGrp)
  call initGrps(nGrp,irrep(GrpSym(:))%Dim)

  if(allocated(eguess))deallocate(eguess)
  allocate(eguess(nGrp))
  if(use_eguess)then
    eguess=e_guess(1:nGrp)/AU2CM1
    print *,"Guess energies: ", eguess
  else
    do i=1,nGrp
     eguess(i)=i*10000/AU2CM1
    end do
  end if

  if(printlvl>0)print *,"    Reading job specific inputs."
  !------------ JOBTYPE SECTION ------------!

  select case(jtype)

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
   !--------------------------------------------------
   ! Find extrema on the surface
   !--------------------------------------------------
   case(2)
    maxstep=0.2
    minguess=dble(0)
    mexguess=dble(0)
    updatehess=0
    surface=0
    ci=0
    minstates=0
    mexstates=0
    nmin = 0
    nmex = 0
    optiter = 3
    opttoler = 1D-3
    degtoler = 1D-2
    xscale   = 1D0
    sscale   = 1D0
    h_update(1) = .false.
    h_recal = 1
    enforcepd = .false.

    read(unit=INPUTFILE,NML=MINMEX)
    if(nmin.gt.0)allocate(minstart(3*natoms,nmin))
    if(nmex.gt.0)allocate(mexstart(3*natoms,nmex))

    do i = 1,nmin
     if(surface(i).gt.0.and.surface(i).le.nstates)minstates(i) = surface(i)
     if(updatehess(i).eq.0)then
      h_update(i) = .false.
     else
      h_update(i) = .true.
     endif
     do j = 1,3*natoms
      minstart(j,i) = minguess((i-1)*3*natoms+j)
     enddo
    enddo

    do i = 1,nmex
     mexstates(i,:) = ci(i*2-1:i*2)
     do j = 1,3*natoms
      mexstart(j,i) = mexguess((i-1)*3*natoms+j)
     enddo
    enddo

   !-------------------------------------------------------
   ! Load geometries, print out geom/frequency information
   !-------------------------------------------------------
   case(3)
    isloop   = .false.
    loopst   = 1
    ngeoms   = 1
    geomfl   = ''
    calchess = .true.
    outputdiab=.false.
    outputdir= ''
    read(unit=INPUTFILE,NML=LOADGEOM)
    if(ngeoms<1)ngeoms=1
    if(.not.isloop)then
      loopst=1
      ngeoms =1
    end if
  end select

  close(unit=INPUTFILE)


  if(printlvl>0)print *,"Exiting readinput..."
  return
end SUBROUTINE readinput
