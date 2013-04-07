!------------------------------------------------------------------------
! testsurfgen
!------------------------------------------------------------------------
! A testing program for surfgen.   
!
! Surfgen is a collection of a set of programs that use quasi-diabatic
! Hamiltonians to reproduce high quality ab initio data for non-adiabatic
! chemistry simulations.!
!------------------------------------------------------------------------
! Purpose
! =======
! This program perform a series of tests for the surfgen subroutines
! The purpose of this program is to ensure the completeness and reliability
! of the program, and to check any inconsistency between versions
!
! Xiaolei Zhu,  Yarkony group, Johns Hopkins University
! April 2013
program testsurfgen
    implicit none

    ! print information
    call PrintHeader()

    ! initialize Hd
    print *,"Initializing surfgen"
    call readginput()
    call initialize(0)

    call testCoord(1d-5)
end program testsurfgen
!---------------------------------------------
!  printHeader
!---------------------------------------------
SUBROUTINE printHeader()
    IMPLICIT NONE
END SUBROUTINE printHeader


!---------------------------------------------
!  testCoord
!---------------------------------------------
!  Testing internal coordinate definitions for each coordinates 
!  that are defined in coord.in
SUBROUTINE testCoord(step)
    use hddata, only : ncoord
    use progdata, only:natoms,coordmap,CoordSet
    IMPLICIT NONE

    double precision,intent(in) :: step
    DOUBLE PRECISION,dimension(3*natoms)::geom
    INTEGER :: I,j,k,m,n
    double precision :: dgeom(3*natoms),cgrad, bm2(ncoord,3*natoms), &
                        bmat(ncoord,3*natoms),igeom(ncoord)
    integer,parameter :: maxstep=2
    double precision,parameter :: fdcoef(maxstep)=[2d0/3,-1d0/12]

    print *,"Testing Coordinate Gradients."

    call random_number(geom)
    geom = geom*5
    print *,"initial Cartesian geometry: "
    j=1
    do i=1,natoms
        print "(3F12.5)",geom(j:j+2)
        j=j+3
    end do

    call buildWBmat(geom,igeom,bmat)
    do i=1,ncoord
    m=coordmap(i,1) !index of set
    n=coordmap(i,2) !index in set
    print *,""
    print *,"Testing Coordinate ",i," (set",m,",ind",n,")"
    print *,"  Coord Type=",CoordSet(m)%Type,"  , Scaling Mode=",CoordSet(m)%Scaling
    print *,"  Atoms: ",CoordSet(m)%coord(:,n)
    print "(A)","   Cart# Analytical Numerical  Difference "
    do j=1,3*natoms
    ! calculating numerical gradient of coord i with respect to cartesian j
    cgrad = 0d0
    do k=1,maxstep
    dgeom = geom
    dgeom(j)=dgeom(j)+step*k
    call buildWBmat(dgeom,igeom,bm2)
    cgrad = cgrad+fdcoef(k)*igeom(i)
    dgeom = geom
    dgeom(j)=dgeom(j)-step*k
    call buildWBmat(dgeom,igeom,bm2)
    cgrad = cgrad-fdcoef(k)*igeom(i)
    end do!k
    cgrad=cgrad/step
    if(abs(cgrad-bmat(i,j))>1d-8)  &
    print "(3x,I4,2x,3E11.4)",j,bmat(i,j),cgrad,cgrad-bmat(i,j)
    end do!j=1,3*natoms
    end do!i=1,ncoord
END SUBROUTINE testCoord


!-----------------------------------------------------------------------------------
! read general input from input file (surfgen.in) 
SUBROUTINE readginput()
    use progdata
    use hddata, only: initGrps,ncoord,nstates,order,getFLUnit
    use CNPI, only: irrep,GrpPrty,GrpSym,nSymLineUps
    IMPLICIT NONE
    INTEGER                         :: i,j,ios,k
    INTEGER                         :: nGrp,jobtype
    INTEGER,dimension(10)           :: surface,updatehess
    INTEGER,dimension(50)           :: ci,atmgrp
    INTEGER,DIMENSION(20,MAX_ALLOWED_SYM)        :: groupSym,groupPrty

    NAMELIST /GENERAL/      jobtype,natoms,order,nGrp,groupsym,groupprty,&
                            printlvl,inputfl,atmgrp,nSymLineUps,cntfl

    nSymLineUps = 1
    natoms     = 2
    printlvl   = 1
    inputfl    = ''
    print *,"Entering readginput()."
    !----------- READ GENERAL INPUT ----------------!
    open(unit=INPUTFILE,file='surfgen.in',access='sequential',form='formatted',&
    IOSTAT=ios,POSITION='REWIND',ACTION='READ',STATUS='OLD')
    if(ios/=0) print *,"readinput: cannot open file surfgen.in.  IOSTAT=",ios

    read(unit=INPUTFILE,NML=GENERAL)
    if(printlvl>0)print *,"    readinput():  Control parameters read from surfgen.in"

    call genAtomList(atmgrp)

    call readIrreps()
    if(allocated(GrpSym))deallocate(GrpSym)
    if(allocated(GrpPrty))deallocate(GrpPrty)
    allocate(GrpSym(nGrp,nSymLineUps))
    allocate(GrpPrty(nGrp,nSymLineUps))
    GrpSym=groupsym(1:nGrp,1:nSymLineUps)
    GrpPrty=groupprty(1:nGrp,1:nSymLineUps)
    call initGrps(nGrp,irrep(GrpSym(:,1))%Dim)

end SUBROUTINE readginput
