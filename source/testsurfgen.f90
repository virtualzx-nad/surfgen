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
    use makesurfdata, only: initMakesurf 
    implicit none

    ! print information
    call PrintHeader()

    ! initialize Hd
    print *,"Initializing surfgen"
    call readginput()
    call initialize(1)

    ! test coordinate gradients
    call testCoord(1d-5)

    ! test gradients of Hd and Lagrangian
    call readdisps()
    call initMakesurf

    call testHd(5 ,1D-5)

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
    logical    :: problem, pcoord

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
    problem  = .false.
    do i=1,ncoord
      pcoord = .false.
      m=coordmap(i,1) !index of set
      n=coordmap(i,2) !index in set
      print *,""
      print "(3(A,I5),A)","Testing Coordinate ",i," (set",m,",ind",n,")"
      print "(2(A,I5))","  Coord Type=",CoordSet(m)%Type,"  , Scaling Mode=",CoordSet(m)%Scaling
      print "(A,4I5)","  Atoms: ",CoordSet(m)%coord(:,n)
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
        if(abs(cgrad-bmat(i,j))>1d-8)then
            print "(A)","   Cart# Analytical Numerical  Difference "
            print "(3x,I4,2x,3E11.4)",j,bmat(i,j),cgrad,cgrad-bmat(i,j)
            problem = .true.
            pcoord  = .true.
        end if 
      end do!j=1,3*natoms
      if(.not. pcoord)print "(A)", "  No problem."
    end do!i=1,ncoord
    if(problem)then
         stop "Coordinate test failed! Values and gradients inconsistent."
    else
        print *,"Coordinate test finished without problem. "
        print *,""
    end if
END SUBROUTINE testCoord


!-----------------------------------------------------------------------------------
! read general input from input file (surfgen.in) 
SUBROUTINE readginput()
    use progdata
    use hddata, only: initGrps,ncoord,nstates,order,getFLUnit,CpOrder
    use CNPI, only: irrep,GrpPrty,GrpSym,nSymLineUps
    IMPLICIT NONE
    INTEGER                         :: i,j,ios,k
    INTEGER                         :: nGrp,jobtype
    INTEGER,dimension(10)           :: surface,updatehess
    INTEGER,dimension(50)           :: ci,atmgrp
    INTEGER,DIMENSION(20,MAX_ALLOWED_SYM)        :: groupSym,groupPrty

    NAMELIST /GENERAL/      jobtype,natoms,order,nGrp,groupsym,groupprty,&
                            printlvl,inputfl,atmgrp,nSymLineUps,cntfl,CpOrder,&
                            basisfl

    nSymLineUps = 1
    CpOrder=-1
    natoms     = 2
    printlvl   = 1
    inputfl    = ''
    basisfl    = ''
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

    if(CpOrder==-1)CpOrder=order
PRINT *,"Maximum order of off-diagonal coupling blocks:",CpOrder
    CALL readMAKESURF(INPUTFILE)

end SUBROUTINE readginput

!------------------------------------------------------------------------------------------------------
! subroutine for testing gradients for Hd and/or gradients for Lagrangian
! An eight point 8th order finite difference scheme is used to evaluate the gradients
!
SUBROUTINE testHd(ntest,disp)
    USE progdata, only : natoms, printlvl
    USE hddata, only : nstates, getHdvec,updateHd,linearizeHd
    USE makesurfdata
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN)                  :: disp
    INTEGER,INTENT(IN)                            :: ntest
    double precision  :: hvec(ncons+nex),hvec_raw(ncon_total)
    integer           ::  i,  n, j, k, l, prtl
    double precision  ::  Ea(nstates), dHa(3*natoms,nstates,nstates)   ,&
    hmat(nstates,nstates), dcgrads(3*natoms,nstates,nstates)
    double precision,dimension(3*natoms) ::  cgeom,dgeom,cgeom0
    double precision                     :: hnorm, ind, maxratio,maxrLag
    double precision, dimension(ncons+nex)   ::  asol0, asol, dsol
    double precision,external            :: dnrm2
    double precision  :: DNumHd(nstates,nstates), DAnaHd(nstates,nstates),& ! Numerical and analytical gradients
                                        DNumHa(nstates),DAnaHa(nstates),  & ! for Hd in both representations
                                                            DNumL, DAnaL    ! for Hd and Lagrangian
    double precision  :: lag, dCi(ncons), dL(nex)
    double precision, dimension(nex,ncons)                      :: jaco
    integer, parameter:: ndisp=4
    double precision  :: coef(ndisp)=(/ 4/5D0 , -1/5D0, 4/105D0, -1/280D0 /)
    double precision  :: lagval(-ndisp:ndisp)
!    double precision  :: h0, ckl_sum(npoints,nstates,nstates), ckl_store(npoints,nstates,nstates),dij(nstates,nstates)
    call getHdvec(hvec_raw,coefmap,ncon_total)
    call tranHd('F',hvec_raw,hvec)   ! convert primitive h expansion to block orthogonal basis expansions

    prtl = printlvl
    printlvl = printlvl-2
    hnorm = dnrm2(ncons,hvec,int(1))
    print *,"Testing gradients of Hd along 3 random pathes"
    print *,"  Cartesian displacement size :", disp
    call init_random_seed()
    maxratio = 0d0
    maxrLag  = 0d0
    do i=1,ntest
        print *,"  displacement",i," :  evaluating analytical gradients"
        ! initialize random Hd, geometry and displacements
        ! Hd is initialized as a random offset from existing Hd
        call random_number(asol0)
        asol0 = asol0-5d-1
        asol0 = asol0/dnrm2(ncons,asol0,int(1))*hnorm*5D-3+hvec
        ! geometry is taken as a random offset from a random data point
        call random_number(cgeom0)
        cgeom0=(cgeom0-5d-1)*1d-1
        call random_number(ind)
        cgeom0=cgeom0+dispgeoms(int(ind*npoints+1))%cgeom
        ! generate a random direction for geometry displacement
        call random_number(dgeom)
        dgeom = (dgeom-5d-1)
        dgeom = dgeom/dnrm2(3*natoms,dgeom,int(1))*disp

        ! generate a random displacement in coefficient space
        call random_number(dsol)
        dsol = dsol - 5D-1
        dsol = dsol/dnrm2(ncons,dsol,int(1))*disp

        call tranHd('B',asol0,hvec_raw)                ! convert hd to nascent basis
        call updateHd(hvec_raw,CoefMap,ncon_total)     ! convert vector hd to list form
        call linearizeHd()

        ! analytical evaluation of Hd gradients
        call getCartHd(cgeom0,Ea,dHa,hmat,dcgrads)
        DAnaHd = 0d0
        DAnaHa = 0d0
        do j=1,nstates
            DAnaHa(j) = dot_product(dHa(:,j,j),dgeom)
            DAnaHd(j,j) = dot_product(dcgrads(:,j,j),dgeom)
            do k=j+1,nstates
              DAnaHd(j,k) = dot_product(dcgrads(:,j,k),dgeom)
              DAnaHd(k,j) = DAnaHd(j,k)
            end do !k=j+1,nstates
        end do!j=1,nstates
        DAnaHd = DAnaHd/disp
        DAnaHa = DAnaHa/disp

        ! analytical evaluation of Lag gradients
        call updateEigenVec(asol0,.false.)
        CALL getCGrad(asol0,dCi,dL,lag,jaco)
        lagval(0) = lag
        DAnaL = dot_product(dsol(1:ncons),dCi)
        DAnaL = DAnaL/disp

        print *,""
        print *,"Generating Hd Numerical Gradients (Ana-Num=Diff)"
        ! numerical evaluation of Hd gradients
        DNumHa = 0d0
        DNumHd = 0d0
        do n = 1,ndisp
            cgeom = cgeom0 + dgeom*n
            call getCartHd(cgeom,Ea,dHa,hmat,dcgrads)
            DNumHa = DNumHa+coef(n)*Ea
            DNumHd = DNumHd+coef(n)*hmat
            cgeom = cgeom0 - dgeom*n
            call getCartHd(cgeom,Ea,dHa,hmat,dcgrads)
            DNumHa = DNumHa-coef(n)*Ea
            DNumHd = DNumHd-coef(n)*hmat
        end do
        DNumHa = DNumHa/disp
        DNumHd = DNumHd/disp
        !print out comparisons
        print *,"Diabatic Representation  : "
        do j=1,nstates
            do k=j,nstates
                print "(2(A,I4),3(A,E14.7))","  Block(",j,",",k,"):  ",DAnaHd(j,k),&
                            " - ",DNumHd(j,k)," = ",DAnaHd(j,k)-DNumHd(j,k)
                maxratio = max(maxratio,abs((DAnaHd(j,k)-DNumHd(j,k))/DNumHd(j,k)))
            end do!k=1,nstates
        end do!j=1,nstates
        print *,"Adiabatic Representation  : "
        do j=1,nstates
            print "(A,I4,3(A,E14.7))","  State ",j,":   ",DAnaHa(j),&
                            " - ",DNumHa(j)," = ",DAnaHa(j)-DNumHa(j)
            maxratio = max(maxratio,abs((DAnaHa(j)-DNumHa(j))/DNumHa(j)))
        end do!j

        ! numerical evaluation of Lag gradients
        print *,""
        Print *,"Generating Lagrangian Gradients (Ana-Num=Diff)"
        DNumL = 0d0
        do n = 1,ndisp
            asol = asol0 + dsol*n
            call updateEigenVec(asol,.true.)
            CALL getCGrad(asol,dCi,dL,lag,jaco)
            lagval(n) = lag
            DNumL = DNumL + coef(n)*lag
            asol = asol0 - dsol*n
            call updateEigenVec(asol,.true.)
            CALL getCGrad(asol,dCi,dL,lag,jaco)
            lagval(-n) = lag
            DNumL = DNumL - coef(n)*lag
        end do!n=1,ndisp
        DNumL = DNumL/disp

        ! Comparisons with previous largest error
        print "(3(A,E14.7))","   ",DAnaL," - ",DNumL," = ",DNumL-DAnaL
        maxrLag = max(maxrLag,abs((DNumL-DAnaL)/DNumL))
        if(abs(DNumL-DAnaL)/DNumL>1d-4)then
            print *,"Lagrangian has large error."
            DNumL = -lagval(-3)/60+lagval(-2)*3/20-lagval(-1)*3/4 + &
                     lagval(3)/60 -lagval(2) *3/20+lagval(1) *3/4
            DNumL = DNumL/disp
            print "(A,E14.7)","6 instead of 8 order numerical gradients:",DNumL
        end if
    end do!i=1,ntest
    print "(A,F10.5,A)"," Maximum relative error for Hd gradients: ",maxratio*100,"%"
    print "(A,F10.5,A)"," Maximum relative error for Lagrangian: ",maxrLag*100,"%"
    print *,""
    if(maxratio>1d-6)stop"Hd gradient test failed."
    if(maxrLag>1d-3) stop"Lagrangian test failed."
    print *,"Hd gradient and Lagrangian tests finished."
    printlvl = prtl
END SUBROUTINE testHd
