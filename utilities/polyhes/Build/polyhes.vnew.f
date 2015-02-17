c========================================================
c program to find conical intersections                 =
c written by David R. Yarkony                           =
c modified by cthree-40                                 =
c                                                       =
c Yarkony Group                                         =
c The Johns Hopkins University                          =
c========================================================
      program polyhes
      
      implicit none                                               
c     data  = size of core memory to use
c     fdate = date of program execution
c     ifail = error variable
      integer        ncor, ifail, iap, npflg
      character*24   fdate
      character*100  workdir, cpgeomcmd
      data           ncor/1500000/ 
      data           ifail/0/

c     clm: unsure what these variables do. set at DRY suggestion (02-13-15)
      iap=1
      npflg=2 

c     set up work directory
c      call setupjob( workdir )
c      call chdir( trim(workdir) )
c      write ( cpgeomcmd, "(a)" ) "cp ../geom.polyhes ../"
c      call system( cpgeomcmd )
c      call setupfiles()
      
      write(6,*)' Polyhes clock starts: ',fdate()
c     print allocation information
      write(6,"(1x,A,I10,A)") " Allocating ", 2*ncor, " bytes"
c     call main polyhes driver
      call vmain(ncor)
      write(6,*)' Polyhes clock ends:   ',fdate() 
      
      end program
c----------------------------------------------------------------------
c----------------------------------------------------------------------
      subroutine vmain( ncor )

c     main polyhes subroutine.
      implicit none
c     Input:
c      ncor  = bytes of memory to be allocated in ia and a
      integer         ncor
      real*8          a( ncor ), ia( ncor )
      integer         natoms, mfixdir, matoms,     m3,     n3,  ndim 
      integer         nd1,        n31,    n3s,  nprec,  nacme,   scr 
      integer         hesss,   egradh, ggradh, xgradh,   hess,  grad   
      integer         eval,      root,   vect,  macme,  lacme, rvect  
      integer         bvect,     scr1,   scr2,  tvect,         vvect
      integer         lambda,  kambda,  ichar, isymtr,  icode,  itop
      integer         intowp, maxiter,    igo
      integer         i, cmdios, close_stat
      character*10    progfli
      character*100   cpprogcmd, copygmcmd
c     Files:
c      i74  = input   geometry file
c      i75  = output  geometry file
c      i78  = scratch geometry file
      integer         i74, i75, i78
      parameter       ( i74 = 74 )
      parameter       ( i75 = 75 )
      parameter       ( i78 = 78 )


c     rewind unit 74 (geometry input file)
      rewind i74
c     get number of atoms from unit 74
      call numatm( i74, 0, natoms)
c     if too many atoms, stop
      if ( natoms .gt. 10 ) stop ' natoms > 10 '
c     mfixdir is
      mfixdir = 30
c     matoms is maximum of natoms and 10
      matoms = max0( natoms, 10)
      m3     = natoms * 3
      n3     = matoms * 3
      ndim   = n3 + mfixdir
      nd1    = ndim * (ndim + 2)
      n31    = (n3 + 1) * n3
      n3s    = n3 * n3
      nprec  = intowp(1)

      write(6,1001) matoms, n3, n31, n3s, nprec
c set pointers for ia(ncor) array
      nacme  = 1
      scr    = nacme  + (n3      * nprec)
      hesss  = scr    + (nd1     * nprec)
      egradh = hesss  + (nd1 / 2 * nprec) 
      ggradh = egradh + (n31 * 2 * nprec)
      xgradh = ggradh + (n31 * 2 * nprec)
      hess   = xgradh + (n31 * 2 * nprec)
      grad   = hess   + (nd1     * nprec)
      eval   = grad   + (ndim    * nprec) 
      root   = eval   + (n3      * nprec)
      vect   = root   + (n3      * nprec)
      macme  = vect   + (n3s     * nprec)
      lacme  = macme  + (n3      * nprec) 
      rvect  = lacme  + (n3      * nprec)
      bvect  = rvect  + (n3  * 3 * nprec)
      scr1   = bvect  + (n3s     * nprec)
      scr2   = scr1   + (nd1     * nprec)
      tvect  = scr2   + (nd1     * nprec)
      vvect  = tvect  + (n31     * nprec)
      lambda = vvect  + (n31     * nprec)
      kambda = lambda + (mfixdir * 2 * nprec)
      ichar  = kambda + (mfixdir)
      isymtr = ichar  + (n3 * 11 * nprec)
      itop   = isymtr + (n3s     * nprec) 

      write(6,1000) nacme, scr, hesss, egradh, ggradh, xgradh, hess,
     1              grad, eval, root, vect, macme, lacme, rvect,
     3              bvect, scr1, scr2, tvect, vvect, lambda, kambda,
     2              ichar, isymtr, itop

c exit if not enough memory was allocated
      if ( itop .gt. ncor ) then
         write(6,"(1x,a)") ' Insufficient memory in polyhes'
         write(6,"(1x,a)") ' Please adjust ncor and recompile.'
         stop ' Exiting...'
      end if

c loop for intersection searching
      maxiter = 2
      i   = 0
      igo = 0

 10   if ( igo .eq. 0 .and. i .lt. maxiter ) then
         i = i + 1
         write(6,"(1x,' Iteration ',i3)") i
         call polyhsd(     ia(nacme),  ia(nacme),    ia(scr), ia(hesss),
     &        ia(egradh), ia(ggradh), ia(xgradh),   ia(hess),  ia(grad),
     &        ia(eval),     ia(root),   ia(vect),  ia(macme), ia(macme),
     &        ia(lacme),   ia(lacme),  ia(rvect),  ia(bvect),  ia(scr1),
     &        ia(scr2),    ia(tvect),  ia(vvect), ia(lambda),ia(kambda),
     &        ia(ichar),  ia(isymtr),   n3,   m3,     matoms,   mfixdir,
     &        i, igo )
         print *, "igo = ", igo
c     close geometry files
         close(unit = i74, iostat = close_stat )
         if ( close_stat .ne. 0 ) stop " ** Could not close file 74. **"
         close(unit = i75, iostat = close_stat )
         if ( close_stat .ne. 0 ) stop " ** Could not close file 75. **"
c     copy geometry files over
         write(copygmcmd, "(a)") "cp fort.75 fort.74"
         write(6, "(a)") copygmcmd

         call system( copygmcmd )
         rewind i74
         go to 10
      else if ( igo .eq. 1 ) then
c        calculation converged
         write(6,"(1x,A)") "Calculation converged. Intersection found."
c        write out final geometry in Columbus format
         call wrtcolgm( i75, natoms )
      else if ( igo .eq. -1 ) then
c        calculation failed
         write(6,"(1x,A)") "Calculation failed."
      else if ( igo .eq. 0  .and. i .ge. maxiter ) then
c        calcuation did not converge
         write(6,"(1x,A)") "Calcuation did not converge."
      else
         write(6,"(1x,A,i3,A,i5)") 
     &   "What's this scenario?? igo=",igo,", i=",i
      end if


      return
 1000 format(/2X,'NACME SCR HESSS EGRADH GGRADH XGRADH HESS GRAD',/8I6,
     1       /2X,'EVAL ROOT VECT MACME LACME RVECT BVECT   SCR1',/8I6, 
     2       /2X,'SCR2 TVECT VVECT LAMBDA KAMBDA ICHAR ISYMTR ITOP',
     3       /8I6)
 1001 format(/2X,'MEMORY LOCATION PARAMETERS:MATOMS=',I3,' N3=',I5,
     1      ' N31=',I5,' N3S=',I5,' NPREC=',I1)

      end subroutine vmain
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine polyhsd( nacme,  rnacme,   scr,  hesss, egradh, ggradh,
     &          xgradh,   hess,     grad,  eval,   root,   vect,  macme,
     &                    rmacme,  lacme, rlacme, rvect,  bvect,   scr1,
     &                    scr2,    tvect,  vvect,lambda,  klass,   char,
     &                    symtrn,     n3,     m3,matoms, mfixdir, iter,
     &                    igo )
      
      implicit real*8 (a-h,o-z)
      character*6    verb
      character*19   kmes(3)
      logical        same
      integer        simult,   across, gmfile,    restart, accel, saddle 
      integer        atom(10), state1, state2, killer(30),  cdif, origin
      integer        fixdir(6,30),     klass(mfixdir),      actual
      integer        igo, iter
      real*8         mass(10),       gm(3,10), disp(30),  dampv(30) 
      real*8         xgradr(100),      scr(2), hesss(2), egradh(N3,2)
      real*8         ggradh(n3,2), dipolr(50),  hess(2),   grad(n3)
      real*8         eval(n3),       root(n3),  vect(2),   norm
      real*8         macme(3,matoms), rmacme(n3), lambda(mfixdir,2)
      real*8         rvect(2), bvect(2), scr1(2), scr2(2), valdir(30)
      real*8         tvect(2), vvect(2), nacme(3,matoms), rnacme(n3)
      real*8         lacme(3,matoms), rlacme(n3), xgradh(n3,2)
      real*8         symtrn(m3,m3), char(m3,11)
c    CLM: variables for interfacing with surfJay
c         polyhes has inequality test NATOMS < 10, 
c         so JGRADS has dimension (3*natoms,nstates,nstates)
c         We have chosen to assume 3 states. Other arrays won't be
c         used: SCGD1 SCENG
      real*8         jgrads(15,3,3), scgd1(15,3,3), sceng(5)
      dimension      grad1(3), grad2(3), grad3(3), grd(3,3), ipflag(5)
      equivalence (  grd(1,1), grad1 )
      equivalence (  grd(1,2), grad2 )
      equivalence (  grd(1,3), grad3 )
      common   /vibord /  ivbord(30), ivbord2(30)
      common   /diatat /  ndiat(10)
      common   /chess  /  norr
      common   /pscent /  maxnew, psc(30,10), dpsc(30,10), ipsc(30,10)

      namelist /path   /  gm, natoms
      namelist /nacint /  ncalc, nacout, mass, gm, idisp, disp, simult,
     &                    newton, valdir, idbg, zeig, scale, restart,
     &                    accel, ihess, damp, dampv, fixdir, conv,
     &                    maxit, mhess, atom, kz, cdif, ipflg, iheseq1,
     &                    tolrot, actual, nsfig, kscale, norr, lineq,
     &                    dpsc, ipsc, ipflag
      namelist /cross  /  state1, state2, motion, across, actual, 
     &                    derivc
c         newton-raphson convergence criteria
      data  conv/1.0d-06/, maxit/20/, nsfig/5/, ipflag/5*0/

      data  atom/1,2,3,4,5,6,7,8,9,10/, fixdir/180*0/, valdir/30*0.D0/
      data  ncalc/1/, nacout/2/, simult/0/, intc/1/, iorg/2/
      data  ipflg/-1/, idbg/-1/, lineq/0/
      data  gmfile/1/, restart/0/, accel/0/, damp/0.0D+0/
      data  idisp/0/, zeig/5.D-6/, scale/1.D0/, eigrot/1.0d-3/
      data  tolrot/1.0d-3/
      data  n32/32/, ihess/0/, saddle/0/, zero/0.D0/, mhess/1/
      data  iciu1, iciu2/31,41/
      data  kz/3/, cdif/0/, iheseq1/0/, kscale/0/
      data  kmes/' entire ', 'entire - geometric', ' energy part of '/
      data  grd/9*1/, nprop/4/
C
C      ITEST  NE 0  NO TAPE10
C      IDISP  LT 0  NEW GEOM ONLY
C             EQ 0  NO NEW GEOM INTERNAL NACME ONLY
C             GT 0  NEW GEOM AND INTERNAL NACME
C      DAMPV        0 < DAMPV(I) < 1  RETAIN DAMPV(I) OF GRADIENT IN 
C                   DIRECTION OF ITH EIGENVALUE
C      DAMP         SET ALL ZERO ENTRIES OF DAMPV TO DAMP
C      SIMULT EQ 1  ALL DISPLACEMENTS AT ONCE
C             EQ 0  INDIVIDUAL DISPLACEMENTS
C                   IDISP=-N NAMELISTED 
C      NEWTON EQ 1  USE HESS AND GRAD TO DETERMINE DISP
C                   FOR NEW GEOM ( ITEST=0 SIMULT=1 NEWTON=1)
C      ACROSS EQ 1  GRADIENT IS (E(J)-E(I))*GRAD(K)
C                   GET E(J) AND E(I) FROM SAXE TAPE12
C      RESTART = 1  NEWTON=1 GET GM AND GRAD FROM PROGRESS  TAPE
C                   TO GENERATE NEXT STEP VIA NEWTON
C                   IHESS=1 GET GM AND IDISP FROM PROGRESS  TAPE
C                   TO GENERATE NEXT STEP FOR HESSIAN CONSTRUCTION
C      ACCEL  EQ 1  ACCELERATE SEARCH USING SCALE
C                   
C      IHESS  EQ 1  DRIVE HESSIAN CONSTRUCTION
C
C      MHESS  EQ 1  UPPER TRIANGLE DIAGONALIZED(DEFAULT)
C               -1  LOWER TRIANGLE DIAGONALIZED
C                0  AVERAGE OF OFF DIAGONAL DIAGONALIZED
C      LAMBDA(1,*)  LAGRANGE MULTIPLIER USED TO CONSTRUCT RHS WITH 
C      LAMBDA(2,*)  USED TO CONSTRAIN INTERNAL COORDINATE
C     
C      FIXDIR       FIXDIR(1,*)=I FIXDIR(2,*)=J - R(I,J) FIXED AT VALDIR(I)
C                   FIXDIR(1,*)=I FIXDIR(2,*)=J  FIXDIR(3,*)=K - 
C                   ANGLE(I,J,K ) FIXED AT VALDIR(I)
C                   FIXDIR(I,J,K,0,0,L) ANGLE(I,J,K,L)
C                   TO INTRODUCE THIS OPTION WITH THE HESSIAN USER MUST 
C                   REGENERATE WORKING HESSIAN USING RESTART
C      IHESEQ1      =1 USE A UNIT HESSIAN ( W PORTION ONLY)
C
C      ACTUAL       =1 FOR AN ACTUAL CROSSING 
C      KSCALE       =  0  SCALE ENTIRE GRADIENT
C                   =  1  SCALE GRADIENT WITHOUT GEOMETRIC CONTRAINTS
C                   >  1  SCALE ENERGY PART OF GRADIENT ONLY
C
C                      FILE USAGE
C      I70       PROPERTY UNIT                  70 
C      I73       DEPENDENT INPUT UNIT           73
C      I74       CURRENT GEOMETRY FILE          74
C      I75       NEXT GEOMETRY FILE             75
C      I76       SYMMETRY EQUIVALENT CENTER     76
C      I77       SEC NEXT GEOMETRY              77
C      I99       TEMPORARY UNIT                 99
C      N32       PROGRESS FILE                  32
C
C      ICIU1     ALCHEMY CI HEADER FILE         31
C      ICIU2     ALCHEMY CI HEADER FILE         41
C      ITAP12    GUGA CI FILE                   12
C      ITAP13    GUGA CI FILE                   13
C      ITAP10    SAXE TAPE 10    E(A)           10
C      ITAP11    SAXE TAPE 11    G(A)           11
C      ITAP14    SAXE TAPE 14    H(A)           14
C
C         UNITS
      ITAP10=10
      ITAP11=11
      ITAP14=14
      INP=5
      IOUT=6
      I70=70
      I73=73
      I74=74
      I75=75
      I76=76
      I77=77
      I99=99
C
      PI=ACOS(-1.D0)
      NORR=1
      DO 6 I=1,MFIXDIR
      LAMBDA(I,1)=ZERO
    6 LAMBDA(I,2)=ZERO
      RMAG=ZERO
C
      WRITE(IOUT,1001)
      REWIND I74
      CALL NUMATM(I74,0,NATOMS)
CCLM: this print is for debugging
      PRINT *, "NATOMS = ", NATOMS
      DO 2 I=1,NATOMS*3
      DAMPV(I)=ZERO
      DISP(I)=0.1D-2
    2 IVBORD(I)=I
C
      REWIND I73
CCLM: this print is for debugging
      PRINT *, "CALLING LOCATOR..."
      CALL LOCATOR(I73,10H**NACINT  ,IERR1)
      IF(IERR1.EQ.0)READ(I73,NACINT,ERR=920)
C
      REWIND INP
CCLM: this print is for debugging
      PRINT *, "CALLING LOCATOR..."
      CALL LOCATOR(INP,10H**NACINT  ,IERRI)
      IF(IERRI*IERR1.NE.0)THEN
      WRITE(IOUT,1030)INP,I73
      CALL EXIT(100)
      ENDIF
      IF(IERRI.EQ.0)READ(INP,NACINT,ERR=921)
      IF(IERRI.NE.0)WRITE(IOUT,1029)INP
      IF(IERR1.NE.0)WRITE(IOUT,1029)I73
      SCLNAC=SCALE
C
      REWIND I73
CCLM: this print is for debugging
      PRINT *, "CALLING REPOSN..."
      CALL REPOSN(I73,8H***CROSS,IOUT)
      READ(I73,CROSS)
C
      SAME=ACTUAL.NE.0
      NCGNT=1
      IF(SAME)NCGNT=2
C
      IF(SAME)WRITE(IOUT,1044)
      WRITE(IOUT,1041)MOTION
      NTDEGF=NATOMS*3
      NDEGF=NATOMS*3-3
      NDEGFI=NDEGF-3
      WRITE(IOUT,1057)KMES(KSCALE+1)
C
      IF(MOTION.EQ.2)THEN
      NDEGF=NATOMS*2-2
      NDEGFI=NDEGF-1
      ENDIF
      IF(MOTION.EQ.1)THEN
      NDEGF=NATOMS-1
      NDEGFI=NDEGF
      ENDIF
C
      NROT=  NDEGF-NDEGFI
C
      maxnew=0
      NFIXDIR=0
      DO 3 I=1,30
      IF(FIXDIR(1,I).EQ.0)GO TO 3
      NFIXDIR=I
      do 44 j=1,6
   44 if(fixdir(j,i).gt.natoms)maxnew=max0(fixdir(j,i),maxnew)
    3 CONTINUE
      maxnew=maxnew-natoms
      leq=0
      if((nfixdir.eq.1).and.(ndegfi.eq.3).and.(lineq.eq.1))leq=1
      if(maxnew.gt.0)call mkpscent(natoms,maxnew,psc,dpsc,ipsc,mass,
     x    iout)
C
      NINTI=    NDEGFI+NFIXDIR+NCGNT
      NINT=     NINTI+NROT
      NDEGFL=   NDEGF+NCGNT
      NDEGFIL=  NDEGFI+NCGNT
C
      NINT2=    NINT*(NINT+1)/2
      NINTI2=   NINTI*(NINTI+1)/2
      NDEGFL2=  NDEGFL*(NDEGFL+1)/2
      NDEGFIL2= NDEGFIL*(NDEGFIL+1)/2
      NDEGFI2=  NDEGFI*(NDEGFI+1)/2
      NDEGF2=   NDEGF*(NDEGF+1)/2
C
      WRITE(IOUT,1043)NDEGF,NROT,NFIXDIR,NINTI,NINT
      IF(NFIXDIR.EQ.0)GO TO 5
      DO 4 I=1,NFIXDIR
      KLASS(I)=1
      IF(FIXDIR(3,I).NE.0)GO TO 14
      WRITE(IOUT,1045)I,FIXDIR(1,I),FIXDIR(2,I),VALDIR(I)
      GO TO 4
   14 CONTINUE
      IF(FIXDIR(4,I).NE.0)GO TO 15
        IF(FIXDIR(6,I).NE.0)THEN
        KLASS(I)=5
        WRITE(IOUT,1062)I,FIXDIR(1,I),FIXDIR(2,I),FIXDIR(3,I),
     X               FIXDIR(6,I),VALDIR(I)
        GO TO 4
        ENDIF
      KLASS(I)=2
      WRITE(IOUT,1046)I,FIXDIR(1,I),FIXDIR(2,I),FIXDIR(3,I),VALDIR(I)
      GO TO 4
   15 CONTINUE
      IF(FIXDIR(6,I).NE.0)GO TO 17
      KLASS(I)=3
      WRITE(IOUT,1047)I,FIXDIR(1,I),FIXDIR(2,I),FIXDIR(3,I),FIXDIR(4,I)
      GO TO 4
   17 CONTINUE
      KLASS(I)=4
      WRITE(IOUT,1048)I,FIXDIR(1,I),FIXDIR(2,I),FIXDIR(3,I),
     X    FIXDIR(4,I),FIXDIR(5,I),FIXDIR(6,I)
    4 CONTINUE
    5 CONTINUE
C
      DO 7 I=1,NATOMS
    7 NDIAT(I)=ATOM(I)
      WRITE(IOUT,1012)(I,MASS(I),I=1,NATOMS)
      WRITE(IOUT,1031)(NDIAT(I),I=1,NATOMS)
      CALL IIDOF(NATOMS,MOTION,NDIAT,KILLER,KZ,IOUT,ORIGIN)
C
      MEWTON=NEWTON
C
      IF(RESTART.EQ.0)GO TO 130
      IF(NEWTON.EQ.0)GO TO 135
      WRITE(IOUT,1027)N32
      CALL GETPRG(N32,NATOMS,NDEGF,GM,RNACME,RMACME,RLACME,
     X           EJMEIA,ES1,LAMBDA,IOUT,NFIXDIR)
      CRTERN=EJMEIA
      RLAG=LAMBDA(NFIXDIR+1,1)
      RMAG=LAMBDA(NFIXDIR+2,1)
      WRITE(IOUT,1021)(J,(GM(I,J),I=1,3),J=1,NATOMS)
      CALL GETPROP1(N32,DIPOLR,XGRADR,NDEGF,IOUT)
      GO TO 450
  135 CONTINUE
      CALL GETHPRG(N32,NATOMS,GM,IGO,IOUT)
      WRITE(IOUT,1034)N32,IGO
      MAXDO=NDEGF+1
      IF(CDIF.GT.0)MAXDO=NDEGF*2+1
      IF(IGO.LT.MAXDO)GO TO 430
      CALL GETDDH(N32,NDEGF,EGRADH,GGRADH,XGRADH,RNACME,RMACME,
     1            RLACME,IOUT,CDIF,LAMBDA,NFIXDIR,SAME) 
      RLAG=LAMBDA(NFIXDIR+1,1)
      RMAG=LAMBDA(NFIXDIR+2,1)
      GO TO 410
  130 CONTINUE
c
c     rewind unit 74 (input geometry file)
      rewind i74
c     read in geometry
      read(i74,path)
      write(iout,1022) i74
      go to 160
 150  continue
      write(iout,1023)
      write(iout,1021) (j,(gm(i,j),i=1,3),j=1,natoms)
 160  continue
c     initialize surfgen potential
      if ( iter .eq. 1 ) call initPotential()
c     evaluate surface
      jgrads=0d0
      call EvaluateSurfgen77(NATOMS,3,GM,EVAL,JGRADS,
     1     SCGD1,SCENG)
      print *, "surfjay call successful"
C      CALL RDT30(ICIU1,EVAL,NRTU1,IOUT,IERR)
      NCIU=ICIU1
      IF(IERR.NE.0)GO TO 900
      print *, "Evaluating intersection between states: ", STATE1
     & , " and ", state2
      EJMEIA=EVAL(STATE1)
      ES1=EVAL(STATE1)
C      IF(.NOT.SAME)CALL RDT30(ICIU2,EVAL,NRTU2,IOUT,IERR)
      NCIU=ICIU2
      IF(IERR.NE.0)GO TO 900
      EJMEIA=-(EJMEIA-EVAL(STATE2))
      ES2=EVAL(STATE2) 
      WRITE(IOUT,1026)ICIU1,ICIU2,EJMEIA
C
C      REWIND I74
C      READ(I74,PATH)
C      WRITE(IOUT,1022)I74
C CLM: removed go to 02-09-2015
C      GO TO 160
C  150 CONTINUE
C      WRITE(IOUT,1023)
C      WRITE(IOUT,1021)(J,(GM(I,J),I=1,3),J=1,NATOMS)
C  160 CONTINUE
C
C CLM: These calls unnecessary.
C          ENERGY  GRADIENTS
C      CALL RFILE(ITAP10)
C         GET GEOM AND GRADIENT POINTERS
C      CALL WREADW(ITAP10,ISTART,1,500+NCALC,JUNK)
C      IPNTG=ISTART+80
C      IPNTEN=IPNTG+INTOWP(NATOMS*3)
C      CALL WREADW(ITAP10,ENUC,INTOWP(1),IPNTEN,JUNK)
C      WRITE(IOUT,1035)ES1,ES2,ENUC
C      IPNTN=ISTART+60+NACOUT-1
C      CALL WREADW(ITAP10,IPNTN,1,IPNTN,JUNK)
C      WRITE(IOUT,1010)IPNTG,IPNTN
C         GET RAW GRADIENTS
C      IF(IPNTN.EQ.0)THEN
C      WRITE(IOUT,1024)ITAP10
C      CALL EXIT(100)
C      ENDIF
C      CALL WREADW(ITAP10,NACME,INTOWP(3*NATOMS),IPNTN,JUNK)
C CLM: rnacme is the energy gradient
      DO I=1,3*NATOMS
         RNACME(I) = JGRADS(I,STATE2,STATE2)
         print *, rnacme(i)
      END DO
      WRITE(IOUT,1003)'ENERGY',(RNACME(I),I=1,3*NATOMS)

C CLM: energy gradients
C      do i=1,3*NATOMS
C         rnacme(i) = JGRADS(i)
C      end do
C          ENERGY DIFFERENCE GRADIENTS
C      CALL RFILE(ITAP11)
C         GET GEOM AND GRADIENT POINTERS
C      CALL WREADW(ITAP11,ISTART,1,500+NCALC,JUNK)
C      IPNTG=ISTART+80
C      IPNTN=ISTART+60+NACOUT-1
C      CALL WREADW(ITAP11,IPNTN,1,IPNTN,JUNK)
C      WRITE(IOUT,1010)IPNTG,IPNTN
C         GET RAW GRADIENTS
C      IF(IPNTN.EQ.0)THEN
C      WRITE(IOUT,1024)ITAP11
C      CALL EXIT(100)
C      ENDIF
C      CALL WREADW(ITAP11,MACME,INTOWP(3*NATOMS),IPNTN,JUNK)
C CLM: compute energy difference gradients
      DO I=1,3*NATOMS
         RMACME(I) = JGRADS(I,STATE2,STATE2) - JGRADS(I,STATE1,STATE1)
      END DO
      WRITE(IOUT,1003)'E-DIF ',(RMACME(I),I=1,3*NATOMS)
C
C      CALL RZERO(LACME,3*NATOMS)
C      IF(SAME)THEN
C          H(1,2) GRADIENTS
C      CALL RFILE(ITAP14)
C         GET GEOM AND GRADIENT POINTERS
C      CALL WREADW(ITAP14,ISTART,1,500+NCALC,JUNK)
C      IPNTG=ISTART+80
C      IPNTN=ISTART+60+NACOUT-1
C      CALL WREADW(ITAP14,IPNTN,1,IPNTN,JUNK)
C      WRITE(IOUT,1010)IPNTG,IPNTN
C         GET RAW GRADIENTS
C      IF(IPNTN.EQ.0)THEN
C      WRITE(IOUT,1024)ITAP14
C      CALL EXIT(100)
C      ENDIF
C      CALL WREADW(ITAP14,LACME,INTOWP(3*NATOMS),IPNTN,JUNK)
C CLM: couplings 
      DO I=1,3*NATOMS
         RLACME(I) = JGRADS(I,STATE1,STATE2)
      END DO
      WRITE(IOUT,1003)'H(IJ) ',(RLACME(I),I=1,3*NATOMS)
C      ENDIF
C
   50 CONTINUE
      CALL GETSYMT(NATOMS,NATOMS*3,SYMTRN,CHAR,SCR,I76,IOUT,KSEC,
     X             ATOM(1),NTSYM,NTSYMNT,NOPS)
      REWIND I99
C
      CALL SECNT(NATOMS,NACME,I76,IOUT)
      WRITE(I99)(RNACME(I),I=1,NATOMS*3)
C
      CALL SECNT(NATOMS,MACME,I76,IOUT)
      WRITE(I99)(RMACME(I),I=1,NATOMS*3)
C
      IF(SAME)THEN
      CALL SECNT(NATOMS,LACME,I76,IOUT)
      WRITE(I99)(RLACME(I),I=1,NATOMS*3)
      ENDIF

      WRITE(IOUT,1002)(J,(GM(I,J),I=1,3),J=1,NATOMS)
      CALL DIST(NATOMS,GM,IOUT)
C
      IF(SAME)THEN
      DO 49 I=1,NATOMS*3
  49  SCR(I)=RNACME(I)+RMACME(I)/2.d0
      RNL=SQRT(DOT(LACME,LACME,NATOMS*3))
      RNM=SQRT(DOT(MACME,MACME,NATOMS*3))
      DLX=(RNL**2 - RNM**2)/(RNL**2 + RNM**2)
      SCRN=SQRT(DOT(SCR,SCR,NATOMS*3))
      DP=DOT(MACME,LACME,NATOMS*3)/(RNL*RNM)
      DP=ACOS(DP)
      GXHD2=RNL*RNM/2*SIN(DP)
      DP=DP/PI*180.
      RAVG1=DOT(SCR,LACME,NATOMS*3)/RNL
      RAVG2=DOT(SCR,MACME,NATOMS*3)/RNM
      PERC=sqrt(ravg1**2 +ravg2**2)/scrn*100
      RAVG1=RAVG1/(RNL*SCRN)
      RAVG2=RAVG2/(SCRN*RNM)
      WRITE(IOUT,1060)RNL,RNM/2.,DP,DLX,GXHD2,SCRN,PERC,
     x                RAVG2,RAVG1
      IF(NATOMS.EQ.3)THEN
      WRITE(IOUT,1059)
      CALL TOINTRL(MASS,GM,MACME,GRAD1,IOUT,'GIJ')
      CALL TOINTRL(MASS,GM,LACME,GRAD2,IOUT,'HIJ')
      CALL SCHMO(3,GRD,IOUT,IFAIL)
      WRITE(IOUT,1058)GRD
      ENDIF
      ENDIF
C
      CALL IDOF(NATOMS,NACME,SCR,KILLER)
      CALL IDOF(NATOMS,MACME,SCR,KILLER)
      IF(SAME)THEN
      CALL IDOF(NATOMS,LACME,SCR,KILLER)
      CALL SGNHIJ(N32,NDEGF,LACME,IOUT,CHANGE)
      ENDIF
C
      CRTERN=EJMEIA
C
      ICODE=0
      IF(NEWTON)71,81,91
   71 CONTINUE
      STOP ' ONLY NR SEARCHES ALLOWED'
   70 CONTINUE
c remove Exit(100)...make return statement
      WRITE(IOUT,1025)'FAILED'
      CALL EXIT(100)     
   80 CONTINUE
      WRITE(IOUT,1025)'NEW DIRECTION'
      CALL GETGM(N32,NATOMS,GM,IOUT)
      CALL GMWRIT(NATOMS,GM,IOUT,NSFIG,I75)
      CALL EXIT(50)
   75 CONTINUE
      WRITE(IOUT,1025)'CONTINUES'
      ICODE=1
      GO TO 81
   91 CONTINUE
C         NEWTON-RAPHSON
      DO 85 I=1,NDEGF
   85 GRAD(I)=RNACME(I)
      IF(NFIXDIR.NE.0)THEN
C        CORRECT GRADIENT FOR CONSTRAINTS
      CALL GETLM(N32,NFIXDIR,LAMBDA,IOUT,RESTART)
      CALL FREEZG(NFIXDIR,FIXDIR,VALDIR,NDEGF,GRAD,NATOMS,
     X            GM,KILLER,SCR,LAMBDA,KLASS,IOUT)
      ENDIF
      GRAD(NDEGF+NFIXDIR+1)=CRTERN
      GRAD(NDEGF+NFIXDIR+2)=ZERO
      CALL CHKNRP(NATOMS,NDEGF,N32,IOUT,CRTERN,GM,GRAD,MACME,
     1            LACME,SCALE,LAMBDA,IGO,MAXIT,CONV,ES1,
     2            NFIXDIR,NCGNT)
      RLAG=LAMBDA(NFIXDIR+1,1)
      RMAG=LAMBDA(NFIXDIR+2,1)
      IF(IGO)70,75,92
   92 CONTINUE
C         CALCULATION CONVERGED OR MAXIT REACHED
      WRITE(IOUT,1025)'HALTED'
      ICODE=100
   81 CONTINUE
C
C          CONSTRUCT HESSIAN
      IF(IHESS.EQ.0)GO TO 400
      SIMULT=1
      IF(KSEC.EQ.0)THEN
      CALL BLDHES(N32,NATOMS,NDEGF,GM,NACME,MACME,LACME,
     1            EGRADH,GGRADH,XGRADH,DISP,CRTERN,IGO,CDIF,
     2            LAMBDA,NFIXDIR,ES1,SAME)
      ELSE
      CALL BLDSECH(N32,I99,NATOMS,NATOMS*3,NTSYMNT,ORIGIN,GM,
     1             NACME,MACME,LACME,EGRADH,GGRADH,XGRADH,
     2             SYMTRN,DISP,CRTERN,IGO,CDIF,LAMBDA,NFIXDIR,ES1,SAME)
      CALL IDOF(NATOMS,NACME,SCR,KILLER)
      CALL IDOF(NATOMS,MACME,SCR,KILLER)
      CALL IDOF(NATOMS,LACME,SCR,KILLER)
      ENDIF
      IF(IGO.GT.0)WRITE(IOUT,1011)IGO
      IF(IGO) 410,420,430
  420 CONTINUE
      STOP ' ILLEGAL PATH '
C
  410 CONTINUE
      RLAG=LAMBDA(NFIXDIR+1,1)
      RMAG=LAMBDA(NFIXDIR+2,1)
C         ANALYZE W-MATRIX
      CALL ANALW(NDEGF,LAMBDA,RVECT,RMACME,RLACME,EGRADH,GGRADH,
     1          XGRADH,SCR,SCR1,SCR2,IOUT,N32,CDIF,MHESS,DISP,
     2          HESS,NROT,NFIXDIR,FIXDIR,KILLER,VALDIR,KLASS,
     3          NATOMS,GM,KZ,ORIGIN,SAME,KSEC,NOPS,SYMTRN,CHAR,
     4          NATOMS*3,NTSYM,NTSYMNT)
      IF(KSEC.EQ.0)THEN
      CALL BLDHESA(NDEGF,EGRADH,GGRADH,XGRADH,HESS,DISP,IOUT,
     1             RLAG,RMAG,MHESS,CDIF,SAME)
      ELSE
      CALL GETEGG(N32,NATOMS*3,EGRADH,GGRADH,XGRADH,IOUT,ICDF)
      CALL BLDHESE(NATOMS,NTSYMNT,NTSYM,NATOMS*3,NOPS,SYMTRN,CHAR,
     1             EGRADH,GGRADH,XGRADH,HESS,SCR,SCR1,ORIGIN,DISP,IOUT,
     2             RLAG,RMAG,CDIF,SAME,MHESS)
      CALL IDOF(NATOMS,EGRADH,SCR,KILLER)
      CALL IDOF(NATOMS,GGRADH,SCR,KILLER)
      CALL IDOF(NATOMS,XGRADH,SCR,KILLER)
      ENDIF
      LAMBDA(NFIXDIR+1,1)=RLAG
      LAMBDA(NFIXDIR+2,1)=RMAG
C         MAKE AN ON INTERNAL BASIS
C         BUILD THE ROTATION VECTOR
      LOOPON=0
      IF(NROT.EQ.0)GO TO 415
      CALL GETIGM(N32,NATOMS,GM,IOUT)
      CALL MROTV(NATOMS,ORIGIN,RVECT,GM,KZ,KILLER)
      CALL CHKR(NDEGF,EGRADH,RVECT,IOUT,NROT,TOLROT,1)
      CALL CHKR(NDEGF,GGRADH,RVECT,IOUT,NROT,TOLROT,2)
      IF(SAME)CALL CHKR(NDEGF,XGRADH,RVECT,IOUT,NROT,TOLROT,3)
  415 CONTINUE
      CALL MVIBV(NATOMS,NDEGFI,NDEGF,TVECT,VVECT,I76,KILLER)
C         O.N. INTERNAL MOTION 
      CALL IBASIS(NROT,RVECT,VVECT,NDEGF,BVECT,SCR,IOUT,IPFLG,IFAIL,0)
      IF(IFAIL.EQ.0)GO TO 416
      LOOPON=LOOPON+1
      IF(LOOPON.GT.NROT)GO TO 910
      GO TO 415
  416 CONTINUE
      CALL RHESS(NDEGF,NDEGFI,BVECT,HESS,SCR,SCR1,IOUT,IDBG,
     1          'HESSIAN ')
      ICODE=100
      CALL PUTHES(N32,NATOMS,NDEGFI,HESS,IOUT,NFIXDIR)
      NEWTON=1
C
C         MAKE AN ON INTERNAL BASIS COORDINATE
C         BUILD THE ROTATION VECTOR
      LOOPON=0
      CALL GETIGM(N32,NATOMS,GM,IOUT)
      CALL MROTV(NATOMS,ORIGIN,RVECT,GM,KZ,KILLER)
      CALL CHKR(NDEGF,EGRADH,RVECT,IOUT,NROT,TOLROT,1)
      CALL CHKR(NDEGF,GGRADH,RVECT,IOUT,NROT,TOLROT,2)
      IF(SAME)CALL CHKR(NDEGF,XGRADH,RVECT,IOUT,NROT,TOLROT,3)
  451 CONTINUE
      CALL MVIBV(NATOMS,NDEGFI,NDEGF,TVECT,VVECT,I76,KILLER)
C         O.N. INTERNAL MOTION BASIS
      CALL IBASIS(NROT,RVECT,VVECT,NDEGF,BVECT,SCR,IOUT,IPFLG,IFAIL,
     X            NFIXDIR+NCGNT)
      IF(IFAIL.EQ.0)GO TO 452
      LOOPON=LOOPON+1
      IF(LOOPON.GT.NROT)GO TO 910
      GO TO 451
  452 CONTINUE
      GO TO 400
C
  430 CONTINUE
C       NEXT DISPLACEMENT FOR HESSIAN
      print *, "ksec = ", ksec
      IF(KSEC.EQ.0)THEN
      DO 431 I=1,NDEGF
  431 SCR(I)=ZERO
      IF(IGO.LT.NDEGF+1)GO TO 432
      SCR(IGO-NDEGF)=-DISP(IGO-NDEGF)
      GO TO 433
  432 CONTINUE
      SCR(IGO)=DISP(IGO)
  433 CONTINUE
      CALL INTMO(NATOMS,GM,SCR,KILLER)
      print *, "Calling gmwrit 1"
      CALL GMWRIT(NATOMS,GM,IOUT,NSFIG,I75)
      print *, "Calling exit(50)"
      CALL EXIT(50)
      ELSE
C           SYMMETRY EQUIVALENT CENTER PATH
      IF(IGO.LT.NDEGF+1)THEN
      XDSP=DISP(IGO)
      ELSE
      XDSP=-DISP(IGO-NDEGF)
      ENDIF
      CALL JNTMO(NATOMS,GM,SYMTRN(1,IGO),XDSP)
      CALL GMWRIT(NATOMS,GM,IOUT,NSFIG,I75)
C
      CALL SECGM(NATOMS,GM,I76,IOUT,NATUN)
      CALL GMWRIT(NATUN,GM,IOUT,NSFIG,I77)
      CALL EXIT(50)
      ENDIF 
  400 CONTINUE
C
      IF(IDISP.LT.0)GO TO 300
C
  450 CONTINUE
      if(leq.ne.1)go to 460
      ndegff=ndegf+1
      call cnstrv(NFIXDIR,FIXDIR,VALDIR,NDEGF,hess,NATOMS,
     1                  GM,KILLER,KLASS,IOUT,rhs)
      CALL GETGM(N32,NATOMS,GM,IOUT)
      CALL MROTV(NATOMS,ORIGIN,hess(ndegff+1),GM,KZ,KILLER)
      do 453 i=1,ndegf
      hess(2*ndegff+i)=rmacme(i)
      hess(3*ndegff+i)=rlacme(i)
  453 continue
      ix=0
      iy=0
      XSCALE1=1.D0
      XSCALE2=1.D0
      IF((ACCEL.NE.0).AND.(KSCALE.LT.2))XSCALE2=SCALE
      IF((ACCEL.NE.0).AND.(KSCALE.LT.1))XSCALE1=SCALE
      hess(ndegff)=-rhs*XSCALE1
      hess(2*ndegff)=-0.d0
      hess(3*ndegff)=-crtern*XSCALE2
      hess(4*ndegff)=-0.d0
      ioff=ndegf**2
      do 454 i=1,ndegf
      disp(i)=0.d0
      do 455 j=1,ndegf
      ix=ix+1
      iy=iy+1
      hesss(ix)=hess(iy)
  455 continue
      iy=iy+1
      hesss(ioff+i)=hess(iy)
  454 continue
      write(iout,1061)(hess(i),i=1,ndegff*ndegf)
c        transpose the matrix
      do 458 i=2,ndegf
      do 458 j=1,i-1
      ij=(i-1)*ndegff+J
      ji=(j-1)*ndegff+i
      xx=hess(ij)
      hess(ij)=hess(ji)
      hess(ji)=xx
  458 continue
      call  GAUSS(hess,ndegf)
      x=0.d0
      do 456 i=1,ndegf
      x=x+hess(i)**2
  456 disp(i)=hess(i)
      x=dsqrt(x)
      if(x.lt.conv)icode=100
      WRITE(IOUT,1037)(DISP(I),I=1,NDEGF)
      CALL CHKJ(NDEGF,HESSS(ioff+1),DISP,hesss,IOUT)
      go to 370
  460 continue
      IF(RLAG.EQ.0.D0)THEN
      CALL RESTLAG(NDEGF,RNACME,RMACME,RLAG,IOUT)
C CLM: This print is for debugging
      PRINT *, "RESTLAG SUCCESSFUL"
      LAMBDA(NFIXDIR+1,1)=RLAG
      LAMBDA(NFIXDIR+2,1)=RMAG
      ENDIF
C
      IF(IPFLG.GT.2)WRITE(IOUT,1052)
      DO 250 I=1,NDEGF
      IF(IPFLG.GT.2)WRITE(IOUT,1053)I,RNACME(I),RMACME(I),
     X                              RLAG,RLACME(I),RMAG
  250 GRAD(I)=RNACME(I)+RLAG*RMACME(I)+RMAG*RLACME(I)
      IF(NFIXDIR.NE.0)THEN
C        CORRECT GRADIENT FOR CONSTRAINTS
      CALL GETLM(N32,NFIXDIR,LAMBDA,IOUT,RESTART)
      CALL FREEZG(NFIXDIR,FIXDIR,VALDIR,NDEGF,GRAD,NATOMS,
     X            GM,KILLER,SCR,LAMBDA,KLASS,IOUT)
      ENDIF
      IF(SAME)THEN
      GRAD(NDEGFL+NFIXDIR-1)=CRTERN
      GRAD(NDEGFL+NFIXDIR)=ZERO
      ELSE
      GRAD(NDEGFL+NFIXDIR)=CRTERN
      ENDIF
C         CONVERT GRADIENT TO INTERNAL COORDINATES
      IF(IHESS.NE.0)GO TO 253
      LOOPON=0
C     CLM: This print is for debugging
      IF(NROT.EQ.0)GO TO 252
      CALL GETGM(N32,NATOMS,GM,IOUT)
      CALL MROTV(NATOMS,ORIGIN,RVECT,GM,KZ,KILLER)
      CALL CHKR(NDEGF,RNACME,RVECT,IOUT,NROT,TOLROT,1)
      CALL CHKR(NDEGF,RMACME,RVECT,IOUT,NROT,TOLROT,2)
      IF(SAME)CALL CHKR(NDEGF,RLACME,RVECT,IOUT,NROT,TOLROT,3)
      CALL CHKR(NDEGF,GRAD,  RVECT,IOUT,NROT,TOLROT,4)
  252 CONTINUE
      CALL MVIBV(NATOMS,NDEGFI,NDEGF,TVECT,VVECT,I76,KILLER)
C CLM: This print is for debugging
      print *, " Calling ibasis... "
      CALL IBASIS(NROT,RVECT,VVECT,NDEGF,BVECT,SCR,IOUT,IPFLG,IFAIL,
     X            NFIXDIR+NCGNT)
      IF(IFAIL.EQ.0)GO TO 253
      LOOPON=LOOPON+1
      IF(LOOPON.GT.NROT)GO TO 910
      GO TO 252
  253 CONTINUE
      WRITE(IOUT,1051)(GRAD(I),I=1,NINT)
      CALL EBTC(SCR,BVECT,GRAD,NINTI,NINT,1)
      GNEW=ZERO
      GOLD=ZERO
      DO 255 I=1,NINTI
      GOLD=GOLD+GRAD(I)**2
      GRAD(I)=SCR(I)
      GNEW=GNEW+GRAD(I)**2
  255 CONTINUE
      DO 257 I=NINTI+1,NINT
  257 GOLD=GOLD+GRAD(I)**2
      GNEW=DSQRT(GNEW)
      GOLD=DSQRT(GOLD)
      WRITE(IOUT,1050)GOLD,GNEW
C
  300 CONTINUE
      IHTOP=NINT2-NINT
      WRITE(IOUT,1015)
      IF(IHESEQ1.EQ.0)THEN
      CALL GETHES(N32,NDEGFI,HESSS,IOUT)
      IF(IDBG.GT.0)THEN
      WRITE(IOUT,1049)' R-R IC ',' HESSIAN '
      CALL PRVOM(HESSS,NDEGFI,2,IOUT)
      ENDIF
      DO 332 I=1,NINT2
  332 HESS(I)=ZERO
      IF(NFIXDIR.GT.0)THEN
C        BUILD LAGRANGE MULTIPLIER CONTRIBUTIONS TO  HESSIAN
      CALL GETLM(N32,NFIXDIR,LAMBDA,IOUT,RESTART)
      CALL BLDCHES(NFIXDIR,FIXDIR,VALDIR,NDEGF,HESS,NATOMS,
     X            GM,KILLER,SCR,SCR,LAMBDA,KLASS,IOUT,ipflag(1))
      ENDIF
      ISTART=IHTOP
      IF(SAME)ISTART=IHTOP-(NINT-1)
      DO 334 I=1,NDEGF
  334 HESS(ISTART+I)=RMACME(I)
      IF(SAME)THEN
      DO 335 I=1,NDEGF
  335 HESS(IHTOP+I)=RLACME(I)
      ENDIF
      CALL RHESS(NINT,NINTI,BVECT,HESS,SCR,SCR1,IOUT,IDBG,
     1      'HESSIAN ')
C         COPY INTERNAL HESSIAN INTO CORRECT SLOTS
      NFXD1=NFIXDIR+1
      IF(SAME)NFXD1=NFIXDIR+2
      IX=0
      DO 333 I=1,NDEGFI
      IY=(I+NFXD1)*(I+NFXD1-1)/2+NFXD1
      DO 333 J=1,I
      IX=IX+1
      IJ=IY+J
  333 HESS(IJ)=HESS(IJ)+HESSS(IX)
C
      ELSE
C
      DO 304 I=1,NINT2
  304 HESS(I)=ZERO
      WRITE(IOUT,1038)
      IJ=0
      DO 303 I=1,NDEGF
      IJ=IJ+I
  303 HESS(IJ)=1.D0+DFLOAT(I)*1.0D-04
C
      IF(NFIXDIR.GT.0)THEN
C        FIX LAGRANGE MULTIPLIER PORTION OF HESSIAN
      CALL GETLM(N32,NFIXDIR,LAMBDA,IOUT,RESTART)
      CALL BLDCHES(NFIXDIR,FIXDIR,VALDIR,NDEGF,HESS,NATOMS,
     X            GM,KILLER,SCR,SCR,LAMBDA,KLASS,IOUT,ipflag(1))
      ENDIF
      ISTART=IHTOP
      IF(SAME)ISTART=IHTOP-(NINT-1)
      DO 305 I=1,NDEGF
  305 HESS(ISTART+I)=RMACME(I)
      IF(SAME)THEN
      DO 306 I=1,NDEGF
  306 HESS(IHTOP+I)=RLACME(I)
      ENDIF
      CALL RHESS(NINT,NINTI,BVECT,HESS,SCR,SCR1,IOUT,IDBG,
     1          'HESSIAN ')
C
      ENDIF
C
C         OVERWRITE SCALE
      REWIND INP
      IF(IERRI.EQ.0)READ(INP,NACINT)
      SCALE=SCLNAC
      IF(ACCEL.EQ.0)SCALE=1.D0
      GNORM=0.D0
      DO 311 I=1,NINTI
  311 GNORM=GNORM+GRAD(I)**2
      GNORM=DSQRT(GNORM)
      WRITE(IOUT,1020)SCALE,GNORM,(GRAD(I),I=1,NINTI)
C         SCALE AND CHANGE SIGN OF GRADIENTS
      DO 307 I=1,NINTI
  307 GRAD(I)=-GRAD(I)
      LOW=NFIXDIR+NCGNT+1
      IF(KSCALE.EQ.0)LOW=1
      IF(KSCALE.EQ.1)LOW=NFIXDIR+1
      DO 308 I=LOW,NINTI
  308 GRAD(I)=GRAD(I)*SCALE
C
      DO 301 I=1,NINTI2
  301 HESSS(I)=HESS(I)
      WRITE(IOUT,1049)' NR ',' HESSIAN '
      CALL PRVOM(HESS,NINTI,2,IOUT)
      CALL GIVENS (NINTI,NINTI,NINTI,HESSS,SCR,ROOT,VECT)
      IHI=0
      WRITE(IOUT,1018)
      DO 322 I=1,NINTI
      LOW=IHI+1
      IHI=IHI+NINTI
  322 WRITE(IOUT,1036)I,ROOT(I),(VECT(J),J=LOW,IHI)
      CALL CHKONB(NINTI,VECT,IOUT)
      DO 321 I=1,NINTI
      IF(DAMPV(I).EQ.0.D0)DAMPV(I)=DAMP
  321 CONTINUE
      WRITE(IOUT,1039)(DAMPV(I),I=1,NINTI)
c     this print is for debugging
      CALL SOLNR(NINTI,GRAD,VECT,ROOT,ZEIG,DAMPV,DISP,SCR,IOUT,
     X           SADDLE,NFIXDIR+NCGNT)
      CALL CHKI(NINTI,HESS,GRAD,DISP,SCR,IOUT)
      WRITE(IOUT,1032)(DISP(I),I=1,NINTI)
      CALL TBAK(NINT,NINTI,BVECT,DISP,SCR,0)
      WRITE(IOUT,1037)(DISP(I),I=1,NDEGF)
      DO 323 I=1,NFIXDIR+NCGNT
  323 LAMBDA(I,2)=DISP(NDEGF+I)+LAMBDA(I,1)
      WRITE(IOUT,1033)'OLD',(LAMBDA(I,1),I=1,NFIXDIR+NCGNT)
      WRITE(IOUT,1033)'NEW',(LAMBDA(I,2),I=1,NFIXDIR+NCGNT)
C
C        THE FOLLOWING CODE IS USED TO HELP DEAL WITH OSCILLATORY BEHAVIOUR
C        OF SAME SYMMETRY MEX
      IF(.NOT.SAME)GO TO 370
      WRITE(IOUT,1056)
      RLAG=LAMBDA(NFIXDIR+1,2)
      RMAG=LAMBDA(NFIXDIR+2,2)
      DO 360 I=1,NDEGF
  360 GRAD(I)=RNACME(I)+RLAG*RMACME(I)+RMAG*RLACME(I)
      IF(NFIXDIR.NE.0)THEN
C        CORRECT GRADIENT FOR CONSTRAINTS
      CALL FREEZG(NFIXDIR,FIXDIR,VALDIR,NDEGF,GRAD,NATOMS,
     X            GM,KILLER,SCR,LAMBDA(1,2),KLASS,IOUT)
      ENDIF
      IF(IPFLG.GT.2)THEN
      DO 361 I=1,NDEGF
  361 WRITE(IOUT,1053)I,RNACME(I),RMACME(I),RLAG,RLACME(I),RMAG
      ENDIF
      X=CRTERN**2
      DO 362 I=1,NDEGF+NFIXDIR
  362 X=X+GRAD(I)**2
      X=DSQRT(X)
      WRITE(IOUT,1054)X,GNORM,CONV,(GRAD(I),I=1,NDEGF+NFIXDIR)
      IF(X.GT.CONV)GO TO 370
      WRITE(IOUT,1055)
      ICODE=100
  370 CONTINUE
C
C          PUT GRADIENTS NEW LAGRANGE MULTIPLIER ON PROGRESS FILE
      CALL PUTGRD(N32,NATOMS,NDEGF,CRTERN,GM,NACME,MACME,LACME,
     1            LAMBDA(1,2),LAMBDA(1,1),IOUT,GNORM,ES1,
     2            NFIXDIR+NCGNT)
C          APPEND DIPOLE DATA
      IF(RESTART.EQ.0)THEN
      CALL PUTPROP(I70,IOUT,NPROP,N32,LACME,NDEGF,CHANGE) 
      ELSE
      CALL PUTPROP1(N32,DIPOLR,XGRADR,NDEGF)
      ENDIF
      CALL INTMO(NATOMS,GM,DISP,KILLER)
      CALL GMWRIT(NATOMS,GM,IOUT,NSFIG,I75)
      CALL DIST(NATOMS,GM,IOUT)
      CALL FANGL(NATOMS,GM,IOUT)
      if(natoms.gt.3)call FDANGL(NATOMS,GM,IOUT)
      if(maxnew.ne.0)call  rptcon(NFIXDIR,FIXDIR,VALDIR,NDEGF,
     1              NATOMS,GM,KILLER,KLASS,IOUT)
C
      CALL SECGM(NATOMS,GM,I76,IOUT,NATUN)
      CALL GMWRIT(NATUN,GM,IOUT,NSFIG,I77)
      call mkcp(killer,natoms,origin,n32,iout)
C
      print *, "icode = ", icode
      return
  900 CONTINUE
      WRITE(IOUT,1019)NCIU
      CALL EXIT(100)
  910 CONTINUE
      WRITE(IOUTU,1042)LOOPON
      CALL EXIT(100)
  920 WRITE(IOUTU,*)'BAD NAMELIST FILE=',I73
      CALL EXIT(100)
  921 WRITE(IOUTU,*)'BAD NAMELIST FILE=',INP
      CALL EXIT(100)
 1000 FORMAT(5X,'NO NAMELIST INPUT- DEFAULT MASSES USED ' )
 1001 FORMAT(/15X,'SEARCH ALGORITHM FOR MINIMUM ENERGY CROSSING ' )
 1002 FORMAT(/5X,' CURRENT GEOMETRY ',/(1X,I2,1X,3F15.8) )
 1003 FORMAT(/5X,'RAW ',A6,' GRAD',/7X,'X',14X,'Y',14X,'Z'/(1X,3F15.9))
 1010 FORMAT(' POINTERS TO GEOMETRY AND NACMES:'2I6 )
 1011 FORMAT('  HESSIAN CONSTRUCTION MODE: NEW DIRECTION=',I3)
 1012 FORMAT(' MASS COMBINATION:',5(' M(',I2,')=',F12.6),
     1       (/20X,5(' M(',I2,')=',F12.6) ))
 1015 FORMAT(5X,'DISPLACEMENT VECTOR GENERATED FROM H ! D> = -G> ' )
 1016 FORMAT('  CI-ENERGIES FROM UNIT ',I4,
     1       /5X,' E1=',E20.12,' E2=',E20.12,' E2-E1=',E20.12)
 1017 FORMAT(/5X,'DEL-E SCALED GRADIENTS',/7X,'X',14X,'Y',/(1X,2E16.9) )
 1018 FORMAT(/5X,'EIGENVALUES AND EIGENVECTORS OF HESSIAN' )
 1019 FORMAT(/5X,'CANT READ CI ENERGY FROM HEADER:UNIT=',I3)
 1020 FORMAT(/5X,'SCALE FACTOR FOR GRADIENT',E15.8,5X,' NORM=',E15.8,
     1       /5X,' UNSCALED',' GRADIENTS ',(/1X,5E15.8))
 1021 FORMAT(5X,'CURRENT GEOMETRY:',/(2X,I3,2X,3F12.6) )
 1022 FORMAT(5X,'**** CURRENT GEOMETRY TAKEN FROM UNIT:',I4 )
 1023 FORMAT(5X,'**** NO CURRENT GEOMETRY FILE ****' )
 1024 FORMAT(/5X,'**** NO GRADIENT POINTER ON TAPE=',I2,' ******' )
 1025 FORMAT(/5X,'ITERATION STATUS  ',A20)
 1026 FORMAT('  CI-ENERGIES FROM UNITS ',2I4,' E2-E1=',E20.12)
 1027 FORMAT(5X,'RESTART SEARCH FROM PROGRESS FILE=',I3)
 1029 FORMAT(5X,' NO INPUT FOUND ON INPUT UNIT=',I5)
 1030 FORMAT(5X,' NO INPUT AT ALL ON UNITS=',2I5,' ABORT' )
 1031 FORMAT(5X,'ATOM MAP=',10I3)
 1032 FORMAT(5X,'SOLUTION IN INTERNAL BASIS:',/,(2X,6F12.6))
 1033 FORMAT(5X,'LAGRANGE MULIPLIER:',A3,'=',30F12.6)
 1034 FORMAT(5X,'RESTART HESSIAN CONSTRUCT FROM PROGRESS FILE='
     1     ,I3,' IDISP=',I3)
 1035 FORMAT(5X,'STATE ENERGIES: ECI=',2E20.12,'  ENUC=',E20.12)
 1036 FORMAT(/,' ROOT=',I3,' EVAL=',E15.8,
     1       /,' EVECT=',6F12.6,/,(7X,6F12.6))
 1037 FORMAT(5X,'NUCLEAR DISPLACEMENTS:',/,(2X,6F12.6))
 1038 FORMAT(5X,'**** UNIT HESSIAN USED ****' )
 1039 FORMAT(5X,' DAMPING VECTOR:',/,(2X,6E12.5))
 1040 FORMAT(5X,' NORM OF GGRAD-OC=',E15.8,' IC=',E15.8)
 1041 FORMAT(5X,' COORDINATES IN ', I2,' DIMENSIONS SEARCHED')
 1042 FORMAT(5X,'SCHMO FAILURE IN IBASIS-MVIBV',I3)
 1043 FORMAT(5X,'DIMENSIONS:NDEGF=',I3,' NROT=',I3,' NFIX=',I3,' IVAR=',
     1       I3,' IVAR+ROT=',I3)
 1044 FORMAT(15X,'ACTUAL CROSSING MINIMIZATION')
 1045 FORMAT(5X,'CONSTRAINT CONDITION ',I2,' R(',I2,',',I2,
     X       ')=     ',F12.6)
 1046 FORMAT(5X,'CONSTRAINT CONDITION ',I2,' ANG(',I2,',',I2,
     X       ',',I2,')=',F12.6)
 1047 FORMAT(5X,'CONSTRAINT CONDITION ',I2,' R(',I2,',',I2,
     X       ')=         R(',I2,',',I2,')')
 1048 FORMAT(5X,'CONSTRAINT CONDITION ',I2,' ANG(',I2,',',I2,
     X       ',',I2,')=ANG(',I2,',',I2,',',I2,')')
 1049 FORMAT(3X,A8,1X,A8)
 1050 FORMAT(5X,'NORM OF RHS: ORIG BASIS ',E15.8,' INT BASIS ',E15.8)
 1051 FORMAT(5X,'GRADIENT: ORIGINAL BASIS',(/5X,5E15.8))
 1052 FORMAT(15X,'GRADIENT DECOMPOSITION')
 1053 FORMAT(3X,'I=',I2,' E=',F11.5,' G=',F11.5,' L1=',F11.5,
     X                  ' H=',F11.5,' L2=',F11.5)
 1054 FORMAT(/5X,'CURRENT GRADIENT FROM NEW LM - NORM=',
     X       E15.8,' OLD NORM=',E15.8,' CONV=',E15.8,(/5X,6E12.5))
 1055 FORMAT(/5X,'ITERATION STATUS  HALTED - NEW GRADIENT CRITERION' )
 1056 FORMAT(/15X,'RE-ANALYZE GRADIENT')
 1057 FORMAT(/15X,'SCALE RHS METHOD: ',A20,' GRADIENT')
 1058 FORMAT(/15X,'CANONICAL COORDINATES',/13X,'r',12X,'R',12X,'GMA',
     X       /3X,'g   ',3F12.6,/3X,'h   ',3F12.6,/3X,'seam',3F12.6)
 1059 FORMAT(/15X,'JACOBI ANALYSIS OF G-H PLANE ')
 1060 FORMAT(/15X,' ANALYSIS OF G-H PLANE ',/5X,
     X      ' NORMH=  ',F10.5, ' NORMG/2= ',F10.5,' ANG=   ',F9.3,
     X       /5X,' DEL=  ', F7.3,5X,' gxh/2=',E15.6, /5X,
     X      ' NORM-AV=', F10.5,' PERCENT-GH= ',F8.3,
     X      ' C-AVG-G=',F10.5,' C-AVG-H=',F10.5,/)
 1061  format(3x,'Linear equation matrix:',/(5x,4f10.5,3x,f12.7))
 1062  FORMAT(5X,'CONSTRAINT CONDITION ',I2,' ANG(',I2,',',I2,
     X       ',',I2,',',i2,' )= ',f12.6)
      END
      SUBROUTINE ANALW(NDEGF,LAMBDA,RVECT,RMACME,RLACME,EGRADH,
     1           GGRADH,XGRADH,SCR,BF,SCR1,IOUT,N32,ICDF,MHESS,
     2           DISP,HESS,NROT,NFIXDIR,FIXDIR,KILLER,VALDIR,
     3           KLASS,NATOMS,GM,KZ,ORIGIN,SAME,KSEC,NOPS,SYMTRN,
     4           CHAR,N3,NTSYM,NTSYMNT)
      IMPLICIT REAL*8(A-H,O-Z)
      DOUBLE PRECISION LAMBDA(2)
      LOGICAL SAME
      INTEGER FIXDIR(6,2),ORIGIN
      DIMENSION RVECT(NDEGF,1),EGRADH(NDEGF,1),GGRADH(NDEGF,1),
     1          SCR(2),RMACME(NDEGF),DISP(2),HESS(2),GM(2),
     2          BF(2),SCR1(2),KILLER(2),KLASS(2),VALDIR(2),
     3          XGRADH(NDEGF,1),RLACME(NDEGF),SYMTRN(N3,N3),
     4          CHAR(N3,2)
C
      WRITE(IOUT,1004)
      NELIM=NROT+NFIXDIR+1 
      IF(SAME)NELIM=NELIM+1
      RLAG=LAMBDA(NFIXDIR+1)
      RMAG=LAMBDA(NFIXDIR+2)
      IF(KSEC.EQ.0)THEN
      CALL GETEGG(N32,NDEGF,EGRADH,GGRADH,XGRADH,IOUT,ICDF)
      CALL BLDHESA(NDEGF,EGRADH,GGRADH,XGRADH,HESS,DISP,
     1             IOUT,RLAG,RMAG,MHESS,CDIF,SAME)
      ELSE
      CALL GETEGG(N32,N3,EGRADH,GGRADH,XGRADH,IOUT,ICDF)
      CALL BLDHESE(NATOMS,NTSYMNT,NTSYM,NATOMS*3,NOPS,SYMTRN,CHAR,
     1             EGRADH,GGRADH,XGRADH,HESS,SCR,SCR1,ORIGIN,DISP,IOUT,
     2             RLAG,RMAG,CDIF,SAME,MHESS)
      CALL IDOF(NATOMS,EGRADH,SCR,KILLER)
      CALL IDOF(NATOMS,GGRADH,SCR,KILLER)
      CALL IDOF(NATOMS,XGRADH,SCR,KILLER)
      ENDIF
      LAMBDA(NFIXDIR+1)=RLAG
      LAMBDA(NFIXDIR+2)=RMAG
      IF(NROT.EQ.0)GO TO 415
      CALL GETIGM(N32,NATOMS,GM,IOUT)
      CALL MROTV(NATOMS,ORIGIN,RVECT,GM,KZ,KILLER)
  415 CONTINUE
      DO 10 I=1,NDEGF
   10 RVECT(I,NROT+1)=RMACME(I)
      IF(SAME)THEN
      DO 11 I=1,NDEGF
   11 RVECT(I,NROT+2)=RLACME(I)
      ENDIF
C
      IF(NFIXDIR.NE.0)THEN
      CALL GETLM(N32,NFIXDIR,LAMBDA,IOUT,0)
      CALL BLDCHES(NFIXDIR,FIXDIR,VALDIR,NDEGF,HESS,NATOMS,
     X            GM,KILLER,SCR,SCR,LAMBDA,KLASS,IOUT)
      IOFF=NDEGF*(NDEGF+1)/2
      JOFF=1
      IF(SAME)JOFF=2
      DO 30 I=1,NFIXDIR
      DO 35 K=1,NDEGF
   35 RVECT(K,NROT+JOFF+I)=HESS(IOFF+K)
      IOFF=IOFF+NDEGF+I
   30 CONTINUE
      ENDIF
C
      CALL RBASIS(NELIM,RVECT,NDEGF,SCR,SCR1,BF,IOUT,1)
      NDEGR=NDEGF-NELIM
      WRITE(IOUT,1002)NDEGF,NDEGR
      CALL CHKR(NDEGF,EGRADH,BF,IOUT,NDEGR,10.D0,1)
      CALL CHKR(NDEGF,GGRADH,BF,IOUT,NDEGR,10.D0,2)
      IF(SAME)CALL CHKR(NDEGF,XGRADH,BF,IOUT,NDEGR,10.D0,3)
      CALL RHESS(NDEGF,NDEGR,BF,HESS,SCR,SCR1,IOUT,1,
     1          'W-MATRIX')
      IOFF=NDEGF*NDEGR
      CALL GIVENS (NDEGR,NDEGR,NDEGR,HESS,SCR,BF(IOFF+1),SCR1)
      IHI=0
      WRITE(IOUT,1001)
      DO 20 I=1,NDEGR
      LOW=IHI+1
      IHI=IHI+NDEGR
      WRITE(IOUT,1000)I,BF(I+IOFF),(SCR1(J),J=LOW,IHI)
      CALL TBAK(NDEGF,NDEGR,BF,SCR1(LOW),SCR,1)
      WRITE(IOUT,1003)(SCR(J),J=1,NDEGF)
   20 CONTINUE
      RETURN
 1000 FORMAT(/,' ROOT=',I3,' EVAL=',E15.8,
     1       /,' EVECT=',6F12.6,/,(7X,6F12.6))
 1001 FORMAT(/5X,'EIGENVALUES AND EIGENVECTORS OF W-MATRIX' )
 1002 FORMAT(/15X,'ANALYSIS OF W MATRIX',/5X,'NDEGF=',I3,' NDEGFS=',
     1        I3,/5X,'GRADIENTS IN SEAM SUBSPACE FOLLOW: ',
     2        '1=EGRAD 2=GGRAD 3=XGRAD',/)
 1003 FORMAT(' EV-INT',6F12.6,/,(7X,6F12.6))
 1004 FORMAT(/15X,'CONSTRUCT  W MATRIX',/)
      END
      SUBROUTINE ANGL(GM,II,JJ,KK,ANG,CA,R13,R12,R23)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION GM(3,2)
      DATA ZERO/0.D0/,TWO/2.D0/,PI/3.141592654D0/
C         LAW OF COSINES
      R13=ZERO
      R12=ZERO
      R23=ZERO
      DO 10 K=1,3
      R12=R12+(GM(K,II)-GM(K,JJ))**2
      R23=R23+(GM(K,JJ)-GM(K,KK))**2
      R13=R13+(GM(K,II)-GM(K,KK))**2
   10 CONTINUE
      CA=-(R13-R12-R23)/(TWO*DSQRT(R12*R23))
      ANG=DACOS(CA)*180.D0/PI
      R12=DSQRT(R12)
      R23=DSQRT(R23)
      R13=DSQRT(R13)
      RETURN
      END
      SUBROUTINE BLDHES(IPRU,NATOMS,NDEGF,GMP,EGRADPP,GGRADPP,
     1                  XGRADPP,EGRADP,GGRADP,XGRADP,DISPP,
     2                  EDIFP,IGO,ICDF,LAMBDAC,NFIXDIR,ES1,SAME)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL SAME
      REAL*8 LAMBDA(30),LAMBDAC(2),LAMBD0(2),LAMBDAR(30)
      DIMENSION EGRADPP(NDEGF),GGRADPP(NDEGF),EGRADP(NDEGF,2),
     1          GGRADP(NDEGF,2),HESSP(2),DISPP(NDEGF),GMP(2),
     2          XGRADP(NDEGF,2),XGRADPP(NDEGF)
      DIMENSION EGRADH(930),GGRADH(930),GM(30),EGRAD(30),GGRAD(30),
     1          XGRADH(930),HESS(465),DISP(30),XGRAD(30)
      NAMELIST/PROGRSH/GGRADH,EGRADH,IDISP,GM,EDIF,DISP,LAMBD0,EREF
     1                ,XGRADH
      NAMELIST/PROGRS/EGRAD,GGRAD,XGRAD,HESS,NEWDIR,SCALE,
     1                GM,EDIF,LAMBDA,RHSN,ESTATE1,LAMBDAR
      NAMELIST/STATUS/ISTATUS
      DATA SCALE/1.D0/,NEWDIR/1/,ISTATUS/-1/
C
      MAXDO=NDEGF+1 
      IF(ICDF.GT.0)MAXDO=NDEGF*2+1
C
      REWIND IPRU
      READ(IPRU,STATUS,END=5)
    5 CONTINUE
      REWIND IPRU
      READ(IPRU,PROGRSH,END=95)
  100 CONTINUE
      IF(IDISP.NE.0)GO TO 36
      X=EDIF
      READ(IPRU,PROGRS,END=36)
      LAMBD0(1)=LAMBDA(NFIXDIR+1)
      LAMBD0(2)=LAMBDA(NFIXDIR+2)
      DO 38 I=1,30
   38 GM(I)=GMP(I)
      EDIF=X
   36 CONTINUE
      IDISP=IDISP+1
      DO 20 I=1,NDEGF
      IF(IDISP.LE.1)DISP(I)=DISPP(I)
      DISPP(I)=DISP(I)
      EGRADH(I+(IDISP-1)*NDEGF)=EGRADPP(I)
      XGRADH(I+(IDISP-1)*NDEGF)=XGRADPP(I)
   20 GGRADH(I+(IDISP-1)*NDEGF)=GGRADPP(I)
      IF(IDISP.GT.1)GO TO 30
      EDIF=EDIFP
      EREF=ES1
   30 CONTINUE
      REWIND IPRU
      WRITE(IPRU,STATUS)
      CALL WRTPGSH(IPRU,NATOMS,NDEGF,EGRADH,GGRADH,
     1             XGRADH,GM,DISPP,EDIF,IDISP,LAMBD0,EREF)
      IGO=IDISP
      DO 26 I=1,NATOMS*3
   26 GMP(I)=GM(I)
      IF(IDISP.LT.MAXDO)RETURN
      EDIFP=EDIF
      ES1=EREF
      IX=0
      DO 25 J=1,MAXDO
      IF(J.LE.NDEGF)THEN
      XGRADPP(J)=XGRADH(J)
      GGRADPP(J)=GGRADH(J)
      EGRADPP(J)=EGRADH(J)
      ENDIF
      DO 25 I=1,NDEGF
      IX=IX+1
      XGRADP(I,J)=XGRADH(IX)
      GGRADP(I,J)=GGRADH(IX)
      EGRADP(I,J)=EGRADH(IX)
   25 CONTINUE
      LAMBDAC(NFIXDIR+1)=LAMBD0(1)
      LAMBDAC(NFIXDIR+2)=LAMBD0(2)
      IGO=-1
      RETURN
   95 CONTINUE
      DO 27 I=1,NATOMS*3
   27 GM(I)=GMP(I)
      IDISP=0
      GO TO 100
  900 CONTINUE
      WRITE(IOUT,1000)
      CALL EXIT(100)
 1000 FORMAT(5X,'CANT FIND DD-HESSIAN DATA - ABORT' )
      END
      SUBROUTINE BLDCHES(NFIXDIR,FIXDIR,VALDIR,NDEGF,HESS,NATOMS,
     1                  GM,KILLER,SCR,ISCR,RLAMBDA,KLASS,IOUT,iprint)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL LIK,LJK,LKK,LIL,LJL,LKL
      INTEGER FIXDIR(6,2),THESS,ISCR(2),KLASS(2)
      DIMENSION KILLER(2),VALDIR(2),HESS(2),SCR(2),RLAMBDA(2)
      DIMENSION GM(3,2)
      COMMON/CHESS/NORR
      common/pscent/maxnew,psc(30,10),dpsc(30,10),ipsc(30,10)
      DATA TWO/2.D0/,ZERO/0.D0/,PI/3.141592654D0/
C
C        STATEMENT FUNCTION GRADIENTS
      C1(X12,X23,R12,R23)=-(X23/R23+CA*X12/R12)/R12
      C2(X12,X23,R12,R23)=(X23-X12)/(R12*R23)
     1   +CA*(X12/R12**2-(X23/R23**2))
      C3(X12,X23,R12,R23)=(X12/R12+CA*X23/R23)/R23
C
      rl=1.d0
      if(maxnew.eq.0)go to 99
      do 40 i=1,maxnew
      do 41 k=1,3
   41 gm(k,natoms+i)=0
      do 42 k=1,3
      do 42 j=1,natoms
   42 gm(k,natoms+i)=gm(k,natoms+i)+gm(k,j)*psc(j,i)
   40 continue
   99 continue
C
      IF(NFIXDIR.EQ.0)GO TO 900
      N3=NATOMS*3
      THESS=NDEGF*(NDEGF+1)/2
C
      if(iprint.lt.2)go to 52
      WRITE(IOUT,1027)'WITHOUT','PRIMITIVE',NDEGF
      IHI=0
      DO 51 I=1,NDEGF
      LOW=IHI+1
      IHI=IHI+I
      WRITE(IOUT,1028)I,(HESS(K),K=LOW,IHI)
   51 CONTINUE
   52 continue
C
C            R-LAMBDA BLOCK
      DO 699 KFIXDIR=1,NFIXDIR
      DO 101 I=1,N3
  101 SCR(I)=ZERO
      II=FIXDIR(1,KFIXDIR)
      JJ=FIXDIR(2,KFIXDIR)
      GO TO (400,600,400,600,399),KLASS(KFIXDIR)
      STOP ' ILLEGAL KLASS '
  400 CONTINUE
C         DISTANCE CONSTRAINT
      call concon0(natoms,gm,scr,ii,jj,rl,zz)
      IF(KLASS(KFIXDIR).EQ.1)GO TO 200
C         DISTANCE DIFFERNECE CONSTRAINT
      II=FIXDIR(3,KFIXDIR)
      JJ=FIXDIR(4,KFIXDIR)
      call concon0(natoms,gm,scr,ii,jj,-rl,zz)
      GO TO 200
C
 399  CONTINUE
C         dihedral ANGLE CONSTRAINT
      KK=FIXDIR(3,KFIXDIR)
      LL=FIXDIR(6,KFIXDIR)
      call concon2(NATOMS,GM,SCR,II,JJ,KK,LL,RL)
      GO TO 200
  600 CONTINUE
C         ANGLE CONSTRAINT
      KK=FIXDIR(3,KFIXDIR)
      CALL ANGL(GM,II,JJ,KK,ANG,CA,R13,R12,R23)
      call concon1(NATOMS,GM,SCR,II,JJ,KK,RL)
      IF(KLASS(KFIXDIR).EQ.2)GO TO 200
C         ANGLE DIFFERENCE CONSTRAINT
      II=FIXDIR(4,KFIXDIR)
      JJ=FIXDIR(5,KFIXDIR)
      KK=FIXDIR(6,KFIXDIR)
      CALL ANGL(GM,II,JJ,KK,ANG,CA,R13,R12,R23)
      call concon1(NATOMS,GM,SCR,II,JJ,KK,-RL)
C
  200 CONTINUE
      KOUNT=0
      DO 140 I=1,N3
      IF(KILLER(I).EQ.1)GO TO 140
      KOUNT=KOUNT+1
      HESS(THESS+KOUNT)=SCR(I)
  140 CONTINUE
      DO 141 I=KOUNT+1,KOUNT+KFIXDIR
      HESS(THESS+I)=ZERO
  141 CONTINUE
      THESS=THESS+KOUNT+KFIXDIR
  699 CONTINUE
C
      if(maxnew.ne.0)go to 901
C            R-R BLOCK
      KOUNT=0
      DO 250 I=1,N3
      ISCR(I)=0
      IF(KILLER(I).EQ.1)GO TO 250
      KOUNT=KOUNT+1
      ISCR(I)=KOUNT
  250 CONTINUE
C
      DO 899 KFIXDIR=1,NFIXDIR
      II=FIXDIR(1,KFIXDIR)
      JJ=FIXDIR(2,KFIXDIR)
      GO TO (700,800,700,800,899),KLASS(KFIXDIR)
      STOP ' ILLEGAL KLASS '
  700 CONTINUE
C         DISTANCE CONSTRAINT
      DO 260 K=1,3
      IIK=ISCR((II-1)*3+K)
      JJK=ISCR((JJ-1)*3+K)
      LIK=IIK.EQ.0
      LJK=JJK.EQ.0
      IJ=IIK*(IIK-1)/2 + JJK
      JI=JJK*(JJK-1)/2 + IIK
      IF(LIK.OR.LJK)GO TO 243
      IF(II.GT.JJ)HESS(IJ)=HESS(IJ)-RLAMBDA(KFIXDIR)*TWO
      IF(JJ.GE.II)HESS(JI)=HESS(JI)-RLAMBDA(KFIXDIR)*TWO
  243 CONTINUE
      IF(LIK)GO TO 244
      IJ=IIK*(IIK+1)/2
      HESS(IJ)=HESS(IJ)+RLAMBDA(KFIXDIR)*TWO
  244 CONTINUE
      IF(LJK)GO TO 260
      IJ=JJK*(JJK+1)/2 
      HESS(IJ)=HESS(IJ)+RLAMBDA(KFIXDIR)*TWO
  260 CONTINUE
      IF(KLASS(KFIXDIR).EQ.1)GO TO 899
C         DISTANCE DIFFERNECE CONSTRAINT
      II=FIXDIR(3,KFIXDIR)
      JJ=FIXDIR(4,KFIXDIR)
      DO 270 K=1,3
      IIK=ISCR((II-1)*3+K)
      JJK=ISCR((JJ-1)*3+K)
      LIK=IIK.EQ.0
      LJK=JJK.EQ.0
      IJ=IIK*(IIK-1)/2 + JJK
      JI=JJK*(JJK-1)/2 + IIK
      IF(LIK.OR.LJK)GO TO 283
      IF(II.GT.JJ)HESS(IJ)=HESS(IJ)+RLAMBDA(KFIXDIR)*TWO
      IF(JJ.GE.II)HESS(JI)=HESS(JI)+RLAMBDA(KFIXDIR)*TWO
  283 CONTINUE
      IF(LIK)GO TO 284
      IJ=IIK*(IIK+1)/2
      HESS(IJ)=HESS(IJ)-RLAMBDA(KFIXDIR)*TWO
  284 CONTINUE
      IF(LJK)GO TO 270
      IJ=JJK*(JJK+1)/2 
      HESS(IJ)=HESS(IJ)-RLAMBDA(KFIXDIR)*TWO
  270 CONTINUE
C
      GO TO 899
  800 CONTINUE
C         ANGLE CONSTRAINT
C        NEXT LINE PERMITS SKIPPING OF RR ANGLE CONTRIBUTION
      IF(NORR.EQ.1)GO TO 899
      KK=FIXDIR(3,KFIXDIR)
      CALL ANGL(GM,II,JJ,KK,ANG,CA,R13,R12,R23)
      DO 370 K=1,3
      IIK=ISCR((II-1)*3+K)
      JJK=ISCR((JJ-1)*3+K)
      KKK=ISCR((KK-1)*3+K)
      LIK=IIK.NE.0
      LJK=JJK.NE.0
      LKK=KKK.NE.0
      X12=GM(K,II)-GM(K,JJ)
      X23=GM(K,JJ)-GM(K,KK)
C         X(I)X(I) Y(I)Y(I) Z(I) Z(I) DERIVATIVES
C         X(J)X(J) Y(J)Y(J) Z(J) Z(J) DERIVATIVES
C         X(K)X(K) Y(K)Y(K) Z(K) Z(K) DERIVATIVES
      IF(LIK)THEN
      IJ=IIK*(IIK+1)/2
      ZZ=-C1(X12,X23,R12,R23)*TWO*X12/R12**2 - CA/R12**2  
     1    + CA*(X12/R12**2)**2
      HESS(IJ)=HESS(IJ)+ZZ*RLAMBDA(KFIXDIR)
      ENDIF
      IF(LKK)THEN
      IJ=KKK*(KKK+1)/2
      ZZ=-C3(X12,X23,R12,R23)*TWO*X23/R23**2 - CA/R23**2 
     1    + CA*(X23/R23**2)**2
      HESS(IJ)=HESS(IJ)+ZZ*RLAMBDA(KFIXDIR)
      ENDIF
      IF(LJK)THEN
      IJ=JJK*(JJK+1)/2
      ZZ=TWO/(R12*R23)+(X12-X23)*(-X12/R12**2+X23/R23**2)/(R23*R12)
     1  -C2(X12,X23,R12,R23)*(-X12/R12**2+X23/R23**2)+CA*
     2  (TWO*X12**2/R12**4+TWO*X23**2/R23**4-ONE/R12**2-ONE/R23**2)
      HESS(IJ)=HESS(IJ)+ZZ*RLAMBDA(KFIXDIR)
      ENDIF
C
C          X(I)X(J)  Y(I)Y(J)  Z(I)Z(J) DERIVATIVES
C          X(J)X(K)  Y(J)Y(K)  Z(J)Z(K) DERIVATIVES
C          X(I)X(K)  Y(I)Y(K)  Z(I)Z(K) DERIVATIVES
C
      IF(LIK.AND.LJK)THEN
      IF(IIK.GT.JJK)THEN
      IJ=IIK*(IIK-1)/2+JJK
      ELSE
      IJ=JJK*(JJK-1)/2+IIK
      ENDIF
      ZZ=-ONE/(R12*R23)-(-X12+X23)*X12/(R23*R12**3)-
     1  CA*(-ONE/R12**2-X12**2/R12**4)-C1(X12,X23,R12,R23)*(-X12/R12**2
     2  +X23/R23**2)
      HESS(IJ)=HESS(IJ)+ZZ*RLAMBDA(KFIXDIR)
      ENDIF
      IF(LJK.AND.LKK)THEN
      IF(JJK.GT.KKK)THEN
      IJ=JJK*(JJK-1)/2+KKK
      ELSE
      IJ=KKK*(KKK-1)/2+JJK
      ENDIF
      ZZ=-ONE/(R12*R23)+(-X12+X23)*X23/(R12*R23**3)-
     1  CA*(-ONE/R23**2+X23**2/R23**4)-C3(X12,X23,R12,R23)*(-X12/R12**2
     2  +X23/R23**2)
      HESS(IJ)=HESS(IJ)+ZZ*RLAMBDA(KFIXDIR)
      ENDIF
      IF(LIK.AND.LKK)THEN
      IF(IIK.GT.KKK)THEN
      IJ=IIK*(IIK-1)/2+KKK
      ELSE
      IJ=KKK*(KKK-1)/2+IIK
      ENDIF
      ZZ=ONE/(R12*R23)-X12**2/(R23*R12**3)
     1    +C1(X12,X23,R12,R23)*(X23/R23**2)
      HESS(IJ)=HESS(IJ)+ZZ*RLAMBDA(KFIXDIR)
      ENDIF
      DO 371 L=1,3
      IF(K.EQ.L)GO TO 371
      IIL=ISCR((II-1)*3+L)
      JJL=ISCR((JJ-1)*3+L)
      KKL=ISCR((KK-1)*3+L)
      LIL=IIL.NE.0
      LJL=JJL.NE.0
      LKL=KKL.NE.0
      X12L=GM(L,II)-GM(L,JJ)
      X23L=GM(L,JJ)-GM(L,KK)
C
C          X(I)Y(J)   X(J)Y(K)   X(I)Y(K)  DERIVATIVES
C
      IF(LJL.AND.LIK)THEN
C  K(1)  L(2)     X1 W2
      IF(JJL.GT.IIK)THEN
      IJ=JJL*(JJL-1)/2+IIK
      ELSE
      IJ=IIK*(IIK-1)/2+JJL
      ENDIF
      ZZ=-(-X12L+X23L)*X12/(R23*R12**3)+
     1    C1(X12,X23,R12,R23)*(-X12L/R12**2 + X23L/R23**2)-
     2    CA*TWO*X12*X12L/R12**4
      HESS(IJ)=HESS(IJ)+ZZ*RLAMBDA(KFIXDIR)
      ENDIF
      IF(LJL.AND.LKK)THEN
C  K(3)  L(2)  X3 W2 
      IF(JJL.GT.KKK)THEN
      IJ=JJL*(JJL-1)/2+KKK
      ELSE
      IJ=KKK*(KKK-1)/2+JJL
      ENDIF
      ZZ=-(-X12L+X23L)*X23/(R12*R23**3)+
     1    C3(X12,X23,R12,R23)*(-X12L/R12**2 + X23L/R23**2)+
     2    CA*TWO*X23*X23L/R23**4
      HESS(IJ)=HESS(IJ)+ZZ*RLAMBDA(KFIXDIR)
      ENDIF
      IF(LKL.AND.LIK)THEN
C   K(1)  L(3)   X1 W3
      IF(KKL.GT.IIK)THEN
      IJ=KKL*(KKL-1)/2+IIK
      ELSE
      IJ=IIK*(IIK-1)/2+KKL
      ENDIF
      ZZ=-X12*X12L/(R23*(R12**3))+C1(X12,X23,R12,R23)*(X23L/R23**2)
     1    -(X23/(R12*R23)-CA*X12/R12**2)*X23L/R23**2
      HESS(IJ)=HESS(IJ)+ZZ*RLAMBDA(KFIXDIR)
      ENDIF
      IF(L.GT.K)GO TO 371 
C
C          X(I)Y(I) X(J)Y(J) X(K)Y(K) DERIVATIVES
C
      IF(LIK.AND.LIL)THEN
C  K(1)  L(1)   W1  X1
      IJ=IIK*(IIK-1)/2+IIL
      ZZ=X23*X12L/(R23*R12**3)-C1(X12L,X23L,R12,R23)*X12/R12**2
     1   +TWO*CA*X12*X12L/R12**4
      HESS(IJ)=HESS(IJ)+ZZ*RLAMBDA(KFIXDIR)
      ENDIF
      IF(LJK.AND.LJL)THEN
C  K(2)  L(2)   W2  X2
      IJ=JJK*(JJK-1)/2+JJL
      ZZ=-(-X12+X23)*(X23L/R23**2-X12L/R12**2)/(R12*R23)+
     1    C2(X12L,X23L,R12,R23)*(X12/R12**2-X23/R23**2)+
     2    CA*TWO*(X12*X12L/R12**4+X23*X23L/R23**4)
      HESS(IJ)=HESS(IJ)+ZZ*RLAMBDA(KFIXDIR)      
      ENDIF
      IF(LKK.AND.LKL)THEN
C  K(3)  L(3)   W3  X3
      IJ=KKK*(KKK-1)/2+KKL
      ZZ=X12*X23L/(R12*R23**3)+C3(X12L,X23L,R12,R23)*X23/R23**2
     1   +TWO*CA*X23*X23L/R23**4
      HESS(IJ)=HESS(IJ)+ZZ*RLAMBDA(KFIXDIR)
      ENDIF
  371 CONTINUE
  370 CONTINUE
      IF(KLASS(KFIXDIR).EQ.2)GO TO 899
C         ANGLE DIFFERENCE CONSTRAINT
      II=FIXDIR(4,KFIXDIR)
      JJ=FIXDIR(5,KFIXDIR)
      KK=FIXDIR(6,KFIXDIR)
      CALL ANGL(GM,II,JJ,KK,ANG,CA,R13,R12,R23)
      DO 380 K=1,3
      IIK=ISCR((II-1)*3+K)
      JJK=ISCR((JJ-1)*3+K)
      KKK=ISCR((KK-1)*3+K)
      LIK=IIK.NE.0
      LJK=JJK.NE.0
      LKK=KKK.NE.0
      X12=GM(K,II)-GM(K,JJ)
      X23=GM(K,JJ)-GM(K,KK)
C         X(I)X(I) Y(I)Y(I) Z(I) Z(I) DERIVATIVES
C         X(J)X(J) Y(J)Y(J) Z(J) Z(J) DERIVATIVES
C         X(K)X(K) Y(K)Y(K) Z(K) Z(K) DERIVATIVES
      IF(LIK)THEN
      IJ=IIK*(IIK+1)/2
      ZZ=-C1(X12,X23,R12,R23)*TWO*X12/R12**2 - CA/R12**2  
     1    + CA*(X12/R12**2)**2
      HESS(IJ)=HESS(IJ)-ZZ*RLAMBDA(KFIXDIR)
      ENDIF
      IF(LKK)THEN
      IJ=KKK*(KKK+1)/2
      ZZ=-C3(X12,X23,R12,R23)*TWO*X23/R23**2 - CA/R23**2 
     1    + CA*(X23/R23**2)**2
      HESS(IJ)=HESS(IJ)-ZZ*RLAMBDA(KFIXDIR)
      ENDIF
      IF(LJK)THEN
      IJ=JJK*(JJK+1)/2
      ZZ=TWO/(R12*R23)+(X12-X23)*(-X12/R12**2+X23/R23**2)/(R23*R12)
     1  -C2(X12,X23,R12,R23)*(-X12/R12**2+X23/R23**2)+CA*
     2  (TWO*X12**2/R12**4+TWO*X23**2/R23**4-ONE/R12**2-ONE/R23**2)
      HESS(IJ)=HESS(IJ)-ZZ*RLAMBDA(KFIXDIR)
      ENDIF
C
C          X(I)X(J)  Y(I)Y(J)  Z(I)Z(J) DERIVATIVES
C          X(J)X(K)  Y(J)Y(K)  Z(J)Z(K) DERIVATIVES
C          X(I)X(K)  Y(I)Y(K)  Z(I)Z(K) DERIVATIVES
C
      IF(LIK.AND.LJK)THEN
      IF(IIK.GT.JJK)THEN
      IJ=IIK*(IIK-1)/2+JJK
      ELSE
      IJ=JJK*(JJK-1)/2+IIK
      ENDIF
      ZZ=-ONE/(R12*R23)-(-X12+X23)*X12/(R23*R12**3)-
     1  CA*(-ONE/R12**2-X12**2/R12**4)-C1(X12,X23,R12,R23)*(-X12/R12**2
     2  +X23/R23**2)
      HESS(IJ)=HESS(IJ)-ZZ*RLAMBDA(KFIXDIR)
      ENDIF
      IF(LJK.AND.LKK)THEN
      IF(JJK.GT.KKK)THEN
      IJ=JJK*(JJK-1)/2+KKK
      ELSE
      IJ=KKK*(KKK-1)/2+JJK
      ENDIF
      ZZ=-ONE/(R12*R23)+(-X12+X23)*X23/(R12*R23**3)-
     1  CA*(-ONE/R23**2+X23**2/R23**4)-C3(X12,X23,R12,R23)*(-X12/R12**2
     2  +X23/R23**2)
      HESS(IJ)=HESS(IJ)-ZZ*RLAMBDA(KFIXDIR)
      ENDIF
      IF(LIK.AND.LKK)THEN
      IF(IIK.GT.KKK)THEN
      IJ=IIK*(IIK-1)/2+KKK
      ELSE
      IJ=KKK*(KKK-1)/2+IIK
      ENDIF
      ZZ=ONE/(R12*R23)-X12**2/(R23*R12**3)
     1    +C1(X12,X23,R12,R23)*(X23/R23**2)
      HESS(IJ)=HESS(IJ)-ZZ*RLAMBDA(KFIXDIR)
      ENDIF
      DO 381 L=1,3
      IF(K.EQ.L)GO TO 381
      IIL=ISCR((II-1)*3+L)
      JJL=ISCR((JJ-1)*3+L)
      KKL=ISCR((KK-1)*3+L)
      LIL=IIL.NE.0
      LJL=JJL.NE.0
      LKL=KKL.NE.0
      X12L=GM(L,II)-GM(L,JJ)
      X23L=GM(L,JJ)-GM(L,KK)
C
C          X(I)Y(J)   X(J)Y(K)   X(I)Y(K)  DERIVATIVES
C
      IF(LJL.AND.LIK)THEN
C  K(1)  L(2)     X1 W2
      IF(JJL.GT.IIK)THEN
      IJ=JJL*(JJL-1)/2+IIK
      ELSE
      IJ=IIK*(IIK-1)/2+JJL
      ENDIF
      ZZ=-(-X12L+X23L)*X12/(R23*R12**3)-
     1    C1(X12,X23,R12,R23)*(-X12L/R12**2 + X23L/R23**2)-
     2    CA*TWO*X12*X12L/R12**4
      HESS(IJ)=HESS(IJ)-ZZ*RLAMBDA(KFIXDIR)
      ENDIF
      IF(LJL.AND.LKK)THEN
C  K(3)  L(2)  X3 W2 
      IF(JJL.GT.KKK)THEN
      IJ=JJL*(JJL-1)/2+KKK
      ELSE
      IJ=KKK*(KKK-1)/2+JJL
      ENDIF
      ZZ=(-X12L+X23L)*X23/(R12*R23**3)-
     1    C3(X12,X23,R12,R23)*(-X12L/R12**2 + X23L/R23**2)+
     2    CA*TWO*X23*X23L/R23**4
      HESS(IJ)=HESS(IJ)-ZZ*RLAMBDA(KFIXDIR)
      ENDIF
      IF(LKL.AND.LIK)THEN
C   K(1)  L(3)   X1 W3
      IF(KKL.GT.IIK)THEN
      IJ=KKL*(KKL-1)/2+IIK
      ELSE
      IJ=IIK*(IIK-1)/2+KKL
      ENDIF
      ZZ=-X12*X12L/(R23*(R12**3))+C1(X12,X23,R12,R23)*(X23L/R23**2)
     1    -(X23/(R12*R23)-CA*X12/R12**2)*X23L/R23**2
      HESS(IJ)=HESS(IJ)-ZZ*RLAMBDA(KFIXDIR)
      ENDIF
      IF(L.GT.K)GO TO 381 
C
C          X(I)Y(I) X(J)Y(J) X(K)Y(K) DERIVATIVES
C
      IF(LIK.AND.LIL)THEN
C  K(1)  L(1)   W1  X1
      IJ=IIK*(IIK-1)/2+IIL
      ZZ=X23*X12L/(R23*R12**3)-C1(X12L,X23L,R12,R23)*X12/R12**2
     1   +TWO*CA*X12*X12L/R12**4
      HESS(IJ)=HESS(IJ)-ZZ*RLAMBDA(KFIXDIR)
      ENDIF
      IF(LJK.AND.LJL)THEN
C  K(2)  L(2)   W2  X2
      IJ=JJK*(JJK-1)/2+JJL
      ZZ=-(-X12+X23)*(X23L/R23**2-X12L/R12**2)/(R12*R23)+
     1    C2(X12L,X23L,R12,R23)*(X12/R12**2-X23/R23**2)+
     2    CA*TWO*(X12*X12L/R12**4+X23*X23L/R23**4)
      HESS(IJ)=HESS(IJ)-ZZ*RLAMBDA(KFIXDIR)      
      ENDIF
      IF(LKK.AND.LKL)THEN
C  K(3)  L(3)   W3  X3
      IJ=KKK*(KKK-1)/2+KKL
      ZZ=X12*X23L/(R12*R23**3)+C3(X12L,X23L,R12,R23)*X23/R23**2
     1   +TWO*CA*X23*X23L/R23**4
      HESS(IJ)=HESS(IJ)-ZZ*RLAMBDA(KFIXDIR)
      ENDIF
  381 CONTINUE
  380 CONTINUE
  899 CONTINUE
C
  901 continue
      if(iprint.lt.1)go to 851
      L=NDEGF+NFIXDIR
      WRITE(IOUT,1027)'WITH','ORIGINAL ',L
      IHI=0
      DO 850 I=1,L
      LOW=IHI+1
      IHI=IHI+I
      WRITE(IOUT,1028)I,(HESS(K),K=LOW,IHI)
  850 CONTINUE
  851 continue
C
  900 CONTINUE
      IF(KOUNT.NE.NDEGF)GO TO 999
      RETURN
  999 CONTINUE
      WRITE(IOUT,1024)KOUNT,NDEGF
      STOP ' ERROR BLDCHES '
 1024 FORMAT(5X,'PRIMITIVE NUCLEAR DEGREE OF FREEDOM ERROR:LOCAL=',
     1      I3,' PASSED=',I3)
 1025 FORMAT(2X,I3,(8F12.6),(/5X,8F12.6))
 1026 FORMAT(5X,'R-LAMBDA HESSIAN BLOCK: PRIMITIVE BASIS DIM=',I3 )
 1027 FORMAT(5X,'HESSIAN ',A7,' CONSTRAINTS:',A9,
     1    ' BASIS DIMENSION=',I3)
 1028 FORMAT(5X,I3,7E15.8,(/8X,7E15.8))
      END
      SUBROUTINE BLDHESA(NDIM,EGRADH,GGRADH,XGRADH,HESS,DISP,
     1                   IOUT,LAMBDA,LAMBDAX,MHESS,ICDF,SAME)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL SAME
      REAL*8 LAMBDA,LAMBDAX
      CHARACTER*8 MOT(3)
C
      DIMENSION GGRADH(NDIM,2),EGRADH(NDIM,2),XGRADH(NDIM,2),
     X          HESS(2),DISP(NDIM)
      DATA MOT/' LOWER  ','AVERAGE ',' UPPER  '/
      DATA TWO/2.D0/
      WRITE(IOUT,1002)MOT(MHESS+2)
C
      MAXDO=NDIM+1
      IF(ICDF.NE.0)THEN
      MAXDO=2*NDIM+1
      WRITE(IOUT,1004)
      ENDIF
C
      WRITE(IOUT,1005)'ENERGY-GRAD'
      DO 50 I=1,MAXDO
      WRITE(IOUT,1001)I,(EGRADH(J,I),J=1,NDIM)
   50 CONTINUE
      WRITE(IOUT,1005)'EDIF  -GRAD'
      DO 51 I=1,MAXDO
      WRITE(IOUT,1001)I,(GGRADH(J,I),J=1,NDIM)
   51 CONTINUE
      IF(SAME)THEN
      WRITE(IOUT,1005)'H(IJ) GRAD'
      DO 52 I=1,MAXDO
      WRITE(IOUT,1001)I,(XGRADH(J,I),J=1,NDIM)
   52 CONTINUE
      ENDIF
      IF(LAMBDA.EQ.0.D0)CALL RESTLAG(NDIM,EGRADH,GGRADH,LAMBDA,IOUT)
      IJ=0
      DO 100 I=1,NDIM
      DO 100 J=1,I
      IJ=IJ+1
      IF(ICDF.GT.0)GO TO 250
      HESSIJ=EGRADH(J,I+1)-EGRADH(J,1)+
     1       LAMBDA*(GGRADH(J,I+1)-GGRADH(J,1))+
     1       LAMBDAX*(XGRADH(J,I+1)-XGRADH(J,1))
      HESSIJ=HESSIJ/DISP(I)
      HESSJI=EGRADH(I,J+1)-EGRADH(I,1)+
     1       LAMBDA*(GGRADH(I,J+1)-GGRADH(I,1))+
     1       LAMBDAX*(XGRADH(I,J+1)-XGRADH(I,1))
      HESSJI=HESSJI/DISP(J)
      HESS(IJ)=HESSIJ
      GO TO 275
  250 CONTINUE
      HESSIJ=EGRADH(J,I+1)-EGRADH(J,I+NDIM+1)+
     1       LAMBDA*(GGRADH(J,I+1)-GGRADH(J,I+NDIM+1))+
     1       LAMBDAX*(XGRADH(J,I+1)-XGRADH(J,I+NDIM+1))
      HESSIJ=HESSIJ/(TWO*DISP(I))
      HESSJI=EGRADH(I,J+1)-EGRADH(I,J+NDIM+1)+
     1       LAMBDA*(GGRADH(I,J+1)-GGRADH(I,J+NDIM+1))+
     1       LAMBDAX*(XGRADH(I,J+1)-XGRADH(I,J+NDIM+1))
      HESSJI=HESSJI/(TWO*DISP(J))
      HESS(IJ)=HESSIJ
  275 CONTINUE
      IF(MHESS.EQ.-1)HESS(IJ)=HESSJI
      IF(MHESS.EQ.0)HESS(IJ)=(HESSJI+HESSIJ)/2.D0
      WRITE(IOUT,1000)I,J,IJ,HESSIJ,HESSJI
  100 CONTINUE
C
      DO 150 I=1,NDIM
      IJ=IJ+1
      HESS(IJ)=GGRADH(I,1)
      WRITE(IOUT,1000)NDIM+1,I,IJ,HESS(IJ)
  150 CONTINUE
      IJ=IJ+1
      HESS(IJ)=0.D0
      WRITE(IOUT,1000)NDIM+1,NDIM+1,IJ,HESS(IJ)
C
      IF(SAME)THEN
      DO 160 I=1,NDIM
      IJ=IJ+1
      HESS(IJ)=XGRADH(I,1)
      WRITE(IOUT,1000)NDIM+1,I,IJ,HESS(IJ)
  160 CONTINUE
      IJ=IJ+1
      HESS(IJ)=0.D0
      IJ=IJ+1
      HESS(IJ)=0.D0
      WRITE(IOUT,1000)NDIM+2,NDIM+2,IJ,HESS(IJ)
      ENDIF
C
      RETURN
 1000 FORMAT(5X,'I=',I2,' J=',I2,' HESS(',I2,')=',2E15.8)
 1001 FORMAT(2X,I3,7E15.8,(/5X,7E15.8))
 1002 FORMAT(/5X,'BUILD HESSIAN FROM DIVIDED DIFFERENCE OF GRADIENTS',
     1       /15X,'METHOD=',A8)
 1003 FORMAT(5X,'I=NDIM+1',' J=',I2,' HESS(',I2,')=',E15.8)
 1004 FORMAT(5X,'*** CENTERED DIFFERENCE USED ***',/)
 1005 FORMAT(/5X,A12)
      END
      SUBROUTINE BLDHESE(NATOMS,NTSYMNT,NTSYM,N3,NOPS,SYMTRN,CHAR,
     1                  EGRADH,GGRADH,XGRADH,HESS,TEM,SCR,
     2                  IORIGIN,DISP,IOUT,LAMBDA,LAMBDAX,ICDF,
     3                  SAME,MHESS)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL SAME
      REAL*8 LAMBDA,LAMBDAX
C
      DIMENSION EGRADH(N3,2),GGRADH(N3,2),XGRADH(N3,2)
     X         ,HESS(N3,N3),DISP(2),SYMTRN(N3,N3)
     X         ,TEM(N3),CHAR(N3,2),SCR(2)
      DATA TWO/2.D0/
C
      DO 4 I=1,N3
      DO 5 J=1,N3
    5 HESS(J,I)=0.D0
    4 IF(I.GT.NTSYM)HESS(I,I)=40.D0
C
      MAXDO=NTSYMNT+1
      IF(ICDF.NE.0)THEN
      MAXDO=2*NTSYMNT+1
      WRITE(IOUT,1004)
      ENDIF
C
      WRITE(IOUT,1005)'ENERGY-GRAD'
      DO 50 I=1,MAXDO
      WRITE(IOUT,1001)I,(EGRADH(J,I),J=1,N3)
   50 CONTINUE
      WRITE(IOUT,1005)'EDIF  -GRAD'
      DO 51 I=1,MAXDO
      WRITE(IOUT,1001)I,(GGRADH(J,I),J=1,N3)
   51 CONTINUE
      IF(SAME)THEN
      WRITE(IOUT,1005)'H(IJ) GRAD'
      DO 52 I=1,MAXDO
      WRITE(IOUT,1001)I,(XGRADH(J,I),J=1,N3)
   52 CONTINUE
      ENDIF
C
      IF(LAMBDA.EQ.0.D0)CALL RESTLAG(N3,EGRADH,GGRADH,LAMBDA,IOUT)
      WRITE(IOUT,1009)LAMBDA,LAMBDAX
C
      DO 100 I=1,NTSYMNT
      DO 100 J=1,N3
      IF(ICDF.GT.0)GO TO 55
      HESSIJ=EGRADH(J,I+1)-EGRADH(J,1)+
     1       LAMBDA*(GGRADH(J,I+1)-GGRADH(J,1))
      IF(SAME)HESSIJ=HESSIJ+ LAMBDAX*(XGRADH(J,I+1)-XGRADH(J,1))
      HESSIJ=HESSIJ/DISP(I)
      HESS(J,I)=HESSIJ
      GO TO 75
   55 CONTINUE
      HESSIJ=EGRADH(J,I+1)-EGRADH(J,I+NTSYMNT+1)+
     1       LAMBDA*(GGRADH(J,I+1)-GGRADH(J,I+NTSYMNT+1))
      IF(SAME)HESSIJ=HESSIJ+LAMBDAX*(XGRADH(J,I+1)-
     1       XGRADH(J,I+NTSYMNT+1))
      HESSIJ=HESSIJ/(TWO*DISP(I))
      HESS(J,I)=HESSIJ
   75 CONTINUE
      WRITE(IOUT,1000)I,J,HESSIJ
  100 CONTINUE
C
C        TRANSLATIONAL INVARIANCE
      DO 200 I=1,N3
      DO 205 K=1,3
      X=0.D0
      DO 210 J=1,NTSYMNT
      X=HESS(I,J)*CHAR(J,NOPS+K)+X
  210 CONTINUE
      DO 206 L=NTSYMNT+1,NTSYM
      IF(DABS(CHAR(L,NOPS+K)).LT.1.0D-04)GO TO 206
      HESS(I,L)=-X/CHAR(L,NOPS+K)
      GO TO 205
  206 CONTINUE
      IF(DABS(X).LT.1.0D-03)GO TO 205
      WRITE(IOUT,1007)I,K,X
      STOP ' TI ERROR IN BLDHESE '
  205 CONTINUE
  200 CONTINUE
      WRITE(IOUT,1006)' TI ',N3,NTSYM
      CALL MATOUT(HESS,N3,NTSYM,N3,NTSYM,IOUT)
      WRITE(IOUT,1008)N3
      CALL MATOUT(SYMTRN,N3,N3,N3,N3,IOUT)
C
C        TRANSFORM TO SYMMETRY COORDINATES TOTALLY SYMMETRIC BLOCK ONLY
      DO 145 I=1,N3
  145 TEM(I)=0.D0
      DO 150 I=1,NTSYM
      DO 140 K=1,NTSYM
      X=0.D0
      DO 160 J=1,N3
  160 X=X+HESS(J,I)*SYMTRN(J,K)
  140 TEM(K)=X
      DO 141 K=1,N3
  141 HESS(K,I)=TEM(K)
  150 CONTINUE
      WRITE(IOUT,1006)' SC ',N3,NTSYM
      CALL MATOUT(HESS,N3,NTSYM,N3,NTSYM,IOUT)
      IF(MHESS.NE.0)GO TO 229
      WRITE(IOUT,1011)
      DO 225 I=1,NTSYM
      DO 225 J=1,I
      X=(HESS(I,J)+HESS(J,I))/2.D0
      HESS(I,J)=X
      HESS(J,I)=X
  225 CONTINUE
  229 CONTINUE
      IX=0
      DO 230 I=1,NTSYM
      DO 230 J=1,I
      IX=IX+1
      SCR(IX)=HESS(J,I)
  230 CONTINUE
      CALL GIVENS(NTSYM,NTSYM,NTSYM,SCR,SCR(IX+1),TEM,TEM(NTSYM+1))
      WRITE(IOUT,1010)(TEM(I),I=1,NTSYM)
C        TRANSFORM TO ATOMIC COORDINATES 
      DO 250 I=1,N3
      DO 240 K=1,N3
      X=0.D0
      DO 260 J=1,N3
  260 X=X+HESS(J,I)*SYMTRN(K,J)
  240 TEM(K)=X
      DO 241 K=1,N3
  241 HESS(K,I)=TEM(K)
  250 CONTINUE
      DO 350 I=1,N3
      DO 340 K=1,N3
      X=0.D0
      DO 360 J=1,N3
  360 X=X+HESS(I,J)*SYMTRN(K,J)
  340 TEM(K)=X
      DO 341 K=1,N3
  341 HESS(I,K)=TEM(K)
  350 CONTINUE
C
      CALL CMPRS(NATOMS,HESS,IORIGIN,IOUT)
      RETURN
 1000 FORMAT(5X,' HESS(I=',I2,',J=',I2,') = ',2E15.8)
 1001 FORMAT(2X,I3,7E15.8,(/5X,7E15.8))
 1002 FORMAT(/5X,'BUILD HESSIAN FROM DIVIDED DIFFERENCE OF GRADIENTS',
     1       /15X,'METHOD=',A8)
 1003 FORMAT(5X,'I=NDIM+1',' J=',I2,' HESS(',I2,')=',E15.8)
 1004 FORMAT(5X,'*** CENTERED DIFFERENCE USED ***',/)
 1005 FORMAT(/5X,A12)
 1006 FORMAT(/5X, A8,' HESSIAN FROM BLDHESE: DIMENSION =',2I3)
 1007 FORMAT(/5X,'TI FAILURE BLDHESE:ROW=',I3,' XYZ=',I3,
     1      ' ROW SUM=',F12.7)
 1008 FORMAT(/5X,'SYMTRN  FROM BLDHESE:DIMENSION=',I3)
 1009 FORMAT(/5X,'LAGRANGE MULTIPLIERS FOR Q BLOCK ',2F12.6)
 1010 FORMAT(/5X,'IN BLDHESE: EIGENVALUES OF TS BLOCK OF HESSIAN',
     1       /(6(1X,F12.6)))
 1011 FORMAT(/5X,'AVERAGE TS BLOCK METHOD USED',/)
      END
      SUBROUTINE BLDSECH(IPRU,I99,NATOMS,N3,NTSYMNT,IORIGIN,GMP,
     1           EGRADPP,GGRADPP,XGRADPP,EGRADP,GGRADP,XGRADP,
     2           SYMTRN,DISPP,EDIFP,IGO,ICDF,LAMBDAC,NFIXDIR,ES1,SAME)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL SAME
      REAL*8 LAMBDA(30),LAMBDAC(2),LAMBD0(2),LAMBDAR(30)
      DIMENSION EGRADPP(N3),EGRADP(N3,2),SYMTRN(N3,N3),
     1         HESSP(2),DISPP(2),GMP(2),GGRADPP(N3),GGRADP(N3,2),
     2         XGRADPP(N3),XGRADP(N3,2)
      DIMENSION EGRADH(930),GGRADH(930),GM(30),EGRAD(30),GGRAD(30),
     1          XGRADH(930),HESS(465),DISP(30),XGRAD(30)
      NAMELIST/PROGRSH/GGRADH,EGRADH,IDISP,GM,EDIF,DISP,LAMBD0,EREF
     1                ,XGRADH
      NAMELIST/PROGRS/EGRAD,GGRAD,XGRAD,HESS,NEWDIR,SCALE,
     1                GM,EDIF,LAMBDA,RHSN,ESTATE1,LAMBDAR
      NAMELIST/STATUS/ISTATUS
      DATA SCALE/1.D0/,NEWDIR/1/,ISTATUS/-1/
C
      MAXDO=NTSYMNT+1 
      IF(ICDF.GT.0)MAXDO=NTSYMNT*2+1
      REWIND I99
      READ(I99)EGRADPP
      READ(I99)GGRADPP
      IF(SAME)READ(I99)XGRADPP
C
      REWIND IPRU
      READ(IPRU,STATUS,END=5)
    5 CONTINUE
      REWIND IPRU
      READ(IPRU,PROGRSH,END=95)
      GO TO 100
   95 CONTINUE
      DO 27 I=1,N3
   27 GM(I)=GMP(I)
      IDISP=0
  100 CONTINUE
      IF(IDISP.NE.0)GO TO 36
      X=EDIF
      READ(IPRU,PROGRS,END=36)
      LAMBD0(1)=LAMBDA(NFIXDIR+1)
      LAMBD0(2)=LAMBDA(NFIXDIR+2)
      DO 38 I=1,30
   38 GM(I)=GMP(I)
      EDIF=X
   36 CONTINUE
      IDISP=IDISP+1
      DO 20 I=1,N3
      EGRADH(I+(IDISP-1)*N3)=EGRADPP(I)
      GGRADH(I+(IDISP-1)*N3)=GGRADPP(I)
      XGRADH(I+(IDISP-1)*N3)=XGRADPP(I)
   20 CONTINUE
      DO 21 I=1,N3-3
      IF(IDISP.LE.1)DISP(I)=DISPP(I)
      DISPP(I)=DISP(I)
   21 CONTINUE
      IF(IDISP.GT.1)GO TO 30
      EDIF=EDIFP
      EREF=ES1
   30 CONTINUE
      REWIND IPRU
      WRITE(IPRU,STATUS)
      CALL WRTPGSH(IPRU,NATOMS,N3,EGRADH,GGRADH,XGRADH
     1            ,GM,DISPP,EDIF,IDISP,LAMBD0,EREF)
      IGO=IDISP
      DO 26 I=1,NATOMS*3
   26 GMP(I)=GM(I)
      IF(IDISP.LT.MAXDO)RETURN
      EDIFP=EDIF
      ES1=EREF
      IX=0
      DO 25 J=1,MAXDO
      DO 25 I=1,N3
      IX=IX+1
      EGRADP(I,J)=EGRADH(IX)
      GGRADP(I,J)=GGRADH(IX)
      XGRADP(I,J)=XGRADH(IX)
   25 CONTINUE
      DO 37 J=1,N3
      EGRADPP(J)=EGRADH(J)
      GGRADPP(J)=GGRADH(J)
      XGRADPP(J)=XGRADH(J)
   37 CONTINUE
      LAMBDAC(NFIXDIR+1)=LAMBD0(1)
      LAMBDAC(NFIXDIR+2)=LAMBD0(2)
      IGO=-1
      RETURN
      END
      SUBROUTINE CHKI(NDIM,S,RHS,LHS,WMAT,IOUT)
      IMPLICIT REAL*8(A-H,O-Z)
      DOUBLE PRECISION LHS
      DIMENSION RHS(NDIM),LHS(NDIM),WMAT(NDIM,NDIM),S(2)
      IX=0
      DO 10 I=1,NDIM
      DO 10 J=1,I
      IX=IX+1
      WMAT(I,J)=S(IX)
      WMAT(J,I)=S(IX)
   10 CONTINUE
      DO 20 I=1,NDIM
      DP=DOT(WMAT(1,I),LHS,NDIM)
      IF(DABS(DP-RHS(I)).LT.1.0D-10)GO TO 20
      WRITE(IOUT,1000)I,DP,RHS(I)
   20 CONTINUE
      RETURN
 1000 FORMAT(5X,'CHKI FAILURE: ROW=',I2,' HX=',E15.8,' RHS=',E15.8)
      END
      SUBROUTINE CHKJ(NDIM,RHS,LHS,WMAT,IOUT)
      IMPLICIT REAL*8(A-H,O-Z)
      DOUBLE PRECISION LHS
      DIMENSION RHS(NDIM),LHS(NDIM),WMAT(NDIM,NDIM),S(2)
      DO 20 I=1,NDIM
      dp=0.
      do 21 j=1,ndim
   21 DP=DP+WMAT(i,J)*LHS(J)
      IF(DABS(DP-RHS(I)).LT.1.0D-10)GO TO 20
      WRITE(IOUT,1000)I,DP,RHS(I)
   20 CONTINUE
      RETURN
 1000 FORMAT(5X,'CHKI FAILURE: ROW=',I2,' HX=',E15.8,' RHS=',E15.8)
      END
      SUBROUTINE BCOM(MASS,B)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 MASS(3),B(3,3),MD,MT
      EQUIVALENCE(N1,NDIAT(1)),(N2,NDIAT(2)),(N3,NDIAT(3))
      COMMON/DIATAT/NDIAT(3)
      MD=MASS(N2)+MASS(N3)
C         INVERSE TRANSFORMATION FROM ATOM CENTERED COORDINATES
C         TO INTERNAL COORDINATES RELATIVE TO
C         FRAME AT COM OF TRIATOM
      MT=MD+MASS(N1)
      B(N1,1)=MD/MT
      B(N2,1)=-MASS(N1)/MT
      B(N3,1)=-MASS(N1)/MT
      B(N1,2)=0.D0
      B(N2,2)=-MASS(N3)/MD
      B(N3,2)= MASS(N2)/MD
      B(N1,3)=1.D0
      B(N2,3)=1.D0
      B(N3,3)=1.D0
      RETURN
      END
      SUBROUTINE CHKONB(NDIM,WMAT,IOUT)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION WMAT(NDIM,NDIM)
      DATA ZERO,ONE/0.D0,1.D0/
      DO 10 I=1,NDIM
      TEST=ZERO
      DO 10 J=1,I
      IF(J.EQ.I)TEST=ONE
      DP=DOT(WMAT(1,I),WMAT(1,J),NDIM)
      IF(DABS(TEST-DP).LT.1.0D-08)GO TO 10
      WRITE(IOUT,1000)I,J,DP
   10 CONTINUE
      RETURN
 1000 FORMAT(5X,'CHKONB FAILURE: (',I2,' ! ',I2,')=',F15.8)
      END
      SUBROUTINE CHKNRP(NATOMS,NDEGF,IPRU,IOUT,EDIFC,GMC,EGRADC,
     1             GGRADC,XGRADC,SCALEC,LAMBDAC,
     2             IPATH,MAXIT,CONV,ES1,NFIXDIR,NCGNT)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 LAMBDA(30),LAMBDAC(2),LAMBD0(2),LAMBDAR(30)
      DIMENSION GMC(2),GGRADC(2),HESSC(2),EGRADC(2),XGRADC(2)
      DIMENSION GM(30),GGRAD(30),HESS(465),EGRAD(30),XGRAD(30)
      DIMENSION GGRADH(930),EGRADH(930),XGRADH(930),DISP(30)
      DATA RHSN/100.D0/,ESTATE/0.0D0/
      NAMELIST/STATUS/ISTATUS
      NAMELIST/PROGRSH/GGRADH,EGRADH,XGRADH,IDISP,GM,EDIF,DISP,
     1                 LAMBD0,EREF
      NAMELIST/PROGRS/GM,EDIF,NEWDIR,EGRAD,GGRAD,XGRAD,
     1                SCALE,HESS,LAMBDA,RHSN,ESTATE1,LAMBDAR
C
      IPATH=0
      SCALE=SCALEC
      NEWDIR=-1
      REWIND IPRU
      READ(IPRU,PROGRS,END=150)
      GO TO 175
  150 CONTINUE
      CALL RESTLAG(NDEGF,EGRADC,GGRADC,RLAG,IOUT)
      LAMBDA(NFIXDIR+1)=RLAG
      LAMBDA(NFIXDIR+2)=0.D0
      NEWDIR=0
      IPATH=0
      GO TO 201
  175 CONTINUE
C        CRITERIA ARE: (1) ENERGY DIFFERNECE WHICH CAN BE NEGATIVE
C                      (2) ESTATE1
C                      (3) NORM OF RHS
      IF(DABS(EDIFC).LT.DABS(EDIF))GO TO 200
C        MAGNITUDE OF ENERGY DIFFERENCE INCREASED
      IPATH=-1
  200 CONTINUE
      IF(ES1.LT.ESTATE1)IPATH=0
  201 CONTINUE
      LAMBDAC(NFIXDIR+1)=LAMBDA(NFIXDIR+1)
      LAMBDAC(NFIXDIR+2)=LAMBDA(NFIXDIR+2)
      SCALEC=SCALE
      GNORM=EDIFC**2
      DO 202 I=1,NDEGF
      GNORM=GNORM+(EGRADC(I)+LAMBDA(NFIXDIR+1)*GGRADC(I)+
     X             LAMBDA(NFIXDIR+2)*XGRADC(I))**2
  202 CONTINUE
      DO 207 I=NDEGF+1,NDEGF+1+NFIXDIR
  207 GNORM=GNORM+EGRADC(I)**2
      GNORM=DSQRT(GNORM)
      IF((IPATH.NE.-1) .OR. (GNORM.LT.RHSN))GO TO 203
      WRITE(IOUT,1004)'DIVERGED',NEWDIR,EDIFC,EDIF,GNORM,RHSN,
     X                 ES1,ESTATE1
      IF(NFIXDIR.GT.0)WRITE(IOUT,1007)(LAMBDA(I),I=1,NFIXDIR)
      IPATH=-1
      RETURN
  203 CONTINUE
      NEWDIR=NEWDIR+1
      WRITE(IOUT,1004)'CONTINUE',NEWDIR,EDIFC,EDIF,GNORM,RHSN,
     X                 ES1,ESTATE1
      REWIND IPRU
      READ(IPRU,PROGRSH,END=204)
      CALL WRTPGS(IPRU,NATOMS,NDEGF,GMC,EDIFC,NEWDIR,
     1    EGRAD,GGRAD,XGRAD,SCALE,HESS,LAMBDA,LAMBDAR,GNORM,
     2    ES1,NFIXDIR+2)
      IPATH=0
      GO TO 205
  204 CONTINUE
      REWIND IPRU
      READ(IPRU,STATUS)
      CALL WRTPGS(IPRU,NATOMS,NDEGF,GMC,EDIFC,NEWDIR,
     1    EGRADC,GGRADC,XGRADC,SCALE,HESS,LAMBDA,LAMBDAR,GNORM,
     2    ES1,NFIXDIR+2)
      IPATH=0
  205 CONTINUE
      IF(GNORM.GT.CONV)GO TO 206
      WRITE(IOUT,1003)'   ',NEWDIR,GNORM,EDIFC,CONV,MAXIT
      WRITE(IOUT,1007)(LAMBDA(I),I=1,NFIXDIR+NCGNT)
      IPATH=1
      RETURN
  206 CONTINUE
C      IF(NEWDIR.LT.MAXIT)RETURN
      IPATH=-1
      IF(NEWDIR.LE.MAXIT)IPATH=0
      WRITE(IOUT,1003)'NOT',NEWDIR,GNORM,EDIFC,CONV,MAXIT
      WRITE(IOUT,1007)(LAMBDA(I),I=1,NFIXDIR+NCGNT)
C      IPATH=1
      RETURN
C
 1000 FORMAT(5X,'NO PROGRESS FILE-INITIALZING FILE=',I3)  
 1001 FORMAT(5X,'NO PROGRESS FILE- ABORT' )
 1002 FORMAT(5X,'NO PROGRESS FILE-FOR HESSIAN RETRIVAL' )
 1003 FORMAT(5X,'NEWTON RAPHSON ITERATIONS ',A3,' CONVERGED',
     1          ' AFTER ',I3,' ITERATIONS ',
     2      /15X,'G-NORM=',E15.8,' VALUE=',E15.8,
     3      /15X,'CONV  =',E15.8,' MAXIT=',I3 )
 1004 FORMAT(5X,'NEWTON RAPHSON ITERATIONS ',A8, ': ',/5X,
     1       'ITER=',I4,4X,' EDIF: CURRENT=', E15.8,' LAST=',E15.8,
     2       /15X,'NORM RHS: CURRENT=', E15.8,' LAST=',E15.8,
     3       /15X,'E(1)    : CURRENT=', E15.8,' LAST=',E15.8)
 1005 FORMAT(5X,'NORM OF RHS=',E15.8)
 1007 FORMAT(15X,'CURRENT LAGRANGE MULTIPLIERS:',6F12.6,(/27X,6F12.6))
      END
      SUBROUTINE CHKR(N,A,V,IOUT,NV,TOL,JFAIL)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(2),V(N,NV)
C
      DP1=DSQRT(DOT(A,A,N))
      DO 100 I=1,NV
      IFAIL=I
      DP2=DSQRT(DOT(V(1,I),V(1,I),N))
      DP=DOT(A,V(1,I),N)
      WRITE(IOUT,1001)JFAIL,I,DP,DP1,DP2
      IF(DP.GT.TOL)GO TO 900
  100 CONTINUE
      RETURN
  900 CONTINUE
      WRITE(IOUT,1002)
      WRITE(IOUT,1000)JFAIL,IFAIL,TOL,DP,DP1,DP2
      WRITE(IOUT,1003)'ROT ',(V(K,IFAIL),K=1,N)
      WRITE(IOUT,1003)'GRAD',(A(K),K=1,N)
 1000 FORMAT(5X,'ROTATION VECTOR FAILURE:  TEST=',I1,
     1          ' DIREC=',I1,' TOLERANCE=',E15.8,/5X,
     2  ' DP=',E15.8,' NORMG=', E15.8,' NORMR=',E15.8)
 1001 FORMAT(5X,'TEST=',I2,' DIRECTION=',I2,' DP=',E15.8,' NORMG=',
     X       E15.8,' NORMR=',E15.8)
 1002 FORMAT(15X,'BAD ROTATION VECTOR OR GRADIENT:')
 1003 FORMAT(2X,A4,/,(5F12.6))
      CALL EXIT(100)
      END
      SUBROUTINE CMPRS(NATOMS,HESS,II,IOUT)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION HESS(2) 
      N3=NATOMS*3
      WRITE(IOUT,1000)'ORGINAL',N3
      CALL MATOUT(HESS,N3,N3,N3,N3,IOUT)
      IX=0
      IY=0
      DO 10 I=1,N3
      IF((I-1)/3+1.NE.II)GO TO 11
      IY=IY+N3
      GO TO 10
   11 CONTINUE
      DO 20 J=1,N3
      IY=IY+1
      IF((J-1)/3+1.EQ.II)GO TO 20
      IX=IX+1
      HESS(IX)=HESS(IY)
   20 CONTINUE
   10 CONTINUE
      NDO=N3-3
      WRITE(IOUT,1000)'NDEF',NDO
      CALL MATOUT(HESS,NDO,NDO,NDO,NDO,IOUT)
C
      IX=0
      DO 100 I=1,NDO
      DO 100 J=1,I
      IX=IX+1
      IJ=(I-1)*NDO+J
      HESS(IX)=HESS(IJ)
  100 CONTINUE
      RETURN
 1000 FORMAT(/5X,'SQUARE ', A8,' HESSIAN FROM CMPRS: DIMENSION =',I3)
      END
      SUBROUTINE CNSTRV(NFIXDIR,FIXDIR,VALDIR,NDEGF,GRAD,NATOMS,
     1                  GM,KILLER,KLASS,IOUT,rhs)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER FIXDIR(6,2),KLASS(2)
      DIMENSION KILLER(2),VALDIR(2),GRAD(NDEGF)
      DIMENSION GM(3,2)
      common/pscent/maxnew,psc(30,10),dpsc(30,10),ipsc(30,10)
      DATA TWO/2.D0/,ZERO/0.D0/,PI/3.141592654D0/
      rl=1.d0
C
      if(maxnew.eq.0)go to 99
      do 40 i=1,maxnew
      do 41 k=1,3
   41 gm(k,natoms+i)=0
      do 42 k=1,3
      do 42 j=1,natoms
   42 gm(k,natoms+i)=gm(k,natoms+i)+gm(k,j)*psc(j,i)
   40 continue
   99 continue
C
      N3=NATOMS*3
      DO 101 I=1,NDEGF
  101 grad(I)=ZERO
C
      MYCLASS=0
      ZNORM=ZERO
      IF(NFIXDIR.NE.1)STOP ' CONSTRV ERROR '
      DO 200 I=1,NFIXDIR
      II=FIXDIR(1,I)
      JJ=FIXDIR(2,I)
      IIU=(II-1)*3
      JJU=(JJ-1)*3
      GO TO (400,600,400,600,700),KLASS(I)
      STOP ' ILLEGAL KLASS '
  400 CONTINUE
C         DISTANCE CONTRAINT
      MYCLASS=1
      call concon0(natoms,gm,grad,ii,jj,rl,x)
      IF(KLASS(I).EQ.3)GO TO 150
      WRITE(IOUT,1028)II,JJ,DSQRT(X),VALDIR(I)
      RHS=X-VALDIR(I)**2
      ZNORM=ZNORM+RHS**2
      GO TO 200
  150 CONTINUE
C         DISTANCE DIFFERENCE CONSTRAINT
      KK=FIXDIR(3,I)
      LL=FIXDIR(4,I)
      IIU=(KK-1)*3
      JJU=(LL-1)*3
      call concon0(natoms,gm,grad,kk,ll,-rl,y)
      WRITE(IOUT,1027)II,JJ,DSQRT(X),KK,LL,DSQRT(Y)
      RHS=X-Y
      ZNORM=ZNORM+RHS**2
      GO TO 200
  700 CONTINUE
C         dihedral ANGLE CONSTRAINT
      KK=FIXDIR(3,I)
      LL=FIXDIR(6,I)
      KKU=(KK-1)*3
      LLU=(LL-1)*3
      CALL DIHED1(GM(1,II),GM(1,JJ),GM(1,KK),GM(1,LL),DIH1)
      call concon2(NATOMS,GM,SCR,II,JJ,KK,LL,RL)
      WRITE(IOUT,1030)II,JJ,KK,LL,acos(DIH1)*180.d0/pi,VALDIR(I)
      XX=VALDIR(I)*PI/180.D0
      RHS=DIH1-DCOS(XX)
      GO TO 200
  600 CONTINUE
C         ANGLE CONSTRAINT
      KK=FIXDIR(3,I)
      KKU=(KK-1)*3
      CALL ANGL(GM,II,JJ,KK,ANG,CA,R13,R12,R23)
      call concon1(NATOMS,GM,SCR,II,JJ,KK,RL)
      IF(KLASS(I).EQ.4)GO TO 170
      WRITE(IOUT,1029)II,JJ,KK,ANG,VALDIR(I)
      XX=VALDIR(I)*PI/180.D0
      RHS=CA-DCOS(XX)
      GO TO 200
  170 CONTINUE
C         ANGLE DIFFERENCE
      IIO=II
      JJO=JJ
      KKO=KK
      II=FIXDIR(4,I)
      JJ=FIXDIR(5,I)
      KK=FIXDIR(6,I)
      IIU=(II-1)*3
      JJU=(JJ-1)*3
      KKU=(KK-1)*3
      CALL ANGL(GM,II,JJ,KK,ANGB,CB,R13,R12,R23)
      CAT=CA
      CA=CB
      call concon1(NATOMS,GM,SCR,II,JJ,KK,-RL)
      WRITE(IOUT,1026)IIO,JJO,KKO,ANG,II,JJ,KK,ANGB
      RHS=CAT-CB
  200 CONTINUE
      ZNORM=DSQRT(ZNORM)
C
      KOUNT=0
      DO 800 I=1,N3
      IF(KILLER(I).EQ.1)GO TO 800
      KOUNT=KOUNT+1
      GRAD(KOUNT)=grad(I)
  800 CONTINUE
C
      WRITE(IOUT,1000)Ndegf,(GRAD(I),I=1,Ndegf)
      IF(MYCLASS.EQ.1)WRITE(IOUT,1001)ZNORM
  900 CONTINUE
      RETURN
 1000 FORMAT(5X,' CONSTRAINTS:DIMENSION=',I3,
     X     (/5X,8F12.6))
 1001 FORMAT(/5X,'DISTANCE NORM=',E15.8)
 1026 FORMAT(5X,'ANG(',I2,',',I2,',',I2,')  ='
     X       ,F12.6,' ANG(',I2,',',I2,',',I2,')  =',F12.6)
 1027 FORMAT(5X,'|R(',I2,') - R(',I2,')|=',F12.6,
     X ' |R(',I2,') - R(',I2,')|=',F12.6)
 1028 FORMAT(5X,'|R(',I2,') - R(',I2,')|=',F12.6,
     X ' CONSTRAINT VALUE=',F12.6)
 1029 FORMAT(5X,'ANG(',I2,',',I2,
     X       ',',I2,')  =',F12.6,' CONSTRAINT VALUE=',F12.6)
 1030 FORMAT(5X,'DIH ANG(',I2,',',I2,
     X       ',',I2,',',i2,')  =',F12.6,' CONSTRAINT VALUE=',F12.6)
      END
      SUBROUTINE CONCON0(NATOMS,GM,SCR,II,JJ,RLAMBDA,X)
      implicit real*8(a-h,o-z)
      dimension gm(3,10),scr(2)
      common/pscent/maxnew,psc(30,10),dpsc(30,10),ipsc(30,10)
      data two/2.d0/
C
      X = (GM(1,II)-GM(1,JJ))**2 + (GM(2,II)-GM(2,JJ))**2+
     1         (GM(3,II)-GM(3,JJ))**2 
      DO 128 K=1,3
      if(ii.le.natoms)then
      iiu=(ii-1)*3
      SCR(IIU+K)=SCR(IIU+K)+RLAMBDA*TWO*(GM(K,II)-GM(K,JJ))
      else
      i=ii-natoms
      do 111 j=1,natoms
      iiu=(j-1)*3
 111  SCR(IIU+K)=SCR(IIU+K)+RLAMBDA*TWO*(GM(K,II)-GM(K,JJ))*psc(j,i)
      endif
C    
      if(jj.le.natoms)then
      jju=(jj-1)*3
      SCR(JJU+K)=SCR(JJU+K)-RLAMBDA*TWO*(GM(K,II)-GM(K,JJ))
      else
      i=jj-natoms
      do 112 j=1,natoms
      jju=(j-1)*3
  112 SCR(JJU+K)=SCR(JJU+K)-RLAMBDA*TWO*(GM(K,II)-GM(K,JJ))*psc(j,i)
      endif
  128 CONTINUE
      return
      end     
      SUBROUTINE CONCON1(NATOMS,GM,SCR,II,JJ,KK,RLAMBDA)
      implicit real*8(a-h,o-z)
      dimension gm(3,10),scr(2)
      common/pscent/maxnew,psc(30,10),dpsc(30,10),ipsc(30,10)
      data two/2.d0/
C        STATEMENT FUNCTION GRADIENTS
      C1(X12,X23,R12,R23)=-(X23/R23+CA*X12/R12)/R12
      C2(X12,X23,R12,R23)=(X23-X12)/(R12*R23)
     1   +CA*(X12/R12**2-(X23/R23**2))
      C3(X12,X23,R12,R23)=(X12/R12+CA*X23/R23)/R23
C
      CALL ANGL(GM,II,JJ,KK,ANG,CA,R13,R12,R23)
      DO 138 K=1,3
C
      X12=GM(K,II)-GM(K,JJ)
      X23=GM(K,JJ)-GM(K,KK)
      if(ii.le.natoms)then
      iiu=(ii-1)*3
      SCR(IIU+K)=SCR(IIU+K)+RLAMBDA*C1(X12,X23,R12,R23)
      else
      i=ii-natoms
      do 111 j=1,natoms
      iiu=(j-1)*3
 111  SCR(IIU+K)=SCR(IIU+K)+RLAMBDA*C1(X12,X23,R12,R23)*psc(j,i)
      endif
C
      if(kk.le.natoms)then
      KKU=(KK-1)*3
      SCR(KKU+K)=SCR(KKU+K)+RLAMBDA*C3(X12,X23,R12,R23)
      else
      i=kk-natoms
      do 112 j=1,natoms
      kku=(j-1)*3
 112  SCR(KKU+K)=SCR(KKU+K)+RLAMBDA*C3(X12,X23,R12,R23)*psc(j,i)
      endif
C
      if(jj.le.natoms)then
      jju=(jj-1)*3
      SCR(JJU+K)=SCR(JJU+K)+RLAMBDA*C2(X12,X23,R12,R23)
      else
      i=jj-natoms
      do 113 j=1,natoms
      jju=(j-1)*3
 113  SCR(JJU+K)=SCR(JJU+K)+RLAMBDA*C2(X12,X23,R12,R23)*psc(j,i)
      endif
C
  138 CONTINUE
C
      return
      end     
      SUBROUTINE CONCON2(NATOMS,GM,SCR,II,JJ,KK,LL,RLAMBDA)
      implicit real*8(a-h,o-z)
      dimension gm(3,10),scr(2)
      common/pscent/maxnew,psc(30,10),dpsc(30,10),ipsc(30,10)
      data two/2.d0/
C

      DO 138 K=1,3
C
      call gdih(gm(1,II),gm(1,JJ),gm(1,kk),gm(1,ll),1,k,gradiik)
      if(ii.le.natoms)then
      iiu=(ii-1)*3
      SCR(IIU+K)=SCR(IIU+K)+RLAMBDA*gradiik
      else
      i=ii-natoms
      do 111 j=1,natoms
      iiu=(j-1)*3
 111  SCR(IIU+K)=SCR(IIU+K)+RLAMBDA*gradiik*psc(j,i)
      endif
C
      CALL GDIH(GM(1,II),GM(1,JJ),GM(1,KK),GM(1,LL),3,K,GRADIIK)
      if(kk.le.natoms)then
      KKU=(KK-1)*3
      SCR(KKU+K)=SCR(KKU+K)+RLAMBDA*gradIIk
      else
      i=kk-natoms
      do 112 j=1,natoms
      kku=(j-1)*3
 112  SCR(KKU+K)=SCR(KKU+K)+RLAMBDA*gradiik*psc(j,i)
      endif
C
      call gdih(gm(1,II),gm(1,JJ),gm(1,kk),gm(1,ll),2,k,gradiik)
      if(jj.le.natoms)then
      jju=(jj-1)*3
      SCR(JJU+K)=SCR(JJU+K)+RLAMBDA*gradiik
      else
      i=jj-natoms
      do 113 j=1,natoms
      jju=(j-1)*3
 113  SCR(JJU+K)=SCR(JJU+K)+RLAMBDA*gradiik*psc(j,i)
      endif
C
      call gdih(gm(1,II),gm(1,JJ),gm(1,kk),gm(1,ll),4,k,gradiik)
      if(ll.le.natoms)then
      llu=(ll-1)*3
      SCR(LLU+K)=SCR(LLU+K)+RLAMBDA*gradiik
      else
      i=ll-natoms
      do 114 l=1,natoms
      llu=(l-1)*3
 114  SCR(LLU+K)=SCR(LLU+L)+RLAMBDA*gradiik*psc(j,i)
      endif
C
  138 CONTINUE
C
      return
      end     

      subroutine  gDIH(C1,C2,C3,C4,ii,kk,grad)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION C1(3),C2(3),C3(3),C4(3)
     1          ,RGM(3,4),GM(3,4),RM1(3,3)
      DATA ZERO/0.D0/,ONE/1.D0/,PI/3.141592654D0/
      data h/0.001/
      do 10 i=1,3
      gm(i,1)=c1(i)
      gm(i,2)=c2(i)
      gm(i,3)=c3(i)
      gm(i,4)=c4(i)
 10   continue
      gm(kk,ii)=gm(kk,ii)+h
      call dihed1(gm(1,1),gm(1,2),gm(1,3),gm(1,4),dihp)
      gm(kk,ii)=gm(kk,ii)-h*2.d0
      call dihed1(gm(1,1),gm(1,2),gm(1,3),gm(1,4),dihm)
      grad=(dihp-dihm)/(2*h)
      return 
      end

      SUBROUTINE DIHEDx(C1,C2,C3,C4,DIH)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION C1(3),C2(3),C3(3),C4(3)
     1          ,RGM(3,4),GM(3,4),RM1(3,3)
      DATA ZERO/0.D0/,ONE/1.D0/,PI/3.141592654D0/
      R=ZERO
      DO 10 I=1,3
      GM(I,1)=C1(I)-C2(I)
      GM(I,2)=ZERO
      GM(I,3)=C3(I)-C2(I)
      GM(I,4)=C4(I)-C2(I)
      R=R+GM(I,3)**2
   10 CONTINUE
      R=SQRT(R)
      THETA=ACOS(GM(3,3)/R)
      if(abs(gm(1,3)).gt.1.0e-08)then
      PHI=ATAN(GM(2,3)/GM(1,3))
      if(gm(1,3).lt.0)phi=phi+pi
      else
      phi=pi/2.d0
      if(gm(2,3).lt.0)phi=-phi
      endif
      DO 11 I=1,2
      RM1(3,I)=ZERO
  11  RM1(I,3)=ZERO
      RM1(3,3)=ONE
      RM1(1,1)=COS(PHI)
      RM1(2,2)=COS(PHI)
      RM1(1,2)=SIN(PHI)
      RM1(2,1)=-SIN(PHI)
      DO 15 I=1,4
   15 CALL EBC(RGM(1,I),RM1,GM(1,I),3,3,1)
      DO 21 I=1,3
      RM1(2,I)=ZERO
   21 RM1(I,2)=ZERO
      RM1(2,2)=ONE
      RM1(1,1)=COS(THETA)
      RM1(3,3)=COS(THETA)
      RM1(1,3)=-SIN(THETA)
      RM1(3,1)= SIN(THETA)
      DO 25 I=1,4
   25 CALL EBC(GM(1,I),RM1,RGM(1,I),3,3,1)
      if(abs(gm(1,3))/r.gt.1.0e-08)go to 900
      if(abs(gm(2,3))/r.gt.1.0e-08)go to 900
      A1=ATAN(GM(2,1)/GM(1,1))
      if(gm(1,1).lt.0)a1=a1+pi
      A4=ATAN(GM(2,4)/GM(1,4))
      if(gm(1,4).lt.0)a4=a4+pi
      DIH=(A4-A1)/PI*180.
      RETURN
  900 continue
      WRITE(6,1000)gm(1,3),gm(2,3),gm(3,3),theta/pi*180.,phi/pi*180.
      DIH=-999
 1000 format('DIHEDRAL ANGLE ALGORITHM ERROR',/5x,5f12.6)
      END
      SUBROUTINE DIHED1(C1,C2,C3,C4,DIH)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION C1(3),C2(3),C3(3),C4(3)
     1          ,RGM(3,4),GM(3,4),RM1(3,3)
      DATA ZERO/0.D0/,ONE/1.D0/,PI/3.141592654D0/
      R=ZERO
      rl1=0
      rl2=0
      dot=0
      DO 10 I=1,3
      GM(I,1)=C1(I)-C2(I)
      GM(I,2)=ZERO
      GM(I,3)=C3(I)-C2(I)
      GM(I,4)=C4(I)-C2(I)
      dot=dot+gm(i,1)*(gm(i,4)-gm(i,3))
      rl1=rl1+gm(i,1)**2
      rl2=rl2+(gm(i,4)-gm(i,3))**2
      R=R+GM(I,3)**2
   10 CONTINUE
      R=SQRT(R)
      THETA=ACOS(GM(3,3)/R)
      if(abs(gm(1,3)).gt.1.0e-08)then
      PHI=ATAN(GM(2,3)/GM(1,3))
      if(gm(1,3).lt.0)phi=phi+pi
      else
      phi=pi/2.d0
      if(gm(2,3).lt.0)phi=-phi
      endif
      DO 11 I=1,2
      RM1(3,I)=ZERO
  11  RM1(I,3)=ZERO
      RM1(3,3)=ONE
      RM1(1,1)=COS(PHI)
      RM1(2,2)=COS(PHI)
      RM1(1,2)=SIN(PHI)
      RM1(2,1)=-SIN(PHI)
      DO 15 I=1,4
   15 CALL EBC(RGM(1,I),RM1,GM(1,I),3,3,1)
      DO 21 I=1,3
      RM1(2,I)=ZERO
   21 RM1(I,2)=ZERO
      RM1(2,2)=ONE
      RM1(1,1)=COS(THETA)
      RM1(3,3)=COS(THETA)
      RM1(1,3)=-SIN(THETA)
      RM1(3,1)= SIN(THETA)
      DO 25 I=1,4
   25 CALL EBC(GM(1,I),RM1,RGM(1,I),3,3,1)
      if(abs(gm(1,3))/r.gt.1.0e-08)go to 900
      if(abs(gm(2,3))/r.gt.1.0e-08)go to 900
      do 31 i=3,3
      rl1=rl1-gm(i,1)**2
      rl2=rl2-(gm(i,4)-gm(i,3))**2
      dot=dot-gm(i,1)*(gm(i,4)-gm(i,3))
   31 continue
c  (-) for consistency with other dih method
C      dih=-dot/sqrt(rl1*rl2)
c8/1/01      dih=-(dot/sqrt(rl1))/sqrt(rl2)
      dih=0.
      if(abs(dot).gt.1.0e-07)dih=-dot/sqrt(rl1*rl2)
      RETURN
  900 continue
      WRITE(6,1000)gm(1,3),gm(2,3),gm(3,3),theta/pi*180.,phi/pi*180.
      DIH=-999
 1000 format('DIHEDRAL ANGLE ALGORITHM ERROR',/5x,5f12.6)
      END
      SUBROUTINE DIHED(C1,C2,C3,C4,DIH)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION C1(3),C2(3),C3(3),C4(3)
     1          ,RGM(3,4),GM(3,4),RM1(3,3)
      DATA ZERO/0.D0/,ONE/1.D0/,PI/3.141592654D0/
      R=ZERO
      rl1=0
      rl2=0
      dot=0
      DO 10 I=1,3
      GM(I,1)=C1(I)-C2(I)
      GM(I,2)=ZERO
      GM(I,3)=C3(I)-C2(I)
      GM(I,4)=C4(I)-C2(I)
      dot=dot+gm(i,1)*(gm(i,4)-gm(i,3))
      rl1=rl1+gm(i,1)**2
      rl2=rl2+(gm(i,4)-gm(i,3))**2
      R=R+GM(I,3)**2
   10 CONTINUE
      R=SQRT(R)
      THETA=ACOS(GM(3,3)/R)
      if(abs(gm(1,3)).gt.1.0e-08)then
      PHI=ATAN(GM(2,3)/GM(1,3))
      if(gm(1,3).lt.0)phi=phi+pi
      else
      phi=pi/2.d0
      if(gm(2,3).lt.0)phi=-phi
      endif
      DO 11 I=1,2
      RM1(3,I)=ZERO
  11  RM1(I,3)=ZERO
      RM1(3,3)=ONE
      RM1(1,1)=COS(PHI)
      RM1(2,2)=COS(PHI)
      RM1(1,2)=SIN(PHI)
      RM1(2,1)=-SIN(PHI)
      DO 15 I=1,4
   15 CALL EBC(RGM(1,I),RM1,GM(1,I),3,3,1)
      DO 21 I=1,3
      RM1(2,I)=ZERO
   21 RM1(I,2)=ZERO
      RM1(2,2)=ONE
      RM1(1,1)=COS(THETA)
      RM1(3,3)=COS(THETA)
      RM1(1,3)=-SIN(THETA)
      RM1(3,1)= SIN(THETA)
      DO 25 I=1,4
   25 CALL EBC(GM(1,I),RM1,RGM(1,I),3,3,1)
      if(abs(gm(1,3))/r.gt.1.0e-08)go to 900
      if(abs(gm(2,3))/r.gt.1.0e-08)go to 900
      do 31 i=3,3
      rl1=rl1-gm(i,1)**2
      rl2=rl2-(gm(i,4)-gm(i,3))**2
      dot=dot-gm(i,1)*(gm(i,4)-gm(i,3))
   31 continue
c  (-) for consistency with other dih method
C      dih=-dot/sqrt(rl1*rl2)
      if(abs(dot).gt.1.0e-07)then
      dih=-dot/sqrt(rl1*rl2)
      dih=180.d0/pi*acos(dih)
      else
      dih=90.
      endif
      RETURN
  900 continue
      WRITE(6,1000)gm(1,3),gm(2,3),gm(3,3),theta/pi*180.,phi/pi*180.
      DIH=-999
 1000 format('DIHEDRAL ANGLE ALGORITHM ERROR',/5x,5f12.6)
      END

      SUBROUTINE DIST(NATOMS,GM,IOUT)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION GM(3,NATOMS)
      DATA TOA/0.52917715D0/
      DIJ(I,J)=DSQRT((GM(1,I)-GM(1,J))**2+(GM(2,I)-GM(2,J))**2+
     1         (GM(3,I)-GM(3,J))**2)
      WRITE(IOUT,1001)'AU      '
      DO 10 I=2,NATOMS
      WRITE(IOUT,1000)I,(DIJ(I,J),J=1,I-1)
   10 CONTINUE
      WRITE(IOUT,1001)'ANGSTROM'
      DO 20 I=2,NATOMS
      WRITE(IOUT,1000)I,(DIJ(I,J)*TOA,J=1,I-1)
   20 CONTINUE
      RETURN
 1000 FORMAT(2X,I3,2X,6F12.6,/,(7X,6F12.6))
 1001 FORMAT(/5X,'INTERNUCLEAR DISTANCES IN ',A8)
      END
C
      SUBROUTINE FANGL(NATOMS,GM,IOUT)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION GM(3,2)
      WRITE(IOUT,1000)
      DO 10 I=1,NATOMS-2
      DO 11 J=I+1,NATOMS-1
      DO 11 K=J+1,NATOMS
      CALL ANGL(GM,I,J,K,ANG1,W,X,Y,Z)
      CALL ANGL(GM,J,K,I,ANG2,W,X,Y,Z)
      CALL ANGL(GM,K,I,J,ANG3,W,X,Y,Z)
      WRITE(IOUT,1001)I,J,K,ANG1, J,K,I,ANG2, K,I,J,ANG3
   11 CONTINUE
   10 CONTINUE
      RETURN
 1000 FORMAT(/5X,'INTERNUCLEAR ANGLES IN DEGREES',/)
 1001 FORMAT(2X,3(3(I2,1X),F6.1,2X))
      END
      SUBROUTINE FDANGL(NATOMS,GM,IOUT)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION GM(3,2)
      WRITE(IOUT,1000)
      IEVE=0
      DO 10 I=1,NATOMS-3
      DO 11 J=I+1,NATOMS-2
      DO 11 K=J+1,NATOMS-1
      DO 11 L=K+1,NATOMS
      I0=I
      J0=J
      K0=K
      L0=L
      CALL DIHED(GM(1,I),GM(1,J),GM(1,K),GM(1,L),DIH)
      IEVE=IEVE+1
      IF(IAND(IEVE,1).NE.0)THEN
      I1=I
      J1=J
      K1=K
      L1=L
      DIH1=DIH
      GO TO 11
      ENDIF
      WRITE(IOUT,1001)I1,J1,K1,L1,DIH1,I,J,K,L,DIH
   11 CONTINUE
      IF(IAND(IEVE,1).NE.0)WRITE(IOUT,1001)I0,J0,K0,L0,DIH
   10 CONTINUE
      RETURN
 1000 FORMAT(/5X,'INTERNUCLEAR DIHEDRAL IN DEGREES',/)
 1001 FORMAT(2X,3(4(I2,1X),F6.1,2X))
      END
      SUBROUTINE FREEZG(NFIXDIR,FIXDIR,VALDIR,NDEGF,GRAD,NATOMS,
     1                  GM,KILLER,SCR,RLAMBDA,KLASS,IOUT)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER FIXDIR(6,2),KLASS(2)
      DIMENSION KILLER(2),VALDIR(2),GRAD(NDEGF),SCR(2),RLAMBDA(2)
      DIMENSION GM(3,2)
      common/pscent/maxnew,psc(30,10),dpsc(30,10),ipsc(30,10)
      DATA TWO/2.D0/,ZERO/0.D0/,PI/3.141592654D0/
      if(maxnew.eq.0)go to 99
      do 40 i=1,maxnew
      do 41 k=1,3
   41 gm(k,natoms+i)=0
      do 42 k=1,3
      do 42 j=1,natoms
   42 gm(k,natoms+i)=gm(k,natoms+i)+gm(k,j)*psc(j,i)
   40 continue
   99 continue
C
      IF(NFIXDIR.EQ.0)GO TO 900
      NF=NDEGF+NFIXDIR
      N3=NATOMS*3
      DO 101 I=1,N3
  101 SCR(I)=ZERO
C
      MYCLASS=0
      ZNORM=ZERO
      DO 200 I=1,NFIXDIR
      II=FIXDIR(1,I)
      JJ=FIXDIR(2,I)
      IIU=(II-1)*3
      JJU=(JJ-1)*3
      GO TO (400,600,400,600,700),KLASS(I)
      STOP ' ILLEGAL KLASS '
  400 CONTINUE
C         DISTANCE CONTRAINT
      MYCLASS=1
      call concon0(natoms,gm,scr,ii,jj,rlambda(i),x)
      IF(KLASS(I).EQ.3)GO TO 150
      WRITE(IOUT,1028)II,JJ,DSQRT(X),VALDIR(I)
      GRAD(NDEGF+I)=X-VALDIR(I)**2
      ZNORM=ZNORM+GRAD(NDEGF+I)**2
      GO TO 200
  150 CONTINUE
C         DISTANCE DIFFERENCE CONSTRAINT
      KK=FIXDIR(3,I)
      LL=FIXDIR(4,I)
      call concon0(natoms,gm,scr,kk,ll,-rlambda(i),y)
      WRITE(IOUT,1027)II,JJ,DSQRT(X),KK,LL,DSQRT(Y)
      GRAD(NDEGF+I)=X-Y
      ZNORM=ZNORM+GRAD(NDEGF+I)**2
      GO TO 200
  700 CONTINUE
C         dihedral ANGLE CONSTRAINT
      KK=FIXDIR(3,I)
      LL=FIXDIR(6,I)
      KKU=(KK-1)*3
      LLU=(LL-1)*3
      CALL DIHED1(GM(1,II),GM(1,JJ),GM(1,KK),GM(1,LL),DIH1)
      call concon2(NATOMS,GM,SCR,II,JJ,KK,LL,RLAMBDA(I))
      WRITE(IOUT,1030)II,JJ,KK,LL,acos(DIH1)*180.d0/pi,VALDIR(I)
      XX=VALDIR(I)*PI/180.D0
      GRAD(NDEGF+I)=dih1-DCOS(XX)
      GO TO 200
  600 CONTINUE
C         ANGLE CONSTRAINT
      KK=FIXDIR(3,I)
      KKU=(KK-1)*3
      CALL ANGL(GM,II,JJ,KK,ANG,CA,R13,R12,R23)
      call concon1(NATOMS,GM,SCR,II,JJ,KK,RLAMBDA(I))
      IF(KLASS(I).EQ.4)GO TO 170
      WRITE(IOUT,1029)II,JJ,KK,ANG,VALDIR(I)
      XX=VALDIR(I)*PI/180.D0
      GRAD(NDEGF+I)=CA-DCOS(XX)
      GO TO 200
 170  continue
C         ANGLE DIFFERENCE
      IIO=II
      JJO=JJ
      KKO=KK
      II=FIXDIR(4,I)
      JJ=FIXDIR(5,I)
      KK=FIXDIR(6,I)
      IIU=(II-1)*3
      JJU=(JJ-1)*3
      KKU=(KK-1)*3
      CALL ANGL(GM,II,JJ,KK,ANGB,CB,R13,R12,R23)
      CAT=CA
      CA=CB
      call concon1(NATOMS,GM,SCR,II,JJ,KK,-RLAMBDA(I))
      WRITE(IOUT,1026)IIO,JJO,KKO,ANG,II,JJ,KK,ANGB
      GRAD(NDEGF+I)=CAT-CB
  200 CONTINUE
      ZNORM=DSQRT(ZNORM)
C
      WRITE(IOUT,1000)'WITHOUT',NDEGF,(GRAD(I),I=1,NDEGF)
      KOUNT=0
      DO 800 I=1,N3
      IF(KILLER(I).EQ.1)GO TO 800
      KOUNT=KOUNT+1
      GRAD(KOUNT)=SCR(I)+GRAD(KOUNT)
  800 CONTINUE
C
      WRITE(IOUT,1000)'WITH',NF,(GRAD(I),I=1,NF)
      IF(MYCLASS.EQ.1)WRITE(IOUT,1001)ZNORM
  900 CONTINUE
      RETURN
 1000 FORMAT(5X,'GRADIENT ',A7,' CONSTRAINTS:DIMENSION=',I3,
     X     (/5X,8F12.6))
 1001 FORMAT(/5X,'DISTANCE NORM=',E15.8)
 1026 FORMAT(5X,'ANG(',I2,',',I2,',',I2,')  ='
     X       ,F12.6,' ANG(',I2,',',I2,',',I2,')  =',F12.6)
 1027 FORMAT(5X,'|R(',I2,') - R(',I2,')|=',F12.6,
     X ' |R(',I2,') - R(',I2,')|=',F12.6)
 1028 FORMAT(5X,'|R(',I2,') - R(',I2,')|=',F12.6,
     X ' CONSTRAINT VALUE=',F12.6)
 1029 FORMAT(5X,'ANG(',I2,',',I2,
     X       ',',I2,')  =',F12.6,' CONSTRAINT VALUE=',F12.6)
 1030 FORMAT(5X,'DIH ANG(',I2,',',I2,
     X       ',',I2,',',i2,')  =',F12.6,' CONSTRAINT VALUE=',F12.6)
      END
      SUBROUTINE GAUSS(B,NMIX)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION B(1)
      DATA ZERO/1.0D-12/
      IND(I)=(I-1)*NMIX1

C           GAUSSIAN ELIMINATION
      NMIX1=NMIX+1
C
C           BEGIN THE ELIMINATION
C
      DO 140 I=1,NMIX
      M=IND(I)
      IF(DABS(B(M+I)).LT.ZERO)GO TO 70
      FX=1./B(M+I)
      GO TO 95
  70  CONTINUE
      IF(I.EQ.NMIX)GO TO 1000
      I1=I+1
C           PIVOT SECTION
      DO 90 J=I1,NMIX
      INJ=IND(J)
      IF(DABS(B(INJ+I)).LT.ZERO)GO TO 90
      FX=1./B(INJ+I)
      DO 85 L=I,NMIX1
      TEMP=B(INJ+L)
      B(INJ+L)=B(M+L)
      B(M+L)=TEMP
  85  CONTINUE
      GO TO 95
  90  CONTINUE
      GO TO 1000
  95  CONTINUE
      DO 100 J=I,NMIX1
      B(M+J)=B(M+J)*FX
 100  CONTINUE
C         END PIVOT
      DO 130 J=1,NMIX
      IF(I.EQ.J)GO TO 130
      L=IND(J)
      Y=B(L+I)
      DO 120 K=I,NMIX1
      B(L+K)=B(L+K)-Y*B(M+K)
 120  CONTINUE
 130  CONTINUE
 140  CONTINUE
C
C       MOVE THE SOLUTION THE TO BEGINNING OF CORE
C
      DO 150 I=1,NMIX
      B(I)=B(IND(I)+NMIX1)
 150  CONTINUE
      RETURN
1000  CONTINUE
      WRITE(6,1010)
1010  FORMAT('  ABORT GAUSS  SINGULAR MATRIX ')
      STOP
      END
      SUBROUTINE PUTPROP(IPRPU,IOUTU,NPROP,IPRU,HIJ,NDEGF,change)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER FIXDIR
      DIMENSION IHDG(17),JHDG(28),TITLE(2,17),PR(11)
      REAL*8 LAMBDA(30),LAMBDAC(2),LAMBDAR(30)
      DIMENSION GMC(2),GGRADC(2),HESSC(2),EGRADC(2),XGRADC(2)
      DIMENSION GM(30),GGRAD(30),HESS(465),EGRAD(30),XGRAD(30)
      DIMENSION DIPOLCI(50),HIJ(2),dipolr(2),xgradr(2)
      DATA RHSN/100.D0/,ESTATE/0.0D0/
      NAMELIST/PROGRS/GM,EDIF,NEWDIR,EGRAD,GGRAD,XGRAD,
     1                SCALE,HESS,LAMBDA,LAMBDAR,RHSN,ESTATE1
      NAMELIST/PRPRTY/DIPOLCI,XGRAD
C
      DATA TITLE/'POTENTIA','L       ','ELECTIC ','FIELD   ','FIELD GR',
     1'ADIENT  ','DIPOLE M','OMENT   ','QUAD. MO','MENT    ','THIRD MO',
     2'MENT    ','PLANAR D','ENSITY  ','LINE DEN','SITY    ','CHARGE D',
     3'ENSITY  ','OVERLAP ','        ','SECOND M','OMENT   ','DIAMAGNE',
     4'TIC SHLD','KINETIC ','        ','MODEL PO','TENTIAL ','NET POPU',
     5'LATION  ','GROSS PO','PULATION','OVLP POP','ULATION '/
      DATA JHDG/'1/R ','OVLP','DX  ','DXY ','DEN ','X   ','Y   ','Z   ',
     1'XX  ','YY  ','ZZ  ','XY  ','XZ  ','YZ  ','RSQ ','XXX ','YYY ',
     2'ZZZ ','XXY ','XXZ ','XYY ','YYZ ','XZZ ','YZZ ','XYZ ',
     3'    ','A   ','B   '/
      DATA IHDG/1,6,9,6,9,16,3,4,5,2,9,9,26,27,26,26,26/
      DO 10 I=1,50
      IF(I.LE.30)XGRAD(I)=0.D0
   10 DIPOLCI(I)=0.D0
C
      OPEN(IPRPU,FORM='UNFORMATTED')
      REWIND IPRPU
      READ(IPRPU,END=200)NP,IST,JST,OVIJ
      DO 30 N=1,NP
      READ(IPRPU)NC,NT
      IF(NC.EQ.NPROP)GO TO 31
   30 CONTINUE
      GO TO 200
   31 CONTINUE
C
      REWIND IPRPU
      IA=IHDG(NPROP)
      IB=IA+NT-1
      WRITE (IOUTU,400) (TITLE(I,NPROP),I=1,2)
      WRITE (IOUTU,420) (JHDG(I),I=IA,IB)
C
      JB=0
  100 CONTINUE
      READ(IPRPU,END=200)NP,IST,JST,OVIJ
      DO 110 N=1,NP
      READ(IPRPU)NC,NT,NX,(PR(J),J=1,NT)
C
      IF(NC.NE.NPROP)GO TO 110
      JA=JB+1
      JB=JB+NT
      DO 111 I=1,NT
  111 DIPOLCI(I+JA-1)=PR(I)
      WRITE (IOUTU,434)IST,JST,(DIPOLCI(J),J=JA,JB)
  110 CONTINUE
      GO TO 100
  200 CONTINUE
C            REPOSITION PROGRESS FILE
      REWIND IPRU
      READ(IPRU,PROGRS,END=900)
C          RESTORE ORIGINAL PHASE TO HIJ
      DO 11 I=1,NDEGF
   11 XGRAD(I)=HIJ(I)*CHANGE
C
      WRITE(IPRU,PRPRTY)
      RETURN
  900 CONTINUE
      WRITE(IOUTU,1000)
      RETURN
C
      ENTRY GETPROP1(IPRU,DIPOLR,XGRADR,NDEGF,IOUTU)
      REWIND IPRU
      READ(IPRU,PRPRTY,END=901)
      DO 803 I=1,50
  803 DIPOLR(I)=DIPOLCI(I)
      DO 804 I=1,NDEGF
  804 XGRADR(I)=XGRAD(I)
      RETURN
  901 CONTINUE
      WRITE(IOUTU,1001)
      RETURN
C
      ENTRY PUTPROP1(IPRU,DIPOLR,XGRADR,NDEGF)
      REWIND IPRU
      READ(IPRU,PROGRS,END=900)
      DO 801 I=1,50
  801 DIPOLCI(I)=DIPOLR(I)
      DO 802 I=1,NDEGF
  802 XGRAD(I)=XGRADR(I)
      WRITE(IPRU,PRPRTY)
      RETURN
C
  400 FORMAT(20X,2A8)
  420 FORMAT (3X,'STATES',2X,10(8X,A4,2X))
  434 FORMAT (3X,2I3,2X,5F14.8/(11X,5F14.8))
 1000 FORMAT(/5X,' NO PROGRESS FILE PROPERTIES NOT APPENDED')
 1001 FORMAT(/5X,' NO PRPRTY FILE PROPERTIES NOT APPENDED')
      END
      SUBROUTINE GETDDH(IPRU,NDEGF,EGRADP,GGRADP,XGRADP,
     1        EGRADPP,GGRADPP,XGRADPP,IOUT,ICDF,LAMBDAC,NFIXDIR)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 LAMBDAC(2),LAMBD0(2)
      DIMENSION EGRADPP(NDEGF),GGRADPP(NDEGF),EGRADP(NDEGF,2),
     1          GGRADP(NDEGF,2),XGRADP(NDEGF,2),XGRADPP(NDEGF)
      DIMENSION EGRADH(930),GGRADH(930),GM(30),DISP(30),XGRADH(930)
      NAMELIST/PROGRSH/GGRADH,EGRADH,XGRADH,IDISP,GM,EDIF,DISP,
     X                 LAMBD0,EREF
      MAXDO=NDEGF+1
      IF(ICDF.GT.0)MAXDO=NDEGF*2+1
      REWIND IPRU
      READ(IPRU,PROGRSH,END=900)
      IX=0
      DO 125 J=1,MAXDO
      IF(J.LE.NDEGF)THEN
      XGRADPP(J)=XGRADH(J)
      GGRADPP(J)=GGRADH(J)
      EGRADPP(J)=EGRADH(J)
      ENDIF
      DO 125 I=1,NDEGF
      IX=IX+1
      XGRADP(I,J)=XGRADH(IX)
      GGRADP(I,J)=GGRADH(IX)
      EGRADP(I,J)=EGRADH(IX)
  125 CONTINUE
      LAMBDAC(NFIXDIR+1)=LAMBD0(1)
      LAMBDAC(NFIXDIR+2)=LAMBD0(2)
      RETURN
  900 CONTINUE
      WRITE(IOUT,1000)
      CALL EXIT(100)
 1000 FORMAT(5X,'CANT FIND DD-HESSIAN DATA - ABORT' )
      END
      SUBROUTINE GETEGG(IPRU,NDEGF,EGRADP,GGRADP,XGRADP,IOUT,ICDF)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 LAMBD0(2)
      DIMENSION EGRADP(NDEGF,2),GGRADP(NDEGF,2),XGRADP(NDEGF,2)
      DIMENSION EGRADH(930),XGRADH(930),GGRADH(930),GM(30),DISP(30)
      NAMELIST/PROGRSH/GGRADH,EGRADH,IDISP,XGRADH,
     1                 GM,EDIF,DISP,LAMBD0,EREF
      MAXDO=NDEGF+1
      IF(ICDF.GT.0)MAXDO=NDEGF*2+1
      REWIND IPRU
      READ(IPRU,PROGRSH,END=900)
      IX=0
      DO 225 J=1,MAXDO
      DO 225 I=1,NDEGF
      IX=IX+1
      XGRADP(I,J)=XGRADH(IX)
      GGRADP(I,J)=GGRADH(IX)
      EGRADP(I,J)=EGRADH(IX)
  225 CONTINUE
      RETURN
  900 CONTINUE
      WRITE(IOUT,1000)
      CALL EXIT(100)
 1000 FORMAT(5X,'CANT FIND DD-HESSIAN DATA - ABORT' )
      END
      SUBROUTINE GETGM(IPRU,NATOMS,GMC,IOUTU)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER FIXDIR
      REAL*8 LAMBDA(30),LAMBDAR(30)
      DIMENSION GMC(2)
      DIMENSION GM(30),GGRAD(30),HESS(465),EGRAD(30),XGRAD(30) 
      DATA RHSN/100.D0/,ESTATE/0.0D0/
      NAMELIST/PROGRS/GM,EDIF,NEWDIR,EGRAD,GGRAD,SCALE,HESS,
     1                LAMBDA,RHSN,ESTATE1,LAMBDAR,XGRAD
C
      WRITE(IOUT,"(1x,A)") "REWINDING PROGRESS FILE"
      REWIND IPRU
      WRITE(IOUT,"(1x,A)") "READING PROGRS NAMELIST..."
      READ(IPRU,PROGRS,END=950)
      WRITE(IOUT,"(1x,A)") "READ COMPLETE."
      DO 105 I=1,3*NATOMS
      GMC(I)=GM(I)
  105 CONTINUE
      RETURN
  950 CONTINUE
      WRITE(IOUT,1001)
      CALL EXIT(100)
 1001 FORMAT(5X,'NO PROGRESS FILE- ABORT' )
      END
C
      SUBROUTINE GETIGM(IPRU,NATOMS,GMC,IOUTU)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 LAMBD0(2)
      DIMENSION GMC(2),GGRADC(2),HESSC(2),EGRADC(2)
      DIMENSION GM(30),GGRAD(30),EGRAD(30)
      DIMENSION GGRADH(930),EGRADH(930),DISP(30),XGRADH(930)
      DATA RHSN/100.D0/,ESTATE/0.0D0/
      NAMELIST/PROGRSH/GGRADH,EGRADH,IDISP,GM,EDIF,DISP,LAMBD0,EREF
     X                ,XGRADH
C
      REWIND IPRU
      READ(IPRU,PROGRSH,END=950)
      DO 109 I=1,3*NATOMS
      GMC(I)=GM(I)
  109 CONTINUE
      RETURN
  950 CONTINUE
      WRITE(IOUT,1001)
      CALL EXIT(100)
 1001 FORMAT(5X,'NO PROGRESS FILE- ABORT' )
      END
C
      SUBROUTINE GETLM(IPRU,NFIXDIR,LAMBDAC,IOUT,IRSTRT)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 LAMBDA(30),LAMBDAC(2),LAMBDAR(30)
      DIMENSION GM(30),GGRAD(30),HESS(465),EGRAD(30),XGRAD(30)
      NAMELIST/PROGRS/GM,EDIF,NEWDIR,EGRAD,GGRAD,SCALE,HESS,
     1                LAMBDA,RHSN,ESTATE1,LAMBDAR,XGRAD
      REWIND IPRU
      READ(IPRU,PROGRS,END=950)
      DO 105 I=1,NFIXDIR
      LAMBDAC(I)=LAMBDA(I)
      IF(IRSTRT.EQ.0)GO TO 105
      LAMBDAC(I)=LAMBDAR(I)
  105 CONTINUE
      RETURN
  950 CONTINUE
      WRITE(IOUT,1001)
      DO 106 I=1,NFIXDIR
      LAMBDAC(I)=0.D0
  106 CONTINUE
      RETURN
 1001 FORMAT(5X,'NO PROGRESS FILE- LM SET TO 0' )
      END
      SUBROUTINE GETHPRG(IPRU,NATOMS,GMC,IDISPC,IOUT)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 LAMBD0(2)
      DIMENSION GMC(2),GGRADC(2),HESSC(2),EGRADC(2)
      DIMENSION GM(30),GGRAD(30),HESS(465),EGRAD(30)
      DIMENSION GGRADH(930),EGRADH(930),DISP(30),XGRADH(930)
      DATA RHSN/100.D0/,ESTATE/0.0D0/
      NAMELIST/PROGRSH/GGRADH,EGRADH,IDISP,GM,EDIF,DISP,LAMBD0,EREF
     X                ,XGRADH
C
      REWIND IPRU
      READ(IPRU,PROGRSH,END=900)
      DO 110 I=1,3*NATOMS
  110 GMC(I)=GM(I)
      IDISPC=IDISP
      RETURN
  900 CONTINUE
      WRITE(IOUT,1001)
      CALL EXIT(100)
 1001 FORMAT(5X,'NO PROGRESS FILE- ABORT' )
      END
C
      SUBROUTINE GETPRG(IPRU,NATOMS,NDEGF,GMC,EGRADC,GGRADC,XGRADC,
     1                  EDIFC,ES1,LAMBDAC,IOUT,NFIXDIR)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER FIXDIR
      REAL*8 LAMBDA(30),LAMBDAC(2),LAMBDAR(30)
      DIMENSION GMC(2),GGRADC(2),HESSC(2),EGRADC(2),XGRADC(2)
      DIMENSION GM(30),GGRAD(30),HESS(465),EGRAD(30),XGRAD(30)
      DATA RHSN/100.D0/,ESTATE/0.0D0/
      NAMELIST/PROGRS/GM,EDIF,NEWDIR,EGRAD,GGRAD,XGRAD,
     1                SCALE,HESS,LAMBDA,LAMBDAR,RHSN,ESTATE1
C
C            RESTART PATH USE OLD LAGRANGE MULTIPLIERS
      REWIND IPRU
      READ(IPRU,PROGRS,END=900)
      ES1=ESTATE1
      DO 106 I=1,NDEGF
      XGRADC(I)=XGRAD(I)
      EGRADC(I)=EGRAD(I)
  106 GGRADC(I)=GGRAD(I)
      DO 107 I=1,3*NATOMS
  107 GMC(I)=GM(I)
      DO 108 I=1,NFIXDIR+2
  108 LAMBDAC(I)=LAMBDAR(I)
      EDIFC=EDIF
      RETURN
  900 CONTINUE
      WRITE(IOUT,1001)
      CALL EXIT(100)
 1001 FORMAT(5X,'NO PROGRESS FILE- ABORT' )
      END
C
      SUBROUTINE GETHES(IPRU,NDEGF,HESSC,IOUT)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 LAMBDA(30),LAMBDAR(30)
      DIMENSION HESSC(2)
      DIMENSION GM(30),GGRAD(30),HESS(465),EGRAD(30),XGRAD(30)
      DATA RHSN/100.D0/,ESTATE/0.0D0/
      NAMELIST/PROGRS/GM,EDIF,NEWDIR,EGRAD,GGRAD,XGRAD,SCALE,HESS,
     1                LAMBDA,RHSN,ESTATE1,LAMBDAR
C
      REWIND IPRU
      READ(IPRU,PROGRS,END=950)
      NDEGF2=NDEGF*(NDEGF+1)/2 
      DO 108 IJ=1,NDEGF2
      HESSC(IJ)=HESS(IJ)
  108 CONTINUE
      RETURN
  950 CONTINUE
      WRITE(IOUT,1001)
      CALL EXIT(100)
 1001 FORMAT(5X,'NO PROGRESS FILE- ABORT' )
      END
      SUBROUTINE GETSYMT(NATOM,N3,SYMTRN,CHAR,SCR,I76,IOUT,KSEC,
     X                  IORIGIN,NTSYM,NTSYMNT,NOP)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION JSEC(50),SYMTRN(N3,N3),CHAR(N3,11),SCR(2)
      INTEGER RDI
      MASK=2*16-1
C
      OPEN(I76,FORM='FORMATTED')
      REWIND I76
      READ(I76,1001,END=500)MDF,NATOMS
      READ(I76,1001)(IZ,I=1,3*NATOMS)
      READ(I76,1001)(IZ,I=1,MDF)
      READ(I76,1001)(IZ,I=1,MDF)
      IF(NATOM.NE.NATOMS)GO TO 900
C
C         SYMMETRY EQUIVALENT CENTER INSTRUCTIONS
      KSEC=0
      DO 310 I=1,NATOMS
      READ(I76,1001)NDO
      KSEC=KSEC+NDO
      IF(NDO.EQ.0)GO TO 310
      READ(I76,1002)(JSEC(J),J=1,NDO)
  310 CONTINUE
      IF(KSEC.EQ.0)GO TO 500
C        GET SYMTRN AND CHARACTERS
C        REORDER SYMMETRIC BASIS FUNCTIONS TO PUT FIXED CENTER LAST
      READ(I76,1003)NATM3,NOP,NTSYM,NNTTSYM
      NTSYMNT=0
      IX=(IORIGIN-1)*3
      ISA=0
      ISB=NTSYM
      DO 10 I=1,NATM3
C         READ CHARACTERS
C      READ(I76,1006)I,(SCR(J),J=1,NOP+3)
      READ(I76,1006)RDI,(SCR(J),J=1,NOP+3)
      DO 311 J=1,NOP
  311 IF(DABS(SCR(J)).GT.1.0D-04)GO TO 350
      ISA=ISA+1
      DO 312 J=1,NOP+3
  312 CHAR(ISA,J)=SCR(J)
      READ(I76,1004)(SYMTRN(J,ISA),J=1,NATM3)
      DO 29 K=1,3
      IF(SYMTRN(IX+K,ISA).NE.0.)GO TO 250
   29 CONTINUE
      GO TO 10
  250 CONTINUE
C        KILLED CENTER PUT IT AT END OF TOTALLY SYMMETRIC GROUP
      DO 240 J=1,N3
  240 SYMTRN(J,NTSYM-NTSYMNT)=SYMTRN(J,ISA)
      DO 241 J=1,NOP+3
  241 CHAR(NTSYM-NTSYMNT,J)=CHAR(ISA,J)
      NTSYMNT=NTSYMNT+1
      ISA=ISA-1
      GO TO 10
  350 CONTINUE
C        NON TOTAL SYMMETRIC BASIS FUNCTION
      ISB=ISB+1
      DO 351 J=1,NOP+3
  351 CHAR(ISB,J)=SCR(J)
      READ(I76,1004)(SYMTRN(J,ISB),J=1,NATM3)
   10 CONTINUE
      IF(ISA+NTSYMNT.NE.NTSYM)STOP ' COUNTING PROBLEM IN GETSYMT '
      NTSYMNT=NTSYM-NTSYMNT
C
  500 CONTINUE
      WRITE(IOUT,1005)NATM3,NOP,NTSYM,NNTTSYM,NTSYMNT
      RETURN
C
  900 CONTINUE
      WRITE(IOUT,1000)I76,NATOM,NATOMS
 1000 FORMAT(5X,'UNIT=',I3,' ERROR NATOM=',I3,' NATOMS=',I3)
 1001 FORMAT(16I4)
 1002 FORMAT(8(1X,Z8))
 1003 FORMAT(14I5)
 1004 FORMAT(7F10.6)
 1005 FORMAT(/15X,'FROM GETSYMT:N3=',I3,' NOP=',I2,' NTSYM=',I3,
     X        ' NNTTSYM=',I3,' NTSYMNT=',I3)
 1006 FORMAT(I3,(7F10.6))
      END
      SUBROUTINE GMWRIT(NATOMS,COORD,IOUT,NSFIG,I75)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION COORD(3,2),GM(3,30)
      NAMELIST/PATH/GM
      DO 10 J=1,NATOMS
      DO 20 I=1,3
   20 GM(I,J)=COORD(I,J)
   10 CONTINUE
      WRITE(IOUT,1000)
      WRITE(IOUT,1001)((GM(J,I),J=1,3),I=1,NATOMS)
      OPEN(I75,FORM='FORMATTED')
      REWIND I75
      WRITE(I75,1003)
      WRITE(I75,1004)'PATH    '
      IF(NSFIG.EQ.3)WRITE(I75,1009)
     X'      GM',((GM(J,I),J=1,3),I=1,NATOMS)
      IF(NSFIG.EQ.4)WRITE(I75,1008)
     X'      GM',((GM(J,I),J=1,3),I=1,NATOMS)
      IF(NSFIG.EQ.5)WRITE(I75,1006)
     X'      GM',((GM(J,I),J=1,3),I=1,NATOMS)
      IF(NSFIG.EQ.6)WRITE(I75,1007)
     X'      GM',((GM(J,I),J=1,3),I=1,NATOMS)
       IF((NSFIG.GT.6).OR.(NSFIG.EQ.0))WRITE(I75,1010)
     X'      GM',((GM(J,I),J=1,3),I=1,NATOMS)
      WRITE(I75,1005)
      WRITE(IOUT,1002)I75,NSFIG
      RETURN
 1000 FORMAT(/10X,'NAMELIST INPUT FOR CURRENT GEOMETRY',/)
 1001 FORMAT(3X,'GM=',3(F15.8,','),/(6X,3(F15.8,',')))
 1002 FORMAT(5X,'PATH NAMELIST WRITTEN TO UNIT=',I3,' NSFIG=',I2)
 1003 FORMAT('****PATH')
 1004 FORMAT(1X,'&',A8)
 1005 FORMAT(1X,'&END')
CDRY 1006 FORMAT(1X,A8,'=',5(E20.10,',')/(11X,(5(E20.10,','))))
 1006 FORMAT(1X,A8,'=',3(F12.5,',')/(11X,(3(F12.5,','))))
 1007 FORMAT(1X,A8,'=',3(F12.6,',')/(11X,(3(F12.6,','))))
 1008 FORMAT(1X,A8,'=',3(F12.4,',')/(11X,(3(F12.4,','))))
 1009 FORMAT(1X,A8,'=',3(F12.3,',')/(11X,(3(F12.3,','))))
 1010 FORMAT(1X,A8,'=',3(E20.10,',')/(11X,(3(E20.10,',')))) 
      END
      SUBROUTINE IBASIS(NR,R,V,N,BF,SCR,IOUT,IPFLG,IFAIL,NFIXDIR)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION R(N,NR),SCR(2),BF(2),V(N,2)
      COMMON/VIBORD/IVBORD1(30),IVBORD2(30)
      NF=N+NFIXDIR
      NDO=NF**2 
      DO 5 I=1,NDO
    5 SCR(I)=0.D0
      IJ=0
      IF(NR.EQ.0)GO TO 11
      DO 10 I=1,NR
      DO 20 J=1,N
      SCR(IJ+J)=R(J,I)
   20 CONTINUE
      IJ=IJ+NF
   10 CONTINUE
   11 CONTINUE
C         LAGRANGE MULTIPLIER SLOTS
      IF(NFIXDIR.NE.0)THEN
      DO 27 KFIXDIR=1,NFIXDIR
      SCR(IJ+N+KFIXDIR)=1.0
      IJ=IJ+NF
   27 CONTINUE
      ENDIF
C         VIBRATION SLOTS
      NDO=N-NR
      DO 40 I=1,NDO
      DO 30 J=1,N
      SCR(IJ+J)=V(J,I)
   30 CONTINUE
      IJ=IJ+NF
   40 CONTINUE
C
      IF(IPFLG.LE.0)GO TO 46
      IHI=0
      WRITE(IOUT,1002)
      DO 45 I=1,NF
      LOW=IHI+1
      IHI=IHI+NF
      WRITE(IOUT,1000)I,(SCR(K),K=LOW,IHI)
   45 CONTINUE
   46 CONTINUE
C
C CLM: this print is for debugging
      print *, "NF=", nf
      do i=1, NF
         write(*,"(1x,i3,15f15.8)") i, SCR((i-1)*NF:i*NF)
      end do
      CALL SCHMO(NF,SCR,IOUT,IFAIL)
C
      IF(IFAIL.EQ.0)GO TO 70
      IF(IFAIL.LE.(NR+NFIXDIR))GO TO 900
      IPUT=IFAIL-(NR+NFIXDIR)
      WRITE(IOUT,1003)IVBORD2(IPUT),IVBORD1(IVBORD2(IPUT))
      IVBORD1(IVBORD2(IPUT))=-1
      RETURN
   70 CONTINUE
      WRITE(IOUT,1004)((IVBORD2(I),IVBORD1(IVBORD2(I))),I=1,N-NR)
C
      IX=0
      IY=NR*NF
      NDO=NF-NR
      DO 50 I=1,NDO+1
      DO 50 J=1,NF
      IX=IX+1
      IY=IY+1
   50 BF(IX)=SCR(IY)
C
      IHI=0
      WRITE(IOUT,1001)
      DO 60 I=1,NDO
      LOW=IHI+1
      IHI=IHI+NF
      WRITE(IOUT,1000)I,(BF(K),K=LOW,IHI)
   60 CONTINUE
      RETURN
 900  CONTINUE
      WRITE(IOUT,1003)IFAILED
      STOP ' SCHMO FAILURE IN IBASIS '
 1000 FORMAT(2X,I3,(8F10.6))
 1001 FORMAT(' ON - BASIS FOR HESSIAN ')
 1002 FORMAT(' ORIGINAL - BASIS FOR HESSIAN ')
 1003 FORMAT(5X,'SCHMO FAILURE IN IBASIS:BASIS FUNCTION=',2I3)
CDRY I2 FORMAT FAILED 3/25/92 IN 1004
C    VALUES LOOK OK - SCARY
 1004 FORMAT(/5X,'INTERNAL STATE BASIS',/,(5(F3.0,' (',F3.0,') ') ))
C 1004 FORMAT(/5X,'INTERNAL STATE BASIS',/,5(I2,' (',I2,') ') )
       END      
      SUBROUTINE INTMO(NATOMS,GM,DISP,KILLER)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION GM(2),DISP(2),KILLER(2)
      NDO=NATOMS*3
      IX=0
      DO 10 I=1,NDO
      IF(KILLER(I).EQ.1)GO TO 10
      IX=IX+1
      GM(I)=GM(I)+DISP(IX)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE IIDOF(NATOMS,MOTION,NDIAT,KILLER,Z,IOUT,IONE)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION KILLER(*),NDIAT(2),A(2),SCR(2)
      INTEGER Z
C         ELIMINATE FORMAL DOF (1,2,3)
C         AND Z-AXIS IF CS SYMMETRY
C         OR X-Y MOTION FOR LINEAR SYSTEM
      DO 10 I=1,NATOMS
      IF(NDIAT(I).EQ.1)IONE=I
      DO 10 J=1,3
      IJ=(I-1)*3+J
      IPUT=0
      IF((MOTION.EQ.2).AND.(J.EQ.Z))IPUT=1
      IF((MOTION.EQ.1).AND.(J.NE.Z))IPUT=1
   10 KILLER(IJ)=IPUT
      DO 11 I=1,3
      IPUT=(IONE-1)*3+I
   11 KILLER(IPUT)=1
C
      WRITE(IOUT,1000)(KILLER(I),I=1,NATOMS*3)
      RETURN
C
      ENTRY IDOF(NATOMS,A,SCR,KILLER)
C
      F=-1.D0/DFLOAT(NATOMS)
      DO 200 I=1,3*NATOMS
  200 SCR(I)=0.D0
C
      DO 230 L=1,3
      DO 220 I=1,NATOMS
      IXYZ=(I-1)*3+L
      SCR(IXYZ)=SCR(IXYZ)+A(IXYZ)
      DO 220 J=1,NATOMS
      JXYZ=(J-1)*3+L
      SCR(JXYZ)=SCR(JXYZ)+A(IXYZ)*F
  220 CONTINUE
  230 CONTINUE
C
      IX=0
      DO 210 I=1,NATOMS*3
      IF(KILLER(I).EQ.1)GO TO 210
      IX=IX+1
      A(IX)=SCR(I)
  210 CONTINUE
      RETURN
 1000 FORMAT(5X,'CARTESIAN TO INTERNAL DOF MAPPING',/(10(2X,3I2)))
      END
      SUBROUTINE JNTMO(NATOMS,GM,VECT,XDSP)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION GM(2),VECT(2)
      DO 10 I=1,NATOMS*3
      GM(I)=GM(I)+XDSP*VECT(I)
   10 CONTINUE
      RETURN
      END
      subroutine mkcp(killer,natoms,iorg,ipru,iout)
      implicit real*8(a-h,o-z)
      dimension gvs2(30),gv(30),hv(30),sv(30),killer(3,10),
     x          gvn(30),hvn(30)
      pi=acos(-1.d0)
      nat3=natoms*3
      call gtfgh(IPRU,IOUT,GVs2,GV,HV,SV,IORG,EREF,
     X                  KILLER,NATOMS)
      call mkRGH(SV,GVS2,HV,GVN,HVN,NAT3,beta)
      DPGHN=DOT(GVN,HVN,NAT3)
      GVND=SQRT(DOT(GVN,GVN,NAT3))
      HVND=SQRT(DOT(HVN,HVN,NAT3))
      HXR=DOT(HVN,GVN,NAT3)/GVND
      HYR=DOT(HVN,HVN,NAT3)/HVND
      GXR=DOT(GVN,GVN,NAT3)/GVND
      GYR=DOT(GVN,HVN,NAT3)/HVND
      SXR=DOT(SV,GVN,NAT3)/GVND
      SYR=DOT(SV,HVN,NAT3)/HVND
      WRITE(iout,1044)BETA/PI*180.,DPGHN,GVND,HVND, GVND*HVND
      WRITE(iout,1014)'ROTD',SXR,SYR,0.,GXR,GYR,HXR,HYR
      WRITE(iout,1027)'HR',(HVN(J),J=1,NAT3)
      WRITE(iout,1027)'GR',(GVN(J),J=1,NAT3)
 1014 FORMAT(2X,A6,' SX=',F7.4,' SY=',F7.4,' SZ=',F7.4,
     X    /8X, ' GX=',F9.6,' GY=',F9.6,' HX=',F9.6,' HY=',F9.6)
 1027 FORMAT(/3X,A3,'=',3(F15.8,','),/,(6X,3(F15.8,',')))
 1044 FORMAT(/15X,'ROTATION TO ON  G AND H:  BETA=',F12.5,
     X    /5X, ' <G|H> NEW =',E12.5, 
     X    /5X, ' |G|/2 NEW =',E10.5, 5X,  ' |H|   NEW =',E10.5,
     X    /5X, ' GXH/2 NEW =',E10.5)
      return
      end
      SUBROUTINE gtfgh(IPRU,IOUT,EGV,GV,HV,SV,IORG,EREF,
     X                  KILLER,NATOMS)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 LAMBDA(30),LAMBD0(2),LAMBDAR(30),DIPOLCIP(2)
      DIMENSION EGV(30),GV(30),HV(30),GMP(30),DIPOLCI(50),KILLER(3,10)
      DIMENSION GM(30),EGRAD(30),GGRAD(30),HESS(465),XGRAD(30),
     1          X1(3),X2(3),X3(3),SV(30)
      NAMELIST/PROGRS/EGRAD,GGRAD,HESS,NEWDIR,SCALE,GM,EDIF,
     1                LAMBDA,RHSN,ESTATE1,LAMBDAR,XGRAD
      NAMELIST/PRPRTY/DIPOLCI,XGRAD
      rewind ipru
      READ(IPRU,PROGRS,END=100)
      READ(IPRU,PRPRTY,END=102)
  102 CONTINUE
      EREF=ESTATE1
      NAT3=NATOMS*3
      DO 5 I=1,NAT3
      GV(I)=0.0
      HV(I)=0.0
   5  CONTINUE
      DO 19 I=1,3
      X1(I)=0
      X2(I)=0
      X3(I)=0
  19  CONTINUE
      IPT=0
C        I=1 X    I=2 Y   I=3 Z
      DO 20 J=1,NATOMS
      DO 20 I=1,3
C         J=ATOMS
      IF(J.EQ.IORG)GO TO 20
      IF(KILLER(I,J).EQ.1)GO TO 20
      IPT=IPT+1
      X1(I)=X1(I)+EGRAD(IPT)
      X2(I)=X2(I)+GGRAD(IPT)
      X3(I)=X3(I)+XGRAD(IPT)
      IPUT=(J-1)*3+I
      EGV(IPUT)=EGRAD(IPT)
      GV(IPUT)=GGRAD(IPT)
      HV(IPUT)=XGRAD(IPT)
      SV(IPUT)=EGV(IPUT)+GV(IPUT)/2
   20 CONTINUE
      DO 10 I=1,3
      IPUT=(IORG-1)*3+I
      EGV(IPUT)=-X1(I)
      GV(IPUT)=-X2(I)
      HV(IPUT)=-X3(I)
      SV(IPUT)=EGV(IPUT)+GV(IPUT)/2
   10 CONTINUE
      do 11 i=1,nat3
   11 egv(i)=gv(i)/2.
      WRITE(IOUT,1001)
      RETURN
  100 CONTINUE
      WRITE(IOUT,1000)
 1000 FORMAT(5X,'NO PROGRESS FILE')
 1001 FORMAT(5X,'**** EGV GV AND HV FROM PROGRESS FILE ****')
      RETURN
      END
      SUBROUTINE mkrGH(SV,GVS2,HV,GVN,HVN,NAT3,beta0)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION GVS2(2),HV(2),GVN(2),HVN(2),sv(2),scr(100)
      do 10 i=1,nat3
      scr(i)=gvs2(i)
   10 scr(i+nat3)=hv(i)
      CALL SCHMOM(nat3,2,SCR,6,IFAIL)
      HX=DOT(SCR(1),HV,NAT3)
      HY=DOT(SCR(NAT3+1),HV,NAT3)
      Gx=DOT(SCR(1),GVS2,NAT3)
      GY=DOT(SCR(NAT3+1),GVS2,NAT3)
      SX=DOT(SCR(1),SV,NAT3)
      SY=DOT(SCR(NAT3+1),SV,NAT3)
C     
      DENOB=DOT(HV,HV,NAT3)-DOT(GVS2,GVS2,NAT3)
      BETA0=0.25*ATAN(2.*DOT(GVS2,HV,NAT3)/DENOB)
      CALL RBB(NAT3,BETA0,GVS2,HV,GVN,HVN)
 200  continue
      return
      end
      subroutine mkpscent(natoms,maxnew,psc,dpsc,ipsc,mass,iout)
      IMPLICIT REAL*8(A-H,O-Z)
      dimension psc(30,2),dpsc(30,2),ipsc(30,2)
      real*8 mass(2),masst
      write(iout,1001)
      do 20 i=1,maxnew
      masst=0
      do 30 j=1,natoms
   30 psc(j,i)=0.d0
      do 31 j=1,natoms
      jj=ipsc(j,i)
      if(jj.eq.0)go to 31
      if(jj.gt.0)then
      cf=dpsc(j,i)
      else
      cf=mass(iabs(jj))
      masst=masst+cf
      endif
      psc(iabs(jj),i)=cf
   31 continue
      if(masst.lt.0.01)masst=1
      do 32 j=1,natoms
   32 psc(j,i)=psc(j,i)/masst
      write(iout,1000)i,(j,psc(j,i),j=1,natoms)
   20 continue
      return
 1000 format(2x,i2,6(1x,i2,f9.5))
 1001 format(/15x,'PSEUDO CENTERS')
      end
  
      SUBROUTINE MROTV(NATOMS,IONE,V,GM,KZ,KILLER)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER I
      DIMENSION GM(2),V(2),ICROSS(3,3),LOOP(3),KILLER(2)
      DATA ICROSS/0,-3,2, 3,0,-1, 2,-1,0/
      LOOP(1)=KZ
      IX=1
      DO 1 I=1,3
      IF(KZ.EQ.I)GO TO 1
      IX=IX+1
      LOOP(IX)=I
    1 CONTINUE
      IX=0
      DO 20 JJ=1,3
      J=LOOP(JJ)
      IK=0
      DO 10 I=1,NATOMS
      DO 10 K=1,3
      IK=IK+1
      IF(KILLER(IK).EQ.1)GO TO 10
      IX=IX+1
      KK=ICROSS(K,J)
      V(IX)=0.D0
      IF(KK)8,10,9
    9 V(IX)=GM((I-1)*3+KK)-GM((IONE-1)*3+KK)
      GO TO 10
    8 V(IX)=-(GM((I-1)*3-KK)-GM((IONE-1)*3-KK))
   10 CONTINUE
   20 CONTINUE
      RETURN
      END
      SUBROUTINE MVIBV(NATOM,NVIB,NDEGF,V,VR,I76,KILLER)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER I
      DIMENSION VR(2),KILLER(2),JSEC(50),V(2)
      COMMON/VIBORD/IVBORD(30),IVBORD2(30)
      MASK=2*16-1
C
      M76=0
      OPEN(I76,FORM='FORMATTED')
      REWIND I76
      READ(I76,1001,END=500)MDF,NATOMS
      READ(I76,1001)(IZ,I=1,3*NATOMS)
      READ(I76,1001)(IZ,I=1,MDF)
      READ(I76,1001)(IZ,I=1,MDF)
      IF(NATOM.NE.NATOMS)GO TO 900
      M76=1
  500 CONTINUE
      NATOMS=NATOM
      NATOM3=NATOMS*3
C
      NDO=0
      MVIB=0
      IK=0
      DO 10 I=1,NATOMS
C         SYMMETRY EQUIVALENT CENTER INSTRUCTIONS
      IF(M76.NE.0)READ(I76,1001)NDO
      IF(NDO.NE.0)READ(I76,1002)(JSEC(J),J=1,NDO)
C
      DO 10 K=1,3
      IK=IK+1
      IF(KILLER(IK).EQ.1)GO TO 10
      MVIB=MVIB+1
      ISKIP=(MVIB-1)*NATOM3
      IX=(I-1)*3+K + ISKIP
      V(IX)=1.0
C
      IF(NDO.EQ.0)GO TO 50
      DO 40 J=1,NDO
      NCENT=IAND(MASK,JSEC(J))
      PHASE=1.D0
      IF(NCENT.LT.I)PHASE=-1.D0
      IP=IAND(1,ISHFT(JSEC(J),-(16+K-1)))
      IX=(NCENT-1)*3+K+ISKIP
      V(IX)=PHASE
      IF(IP.NE.0)V(IX)=-PHASE
   40 CONTINUE
   50 CONTINUE
C
   10 CONTINUE
C
  100 CONTINUE
      IX=0
      LVIB=0
      DO 110 KVIB=1,NDEGF
      IVIB=IVBORD(KVIB)
      IF(IVIB.LE.0)GO TO 110
      LVIB=LVIB+1
      ISKIP=(IVIB-1)*NATOM3
      IVBORD2(LVIB)=KVIB
      DO 120 IK=1,NATOM3
      IF(KILLER(IK).EQ.1)GO TO 120
      IX=IX+1
      VR(IX)=V(IK+ISKIP)
  120 CONTINUE
      IF(LVIB.EQ.NVIB)GO TO 111
  110 CONTINUE
  111 CONTINUE
      IF(LVIB.NE.NVIB)GO TO 910
      RETURN
  900 CONTINUE
      WRITE(IOUT,1000)I76,NATOM,NATOMS
      STOP ' ERROR IN VIBV '
  910 CONTINUE
      WRITE(IOUT,1003)LVIB,NVIB
      STOP ' ERROR MAKING VIB '
 1000 FORMAT(5X,'UNIT=',I3,' ERROR NATOM=',I3,' NATOMS=',I3)
 1001 FORMAT(16I4)
 1002 FORMAT(8(1X,Z8))
 1003 FORMAT(5X,'ERROR MAKING VIB BASIS: HAVE=',I3,' NEED=',I3)
      END
      SUBROUTINE NUMATM(IUNIT,IP,NATOMS)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION GM(3,10)
      DATA GM/30*123456789./
      NAMELIST/PATH/GM
C
      READ(IUNIT,PATH,END=250)
      DO 252 I=1,10
      IF(DABS(GM(1,I)-123456789.).LT.0.1)GO TO 253
      NATOMS=I
  252 CONTINUE
  253 CONTINUE
      IF(IP.LE.0)GO TO 260
      WRITE(IP,1004)IUNIT
      WRITE(IP,1005)(J,(GM(I,J),I=1,3),J=1,NATOMS)
      GO TO 260
  250 CONTINUE
      WRITE(IP,1003)
      CALL EXIT(100)
  260 CONTINUE
      RETURN
 1003 FORMAT(5X,'**** NO CURRENT GEOMETRY FILE ****' )
 1004 FORMAT(5X,'**** CURRENT GEOMETRY TAKEN FROM UNIT:',I4 )
 1005 FORMAT(5X,'CURRENT GEOMETRY:',/(2X,I3,2X,3F12.6) )
      END      
      SUBROUTINE PUTGRD(IPRU,NATOMS,NDEGF,EDIFC,GMC,EGRADC,GGRADC,
     1             XGRADC,LAMBDAC,LAMBDR,IOUT,GNORMP,ES1P,NFIXDIR)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 LAMBDA(30),LAMBDAC(2),LAMBD0(2),LAMBDAR(30),LAMBDR(2)
      DIMENSION GMC(2),GGRADC(2),HESSC(2),EGRADC(2),XGRADC(2) 
      DIMENSION GM(30),GGRAD(30),HESS(465),EGRAD(30),XGRAD(30)
      DIMENSION GGRADH(930),EGRADH(930),DISP(30),XGRADH(930)
      DATA RHSN/100.D0/,ESTATE/0.0D0/
      NAMELIST/STATUS/ISTATUS
      NAMELIST/PROGRSH/GGRADH,EGRADH,IDISP,GM,EDIF,DISP,LAMBD0,EREF
     X                ,XGRADH
      NAMELIST/PROGRS/GM,EDIF,NEWDIR,EGRAD,GGRAD,XGRAD,SCALE,HESS,
     1                LAMBDA,RHSN,ESTATE1,LAMBDAR
C
      REWIND IPRU
      READ(IPRU,PROGRS,END=103)
  103 CONTINUE
      REWIND IPRU
      READ(IPRU,PROGRSH,END=104)
      CALL WRTPGS(IPRU,NATOMS,NDEGF,GMC,EDIFC,NEWDIR,
     1  EGRADC,GGRADC,XGRADC,SCALE,HESS,LAMBDAC,LAMBDR,
     2  GNORMP,ES1P,NFIXDIR)
      RETURN
  104 CONTINUE
      REWIND IPRU
      READ(IPRU,STATUS)
      CALL WRTPGS(IPRU,NATOMS,NDEGF,GMC,EDIFC,NEWDIR,
     1  EGRADC,GGRADC,XGRADC,SCALE,HESS,LAMBDAC,LAMBDR,
     2  GNORMP,ES1P,NFIXDIR)
      RETURN
      END
C
      SUBROUTINE PUTHES(IPRU,NATOMS,NDEGF,HESSP,IOUT,NFIXDIR)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 LAMBDA(30),LAMBD0(2),LAMBDAR(30)
      DIMENSION HESSP(2)
      DIMENSION EGRADH(930),GGRADH(930),GM(30),EGRAD(30),GGRAD(30),
     1          HESS(465),DISP(30),XGRAD(30),XGRADH(930)
      NAMELIST/PROGRSH/GGRADH,EGRADH,IDISP,GM,EDIF,DISP,LAMBD0,EREF,
     X                 XGRADH
      NAMELIST/PROGRS/EGRAD,GGRAD,HESS,NEWDIR,SCALE,GM,EDIF,
     1                LAMBDA,RHSN,ESTATE1,LAMBDAR,XGRAD
      NAMELIST/STATUS/ISTATUS
      NDEGF2=(NDEGF+1)*NDEGF/2
      REWIND IPRU
      READ(IPRU,PROGRS,END=64)
   64 CONTINUE
      REWIND IPRU
      READ(IPRU,STATUS)
      READ(IPRU,PROGRSH)
      DO 63 I=1,NDEGF2
   63 HESS(I)=HESSP(I)
      CALL WRTPGS(IPRU,NATOMS,NDEGF,GM,EDIF,NEWDIR,
     1  EGRAD,GGRAD,XGRAD,SCALE,HESSP,LAMBDA,LAMBDAR,RHSN,
     2  ES1,NFIXDIR)
      RETURN
      END
      SUBROUTINE RBB(N,BETA,A,B,AN,BN)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(N),B(N),AN(N),BN(N)
      C2B=COS(2*BETA)
      S2B=SIN(2*BETA)
      DO 10 I=1,N
      AN(I)=C2B*A(I)-S2B*B(I)
 10   BN(I)=S2B*A(I)+C2B*B(I)
      RETURN
      END
      SUBROUTINE RBASIS(NR,R,N,SCR,ISCR,BF,IOUT,IPFLG)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION R(N,NR),SCR(2),BF(2),ISCR(2,2)
C         ISCR(1,*) LIST OF AVAILABLE BF
C         ISCR(2,*) CURRENT BF ( NR BF NOT COUNTED )
      DO 2 I=1,N
    2 ISCR(1,I)=I
      LOOP=0
  500 CONTINUE
      NDO=N**2
      DO 5 I=1,NDO
    5 SCR(I)=0.D0
      IJ=0
      IF(NR.EQ.0)GO TO 11
      DO 10 I=1,NR
      DO 10 J=1,N
      IJ=IJ+1
      SCR(IJ)=R(J,I)
   10 CONTINUE
   11 CONTINUE
      IOFF=NR*N
      KOFF=0
      DO 40 I=1,N
      IF(ISCR(1,I).LE.0)GO TO 40
      SCR(ISCR(1,I)+IOFF)=1.D0
      KOFF=KOFF+1
      ISCR(2,KOFF)=I
      IOFF=IOFF+N
      IF(KOFF.EQ.(N-NR))GO TO 41      
   40 CONTINUE 
   41 CONTINUE
      NDO=N-NR
C
      IF(IPFLG.LE.0)GO TO 46
      IHI=0
      WRITE(IOUT,1002)
      DO 45 I=1,N
      LOW=IHI+1
      IHI=IHI+N
      WRITE(IOUT,1000)I,(SCR(K),K=LOW,IHI)
   45 CONTINUE
   46 CONTINUE
C
      CALL SCHMO(N,SCR,IOUT,IFAIL)
      IF(IFAIL.EQ.0)GO TO 47
      IF(IFAIL.LE.NR)GO TO 900
      IPUT=IFAIL-NR
      ISCR(1,ISCR(2,IPUT))=-1
      LOOP=LOOP+1
      WRITE(IOUT,1004)LOOP,ISCR(2,IPUT)
      IF(LOOP.GT.NR)GO TO 900
      GO TO 500
   47 CONTINUE
      IX=0
      IY=NR*N
      DO 50 I=1,NDO
      DO 50 J=1,N
      IX=IX+1
      IY=IY+1
   50 BF(IX)=SCR(IY)
C
      IHI=0
      WRITE(IOUT,1001)
      DO 60 I=1,N-NR
      LOW=IHI+1
      IHI=IHI+N
      WRITE(IOUT,1000)I,(BF(K),K=LOW,IHI)
   60 CONTINUE
      RETURN 
 900  CONTINUE
      WRITE(IOUT,1003)IFAILED
      STOP ' SCHMO FAILURE IN RBASIS '
 1000 FORMAT(2X,I3,(8F10.6))
 1001 FORMAT(' ON - BASIS FOR REDUCED W-MATRIX ')
 1002 FORMAT(' ORIGINAL - BASIS FOR W-MATRIX ')
 1003 FORMAT(5X,'SCHMO FAILURE IN RBASIS:BASIS FUNCTION=',I3)
 1004 FORMAT(5X,'SCHMO FAILURE:LOOP=',I3,' BF=',I3)
      END      
      SUBROUTINE RHESS(NL,NS,TRAN,H,SCR,SCR1,IOUT,IPFLG,NAME)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION TRAN(NL,NS),H(2),SCR(2),SCR1(2)
      CHARACTER*8 NAME
C
      IX=0
      DO 10 I=1,NL
      DO 10 J=1,I
      IX=IX+1
      IJ=(I-1)*NL+J
      SCR(IJ)=H(IX)
      IJ=(J-1)*NL+I
      SCR(IJ)=H(IX)
   10 CONTINUE
      IF(IPFLG.GT.0)THEN
      WRITE(IOUT,1000)'ORIGINAL',NAME
      CALL PRVOM(H,NL,2,IOUT)
      ENDIF
C
      CALL EBC(SCR1,SCR,TRAN,NL,NL,NS)
      CALL EBTC(SCR,TRAN,SCR1,NS,NL,NS)
C
      IX=0
      DO 20 I=1,NS
      DO 20 J=1,I
      IX=IX+1
      IJ=(I-1)*NS+J
      H(IX)=SCR(IJ)
      JI=(J-1)*NS+I
      IF(DABS(SCR(JI)-SCR(IJ)).GT.1.0D-08)GO TO 900
   20 CONTINUE
C
      IF(IPFLG.GT.0)THEN
      WRITE(IOUT,1000)'REDUCED',NAME
      CALL PRVOM(H,NS,2,IOUT)
      ENDIF
C
      RETURN
  900 CONTINUE
      STOP ' ERROR IN RHESS '
 1000 FORMAT(3X,A8,1X,A8)
      END
      SUBROUTINE RESTLAG(N,A,B,RLAG,IOUT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(N),B(N)
      AB=0.D0
      BB=0.D0
      DO 10 I=1,N
      AB=AB+A(I)*B(I)
      BB=BB+B(I)*B(I)
   10 CONTINUE
      RLAG=-AB/BB
      WRITE(IOUT,1000)RLAG
      RETURN
 1000 FORMAT(5X,'LAGRANG MULT=',E15.8,' FROM RESTLAG')
      END
      SUBROUTINE RDT30(N30,E,NSTATE,IOUTU,IERR)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION E(2)
      IERR=0
      REWIND N30
      READ(N30,END=60)N,NSTATE,E0
   10 READ(N30,END=20)IT,IUNIT,(E(I),I=1,NSTATE)
      GO TO 10
   20 CONTINUE
C      DO 21 I=1,NSTATE
C   21 E(I)=E(I)+E0
      WRITE(IOUTU,1001)(E(I),I=1,NSTATE)
      RETURN
   60 CONTINUE
      WRITE(IOUTU,1000)N30
      IERR=1
      RETURN
 1000 FORMAT(5X,'UNIT ',I3,' EMPTY ')
 1001 FORMAT(5X,'CI ENERGIES FROM ALCHEMY TAPE',(/2X,5E20.13))
      END
      SUBROUTINE rptcon(NFIXDIR,FIXDIR,VALDIR,NDEGF,NATOMS,
     1                  GM,KILLER,KLASS,IOUT)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER FIXDIR(6,2),KLASS(2)
      DIMENSION KILLER(2),VALDIR(2)
      DIMENSION GM(3,2)
      common/pscent/maxnew,psc(30,10),dpsc(30,10),ipsc(30,10)
      DATA TWO/2.D0/,ZERO/0.D0/,PI/3.141592654D0/
      if(maxnew.eq.0)go to 99
      do 40 i=1,maxnew
      do 41 k=1,3
   41 gm(k,natoms+i)=0
      do 42 k=1,3
      do 42 j=1,natoms
   42 gm(k,natoms+i)=gm(k,natoms+i)+gm(k,j)*psc(j,i)
   40 continue
   99 continue
C
      IF(NFIXDIR.EQ.0)GO TO 900
      NF=NDEGF+NFIXDIR
      N3=NATOMS*3
C
      MYCLASS=0
      ZNORM=ZERO
      DO 200 I=1,NFIXDIR
      II=FIXDIR(1,I)
      JJ=FIXDIR(2,I)
      GO TO (400,600,400,600,700),KLASS(I)
      STOP ' ILLEGAL KLASS '
  400 CONTINUE
C         DISTANCE CONTRAINT
      MYCLASS=1
      X = (GM(1,II)-GM(1,JJ))**2 + (GM(2,II)-GM(2,JJ))**2+
     1         (GM(3,II)-GM(3,JJ))**2 
      IF(KLASS(I).EQ.3)GO TO 150
      WRITE(IOUT,1028)II,JJ,DSQRT(X),VALDIR(I)
      ZNORM=ZNORM+(X-VALDIR(I)**2)**2
      GO TO 200
  150 CONTINUE
C         DISTANCE DIFFERENCE CONSTRAINT
      kk=FIXDIR(3,I)
      LL=FIXDIR(4,I)
      y = (GM(1,kk)-GM(1,LL))**2 + (GM(2,kk)-GM(2,LL))**2+
     1         (GM(3,kk)-GM(3,LL))**2 
      WRITE(IOUT,1027)II,JJ,DSQRT(X),KK,LL,DSQRT(Y)
      GO TO 200
  700 CONTINUE
C         dihedral ANGLE CONSTRAINT
      KK=FIXDIR(3,I)
      LL=FIXDIR(6,I)
      KKU=(KK-1)*3
      LLU=(LL-1)*3
      CALL DIHED1(GM(1,II),GM(1,JJ),GM(1,KK),GM(1,LL),DIH1)
      WRITE(IOUT,1030)II,JJ,KK,LL,acos(DIH1)*180.d0/pi,VALDIR(I)
      XX=VALDIR(I)*PI/180.D0
      GO TO 200
  600 CONTINUE
C         ANGLE CONSTRAINT
      KK=FIXDIR(3,I)
      CALL ANGL(GM,II,JJ,KK,ANG,CA,R13,R12,R23)
      IF(KLASS(I).EQ.4)GO TO 170
      WRITE(IOUT,1029)II,JJ,KK,ANG,VALDIR(I)
      XX=VALDIR(I)*PI/180.D0
      GO TO 200
  170 CONTINUE
C         ANGLE DIFFERENCE
      IIO=II
      JJO=JJ
      KKO=KK
      II=FIXDIR(4,I)
      JJ=FIXDIR(5,I)
      KK=FIXDIR(6,I)
      CALL ANGL(GM,II,JJ,KK,ANGB,CB,R13,R12,R23)
      CAT=CA
      CA=CB
      WRITE(IOUT,1026)IIO,JJO,KKO,ANG,II,JJ,KK,ANGB
  200 CONTINUE
      ZNORM=DSQRT(ZNORM)
      IF(MYCLASS.EQ.1)WRITE(IOUT,1001)ZNORM
C
  900 CONTINUE
      RETURN
 1000 FORMAT(5X,'GRADIENT ',A7,' CONSTRAINTS:DIMENSION=',I3,
     X     (/5X,8F12.6))
 1001 FORMAT(/5X,'DISTANCE NORM=',E15.8)
 1026 FORMAT(5X,'ANG(',I2,',',I2,',',I2,')  ='
     X       ,F12.6,' ANG(',I2,',',I2,',',I2,')  =',F12.6)
 1027 FORMAT(5X,'|R(',I2,') - R(',I2,')|=',F12.6,
     X ' |R(',I2,') - R(',I2,')|=',F12.6)
 1028 FORMAT(5X,'|R(',I2,') - R(',I2,')|=',F12.6,
     X ' CONSTRAINT VALUE=',F12.6)
 1029 FORMAT(5X,'ANG(',I2,',',I2,
     X       ',',I2,')  =',F12.6,' CONSTRAINT VALUE=',F12.6)
 1030 FORMAT(5X,'DIH ANG(',I2,',',I2,
     X       ',',I2,',',i2,')  =',F12.6,' CONSTRAINT VALUE=',F12.6)
      END
      SUBROUTINE SCHMO(N,V,IOUT,IFAIL)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION V(N,N)
      IFAIL=0
      DO 100 I=1,N
      Z=0.D0
      DO 90 J=1,N
   90 Z=Z+V(J,I)**2 
      Z=DSQRT(Z)
      DO 80 J=1,N
   80 V(J,I)=V(J,I)/Z
  100 CONTINUE
C
      DO 500 I=2,N
      DO 300 J=1,I-1
      DP=DOT(V(1,I),V(1,J),N)
      DO 250 K=1,N
  250 V(K,I)=V(K,I)-DP*V(K,J)
  300 CONTINUE
      Z=0.D0
      DO 260 K=1,N
  260 Z=Z+V(K,I)**2
      Z=DSQRT(Z)
      IF(Z.LT.1.0D-06)THEN
      WRITE(IOUT,1000)I,Z
      IFAIL=I
      RETURN
      ENDIF
      DO 270 K=1,N
  270 V(K,I)=V(K,I)/Z
  500 CONTINUE
      RETURN
 1000 FORMAT(5X,'BAD VECTOR IN SCHMO:I=',I3,' NORM=',E15.6)
      END
      SUBROUTINE SCHMOM(N,M,V,IOUT,IFAIL)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION V(N,M)
      IFAIL=0
      DO 100 I=1,M
      Z=0.D0
      DO 90 J=1,N
   90 Z=Z+V(J,I)**2 
      Z=DSQRT(Z)
      DO 80 J=1,N
   80 V(J,I)=V(J,I)/Z
  100 CONTINUE
C
      DO 500 I=2,M
      DO 300 J=1,I-1
      DP=DOT(V(1,I),V(1,J),N)
      DO 250 K=1,N
  250 V(K,I)=V(K,I)-DP*V(K,J)
  300 CONTINUE
      Z=0.D0
      DO 260 K=1,N
  260 Z=Z+V(K,I)**2
      Z=DSQRT(Z)
      IF(Z.LT.1.0D-06)THEN
      WRITE(IOUT,1000)I,Z
      IFAIL=I
      RETURN
      ENDIF
      DO 270 K=1,N
  270 V(K,I)=V(K,I)/Z
  500 CONTINUE
      RETURN
 1000 FORMAT(5X,'BAD VECTOR IN SCHMO:I=',I3,' NORM=',E15.6)
      END
      SUBROUTINE SOLNR(N,RHS,VEC,EVAL,TOL,DAMP,SOL,SCR,ITP19,
     X                 ISADDLE,NFIXDIR)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL FLAG
      DIMENSION RHS(N),VEC(N,N),EVAL(N),DAMP(N),SOL(N),SCR(N)
      DATA ZERO/0.D0/
      FLAG=.FALSE.
C
      WRITE(ITP19,1001)
      RHSNRM1=0.D0
      DO 10 I=1,N
   10 RHSNRM1=RHSNRM1+RHS(I)**2
      RHSNRM1=DSQRT(RHSNRM1)
      DO 100 J=1,N
      DP=ZERO
      CNORM=ZERO
      PNORM=ZERO
      DO 30 K=1,N
      CNORM=CNORM+RHS(K)**2
   30 DP=DP+RHS(K)*VEC(K,J)
      SCR(J)=DP
      IF(DABS(EVAL(J)).GT.TOL)GO TO 100
      FLAG=.TRUE.
CDRY      DMPF=1.D0-DABS(EVAL(J)/DAMP(J))
CDRY 7/95
      DMPF=DAMP(J)
      WRITE(ITP19,1004)J,EVAL(J),TOL,DMPF
      DO 40 K=1,N
      RHS(K)=RHS(K)-DP*DMPF*VEC(K,J)
   40 PNORM=PNORM+RHS(K)**2
      CNORM=DSQRT(CNORM)
      PNORM=DSQRT(PNORM)
CDRY 7/20/98
      DP=0.0
      DO 31 K=1,N
   31 DP=DP+RHS(K)*VEC(K,J)
      SCR(J)=DP
C
      WRITE(ITP19,1002)J,DP,DMPF,CNORM,PNORM,SCR(J)
  100 CONTINUE
      IF(FLAG)WRITE(ITP19,1000)RHS
      WRITE(ITP19,1003)(SCR(I),I=1,N)
      RHSNRM2=0.D0
      DO 11 I=1,N
   11 RHSNRM2=RHSNRM2+SCR(I)**2
      RHSNRM2=DSQRT(RHSNRM2)
      WRITE(ITP19,1005)RHSNRM1,RHSNRM2
C
      DO 210 J=1,N
      X=ZERO
      DO 200 I=1,N
CDRY      C=SCR(I)/(EVAL(I)+DAMP(I))
      C=SCR(I)/EVAL(I)
      IF(DABS(EVAL(I)).LT.TOL)C=SCR(I)/(EVAL(I)+DAMP(I))
      IF((I.GT.NFIXDIR).AND.(ISADDLE.EQ.0).AND.(EVAL(I).LT.0.D0))
     1          C=SCR(I)
      X=X+VEC(J,I)*C
  200 CONTINUE
      SOL(J)=X
  210 CONTINUE
C
      RETURN
 1000 FORMAT(5X,'PROJECTED GRADIENT:',4E15.8,(/2X,5E15.8))
 1001 FORMAT(/15X,' GRADIENT PROJECTION ROUTINE',/)
 1002 FORMAT(5X,'< ',I1,' ! GRAD >=',E15.8,' DAMPING=',E15.8
     1,/5X,' NORM BEFORE(AFTER) PROJECTION=',E15.8,'(',E15.8,')',
     2 ' PROJECTED RHS:',E15.8) 
 1003 FORMAT(5X,'GRADIENTS IN EIGENSTATE BASIS ',4E15.8,(/2X,5E15.8))
 1004 FORMAT(5X,'DAMPING FOR ROOT=',I2,' EVAL=',E11.4,' TOL=',E11.4,
     1          'DAMP FACT=',E11.4)
 1005 FORMAT(5X,'FROM SOLNR: NORM RHS- OB=',F15.8,' ESB=',F15.8)
      END
      SUBROUTINE SGNHIJ(IPRU,NDEGF,HIJP,IOUT,CHANGE)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 LAMBDA(30),LAMBD0(2),LAMBDAR(30)
      DIMENSION HIJP(2)
      DIMENSION GM(30),EGRAD(30),GGRAD(30),HESS(465),XGRAD(30)
      NAMELIST/PROGRS/EGRAD,GGRAD,HESS,NEWDIR,SCALE,GM,EDIF,
     1                LAMBDA,RHSN,ESTATE1,LAMBDAR,XGRAD
      NAMELIST/STATUS/ISTATUS
C
      CHANGE=1.D0
      DOT=0.D0
      RNORM1=0.D0
      RNORM2=0.D0
      REWIND IPRU
      READ(IPRU,PROGRS,END=800)
      DO 100 I=1,NDEGF
      RNORM1=RNORM1+XGRAD(I)**2
      RNORM2=RNORM2+HIJP(I)**2
      DOT=DOT+HIJP(I)*XGRAD(I)
  100 CONTINUE
      RNORM1=DSQRT(RNORM1)
      RNORM2=DSQRT(RNORM2)
      IF(DOT.GE.0.0D0)GO TO 800
      WRITE(IOUT,1000)NDEGF,DOT,RNORM1,RNORM2
      CHANGE=-1.D0
      DO 101 I=1,NDEGF
  101 HIJP(I)=-HIJP(I)
  800 CONTINUE
      REWIND IPRU
      RETURN
 1000 FORMAT(5X,'CHANGING SIGN OF HIJ: LENGTH=',I3,
     X      /5X,'OVERLAP=',E12.6,' NORMS=',2E12.6)
      END
      SUBROUTINE SECGM(NATOM,GM,I76,IOUT,NATUN)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION JSEC(50),GM(3,2),GML(3,30)
      MASK=2*16-1
C
      OPEN(I76,FORM='FORMATTED')
      REWIND I76
      READ(I76,1001,END=500)MDF,NATOMS
      READ(I76,1001)(IZ,I=1,3*NATOMS)
      READ(I76,1001)(IZ,I=1,MDF)
      READ(I76,1001)(IZ,I=1,MDF)
      IF(NATOM.NE.NATOMS)GO TO 900
C
C         SYMMETRY EQUIVALENT CENTER INSTRUCTIONS
      NATUN=0
      KSEC=0
      DO 310 I=1,NATOMS
      READ(I76,1001)NDO
      KSEC=KSEC+NDO
      IF(NDO.EQ.0)THEN
      NATUN=NATUN+1
      DO 311 J=1,3
  311 GML(J,NATUN)=GM(J,I)
      GO TO 310
      ENDIF
      READ(I76,1002)(JSEC(J),J=1,NDO)
      DO 320 J=1,NDO
      NCENT=IAND(MASK,JSEC(J))
      IF(NCENT.LT.I)GO TO 310
  320 CONTINUE
      NATUN=NATUN+1
      DO 321 J=1,3
  321 GML(J,NATUN)=GM(J,I)
  310 CONTINUE
      IF(KSEC.EQ.0)GO TO 500
      DO 350 I=1,NATUN
      DO 350 J=1,3
  350 GM(J,I)=GML(J,I)
      WRITE(IOUT,1003)NATOMS,NATUN,((GM(K,J),K=1,3),J=1,NATUN)
C
  500 CONTINUE
      RETURN
C
  900 CONTINUE
      WRITE(IOUT,1000)I76,NATOM,NATOMS
 1000 FORMAT(5X,'UNIT=',I3,' ERROR NATOM=',I3,' NATOMS=',I3)
 1001 FORMAT(16I4)
 1002 FORMAT(8(1X,Z8))
 1003 FORMAT(5X,'UNIQ ATOM GEOM: NATOMS NATUN',2I3,/(2X,3F16.9))
      END
      SUBROUTINE SECNT(NATOM,GRAD,I76,IOUT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION JSEC(50),GRAD(3,2)
      MASK=2*16-1
C
      OPEN(I76,FORM='FORMATTED')
      REWIND I76
      READ(I76,1001,END=500)MDF,NATOMS
      READ(I76,1001)(IZ,I=1,3*NATOMS)
      READ(I76,1001)(IZ,I=1,MDF)
      READ(I76,1001)(IZ,I=1,MDF)
      IF(NATOM.NE.NATOMS)GO TO 900
C
      WRITE(IOUT,1003)'BEFORE',((GRAD(K,J),K=1,3),J=1,NATOMS)
C         SYMMETRY EQUIVALENT CENTER INSTRUCTIONS
      KSEC=0
      DO 310 I=1,NATOMS
      READ(I76,1001)NDO
      KSEC=KSEC+NDO
      IF(NDO.EQ.0)GO TO 310
      READ(I76,1002)(JSEC(J),J=1,NDO)
      DO 320 J=1,NDO
      NCENT=IAND(MASK,JSEC(J))
      IF(NCENT.LT.I)GO TO 320
      DO 330 K=1,3
      IP=IAND(1,ISHFT(JSEC(J),-(16+K-1)))
      GP=GRAD(K,I)
      IF(IP.NE.0)GP=-GP
      GRAD(K,NCENT)=GP
  330 CONTINUE
  320 CONTINUE
  310 CONTINUE
      IF(KSEC.GT.0)WRITE(IOUT,1003)'AFTER',
     1                  ((GRAD(K,J),K=1,3),J=1,NATOMS)
C
  500 CONTINUE
      RETURN
C
  900 CONTINUE
      WRITE(IOUT,1000)I76,NATOM,NATOMS
 1000 FORMAT(5X,'UNIT=',I3,' ERROR NATOM=',I3,' NATOMS=',I3)
 1001 FORMAT(16I4)
 1002 FORMAT(8(1X,Z8))
 1003 FORMAT(5X,'GRADIENTS FROM SECNT: ',A6,/(2X,3F16.9))
      END
      SUBROUTINE TBAK(NL,NS,T,X,Y,IGO)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION T(NL,NS),X(2),Y(2)
      DO 10 I=1,NL
   10 Y(I)=0.D0
      DO 20 I=1,NS
      DO 20 J=1,NL
      Y(J)=Y(J)+T(J,I)*X(I)
   20 CONTINUE
      IF(IGO.EQ.1)RETURN
      DO 30 I=1,NL
   30 X(I)=Y(I)
      RETURN
      END      
      SUBROUTINE TOINTRL(MASS,GM,NACME,GRAD,IOUT,NAME)
      IMPLICIT REAL*8(A-H,O-Z)
C       NACME(XYZ,ATOM)
      REAL*8 MASS(3),NACME(3,3),GM(3,3),CCORDI(2,2),PCORDI(2,2)
     1      ,NACINTC(2,3),NACINTP(2,2),COMD(2),BI(3,3),GRAD(3)
      EQUIVALENCE(N1,NDIAT(1)),(N2,NDIAT(2)),(N3,NDIAT(3))
      CHARACTER*3 NAME
      COMMON/DIATAT/NDIAT(3)
C
      PI=DACOS (-1.D0)
      RTD=180.D0/PI
      WRITE(IOUT,1003)NAME,NACME
      TMASSD=MASS(N2)+MASS(N3)
C         RELATIVE TO COM OF TRIATOM
C         BUILD BI MATRIX  TO DETERMINE INTERNAL NACME'S
      CALL BCOM(MASS,BI)
      DO 110 I=1,2
      COMD(I)=( MASS(N2)*GM(I,N2) + MASS(N3)*GM(I,N3) )/TMASSD
      CCORDI(I,2)=GM(I,N3)-GM(I,N2)
      CCORDI(I,1)=GM(I,N1)-COMD(I)
      X=0.D0
      Y=0.D0
      Z=0.D0
      DO 115 J=1,3
      X=X + NACME(I,J)*BI(J,1)
      Y=Y + NACME(I,J)*BI(J,2)
      Z=Z + NACME(I,J)*BI(J,3)
  115 CONTINUE
      NACINTC(I,1)=X
      NACINTC(I,2)=Y
      NACINTC(I,3)=Z
  110 CONTINUE
C      WRITE(IOUT,1009)((NACINTC(I,J),J=1,3),I=1,2)
      DO 120 I=1,2
      PCORDI(1,I)=DSQRT( CCORDI(1,I)**2 + CCORDI(2,I)**2 )
      X=DACOS ( CCORDI(1,I)/PCORDI(1,I) )
      IF(CCORDI(2,I).LT.0.D0)X=2.D0*PI-X
      PCORDI(2,I)=X*RTD
      NACINTP(1,I)=(CCORDI(1,I)*NACINTC(1,I)+CCORDI(2,I)*NACINTC(2,I))
     1             /PCORDI(1,I)
      NACINTP(2,I)=CCORDI(1,I)*NACINTC(2,I)-CCORDI(2,I)*NACINTC(1,I)
  120 CONTINUE
      WRITE(IOUT,1008)PCORDI(1,1),PCORDI(1,2),
     X               PCORDI(2,1)-PCORDI(2,2)
C      WRITE(IOUT,1005)NACINTP
      DBR=NACINTP(1,1)
      DLR=NACINTP(1,2)
      DGAMMA=(NACINTP(2,1)-NACINTP(2,2) )/2.D0
      DG =(NACINTP(2,2)-NACINTP(2,1))/( 2.D0*PCORDI(1,1) )
      DGAMMA=(NACINTP(2,1)-NACINTP(2,2) )/2.D0
      DOMEGA=(NACINTP(2,2)+NACINTP(2,1))/2.D0
C      WRITE(IOUT,1013)' TRIATOM'
      GRAD(1)=DLR
      GRAD(2)=DBR
      GRAD(3)=DGAMMA
      WRITE(IOUT,1004)NAME,GRAD
      WRITE(IOUT,1014)NAME,DOMEGA,NACINTC(1,3),NACINTC(2,3)
C      WRITE(IOUT,1011)DG,PCORDI(1,1)
      RETURN
 1003 FORMAT(/5X,'CARTESIAN ',A3,/7X,'X',14X,'Y',/(1X,3F15.9) )
C 1004 FORMAT(/5X,A3,' IN JACOBI r  R gma:', 3(1X,,E16.9))
 1004 FORMAT(/5X,A3,' IN JACOBI r  R gma:', 3(1X,E16.9))
 1005 FORMAT(/5X,'  POLAR NACMES',/7X,'R',12X,'THETA',
     2       (2(/1X,2E16.9)) )
 1006 FORMAT(/5X,' INTERNAL POLAR COORDINATES',/7X,'R',12X,'THETA',
     2       (2(/1X,2E16.9)) )
C1007 FORMAT(' BI MATRIX:',/(2X,3F12.6) )
 1008 format(/,' R  r  gamma  ',3F15.9)
 1009 FORMAT(/,' INTERNAL CARTESIAN NACMES',/' X',3E16.9,/' Y',3E16.9)
 1011 FORMAT(' REVISED VERSION: - 1/R DGM=',E16.9,' R =',F12.6)
 1013 FORMAT(/5X,' FIXED ORIGIN IS ',A8 )
 1014 FORMAT(5X,A3,' TR-ROT: ANG=',F12.6,' X3=',F12.6,
     1           ' Y3= ',F12.6 )
      END
      SUBROUTINE WRTPGS(IPRU,NATOMS,NDEGF,GM,EDIF,NEWDIR,
     1                  EGRAD,GGRAD,XGRAD,SCALE,HESS,
     2                  LAMBDA,LAMBDAR,RHSN,ESTATE1,NFIXDIR)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION GM(2),EGRAD(2),GGRAD(2),XGRAD(2),HESS(2),
     1          GGRADH(2),EGRADH(2),XGRADH(2),DISP(2)
      REAL*8 LAMBDA(2),LAMBD0(2),LAMBDAR(2)
      NDEGF2=NDEGF*(NDEGF+1)/2
      WRITE(IPRU,1000)'PROGRS  '
      WRITE(IPRU,1002)'      GM',(GM(I),I=1,3*NATOMS)
      WRITE(IPRU,1002)'   EGRAD',(EGRAD(I),I=1,NDEGF)
      WRITE(IPRU,1002)'   GGRAD',(GGRAD(I),I=1,NDEGF)
      WRITE(IPRU,1002)'   XGRAD',(XGRAD(I),I=1,NDEGF)
      WRITE(IPRU,1002)'    HESS',(HESS(I),I=1,NDEGF2)
      WRITE(IPRU,1002)'  LAMBDA',(LAMBDA(I),I=1,NFIXDIR)
      WRITE(IPRU,1002)' LAMBDAR',(LAMBDAR(I),I=1,NFIXDIR)
      WRITE(IPRU,1002)'   SCALE',SCALE
      WRITE(IPRU,1002)'    EDIF',EDIF
      WRITE(IPRU,1002)'    RHSN',RHSN
      WRITE(IPRU,1002)' ESTATE1',ESTATE1
      WRITE(IPRU,1003)'  NEWDIR',NEWDIR
      WRITE(IPRU,1001)
      RETURN
      ENTRY WRTPGSH(IPRU,NATOMS,NDEGF,EGRADH,GGRADH,XGRADH,GM,DISP,
     1              EDIF,IDISP,LAMBD0,EREF)
      NWRIT=(2*NDEGF+1)*NDEGF
      WRITE(IPRU,1000)'PROGRSH '
      WRITE(IPRU,1002)'      GM',(GM(I),I=1,3*NATOMS)
      WRITE(IPRU,1002)'  EGRADH',(EGRADH(I),I=1,NWRIT)
      WRITE(IPRU,1002)'  GGRADH',(GGRADH(I),I=1,NWRIT)
      WRITE(IPRU,1002)'  XGRADH',(XGRADH(I),I=1,NWRIT)
      WRITE(IPRU,1002)'    DISP',(DISP(I),I=1,NDEGF)
      WRITE(IPRU,1002)'  LAMBD0',(LAMBD0(I),I=1,2)
      WRITE(IPRU,1002)'    EDIF',EDIF
      WRITE(IPRU,1002)'    EREF',EREF
      WRITE(IPRU,1003)'  IDISP',IDISP
      WRITE(IPRU,1001)
      RETURN
 1000 FORMAT(1X,'&',A8)
 1001 FORMAT(1X,'&END')
 1002 FORMAT(1X,A8,'=',5(E20.10,',')/(11X,(5(E20.10,','))))
 1003 FORMAT(1X,A8,'=',15(I5,',')/(11X,(15(I5,','))))
      END
c---------------------------------------------------------------------
c---------------------------------------------------------------------
      subroutine wrtcolgm( gmunit, natoms )

c     Subtoutine to write out geometry in $unit in columbus format
      implicit none
      integer             gmunit, natoms
      integer             i, ios
c     namelis for reading in geometry
      real*8              gm( 30 )
      namelist /path/     gm
c     i73 is unit of dependent input file (unit=73)
      integer             i73
      parameter           ( i73 = 73 )
c     ifnlgm is unit of final geometry file (unit=79)
      integer             ifnlgm
      parameter           ( ifnlgm = 79 )
c     namelist for printout
      real*8              masses( 10 ), numbers( 10 )
      character*1         names ( 10 )
      namelist /printout/ names, masses, numbers
      
c     open unit 73 and read printout namelist
      open ( unit = i73, status = "old", action = "read", iostat = ios )
      if ( ios .ne. 0 ) then
         write( 6, "(1x,a)" ) "Could not open unit 73."
         return
      end if
      read ( unit = i73, nml = printout )
c     close file
      close ( unit = i73 )

c     open unit $unit and read geometry in file
      open ( unit = gmunit, status = "old", action = "read",
     &       iostat = ios )
      if ( ios .ne. 0 ) then
         write ( 6, "(1x,a)" ) "Could not open old geometry file."
         return
      end if
      read ( unit = gmunit, nml = path )
c     close file
      close ( unit = gmunit )

c     print final geometry
      open ( unit = ifnlgm, status = "new", action = "write", 
     &       name = "geom.final", iostat = ios )
      if ( ios .ne. 0 ) then
         write( 6, "(1x,a)" ) "Could not open new geometry file."
         return
      end if
      do 10 i = 1, natoms
         write( ifnlgm, 1000 ) names(i), numbers(i),
     &   gm( (i-1)*3 + 1 ), gm( (i-1)*3 + 2 ), gm( (i-1)*3 + 3 )
 10   continue
c     close file
      close ( ifnlgm )

      return

c     Columbus geometry format
 1000 format(1x,a2,2x,f5.1,4f14.8)

      end subroutine wrtcolgm
c---------------------------------------------------------------
c---------------------------------------------------------------
      subroutine setupjob( workdir )
c     Set up polyhes calculation for surfgen job
c     Creates the following directories:
c       ./Polyhes, ./Polyhes/ProgDir
c     Returns: workdir
      implicit none
      character*100       jobdir,  workdir, progdir
      character*100       inputfl,  progfl,  geomfl 
      character*100       mkwkcmd, mkpgcmd
      character*100       enterwk
      integer             ios
      integer             i

c     get current directory
      call getcwd( jobdir )
      write( 6,* ) " Current directory: ", trim( jobdir )

c     make work directory
      workdir = trim( jobdir ) // '/Polyhes'
      write ( mkwkcmd, "(2a)" ) "mkdir ", workdir
      call system ( mkwkcmd, ios )
      if ( ios .ne. 0 ) stop "Could not make work directory."
  
c     make progress file directory
      progdir = trim( workdir ) // '/ProgDir'
      write ( mkpgcmd, "(2a)" ) "mkdir ", progdir
      call system ( mkpgcmd, ios )
      if ( ios .ne. 0 ) stop "Could not make progress directory."

      return
      
      end subroutine setupjob
c-------------------------------------------------------------
c-------------------------------------------------------------
      subroutine setupfiles( )
c     Set up files for polyhes run.
c     generates: fort.32 (progress file)
c                fort.74 (input geometry file)
c                fort.73 (input file)
      implicit none
      return
      end subroutine setupfiles
