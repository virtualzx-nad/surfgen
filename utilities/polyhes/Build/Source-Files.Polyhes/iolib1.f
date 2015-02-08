      SUBROUTINE ASREAD(FILE,BUFFER,NWORDS)
C
      IMPLICIT INTEGER(A-Z)
C
      DIMENSION BUFFER(NWORDS)
C
      CALL RREAD(FILE,BUFFER,NWORDS)
C
      RETURN
      END
      SUBROUTINE ASWRIT(FILE,BUFFER,NWORDS)
C
      IMPLICIT INTEGER(A-Z)
C
      DIMENSION BUFFER(NWORDS)
C
      CALL RWRIT(FILE,BUFFER,NWORDS)
C
      RETURN
      END
      SUBROUTINE AWAIT(FILE)
C
      IMPLICIT INTEGER(A-Z)
C
C       A DO NOTHING ROUTINE
      I=5
      RETURN
      END
      SUBROUTINE GETWA(FILE,ADDR)
C
      IMPLICIT INTEGER(A-Z)
C
      COMMON /POINTR/ WPTR(128),TPTR(128)
C
C     FETCH THE ADDRESS OF THE WORD POINTER
      ADDR=WPTR(FILE)
C
      RETURN
      END
      SUBROUTINE SETWA(FILE,ADDR)
C
      IMPLICIT INTEGER(A-Z)
C
      COMMON /POINTR/ WPTR(128),TPTR(128)
C
C     SET THE WORD POINTER
      WPTR(FILE)=ADDR
C
      RETURN
      END
       subroutine rfile(file)
cdry modified 6/21/93
c      SUBROUTINE SFILE(FILE,FCB)
C
      IMPLICIT INTEGER(A-Z)
      character*1 nc(3)
      character*8 fyle
      data fyle/'fort.   '/
C
      COMMON /POINTR/ WPTR(128),TPTR(128)
C
      DIMENSION FCB(16)
C
      WPTR(FILE)=1
      TPTR(FILE)=1
C
C     -------- OPEN THE FILE -------
C
      call tochar(file,nc)
      fyle(6:6)=nc(1)
      fyle(7:7)=nc(2)
      fyle(8:8)=nc(3)
C      OPEN  (UNIT=FILE,file='fyle',
      OPEN  (UNIT=FILE,
     #       ACCESS='DIRECT',
     #       RECL=4*1024,
     #       IOSTAT=IERR)
C
      IF(IERR.NE.0) THEN
         WRITE(6,1) IERR,FILE,fyle
    1    FORMAT(//10X,'SFILE, ERROR CODE ',I5,' ENCOUNTERED ',
     #          'OPENING FILE ',I3,2x,a8)
         CALL PERROR(' ERROR rfile ' )
         CALL EXIT
      ELSE
CDRY         WRITE(6,2) FILE
    2    FORMAT(//10X,'FILE',I4,' IS OPEN')
      ENDIF
C
      RETURN
      END
       subroutine rread(file,buffer,nwords,iadr)
c      SUBROUTINE SREAD(FILE,BUFFER,NWORDS)
C
      IMPLICIT INTEGER(A-Z)
C
      COMMON /POINTR/ WPTR(128),TPTR(128)
C
      DIMENSION BUFFER(NWORDS)
C
CDBG
C      WRITE(6,1000) FILE,NWORDS,IADR
C1000 FORMAT(' RREAD:    FILE',I3,I12,' WORDS','  IADR ',i12)
CDBG
C
C     DETERMINE FWORD AND THEN CALL WREAD
C
c      FWORD=WPTR(FILE)
c
      reclen=1024
      fword=(iadr-1)*reclen+1
      CALL WREADW(FILE,BUFFER,NWORDS,FWORD,END)
      icheck=mround(end-1)+1
      wptr(file)=icheck
C
      RETURN
      END
      subroutine sread(file,buffer,nwords)
      implicit integer(a-z)
      common/pointr/wptr(128),tptr(128)
      dimension buffer(nwords)
      fword = wptr(file)
       icheck=mround(fword-1)+1
        if(icheck.ne.fword) then
         write(6,*)' sread: fword is not on a sector boundary round-up '
         fword=icheck
        endif
c
      call wreadw(file,buffer,nwords,fword,end)
c
      icheck=mround(end-1)+1
      wptr(file)=icheck
c
      return
      end
C
      SUBROUTINE SREW(FILE)
C
      IMPLICIT INTEGER(A-Z)
C
      COMMON /POINTR/ WPTR(128),TPTR(128)
C
C     --- SET THE WORD POINTER BACK TO 1 ---
C
      WPTR(FILE)=1
C
      RETURN
      END
        subroutine rwfil(file,i,j,k,l,m)
c     dummy IBM rewind file
      return
      end
      SUBROUTINE SWEOF(FILE)
C
      IMPLICIT INTEGER(A-Z)
C     PRETEND YOU'RE WRITING AN
C     END OF FILE MARKER
      RETURN
      END
      subroutine rwrit(file,buffer,nwords,iadr)
c      SUBROUTINE SWRIT(FILE,BUFFER,NWORDS)
C
      IMPLICIT INTEGER(A-Z)
C
      COMMON /POINTR/ WPTR(128),TPTR(128)
C
      DIMENSION BUFFER(NWORDS)
C
C     DETERMINE FWORD AND THEN CALL SWRIT
C
c      FWORD=WPTR(FILE)
c
      reclen=1024
      fword=(iadr-1)*reclen+1
CDBG
C     WRITE(6,12) FILE,NWORDS,IADR
C  12 FORMAT(' RWRIT:  FILE',I3,I12,' WORDS','  IADR ',i12)
CDBG
      CALL WWRITW(FILE,BUFFER,NWORDS,FWORD,END)
      icheck=mround(end-1)+1
      wptr(file)=icheck
C
      RETURN
      END
      subroutine swrit(file, buffer, nwords)
      implicit integer (a-z)
      common /pointr/wptr(128), tptr(128)
      dimension buffer(nwords)
      fword = wptr(file)
       icheck=mround(fword-1)+1
       if(icheck.ne.fword) then
         write(6,*)' swrit: fword is not at a sector boundary round up '
         fword=icheck
       endif
c
      call wwritw(file, buffer, nwords, fword, end)
c
c     put word pointer on a sector boundary
c
      icheck=mround(end-1)+1
      wptr(file)=icheck
c
      return
      end
      subroutine wreadw(file,buffer,nwords,fword,nxtwrd)
c     SUBROUTINE WREAD(FILE,BUFFER,NWORDS,FWORD,NXTWRD)
C
      IMPLICIT INTEGER(A-Z)
C
      COMMON /POINTR/ WPTR(128),TPTR(128)
C
      DIMENSION BUFFER(NWORDS)
      DIMENSION RECBUF(1024)
C
CDBG
C     WRITE(6,1000) FILE,NWORDS,FWORD
C1000 FORMAT(' WREAD:    FILE',I3,I12,' WORDS     FIRST WORD',I8)
CDBG
C
C     ---- DETERMINE WHICH RECORDS NEED TO BE READ ----
C
      RECLEN=1024
      FREC=(FWORD-1)/RECLEN + 1
      LWORD=FWORD+NWORDS-1
      LREC=(LWORD-1)/RECLEN + 1
C
C      ---- LOOP OVER THE RECORDS ----
C
      DO 10 IREC=FREC,LREC
C
         READ(FILE,REC=IREC,IOSTAT=IERR) RECBUF
C
C
      IF(IERR.NE.0) THEN
      WRITE(6,3000) IERR,FILE,IREC
 3000 FORMAT(/'  ERROR CODE IN FORTRAN READ (WREAD): ',I4,' FILE ',I4,
     1 ' IREC ',I4)
      CALL EXIT
      ENDIF
C           -- TRANSFER THE APPROPRIATE PORTION --
C           --      OF RECBUF TO BUFFER         --
C
            IF(IREC.EQ.FREC) THEN
C
               OFFSET=FWORD-(IREC-1)*RECLEN - 1
               COUNT=RECLEN-OFFSET
               IF(COUNT.GT.NWORDS) COUNT=NWORDS
               DO 20 I=1,COUNT
                  BUFFER(I)=RECBUF(I+OFFSET)
   20          CONTINUE
C
            ELSE IF(IREC.EQ.LREC) THEN
C
               COUNT=LWORD-(IREC-1)*RECLEN
               OFFSET=NWORDS-COUNT
               DO 30 I=1,COUNT
                  BUFFER(I+OFFSET)=RECBUF(I)
   30          CONTINUE
C
            ELSE
C
              OFFSET=RECLEN*(IREC-1)+1-FWORD
              DO 40 I=1,RECLEN
              buffer(i+offset) = recbuf(i)
   40         CONTINUE
C
            END IF
C
C
   10 CONTINUE
C
C      ---- UPDATE POINTER ----
C
C
C     WPTR IS THE NEXT WORD AFTER THE LAST WORD TRANSFERRED
      WPTR(FILE)=LWORD+1
      NXTWRD=WPTR(FILE)
C
      RETURN
      END
      subroutine wwritw(file,buffer,nwords,fword,nxtwrd)
c     SUBROUTINE WWRIT(FILE,BUFFER,NWORDS,FWORD,NXTWRD)
C
      IMPLICIT INTEGER(A-Z)
C
      COMMON /POINTR/ WPTR(128),TPTR(128)
C
      DIMENSION BUFFER(NWORDS)
      DIMENSION RECBUF(1024)
C
CDBG
C     WRITE(6,1000) FILE,NWORDS,FWORD
C1000 FORMAT(' WWRIT:    FILE',I3,I12,' WORDS   FIRST WORD:',I8)
CDBG
C     ---- DETERMINE WHICH RECORDS NEED TO BE WRITTEN ----
C
      RECLEN=1024
      FREC=(FWORD-1)/RECLEN + 1
      LWORD=FWORD+NWORDS-1
      LREC=(LWORD-1)/RECLEN + 1
C
C      ---- LOOP OVER THE RECORDS ----
C
      DO 10 IREC=FREC,LREC
C
C           DOES IREC HAVE TO BE READ?
C           IF(IREC.LE.TPTR(FILE)) THEN
c              READ(FILE,REC=IREC,ERR=300) RECBUF
c*needed only on Cray, but the end equals is needed to make it work
c* on alliant the way is does on the Cray
             READ(FILE,REC=IREC,ERR=300,end=300) RECBUF
C           ENDIF
            GO TO 310
  300       CONTINUE
C           WRITE(6,1300) IREC
C1300       FORMAT(/'  IN WWRIT, FIRST ACCESS OF RECORD',I3)
  310       CONTINUE
C
C           -- UPDATE RECBUF --
C
            IF(IREC.EQ.FREC) THEN
C
               OFFSET=FWORD-(IREC-1)*RECLEN - 1
               COUNT=RECLEN-OFFSET
               IF(COUNT.GT.NWORDS) COUNT=NWORDS
               DO 20 I=1,COUNT
                  RECBUF(I+OFFSET)=BUFFER(I)
   20          CONTINUE
C
            ELSE IF(IREC.EQ.LREC) THEN
C
               COUNT=LWORD-(IREC-1)*RECLEN
               OFFSET=NWORDS-COUNT
               DO 30 I=1,COUNT
                  RECBUF(I)=BUFFER(I+OFFSET)
   30          CONTINUE
C
            ELSE
C
              OFFSET=RECLEN*(IREC-1)+1-FWORD
              DO 40 I=1,RECLEN
                 RECBUF(I)=BUFFER(I+OFFSET)
   40         CONTINUE
C
            END IF
C
C      ---- WRITE RECBUF ONTO RECORD IREC ----
C
            WRITE(FILE,REC=IREC,IOSTAT=IERR) RECBUF
      IF(IERR.NE.0) THEN
      WRITE(6,14) IERR,FILE
   14 FORMAT(//' ERROR CODE FOR WRITE OF RECBUF ',I4,' FILE=',I3)
      CALL PERROR('FROM WWRITW' )
      CALL EXIT
      ENDIF
C
   10 CONTINUE
C
C      ---- UPDATE POINTERS ----
C
C     TPTR IS THE HIGHEST RECORD WRITTEN
      IF(LREC.GT.TPTR(FILE)) TPTR(FILE)=LREC
C
C     WPTR IS THE NEXT WORD AFTER THE LAST WORD TRANSFERRED
      WPTR(FILE)=LWORD+1
      NXTWRD=WPTR(FILE)
C
C
      RETURN
      END


C********************************************************************
C*  THESE ARE  IO, MATH, AND ROOTS ROUTINES USED BY Y.Y.            *
C*  LAST UPDATED ON SEPTEMBER 28, 1984 BY YUKIO YAMAGUCHI           *
C********************************************************************
C*  LAST UPDATED ON MARCH 1, 1985  BY RICHARD REMINGTON             *
C*  MODIFIED TSTART TO SEND TIMING INFO TO JOBCNTRL FILES           *
C********************************************************************
C*  LAST UPDATED ON MARCH 7, 1985  BY RICHARD REMINGTON             *
C*SUBROUTINE ZFILE MODIFIED,WFILE NOT REMOVED (USED IN SCF & GUGACI)*
C********************************************************************
C*  LAST UPDATED ON MARCH 7, 1985     BY RICHARD REMINGTON AND Y.Y. *
C*  SUBROUTINE  "MWRIT" CHANGED AS PER YUKIO.                       *
C********************************************************************
C
       SUBROUTINE TSTART(ITAPE)
C        THIS SUBROUTINE USES THE IBM UTILITY TXTLIB TO PRINTOUT
C        DATE AND TIMING INFORMATION TO FILE (ITAPE)
C
c      INTEGER IERROR
c      REAL*4 ECPU,ETIME,ETCPU,EECPU,EETIME,EETCPU,DECPU,DETCPU,DETIME
c      CHARACTER*32  DATTIM
C
c   1  FORMAT(4X,' TIMING INFORMATION :  TSTART CALLED         ',A23,/)
c   2  FORMAT(9X,' VIRTUAL CPU TIME (LOGON)    =',F13.4,' SEC',5X,
c    + ' OR ',I9,' HOUR(S) ',I3,' MINUTE(S) ',F10.4,' SECONDS ')
c   3  FORMAT(9X,' TOTAL CPU TIME (LOGON)      =',F13.4,' SEC',5X,
c    + ' OR ',I9,' HOUR(S) ',I3,' MINUTE(S) ',F10.4,' SECONDS ')
c   4  FORMAT(9X,' TOTAL CONNECT TIME (LOGON)  =',F13.4,' SEC',5X,
c    + ' OR ',I9,' HOUR(S) ',I3,' MINUTE(S) ',F10.4,' SECONDS ')
c   5  FORMAT(4X,' TIMING INFORMATION :  TSTOP CALLED          ',A23,//
c    + 4X,' ELAPSED TIME SINCE LAST CALL TO TSTART ')
c   6  FORMAT(9X,' PROGRAM VIRTUAL CPU TIME    =',F13.4,' SEC',5X,
c    + ' OR ',I9,' HOUR(S) ',I3,' MINUTE(S) ',F10.4,' SECONDS ')
c   7  FORMAT(9X,' PROGRAM TOTAL CPU TIME      =',F13.4,' SEC',5X,
c    + ' OR ',I9,' HOUR(S) ',I3,' MINUTE(S) ',F10.4,' SECONDS ')
c   8  FORMAT(9X,' PROGRAM TOTAL CONNECT TIME  =',F13.4,' SEC',5X,
c    + ' OR ',I9,' HOUR(S) ',I3,' MINUTE(S) ',F10.4,' SECONDS ')
c  10  FORMAT(/,130('*'))
c  20  FORMAT(33X,F13.5)
c  26  FORMAT(9X,' TOTAL JOB VIRTUAL CPU TIME  =',F13.4,' SEC',5X,
c    + ' OR ',I9,' HOUR(S) ',I3,' MINUTE(S) ',F10.4,' SECONDS ')
c  27  FORMAT(9X,' TOTAL JOB TOTAL CPU TIME    =',F13.4,' SEC',5X,
c    + ' OR ',I9,' HOUR(S) ',I3,' MINUTE(S) ',F10.4,' SECONDS ')
c  28  FORMAT(9X,' TOTAL JOB CONNECT TIME      =',F13.4,' SEC',5X,
c    + ' OR ',I9,' HOUR(S) ',I3,' MINUTE(S) ',F10.4,' SECONDS ')
C
c      CALL DATETM(DATTIM,ECPU,ETIME,ETCPU)
C
c      IVHRS   = ECPU/3600
c      IVMINS  = (ECPU - IVHRS*3600)/60
c      VSECS   = AMOD(ECPU,60.0)
C
c      ITHRS   = ETCPU/3600
c      ITMINS  = (ETCPU - ITHRS*3600)/60
c      TSECS   = AMOD(ETCPU,60.0)
C
c      ICHRS   = ETIME/3600
c      ICMINS  = (ETIME - ICHRS*3600)/60
c      CSECS   = AMOD(ETIME,60.0)
C
c      IF(ITAPE .NE. 29) GOTO 291
C
c      WRITE(ITAPE,2) ECPU,IVHRS,IVMINS,VSECS
c      WRITE(ITAPE,3) ETCPU,ITHRS,ITMINS,TSECS
c      WRITE(ITAPE,4) ETIME,ICHRS,ICMINS,CSECS
c     RETURN
C
c 291  CONTINUE
C
c      WRITE(ITAPE,10)
c      WRITE(ITAPE,1) DATTIM
C
      RETURN
C
       ENTRY TSTOP(ITAPE)
C
C      THIS COMPUTES THE DELTA TIMES BETWEEN TSTOP AND TSTART
C
c      CALL DATETM(DATTIM,EECPU,EETIME,EETCPU)
C
c      IVHRS   = EECPU/3600
c      IVMINS  = (EECPU - IVHRS*3600)/60
c      VSECS   = AMOD(EECPU,60.0)
C
c      ITHRS   = EETCPU/3600
c      ITMINS  = (EETCPU - ITHRS*3600)/60
c      TSECS   = AMOD(EETCPU,60.0)
C
c      ICHRS   = EETIME/3600
c      ICMINS  = (EETIME - ICHRS*3600)/60
c      CSECS   = AMOD(EETIME,60.0)
C
C      DELTA TIMES    BETWEEN THOSE OF TSTOP AND TSTART
C
c      IF(ITAPE .NE. 29) GOTO 295
C
c      WRITE(ITAPE,2) EECPU,IVHRS,IVMINS,VSECS
c      WRITE(ITAPE,3) EETCPU,ITHRS,ITMINS,TSECS
c      WRITE(ITAPE,4) EETIME,ICHRS,ICMINS,CSECS
      RETURN
C
c 295 CONTINUE
C
c      DECPU   = EECPU - ECPU
c      DETCPU  = EETCPU - ETCPU
c      DETIME  = EETIME - ETIME
C
c      IDVHRS   = DECPU/3600
c      IDVMIN   = (DECPU - IDVHRS*3600)/60
c      DVSECS   = AMOD(DECPU,60.0)
C
c      IDTHRS   = DETCPU/3600
c      IDTMIN   = (DETCPU - IDTHRS*3600)/60
c      DTSECS   = AMOD(DETCPU,60.0)
C
c      IDCHRS   = DETIME/3600
c      IDCMIN   = (DETIME - IDCHRS*3600)/60
c      DCSECS   = AMOD(DETIME,60.0)
C
c      WRITE(ITAPE,10)
c      WRITE(ITAPE,5) DATTIM
C
c      WRITE(ITAPE,6) DECPU,IDVHRS,IDVMIN,DVSECS
c      WRITE(ITAPE,7) DETCPU,IDTHRS,IDTMIN,DTSECS
c      WRITE(ITAPE,8) DETIME,IDCHRS,IDCMIN,DCSECS
C
      RETURN
      END
C
C
       subroutine iwfile(itape,isect)
c      SUBROUTINE WFILE(ITAPE,ISECT)
C
C THIS ROUTINE ZEROES OUT A FILE ON DISK
C
      IMPLICIT INTEGER (A-Z)
      INTEGER BUF(224)
C
      CALL IZERO(BUF,224)
      NSEC = ISECT / 2
      IF(2*NSEC.NE.ISECT) NSEC = NSEC + 1
      DO 1  I=1,NSEC
   1  CALL RWRIT(ITAPE,BUF,224)
C
      RETURN
      END
       subroutine izfile(itape,isect)
c      SUBROUTINE ZFILE(ITAPE,ISECT)
C
C THIS ROUTINE ZEROES OUT A FILE ON DISK
C
      IMPLICIT INTEGER (A-Z)
      COMMON /IOBUFF/ BUFF(1024)
CTJL  INTEGER BUF(1024)
C         IN BLOCK DATA IOBUF IS 1024 INTEGER WORD ZEROES
CTJL  CALL IZERO(BUF,1024)
C     NSEC = ISECT * 2
C     IF(2*NSEC.NE.ISECT) NSEC = NSEC + 1
      DO 1  I=1,ISECT
   1  CALL RWRIT(ITAPE,BUFF,1024)
C
      RETURN
      END
      subroutine itc1(file,n,string)
         integer file
      character*1 c1,string
      c1=char(file)
      string=c1
      return
      end
      subroutine itc2(file,n,string)
      integer file
       character*2 c2,string
       c2=char(file)
       string=c2
       return
       end
       subroutine irfile(itape)
c      SUBROUTINE RFILE(ITAPE)
C
C USE THIS FOR ALL UNFORMATTED FILES
C
      IMPLICIT INTEGER (A-Z)
C
      COMMON /POINTR/ POS(99)
      COMMON /SECT/ SECTOR
C
       INTEGER*4  IERR
C
      CHARACTER*20 FILEID
      CHARACTER*10 FFTYPE,STRING,BLANKS
      CHARACTER*8  FNAME,FTYPE
      CHARACTER*2  FMODE,CHAR2,RECF
      CHARACTER*1  BLANK,CHAR1
      ITAP3 = 3
      ITAP6 = 6
C
C  READ IN THE "FILENAME" FROM LABEL LINE OF INPUT OF INPUT DECK
C  AND OBTAIN POINTER TO INPUT FILE FOR RESTORING POINTER AT RETURN
C
      INQUIRE (UNIT=5,NEXTREC=NXR5)
CRBR  WRITE(ITAP3,*) ' INPUT FILE POINTER = ',NXR5
      IERR = 0
      CALL LOCATOR(5,'# FILES ##',IERR)
      IF(IERR.NE.0) GOTO 100
C
      READ(5,1) FNAME
   1  FORMAT(A8)
      GO TO 101
 100  WRITE(ITAP3,*) ' FILENAME DEFAULTING TO "FI", IERR=',IERR
      FNAME = 'FI'
 101  CONTINUE
C
      IF (ITAPE .GT. 9) GOTO 110
      CALL ITC1(ITAPE,1,CHAR1)
      STRING = CHAR1
      GOTO 112
 110  CONTINUE
      CALL ITC2(ITAPE,2,CHAR2)
      STRING = CHAR2
 112  CONTINUE
C
      FTYPE = 'FILE'//STRING
      FFTYPE = 'FILE'//STRING
CDBG  WRITE(ITAP3,*) ' ITAPE=',ITAPE,' STRING=',STRING,' FTYPE=',FTYPE
C
C  READ IN "FILEMODE" FROM THE "FILEMODE DATA" FILE
C
c      IERR = 0
c      CALL LOCATE(28,FFTYPE,IERR)
c      IF(IERR .NE. 0) GO TO 120
C
c      BACKSPACE 28
c      READ(28,2) FMODE
c   2  FORMAT(10X,A2)
c      GO TO 121
c 120  FMODE = 'A4'
c      WRITE(ITAP3,*) ' FILEMODE DEFAULTING TO A4; FMODE=',FMODE
c 121  CONTINUE
C
C A SECTOR LENGTH OF 4096 BYTES IS SET BY THE ROUTINE DABUFF
C
      SECTOR = 1024
      TYPE = 5
      ICODE = 3
      LREC = 4096
      BLK = 4096
c      RECF = 'FB'
CDBG  WRITE(ITAP3,*) ' FNAME=',FNAME,' FTYPE=',FTYPE,' FMODE =',FMODE
c      FILEID = FNAME//' '//FTYPE//' '//FMODE
          fileid = fname//'.'//string
      WRITE(ITAP3,*) ' FILEID =',FILEID
      IERR = 0
c      CALL UOPEN(ITAPE,FILEID,TYPE,ICODE,LREC,BLK,RECF,0,IERR)
c   insert sfile call to open at this point
Cibm
C          open(itape,file=fileid,
C    #             access='direct',
C    #             initialsize = 100,
C    #             recl= lrec,
C    #             recordtype='fixed',
C    #             iostat=ierr )
           open(itape,file=fileid,
     #             access='direct',
     #             recl= lrec,
     #             iostat=ierr )
Cibm
c 
C
      IF (IERR .EQ. 0) GO TO 130
      WRITE(*,*) ' ERROR ENCOUNTERED OPENING FILE ',ITAPE,' IERR =',IERR
      WRITE(*,*) ' FILEID =',FILEID
      WRITE(*,*) ' FNAME =',FNAME,' FTYPE =',FTYPE,' FMODE =',FMODE
      WRITE(*,*) ' TYPE  =',TYPE,' ICODE =',ICODE,' LREC  =',LREC
      WRITE(*,*) ' BLK   =',BLK,' RECF  =',RECF,' IERR  =',IERR
C     STOP
 130  CONTINUE
C
C PLACE THE POINTER OF ITAPE AT THE BEGINNING
C
      CALL SREW(ITAPE)
C
C RETURN
C
      RETURN
      END
c       subroutine irgetsa(itape,iadr)
       SUBROUTINE RGETSA(ITAPE,IADR)
C
C THIS ROUTINE RETURNS THE CURRENT SECTOR ADDRESS FOR FILES USING
C THE SREAD/SWRIT AND/OR RREAD/RWRIT IO ROUTINES.
C
      IMPLICIT INTEGER (A-Z)
C
       common /pointr/wptr(128),tptr(128)
c      COMMON /SECT/ SECTOR
C
c      CALL RDPOS(ITAPE,IPOS,IERR)
c      IREC = (IPOS-1)/SECTOR + 1
c      TEST = SECTOR*(IREC-1) + 1
c      IF(IPOS.NE.TEST.OR.IERR.NE.0) THEN
c      WRITE(*,*) ' ERROR ENCOUNTERED IN RGETSA FOR FILE ',ITAPE
c      WRITE(*,*) ' IERR,IPOS,TEST = ',IERR,IPOS,TEST
        reclen=1024
       fword = wptr(itape)
       irec = (fword-1)/reclen + 1
c      STOP
c      END IF
C
      IADR=IREC
C
      RETURN
      END
c       subroutine irsetsa(itape,iadr)
       SUBROUTINE RSETSA(ITAPE,IADR)
C
C THIS ROUTINE SETS THE POINTER OF A UTILITIES FILE TO THE DESIRED
C SECTOR ADDRESS.
C
      IMPLICIT INTEGER (A-Z)
C
       common/pointr/wptr(128),tptr(128)
c      COMMON /SECT/ SECTOR
C
        IPOS = SEC2I(IADR-1) + 1
c      CALL STPOS(ITAPE,IPOS,IERR)
c      IF(IERR.NE.0) THEN
c      WRITE(*,*) ' ERROR ENCOUNTERED IN RSETSA FOR FILE ',ITAPE
c      WRITE(*,*) ' IERR,IPOS,IADR = ',IERR,IPOS,IADR
c      STOP
c      END IF
       wptr(itape) = ipos
c
C
      RETURN
      END
       subroutine ilocate(input,token,ierror)
c      SUBROUTINE LOCATE(INPUT,TOKEN,IERROR)
C
C     ----- SEARCH THROUGH INPUT FILE FOR TOKEN BEGINNING
C           WITH # TO LOCATE INPUT FOR PROGRAM.  IERROR IS
C           SET TO 0 IF NO ERRORS, 1 IF ANY ERROR OCCURS.
C
C
      CHARACTER*10 TOKEN,LINE
C
      REWIND (UNIT=INPUT,ERR=99)
C
    1 CONTINUE
      READ (UNIT=INPUT,FMT='(A10)',END=99,ERR=99) LINE
      IF (LINE.NE.TOKEN) GO TO 1
C
      IERROR=0
      RETURN
C
C
   99 CONTINUE
      IERROR=1
      RETURN
      END
      INTEGER FUNCTION I2SEC(N)
CC    I2SEC=(N+191)/192
      I2SEC=(N+1023)/1024
      RETURN
      END
      INTEGER FUNCTION SEC2I(N)
CC    SEC2I=192*N
      SEC2I=1024*N
      RETURN
      END
CDRY   ON OS3.1 FOLLOWING FUNCTION CAUSED TROUBLE
CDRY   ON OS3.2 IT SEEMS TO BE NEED
CDRY   10/29/93  iround changed to mround
      INTEGER FUNCTION MROUND(N)
CC    MROUND=(N+191)/192*192
      MROUND=(N+1023)/1024*1024
      RETURN
      END
      INTEGER FUNCTION MTRUNC(N)
CC    ITRUNC=N/192*192
      MTRUNC=N/1024*1024
      RETURN
      END
      INTEGER FUNCTION IADTWP(N)
      IADTWP=(N+3)/2
      RETURN
      END
      INTEGER FUNCTION WPADTI(N)
      WPADTI=N*2-1
      RETURN
      END
       subroutine gwrite(file,buffer,nwords,fword,nxtwrd)
c      subroutine wwritw(file,buffer,nwords,fword,nxtwrd)
c
        implicit integer(a-z)
c
       common/pointr/ wptr(128),tptr(128)
c
       dimension buffer(nwords)
       dimension recbuf(1024)
c
c
       reclen= 1024
       frec = (fword-1)/(reclen) + 1
       lword = fword + nwords - 1
       lrec = (lword - 1)/reclen + 1
c
c        loop over records
c
          do 10  irec = frec, lrec
c
c
c      does irec have to be read
c       if (irec.le. tptr(file) then
c          read ( file,rec=irec,err=300) recbuf
c       needed only on the Cray, but the end equals is needed to make
c       the alliant behave the same way the Cray does.
C         read(file,rec=irec,err = 300, end = 300) recbuf
Cibm
          read(file,rec=irec,err = 300) recbuf
Cibm
c       endif
        go to 310
  300   continue
C       write (6,1300) irec
C 1300  format(/'  in wwrit, first access of record ', i3)
  310   continue
c
c      --- update recbuf ----
c
         if (irec.eq.frec) then
c
              offset = fword -(irec-1)*reclen -1
              count = reclen - offset
              if (count.gt.nwords) count = nwords
              do 20 i=1,ncount
                 recbuf(i+offset) = buffer(i)
   20         continue
c
        else if (irec.eq.lrec) then
c
              count = lword - (irec-1)*reclen
              offset = nwords - count
              do 30 i = 1,ncount
                  recbuf(i) = buffer (i + offset)
   30         continue
c
        else
c
              offset = reclen * (irec-1) + 1 - fword
              do 40 i = 1,reclen
                   recbuf(i) = buffer (i+offset)
   40          continue
c
        end if
c
c      ---  write recbuf onto record irec  -----
c
               write( file,rec=irec,iostat=ierr) recbuf
        if( ierr.ne.0) then
        write (6,14) ierr
   14   format(//'  error code for write of recbuf  ',i4)
        call exit
        endif
c
   10   continue
c
c
c     ----  update pointers ---
c
c      tptr is the highest record written
       if (lrec.gt.tptr(file)) tptr(file) = lrec
c
c      wptr is the next word after the last word transferred
       wptr (file) = lword + 1
       nxtwrd = wptr (file)
c
c
       return
       end
      subroutine iwritw(itape,array,nlen,fword,lword)
c      original IBM version of WWRITW
c      SUBROUTINE WWRITW(ITAPE,ARRAY,NLEN,FWORD,LWORD)
C
      IMPLICIT INTEGER (A-Z)
C
      COMMON /SECT/ SECTOR
C     COMMON / IOBUFF / BUFF(256)
C
      DIMENSION ARRAY(NLEN)
C
CTJL  WRITE(*,*) ' IN WWRIT, FWORD,NLEN,ITAPE = ',FWORD,NLEN,ITAPE
      IOP = 2
      CALL RWFIL(ITAPE,IOP,FWORD,ARRAY,NLEN,IERR)
      IF(IERR.NE.0) THEN
      WRITE(*,*) ' ERROR IN WWRITW FOR FILE ',ITAPE
      WRITE(*,*) ' IERR = ',IERR
      STOP
      END IF
C
C UPDATE LWORD AND THEN RETURN
C
      LWORD = FWORD + NLEN
C
      RETURN
      END
       subroutine iwreadw(itape,array,nlen,fword,lword)
c      SUBROUTINE WREADW(ITAPE,ARRAY,NLEN,FWORD,LWORD)
C
      IMPLICIT INTEGER (A-Z)
C
      COMMON /SECT/ SECTOR
C
      DIMENSION ARRAY(NLEN)
C
CTJL  WRITE(*,*) ' IN WREAD, FWORD,NLEN,ITAPE = ',FWORD,NLEN,ITAPE
      IOP = 1
      CALL RWFIL(ITAPE,IOP,FWORD,ARRAY,NLEN,IERR)
      IF(IERR.NE.0) THEN
      WRITE(*,*) ' ERROR IN WREADW FOR FILE ',ITAPE
      WRITE(*,*) ' IERR = ',IERR
      STOP
      END IF
C
C UPDATE LWORD AND THEN RETURN
C
      LWORD = FWORD + NLEN
C
      RETURN
      END
       subroutine irread(itape,array,nlen,irec)
c      SUBROUTINE RREAD(ITAPE,ARRAY,NLEN,IREC)
C
C THIS ROUTINE READS NLEN INTEGER WORDS INTO ARRAY STARTING AT
C SECTOR IREC.
C
      IMPLICIT INTEGER (A-Z)
C
      COMMON /SECT/ SECTOR
      COMMON /IOBUFF/ BUFF(1024)
C
      DIMENSION ARRAY(NLEN)
C
      IOP = 1
      IPOS = (IREC-1)*SECTOR + 1
CTJL  WRITE(*,*) ' IN RREAD FOR FILE',ITAPE
CTJL  WRITE(*,*) ' IREC,IPOS,NLEN = ',IREC,IPOS,NLEN
      CALL RWFIL(ITAPE,IOP,IPOS,ARRAY,NLEN,IERR)
      IF(IERR.NE.0) THEN
      WRITE(*,*) ' ERROR ENCOUNTERED IN RREAD FOR FILE ',ITAPE
      WRITE(*,*) ' IERR = ',IERR
      STOP
      END IF
C
C WE MUST POSITION THE FILE POINTER AT THE BEGINNING OF THE NEXT SECTOR
C
      EXTRA = MOD(NLEN,SECTOR)
      IF(EXTRA.NE.0) THEN
      IOP = 2
      IPOS = 0
      NW = SECTOR - EXTRA
      CALL RWFIL(ITAPE,IOP,IPOS,BUFF,NW,IERR)
      IF(IERR.NE.0) THEN
      WRITE(*,*) ' ERROR ENCOUNTERED WHILE POSITIONING FILE POINTER '
      WRITE(*,*) ' FOR FILE ',ITAPE,'   IN RREAD,  IERR = ',IERR
      STOP
      END IF
      END IF
C
      RETURN
      END
       subroutine irwrit(itape,array,nlen,irec)
c      SUBROUTINE RWRIT(ITAPE,ARRAY,NLEN,IREC)
C
C THIS ROUTINE WRITES NLEN INTEGER WORDS FROM ARRAY TO FILE ITAPE
C STARTING AT SECTOR IREC.
C
      IMPLICIT INTEGER (A-Z)
C
      COMMON /SECT/ SECTOR
      COMMON /IOBUFF/ BUFF(1024)
C
      DIMENSION ARRAY(NLEN)
C
      IOP = 2
      IPOS = (IREC-1)*SECTOR + 1
CTJL  WRITE(*,*) ' IN RWRIT FOR FILE',ITAPE
CTJL  WRITE(*,*) ' IREC,IPOS,NLEN = ',IREC,IPOS,NLEN
      CALL RWFIL(ITAPE,IOP,IPOS,ARRAY,NLEN,IERR)
      IF(IERR.NE.0) THEN
      WRITE(*,*) ' ERROR ENCOUNTERED IN RWRIT FOR FILE ',ITAPE
      WRITE(*,*) ' IERR = ',IERR
      STOP
      END IF
C
C WE MUST POSITION THE FILE POINTER AT THE BEGINNING OF THE NEXT SECTOR
C
      EXTRA = MOD(NLEN,SECTOR)
      IF(EXTRA.NE.0) THEN
      IOP = 2
      IPOS = 0
      NW = SECTOR - EXTRA
      CALL RWFIL(ITAPE,IOP,IPOS,BUFF,NW,IERR)
      IF(IERR.NE.0) THEN
      WRITE(*,*) ' ERROR ENCOUNTERED WHILE POSITIONING FILE POINTER '
      WRITE(*,*) ' FOR FILE ',ITAPE,'   IN RWRIT,  IERR = ',IERR
      STOP
      END IF
      END IF
C
      RETURN
      END
       subroutine isread(itape,array,nlen)
c      SUBROUTINE SREAD(ITAPE,ARRAY,NLEN)
C
C THIS ROUTINE READS NLEN INTEGER WORDS FROM ITAPE INTO ARRAY
C STARTING AT THE CURRENT POINTER LOCATION.
C
      IMPLICIT INTEGER (A-Z)
C
      COMMON /IOBUFF/ BUFF(1024)
      COMMON /SECT/ SECTOR
C
      DIMENSION ARRAY(NLEN)
C
CTJL  WRITE(*,*) ' IN SREAD,  ITAPE,NLEN = ',ITAPE,NLEN
      IOP = 1
      IPOS = 0
      CALL RWFIL(ITAPE,IOP,IPOS,ARRAY,NLEN,IERR)
      IF(IERR.NE.0) THEN
      WRITE(*,*) ' ERROR ENCOUNTERED IN SREAD FOR FILE ',ITAPE
      WRITE(*,*) ' IERR = ',IERR
      STOP
      END IF
C
C WE MUST POSITION THE FILE POINTER AT THE BEGINNING OF THE NEXT SECTOR
C
      EXTRA = MOD(NLEN,SECTOR)
      IF(EXTRA.NE.0) THEN
      IOP = 2
      IPOS = 0
      NW = SECTOR - EXTRA
      CALL RWFIL(ITAPE,IOP,IPOS,BUFF,NW,IERR)
      IF(IERR.NE.0) THEN
      WRITE(*,*) ' ERROR ENCOUNTERED WHILE POSITIONING FILE POINTER '
      WRITE(*,*) ' FOR FILE ',ITAPE,'   IN SREAD,  IERR = ',IERR
      STOP
      END IF
      END IF
C
      RETURN
      END
       subroutine iswrit(itape,array,nlen)
c      SUBROUTINE SWRIT(ITAPE,ARRAY,NLEN)
C
C THIS ROUTINE WTITES NLEN INTEGER WORDS FROM ARRAY TO FILE ITAPE
C STARTING AT THE CURRENT POINTER LOCATION.
C
      IMPLICIT INTEGER (A-Z)
C
      COMMON /IOBUFF/ BUFF(1024)
      COMMON /SECT/ SECTOR
C
      DIMENSION ARRAY(NLEN)
C
CTJL  WRITE(*,*) ' IN SWRIT,  ITAPE,NLEN = ',ITAPE,NLEN
CTJL  WRITE(*,*) ' N ',ARRAY
      IOP = 2
      IPOS = 0
      CALL RWFIL(ITAPE,IOP,IPOS,ARRAY,NLEN,IERR)
      IF(IERR.NE.0) THEN
      WRITE(*,*) ' ERROR ENCOUNTERED IN SWRIT FOR FILE ',ITAPE
      WRITE(*,*) ' IERR = ',IERR
      STOP
      END IF
C
C WE MUST POSITION THE FILE POINTER AT THE BEGINNING OF THE NEXT SECTOR
C
      EXTRA = MOD(NLEN,SECTOR)
      IF(EXTRA.NE.0) THEN
      IOP = 2
CTJL  CALL RDPOS(ITAPE,IPOS,IERR)
CTJL  WRITE(*,*) ' AFTER SWRIT IPOS,IERR = ',IPOS,IERR
      IPOS = 0
      NW = SECTOR - EXTRA
CTJL  WRITE(*,*) ' NLEN,SECTOR,EXTRA,NW = ',NLEN,SECTOR,EXTRA,NW
      CALL RWFIL(ITAPE,IOP,IPOS,BUFF,NW,IERR)
      IF(IERR.NE.0) THEN
      WRITE(*,*) ' ERROR ENCOUNTERED WHILE POSITIONING FILE POINTER '
      WRITE(*,*) ' FOR FILE ',ITAPE,'   IN SWRIT,  IERR = ',IERR
      STOP
      END IF
      END IF
C
      RETURN
      END
      SUBROUTINE EIGOUT(A,B,C,NAD,NBD,M,N,IOUT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NAD,NBD),B(NAD),C(NAD)
      CHARACTER*6 LINE
      LINE='------'
    1 FORMAT(2X,10(7X,I5))
    2 FORMAT(2X,21A6)
    3 FORMAT(2X,I2,2X,10F12.7)
    4 FORMAT(/,6X,10F12.7)
    5 FORMAT(/)
C
      II=0
      JJ=0
  200 II=II+1
      JJ=JJ+1
      KK=10*JJ
      NN=N
      IF(N.GT.KK) NN=KK
      LL=2*(NN-II+1)+1
      WRITE(IOUT,1) (I,I=II,NN)
      WRITE(IOUT,2) (LINE,I=1,LL)
      DO 101 I=1,M
  101 WRITE(IOUT,3) I,(A(I,J),J=II,NN)
      WRITE(IOUT,4) (B(J),J=II,NN)
      WRITE(IOUT,4) (C(J),J=II,NN)
      IF(N.LE.KK) GO TO 201
      WRITE(IOUT,5)
      II=KK
      GO TO 200
  201 RETURN
      END
      INTEGER FUNCTION NREC(NLEN)
      IMPLICIT INTEGER (A-Z)
CC    IADD=(NLEN+8188)/8189
      IADD=(NLEN+255)/256
      NREC=IADD
      RETURN
      END
      SUBROUTINE INITMF(ISTART)
C   THIS PROGRAM INTIALIZES USEAGE OF THE MASTER FILE
C**********************************************************
C***LAST UPDATED ON FEBRUARY 06, 1985 BY YUKIO YAMAGUCHI***
C**********************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/CALIF/LPARA(1024),APARA(1024)
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4
      COMMON/LOCAT/LOCS(1024)
      COMMON/MFSEC/MFILE,NSECT
C
C   DEFINE CONSTANTS
      MFILE=40
      NSECT=1024
      MAXBAS=52
      MAXBUF=4096
      MAXBF2=MAXBUF*2
      MAXBF4=MAXBUF*4
      IF(ISTART.EQ.0) GO TO 201
C
C   FIND LOCATIONS OF PARAMETERS AND MATRICES IN THE MASTER FILE
      CALL RFILE(MFILE)
      CALL RREAD(MFILE,LOCS,NSECT,1)
      CALL MREAD(LPARA,2)
      CALL MREAD(APARA,3)
C
  201 CONTINUE
      RETURN
      END
      SUBROUTINE MWRIT(IA,IADD)
      IMPLICIT INTEGER (A-Z)
      DIMENSION IA(1)
      COMMON/LOCAT/LOCS(1024)
      COMMON/MFSEC/MFILE,NSECT
C
C  NO.  :  CONTENTS
C
C   TRIA. = LOWER TRIANGULAR MATRIX
C   SQUA. = SQUARE MATRIX
C   RECT. = RECTANGULAR MATRIX
C   1   : LOCATION
C   2   : PARAMETERS
C   3   : ENERGIES AND REAL CONSTANTS
C   4   : MO ORDERING IN SCF
C   5   : MO ORDERING IN DRT
C   6   : MO CODE FOR CI
C   7   : MO CODE FOR MCSCF
C   8   : NUCLEAR CHARGE AND GEOMETRY
C   9   : AO-SO TRANSFORMATION MATRIX (RECT.)
C  10   : EIGENVALUES IN PITZER'S SCF (SQUA.)
C  11   : OCCUPATION IN PITZER'S SCF
C  12   : SO-MO EIGENVECTORS IN PITZER'S SCF (SQUA.)
C  13   : AO-MO EIGENVECTORS IN PITZER'S SCF (RECT.)
C  14   : SO OVERLAP INTEGRALS IN PITZER'S SCF (TRIA.)
C  15   : ONE ELECTRON SO INTEGRALS IN PITZER'S SCF (TRIA.)
C  16   : ONE ELECTRON MO INTEGRALS IN PITZER'S SCF (TRIA.)
C  17   : EIGENVALUES IN SORTED SCF
C  18   : OCCUPATION IN SORTED SCF
C  19   : SO-MO EIGENVECTORS IN SORTED SCF (SQUA.)
C  20   : AO-MO EIGENVECTORS IN SORTED SCF (RECT.)
C  21   : SO OVERLAP INTEGRALS IN SORTED SCF (TRIA.)
C  22   : ONE ELECTRON SO INTEGRALS IN SORTED SCF (TRIA.)
C  23   : ONE ELECTRON MO INTEGRLAS IN SORTED SCF (TRIA.)
C  24   : LAGRANGIAN MATRIX FOR SCF IN AO BASIS (TRIA.)
C  25   : LAGRANGIAN MATRIX FOR SCF IN MO BASIS (TRIA.)
C  26   : K MATRIX FOR HIGH SPIN OPEN-SHELL SCF (TRIA.)
C  27   : FIRST ZETA MATRIX  FOR GENERALIZED OPEN-SHELL SCF (TRIA.)
C  28   : SECOND ZETA MATRIX FOR GENERALIZED OPEN-SHELL SCF (TRIA.)
C  29   : THIRD ZETA MATRIX  FOR GENERALIZED OPEN-SHELL SCF (TRIA.)
C  30   : FOURTH ZETA MATRIX FOR GENERALIZED OPEN-SHELL SCF (TRIA.)
C  31   : FIFTH ZETA MATRIX  FOR GENERALIZED OPEN-SHELL SCF (TRIA.)
C  32   : LAGRANGIAN MATRIX FOR CI IN AO BASIS (TRIA.)
C  33   : LAGRANGIAN MATRIX FOR CI IN MO BASIS (SQUA.)
C  34   : ONE PDM IN AO BASIS FOR CI (TRIA.)
C  35   : ONE PDM IN MO BASIS FOR CI (TRIA.)
C  36   : LAGRANGIAN MATRIX FOR MCSCF IN AO BASIS (TRIA.)
C  37   : LAGRANGIAN MATRIX FOR MCSCF IN MO BASIS (SQUA.)
C  38   : ONE PDM IN AO BASIS FOR MCSCF (TRIA.)
C  39   : ONE PDM IN MO BASIS FOR MCSCF (TRIA.)
C
      IADDR=LOCS(IADD)
      NSECTH=NSECT/2
      NLEN=LOCS(IADD+NSECTH)
      CALL RWRIT(MFILE,IA,NLEN,IADDR)
      RETURN
C
      ENTRY MREAD(IA,IADD)
C
      IADDR=LOCS(IADD)
      NSECTH=NSECT/2
      NLEN=LOCS(IADD+NSECTH)
      CALL RREAD(MFILE,IA,NLEN,IADDR)
      RETURN
      END
      SUBROUTINE EIVOUT(A,B,NAD,NBD,M,N,IOUT)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*6 LINE(21)
      DIMENSION A(NAD,NBD),B(NAD)
      DATA LINE / 21*'------' /
    1 FORMAT(2X,10(7X,I5))
    2 FORMAT(2X,21A6)
    3 FORMAT(2X,I2,2X,10F12.7)
    4 FORMAT(/,6X,10F12.7)
    5 FORMAT(/)
C
      II=0
      JJ=0
  200 II=II+1
      JJ=JJ+1
      KK=10*JJ
      NN=N
      IF(N.GT.KK) NN=KK
      LL=2*(NN-II+1)+1
      WRITE(IOUT,1) (I,I=II,NN)
      WRITE(IOUT,2) (LINE(I),I=1,LL)
      DO 101 I=1,M
  101 WRITE(IOUT,3) I,(A(I,J),J=II,NN)
      WRITE(IOUT,4) (B(J),J=II,NN)
      IF(N.LE.KK) GO TO 201
      WRITE(IOUT,5)
      II=KK
      GO TO 200
  201 RETURN
      END
      SUBROUTINE FRQOUT(A,B,NAD,NBD,M,N,IOUT)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*6 LINE(21)
      DIMENSION A(NAD,NBD),B(NBD)
      DATA LINE / 21*'------' /
    1 FORMAT(2X,10(7X,I5))
    2 FORMAT(6X,10(4X,F8.2))
    3 FORMAT(2X,21A6)
    4 FORMAT(2X,I2,2X,10F12.7)
    5 FORMAT(/)
C
      II=0
      JJ=0
  200 II=II+1
      JJ=JJ+1
      KK=10*JJ
      NN=N
      IF(N.GT.KK) NN=KK
      LL=2*(NN-II+1)+1
      WRITE(IOUT,1) (I,I=II,NN)
      WRITE(IOUT,2) (B(I),I=II,NN)
      WRITE(IOUT,3) (LINE(I),I=1,LL)
      DO 101 I=1,M
      WRITE(IOUT,4) I,(A(I,J),J=II,NN)
  101 CONTINUE
      IF(N.LE.KK) GO TO 201
      WRITE(IOUT,5)
      II=KK
      GO TO 200
  201 RETURN
      END
      SUBROUTINE EBC(A,B,C,NI,NK,NJ)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION A(NI,NJ),B(NI,NK),C(NK,NJ)
C
      CALL ZERO(A,NI*NJ)
c      CALL MXMB(B,1,NI,C,1,NK,A,1,NI,NI,NK,NJ)
C
      DO 3 I=1,NI
           DO 2 J=1,NJ
                T=0.0D+00
                DO 1 K=1,NK
                     T=T+B(I,K)*C(K,J)
    1           CONTINUE
                A(I,J)=T
    2      CONTINUE
    3 CONTINUE
      RETURN
C     END
C
C     DIMENSION A(1),B(1),C(1)
C
C     IJ=1
C     JA=1
C     DO 3 J=1,NJ
C     DO 2 I=1,NI
C     IK=I
C     KJ=JA
C     T=0.0D+00
C     DO 1 K=1,NK
C     T=T+B(IK)*C(KJ)
C     IK=IK+NI
C     KJ=KJ+1
C   1 CONTINUE
C     A(IJ)=T
C     IJ=IJ+1
C   2 CONTINUE
C     JA=JA+NK
C   3 CONTINUE
C     RETURN
      END
      SUBROUTINE EBTC(A,B,C,NI,NK,NJ)
C
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION A(NI,NJ),B(NK,NI),C(NK,NJ)
C
      CALL ZERO(A,NI*NJ)
c     CALL MXMB(B,NK,1,C,1,NK,A,1,NI,NI,NK,NJ)
C
      DO 3 I=1,NI
           DO 2 J=1,NJ
                T=0.0D+00
                DO 1 K=1,NK
                     T=T+B(K,I)*C(K,J)
    1           CONTINUE
                A(I,J)=T
    2      CONTINUE
    3 CONTINUE
      RETURN
C     END
C
C     DIMENSION A(1),B(1),C(1)
C
C     IJ=1
C     JA=1
C     DO 3 J=1,NJ
C     KI=1
C     DO 2 I=1,NI
C     KJ=JA
C     T=0.0D+00
C     DO 1 K=1,NK
C     T=T+B(KI)*C(KJ)
C     KI=KI+1
C     KJ=KJ+1
C   1 CONTINUE
C     A(IJ)=T
C     IJ=IJ+1
C   2 CONTINUE
C     JA=JA+NK
C   3 CONTINUE
C     RETURN
      END
      SUBROUTINE EBCT(A,B,C,NI,NK,NJ)
C
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION A(NI,NJ),B(NI,NK),C(NJ,NK)
C
      CALL ZERO(A,NI*NJ)
c     CALL MXMB(B,1,NI,C,NJ,1,A,1,NI,NI,NK,NJ)
C
      DO 3 I=1,NI
           DO 2 J=1,NJ
                T=0.0D+00
                DO 1 K=1,NK
                     T=T+B(I,K)*C(J,K)
    1           CONTINUE
                A(I,J)=T
    2      CONTINUE
    3 CONTINUE
      RETURN
C     END
C
C     DIMENSION A(1),B(1),C(1)
C
C     IJ=1
C     DO 3 J=1,NJ
C     DO 2 I=1,NI
C     IK=I
C     JK=J
C     T=0.0D+00
C     DO 1 K=1,NK
C     T=T+B(IK)*C(JK)
C     IK=IK+NI
C     JK=JK+NJ
C   1 CONTINUE
C     A(IJ)=T
C     IJ=IJ+1
C   2 CONTINUE
C   3 CONTINUE
C     RETURN
      END
      SUBROUTINE EBTCT(A,B,C,NI,NK,NJ)
C
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION A(NI,NJ),B(NK,NI),C(NJ,NK)
C
      CALL ZERO(A,NI*NJ)
c     CALL MXMB(B,NK,1,C,NJ,1,A,1,NI,NI,NK,NJ)
C
      DO 3 I=1,NI
           DO 2 J=1,NJ
                T=0.0D+00
                DO 1 K=1,NK
                     T=T+B(K,I)*C(J,K)
    1           CONTINUE
                A(I,J)=T
    2      CONTINUE
    3 CONTINUE
      RETURN
C     END
C
C     DIMENSION A(1),B(1),C(1)
C
C     IJ=1
C     DO 3 J=1,NJ
C     KI=1
C     DO 2 I=1,NI
C     JK=J
C     T=0.0D+00
C     DO 1 K=1,NK
C     T=T+B(KI)*C(JK)
C     JK=JK+NJ
C     KI=KI+1
C   1 CONTINUE
C     A(IJ)=T
C     IJ=IJ+1
C   2 CONTINUE
C   3 CONTINUE
C     RETURN
      END
      SUBROUTINE EMBC(A,B,C,NI,NK,NJ)
C
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION A(NI,NJ),B(NI,NK),C(NK,NJ)
C
      DO 3 I=1,NI
           DO 2 J=1,NJ
                T=0.0D+00
                DO 1 K=1,NK
                     T=T-B(I,K)*C(K,J)
    1           CONTINUE
                A(I,J)=T
    2      CONTINUE
    3 CONTINUE
      RETURN
C     END
C
C     DIMENSION A(1),B(1),C(1)
C
C     IJ=1
C     JA=1
C     DO 3 J=1,NJ
C     DO 2 I=1,NI
C     IK=I
C     KJ=JA
C     T=0.0D+00
C     DO 1 K=1,NK
C     T=T-B(IK)*C(KJ)
C     IK=IK+NI
C     KJ=KJ+1
C   1 CONTINUE
C     A(IJ)=T
C     IJ=IJ+1
C   2 CONTINUE
C     JA=JA+NK
C   3 CONTINUE
C     RETURN
      END
      SUBROUTINE EMBTC(A,B,C,NI,NK,NJ)
C
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION A(NI,NJ),B(NK,NI),C(NK,NJ)
C
      DO 3 I=1,NI
           DO 2 J=1,NJ
                T=0.0D+00
                DO 1 K=1,NK
                     T=T-B(K,I)*C(K,J)
    1           CONTINUE
                A(I,J)=T
    2      CONTINUE
    3 CONTINUE
      RETURN
C     END
C
C     DIMENSION A(1),B(1),C(1)
C
C     IJ=1
C     JA=1
C     DO 3 J=1,NJ
C     KI=1
C     DO 2 I=1,NI
C     KJ=JA
C     T=0.0D+00
C     DO 1 K=1,NK
C     T=T-B(KI)*C(KJ)
C     KI=KI+1
C     KJ=KJ+1
C   1 CONTINUE
C     A(IJ)=T
C     IJ=IJ+1
C   2 CONTINUE
C     JA=JA+NK
C   3 CONTINUE
C     RETURN
      END
      SUBROUTINE EMBCT(A,B,C,NI,NK,NJ)
C
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION A(NI,NJ),B(NI,NK),C(NJ,NK)
C
      DO 3 I=1,NI
           DO 2 J=1,NJ
                T=0.0D+00
                DO 1 K=1,NK
                     T=T-B(I,K)*C(J,K)
    1           CONTINUE
                A(I,J)=T
    2      CONTINUE
    3 CONTINUE
      RETURN
C     END
C
C     DIMENSION A(1),B(1),C(1)
C
C     IJ=1
C     DO 3 J=1,NJ
C     DO 2 I=1,NI
C     IK=I
C     JK=J
C     T=0.0D+00
C     DO 1 K=1,NK
C     T=T-B(IK)*C(JK)
C     IK=IK+NI
C     JK=JK+NJ
C   1 CONTINUE
C     A(IJ)=T
C     IJ=IJ+1
C   2 CONTINUE
C   3 CONTINUE
C     RETURN
      END
      SUBROUTINE EMBTCT(A,B,C,NI,NK,NJ)
C
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION A(NI,NJ),B(NK,NI),C(NJ,NK)
C
      DO 3 I=1,NI
           DO 2 J=1,NJ
                T=0.0D+00
                DO 1 K=1,NK
                     T=T-B(K,I)*C(J,K)
    1           CONTINUE
                A(I,J)=T
    2      CONTINUE
    3 CONTINUE
      RETURN
C     END
C
C     DIMENSION A(1),B(1),C(1)
C
C     IJ=1
C     DO 3 J=1,NJ
C     KI=1
C     DO 2 I=1,NI
C     JK=J
C     T=0.0D+00
C     DO 1 K=1,NK
C     T=T-B(KI)*C(JK)
C     JK=JK+NJ
C     KI=KI+1
C   1 CONTINUE
C     A(IJ)=T
C     IJ=IJ+1
C   2 CONTINUE
C   3 CONTINUE
C     RETURN
      END
      SUBROUTINE APBC(A,B,C,NI,NK,NJ)
C
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION A(NI,NJ),B(NI,NK),C(NK,NJ)
C
c     CALL MXMB(B,1,NI,C,1,NK,A,1,NI,NI,NK,NJ)
C
      DO 3 I=1,NI
           DO 2 J=1,NJ
                T=A(I,J)
                DO 1 K=1,NK
                     T=T+B(I,K)*C(K,J)
    1           CONTINUE
                A(I,J)=T
    2      CONTINUE
    3 CONTINUE
      RETURN
C     END
C
C     DIMENSION A(1),B(1),C(1)
C
C     IJ=1
C     JA=1
C     DO 3 J=1,NJ
C     DO 2 I=1,NI
C     IK=I
C     KJ=JA
C     T=A(IJ)
C     DO 1 K=1,NK
C     T=T+B(IK)*C(KJ)
C     IK=IK+NI
C     KJ=KJ+1
C   1 CONTINUE
C     A(IJ)=T
C     IJ=IJ+1
C   2 CONTINUE
C     JA=JA+NK
C   3 CONTINUE
C     RETURN
      END
      SUBROUTINE APBTC(A,B,C,NI,NK,NJ)
C
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION A(NI,NJ),B(NK,NI),C(NK,NJ)
C
c     CALL MXMB(B,NK,1,C,1,NK,A,1,NI,NI,NK,NJ)
C
      DO 3 I=1,NI
           DO 2 J=1,NJ
                T=A(I,J)
                DO 1 K=1,NK
                     T=T+B(K,I)*C(K,J)
    1           CONTINUE
                A(I,J)=T
    2      CONTINUE
    3 CONTINUE
      RETURN
C     END
C
C     DIMENSION A(1),B(1),C(1)
C
C     IJ=1
C     JA=1
C     DO 3 J=1,NJ
C     KI=1
C     DO 2 I=1,NI
C     KJ=JA
C     T=A(IJ)
C     DO 1 K=1,NK
C     T=T+B(KI)*C(KJ)
C     KI=KI+1
C     KJ=KJ+1
C   1 CONTINUE
C     A(IJ)=T
C     IJ=IJ+1
C   2 CONTINUE
C     JA=JA+NK
C   3 CONTINUE
C     RETURN
      END
      SUBROUTINE APBCT(A,B,C,NI,NK,NJ)
C
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION A(NI,NJ),B(NI,NK),C(NJ,NK)
C
c     CALL MXMB(B,1,NI,C,NJ,1,A,1,NI,NI,NK,NJ)
C
      DO 3 I=1,NI
           DO 2 J=1,NJ
                T=A(I,J)
                DO 1 K=1,NK
                     T=T+B(I,K)*C(J,K)
    1           CONTINUE
                A(I,J)=T
    2      CONTINUE
    3 CONTINUE
      RETURN
C     END
C
C     DIMENSION A(1),B(1),C(1)
C
C     IJ=1
C     DO 3 J=1,NJ
C     DO 2 I=1,NI
C     IK=I
C     JK=J
C     T=A(IJ)
C     DO 1 K=1,NK
C     T=T+B(IK)*C(JK)
C     IK=IK+NI
C     JK=JK+NJ
C   1 CONTINUE
C     A(IJ)=T
C     IJ=IJ+1
C   2 CONTINUE
C   3 CONTINUE
C     RETURN
      END
      SUBROUTINE APBTCT(A,B,C,NI,NK,NJ)
C
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION A(NI,NJ),B(NK,NI),C(NJ,NK)
C
c     CALL MXMB(B,NK,1,C,NJ,1,A,1,NI,NI,NK,NJ)
C
      DO 3 I=1,NI
           DO 2 J=1,NJ
                T=A(I,J)
                DO 1 K=1,NK
                     T=T+B(K,I)*C(J,K)
    1           CONTINUE
                A(I,J)=T
    2      CONTINUE
    3 CONTINUE
      RETURN
C     END
C
C     DIMENSION A(1),B(1),C(1)
C
C     IJ=1
C     DO 3 J=1,NJ
C     KI=1
C     DO 2 I=1,NI
C     JK=J
C     T=A(IJ)
C     DO 1 K=1,NK
C     T=T+B(KI)*C(JK)
C     JK=JK+NJ
C     KI=KI+1
C   1 CONTINUE
C     A(IJ)=T
C     IJ=IJ+1
C   2 CONTINUE
C   3 CONTINUE
C     RETURN
      END
      SUBROUTINE AMBC(A,B,C,NI,NK,NJ)
C
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION A(NI,NJ),B(NI,NK),C(NK,NJ)
C
      DO 3 I=1,NI
           DO 2 J=1,NJ
                T=A(I,J)
                DO 1 K=1,NK
                     T=T-B(I,K)*C(K,J)
    1           CONTINUE
                A(I,J)=T
    2      CONTINUE
    3 CONTINUE
      RETURN
C     END
C
C     DIMENSION A(1),B(1),C(1)
C
C     IJ=1
C     JA=1
C     DO 3 J=1,NJ
C     DO 2 I=1,NI
C     IK=I
C     KJ=JA
C     T=A(IJ)
C     DO 1 K=1,NK
C     T=T-B(IK)*C(KJ)
C     IK=IK+NI
C     KJ=KJ+1
C   1 CONTINUE
C     A(IJ)=T
C     IJ=IJ+1
C   2 CONTINUE
C     JA=JA+NK
C   3 CONTINUE
C     RETURN
      END
      SUBROUTINE AMBTC(A,B,C,NI,NK,NJ)
C
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION A(NI,NJ),B(NK,NI),C(NK,NJ)
C
      DO 3 I=1,NI
           DO 2 J=1,NJ
                T=A(I,J)
                DO 1 K=1,NK
                     T=T-B(K,I)*C(K,J)
    1           CONTINUE
                A(I,J)=T
    2      CONTINUE
    3 CONTINUE
      RETURN
C     END
C
C     DIMENSION A(1),B(1),C(1)
C
C     IJ=1
C     JA=1
C     DO 3 J=1,NJ
C     KI=1
C     DO 2 I=1,NI
C     KJ=JA
C     T=A(IJ)
C     DO 1 K=1,NK
C     T=T-B(KI)*C(KJ)
C     KI=KI+1
C     KJ=KJ+1
C   1 CONTINUE
C     A(IJ)=T
C     IJ=IJ+1
C   2 CONTINUE
C     JA=JA+NK
C   3 CONTINUE
C     RETURN
      END
      SUBROUTINE AMBCT(A,B,C,NI,NK,NJ)
C
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION A(NI,NJ),B(NI,NK),C(NJ,NK)
C
      DO 3 I=1,NI
           DO 2 J=1,NJ
                T=A(I,J)
                DO 1 K=1,NK
                     T=T-B(I,K)*C(J,K)
    1           CONTINUE
                A(I,J)=T
    2      CONTINUE
    3 CONTINUE
      RETURN
C     END
C
C     DIMENSION A(1),B(1),C(1)
C
C     IJ=1
C     DO 3 J=1,NJ
C     DO 2 I=1,NI
C     IK=I
C     JK=J
C     T=A(IJ)
C     DO 1 K=1,NK
C     T=T-B(IK)*C(JK)
C     IK=IK+NI
C     JK=JK+NJ
C   1 CONTINUE
C     A(IJ)=T
C     IJ=IJ+1
C   2 CONTINUE
C   3 CONTINUE
C     RETURN
      END
      SUBROUTINE AMBTCT(A,B,C,NI,NK,NJ)
C
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION A(NI,NJ),B(NK,NI),C(NJ,NK)
C
      DO 3 I=1,NI
           DO 2 J=1,NJ
                T=A(I,J)
                DO 1 K=1,NK
                     T=T-B(K,I)*C(J,K)
    1           CONTINUE
                A(I,J)=T
    2      CONTINUE
    3 CONTINUE
      RETURN
C     END
C
C     DIMENSION A(1),B(1),C(1)
C
C     IJ=1
C     DO 3 J=1,NJ
C     KI=1
C     DO 2 I=1,NI
C     JK=J
C     T=A(IJ)
C     DO 1 K=1,NK
C     T=T-B(KI)*C(JK)
C     JK=JK+NJ
C     KI=KI+1
C   1 CONTINUE
C     A(IJ)=T
C     IJ=IJ+1
C   2 CONTINUE
C   3 CONTINUE
C     RETURN
      END
CBHL      FUNCTION DOT(A,B,N)
C     REAL*8 FUNCTION DOT(A,B,N)
Cibm
      DOUBLE PRECISION FUNCTION DOT(A,B,N)
Cibm
C
      REAL*8 A(N),B(N),T
C
      T=0.0D+00
      DO 1 I=1,N
           T=T+A(I)*B(I)
    1 CONTINUE
      DOT=T
      RETURN
      END
      SUBROUTINE ATEBC(A,B,C,NI,NK,NJ)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      REAL*8 A(NJ,NI),B(NI,NK),C(NK,NJ)
C
      CALL ZERO(A,NJ*NI)
C
c     CALL MXMB(B,1,NI,C,1,NK,A,NJ,1,NI,NK,NJ)
C
      DO 3 J=1,NJ
         DO 2 K=1,NK
            CKJ=C(K,J)
            IF (CKJ.EQ.0.0D+00) GO TO 2
            DO 1 I=1,NI
               A(J,I)=A(J,I)+B(I,K)*CKJ
    1       CONTINUE
    2    CONTINUE
    3 CONTINUE
C
      RETURN
      END
      SUBROUTINE ATEMBC(A,B,C,NI,NK,NJ)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      REAL*8 A(NJ,NI),B(NI,NK),C(NK,NJ)
C
      CALL ZERO(A,NJ*NI)
C
      DO 3 J=1,NJ
         DO 2 K=1,NK
            CKJ=C(K,J)
            IF (CKJ.EQ.0.0D+00) GO TO 2
            DO 1 I=1,NI
               A(J,I)=A(J,I)-B(I,K)*CKJ
    1       CONTINUE
    2    CONTINUE
    3 CONTINUE
C
      RETURN
      END
      SUBROUTINE ATPBC(A,B,C,NI,NK,NJ)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      REAL*8 A(NJ,NI),B(NI,NK),C(NK,NJ)
C
c      CALL MXMB(B,1,NI,C,1,NK,A,NJ,1,NI,NK,NJ)
C
      DO 3 J=1,NJ
         DO 2 K=1,NK
            CKJ=C(K,J)
            IF (CKJ.EQ.0.0D+00) GO TO 2
            DO 1 I=1,NI
               A(J,I)=A(J,I)+B(I,K)*CKJ
    1       CONTINUE
    2    CONTINUE
    3 CONTINUE
C
      RETURN
      END
      SUBROUTINE ATMBC(A,B,C,NI,NK,NJ)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      REAL*8 A(NJ,NI),B(NI,NK),C(NK,NJ)
C
      DO 3 J=1,NJ
         DO 2 K=1,NK
            CKJ=C(K,J)
            IF (CKJ.EQ.0.0D+00) GO TO 2
            DO 1 I=1,NI
               A(J,I)=A(J,I)-B(I,K)*CKJ
    1       CONTINUE
    2    CONTINUE
    3 CONTINUE
C
      RETURN
      END
      SUBROUTINE ATPBCT(A,B,C,NI,NK,NJ)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      REAL*8 A(NJ,NI),B(NI,NK),C(NJ,NK)
C
c     CALL MXMB(B,1,NI,C,NJ,1,A,NJ,1,NI,NK,NJ)
C
      DO 3 K=1,NK
         DO 2 J=1,NJ
            CJK=C(J,K)
            IF (CJK.EQ.0.0D+00) GO TO 2
            DO 1 I=1,NI
               A(J,I)=A(J,I)+B(I,K)*CJK
    1       CONTINUE
    2    CONTINUE
    3 CONTINUE
C
      RETURN
      END
      SUBROUTINE APBCTT(A,B,C,NI,NK,NJ)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      REAL*8 A(NJ,NI),B(NK,NI),C(NJ,NK)
C
CTJL  CALL MXMB(B,1,NI,C,NK,1,A,NI,1,NI,NK,NJ)
C
      DO 3 K=1,NK
         DO 2 J=1,NJ
            CJK=C(J,K)
            IF (CJK.EQ.0.0D+00) GO TO 2
            DO 1 I=1,NI
               A(J,I)=A(J,I)+B(K,I)*CJK
    1       CONTINUE
    2    CONTINUE
    3 CONTINUE
C
      RETURN
      END
      SUBROUTINE IZERO(A,N)
C
      INTEGER A(N)
C
      NBYTES = 4 * N
c      CALL VXINIT(A,0,NBYTES)
       DO 1 I=1,N
            A(I)=0
   1   CONTINUE
C
      RETURN
      END


      SUBROUTINE GIVENS (NX,NROOTX,NJX,A,B,ROOT,VECT)
C     DOUBLE PRECISION VERSION BY MEL LEVY 8/72
C 62.3  GIVENS  -EIGENVALUES AND EIGENVECTORS BY THE GIVENS METHOD.
C     BY FRANKLIN PROSSER, INDIANA UNIVERSITY.
C     SEPTEMBER, 1967
C     CALCULATES EIGENVALUES AND EIGENVECTORS OF REAL SYMMETRIC MATRIX
C     STORED IN PACKED UPPER TRIANGULAR FORM.
C
C     THANKS ARE DUE TO F. E. HARRIS (STANFORD UNIVERSITY) AND H. H.
C     MICHELS (UNITED AIRCRAFT RESEARCH LABORATORIES) FOR EXCELLENT
C     WORK ON NUMERICAL DIFFICULTIES WITH EARLIER VERSIONS OF THIS
C     PROGRAM.
C
C     THE PARAMETERS FOR THE ROUTINE ARE...
C         NX     ORDER OF MATRIX
C         NROOTX NUMBER OF ROOTS WANTED.  THE NROOTX SMALLEST (MOST
C                 NEGATIVE) ROOTS WILL BE CALCULATED.  IF NO VECTORS
C                 ARE WANTED, MAKE THIS NUMBER NEGATIVE.
C         NJX    ROW DIMENSION OF VECT ARRAY.  SEE -VECT- BELOW.
C                 NJX MUST BE NOT LESS THAN NX.
C         A      MATRIX STORED BY COLUMNS IN PACKED UPPER TRIANGULAR
C                FORM, I.E. OCCUPYING NX*(NX+1)/2 CONSECUTIVE
C                LOCATIONS.
C         B      SCRATCH ARRAY USED BY GIVENS.  MUST BE AT LEAST
C                 NX*5 CELLS.
C         ROOT   ARRAY TO HOLD THE EIGENVALUES.  MUST BE AT LEAST
C                NROOTX CELLS LONG.  THE NROOTX SMALLEST ROOTS ARE
C                 ORDERED LARGEST FIRST IN THIS ARRAY.
C         VECT   EIGENVECTOR ARRAY.  EACH COLUMN WILL HOLD AN
C                 EIGENVECTOR FOR THE CORRESPONDING ROOT.  MUST BE
C                 DIMENSIONED WITH -NJX- ROWS AND AT LEAST -NROOTX-
C                 COLUMNS, UNLESS NO VECTORS
C                 ARE REQUESTED (NEGATIVE NROOTX).  IN THIS LATTER
C                 CASE, THE ARGUMENT VECT IS JUST A DUMMY, AND THE
C                 STORAGE IS NOT USED.
C                 THE EIGENVECTORS ARE NORMALIZED TO UNITY.
C
C     THE ARRAYS A AND B ARE DESTROYED BY THE COMPUTATION.  THE RESULTS
C     APPEAR IN ROOT AND VECT.
C     FOR PROPER FUNCTIONING OF THIS ROUTINE, THE RESULT OF A FLOATING
C     POINT UNDERFLOW SHOULD BE A ZERO.
C     TO CONVERT THIS ROUTINE TO DOUBLE PRECISION (E.G. ON IBM 360
C     MACHINES), BE SURE THAT ALL REAL VARIABLES AND FUNCTION
C     REFERENCES ARE PROPERLY MADE DOUBLE PRECISION.
C     THE VALUE OF -ETA- (SEE BELOW) SHOULD ALSO BE CHANGED, TO REFLECT
C     THE INCREASED PRECISION.
C
C     THE ORIGINAL REFERENCE TO THE GIVENS TECHNIQUE IS IN OAK RIDGE
C     REPORT NUMBER ORNL 1574 (PHYSICS), BY WALLACE GIVENS.
C     THE METHOD AS PRESENTED IN THIS PROGRAM CONSISTS OF FOUR STEPS,
C     ALL MODIFICATIONS OF THE ORIGINAL METHOD...
C     FIRST, THE INPUT MATRIX IS REDUCED TO TRIDIAGONAL FORM BY THE
C     HOUSEHOLDER TECHNIQUE (J. H. WILKINSON, COMP. J. 3, 23 (1960)).
C     THE ROOTS ARE THEN LOCATED BY THE STURM SEQUENCE METHOD (J. M.
C     ORTEGA (SEE REFERENCE BELOW).  THE VECTORS OF THE TRIDIAGONAL
C     FORM ARE THEN EVALUATED (J. H. WILKINSON, COMP. J. 1, 90 (1958)),
C     AND LAST THE TRIDIAGONAL VECTORS ARE ROTATED TO VECTORS OF THE
C     ORIGINAL ARRAY (FIRST REFERENCE).
C     VECTORS FOR DEGENERATE (OR NEAR-DEGENERATE) ROOTS ARE FORCED
C     TO BE ORTHOGONAL, USING A METHOD SUGGESTED BY B. GARBOW, ARGONNE
C     NATIONAL LABS (PRIVATE COMMUNICATION, 1964).  THE GRAM-SCHMIDT
C     PROCESS IS USED FOR THE ORTHOGONALIZATION.
C
C     AN EXCELLENT PRESENTATION OF THE GIVENS TECHNIQUE IS FOUND IN
C     J. M. ORTEGA-S ARTICLE IN -MATHEMATICS FOR DIGITAL COMPUTERS,-
C     VOLUME 2, ED. BY RALSTON AND WILF, WILEY (1967), PAGE 94.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DOUBLE PRECISION B(NX,5),A(1),ROOT(NROOTX),VECT(NJX,NROOTX)
      DOUBLE PRECISION ETA,THETA,DEL1,DELTA,SMALL,DELBIG,THETA1,TOLER
      DOUBLE PRECISION RPOWER,RPOW1,RAND1,FACTOR,ANORM,ALIMIT,SUM,TEMP
      DOUBLE PRECISION AK,ROOTL,ROOTX,TRIAL,F0,SAVE,AROOT
      DOUBLE PRECISION ELIM1,ELIM2
C
C ** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C **  USERS PLEASE NOTE...
C **  THE FOLLOWING TWO PARAMETERS, ETA AND THETA, SHOULD BE ADJUSTED
C **  BY THE USER FOR HIS PARTICULAR MACHINE.
C **  ETA IS AN INDICATION OF THE PRECISION OF THE FLOATING POINT
C **  REPRESENTATION ON THE COMPUTER BEING USED (ROUGHLY 10**(-M),
C **  WHERE M IS THE NUMBER OF DECIMALS OF PRECISION ).
C **  THETA IS AN INDICATION OF THE RANGE OF NUMBERS THAT CAN BE
C **  EXPRESSED IN THE FLOATING POINT REPRESENTATION (ROUGHLY THE
C **  LARGEST NUMBER).
C **  SOME RECOMMENDED VALUES FOLLOW.
C **  FOR IBM 7094, UNIVAC 1108, ETC. (27-BIT BINARY FRACTION, 8-BIT
C **  BINARY EXPONENT), ETA=1.E-8, THETA=1.E37.
C **  FOR CONTROL DATA 3600 (36-BIT BINARY FRACTION, 11-BIT BINARY
C **  EXPONENT), ETA=1.E-11, THETA=1.E307.
C **  FOR CONTROL DATA 6600 (48-BIT BINARY FRACTION, 11-BIT BINARY
C **  EXPONENT), ETA=1.E-14, THETA=1.E307.
C **  FOR IBM 360/50 AND 360/65 DOUBLE PRECISION (56-BIT HEXADECIMAL
C **  FRACTION, 7-BIT HEXADECIMAL EXPONENT), ETA=1.E-16, THETA=1.E75.
C **
      THETA = 1.D75
      ETA = 1.D-16
C ** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
      DEL1=ETA/100.0D0
      DELTA=ETA**2*100.0D0
      SMALL=ETA**2/100.0D0
      DELBIG=THETA*DELTA/1000.0D0
      THETA1=1000.0D0/THETA
C     TOLER  IS A FACTOR USED TO DETERMINE IF TWO ROOTS ARE CLOSE
C     ENOUGH TO BE CONSIDERED DEGENERATE FOR PURPOSES OF ORTHOGONALI-
C     ZING THEIR VECTORS.  FOR THE MATRIX NORMED TO UNITY, IF THE
C     DIFFERENCE BETWEEN TWO ROOTS IS LESS THAN TOLER, THEN
C     ORTHOGONALIZATION WILL OCCUR.
      TOLER = ETA*100.
C
C     INITIAL VALUE FOR PSEUDORANDOM NUMBER GENERATOR... (2**23)-3
      RPOWER=8388608.0D0
      RPOW1 = RPOWER/2.
      RAND1=RPOWER-3.0D0
C
      N = NX
      NROOT = IABS(NROOTX)
      IF (NROOT.EQ.0) GO TO 1001
      IF (N-1) 1001,1003,105
1003  ROOT(1) = A(1)
      IF(NROOTX .GT. 0)VECT(1,1)=1.0D0
      GO TO 1001
105   CONTINUE
C     NSIZE    NUMBER OF ELEMENTS IN THE PACKED ARRAY
      NSIZE = (N*(N+1))/2
      NM1 = N-1
      NM2 = N-2
C
C     SCALE MATRIX TO EUCLIDEAN NORM OF 1.  SCALE FACTOR IS ANORM.
      FACTOR=0.0D-0
      DO 70 I=1,NSIZE
 70   FACTOR=DMAX1(FACTOR,DABS(A(I)))
      IF(FACTOR .NE. 0.0D0)GO TO 72
C     NULL MATRIX.  FIX UP ROOTS AND VECTORS, THEN EXIT.
      DO 78 I=1,NROOT
           IF (NROOTX.LT.0) GO TO 78
           DO 77 J=1,N
 77        VECT(J,I)=0.0D0
           VECT(I,I)=1.0D0
 78   ROOT(I)=0.0D0
      GO TO 1001
C
 72   ANORM=0.0D0
      J = 1
      K = 1
      DO 80 I=1,NSIZE
           IF (I.NE.J) GO TO 81
           ANORM=ANORM+(A(I)/FACTOR)**2/2.0D0
           K = K+1
           J = J+K
           GO TO 80
81         ANORM = ANORM + (A(I)/FACTOR)**2
80    CONTINUE
      ANORM=DSQRT(ANORM*2.0D0)*FACTOR
      DO 91 I=1,NSIZE
91    A(I) = A(I)/ANORM
      ALIMIT=1.0D-0
C
C     TRIDIA SECTION.
C     TRIDIAGONALIZATION OF SYMMETRIC MATRIX
      ID = 0
      IA = 1
      IF (NM2.EQ.0) GO TO 201
      DO 200  J=1,NM2
C     J       COUNTS ROW  OF A-MATRIX TO BE DIAGONALIZED
C     IA      START OF NON-CODIAGONAL ELEMENTS IN THE ROW
C     ID      INDEX OF CODIAGONAL ELEMENT ON ROW BEING CODIAGONALIZED.
           IA = IA+J+2
           ID = ID + J + 1
           JP2 = J+2
C     SUM SQUARES OF NON-CODIAGONAL ELEMENTS IN ROW J
           II = IA
           SUM=0.0D-0
           DO 100 I=JP2,N
                SUM = SUM + A(II)**2
100        II = II + I
           TEMP = A(ID)
           IF (SUM.GT.SMALL) GO TO 110
C     NO TRANSFORMATION NECESSARY IF ALL THE NON-CODIAGONAL
C     ELEMENTS ARE TINY.
           B(J,1) = TEMP
           A(ID)=0.0D0
           GO TO 200
C     NOW COMPLETE THE SUM OF OFF-DIAGONAL SQUARES
 110       SUM=DSQRT(SUM+TEMP**2)
C     NEW CODIAGONAL ELEMENT
           B(J,1)=-DSIGN(SUM,TEMP)
C     FIRST NON-ZERO ELEMENT OF THIS W-VECTOR
           B(J+1,2)=DSQRT((1.0D0+DABS(TEMP)/SUM)/2.0D0)
C     FORM REST OF THE W-VECTOR ELEMENTS
           TEMP=DSIGN(0.5D0/(B(J+1,2)*SUM),TEMP)
           II = IA
           DO 130 I=JP2,N
                B(I,2) = A(II)*TEMP
130        II = II + I
C     FORM P-VECTOR AND SCALAR.  P-VECTOR = A-MATRIX*W-VECTOR.
C     SCALAR = W-VECTOR*P-VECTOR.
           AK=0.0D0
C     IC      LOCATION OF NEXT DIAGONAL ELEMENT
           IC = ID + 1
           J1 = J + 1
           DO 190  I=J1,N
                JJ = IC
                TEMP=0.0D0
                DO 180  II=J1,N
C     I       RUNS OVER THE NON-ZERO P-ELEMENTS
C     II      RUNS OVER ELEMENTS OF W-VECTOR
                     TEMP = TEMP + B(II,2)*A(JJ)
C     CHANGE INCREMENTING MODE AT THE DIAGONAL ELEMENTS.
                     IF (II.LT.I) GO TO 210
                     JJ = JJ + II
                     GO TO 180
210                  JJ = JJ + 1
180             CONTINUE
C     BUILD UP THE K-SCALAR (AK)
                AK = AK + TEMP*B(I,2)
                B(I,1) = TEMP
C     MOVE IC TO TOP OF NEXT A-MATRIX -ROW-
190        IC = IC + I
C     FORM THE Q-VECTOR
           DO 150  I=J1,N
150        B(I,1) = B(I,1) - AK*B(I,2)
C     TRANSFORM THE REST OF THE A-MATRIX
C     JJ      START-1 OF THE REST OF THE A-MATRIX
           JJ = ID
C     MOVE W-VECTOR INTO THE OLD A-MATRIX LOCATIONS TO SAVE SPACE
C     I       RUNS OVER THE SIGNIFICANT ELEMENTS OF THE W-VECTOR
           DO 160  I=J1,N
                A(JJ) = B(I,2)
                DO 170  II=J1,I
                     JJ = JJ + 1
 170            A(JJ)=A(JJ)-2.0D0*(B(I,1)*B(II,2)+B(I,2)*B(II,1))
160        JJ = JJ + J
200   CONTINUE
C     MOVE LAST CODIAGONAL ELEMENT OUT INTO ITS PROPER PLACE
201   CONTINUE
      B(NM1,1) = A(NSIZE-1)
      A(NSIZE-1)=0.0D-0
C
C     STURM SECTION.
C     STURM SEQUENCE ITERATION TO OBTAIN ROOTS OF TRIDIAGONAL FORM.
C     MOVE DIAGONAL ELEMENTS INTO SECOND N ELEMENTS OF B-VECTOR.
C     THIS IS A MORE CONVENIENT INDEXING POSITION.
C     ALSO, PUT SQUARE OF CODIAGONAL ELEMENTS IN THIRD N ELEMENTS.
      JUMP=1
      DO 320 J=1,N
           B(J,2)=A(JUMP)
           B(J,3) = B(J,1)**2
320   JUMP = JUMP+J+1
      DO 310 I=1,NROOT
310   ROOT(I) = +ALIMIT
      ROOTL = -ALIMIT
C     ISOLATE THE ROOTS.  THE NROOT LOWEST ROOTS ARE FOUND, LOWEST FIRST
      DO 330 I=1,NROOT
C     FIND CURRENT -BEST- UPPER BOUND
           ROOTX = +ALIMIT
           DO 335 J=I,NROOT
 335       ROOTX=DMIN1(ROOTX,ROOT(J))
           ROOT(I) = ROOTX
C     GET IMPROVED TRIAL ROOT
 500       TRIAL=(ROOTL+ROOT(I))*0.5D0
           IF (TRIAL.EQ.ROOTL.OR.TRIAL.EQ.ROOT(I)) GO TO 330
C     FORM STURM SEQUENCE RATIOS, USING ORTEGA-S ALGORITHM (MODIFIED).
C     NOMTCH IS THE NUMBER OF ROOTS LESS THAN THE TRIAL VALUE.
           NOMTCH=N
           J=1
360        F0 = B(J,2) - TRIAL
370        CONTINUE
           IF(DABS(F0) .LT. THETA1)GO TO 380
           IF(F0 .GE. 0.0D0)NOMTCH=NOMTCH-1
           J = J + 1
           IF (J.GT.N) GO TO 390
C     SINCE MATRIX IS NORMED TO UNITY, MAGNITUDE OF B(J,3) IS LESS THAN
C     ONE, AND SINCE F0 IS GREATER THAN THETA1, OVERFLOW IS NOT POSSIBLE
C     AT THE DIVISION STEP.
           F0 = B(J,2) - TRIAL - B(J-1,3)/F0
           GO TO 370
380        J = J + 2
           NOMTCH = NOMTCH - 1
           IF (J.LE.N) GO TO 360
390        CONTINUE
C     FIX NEW BOUNDS ON ROOTS
           IF (NOMTCH.GE.I) GO TO 540
           ROOTL = TRIAL
           GO TO 500
540        ROOT(I) = TRIAL
           NOM = MIN0(NROOT,NOMTCH)
           ROOT(NOM) = TRIAL
           GO TO 500
330   CONTINUE
C     REVERSE THE ORDER OF THE EIGENVALUES, SINCE CUSTOM DICTATES
C     -LARGEST FIRST-.  THIS SECTION MAY BE REMOVED IF DESIRED WITHOUT
C     AFFECTING THE REMAINDER OF THE ROUTINE.
C     NRT = NROOT/2
C     DO 10 I=1,NRT
C     SAVE = ROOT(I)
C     NMIP1 = NROOT - I + 1
CC    ROOT(I) = ROOT(NMIP1)
C10   ROOT(NMIP1) = SAVE
C
C     TRIVEC SECTION.
C     EIGENVECTORS OF CODIAGONAL FORM
C807  CONTINUE
C     QUIT NOW IF NO VECTORS WERE REQUESTED.
      IF (NROOTX.LT.0) GO TO 1002
C     INITIALIZE VECTOR ARRAY.
      DO 15 I=1,N
           DO 15 J=1,NROOT
 15   VECT(I,J)=1.0D-0
      DO 700 I=1,NROOT
           AROOT = ROOT(I)
C     ORTHOGONALIZE IF ROOTS ARE CLOSE.
           IF (I.EQ.1) GO TO 710
C     THE ABSOLUTE VALUE IN THE NEXT TEST IS TO ASSURE THAT THE TRIVEC
C     SECTION IS INDEPENDENT OF THE ORDER OF THE EIGENVALUES.
           IF(DABS(ROOT(I-1)-AROOT) .LT. TOLER)GO TO 720
710        IA = -1
720        IA = IA + 1
           ELIM1 = A(1) - AROOT
           ELIM2 = B(1,1)
           JUMP = 1
           DO 750  J=1,NM1
                JUMP = JUMP+J+1
C     GET THE CORRECT PIVOT EQUATION FOR THIS STEP.
                IF(DABS(ELIM1) .LE. DABS(B(J,1)))GO TO 760
C     FIRST (ELIM1) EQUATION IS THE PIVOT THIS TIME.  CASE 1.
                B(J,2) = ELIM1
                B(J,3) = ELIM2
                B(J,4)=0.0D0
                TEMP = B(J,1)/ELIM1
                ELIM1 = A(JUMP) - AROOT - TEMP*ELIM2
                ELIM2 = B(J+1,1)
                GO TO 755
C     SECOND EQUATION IS THE PIVOT THIS TIME.  CASE 2.
760             B(J,2) = B(J,1)
                B(J,3) = A(JUMP) - AROOT
                B(J,4) = B(J+1,1)
                TEMP=1.0D0
                IF(DABS(B(J,1)) .GT. THETA1)TEMP=ELIM1/B(J,1)
                ELIM1 = ELIM2 - TEMP*B(J,3)
                ELIM2 = -TEMP*B(J+1,1)
C     SAVE FACTOR FOR THE SECOND ITERATION.
755             B(J,5) = TEMP
750        CONTINUE
           B(N,2) = ELIM1
           B(N,3)=0.0D-0
           B(N,4)=0.0D-0
           B(NM1,4)=0.0D-0
           ITER = 1
           IF (IA.NE.0) GO TO 801
C     BACK SUBSTITUTE TO GET THIS VECTOR.
790        L = N + 1
           DO 780 J=1,N
                L = L - 1
786             CONTINUE
                ELIM1=VECT(L,I)-VECT(L+1,I)*B(L,3)-VECT(L+2,I)*B(L,4)
C     IF OVERFLOW IS CONCEIVABLE, SCALE THE VECTOR DOWN.
C     THIS APPROACH IS USED TO AVOID MACHINE-DEPENDENT AND SYSTEM-
C     DEPENDENT CALLS TO OVERFLOW ROUTINES.
                IF(DABS(ELIM1) .GT. DELBIG)GO TO 782
                TEMP = B(L,2)
                IF(DABS(B(L,2)) .LT. DELTA)TEMP=DELTA
                VECT(L,I) = ELIM1/TEMP
                GO TO 780
C     VECTOR IS TOO BIG.  SCALE IT DOWN.
782             DO 784 K=1,N
784             VECT(K,I) = VECT(K,I)/DELBIG
                GO TO 786
780        CONTINUE
           GO TO (820,800), ITER
C     SECOND ITERATION.  (BOTH ITERATIONS FOR REPEATED-ROOT VECTORS).
820        ITER = ITER + 1
890        ELIM1 = VECT(1,I)
           DO 830 J=1,NM1
                IF (B(J,2).EQ.B(J,1)) GO TO 840
C     CASE ONE.
                VECT(J,I) = ELIM1
                ELIM1 = VECT(J+1,I) - ELIM1*B(J,5)
                GO TO 830
C     CASE TWO.
840             VECT(J,I) = VECT(J+1,I)
                ELIM1 = ELIM1 - VECT(J+1,I)*TEMP
830        CONTINUE
           VECT(N,I) = ELIM1
           GO TO 790
C     PRODUCE A RANDOM VECTOR
801        CONTINUE
           DO 802 J=1,N
C     GENERATE PSEUDORANDOM NUMBERS WITH UNIFORM DISTRIBUTION IN (-1,1).
C     THIS RANDOM NUMBER SCHEME IS OF THE FORM...
C     RAND1 = AMOD((2**12+3)*RAND1,2**23)
C     IT HAS A PERIOD OF 2**21 NUMBERS.
                RAND1=DMOD(4099.0D0*RAND1,RPOWER)
 802       VECT(J,I)=RAND1/RPOW1-1.0D0
           GO TO 790
C
C     ORTHOGONALIZE THIS REPEATED-ROOT VECTOR TO OTHERS WITH THIS ROOT.
800        IF (IA.EQ.0) GO TO 885
           DO 860 J1=1,IA
                K = I - J1
                TEMP=0.0D0
                DO 870 J=1,N
870             TEMP = TEMP + VECT(J,I)*VECT(J,K)
                DO 880 J=1,N
880             VECT(J,I) = VECT(J,I) - TEMP*VECT(J,K)
860        CONTINUE
885        GO TO (890,900), ITER
C     NORMALIZE THE VECTOR
 900       ELIM1=0.0D0
           DO 904 J=1,N
 904       ELIM1=DMAX1(DABS(VECT(J,I)),ELIM1)
           TEMP=0.0D-0
           DO 910 J=1,N
                ELIM2=VECT(J,I)/ELIM1
910        TEMP = TEMP + ELIM2**2
           TEMP=1.0D0/(DSQRT(TEMP)*ELIM1)
           DO 920 J=1,N
                VECT(J,I) = VECT(J,I)*TEMP
                IF(DABS(VECT(J,I)) .LT. DEL1)VECT(J,I)=0.0D0
920        CONTINUE
700   CONTINUE
C
C     SIMVEC SECTION.
C     ROTATE CODIAGONAL VECTORS INTO VECTORS OF ORIGINAL ARRAY
C     LOOP OVER ALL THE TRANSFORMATION VECTORS
      IF (NM2.EQ.0) GO TO 1002
      JUMP = NSIZE - (N+1)
      IM = NM1
      DO 950  I=1,NM2
           J1 = JUMP
C     MOVE A TRANSFORMATION VECTOR OUT INTO BETTER INDEXING POSITION.
           DO 955  J=IM,N
                B(J,2) = A(J1)
955        J1 = J1 + J
C     MODIFY ALL REQUESTED VECTORS.
           DO 960  K=1,NROOT
                TEMP=0.0D0
C     FORM SCALAR PRODUCT OF TRANSFORMATION VECTOR WITH EIGENVECTOR
                DO 970  J=IM,N
970             TEMP = TEMP + B(J,2)*VECT(J,K)
                TEMP = TEMP + TEMP
                DO 980  J=IM,N
980             VECT(J,K) = VECT(J,K) - TEMP*B(J,2)
960        CONTINUE
           JUMP = JUMP - IM
950   IM = IM - 1
1002  CONTINUE
C     RESTORE ROOTS TO THEIR PROPER SIZE.
      DO 95 I=1,NROOT
95    ROOT(I) = ROOT(I)*ANORM
1001  RETURN
      END
      SUBROUTINE NCDLBL(ITAPE,NA,NB,NC,ND,NE,RF,RG,RH)
      REAL*8 RF,RG,RH,RN
      DIMENSION NA(26),NB(26),N(112),RN(10)
      EQUIVALENCE(N(57),RN(1))
      L = 26
      DO 10 I=1,26
           L = L + 1
           N(I) = NA(I)
 10   N(L) = NB(I)
      N(53) = NC
      N(54) = ND
      N(55) = NE
      RN(1) = RF
      RN(2) = RG
      RN(3) = RH
      CALL SWRIT(ITAPE,N,112)
      RETURN
C
      ENTRY DCDLBL(ITAPE,NA,NB,NC,ND,NE,RF,RG,RH)
      CALL SREAD(ITAPE,N,112)
      L = 26
      DO 20 I=1,26
           L = L + 1
           NA(I) = N(I)
 20   NB(I) = N(L)
      NC = N(53)
      ND = N(54)
      NE = N(55)
      RF = RN(1)
      RG = RN(2)
      RH = RN(3)
      RETURN
      END
      SUBROUTINE RSP(NM,N,NV,A,W,MATZ,Z,FV1,FV2)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C
      INTEGER I,J,N,NM,NV,IERR,MATZ
      REAL*8 A(NV),W(N),Z(NM,N),FV1(N),FV2(N)
      DATA ZERO,ONE / 0.0D+00,1.0D+00 /
    1 FORMAT(//,2X,' IERR = ',I5//)
C
C     THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF
C     SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)
C     TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)
C     OF A REAL SYMMETRIC PACKED MATRIX.
C
C     ON INPUT-
C
C        NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
C        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C        DIMENSION STATEMENT,
C
C        N  IS THE ORDER OF THE MATRIX  A,
C
C        NV  IS AN INTEGER VARIABLE SET EQUAL TO THE
C        DIMENSION OF THE ARRAY  A  AS SPECIFIED FOR
C        A  IN THE CALLING PROGRAM.  NV  MUST NOT BE
C        LESS THAN  N*(N+1)/2,
C
C        A  CONTAINS THE LOWER TRIANGLE OF THE REAL SYMMETRIC
C        PACKED MATRIX STORED ROW-WISE,
C
C        MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF
C        ONLY EIGENVALUES ARE DESIRED,  OTHERWISE IT IS SET TO
C        ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.
C
C     ON OUTPUT-
C
C        W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER,
C
C        Z  CONTAINS THE EIGENVECTORS IF MATZ IS NOT ZERO,
C
C        IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN
C        ERROR COMPLETION CODE DESCRIBED IN SECTION 2B OF THE
C        DOCUMENTATION.  THE NORMAL COMPLETION CODE IS ZERO,
C
C        FV1  AND  FV2  ARE TEMPORARY STORAGE ARRAYS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
      IF (N .LE. NM) GO TO 5
      IERR = 10 * N
      GO TO 50
    5 IF (NV .GE. (N * (N + 1)) / 2) GO TO 10
      IERR = 20 * N
      GO TO 50
C
   10 CALL  TRED3(N,NV,A,W,FV1,FV2)
      IF (MATZ .NE. 0) GO TO 20
C     ********** FIND EIGENVALUES ONLY **********
      CALL  TQLRAT(N,W,FV2,IERR)
      GO TO 50
C     ********** FIND BOTH EIGENVALUES AND EIGENVECTORS **********
   20 DO 40 I = 1, N
C
           DO 30 J = 1, N
                Z(J,I) =ZERO
   30      CONTINUE
C
           Z(I,I) =ONE
   40 CONTINUE
C
      CALL  TQL2(NM,N,W,FV1,Z,IERR)
      IF (IERR .NE. 0) GO TO 50
      CALL  TRBAK3(NM,N,NV,A,N,Z)
   50 IF(IERR.NE.0) GO TO 60
      RETURN
   60 WRITE(6,1) IERR
      RETURN
C     ********** LAST CARD OF RSP **********
      END
C
C     ------------------------------------------------------------------
C
      SUBROUTINE TQLRAT(N,D,E2,IERR)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C
      INTEGER I,J,L,M,N,II,L1,MML,IERR
      REAL*8 D(N),E2(N)
      REAL*8 B,C,F,G,H,P,R,S,MACHEP
      DATA ZERO,ONE,TWO /0.0D+00,1.0D+00,2.0D+00/
C     REAL*8 SQRT,ABS,SIGN
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQLRAT,
C     ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC
C     TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD.
C
C     ON INPUT-
C
C        N IS THE ORDER OF THE MATRIX,
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,
C
C        E2 CONTAINS THE SQUARES OF THE SUBDIAGONAL ELEMENTS OF THE
C          INPUT MATRIX IN ITS LAST N-1 POSITIONS.  E2(1) IS ARBITRARY.
C
C      ON OUTPUT-
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
C          THE SMALLEST EIGENVALUES,
C
C        E2 HAS BEEN DESTROYED,
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C
C                **********
      MACHEP = TWO**(-47)
C
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E2(I-1) = E2(I)
C
      F =ZERO
      B =ZERO
      E2(N) =ZERO
C
      DO 290 L = 1, N
           J = 0
           H = MACHEP * (DABS(D(L)) + DSQRT(E2(L)))
           IF (B .GT. H) GO TO 105
           B = H
           C = B * B
C     ********** LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT **********
  105      DO 110 M = L, N
                IF (E2(M) .LE. C) GO TO 120
C     ********** E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP **********
  110      CONTINUE
C
  120      IF (M .EQ. L) GO TO 210
  130      IF (J .EQ. 30) GO TO 1000
           J = J + 1
C     ********** FORM SHIFT **********
           L1 = L + 1
           S = DSQRT(E2(L))
           G = D(L)
           P = (D(L1) - G) / (TWO * S)
           R = DSQRT(P*P+ONE)
           D(L) = S / (P + DSIGN(R,P))
           H = G - D(L)
C
           DO 140 I = L1, N
  140      D(I) = D(I) - H
C
           F = F + H
C     ********** RATIONAL QL TRANSFORMATION **********
           G = D(M)
           IF (G .EQ.ZERO) G = B
           H = G
           S =ZERO
           MML = M - L
C     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********
           DO 200 II = 1, MML
                I = M - II
                P = G * H
                R = P + E2(I)
                E2(I+1) = S * R
                S = E2(I) / R
                D(I+1) = H + S * (H + D(I))
                G = D(I) - E2(I) / G
                IF (G .EQ.ZERO) G = B
                H = G * P / R
  200      CONTINUE
C
           E2(L) = S * G
           D(L) = H
C     ********** GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST **********
           IF (H .EQ.ZERO) GO TO 210
           IF (DABS(E2(L)) .LE. DABS(C/H)) GO TO 210
           E2(L) = H * E2(L)
           IF (E2(L) .NE.ZERO) GO TO 130
  210      P = D(L) + F
C     ********** ORDER EIGENVALUES **********
           IF (L .EQ. 1) GO TO 250
C     ********** FOR I=L STEP -1 UNTIL 2 DO -- **********
           DO 230 II = 2, L
                I = L + 2 - II
                IF (P .GE. D(I-1)) GO TO 270
                D(I) = D(I-1)
  230      CONTINUE
C
  250      I = 1
  270      D(I) = P
  290 CONTINUE
C
      GO TO 1001
C     ********** SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS **********
 1000 IERR = L
 1001 RETURN
C     ********** LAST CARD OF TQLRAT **********
      END
C
C     ------------------------------------------------------------------
C
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C
      INTEGER I,J,K,L,M,N,II,L1,NM,MML,IERR
      REAL*8 D(N),E(N),Z(NM,N)
      REAL*8 B,C,F,G,H,P,R,S,MACHEP
      DATA ZERO,ONE,TWO /0.0D+00,1.0D+00,2.0D+00/
C     REAL*8 SQRT,ABS,SIGN
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2,
C     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND
C     WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
C     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.
C     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO
C     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS
C     FULL MATRIX TO TRIDIAGONAL FORM.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRIX,
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY,
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS
C          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN
C          THE IDENTITY MATRIX.
C
C      ON OUTPUT-
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
C          UNORDERED FOR INDICES 1,2,...,IERR-1,
C
C        E HAS BEEN DESTROYED,
C
C        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC
C          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,
C          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
C          EIGENVALUES,
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
C
C                **********
      MACHEP =TWO**(-47)
C
      IERR = 0
      IF (N .EQ. 1) GO TO 1001
C
      DO 100 I = 2, N
  100 E(I-1) = E(I)
C
      F =ZERO
      B =ZERO
      E(N) = ZERO
C
      DO 240 L = 1, N
           J = 0
           H = MACHEP * (DABS(D(L)) + DABS(E(L)))
           IF (B .LT. H) B = H
C     ********** LOOK FOR SMALL SUB-DIAGONAL ELEMENT **********
           DO 110 M = L, N
                IF (DABS(E(M)) .LE. B) GO TO 120
C     ********** E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP **********
  110      CONTINUE
C
  120      IF (M .EQ. L) GO TO 220
  130      IF (J .EQ. 30) GO TO 1000
           J = J + 1
C     ********** FORM SHIFT **********
           L1 = L + 1
           G = D(L)
           P = (D(L1) - G) / (TWO * E(L))
           R = DSQRT(P*P+ONE)
           D(L) = E(L) / (P + DSIGN(R,P))
           H = G - D(L)
C
           DO 140 I = L1, N
  140      D(I) = D(I) - H
C
           F = F + H
C     ********** QL TRANSFORMATION **********
           P = D(M)
           C = ONE
           S = ZERO
           MML = M - L
C     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********
           DO 200 II = 1, MML
                I = M - II
                G = C * E(I)
                H = C * P
                IF (DABS(P) .LT. DABS(E(I))) GO TO 150
                C = E(I) / P
                R = DSQRT(C*C+ONE)
                E(I+1) = S * P * R
                S = C / R
                C = ONE / R
                GO TO 160
  150           C = P / E(I)
                R = DSQRT(C*C+ONE)
                E(I+1) = S * E(I) * R
                S = ONE / R
                C = C * S
  160           P = C * D(I) - S * G
                D(I+1) = H + S * (C * G + S * D(I))
C     ********** FORM VECTOR **********
                DO 180 K = 1, N
                     H = Z(K,I+1)
                     Z(K,I+1) = S * Z(K,I) + C * H
                     Z(K,I) = C * Z(K,I) - S * H
  180           CONTINUE
C
  200      CONTINUE
C
           E(L) = S * P
           D(L) = C * P
           IF (DABS(E(L)) .GT. B) GO TO 130
  220      D(L) = D(L) + F
  240 CONTINUE
C     ********** ORDER EIGENVALUES AND EIGENVECTORS **********
      DO 300 II = 2, N
           I = II - 1
           K = I
           P = D(I)
C
           DO 260 J = II, N
                IF (D(J) .GE. P) GO TO 260
                K = J
                P = D(J)
  260      CONTINUE
C
           IF (K .EQ. I) GO TO 300
           D(K) = D(I)
           D(I) = P
C
           DO 280 J = 1, N
                P = Z(J,I)
                Z(J,I) = Z(J,K)
                Z(J,K) = P
  280      CONTINUE
C
  300 CONTINUE
C
      GO TO 1001
C     ********** SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS **********
 1000 IERR = L
 1001 RETURN
C     ********** LAST CARD OF TQL2 **********
      END
C
C     ------------------------------------------------------------------
C
      SUBROUTINE TRBAK3(NM,N,NV,A,M,Z)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C
      INTEGER I,J,K,L,M,N,IK,IZ,NM,NV
      REAL*8 A(NV),Z(NM,M)
      REAL*8 H,S
      DATA ZERO /0.0D+00/
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRBAK3,
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A REAL SYMMETRIC
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
C     SYMMETRIC TRIDIAGONAL MATRIX DETERMINED BY  TRED3.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRIX,
C
C        NV MUST BE SET TO THE DIMENSION OF THE ARRAY PARAMETER A
C          AS DECLARED IN THE CALLING PROGRAM DIMENSION STATEMENT,
C
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANSFORMATIONS
C          USED IN THE REDUCTION BY  TRED3  IN ITS FIRST
C          N*(N+1)/2 POSITIONS,
C
C        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED,
C
C        Z CONTAINS THE EIGENVECTORS TO BE BACK TRANSFORMED
C          IN ITS FIRST M COLUMNS.
C
C     ON OUTPUT-
C
C        Z CONTAINS THE TRANSFORMED EIGENVECTORS
C          IN ITS FIRST M COLUMNS.
C
C     NOTE THAT TRBAK3 PRESERVES VECTOR EUCLIDEAN NORMS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
      IF (M .EQ. 0) GO TO 200
      IF (N .EQ. 1) GO TO 200
C
      DO 140 I = 2, N
           L = I - 1
           IZ = (I * L) / 2
           IK = IZ + I
           H = A(IK)
           IF (H .EQ.ZERO) GO TO 140
C
           DO 130 J = 1, M
                S =ZERO
                IK = IZ
C
                DO 110 K = 1, L
                     IK = IK + 1
                     S = S + A(IK) * Z(K,J)
  110           CONTINUE
C     ********** DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW **********
                S = (S / H) / H
                IK = IZ
C
                DO 120 K = 1, L
                     IK = IK + 1
                     Z(K,J) = Z(K,J) - S * A(IK)
  120           CONTINUE
C
  130      CONTINUE
C
  140 CONTINUE
C
  200 RETURN
C     ********** LAST CARD OF TRBAK3 **********
      END
C
C     ------------------------------------------------------------------
C
      SUBROUTINE TRED3(N,NV,A,D,E,E2)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C
      INTEGER I,J,K,L,N,II,IZ,JK,NV
      REAL*8 A(NV),D(N),E(N),E2(N)
      REAL*8 F,G,H,HH,SCALE
      DATA ZERO /0.0D+00/
C     REAL*8 SQRT,ABS,SIGN
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED3,
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX, STORED AS
C     A ONE-DIMENSIONAL ARRAY, TO A SYMMETRIC TRIDIAGONAL MATRIX
C     USING ORTHOGONAL SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT-
C
C        N IS THE ORDER OF THE MATRIX,
C
C        NV MUST BE SET TO THE DIMENSION OF THE ARRAY PARAMETER A
C          AS DECLARED IN THE CALLING PROGRAM DIMENSION STATEMENT,
C
C        A CONTAINS THE LOWER TRIANGLE OF THE REAL SYMMETRIC
C          INPUT MATRIX, STORED ROW-WISE AS A ONE-DIMENSIONAL
C          ARRAY, IN ITS FIRST N*(N+1)/2 POSITIONS.
C
C     ON OUTPUT-
C
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL
C          TRANSFORMATIONS USED IN THE REDUCTION,
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX,
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO,
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     ********** FOR I=N STEP -1 UNTIL 1 DO -- **********
      DO  300 II = 1, N
           I = N + 1 - II
           L = I - 1
           IZ = (I * L) / 2
           H =ZERO
           SCALE = ZERO
           IF (L .LT. 1) GO TO 130
C     ********** SCALE ROW (ALGOL TOL THEN NOT NEEDED) **********
           DO 120 K = 1, L
                IZ = IZ + 1
                D(K) = A(IZ)
                SCALE = SCALE + DABS(D(K))
  120      CONTINUE
C
           IF (SCALE .NE.ZERO) GO TO 140
  130      E(I) =ZERO
           E2(I) = ZERO
           GO TO 290
C
  140      DO 150 K = 1, L
                D(K) = D(K) / SCALE
                H = H + D(K) * D(K)
  150      CONTINUE
C
           E2(I) = SCALE * SCALE * H
           F = D(L)
           G = -DSIGN(DSQRT(H),F)
           E(I) = SCALE * G
           H = H - F * G
           D(L) = F - G
           A(IZ) = SCALE * D(L)
           IF (L .EQ. 1) GO TO 290
           F =ZERO
C
           DO 240 J = 1, L
                G =ZERO
                JK = (J * (J-1)) / 2
C     ********** FORM ELEMENT OF A*U **********
                DO 180 K = 1, L
                     JK = JK + 1
                     IF (K .GT. J) JK = JK + K - 2
                     G = G + A(JK) * D(K)
  180           CONTINUE
C     ********** FORM ELEMENT OF P **********
                E(J) = G / H
                F = F + E(J) * D(J)
  240      CONTINUE
C
           HH = F / (H + H)
           JK = 0
C     ********** FORM REDUCED A **********
           DO 260 J = 1, L
                F = D(J)
                G = E(J) - HH * F
                E(J) = G
C
                DO 260 K = 1, J
                     JK = JK + 1
                     A(JK) = A(JK) - F * E(K) - G * D(K)
  260      CONTINUE
C
  290      D(I) = A(IZ+1)
           A(IZ+1) = SCALE * DSQRT(H)
  300 CONTINUE
C
      RETURN
C     ********** LAST CARD OF TRED3 **********
      END
      SUBROUTINE MATINV(A,B,PIVOT,INDEX,M,N,DETERM,NMAX)
C     MATRIX INVERSION WITH ACCOMPANYING SOLUTION OF LINEAR EQUATIONS
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NMAX,NMAX),B(NMAX,NMAX),PIVOT(NMAX),INDEX(NMAX)
      DATA ZERO,ONE / 0.0D+00 , 1.0D+00 /
C
      DETERM=ONE
      DO 20 I=1,N
      PIVOT(I)=ZERO
   20 CONTINUE
C     PERFORM SUCCESSIVE PIVOT OPERATIONS (GRAND LOOP)
      DO 550 I=1,N
C     SEARCH FOR PIVOT ELEMENT AND EXTEND DETERMINANT PARTIAL PRODUCT
      AMAX=ZERO
      DO 105 J=1,N
      IF (PIVOT(J).NE.ZERO) GO TO 105
      DO 100 K=1,N
      IF (PIVOT(K).NE.ZERO) GO TO 100
      AVAL=A(J,K)
      TEMP=DABS(AVAL)
      IF (TEMP.LT.AMAX) GO TO 100
      IROW=J
      ICOLUM=K
      AMAX=TEMP
  100 CONTINUE
  105 CONTINUE
      INDEX(I)=4096*IROW+ICOLUM
      J=IROW
      AMAX=A(J,ICOLUM)
      DETERM=AMAX*DETERM
C     RETURN IF MATRIX IS SINGULAR (ZERO PIVOT) AFTER COLUMN INTERCHANGE
      IF(DETERM.EQ.ZERO) GO TO 600
      PIVOT(ICOLUM)=AMAX
C     INTERCHANGE ROWS TO PUT PIVOT ELEMENT ON DIAGONAL
      IF (IROW.EQ.ICOLUM) GO TO 260
      DETERM=-DETERM
      DO 200 K=1,N
      SWAP=A(J,K)
      A(J,K)=A(ICOLUM,K)
  200 A(ICOLUM,K)=SWAP
      IF (M.LE.0) GO TO 260
      DO 250 K=1,M
      SWAP=B(J,K)
      B(J,K)=B(ICOLUM,K)
  250 B(ICOLUM,K)=SWAP
C     DIVIDE PIVOT ROW BY PIVOT ELEMENT
  260 K=ICOLUM
      A(ICOLUM,K)=ONE
      DO 350 K=1,N
  350 A(ICOLUM,K)=A(ICOLUM,K)/AMAX
      IF (M.LE.0) GO TO 380
      DO 370 K=1,M
  370 B(ICOLUM,K)=B(ICOLUM,K)/AMAX
C     REDUCE NON-PIVOT ROWS
  380 DO 550 J=1,N
      IF (J.EQ.ICOLUM) GO TO 550
      T=A( J,ICOLUM)
      A( J,ICOLUM)=ZERO
      DO 450 K=1,N
  450 A( J,K)=A( J,K)-A(ICOLUM,K)*T
      IF (M.LE.0) GO TO 550
      DO 500 K=1,M
  500 B( J,K)=B( J,K)-B(ICOLUM,K)*T
  550 CONTINUE
C     INTERCHANGE COLUMNS AFTER ALL PIVOT OPERATIONS HAVE BEEN PERFORMED
  600 DO 710 I=1,N
      I1=N+1-I
      K=INDEX(I1)/4096
      ICOLUM=INDEX(I1)-4096*K
      IF (K.EQ.ICOLUM) GO TO 710
      DO 705 J=1,N
      SWAP=A(J,K)
      A(J,K)=A(J,ICOLUM)
  705 A(J,ICOLUM)=SWAP
  710 CONTINUE
      RETURN
      END
      SUBROUTINE FLIN(A,IDIM,IN,IM,DET)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     LINEAR SIMULTANEOUS EQUATION
C
C     A(IN*IN) * X(IN*IM) = B(IN*IM)
C
C     A & B SHOULD BE STORED ON A(IN*(IN+IM))
C     SOLUTION X WILL BE STORED ON B PART IN DIMENSION A.
C
      DIMENSION A(IDIM,1)
      DATA ZERO,ONE / 0.0D+00 , 1.0D+00 /
C
      N=IN
      NR=IM
      JMAX=N+NR
      SIGN=ONE
C M IS THE STAGE OF ELIMINATION
      DO 49 M=1,N
      TEMP=ZERO
      DO 41 I=M,N
      IF(M.GT.1)A(I,M)=A(I,M)-DOTT(A(I,1),IDIM,A(1,M),1,M-1)
      AVAL=A(I,M)
      IF(DABS(AVAL).LE.TEMP)GOTO 41
      TEMP=DABS(AVAL)
      IMAX=I
 41   CONTINUE
      IF(TEMP.LE.ZERO)GOTO 999
      IF(IMAX.EQ.M)GOTO 45
      SIGN=-SIGN
      DO 44 J=1,JMAX
      STOR=A(M,J)
      A(M,J)=A(IMAX,J)
      A(IMAX,J)=STOR
 44   CONTINUE
 45   CONTINUE
      JJ=M+1
      IF(JJ.GT.JMAX)GOTO 49
      IF(M.GT.1)GOTO 47
      DO 46 J=JJ,JMAX
      A(1,J)=A(1,J)/A(1,1)
 46   CONTINUE
      D=A(1,1)
      GOTO 49
 47   CONTINUE
      DO 48 J=JJ,JMAX
      A(M,J)=(A(M,J)-DOTT(A(M,1),IDIM,A(1,J),1,M-1))/A(M,M)
 48   CONTINUE
      D=D*A(M,M)
 49   CONTINUE
      IF(NR.EQ.0) RETURN
      DO 59 I=1,NR
      NPI=N+I
      DO 58 K=2,N
      J=N+1-K
      A(J,NPI)=A(J,NPI)-DOTT(A(J,J+1),IDIM,A(J+1,NPI),1,K-1)
 58   CONTINUE
 59   CONTINUE
C***  IF(DABS(D).GE.1.0D+36) D=1.0D+36
C***  IF(DABS(D).LE.1.0D-36) D=1.0D-36
      DET=D*SIGN
      RETURN
C ON ZERO PIVOT, SET DET=0.AND RETURN TO CALLING PROGRAM NOV 1972
 999  DET=ZERO
      RETURN
      END
CBHL  FUNCTION DOTT(A,NA,B,NB,N)
C     REAL*8 FUNCTION DOTT(A,NA,B,NB,N)
Cibm
      DOUBLE PRECISION FUNCTION DOTT(A,NA,B,NB,N)
Cibm
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(1),B(1)
      DATA ZERO / 0.0D+00 /
C
      IAPT=1
      IBPT=1
      D=ZERO
      DO 10 I=1,N
      D=D+A(IAPT)*B(IBPT)
      IAPT=IAPT+NA
      IBPT=IBPT+NB
 10   CONTINUE
      DOTT=D
      RETURN
      END
      SUBROUTINE FLINQ(A,IDIM,IN,IM,DET)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 DET,D,SIGN
C
C     LINEAR SIMULTANEOUS EQUATION
C
C     A(IN*IN) * X(IN*IM) = B(IN*IM)
C
C     A & B SHOULD BE STORED ON A(IN*(IN+IM))
C     SOLUTION X WILL BE STORED ON B PART IN DIMENSION A.
C
      DIMENSION A(IDIM,1)
      DATA ZERO,ONE / 0.0D+00 , 1.0D+00 /
C
      N=IN
      NR=IM
      JMAX=N+NR
      SIGN=1.0D+00
C M IS THE STAGE OF ELIMINATION
      DO 49 M=1,N
      TEMP=ZERO
      DO 41 I=M,N
      IF(M.GT.1)A(I,M)=A(I,M)-DOTT(A(I,1),IDIM,A(1,M),1,M-1)
      AVAL=A(I,M)
      IF(DABS(AVAL).LE.TEMP)GOTO 41
      TEMP=DABS(AVAL)
      IMAX=I
 41   CONTINUE
      IF(TEMP.LE.ZERO)GOTO 999
      IF(IMAX.EQ.M)GOTO 45
      SIGN=-SIGN
      DO 44 J=1,JMAX
      STOR=A(M,J)
      A(M,J)=A(IMAX,J)
      A(IMAX,J)=STOR
 44   CONTINUE
 45   CONTINUE
      JJ=M+1
      IF(JJ.GT.JMAX)GOTO 49
      IF(M.GT.1)GOTO 47
      DO 46 J=JJ,JMAX
      A(1,J)=A(1,J)/A(1,1)
 46   CONTINUE
      D=A(1,1)
      GOTO 49
 47   CONTINUE
      DO 48 J=JJ,JMAX
      A(M,J)=(A(M,J)-DOTT(A(M,1),IDIM,A(1,J),1,M-1))/A(M,M)
 48   CONTINUE
      D=D*A(M,M)
 49   CONTINUE
      IF(NR.EQ.0) RETURN
      DO 59 I=1,NR
      NPI=N+I
      DO 58 K=2,N
      J=N+1-K
      A(J,NPI)=A(J,NPI)-DOTT(A(J,J+1),IDIM,A(J+1,NPI),1,K-1)
 58   CONTINUE
 59   CONTINUE
      DET=D*SIGN
      RETURN
C ON ZERO PIVOT, SET DET=0.AND RETURN TO CALLING PROGRAM NOV 1972
 999  DET=0.0D+00
      RETURN
      END
      SUBROUTINE RT123
C             *****   VERSION FEBRUARY 13,1975   *****
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/ROOT/X,U(9),W(9),NROOTS
      EQUIVALENCE (U(1),RT1),(U(2),RT2),(U(3),RT3),(U(4),RT4),(U(5),RT5)
      EQUIVALENCE (W(1),WW1),(W(2),WW2),(W(3),WW3),(W(4),WW4),(W(5),WW5)
      DATA R12,PIE4/2.75255128608411D-01, 7.85398163397448D-01/
      DATA R22,W22/ 2.72474487139158D+00, 9.17517095361369D-02/
      DATA R13/     1.90163509193487D-01/
      DATA R23,W23/ 1.78449274854325D+00, 1.77231492083829D-01/
      DATA R33,W33/ 5.52534374226326D+00, 5.11156880411248D-03/
      IF(X.GT.5.0D+00) GO TO 50
      IF(X.GT.1.0D+00) GO TO 30
      IF(X.GT.3.0D-07) GO TO 20
C     X IS APPROXIMATELY ZERO.         NROOTS=1,2, OR 3
      IF(NROOTS-2) 11,12,13
   11 RT1= 0.5D+00 -X/5.0D+00
      WW1= 1.0D+00 -X/3.0D+00
      RETURN
   12 RT1=1.30693606237085D-01    -2.90430236082028D-02 *X
      RT2=2.86930639376291D+00    -6.37623643058102D-01 *X
      WW1=6.52145154862545D-01    -1.22713621927067D-01 *X
      WW2=3.47854845137453D-01    -2.10619711404725D-01 *X
      RETURN
   13 RT1=6.03769246832797D-02    -9.28875764357368D-03 *X
      RT2=7.76823355931043D-01    -1.19511285527878D-01 *X
      RT3=6.66279971938567D+00    -1.02504611068957D+00 *X
      WW1=4.67913934572691D-01    -5.64876917232519D-02 *X
      WW2=3.60761573048137D-01    -1.49077186455208D-01 *X
      WW3=1.71324492379169D-01    -1.27768455150979D-01 *X
      RETURN
C     X = 0.0 TO 1.0                   NROOTS=1,2, OR 3
   20 IF(NROOTS.EQ.3) GO TO 23
      F1=          ((((((((-8.36313918003957D-08*X+1.21222603512827D-06
     1)*X-1.15662609053481D-05 )*X+9.25197374512647D-05
     2)*X-6.40994113129432D-04 )*X+3.78787044215009D-03
     3)*X-1.85185172458485D-02 )*X+7.14285713298222D-02
     4)*X-1.99999999997023D-01 )*X+3.33333333333318D-01
      WW1=(X+X)*F1+EXP(-X)
      IF(NROOTS.EQ.2) GO TO 22
      RT1=F1/(WW1-F1)
      RETURN
   22 RT1=           (((((((-2.35234358048491D-09*X+2.49173650389842D-08
     1)*X-4.558315364581D-08)*X-2.447252174587D-06)*X+4.743292959463D-05
     2)*X-5.33184749432408D-04 )*X+4.44654947116579D-03
     3)*X-2.90430236084697D-02 )*X+1.30693606237085D-01
      RT2=           (((((((-2.47404902329170D-08*X+2.36809910635906D-07
     1)*X+1.835367736310D-06)*X-2.066168802076D-05)*X-1.345693393936D-04
     2)*X-5.88154362858038D-05 )*X+5.32735082098139D-02
     3)*X-6.37623643056745D-01 )*X+2.86930639376289D+00
      WW2=((F1-WW1)*RT1+F1)*(1.0D+00+RT2)/(RT2-RT1)
      WW1=WW1-WW2
      RETURN
   23 RT1=            ((((((-5.10186691538870D-10*X+2.40134415703450D-08
     1)*X-5.01081057744427D-07 )*X+7.58291285499256D-06
     2)*X-9.55085533670919D-05 )*X+1.02893039315878D-03
     3)*X-9.28875764374337D-03 )*X+6.03769246832810D-02
      RT2=            ((((((-1.29646524960555D-08*X+7.74602292865683D-08
     1)*X+1.56022811158727D-06 )*X-1.58051990661661D-05
     2)*X-3.30447806384059D-04 )*X+9.74266885190267D-03
     3)*X-1.19511285526388D-01 )*X+7.76823355931033D-01
      RT3=            ((((((-9.28536484109606D-09*X-3.02786290067014D-07
     1)*X-2.50734477064200D-06 )*X-7.32728109752881D-06
     2)*X+2.44217481700129D-04 )*X+4.94758452357327D-02
     3)*X-1.02504611065774D+00 )*X+6.66279971938553D+00
      F2=          ((((((((-7.60911486098850D-08*X+1.09552870123182D-06
     1)*X-1.03463270693454D-05 )*X+8.16324851790106D-05
     2)*X-5.55526624875562D-04 )*X+3.20512054753924D-03
     3)*X-1.51515139838540D-02 )*X+5.55555554649585D-02
     4)*X-1.42857142854412D-01 )*X+1.99999999999986D-01
  300 G=EXP(-X)
      F1=((X+X)*F2+G)/3.0D+00
      WW1=(X+X)*F1+G
  301 T1=RT1/(RT1+1.0D+00)
      T2=RT2/(RT2+1.0D+00)
      T3=RT3/(RT3+1.0D+00)
      A2=F2-T1*F1
      A1=F1-T1*WW1
      WW3=(A2-T2*A1)/((T3-T2)*(T3-T1))
      WW2=(T3*A1-A2)/((T3-T2)*(T2-T1))
      WW1=WW1-WW2-WW3
      RETURN
   30 IF(X.GT.3.0D+00) GO TO 40
C     X = 1.0 TO 3.0                   NROOTS=1,2, OR 3
      Y=X-2.0D+00
      IF(NROOTS.EQ.3) GO TO 33
      F1=        ((((((((((-1.61702782425558D-10*Y+1.96215250865776D-09
     1)*Y-2.14234468198419D-08 )*Y+2.17216556336318D-07
     2)*Y-1.98850171329371D-06 )*Y+1.62429321438911D-05
     3)*Y-1.16740298039895D-04 )*Y+7.24888732052332D-04
     4)*Y-3.79490003707156D-03 )*Y+1.61723488664661D-02
     5)*Y-5.29428148329736D-02 )*Y+1.15702180856167D-01
      WW1=(X+X)*F1+EXP(-X)
      IF(NROOTS.EQ.2) GO TO 32
      RT1=F1/(WW1-F1)
      RETURN
   32 RT1=         (((((((((-6.36859636616415D-12*Y+8.47417064776270D-11
     1)*Y-5.152207846962D-10)*Y-3.846389873308D-10)*Y+8.472253388380D-08
     2)*Y-1.85306035634293D-06 )*Y+2.47191693238413D-05
     3)*Y-2.49018321709815D-04 )*Y+2.19173220020161D-03
     4)*Y-1.63329339286794D-02 )*Y+8.68085688285261D-02
      RT2=         ((((((((( 1.45331350488343D-10*Y+2.07111465297976D-09
     1)*Y-1.878920917404D-08)*Y-1.725838516261D-07)*Y+2.247389642339D-06
     2)*Y+9.76783813082564D-06 )*Y-1.93160765581969D-04
     3)*Y-1.58064140671893D-03 )*Y+4.85928174507904D-02
     4)*Y-4.30761584997596D-01 )*Y+1.80400974537950D+00
      WW2=((F1-WW1)*RT1+F1)*(1.0D+00+RT2)/(RT2-RT1)
      WW1=WW1-WW2
      RETURN
   33 RT1=          (((((((( 1.44687969563318D-12*Y+4.85300143926755D-12
     1)*Y-6.55098264095516D-10 )*Y+1.56592951656828D-08
     2)*Y-2.60122498274734D-07 )*Y+3.86118485517386D-06
     3)*Y-5.13430986707889D-05 )*Y+6.03194524398109D-04
     4)*Y-6.11219349825090D-03 )*Y+4.52578254679079D-02
      RT2=           ((((((( 6.95964248788138D-10*Y-5.35281831445517D-09
     1)*Y-6.745205954533D-08)*Y+1.502366784525D-06)*Y+9.923326947376D-07
     2)*Y-3.89147469249594D-04 )*Y+7.51549330892401D-03
     3)*Y-8.48778120363400D-02 )*Y+5.73928229597613D-01
      RT3=          ((((((((-2.81496588401439D-10*Y+3.61058041895031D-09
     1)*Y+4.53631789436255D-08 )*Y-1.40971837780847D-07
     2)*Y-6.05865557561067D-06 )*Y-5.15964042227127D-05
     3)*Y+3.34761560498171D-05 )*Y+5.04871005319119D-02
     4)*Y-8.24708946991557D-01 )*Y+4.81234667357205D+00
      F2=        ((((((((((-1.48044231072140D-10*Y+1.78157031325097D-09
     1)*Y-1.92514145088973D-08 )*Y+1.92804632038796D-07
     2)*Y-1.73806555021045D-06 )*Y+1.39195169625425D-05
     3)*Y-9.74574633246452D-05 )*Y+5.83701488646511D-04
     4)*Y-2.89955494844975D-03 )*Y+1.13847001113810D-02
     5)*Y-3.23446977320647D-02 )*Y+5.29428148329709D-02
      GO TO 300
C     X = 3.0 TO 5.0                   NROOTS =1,2, OR 3
   40 Y=X-4.0D+00
      IF(NROOTS.EQ.3) GO TO 43
      F1=        ((((((((((-2.62453564772299D-11*Y+3.24031041623823D-10
     1)*Y-3.614965656163D-09)*Y+3.760256799971D-08)*Y-3.553558319675D-07
     2)*Y+3.022556449731D-06)*Y-2.290098979647D-05)*Y+1.526537461148D-04
     3)*Y-8.81947375894379D-04 )*Y+4.33207949514611D-03
     4)*Y-1.75257821619926D-02 )*Y+5.28406320615584D-02
      WW1=(X+X)*F1+EXP(-X)
      IF(NROOTS.EQ.2) GO TO 42
      RT1=F1/(WW1-F1)
      RETURN
   42 RT1=          ((((((((-4.11560117487296D-12*Y+7.10910223886747D-11
     1)*Y-1.73508862390291D-09 )*Y+5.93066856324744D-08
     2)*Y-9.76085576741771D-07 )*Y+1.08484384385679D-05
     3)*Y-1.12608004981982D-04 )*Y+1.16210907653515D-03
     4)*Y-9.89572595720351D-03 )*Y+6.12589701086408D-02
      RT2=         (((((((((-1.80555625241001D-10*Y+5.44072475994123D-10
     1)*Y+1.603498045240D-08)*Y-1.497986283037D-07)*Y-7.017002532106D-07
     2)*Y+1.85882653064034D-05 )*Y-2.04685420150802D-05
     3)*Y-2.49327728643089D-03 )*Y+3.56550690684281D-02
     4)*Y-2.60417417692375D-01 )*Y+1.12155283108289D+00
      WW2=((F1-WW1)*RT1+F1)*(1.0D+00+RT2)/(RT2-RT1)
      WW1=WW1-WW2
      RETURN
   43 RT1=           ((((((( 1.44265709189601D-11*Y-4.66622033006074D-10
     1)*Y+7.649155832025D-09)*Y-1.229940017368D-07)*Y+2.026002142457D-06
     2)*Y-2.87048671521677D-05 )*Y+3.70326938096287D-04
     3)*Y-4.21006346373634D-03 )*Y+3.50898470729044D-02
      RT2=          ((((((((-2.65526039155651D-11*Y+1.97549041402552D-10
     1)*Y+2.15971131403034D-09 )*Y-7.95045680685193D-08
     2)*Y+5.15021914287057D-07 )*Y+1.11788717230514D-05
     3)*Y-3.33739312603632D-04 )*Y+5.30601428208358D-03
     4)*Y-5.93483267268959D-02 )*Y+4.31180523260239D-01
      RT3=          ((((((((-3.92833750584041D-10*Y-4.16423229782280D-09
     1)*Y+4.42413039572867D-08 )*Y+6.40574545989551D-07
     2)*Y-3.05512456576552D-06 )*Y-1.05296443527943D-04
     3)*Y-6.14120969315617D-04 )*Y+4.89665802767005D-02
     4)*Y-6.24498381002855D-01 )*Y+3.36412312243724D+00
      F2=        ((((((((((-2.36788772599074D-11*Y+2.89147476459092D-10
     1)*Y-3.18111322308846D-09 )*Y+3.25336816562485D-08
     2)*Y-3.00873821471489D-07 )*Y+2.48749160874431D-06
     3)*Y-1.81353179793672D-05 )*Y+1.14504948737066D-04
     4)*Y-6.10614987696677D-04 )*Y+2.64584212770942D-03
     5)*Y-8.66415899015349D-03 )*Y+1.75257821619922D-02
      GO TO 300
   50 IF(X.GT.15.0D+00) GO TO 70
      G=EXP(-X)
      IF(X.GT.10.0D+00) GO TO 60
C     X = 5.0 TO 10.0                  NROOTS =1,2, OR 3
      WW1=    (((((( 4.6897511375022D-01/X-6.9955602298985D-01)/X
     1+5.3689283271887D-01)/X-3.2883030418398D-01)/X
     2+2.4645596956002D-01)/X-4.9984072848436D-01)/X
     3-3.1501078774085D-06)*G + SQRT(PIE4/X)
      F1=(WW1-G)/(X+X)
      IF(NROOTS-2) 51,52,53
   51 RT1=F1/(WW1-F1)
      RETURN
   52 Y=X-7.5D+00
      RT1=     (((((((((((((-1.43632730148572D-16*Y+2.38198922570405D-16
     1)*Y+1.358319618800D-14)*Y-7.064522786879D-14)*Y-7.719300212748D-13
     2)*Y+7.802544789997D-12)*Y+6.628721099436D-11)*Y-1.775564159743D-09
     3)*Y+1.713828823990D-08)*Y-1.497500187053D-07)*Y+2.283485114279D-06
     4)*Y-3.76953869614706D-05 )*Y+4.74791204651451D-04
     5)*Y-4.60448960876139D-03 )*Y+3.72458587837249D-02
      RT2=      (((((((((((( 2.48791622798900D-14*Y-1.36113510175724D-13
     1)*Y-2.224334349799D-12)*Y+4.190559455515D-11)*Y-2.222722579924D-10
     2)*Y-2.624183464275D-09)*Y+6.128153450169D-08)*Y-4.383376014528D-07
     3)*Y-2.49952200232910D-06 )*Y+1.03236647888320D-04
     4)*Y-1.44614664924989D-03 )*Y+1.35094294917224D-02
     5)*Y-9.53478510453887D-02 )*Y+5.44765245686790D-01
      WW2=((F1-WW1)*RT1+F1)*(1.0D+00+RT2)/(RT2-RT1)
      WW1=WW1-WW2
      RETURN
   53 F2=(F1+F1+F1-G)/(X+X)
      Y=X-7.5D+00
      RT1=       ((((((((((( 5.74429401360115D-16*Y+7.11884203790984D-16
     1)*Y-6.736701449826D-14)*Y-6.264613873998D-13)*Y+1.315418927040D-11
     2)*Y-4.23879635610964D-11 )*Y+1.39032379769474D-09
     3)*Y-4.65449552856856D-08 )*Y+7.34609900170759D-07
     4)*Y-1.08656008854077D-05 )*Y+1.77930381549953D-04
     5)*Y-2.39864911618015D-03 )*Y+2.39112249488821D-02
      RT2=       ((((((((((( 1.13464096209120D-14*Y+6.99375313934242D-15
     1)*Y-8.595618132088D-13)*Y-5.293620408757D-12)*Y-2.492175211635D-11
     2)*Y+2.73681574882729D-09 )*Y-1.06656985608482D-08
     3)*Y-4.40252529648056D-07 )*Y+9.68100917793911D-06
     4)*Y-1.68211091755327D-04 )*Y+2.69443611274173D-03
     5)*Y-3.23845035189063D-02 )*Y+2.75969447451882D-01
      RT3=      (((((((((((( 6.66339416996191D-15*Y+1.84955640200794D-13
     1)*Y-1.985141104444D-12)*Y-2.309293727603D-11)*Y+3.917984522103D-10
     2)*Y+1.663165279876D-09)*Y-6.205591993923D-08)*Y+8.769581622041D-09
     3)*Y+8.97224398620038D-06 )*Y-3.14232666170796D-05
     4)*Y-1.83917335649633D-03 )*Y+3.51246831672571D-02
     5)*Y-3.22335051270860D-01 )*Y+1.73582831755430D+00
      GO TO 301
C     X = 10.0 TO 15.0                 NROOTS=1,2, OR 3
   60 WW1=       (((-1.8784686463512D-01/X+2.2991849164985D-01)/X
     1-4.9893752514047D-01)/X-2.1916512131607D-05)*G + SQRT(PIE4/X)
      F1=(WW1-G)/(X+X)
      IF(NROOTS-2) 61,62,63
   61 RT1=F1/(WW1-F1)
      RETURN
   62 RT1=      ((((-1.01041157064226D-05*X+1.19483054115173D-03)*X
     1 -6.73760231824074D-02)*X+1.25705571069895D+00)*X
     2 +     (((-8.57609422987199D+03/X+5.91005939591842D+03)/X
     3 -1.70807677109425D+03)/X+2.64536689959503D+02)/X
     4 -2.38570496490846D+01)*G + R12/(X-R12)
      RT2=       ((( 3.39024225137123D-04*X-9.34976436343509D-02)*X
     1 -4.22216483306320D+00)*X +       (((-2.08457050986847D+03/X
     2 -1.04999071905664D+03)/X+3.39891508992661D+02)/X
     3 -1.56184800325063D+02)/X+8.00839033297501D+00)*G + R22/(X-R22)
      WW2=((F1-WW1)*RT1+F1)*(1.0D+00+RT2)/(RT2-RT1)
      WW1=WW1-WW2
      RETURN
   63 F2=(F1+F1+F1-G)/(X+X)
      Y=X-12.5D+00
      RT1=       ((((((((((( 4.42133001283090D-16*Y-2.77189767070441D-15
     1)*Y-4.084026087887D-14)*Y+5.379885121517D-13)*Y+1.882093066702D-12
     2)*Y-8.67286219861085D-11 )*Y+7.11372337079797D-10
     3)*Y-3.55578027040563D-09 )*Y+1.29454702851936D-07
     4)*Y-4.14222202791434D-06 )*Y+8.04427643593792D-05
     5)*Y-1.18587782909876D-03 )*Y+1.53435577063174D-02
      RT2=       ((((((((((( 6.85146742119357D-15*Y-1.08257654410279D-14
     1)*Y-8.579165965128D-13)*Y+6.642452485783D-12)*Y+4.798806828724D-11
     2)*Y-1.13413908163831D-09 )*Y+7.08558457182751D-09
     3)*Y-5.59678576054633D-08 )*Y+2.51020389884249D-06
     4)*Y-6.63678914608681D-05 )*Y+1.11888323089714D-03
     5)*Y-1.45361636398178D-02 )*Y+1.65077877454402D-01
      RT3=      (((((((((((( 3.20622388697743D-15*Y-2.73458804864628D-14
     1)*Y-3.157134329361D-13)*Y+8.654129268056D-12)*Y-5.625235879301D-11
     2)*Y-7.718080513708D-10)*Y+2.064664199164D-08)*Y-1.567725007761D-07
     3)*Y-1.57938204115055D-06 )*Y+6.27436306915967D-05
     4)*Y-1.01308723606946D-03 )*Y+1.13901881430697D-02
     5)*Y-1.01449652899450D-01 )*Y+7.77203937334739D-01
      GO TO 301
   70 IF(X.GT.33.0D+00) GO TO 90
C     X = 15.0 TO 33.0                 NROOTS=1,2, OR 3
      G=EXP(-X)
      WW1=        (( 1.9623264149430D-01/X-4.9695241464490D-01)/X
     1-6.0156581186481D-05)*G + SQRT(PIE4/X)
      F1=(WW1-G)/(X+X)
      IF(NROOTS-2) 71,72,73
   71 RT1=F1/(WW1-F1)
      RETURN
   72 RT1=      ((((-1.14906395546354D-06*X+1.76003409708332D-04)*X
     1 -1.71984023644904D-02)*X-1.37292644149838D-01)*X
     2 +       (-4.75742064274859D+01/X+9.21005186542857D+00)/X
     3 -2.31080873898939D-02)*G + R12/(X-R12)
      RT2=       ((( 3.64921633404158D-04*X-9.71850973831558D-02)*X
     1 -4.02886174850252D+00)*X +         (-1.35831002139173D+02/X
     2 -8.66891724287962D+01)/X+2.98011277766958D+00)*G + R22/(X-R22)
      WW2=((F1-WW1)*RT1+F1)*(1.0D+00+RT2)/(RT2-RT1)
      WW1=WW1-WW2
      RETURN
   73 F2=(F1+F1+F1-G)/(X+X)
      IF(X.GT.20.0D+00) GO TO 83
      RT1=    ((((((-2.43270989903742D-06*X+3.57901398988359D-04)*X
     1 -2.34112415981143D-02)*X+7.81425144913975D-01)*X
     2 -1.73209218219175D+01)*X+2.43517435690398D+02)*X
     3 +       (-1.97611541576986D+04/X+9.82441363463929D+03)/X
     4 -2.07970687843258D+03)*G + R13/(X-R13)
      RT2=     (((((-2.62627010965435D-04*X+3.49187925428138D-02)*X
     1 -3.09337618731880D+00)*X+1.07037141010778D+02)*X
     2 -2.36659637247087D+03)*X +        ((-2.91669113681020D+06/X
     3 +1.41129505262758D+06)/X-2.91532335433779D+05)/X
     4 +3.35202872835409D+04)*G + R23/(X-R23)
      RT3=     ((((( 9.31856404738601D-05*X-2.87029400759565D-02)*X
     1 -7.83503697918455D-01)*X-1.84338896480695D+01)*X
     2 +4.04996712650414D+02)*X +         (-1.89829509315154D+05/X
     3 +5.11498390849158D+04)/X-6.88145821789955D+03)*G + R33/(X-R33)
      GO TO 301
   83 RT1=      ((((-4.97561537069643D-04*X-5.00929599665316D-02)*X
     1 +1.31099142238996D+00)*X-1.88336409225481D+01)*X
     2 -6.60344754467191D+02 /X+1.64931462413877D+02)*G + R13/(X-R13)
      RT2=      ((((-4.48218898474906D-03*X-5.17373211334924D-01)*X
     1 +1.13691058739678D+01)*X-1.65426392885291D+02)*X
     2 -6.30909125686731D+03 /X+1.52231757709236D+03)*G + R23/(X-R23)
      RT3=      ((((-1.38368602394293D-02*X-1.77293428863008D+00)*X
     1 +1.73639054044562D+01)*X-3.57615122086961D+02)*X
     2 -1.45734701095912D+04 /X+2.69831813951849D+03)*G + R33/(X-R33)
      GO TO 301
C     X = 33.0 TO INFINITY             NROOTS=1,2, OR 3
   90 WW1=SQRT(PIE4/X)
      IF(NROOTS-2) 91,92,93
   91 RT1=0.5D+00/(X-0.5D+00)
      RETURN
   92 IF(X.GT.40.0D+00) GO TO 102
      G=EXP(-X)
      RT1=(-8.78947307498880D-01*X+1.09243702330261D+01)*G + R12/(X-R12)
      RT2=(-9.28903924275977D+00*X+8.10642367843811D+01)*G + R22/(X-R22)
      WW2=( 4.46857389308400D+00*X-7.79250653461045D+01)*G + W22*WW1
      WW1=WW1-WW2
      RETURN
   93 IF(X.GT.47.0D+00) GO TO 103
      G=EXP(-X)
      RT1=        ((-7.39058467995275D+00*X+3.21318352526305D+02)*X
     1 -3.99433696473658D+03)*G + R13/(X-R13)
      RT2=        ((-7.38726243906513D+01*X+3.13569966333873D+03)*X
     1 -3.86862867311321D+04)*G + R23/(X-R23)
      RT3=        ((-2.63750565461336D+02*X+1.04412168692352D+04)*X
     1 -1.28094577915394D+05)*G + R33/(X-R33)
      WW3=       ((( 1.52258947224714D-01*X-8.30661900042651D+00)*X
     1 +1.92977367967984D+02)*X-1.67787926005344D+03)*G + W33*WW1
      WW2=        (( 6.15072615497811D+01*X-2.91980647450269D+03)*X
     1 +3.80794303087338D+04)*G + W23*WW1
      WW1=WW1-WW2-WW3
      RETURN
  102 RT1=R12/(X-R12)
      RT2=R22/(X-R22)
      WW2=W22*WW1
      WW1=WW1-WW2
      RETURN
  103 RT1=R13/(X-R13)
      RT2=R23/(X-R23)
      RT3=R33/(X-R33)
      WW2=W23*WW1
      WW3=W33*WW1
      WW1=WW1-WW2-WW3
      RETURN
      END
      SUBROUTINE ROOT4
C          *****   VERSION FEBRUARY 16,1975   *****
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/ROOT/X,U(9),W(9),NROOTS
      EQUIVALENCE (U(1),RT1),(U(2),RT2),(U(3),RT3),(U(4),RT4),(U(5),RT5)
      EQUIVALENCE (W(1),WW1),(W(2),WW2),(W(3),WW3),(W(4),WW4),(W(5),WW5)
      DATA R14,PIE4/1.45303521503316D-01, 7.85398163397448D-01/
      DATA R24,W24/ 1.33909728812636D+00, 2.34479815323517D-01/
      DATA R34,W34/ 3.92696350135829D+00, 1.92704402415764D-02/
      DATA R44,W44/ 8.58863568901199D+00, 2.25229076750736D-04/
      IF(X.GT.15.0D+00) GO TO 470
      IF(X.GT.5.0D+00) GO TO 450
      IF(X.GT.1.0D+00) GO TO 430
      IF(X.GT.3.0D-07) GO TO 420
C     X IS APPROXIMATELY ZERO.                   NROOTS = 4
      RT1=3.48198973061471D-02    -4.09645850660395D-03 *X
      RT2=3.81567185080042D-01    -4.48902570656719D-02 *X
      RT3=1.73730726945891D+00    -2.04389090547327D-01 *X
      RT4=1.18463056481549D+01    -1.39368301742312D+00 *X
      WW1=3.62683783378362D-01    -3.13844305713928D-02 *X
      WW2=3.13706645877886D-01    -8.98046242557724D-02 *X
      WW3=2.22381034453372D-01    -1.29314370958973D-01 *X
      WW4=1.01228536290376D-01    -8.28299075414321D-02 *X
      RETURN
C     X=0.0 TO 1.0                               NROOTS = 4
  420 RT1=            ((((((-1.95309614628539D-10*X+5.19765728707592D-09
     1)*X-1.01756452250573D-07 )*X+1.72365935872131D-06
     2)*X-2.61203523522184D-05 )*X+3.52921308769880D-04
     3)*X-4.09645850658433D-03 )*X+3.48198973061469D-02
      RT2=             (((((-1.89554881382342D-08*X+3.07583114342365D-07
     1)*X+1.270981734393D-06)*X-1.417298563884D-04)*X+3.226979163176D-03
     2)*X-4.48902570678178D-02 )*X+3.81567185080039D-01
      RT3=            (((((( 1.77280535300416D-09*X+3.36524958870615D-08
     1)*X-2.58341529013893D-07 )*X-1.13644895662320D-05
     2)*X-7.91549618884063D-05 )*X+1.03825827346828D-02
     3)*X-2.04389090525137D-01 )*X+1.73730726945889D+00
      RT4=             (((((-5.61188882415248D-08*X-2.49480733072460D-07
     1)*X+3.428685057114D-06)*X+1.679007454539D-04)*X+4.722855585715D-02
     2)*X-1.39368301737828D+00 )*X+1.18463056481543D+01
      WW1=            ((((((-1.14649303201279D-08*X+1.88015570196787D-07
     1)*X-2.33305875372323D-06 )*X+2.68880044371597D-05
     2)*X-2.94268428977387D-04 )*X+3.06548909776613D-03
     3)*X-3.13844305680096D-02 )*X+3.62683783378335D-01
      WW2=          ((((((((-4.11720483772634D-09*X+6.54963481852134D-08
     1)*X-7.20045285129626D-07 )*X+6.93779646721723D-06
     2)*X-6.05367572016373D-05 )*X+4.74241566251899D-04
     3)*X-3.26956188125316D-03 )*X+1.91883866626681D-02
     4)*X-8.98046242565811D-02 )*X+3.13706645877886D-01
      WW3=          ((((((((-3.41688436990215D-08*X+5.07238960340773D-07
     1)*X-5.01675628408220D-06 )*X+4.20363420922845D-05
     2)*X-3.08040221166823D-04 )*X+1.94431864731239D-03
     3)*X-1.02477820460278D-02 )*X+4.28670143840073D-02
     4)*X-1.29314370962569D-01 )*X+2.22381034453369D-01
      WW4=         ((((((((( 4.99660550769508D-09*X-7.94585963310120D-08
     1)*X+8.359072409485D-07)*X-7.422369210610D-06)*X+5.763374308160D-05
     2)*X-3.86645606718233D-04 )*X+2.18417516259781D-03
     3)*X-9.99791027771119D-03 )*X+3.48791097377370D-02
     4)*X-8.28299075413889D-02 )*X+1.01228536290376D-01
      RETURN
C     X= 1.0 TO 5.0                              NROOTS = 4
  430 Y=X-3.0D+00
      RT1=         (((((((((-1.48570633747284D-15*Y-1.33273068108777D-13
     1)*Y+4.068543696670D-12)*Y-9.163164161821D-11)*Y+2.046819017845D-09
     2)*Y-4.03076426299031D-08 )*Y+7.29407420660149D-07
     3)*Y-1.23118059980833D-05 )*Y+1.88796581246938D-04
     4)*Y-2.53262912046853D-03 )*Y+2.51198234505021D-02
      RT2=         ((((((((( 1.35830583483312D-13*Y-2.29772605964836D-12
     1)*Y-3.821500128045D-12)*Y+6.844424214735D-10)*Y-1.048063352259D-08
     2)*Y+1.50083186233363D-08 )*Y+3.48848942324454D-06
     3)*Y-1.08694174399193D-04 )*Y+2.08048885251999D-03
     4)*Y-2.91205805373793D-02 )*Y+2.72276489515713D-01
      RT3=         ((((((((( 5.02799392850289D-13*Y+1.07461812944084D-11
     1)*Y-1.482277886411D-10)*Y-2.153585661215D-09)*Y+3.654087802817D-08
     2)*Y+5.15929575830120D-07 )*Y-9.52388379435709D-06
     3)*Y-2.16552440036426D-04 )*Y+9.03551469568320D-03
     4)*Y-1.45505469175613D-01 )*Y+1.21449092319186D+00
      RT4=         (((((((((-1.08510370291979D-12*Y+6.41492397277798D-11
     1)*Y+7.542387436125D-10)*Y-2.213111836647D-09)*Y-1.448228963549D-07
     2)*Y-1.95670833237101D-06 )*Y-1.07481314670844D-05
     3)*Y+1.49335941252765D-04 )*Y+4.87791531990593D-02
     4)*Y-1.10559909038653D+00 )*Y+8.09502028611780D+00
      WW1=        ((((((((((-4.65801912689961D-14*Y+7.58669507106800D-13
     1)*Y-1.186387548048D-11)*Y+1.862334710665D-10)*Y-2.799399389539D-09
     2)*Y+4.148972684255D-08)*Y-5.933568079600D-07)*Y+8.168349266115D-06
     3)*Y-1.08989176177409D-04 )*Y+1.41357961729531D-03
     4)*Y-1.87588361833659D-02 )*Y+2.89898651436026D-01
      WW2=      ((((((((((((-1.46345073267549D-14*Y+2.25644205432182D-13
     1)*Y-3.116258693847D-12)*Y+4.321908756610D-11)*Y-5.673270062669D-10
     2)*Y+7.006295962960D-09)*Y-8.120186517000D-08)*Y+8.775294645770D-07
     3)*Y-8.77829235749024D-06 )*Y+8.04372147732379D-05
     4)*Y-6.64149238804153D-04 )*Y+4.81181506827225D-03
     5)*Y-2.88982669486183D-02 )*Y+1.56247249979288D-01
      WW3=     ((((((((((((( 9.06812118895365D-15*Y-1.40541322766087D-13
     1)*Y+1.919270015269D-12)*Y-2.605135739010D-11)*Y+3.299685839012D-10
     2)*Y-3.86354139348735D-09 )*Y+4.16265847927498D-08
     3)*Y-4.09462835471470D-07 )*Y+3.64018881086111D-06
     4)*Y-2.88665153269386D-05 )*Y+2.00515819789028D-04
     5)*Y-1.18791896897934D-03 )*Y+5.75223633388589D-03
     6)*Y-2.09400418772687D-02 )*Y+4.85368861938873D-02
      WW4=    ((((((((((((((-9.74835552342257D-16*Y+1.57857099317175D-14
     1)*Y-2.249993780112D-13)*Y+3.173422008953D-12)*Y-4.161159459680D-11
     2)*Y+5.021343560166D-10)*Y-5.545047534808D-09)*Y+5.554146993491D-08
     3)*Y-4.99048696190133D-07 )*Y+3.96650392371311D-06
     4)*Y-2.73816413291214D-05 )*Y+1.60106988333186D-04
     5)*Y-7.64560567879592D-04 )*Y+2.81330044426892D-03
     6)*Y-7.16227030134947D-03 )*Y+9.66077262223353D-03
      RETURN
  450 IF(X.GT.10.0D+00) GO TO 460
C     X=5.0 TO 10.0                              NROOTS = 4
      Y=X-7.5D+00
      RT1=         ((((((((( 4.64217329776215D-15*Y-6.27892383644164D-15
     1)*Y+3.462236347446D-13)*Y-2.927229355350D-11)*Y+5.090355371676D-10
     2)*Y-9.97272656345253D-09 )*Y+2.37835295639281D-07
     3)*Y-4.60301761310921D-06 )*Y+8.42824204233222D-05
     4)*Y-1.37983082233081D-03 )*Y+1.66630865869375D-02
      RT2=         ((((((((( 2.93981127919047D-14*Y+8.47635639065744D-13
     1)*Y-1.446314544774D-11)*Y-6.149155555753D-12)*Y+8.484275604612D-10
     2)*Y-6.10898827887652D-08 )*Y+2.39156093611106D-06
     3)*Y-5.35837089462592D-05 )*Y+1.00967602595557D-03
     4)*Y-1.57769317127372D-02 )*Y+1.74853819464285D-01
      RT3=        (((((((((( 2.93523563363000D-14*Y-6.40041776667020D-14
     1)*Y-2.695740446312D-12)*Y+1.027082960169D-10)*Y-5.822038656780D-10
     2)*Y-3.159991002539D-08)*Y+4.327249251331D-07)*Y+4.856768455119D-06
     3)*Y-2.54617989427762D-04 )*Y+5.54843378106589D-03
     4)*Y-7.95013029486684D-02 )*Y+7.20206142703162D-01
      RT4=       (((((((((((-1.62212382394553D-14*Y+7.68943641360593D-13
     1)*Y+5.764015756615D-12)*Y-1.380635298784D-10)*Y-1.476849808675D-09
     2)*Y+1.84347052385605D-08 )*Y+3.34382940759405D-07
     3)*Y-1.39428366421645D-06 )*Y-7.50249313713996D-05
     4)*Y-6.26495899187507D-04 )*Y+4.69716410901162D-02
     5)*Y-6.66871297428209D-01 )*Y+4.11207530217806D+00
      WW1=        ((((((((((-1.65995045235997D-15*Y+6.91838935879598D-14
     1)*Y-9.131223418888D-13)*Y+1.403341829454D-11)*Y-3.672235069444D-10
     2)*Y+6.366962546990D-09)*Y-1.039220021671D-07)*Y+1.959098751715D-06
     3)*Y-3.33474893152939D-05 )*Y+5.72164211151013D-04
     4)*Y-1.05583210553392D-02 )*Y+2.26696066029591D-01
      WW2=      ((((((((((((-3.57248951192047D-16*Y+6.25708409149331D-15
     1)*Y-9.657033089714D-14)*Y+1.507864898748D-12)*Y-2.332522256110D-11
     2)*Y+3.428545616603D-10)*Y-4.698730937661D-09)*Y+6.219977635130D-08
     3)*Y-7.83008889613661D-07 )*Y+9.08621687041567D-06
     4)*Y-9.86368311253873D-05 )*Y+9.69632496710088D-04
     5)*Y-8.14594214284187D-03 )*Y+8.50218447733457D-02
      WW3=     ((((((((((((( 1.64742458534277D-16*Y-2.68512265928410D-15
     1)*Y+3.788890667676D-14)*Y-5.508918529823D-13)*Y+7.555896810069D-12
     2)*Y-9.69039768312637D-11 )*Y+1.16034263529672D-09
     3)*Y-1.28771698573873D-08 )*Y+1.31949431805798D-07
     4)*Y-1.23673915616005D-06 )*Y+1.04189803544936D-05
     5)*Y-7.79566003744742D-05 )*Y+5.03162624754434D-04
     6)*Y-2.55138844587555D-03 )*Y+1.13250730954014D-02
      WW4=    ((((((((((((((-1.55714130075679D-17*Y+2.57193722698891D-16
     1)*Y-3.626606654097D-15)*Y+5.234734676175D-14)*Y-7.067105402134D-13
     2)*Y+8.793512664890D-12)*Y-1.006088923498D-10)*Y+1.050565098393D-09
     3)*Y-9.91517881772662D-09 )*Y+8.35835975882941D-08
     4)*Y-6.19785782240693D-07 )*Y+3.95841149373135D-06
     5)*Y-2.11366761402403D-05 )*Y+9.00474771229507D-05
     6)*Y-2.78777909813289D-04 )*Y+5.26543779837487D-04
      RETURN
C     X=10.0 TO 15.0                             NROOTS = 4
  460 Y=X-12.5D+00
      RT1=       ((((((((((( 4.94869622744119D-17*Y+8.03568805739160D-16
     1)*Y-5.599125915431D-15)*Y-1.378685560217D-13)*Y+7.006511663249D-13
     2)*Y+1.30391406991118D-11 )*Y+8.06987313467541D-11
     3)*Y-5.20644072732933D-09 )*Y+7.72794187755457D-08
     4)*Y-1.61512612564194D-06 )*Y+4.15083811185831D-05
     5)*Y-7.87855975560199D-04 )*Y+1.14189319050009D-02
      RT2=       ((((((((((( 4.89224285522336D-16*Y+1.06390248099712D-14
     1)*Y-5.446260182933D-14)*Y-1.613630106295D-12)*Y+3.910179118937D-12
     2)*Y+1.90712434258806D-10 )*Y+8.78470199094761D-10
     3)*Y-5.97332993206797D-08 )*Y+9.25750831481589D-07
     4)*Y-2.02362185197088D-05 )*Y+4.92341968336776D-04
     5)*Y-8.68438439874703D-03 )*Y+1.15825965127958D-01
      RT3=        (((((((((( 6.12419396208408D-14*Y+1.12328861406073D-13
     1)*Y-9.051094103059D-12)*Y-4.781797525341D-11)*Y+1.660828868694D-09
     2)*Y+4.499058798868D-10)*Y-2.519549641933D-07)*Y+4.977444040180D-06
     3)*Y-1.25858350034589D-04 )*Y+2.70279176970044D-03
     4)*Y-3.99327850801083D-02 )*Y+4.33467200855434D-01
      RT4=       ((((((((((( 4.63414725924048D-14*Y-4.72757262693062D-14
     1)*Y-1.001926833832D-11)*Y+6.074107718414D-11)*Y+1.576976911942D-09
     2)*Y-2.01186401974027D-08 )*Y-1.84530195217118D-07
     3)*Y+5.02333087806827D-06 )*Y+9.66961790843006D-06
     4)*Y-1.58522208889528D-03 )*Y+2.80539673938339D-02
     5)*Y-2.78953904330072D-01 )*Y+1.82835655238235D+00
      WW4=     ((((((((((((( 2.90401781000996D-18*Y-4.63389683098251D-17
     1)*Y+6.274018198326D-16)*Y-8.936002188168D-15)*Y+1.194719074934D-13
     2)*Y-1.45501321259466D-12 )*Y+1.64090830181013D-11
     3)*Y-1.71987745310181D-10 )*Y+1.63738403295718D-09
     4)*Y-1.39237504892842D-08 )*Y+1.06527318142151D-07
     5)*Y-7.27634957230524D-07 )*Y+4.12159381310339D-06
     6)*Y-1.74648169719173D-05 )*Y+8.50290130067818D-05
      WW3=      ((((((((((((-4.19569145459480D-17*Y+5.94344180261644D-16
     1)*Y-1.148797566469D-14)*Y+1.881303962576D-13)*Y-2.413554618391D-12
     2)*Y+3.372127423047D-11)*Y-4.933988617784D-10)*Y+6.116545396281D-09
     3)*Y-6.69965691739299D-08 )*Y+7.52380085447161D-07
     4)*Y-8.08708393262321D-06 )*Y+6.88603417296672D-05
     5)*Y-4.67067112993427D-04 )*Y+5.42313365864597D-03
      WW2=        ((((((((((-6.22272689880615D-15*Y+1.04126809657554D-13
     1)*Y-6.842418230913D-13)*Y+1.576841731919D-11)*Y-4.203948834175D-10
     2)*Y+6.287255934781D-09)*Y-8.307159819228D-08)*Y+1.356478091922D-06
     3)*Y-2.08065576105639D-05 )*Y+2.52396730332340D-04
     4)*Y-2.94484050194539D-03 )*Y+6.01396183129168D-02
      WW1=       (((-1.8784686463512D-01/X+2.2991849164985D-01)/X
     1-4.9893752514047D-01)/X-2.1916512131607D-05)*EXP(-X)
     2 +SQRT(PIE4/X)-WW4-WW3-WW2
      RETURN
  470 WW1=SQRT(PIE4/X)
      IF(X.GT.35.0D+00) GO TO 490
      IF(X.GT.20.0D+00) GO TO 480
C     X=15.0 TO 20.0                             NROOTS = 4
      Y=X-17.5D+00
      RT1=       ((((((((((( 4.36701759531398D-17*Y-1.12860600219889D-16
     1)*Y-6.149849164164D-15)*Y+5.820231579541D-14)*Y+4.396602872143D-13
     2)*Y-1.24330365320172D-11 )*Y+6.71083474044549D-11
     3)*Y+2.43865205376067D-10 )*Y+1.67559587099969D-08
     4)*Y-9.32738632357572D-07 )*Y+2.39030487004977D-05
     5)*Y-4.68648206591515D-04 )*Y+8.34977776583956D-03
      RT2=       ((((((((((( 4.98913142288158D-16*Y-2.60732537093612D-16
     1)*Y-7.775156445127D-14)*Y+5.766105220086D-13)*Y+6.432696729600D-12
     2)*Y-1.39571683725792D-10 )*Y+5.95451479522191D-10
     3)*Y+2.42471442836205D-09 )*Y+2.47485710143120D-07
     4)*Y-1.14710398652091D-05 )*Y+2.71252453754519D-04
     5)*Y-4.96812745851408D-03 )*Y+8.26020602026780D-02
      RT3=       ((((((((((( 1.91498302509009D-15*Y+1.48840394311115D-14
     1)*Y-4.316925145767D-13)*Y+1.186495793471D-12)*Y+4.615806713055D-11
     2)*Y-5.54336148667141D-10 )*Y+3.48789978951367D-10
     3)*Y-2.79188977451042D-09 )*Y+2.09563208958551D-06
     4)*Y-6.76512715080324D-05 )*Y+1.32129867629062D-03
     5)*Y-2.05062147771513D-02 )*Y+2.88068671894324D-01
      RT4=       (((((((((((-5.43697691672942D-15*Y-1.12483395714468D-13
     1)*Y+2.826607936174D-12)*Y-1.266734493280D-11)*Y-4.258722866437D-10
     2)*Y+9.45486578503261D-09 )*Y-5.86635622821309D-08
     3)*Y-1.28835028104639D-06 )*Y+4.41413815691885D-05
     4)*Y-7.61738385590776D-04 )*Y+9.66090902985550D-03
     5)*Y-1.01410568057649D-01 )*Y+9.54714798156712D-01
      WW4=      ((((((((((((-7.56882223582704D-19*Y+7.53541779268175D-18
     1)*Y-1.157318032236D-16)*Y+2.411195002314D-15)*Y-3.601794386996D-14
     2)*Y+4.082150659615D-13)*Y-4.289542980767D-12)*Y+5.086829642731D-11
     3)*Y-6.35435561050807D-10 )*Y+6.82309323251123D-09
     4)*Y-5.63374555753167D-08 )*Y+3.57005361100431D-07
     5)*Y-2.40050045173721D-06 )*Y+4.94171300536397D-05
      WW3=       (((((((((((-5.54451040921657D-17*Y+2.68748367250999D-16
     1)*Y+1.349020069254D-14)*Y-2.507452792892D-13)*Y+1.944339743818D-12
     2)*Y-1.29816917658823D-11 )*Y+3.49977768819641D-10
     3)*Y-8.67270669346398D-09 )*Y+1.31381116840118D-07
     4)*Y-1.36790720600822D-06 )*Y+1.19210697673160D-05
     5)*Y-1.42181943986587D-04 )*Y+4.12615396191829D-03
      WW2=       (((((((((((-1.86506057729700D-16*Y+1.16661114435809D-15
     1)*Y+2.563712856363D-14)*Y-4.498350984631D-13)*Y+1.765194089338D-12
     2)*Y+9.04483676345625D-12 )*Y+4.98930345609785D-10
     3)*Y-2.11964170928181D-08 )*Y+3.98295476005614D-07
     4)*Y-5.49390160829409D-06 )*Y+7.74065155353262D-05
     5)*Y-1.48201933009105D-03 )*Y+4.97836392625268D-02
      WW1=        (( 1.9623264149430D-01/X-4.9695241464490D-01)/X
     1-6.0156581186481D-05)*EXP(-X)+WW1-WW2-WW3-WW4
      RETURN
C     X=20.0 TO 35.0                             NROOTS = 4
  480 G=EXP(-X)
      RT1=    ((((((-4.45711399441838D-05*X+1.27267770241379D-03)*X
     1 -2.36954961381262D-01)*X+1.54330657903756D+01)*X
     2 -5.22799159267808D+02)*X+1.05951216669313D+04)*X
     3 +       (-2.51177235556236D+06/X+8.72975373557709D+05)/X
     4 -1.29194382386499D+05)*G + R14/(X-R14)
      RT2=     (((((-7.85617372254488D-02*X+6.35653573484868D+00)*X
     1 -3.38296938763990D+02)*X+1.25120495802096D+04)*X
     2 -3.16847570511637D+05)*X +        ((-1.02427466127427D+09/X
     3 +3.70104713293016D+08)/X-5.87119005093822D+07)/X
     4 +5.38614211391604D+06)*G + R24/(X-R24)
      RT3=     (((((-2.37900485051067D-01*X+1.84122184400896D+01)*X
     1 -1.00200731304146D+03)*X+3.75151841595736D+04)*X
     2 -9.50626663390130D+05)*X +        ((-2.88139014651985D+09/X
     3 +1.06625915044526D+09)/X-1.72465289687396D+08)/X
     4 +1.60419390230055D+07)*G + R34/(X-R34)
      RT4=    ((((((-6.00691586407385D-04*X-3.64479545338439D-01)*X
     1 +1.57496131755179D+01)*X-6.54944248734901D+02)*X
     2 +1.70830039597097D+04)*X-2.90517939780207D+05)*X
     3 +       (+3.49059698304732D+07/X-1.64944522586065D+07)/X
     4 +2.96817940164703D+06)*G + R44/(X-R44)
      IF(X.LE.25.0D+00)
     1WW4=   ((((((( 2.33766206773151D-07*X-3.81542906607063D-05)*X
     1 +3.51416601267000D-03)*X-1.66538571864728D-01)*X
     2 +4.80006136831847D+00)*X-8.73165934223603D+01)*X
     3 +9.77683627474638D+02)*X +           1.66000945117640D+04/X
     4 -6.14479071209961D+03)*G + W44*WW1
      IF(X.GT.25.0D+00)
     1WW4=    (((((( 5.74245945342286D-06*X-7.58735928102351D-05)*X
     1 +2.35072857922892D-04)*X-3.78812134013125D-03)*X
     2 +3.09871652785805D-01)*X-7.11108633061306D+00)*X
     3 +5.55297573149528D+01)*G + W44*WW1
      WW3=    (((((( 2.36392855180768D-04*X-9.16785337967013D-03)*X
     1 +4.62186525041313D-01)*X-1.96943786006540D+01)*X
     2 +4.99169195295559D+02)*X-6.21419845845090D+03)*X
     3 +      ((+5.21445053212414D+07/X-1.34113464389309D+07)/X
     4 +1.13673298305631D+06)/X-2.81501182042707D+03)*G + W34*WW1
      WW2=    (((((( 7.29841848989391D-04*X-3.53899555749875D-02)*X
     1 +2.07797425718513D+00)*X-1.00464709786287D+02)*X
     2 +3.15206108877819D+03)*X-6.27054715090012D+04)*X
     3 +       (+1.54721246264919D+07/X-5.26074391316381D+06)/X
     4 +7.67135400969617D+05)*G + W24*WW1
      WW1=        (( 1.9623264149430D-01/X-4.9695241464490D-01)/X
     1-6.0156581186481D-05)*G + WW1-WW2-WW3-WW4
      RETURN
  490 IF(X.GT.53.0D+00) GO TO 495
C     X=35.0 TO 53.0                             NROOTS = 4
      G=EXP(-X)*(X*X)**2
      RT4=        ((-2.19135070169653D-03*X-1.19108256987623D-01)*X
     1 -7.50238795695573D-01)*G + R44/(X-R44)
      RT3=        ((-9.65842534508637D-04*X-4.49822013469279D-02)*X
     1 +6.08784033347757D-01)*G + R34/(X-R34)
      RT2=        ((-3.62569791162153D-04*X-9.09231717268466D-03)*X
     1 +1.84336760556262D-01)*G + R24/(X-R24)
      RT1=        ((-4.07557525914600D-05*X-6.88846864931685D-04)*X
     1 +1.74725309199384D-02)*G + R14/(X-R14)
      WW4=        (( 5.76631982000990D-06*X-7.89187283804890D-05)*X
     1 +3.28297971853126D-04)*G + W44*WW1
      WW3=        (( 2.08294969857230D-04*X-3.77489954837361D-03)*X
     1 +2.09857151617436D-02)*G + W34*WW1
      WW2=        (( 6.16374517326469D-04*X-1.26711744680092D-02)*X
     1 +8.14504890732155D-02)*G + W24*WW1
      WW1=WW1-WW2-WW3-WW4
      RETURN
C     X=47.0 TO INFINITY                         NROOTS = 4
  495 RT1=R14/(X-R14)
      RT2=R24/(X-R24)
      RT3=R34/(X-R34)
      RT4=R44/(X-R44)
      WW4=W44*WW1
      WW3=W34*WW1
      WW2=W24*WW1
      WW1=WW1-WW2-WW3-WW4
      RETURN
      END
      SUBROUTINE ROOT5
C          *****   VERSION  FEBRUARY 27,1975   *****
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/ROOT/X,U(9),W(9),NROOTS
      EQUIVALENCE (U(1),RT1),(U(2),RT2),(U(3),RT3),(U(4),RT4),(U(5),RT5)
      EQUIVALENCE (W(1),WW1),(W(2),WW2),(W(3),WW3),(W(4),WW4),(W(5),WW5)
      DATA R15,PIE4/1.17581320211778D-01, 7.85398163397448D-01/
      DATA R25,W25/ 1.07456201243690D+00, 2.70967405960535D-01/
      DATA R35,W35/ 3.08593744371754D+00, 3.82231610015404D-02/
      DATA R45,W45/ 6.41472973366203D+00, 1.51614186862443D-03/
      DATA R55,W55/ 1.18071894899717D+01, 8.62130526143657D-06/
      IF(X.GT.15.0D+00) GO TO 570
      IF(X.GT. 5.0D+00) GO TO 550
      IF(X.GT. 1.0D+00) GO TO 530
      IF(X.GT. 3.0D-07) GO TO 520
C     X IS APPROXIMATELY ZERO.                   NROOTS = 5
      RT1=2.26659266316985D-02    -2.15865967920897D-03 *X
      RT2=2.31271692140903D-01    -2.20258754389745D-02 *X
      RT3=8.57346024118836D-01    -8.16520023025515D-02 *X
      RT4=2.97353038120346D+00    -2.83193369647137D-01 *X
      RT5=1.84151859759051D+01    -1.75382723579439D+00 *X
      WW1=2.95524224714752D-01    -1.96867576909777D-02 *X
      WW2=2.69266719309995D-01    -5.61737590184721D-02 *X
      WW3=2.19086362515981D-01    -9.71152726793658D-02 *X
      WW4=1.49451349150580D-01    -1.02979262193565D-01 *X
      WW5=6.66713443086877D-02    -5.73782817488315D-02 *X
      RETURN
C     X=0.0 TO 1.0                               NROOTS = 5
  520 RT1=            ((((((-4.46679165328413D-11*X+1.21879111988031D-09
     1)*X-2.62975022612104D-08 )*X+5.15106194905897D-07
     2)*X-9.27933625824749D-06 )*X+1.51794097682482D-04
     3)*X-2.15865967920301D-03 )*X+2.26659266316985D-02
      RT2=            (((((( 1.93117331714174D-10*X-4.57267589660699D-09
     1)*X+2.48339908218932D-08 )*X+1.50716729438474D-06
     2)*X-6.07268757707381D-05 )*X+1.37506939145643D-03
     3)*X-2.20258754419939D-02 )*X+2.31271692140905D-01
      RT3=             ((((( 4.84989776180094D-09*X+1.31538893944284D-07
     1)*X-2.766753852879D-06)*X-7.651163510626D-05)*X+4.033058545972D-03
     2)*X-8.16520022916145D-02 )*X+8.57346024118779D-01
      RT4=              ((((-2.48581772214623D-07*X-4.34482635782585D-06
     1)*X-7.46018257987630D-07 )*X+1.01210776517279D-02
     2)*X-2.83193369640005D-01 )*X+2.97353038120345D+00
      RT5=             (((((-8.92432153868554D-09*X+1.77288899268988D-08
     1)*X+3.040754680666D-06)*X+1.058229325071D-04)*X+4.596379534985D-02
     2)*X-1.75382723579114D+00 )*X+1.84151859759049D+01
      WW1=            ((((((-2.03822632771791D-09*X+3.89110229133810D-08
     1)*X-5.84914787904823D-07 )*X+8.30316168666696D-06
     2)*X-1.13218402310546D-04 )*X+1.49128888586790D-03
     3)*X-1.96867576904816D-02 )*X+2.95524224714749D-01
      WW2=           ((((((( 8.62848118397570D-09*X-1.38975551148989D-07
     1)*X+1.602894068228D-06)*X-1.646364300836D-05)*X+1.538445806778D-04
     2)*X-1.28848868034502D-03 )*X+9.38866933338584D-03
     3)*X-5.61737590178812D-02 )*X+2.69266719309991D-01
      WW3=          ((((((((-9.41953204205665D-09*X+1.47452251067755D-07
     1)*X-1.57456991199322D-06 )*X+1.45098401798393D-05
     2)*X-1.18858834181513D-04 )*X+8.53697675984210D-04
     3)*X-5.22877807397165D-03 )*X+2.60854524809786D-02
     4)*X-9.71152726809059D-02 )*X+2.19086362515979D-01
      WW4=          ((((((((-3.84961617022042D-08*X+5.66595396544470D-07
     1)*X-5.52351805403748D-06 )*X+4.53160377546073D-05
     2)*X-3.22542784865557D-04 )*X+1.95682017370967D-03
     3)*X-9.77232537679229D-03 )*X+3.79455945268632D-02
     4)*X-1.02979262192227D-01 )*X+1.49451349150573D-01
      WW5=         ((((((((( 4.09594812521430D-09*X-6.47097874264417D-08
     1)*X+6.743541482689D-07)*X-5.917993920224D-06)*X+4.531969237381D-05
     2)*X-2.99102856679638D-04 )*X+1.65695765202643D-03
     3)*X-7.40671222520653D-03 )*X+2.50889946832192D-02
     4)*X-5.73782817487958D-02 )*X+6.66713443086877D-02
      RETURN
C     X=1.0 TO 5.0                               NROOTS = 5
  530 Y=X-3.0D+00
      RT1=          ((((((((-2.58163897135138D-14*Y+8.14127461488273D-13
     1)*Y-2.11414838976129D-11 )*Y+5.09822003260014D-10
     2)*Y-1.16002134438663D-08 )*Y+2.46810694414540D-07
     3)*Y-4.92556826124502D-06 )*Y+9.02580687971053D-05
     4)*Y-1.45190025120726D-03 )*Y+1.73416786387475D-02
      RT2=         ((((((((( 1.04525287289788D-14*Y+5.44611782010773D-14
     1)*Y-4.831059411392D-12)*Y+1.136643908832D-10)*Y-1.104373076913D-09
     2)*Y-2.35346740649916D-08 )*Y+1.43772622028764D-06
     3)*Y-4.23405023015273D-05 )*Y+9.12034574793379D-04
     4)*Y-1.52479441718739D-02 )*Y+1.76055265928744D-01
      RT3=         (((((((((-6.89693150857911D-14*Y+5.92064260918861D-13
     1)*Y+1.847170956043D-11)*Y-3.390752744265D-10)*Y-2.995532064116D-09
     2)*Y+1.57456141058535D-07 )*Y-3.95859409711346D-07
     3)*Y-9.58924580919747D-05 )*Y+3.23551502557785D-03
     4)*Y-5.97587007636479D-02 )*Y+6.46432853383057D-01
      RT4=          ((((((((-3.61293809667763D-12*Y-2.70803518291085D-11
     1)*Y+8.83758848468769D-10 )*Y+1.59166632851267D-08
     2)*Y-1.32581997983422D-07 )*Y-7.60223407443995D-06
     3)*Y-7.41019244900952D-05 )*Y+9.81432631743423D-03
     4)*Y-2.23055570487771D-01 )*Y+2.21460798080643D+00
      RT5=         ((((((((( 7.12332088345321D-13*Y+3.16578501501894D-12
     1)*Y-8.776668218053D-11)*Y-2.342817613343D-09)*Y-3.496962018025D-08
     2)*Y-3.03172870136802D-07 )*Y+1.50511293969805D-06
     3)*Y+1.37704919387696D-04 )*Y+4.70723869619745D-02
     4)*Y-1.47486623003693D+00 )*Y+1.35704792175847D+01
      WW1=         ((((((((( 1.04348658616398D-13*Y-1.94147461891055D-12
     1)*Y+3.485512360993D-11)*Y-6.277497362235D-10)*Y+1.100758247388D-08
     2)*Y-1.88329804969573D-07 )*Y+3.12338120839468D-06
     3)*Y-5.04404167403568D-05 )*Y+8.00338056610995D-04
     4)*Y-1.30892406559521D-02 )*Y+2.47383140241103D-01
      WW2=       ((((((((((( 3.23496149760478D-14*Y-5.24314473469311D-13
     1)*Y+7.743219385056D-12)*Y-1.146022750992D-10)*Y+1.615238462197D-09
     2)*Y-2.15479017572233D-08 )*Y+2.70933462557631D-07
     3)*Y-3.18750295288531D-06 )*Y+3.47425221210099D-05
     4)*Y-3.45558237388223D-04 )*Y+3.05779768191621D-03
     5)*Y-2.29118251223003D-02 )*Y+1.59834227924213D-01
      WW3=      ((((((((((((-3.42790561802876D-14*Y+5.26475736681542D-13
     1)*Y-7.184330797139D-12)*Y+9.763932908544D-11)*Y-1.244014559219D-09
     2)*Y+1.472744068942D-08)*Y-1.611749975234D-07)*Y+1.616487851917D-06
     3)*Y-1.46852359124154D-05 )*Y+1.18900349101069D-04
     4)*Y-8.37562373221756D-04 )*Y+4.93752683045845D-03
     5)*Y-2.25514728915673D-02 )*Y+6.95211812453929D-02
      WW4=     ((((((((((((( 1.04072340345039D-14*Y-1.60808044529211D-13
     1)*Y+2.183534866798D-12)*Y-2.939403008391D-11)*Y+3.679254029085D-10
     2)*Y-4.23775673047899D-09 )*Y+4.46559231067006D-08
     3)*Y-4.26488836563267D-07 )*Y+3.64721335274973D-06
     4)*Y-2.74868382777722D-05 )*Y+1.78586118867488D-04
     5)*Y-9.68428981886534D-04 )*Y+4.16002324339929D-03
     6)*Y-1.28290192663141D-02 )*Y+2.22353727685016D-02
      WW5=    ((((((((((((((-8.16770412525963D-16*Y+1.31376515047977D-14
     1)*Y-1.856950818865D-13)*Y+2.596836515749D-12)*Y-3.372639523006D-11
     2)*Y+4.025371849467D-10)*Y-4.389453269417D-09)*Y+4.332753856271D-08
     3)*Y-3.82673275931962D-07 )*Y+2.98006900751543D-06
     4)*Y-2.00718990300052D-05 )*Y+1.13876001386361D-04
     5)*Y-5.23627942443563D-04 )*Y+1.83524565118203D-03
     6)*Y-4.37785737450783D-03 )*Y+5.36963805223095D-03
      RETURN
  550 IF(X.GT.10.0D+00) GO TO 560
C     X=5.0 TO 10.0                              NROOTS = 5
      Y=X-7.5D+00
      RT1=          ((((((((-1.13825201010775D-14*Y+1.89737681670375D-13
     1)*Y-4.81561201185876D-12 )*Y+1.56666512163407D-10
     2)*Y-3.73782213255083D-09 )*Y+9.15858355075147D-08
     3)*Y-2.13775073585629D-06 )*Y+4.56547356365536D-05
     4)*Y-8.68003909323740D-04 )*Y+1.22703754069176D-02
      RT2=         (((((((((-3.67160504428358D-15*Y+1.27876280158297D-14
     1)*Y-1.296476623788D-12)*Y+1.477175434354D-11)*Y+5.464102147892D-10
     2)*Y-2.42538340602723D-08 )*Y+8.20460740637617D-07
     3)*Y-2.20379304598661D-05 )*Y+4.90295372978785D-04
     4)*Y-9.14294111576119D-03 )*Y+1.22590403403690D-01
      RT3=         ((((((((( 1.39017367502123D-14*Y-6.96391385426890D-13
     1)*Y+1.176946020731D-12)*Y+1.725627235645D-10)*Y-3.686383856300D-09
     2)*Y+2.87495324207095D-08 )*Y+1.71307311000282D-06
     3)*Y-7.94273603184629D-05 )*Y+2.00938064965897D-03
     4)*Y-3.63329491677178D-02 )*Y+4.34393683888443D-01
      RT4=        ((((((((((-1.27815158195209D-14*Y+1.99910415869821D-14
     1)*Y+3.753542914426D-12)*Y-2.708018219579D-11)*Y-1.190574776587D-09
     2)*Y+1.106696436509D-08)*Y+3.954955671326D-07)*Y-4.398596059588D-06
     3)*Y-2.01087998907735D-04 )*Y+7.89092425542937D-03
     4)*Y-1.42056749162695D-01 )*Y+1.39964149420683D+00
      RT5=        ((((((((((-1.19442341030461D-13*Y-2.34074833275956D-12
     1)*Y+6.861649627426D-12)*Y+6.082671496226D-10)*Y+5.381160105420D-09
     2)*Y-6.253297138700D-08)*Y-2.135966835050D-06)*Y-2.373394341886D-05
     3)*Y+2.88711171412814D-06 )*Y+4.85221195290753D-02
     4)*Y-1.04346091985269D+00 )*Y+7.89901551676692D+00
      WW1=         ((((((((( 7.95526040108997D-15*Y-2.48593096128045D-13
     1)*Y+4.761246208720D-12)*Y-9.535763686605D-11)*Y+2.225273630974D-09
     2)*Y-4.49796778054865D-08 )*Y+9.17812870287386D-07
     3)*Y-1.86764236490502D-05 )*Y+3.76807779068053D-04
     4)*Y-8.10456360143408D-03 )*Y+2.01097936411496D-01
      WW2=       ((((((((((( 1.25678686624734D-15*Y-2.34266248891173D-14
     1)*Y+3.973252415832D-13)*Y-6.830539401049D-12)*Y+1.140771033372D-10
     2)*Y-1.82546185762009D-09 )*Y+2.77209637550134D-08
     3)*Y-4.01726946190383D-07 )*Y+5.48227244014763D-06
     4)*Y-6.95676245982121D-05 )*Y+8.05193921815776D-04
     5)*Y-8.15528438784469D-03 )*Y+9.71769901268114D-02
      WW3=      ((((((((((((-8.20929494859896D-16*Y+1.37356038393016D-14
     1)*Y-2.022863065220D-13)*Y+3.058055403795D-12)*Y-4.387890955243D-11
     2)*Y+5.923946274445D-10)*Y-7.503659964159D-09)*Y+8.851599803902D-08
     3)*Y-9.65561998415038D-07 )*Y+9.60884622778092D-06
     4)*Y-8.56551787594404D-05 )*Y+6.66057194311179D-04
     5)*Y-4.17753183902198D-03 )*Y+2.25443826852447D-02
      WW4=    ((((((((((((((-1.08764612488790D-17*Y+1.85299909689937D-16
     1)*Y-2.730195628655D-15)*Y+4.127368817265D-14)*Y-5.881379088074D-13
     2)*Y+7.805245193391D-12)*Y-9.632707991704D-11)*Y+1.099047050624D-09
     3)*Y-1.15042731790748D-08 )*Y+1.09415155268932D-07
     4)*Y-9.33687124875935D-07 )*Y+7.02338477986218D-06
     5)*Y-4.53759748787756D-05 )*Y+2.41722511389146D-04
     6)*Y-9.75935943447037D-04 )*Y+2.57520532789644D-03
      WW5=   ((((((((((((((( 7.28996979748849D-19*Y-1.26518146195173D-17
     1)*Y+1.886145834486D-16)*Y-2.876728287383D-15)*Y+4.114588668138D-14
     2)*Y-5.44436631413933D-13 )*Y+6.64976446790959D-12
     3)*Y-7.44560069974940D-11 )*Y+7.57553198166848D-10
     4)*Y-6.92956101109829D-09 )*Y+5.62222859033624D-08
     5)*Y-3.97500114084351D-07 )*Y+2.39039126138140D-06
     6)*Y-1.18023950002105D-05 )*Y+4.52254031046244D-05
     7)*Y-1.21113782150370D-04 )*Y+1.75013126731224D-04
      RETURN
C     X=10.0 TO 15.0                             NROOTS = 5
  560 Y=X-12.5D+00
      RT1=        ((((((((((-4.16387977337393D-17*Y+7.20872997373860D-16
     1)*Y+1.395993802064D-14)*Y+3.660484641252D-14)*Y-4.154857548139D-12
     2)*Y+2.301379846544D-11)*Y-1.033307012866D-09)*Y+3.997777641049D-08
     3)*Y-9.35118186333939D-07 )*Y+2.38589932752937D-05
     4)*Y-5.35185183652937D-04 )*Y+8.85218988709735D-03
      RT2=        ((((((((((-4.56279214732217D-16*Y+6.24941647247927D-15
     1)*Y+1.737896339191D-13)*Y+8.964205979517D-14)*Y-3.538906780633D-11
     2)*Y+9.561341254948D-11)*Y-9.772831891310D-09)*Y+4.240340194620D-07
     3)*Y-1.02384302866534D-05 )*Y+2.57987709704822D-04
     4)*Y-5.54735977651677D-03 )*Y+8.68245143991948D-02
      RT3=        ((((((((((-2.52879337929239D-15*Y+2.13925810087833D-14
     1)*Y+7.884307667104D-13)*Y-9.023398159510D-13)*Y-5.814101544957D-11
     2)*Y-1.333480437968D-09)*Y-2.217064940373D-08)*Y+1.643290788086D-06
     3)*Y-4.39602147345028D-05 )*Y+1.08648982748911D-03
     4)*Y-2.13014521653498D-02 )*Y+2.94150684465425D-01
      RT4=        ((((((((((-6.42391438038888D-15*Y+5.37848223438815D-15
     1)*Y+8.960828117859D-13)*Y+5.214153461337D-11)*Y-1.106601744067D-10
     2)*Y-2.007890743962D-08)*Y+1.543764346501D-07)*Y+4.520749076914D-06
     3)*Y-1.88893338587047D-04 )*Y+4.73264487389288D-03
     4)*Y-7.91197893350253D-02 )*Y+8.60057928514554D-01
      RT5=       (((((((((((-2.24366166957225D-14*Y+4.87224967526081D-14
     1)*Y+5.587369053655D-12)*Y-3.045253104617D-12)*Y-1.223983883080D-09
     2)*Y-2.05603889396319D-09 )*Y+2.58604071603561D-07
     3)*Y+1.34240904266268D-06 )*Y-5.72877569731162D-05
     4)*Y-9.56275105032191D-04 )*Y+4.23367010370921D-02
     5)*Y-5.76800927133412D-01 )*Y+3.87328263873381D+00
      WW1=         ((((((((( 8.98007931950169D-15*Y+7.25673623859497D-14
     1)*Y+5.851494250405D-14)*Y-4.234204823846D-11)*Y+3.911507312679D-10
     2)*Y-9.65094802088511D-09 )*Y+3.42197444235714D-07
     3)*Y-7.51821178144509D-06 )*Y+1.94218051498662D-04
     4)*Y-5.38533819142287D-03 )*Y+1.68122596736809D-01
      WW2=        ((((((((((-1.05490525395105D-15*Y+1.96855386549388D-14
     1)*Y-5.500330153548D-13)*Y+1.003849567976D-11)*Y-1.720997242621D-10
     2)*Y+3.533277061402D-09)*Y-6.389171736029D-08)*Y+1.046236652393D-06
     3)*Y-1.73148206795827D-05 )*Y+2.57820531617185D-04
     4)*Y-3.46188265338350D-03 )*Y+7.03302497508176D-02
      WW3=       ((((((((((( 3.60020423754545D-16*Y-6.24245825017148D-15
     1)*Y+9.945311467434D-14)*Y-1.749051512721D-12)*Y+2.768503957853D-11
     2)*Y-4.08688551136506D-10 )*Y+6.04189063303610D-09
     3)*Y-8.23540111024147D-08 )*Y+1.01503783870262D-06
     4)*Y-1.20490761741576D-05 )*Y+1.26928442448148D-04
     5)*Y-1.05539461930597D-03 )*Y+1.15543698537013D-02
      WW4=     ((((((((((((( 2.51163533058925D-18*Y-4.31723745510697D-17
     1)*Y+6.557620865832D-16)*Y-1.016528519495D-14)*Y+1.491302084832D-13
     2)*Y-2.06638666222265D-12 )*Y+2.67958697789258D-11
     3)*Y-3.23322654638336D-10 )*Y+3.63722952167779D-09
     4)*Y-3.75484943783021D-08 )*Y+3.49164261987184D-07
     5)*Y-2.92658670674908D-06 )*Y+2.12937256719543D-05
     6)*Y-1.19434130620929D-04 )*Y+6.45524336158384D-04
      WW5=    ((((((((((((((-1.29043630202811D-19*Y+2.16234952241296D-18
     1)*Y-3.107631557965D-17)*Y+4.570804313173D-16)*Y-6.301348858104D-15
     2)*Y+8.031304476153D-14)*Y-9.446196472547D-13)*Y+1.018245804339D-11
     3)*Y-9.96995451348129D-11 )*Y+8.77489010276305D-10
     4)*Y-6.84655877575364D-09 )*Y+4.64460857084983D-08
     5)*Y-2.66924538268397D-07 )*Y+1.24621276265907D-06
     6)*Y-4.30868944351523D-06 )*Y+9.94307982432868D-06
      RETURN
  570 IF(X.GT.25.0D+00) GO TO 590
      IF(X.GT.20.0D+00) GO TO 580
C     X=15.0 TO 20.0                             NROOTS = 5
      Y=X-17.5D+00
      RT1=        (((((((((( 1.91875764545740D-16*Y+7.8357401095707 D-16
     1)*Y-3.260875931644D-14)*Y-1.186752035569D-13)*Y+4.275180095653D-12
     2)*Y+3.357056136731D-11)*Y-1.123776903884D-09)*Y+1.231203269887D-08
     3)*Y-3.99851421361031D-07 )*Y+1.45418822817771D-05
     4)*Y-3.49912254976317D-04 )*Y+6.67768703938812D-03
      RT2=        (((((((((( 2.02778478673555D-15*Y+1.01640716785099D-14
     1)*Y-3.385363492036D-13)*Y-1.615655871159D-12)*Y+4.527419140333D-11
     2)*Y+3.853670706486D-10)*Y-1.184607130107D-08)*Y+1.347873288827D-07
     3)*Y-4.47788241748377D-06 )*Y+1.54942754358273D-04
     4)*Y-3.55524254280266D-03 )*Y+6.44912219301603D-02
      RT3=        (((((((((( 7.79850771456444D-15*Y+6.00464406395001D-14
     1)*Y-1.249779730869D-12)*Y-1.020720636353D-11)*Y+1.814709816693D-10
     2)*Y+1.766397336977D-09)*Y-4.603559449010D-08)*Y+5.863956443581D-07
     3)*Y-2.03797212506691D-05 )*Y+6.31405161185185D-04
     4)*Y-1.30102750145071D-02 )*Y+2.10244289044705D-01
      RT4=       (((((((((((-2.92397030777912D-15*Y+1.94152129078465D-14
     1)*Y+4.859447665850D-13)*Y-3.217227223463D-12)*Y-7.484522135512D-11
     2)*Y+7.19101516047753D-10 )*Y+6.88409355245582D-09
     3)*Y-1.44374545515769D-07 )*Y+2.74941013315834D-06
     4)*Y-1.02790452049013D-04 )*Y+2.59924221372643D-03
     5)*Y-4.35712368303551D-02 )*Y+5.62170709585029D-01
      RT5=       ((((((((((( 1.17976126840060D-14*Y+1.24156229350669D-13
     1)*Y-3.892741622280D-12)*Y-7.755793199043D-12)*Y+9.492190032313D-10
     2)*Y-4.98680128123353D-09 )*Y-1.81502268782664D-07
     3)*Y+2.69463269394888D-06 )*Y+2.50032154421640D-05
     4)*Y-1.33684303917681D-03 )*Y+2.29121951862538D-02
     5)*Y-2.45653725061323D-01 )*Y+1.89999883453047D+00
      WW1=        (((((((((( 1.74841995087592D-15*Y-6.95671892641256D-16
     1)*Y-3.000659497257D-13)*Y+2.021279817961D-13)*Y+3.853596935400D-11
     2)*Y+1.461418533652D-10)*Y-1.014517563435D-08)*Y+1.132736008979D-07
     3)*Y-2.86605475073259D-06 )*Y+1.21958354908768D-04
     4)*Y-3.86293751153466D-03 )*Y+1.45298342081522D-01
      WW2=        ((((((((((-1.11199320525573D-15*Y+1.85007587796671D-15
     1)*Y+1.220613939709D-13)*Y+1.275068098526D-12)*Y-5.341838883262D-11
     2)*Y+6.161037256669D-10)*Y-1.009147879750D-08)*Y+2.907862965346D-07
     3)*Y-6.12300038720919D-06 )*Y+1.00104454489518D-04
     4)*Y-1.80677298502757D-03 )*Y+5.78009914536630D-02
      WW3=        ((((((((((-9.49816486853687D-16*Y+6.67922080354234D-15
     1)*Y+2.606163540537D-15)*Y+1.983799950150D-12)*Y-5.400548574357D-11
     2)*Y+6.638043374114D-10)*Y-8.799518866802D-09)*Y+1.791418482685D-07
     3)*Y-2.96075397351101D-06 )*Y+3.38028206156144D-05
     4)*Y-3.58426847857878D-04 )*Y+8.39213709428516D-03
      WW4=       ((((((((((( 1.33829971060180D-17*Y-3.44841877844140D-16
     1)*Y+4.745009557656D-15)*Y-6.033814209875D-14)*Y+1.049256040808D-12
     2)*Y-1.70859789556117D-11 )*Y+2.15219425727959D-10
     3)*Y-2.52746574206884D-09 )*Y+3.27761714422960D-08
     4)*Y-3.90387662925193D-07 )*Y+3.46340204593870D-06
     5)*Y-2.43236345136782D-05 )*Y+3.54846978585226D-04
      WW5=     ((((((((((((( 2.69412277020887D-20*Y-4.24837886165685D-19
     1)*Y+6.030500065438D-18)*Y-9.069722758289D-17)*Y+1.246599177672D-15
     2)*Y-1.56872999797549D-14 )*Y+1.87305099552692D-13
     3)*Y-2.09498886675861D-12 )*Y+2.11630022068394D-11
     4)*Y-1.92566242323525D-10 )*Y+1.62012436344069D-09
     5)*Y-1.23621614171556D-08 )*Y+7.72165684563049D-08
     6)*Y-3.59858901591047D-07 )*Y+2.43682618601000D-06
      RETURN
C     X=20.0 TO 25.0                             NROOTS = 5
  580 Y=X-22.5D+00
      RT1=         (((((((((-1.13927848238726D-15*Y+7.39404133595713D-15
     1)*Y+1.445982921243D-13)*Y-2.676703245252D-12)*Y+5.823521627177D-12
     2)*Y+2.17264723874381D-10 )*Y+3.56242145897468D-09
     3)*Y-3.03763737404491D-07 )*Y+9.46859114120901D-06
     4)*Y-2.30896753853196D-04 )*Y+5.24663913001114D-03
      RT2=        (((((((((( 2.89872355524581D-16*Y-1.22296292045864D-14
     1)*Y+6.184065097200D-14)*Y+1.649846591230D-12)*Y-2.729713905266D-11
     2)*Y+3.709913790650D-11)*Y+2.216486288382D-09)*Y+4.616160236414D-08
     3)*Y-3.32380270861364D-06 )*Y+9.84635072633776D-05
     4)*Y-2.30092118015697D-03 )*Y+5.00845183695073D-02
      RT3=        (((((((((( 1.97068646590923D-15*Y-4.89419270626800D-14
     1)*Y+1.136466605916D-13)*Y+7.546203883874D-12)*Y-9.635646767455D-11
     2)*Y-8.295965491209D-11)*Y+7.534109114453D-09)*Y+2.699970652707D-07
     3)*Y-1.42982334217081D-05 )*Y+3.78290946669264D-04
     4)*Y-8.03133015084373D-03 )*Y+1.58689469640791D-01
      RT4=        (((((((((( 1.33642069941389D-14*Y-1.55850612605745D-13
     1)*Y-7.522712577474D-13)*Y+3.209520801187D-11)*Y-2.075594313618D-10
     2)*Y-2.070575894402D-09)*Y+7.323046997451D-09)*Y+1.851491550417D-06
     3)*Y-6.37524802411383D-05 )*Y+1.36795464918785D-03
     4)*Y-2.42051126993146D-02 )*Y+3.97847167557815D-01
      RT5=        ((((((((((-6.07053986130526D-14*Y+1.04447493138843D-12
     1)*Y-4.286617818951D-13)*Y-2.632066100073D-10)*Y+4.804518986559D-09
     2)*Y-1.835675889421D-08)*Y-1.068175391334D-06)*Y+3.292234974141D-05
     3)*Y-5.94805357558251D-04 )*Y+8.29382168612791D-03
     4)*Y-9.93122509049447D-02 )*Y+1.09857804755042D+00
      WW1=         (((((((((-9.10338640266542D-15*Y+1.00438927627833D-13
     1)*Y+7.817349237071D-13)*Y-2.547619474232D-11)*Y+1.479321506529D-10
     2)*Y+1.52314028857627D-09 )*Y+9.20072040917242D-09
     3)*Y-2.19427111221848D-06 )*Y+8.65797782880311D-05
     4)*Y-2.82718629312875D-03 )*Y+1.28718310443295D-01
      WW2=         ((((((((( 5.52380927618760D-15*Y-6.43424400204124D-14
     1)*Y-2.358734508092D-13)*Y+8.261326648131D-12)*Y+9.229645304956D-11
     2)*Y-5.68108973828949D-09 )*Y+1.22477891136278D-07
     3)*Y-2.11919643127927D-06 )*Y+4.23605032368922D-05
     4)*Y-1.14423444576221D-03 )*Y+5.06607252890186D-02
      WW3=         ((((((((( 3.99457454087556D-15*Y-5.11826702824182D-14
     1)*Y-4.157593182747D-14)*Y+4.214670817758D-12)*Y+6.705582751532D-11
     2)*Y-3.36086411698418D-09 )*Y+6.07453633298986D-08
     3)*Y-7.40736211041247D-07 )*Y+8.84176371665149D-06
     4)*Y-1.72559275066834D-04 )*Y+7.16639814253567D-03
      WW4=       (((((((((((-2.14649508112234D-18*Y-2.45525846412281D-18
     1)*Y+6.126212599772D-16)*Y-8.526651626939D-15)*Y+4.826636065733D-14
     2)*Y-3.39554163649740D-13 )*Y+1.67070784862985D-11
     3)*Y-4.42671979311163D-10 )*Y+6.77368055908400D-09
     4)*Y-7.03520999708859D-08 )*Y+6.04993294708874D-07
     5)*Y-7.80555094280483D-06 )*Y+2.85954806605017D-04
      WW5=      ((((((((((((-5.63938733073804D-21*Y+6.92182516324628D-20
     1)*Y-1.586937691507D-18)*Y+3.357639744582D-17)*Y-4.810285046442D-16
     2)*Y+5.386312669975D-15)*Y-6.117895297439D-14)*Y+8.441808227634D-13
     3)*Y-1.18527596836592D-11 )*Y+1.36296870441445D-10
     4)*Y-1.17842611094141D-09 )*Y+7.80430641995926D-09
     5)*Y-5.97767417400540D-08 )*Y+1.65186146094969D-06
      RETURN
  590 WW1=SQRT(PIE4/X)
      IF(X.GT.40.0D+00) GO TO 595
C     X=25.0 TO 40.0                             NROOTS = 5
      G=EXP(-X)
      RT1=  ((((((((-1.73363958895356D-06*X+1.19921331441483D-04)*X
     1 -1.59437614121125D-02)*X+1.13467897349442D+00)*X
     2 -4.47216460864586D+01)*X+1.06251216612604D+03)*X
     3 -1.52073917378512D+04)*X+1.20662887111273D+05)*X
     4 -4.07186366852475D+05)*G + R15/(X-R15)
      RT2=  ((((((((-1.60102542621710D-05*X+1.10331262112395D-03)*X
     1 -1.50043662589017D-01)*X+1.05563640866077D+01)*X
     2 -4.10468817024806D+02)*X+9.62604416506819D+03)*X
     3 -1.35888069838270D+05)*X+1.06107577038340D+06)*X
     4 -3.51190792816119D+06)*G + R25/(X-R25)
      RT3=  ((((((((-4.48880032128422D-05*X+2.69025112122177D-03)*X
     1 -4.01048115525954D-01)*X+2.78360021977405D+01)*X
     2 -1.04891729356965D+03)*X+2.36985942687423D+04)*X
     3 -3.19504627257548D+05)*X+2.34879693563358D+06)*X
     4 -7.16341568174085D+06)*G + R35/(X-R35)
      RT4=  ((((((((-6.38526371092582D-05*X-2.29263585792626D-03)*X
     1 -7.65735935499627D-02)*X+9.12692349152792D+00)*X
     2 -2.32077034386717D+02)*X+2.81839578728845D+02)*X
     3 +9.59529683876419D+04)*X-1.77638956809518D+06)*X
     4 +1.02489759645410D+07)*G + R45/(X-R45)
      RT5=  ((((((((-3.59049364231569D-05*X-2.25963977930044D-02)*X
     1 +1.12594870794668D+00)*X-4.56752462103909D+01)*X
     2 +1.05804526830637D+03)*X-1.16003199605875D+04)*X
     3 -4.07297627297272D+04)*X+2.22215528319857D+06)*X
     4 -1.61196455032613D+07)*G + R55/(X-R55)
      WW5= (((((((((-4.61100906133970D-10*X+1.43069932644286D-07)*X
     1 -1.63960915431080D-05)*X+1.15791154612838D-03)*X
     2 -5.30573476742071D-02)*X+1.61156533367153D+00)*X
     3 -3.23248143316007D+01)*X+4.12007318109157D+02)*X
     4 -3.02260070158372D+03)*X+9.71575094154768D+03)*G + W55*WW1
      WW4= (((((((((-2.40799435809950D-08*X+8.12621667601546D-06)*X
     1 -9.04491430884113D-04)*X+6.37686375770059D-02)*X
     2 -2.96135703135647D+00)*X+9.15142356996330D+01)*X
     3 -1.86971865249111D+03)*X+2.42945528916947D+04)*X
     4 -1.81852473229081D+05)*X+5.96854758661427D+05)*G + W45*WW1
      WW3=  (((((((( 1.83574464457207D-05*X-1.54837969489927D-03)*X
     1 +1.18520453711586D-01)*X-6.69649981309161D+00)*X
     2 +2.44789386487321D+02)*X-5.68832664556359D+03)*X
     3 +8.14507604229357D+04)*X-6.55181056671474D+05)*X
     4 +2.26410896607237D+06)*G + W35*WW1
      WW2=  (((((((( 2.77778345870650D-05*X-2.22835017655890D-03)*X
     1 +1.61077633475573D-01)*X-8.96743743396132D+00)*X
     2 +3.28062687293374D+02)*X-7.65722701219557D+03)*X
     3 +1.10255055017664D+05)*X-8.92528122219324D+05)*X
     4 +3.10638627744347D+06)*G + W25*WW1
      WW1=WW1-0.01962D+00*G -WW2-WW3-WW4-WW5
      RETURN
  595 IF(X.GT.59.0D+00) GO TO 599
C     X=40.0 TO 59.0                             NROOTS = 5
      XXX=X**3
      G=XXX*EXP(-X)
      RT1=       (((-2.43758528330205D-02*X+2.07301567989771D+00)*X
     1 -6.45964225381113D+01)*X+7.14160088655470D+02)*G + R15/(X-R15)
      RT2=       (((-2.28861955413636D-01*X+1.93190784733691D+01)*X
     1 -5.99774730340912D+02)*X+6.61844165304871D+03)*G + R25/(X-R25)
      RT3=       (((-6.95053039285586D-01*X+5.76874090316016D+01)*X
     1 -1.77704143225520D+03)*X+1.95366082947811D+04)*G + R35/(X-R35)
      RT4=       (((-1.58072809087018D+00*X+1.27050801091948D+02)*X
     1 -3.86687350914280D+03)*X+4.23024828121420D+04)*G + R45/(X-R45)
      RT5=       (((-3.33963830405396D+00*X+2.51830424600204D+02)*X
     1 -7.57728527654961D+03)*X+8.21966816595690D+04)*G + R55/(X-R55)
      G=XXX*G
      WW5=        (( 1.35482430510942D-08*X-3.27722199212781D-07)*X
     1 +2.41522703684296D-06)*G + W55*WW1
      WW4=        (( 1.23464092261605D-06*X-3.55224564275590D-05)*X
     1 +3.03274662192286D-04)*G + W45*WW1
      WW3=        (( 1.34547929260279D-05*X-4.19389884772726D-04)*X
     1 +3.87706687610809D-03)*G + W35*WW1
      WW2=        (( 2.09539509123135D-05*X-6.87646614786982D-04)*X
     1 +6.68743788585688D-03)*G + W25*WW1
      WW1=WW1-WW2-WW3-WW4-WW5
      RETURN
C     X=59.0 TO INFINITY                         NROOTS = 5
  599 RT1=R15/(X-R15)
      RT2=R25/(X-R25)
      RT3=R35/(X-R35)
      RT4=R45/(X-R45)
      RT5=R55/(X-R55)
      WW2=W25*WW1
      WW3=W35*WW1
      WW4=W45*WW1
      WW5=W55*WW1
      WW1=WW1-WW2-WW3-WW4-WW5
      RETURN
      END
      SUBROUTINE DROOT
C
C     GENERAL SUBROUTINE TO FIND WEIGHTS AND ROOTS FOR RYS QUADRATURE.
C     CALLS DSMIT, DFUNC AND DNODE
C
C     THIS VERSION USES CHRISTOFFEL FORMULA FOR WEIGHTS.
C
C     IMPLICIT REAL*12 (A-H,O-Z)
C     DOUBLE PRECISION XX,UF,WF
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /ROOT/   XX,UF(9),WF(9),NROOTS
      COMMON /FFM/    FF(19)
      COMMON /DROOTS/  R(9,9),W(9,9)
C
      DIMENSION C(10,10),S(10,10),A(10),RT(10)
C
      DATA A0, A1S2, A1, A4 /0.0D0, 0.5D0, 1.0D0, 4.0D0/
C
C     ITH ROOT OF THE JTH RYS POLYNOMIAL IS RETURNED IN R(I,J) WITH
C     THE CORRESPONDING WEIGHT FACTOR IN W(I,J).   J=1,2,...,N
C
C
C1    WRITE(3,1) N
C1  1 FORMAT(' ##### ENTERED DROOT, NROOTS=',I2,' #####')
C
      N=NROOTS
      X=XX
      IF(N.LT.2) N=2
      N1=N+1
      NN=N+N
      CALL DFUNC(X,NN)
      DO 10 I=1,N1
      DO 10 J=1,N1
   10 S(I,J)=FF(I+J-1)
      CALL DSMIT(C,S,N1)
      DO 20 I=1,N
      DO 20 J=1,I
      W(I,J)=A0
   20 R(I,J)=A0
      WSUM=FF(1)
      W(1,1)=WSUM
      R(1,1)=FF(2)/WSUM
      DUM=DSQRT(C(2,3)**2-A4*C(1,3)*C(3,3))
      R(1,2)=A1S2*(-C(2,3)-DUM)/C(3,3)
      R(2,2)=A1S2*(-C(2,3)+DUM)/C(3,3)
      IF(N.EQ.2) GO TO 70
      DO 25 I=3,N1
   25 RT(I)=A1
      RT(1)=R(1,2)
      RT(2)=R(2,2)
      DO 60 K=3,N
      K1=K+1
      DO 30 I=1,K1
   30 A(I)=C(I,K1)
      CALL DNODE(A,RT,K)
      DO 50 I=1,K
   50 R(I,K)=RT(I)
   60 CONTINUE
   70 DO 150 K=2,N
      JMAX=K-1
      DO 150 I=1,K
      ROOT=R(I,K)
      DUM=A1/FF(1)
      DO 110 J=1,JMAX
      J1=J+1
      POLY=C(J1,J1)
      DO 100 M=1,J
  100 POLY=POLY*ROOT+C(J1-M,J1)
  110 DUM=DUM+POLY*POLY
  150 W(I,K)=A1/DUM
      DO 160 K=1,NROOTS
      DUM=R(K,NROOTS)
      UF(K)=DUM/(A1-DUM)
  160 WF(K)=W(K,NROOTS)
C
C1    WRITE (3,2) (UF(IQ),WF(IQ),IQ=1,NROOTS)
C1  2 FORMAT (' ROOTS AND WEIGHTS',/,(1X,2G25.10))
C1    WRITE (3,3)
C1  3 FORMAT (' ##### EXITING DROOT ###########################')
C
      RETURN
      END
      SUBROUTINE DSMIT(C,S,N)
C
C     RETURNS AN N*N TRIANGULAR MATRIX C SUCH THAT C(TRANSPOSE)SC=I,
C     WHERE I IS A N*N IDENTITY MATRIX
C
C
C     IMPLICIT REAL*12 (A-H,O-Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION C(10,10),S(10,10),V(10),Y(10)
C
      DATA A0, A1 /0.0D0, 1.0D0/
C
      DO 10 I=1,N
      DO 10 J=1,I
   10 C(I,J)=A0
      DO 100 J=1,N
      KMAX=J-1
      FAC=S(J,J)
      IF(KMAX.EQ.0) GO TO 60
      DO 20 K=1,KMAX
      V(K)=A0
   20 Y(K)=S(K,J)
      DO 50 K=1,KMAX
      DOT=A0
      DO 30 I=1,K
   30 DOT=C(I,K)*Y(I)+DOT
      DO 40 I=1,K
   40 V(I)=V(I)-DOT*C(I,K)
   50 FAC=FAC-DOT*DOT
   60 FAC=A1/DSQRT(FAC)
      C(J,J)=FAC
      IF(KMAX.EQ.0) GO TO 100
      DO 70 K=1,KMAX
   70 C(K,J)=FAC*V(K)
  100 CONTINUE
      RETURN
      END
      SUBROUTINE DFUNC(X,N)
C
C
C
C
C
C     IMPLICIT REAL*12 (A-H,O-Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      COMMON /FFM/    FF(19)
C
      DATA TM78, TM29, A1S2, PIE4, A1, A2, A80, A360 /1.0D-78, 1.0D-15,
     1   0.5D0,0.7853981633974483096156608D0,1.0D0,2.0D0,80.0D0,360.0D0/
C
      TOL=TM29
      XX=X+X
      FACMIN=XX
      E=TM78
      QX=-X
      IF(FACMIN.LT.A360) E=DEXP(QX)
      IF (E.EQ.0.0D+00) E=TM78
      IF(FACMIN.GT.A80) GO TO 100
      TERM=A1
      SUM=A1
      FAC=N
      FAC=FAC+A1S2
   10 FAC=FAC+A1
      TERM=TERM*X/FAC
      SUM=SUM+TERM
      IF(FAC.LE.FACMIN) GO TO 10
      T=TERM
      S=SUM
      IF(T.GT.S*TOL) GO TO 10
      FAC=N+N+1
      FF(N+1)=SUM*E/FAC
      M=N-1
      FAC=M+M+1
   20 IF(M.LT.0) RETURN
      FF(M+1)=(E+XX*FF(M+2))/FAC
      M=M-1
      FAC=FAC-A2
      GO TO 20
C
C     USE ASYMPTOTIC EXPANSION FOR LARGE ARGUMENTS.
C
  100 A=DSQRT(PIE4/X)
      TMAX=A*TOL/E
      TERM=A1/XX
      SUM=TERM
      FAC=A1
  110 FAC=FAC-A2
      TERM=FAC*TERM/XX
      SUM=TERM+SUM
      T=TERM
      IF(DABS(T).GT.TMAX) GO TO 110
      FF(1)=A-E*SUM
      FAC=-A1
      M=0
  120 IF(M.EQ.N) RETURN
      M=M+1
      FAC=FAC+A2
      FF(M+1)=(FAC*FF(M)-E)/XX
      GO TO 120
      END
      SUBROUTINE DNODE(A,RT,K)
C
C
C     RETURNS IN RT(I) THE ITH ROOT OF A POLYNOMIAL OF ORDER K WHOSE
C     MTH COEFFICIENT IS STORED IN A(M+1). IT IS ASSUMED THAT THE
C     INITIAL VALUES IN RT BRACKET THE FINAL VALUES.
C
C
C     IMPLICIT REAL*12 (A-H,O-Z)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DIMENSION A(10),RT(10)
C
      DATA A0,TM21,A1S16,A1S4,A3S4 /0.0D0,1.0D-13,6.25D-2,0.25D0,0.75D0/
C
      TOL=TM21
      K1=K+1
      R2=A0
      P2=A(1)
      DO 100 M=1,K
      R1=R2
      P1=P2
      R2=RT(M)
      P2=A(K1)
      DO 10 I=1,K
   10 P2=P2*R2+A(K1-I)
      PROD=P1*P2
      IF(PROD.LT.A0) GO TO 20
      WRITE(6,15) M,K
   15 FORMAT(/12H0ROOT NUMBER,I4,
     1   38H WAS NOT FOUND FOR POLYNOMIAL OF ORDER,I4//)
      STOP
   20 R5=R1
      P5=P1
      R6=R2
      P6=P2
   30 R3=R5
      P3=P5
      R4=R6
      P4=P6
      R =(R3*P4-R4*P3)/(P4-P3)
      DR=R4-R3
      DELTA=DR
      IF(DABS(DELTA).LT.TOL) GO TO 90
      DR=A1S16*DR
      R5=R-DR
      IF(R5.LT.R3) R5=R3
      R6=R+DR
      IF(R6.GT.R4) R6=R4
      P5=A(K1)
      P6=P5
      DO 40 I=1,K
      P5=P5*R5+A(K1-I)
   40 P6=P6*R6+A(K1-I)
   45 PROD=P5*P6
      IF(PROD.LT.A0) GO TO 30
      PROD=P3*P5
      IF(PROD.GT.A0) GO TO 60
      R5=A1S4*R3+A3S4*R5
      P5=A(K1)
      DO 50 I=1,K
   50 P5=P5*R5+A(K1-I)
      GO TO 45
   60 R6=A1S4*R4+A3S4*R6
      P6=A(K1)
      DO 70 I=1,K
   70 P6=P6*R6+A(K1-I)
      GO TO 45
   90 RT(M)=R
  100 CONTINUE
      RETURN
      END
      SUBROUTINE ISREW(ITAPE)
C
      IMPLICIT INTEGER (A-Z)
C
      COMMON /POINTR/ POS(99)
C
C THIS ROUTINE REWINDS FILES USING THE FORTRAN UTILITIES DIRECT ACCESS
C
C PLACE THE POINTER OF ITAPE AT THE BEGINNING
C
      IOP = 2
      IPOS = 1
      IARR = 0
      NW = 0
c      CALL STPOS(ITAPE,IPOS,IERR)
      CALL RWFIL(ITAPE,IOP,IPOS,IARR,NW,IERR)
      IF(IERR.NE.0) THEN
      WRITE(*,*) ' ERROR ENCOUNTERED REWINDING FILE ',ITAPE
      WRITE(*,*) ' IERR = ',IERR
      STOP
      END IF
C
C CHECK CURRENT POSITION OF POINTER
C
c      CALL RDPOS(ITAPE,IPOS,IERR)
      IF(IERR.NE.0.OR.IPOS.NE.1) THEN
      WRITE(*,*) ' ERROR ENCOUNTERED READING POINTER FOR FILE ',ITAPE
      WRITE(*,*) ' IERR,IPOS = ',IERR,IPOS
      STOP
      END IF
C
C UPDATE THE POSITION ARRAY
C
      POS(ITAPE) = IPOS
C
C RETURN
C
      RETURN
      END
      SUBROUTINE RCLOSE(ITAPE,JCODE)
C
      IMPLICIT INTEGER (A-Z)
C
C JCODE = 4     CLOSE AND DELETE FILE
C JCODE = 3     CLOSE AND SAVE FILE
C
CTJL  WRITE(*,*) ' CLOSING FILE',ITAPE,'  JCODE = ',JCODE
c      IF(JCODE.NE.3.AND.JCODE.NE.4) THEN
c      WRITE(*,*) ' INVALID JCODE IN RCLOSE,  JCODE = ',JCODE
c      WRITE(*,*) ' FILE ',ITAPE,'  CLOSED AND SAVED.'
c      JCODE = 3
c      END IF
      IF(ITAPE.EQ.6) STOP ' YOU CANNOT CLOSE A FILE ON UNIT 6'
       if(jcode.eq.3)then
Cibm
C      close(unit=itape,
C    1         disp='save',
C    2        err=300)
C      else if(jcode.eq.4)then
       close(unit=itape, status='KEEP', err=300)
Cibm
       else if(jcode.eq.4)then
Cibm
C      close(unit=itape,
C    1        disp='delete',
C    2        err=300)
       close(unit=itape, status='DELETE', err=300)
Cibm
       else
       write(*,*)'file','itape',' incorrect j code; no close'
       end if
       go to 400
c      CALL UCLOSE(ITAPE,JCODE,IERR)
300   continue
      WRITE(*,*) ' ERROR ENCOUNTERED CLOSING FILE',ITAPE
      WRITE(*,*) ' IERR,JCODE = ',ERR,JCODE
CTJL  STOP
400   continue
      RETURN
      END
      BLOCK DATA IOBUF
C
      IMPLICIT INTEGER (A-Z)
C
      COMMON /IOBUFF/ BUFF(1024)
C
      DATA BUFF / 1024*0 /
C
      END
      SUBROUTINE MABORT
C
      IABORT = 999
      REWIND IABORT
C
      STOP
C
      END

            SUBROUTINE LOCATOR(INPUT,TOKEN,IERROR)
C
C     ----- SEARCH THROUGH INPUT FILE FOR TOKEN BEGINNING
C           WITH # TO LOCATE INPUT FOR PROGRAM.  IERROR IS
C           SET TO 0 IF NO ERRORS, 1 IF ANY ERROR OCCURS.
C
C
      CHARACTER*10 TOKEN,LINE
C
      REWIND (UNIT=INPUT,ERR=99)
C
    1 CONTINUE
      READ (UNIT=INPUT,FMT='(A10)',END=99,ERR=99) LINE
      IF (LINE.NE.TOKEN) GO TO 1
C
      IERROR=0
      RETURN
C
C
   99 CONTINUE
      IERROR=1
      RETURN
      END
      subroutine dump
cbhl
cbhl for IBM compatiblility
cbhl
      call exit(1)
c
      return
      end
