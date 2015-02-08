      SUBROUTINE DEFIL (NDEFL,LDAR,LDSK,NFDA,NALM)                
      DIMENSION LDSKS(10),LDARS(2),IWORK(5)                       
      character*1 nc(3)
      character*8 fyle
      data fyle/'fort.   '/
      COMMON IDSK
      DATA NFLG/0/                                                
      DATA LDSKS/1940,4490,8840,15380,21910,                      
     1    3850,9410,17170,26570,30630/                            
      DATA LDARS/1058,722/                                        
    4 N=NDEFL                                                     
      IF (N.EQ.0) GO TO 15                                        
      IF (N.LE.10) GO TO 10                                       
      WRITE (6,5)                                                 
    5 FORMAT ('0NDEFL OUT OF RANGE')                              
      NALM=1                                                      
      RETURN                                                      
C                                                                 
   10 M=1                                                         
      IF ( N.GT.5) M=2                                            
      LDSK=LDSKS(N)                                               
      LDAR=LDARS(M)                                               
   15 IF (NFLG.NE.0) GO TO 30                                     
CDRY      NFLG=1                                                      
   20 continue
      call tochar(NFDA,nc)
      fyle(6:6)=nc(1)
      fyle(7:7)=nc(2)
      fyle(8:8)=nc(3)
      IDSK = 1
      LRECL = 4 * LDAR
      WRITE(6,1000)NFDA,LDAR,LRECL,FYLE
C      OPEN ( UNIT=NFDA, ACCESS='DIRECT', RECL=4*LDAR, file='fyle')
      OPEN ( UNIT=NFDA, ACCESS='DIRECT', RECL=4*LDAR)
C                                                                  
   30 NALM=0                                                       
      RETURN                                                       
 1000 FORMAT(5X,'OPEN DIRECT ACCESS UNIT=',I3,' LDAR=',I8,' LRECL=',
     X       I8, ' FILE=',A8)
      END                                                          
      SUBROUTINE MATOUT(A,NAD,NBD,M,N,IOUT)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*6 LINE(21)
      DIMENSION A(NAD,NBD)
    1 FORMAT(3X,5(7X,I5))
    2 FORMAT(3X,21A6)
    3 FORMAT(2X,I3,2X,5F12.7)
    4 FORMAT(/)
C
      DATA LINE / 21*'------' /
      II=0
      JJ=0
  200 II=II+1
      JJ=JJ+1
      KK=5*JJ
      NN=N
      IF(N.GT.KK) NN=KK
      LL=2*(NN-II+1)+1
      WRITE (IOUT,1) (I,I=II,NN)
      WRITE (IOUT,2) (LINE(I),I=1,LL)
      DO 101 I=1,M
      WRITE (IOUT,3) I,(A(I,J),J=II,NN)
  101 CONTINUE
      IF(N.LE.KK) GO TO 201
      WRITE (IOUT,4)
      II=KK
      GO TO 200
  201 RETURN
      END
      SUBROUTINE PRINT(A,NAD,M,IOUT)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*6 LINE(21)
      DIMENSION A(NAD)
    1 FORMAT(2X,10(7X,I5))
    2 FORMAT(2X,21A6)
    3 FORMAT(2X,I2,2X,10F12.7)
    4 FORMAT(/)
      DATA LINE / 21*'------' /
      II=0
      JJ=0
  200 II=II+1
      JJ=JJ+1
      KK=10*JJ
      NN=KK+KK*(KK-1)/2
      MM=M
      IF(M.GT.KK) MM=KK
      LL=2*(MM-II+1)+1
      WRITE (IOUT,1) (I,I=II,MM)
      WRITE (IOUT,2) (LINE(I),I=1,LL)
      DO 101 I=II,M
      I1=I*(I-1)/2+II
      I2=I+I*(I-1)/2
      IF(I2.GT.NN) I2=I1+9
      WRITE (IOUT,3) I,(A(J),J=I1,I2)
  101 CONTINUE
      IF(M.LE.KK) GO TO 201
      WRITE (IOUT,4)
      II=KK
      GO TO 200
  201 RETURN
      END
CDRY  5/15/95 IBM PROBLEM  WITH NAME PACK 
C     CHANGE IT TO XPACK
CDRY      SUBROUTINE PACK(I,J,K,L,IX,JA,JAA)
      SUBROUTINE XPACK(I,J,K,L,IX,JA,JAA)
C
      IF(I-L) 202,201,202
C   (11/11)
  201 IX=1
      GO TO 230
  202 IF(J-L) 208,203,211
  203 IF(K-L) 205,204,205
C   (21/11)=(12/11)=(11/21)=(11/12)
  204 IX=4
      GO TO 230
  205 IF(I-K) 207,206,207
C   (21/21)=(21/12)=(12/21)=(12/12)
  206 IX=5
      GO TO 230
C   (31/21)=(31/12)=(13/21)=(13/12)=(21/31)=(21/13)=(12/31)=(12/13)
  207 IX=11
      GO TO 230
  208 IF(K-L) 210,209,210
C   (31/22)=(13/22)=(22/31)=(22/13)
  209 IX=8
      GO TO 230
C   (41/32)=(41/23)=(14/32)=(14/23)=(32/41)=(32/14)=(23/41)=(23/14)
  210 IX=14
      GO TO 230
  211 IF(I-J) 212,212,217
  212 IF(I-K) 213,213,214
C   (22/21)=(22/12)=(21/22)=(12/22)
  213 IX=3
      GO TO 230
  214 IF(K-L) 215,215,216
C   (22/11)=(11/22)
  215 IX=2
      GO TO 230
C   (33/21)=(33/12)=(21/33)=(12/33)
  216 IX=6
      GO TO 230
  217 IF(I-K) 218,218,219
C   (32/31)=(32/13)=(23/31)=(23/13)=(31/32)=(31/23)=(13/32)=(13/23)
  218 IX=9
      GO TO 230
  219 IF(J-K) 220,221,222
C   (42/31)=(42/13)=(24/31)=(24/13)=(31/42)=(31/24)=(13/42)=(13/24)
  220 IX=13
      GO TO 230
C   (32/21)=(32/12)=(23/21)=(23/12)=(21/32)=(21/23)=(12/32)=(12/23)
  221 IX=10
      GO TO 230
  222 IF(K-L) 223,223,224
C   (32/11)=(23/11)=(11/32)=(11/23)
  223 IX=7
      GO TO 230
C   (43/21)=(43/12)=(34/21)=(34/12)=(21/43)=(21/34)=(12/43)=(12/34)
  224 IX=12
C
  230 CONTINUE
      JA=4096*I+16*J+IX
      JAA=256*K+L
      RETURN
C
      ENTRY UNPACK(I,J,K,L,IX,JA,JAA)
      IX=MOD(JA,16)
      JC=JA/16
      J=MOD(JC,256)
      I=JC/256
      L=MOD(JAA,256)
      K=JAA/256
      RETURN
      END
      SUBROUTINE ZERO(A,N)
C
      REAL*8 A(N)
C
      NBYTES = 8 * N
c      CALL VXINIT(A,0,NBYTES)
      DO 1 I=1,N
           A(I)=0.0D+00
    1 CONTINUE
C
      RETURN
      END
      INTEGER FUNCTION INTOWP(N)
      INTOWP=2*N
      RETURN
      END
      SUBROUTINE MATWRTT(A,N,IOUT,TOL)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(N),IVAL(50),VAL(50)
      K=0
      DO 10 I=1,N
      IF(DABS(A(I)).LT.TOL)GO TO 10
      IF(K.GE.50)THEN
      WRITE(IOUT,1000)(IVAL(J),VAL(J),J=1,K)
      K=0
      ENDIF
      K=K+1
      VAL(K)=A(I)
      IVAL(K)=I
   10 CONTINUE
      IF(K.LE.0)RETURN
      WRITE(IOUT,1000)(IVAL(J),VAL(J),J=1,K)
      RETURN
 1000 FORMAT(6(1X,I6,E13.6))
      END
      SUBROUTINE PRTSM(A,N,IUNIT)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(N,N)
      WRITE(IUNIT,1001)N
      DO 10 I=1,N
      WRITE(IUNIT,1000)I,(A(I,J),J=1,N)
   10 CONTINUE
      RETURN
 1000 FORMAT(2X,'ROW=',I3,/(10F12.6))
 1001 FORMAT(/10X,'SQUARE MATRIX: DIMENSION=',I4)
      END
      SUBROUTINE PRVOM(X,LEN,IFLAG,M6)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION X(2)
      IF(IFLAG.NE.1)GO TO 50
C     WRITE(M6,10)
C  10 FORMAT(/,'   VECTOR OUTPUT ',/)
      WRITE(M6,20)(K,X(K),K=1,LEN)
   20 FORMAT(5(4X,I4,2X,D15.7))
      RETURN
   50 CONTINUE
C     WRITE(M6,60)
C  60 FORMAT(/,'  MATRIX OUTPUT ',/)
      IF(IFLAG.EQ.3)GO TO 80
      LX=0
      DO 70 K=1,LEN
      WRITE(M6,75) K,(X(LX+L),L=1,K)
      LX=LX+K
   70 CONTINUE
   75 FORMAT(/,'  ROW ',I4,/,5(2X,D15.7))
      RETURN
   80 CONTINUE
      LX=0
      DO 85 K=1,LEN
      KX=LEN+1-K
      WRITE(M6,75) K,(X(LX+L),L=1,KX)
      LX=LX+KX
   85 CONTINUE
      RETURN
      END
      SUBROUTINE REPOSN(LU,A,ITP6)
      IMPLICIT REAL*8(A-H,O-Z)
   10 CONTINUE
      READ(LU,1000,END=3)B
      IF(A.NE.B)GO TO 10
      RETURN
    3 CONTINUE
      WRITE(ITP6,4)A,LU
    4 FORMAT(5X,A8,18H NOT FOUND ON UNIT     ,I3)
 1000 FORMAT(A8)
      STOP
      END
      SUBROUTINE RZERO(A,N)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(N)
      CALL ZERO(A,N)
      RETURN
      END
      SUBROUTINE SEARCH(LU,A,ITP6)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IB(8),B(4),IA(2)
      EQUIVALENCE (IB(1),B(1)),(IA(1),AL)
      DATA IC/4H****/
      AL=A
    1 READ (LU,END=3) IB
C     IF(B(1).NE.C.OR.B(4).NE.A) GO TO 1
      IF((IC.NE.IB(1)).OR.(IC.NE.IB(2))
     1  .OR.(IB(7).NE.IA(1)).OR.(IB(8).NE.IA(2)))
     2   GO TO 1
      RETURN
    3 CONTINUE
      WRITE(ITP6,4)A,LU
    4 FORMAT (5X,A8,18H NOT FOUND ON UNIT , I3)
      STOP   ' ERROR IN SEARCH '
      END
      SUBROUTINE TOPEN(N,IOPEN,IRECL,IBLK,IUNIT,M6)
      DOUBLE PRECISION WRD(2)
      DIMENSION IBLK(N),IOPEN(N),IRECL(N),IUNIT(N)
      DIMENSION IDEFLT(30)
      DATA WRD/8H PACKED ,8H BINARY  /
 1000 FORMAT(3X,'UNIT=',I3,' OPENED IN',A8,' MODE  LRECL= ',I5 ,
     1 ' BLOCKSIZE= ',I5,' BYTES' )
 1001 FORMAT(/15X, '  FILE DEFINITIONS  ',/ )
 1002 FORMAT(/5X,'P-E DEFAULT RECL AND BLOCKSIZE USED FOR PACKED ',
     1       'UNITS ',30I4)
      KDEFLT=0
      WRITE(M6,1001)
      DO 10 I=1,N
      IF(IOPEN(I)*IUNIT(I).EQ.0)GO TO 10
      NRECL=IRECL(I)*4
      NBLK=IBLK(I)*4
      IWRD=1
      IF(IOPEN(I).LT.0)GO TO 50
      IF(IRECL(I).EQ.0)GO TO 30
      IF(NBLK.EQ.0)NBLK=NRECL+((NRECL-1)/256+1)*4
C     OPEN(UNIT=IUNIT(I),ACCESS='SEQUENTIAL',FORM='PACKED',
C    1     RECL=NRECL,BLOCKSIZE=NBLK )
      OPEN(UNIT=IUNIT(I),ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      GO TO 9
   30 CONTINUE
      IF(NBLK.EQ.0)GO TO 40
C     OPEN(UNIT=IUNIT(I),ACCESS='SEQUENTIAL',FORM='PACKED',
C    1      BLOCKSIZE=NBLK)
      OPEN(UNIT=IUNIT(I),ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      GO TO 9
   40 CONTINUE
      OPEN(UNIT=IUNIT(I),ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      KDEFLT=KDEFLT+1
      IDEFLT(KDEFLT)=IUNIT(I)
      GO TO 10
   50 CONTINUE
      IWRD=2
C     OPEN(UNIT=IUNIT(I),ACCESS='SEQUENTIAL',FORM='BINARY',
C    1     RECL=NRECL,BLOCKSIZE=NRECL)
      OPEN(UNIT=IUNIT(I),ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      NBLK=NRECL
    9 CONTINUE
      WRITE(M6,1000)IUNIT(I),WRD(IWRD),NRECL,NBLK
   10 CONTINUE
      IF(KDEFLT.GT.0)WRITE(M6,1002)(IDEFLT(I),I=1,KDEFLT)
      RETURN
      END
      SUBROUTINE TOCHAR(N,NC)
      CHARACTER*1 NC(3),IX(10)
      DATA IX/'0','1','2','3','4','5','6','7','8','9'/
      I=N/100
      NC(1)=IX(I+1)
      J=N-100*I
      K=J/10
      NC(2)=IX(K+1)
      L=N-100*I-10*K
      NC(3)=IX(L+1)
      RETURN
      END
      SUBROUTINE TRNSQ(C,A,T,S,NAO,NMO)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(NAO,NAO),C(NAO,NMO),T(NMO,NMO),S(2)
      DATA ZERO/0.0D0/
      DO 5 I=1,NMO
      DO 5 J=1,NMO
    5 T(I,J)=ZERO
      DO 10 I=1,NAO
C         CONSTRUCT PK(I)
      DO 20 J=1,NMO
   20 S(J)=ZERO
      DO 30 K=1,NMO
      DO 30 J=1,NAO
   30 S(K)=S(K)+A(J,I)*C(J,K)
C         FORM PK(I)*C(I,L) CONTRI TO T(K,L)
      DO 40 K=1,NMO
      DO 40 L=1,NMO
   40 T(K,L)=T(K,L)+S(K)*C(I,L)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE YREADW(IUNIT,A,LEN)
      CHARACTER*4 A(LEN)
CDRY      WRITE(6,9999)IUNIT,LEN
C9999 FORMAT(5X,'YREAD: IUNIT LEN ',I3,I7)
      READ(IUNIT)A
      RETURN
      END
