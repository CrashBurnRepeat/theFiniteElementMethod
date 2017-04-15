C    *****************************************************************
C    *                         MESH-1D.FOR                           *
C    *                                                               *
C    *              GENERATES MESH AND BOUNDARY CONDITIONS           *
C    *                 FOR 2, 3, OR 4 NODE 1-D ELEMENTS              *
C    *                                                               *
C    *                     WRITTEN FOR THE TEXTBOOK                  *
C    *           "THE FINITE ELEMENT METHOD: BASIC CONCEPTS          *
C    *                        AND APPLICATIONS"                      *
C    *                               BY                              *
C    *             DARRELL W. PEPPER AND JUAN C. HEINRICH            *
C    *                 HEMISPHERE PUBLISHING COMPANY                 *
C    *                                                               *
C    *       (NOTE: MAXIMUM NUMBER OF NODES DIMENSIONED TO 51)       *
C    *                                                               *
C    *                                   NOV. 2, 1991 D. W. PEPPER   *
C    *****************************************************************
C
      CALL INFILE
      CALL GRID
      CALL OUTPUT
      CALL BANDW
      STOP
      END
C -----------------------------------------------------------------
      SUBROUTINE GRID
C    *****************************************************************
C    *                                                               *
C    *   SUBROUTINE GRID:   READ INPUT DATA AND GENERATE MESH        *
C    *                                                               *
C    *****************************************************************
C
      COMMON/A/ X(51),QQ(51),Q(51),U(51),DX(51)
      COMMON/B/ NELEM,NNODE,NUMN,NINT,NTYPE,NSTOP,KPRNT,NVEL
      COMMON/C/ DT,AF,KBC1,KBC2,KBC3,TO,RHO,CP
      COMMON/D/ NET(2),NEQ(2),NEH(2),TEMP(2),TINF(2),H(2)
      COMMON/E/ XMIN,XMAX
      COMMON/F/ NODE(51,4)
      COMMON/G/ TITLE(18)
      DIMENSION XX(51)
C
C     READ INITIAL DATA
C
      READ(5,100)TITLE
  100 FORMAT(18A4)
      READ(5,*)XMIN,XMAX
      READ(5,*)NTYPE,NUMN,NINT,NSTOP,KPRNT,NVEL
      READ(5,*)DT,AF,TO,RHO,CP,DXX,UX
C
      IF(NUMN.LT.2.OR.NUMN.GT.4)THEN
      WRITE(*,101)
  101 FORMAT(10X,'ELEMENT SHAPE NOT PROPERLY DEFINED')
      STOP
      ENDIF
C
      CALL CLEAR(X,51)
      CALL CLEAR(XX,51)
      NELEM=NINT
      XX(1)=XMIN
      X(1)=XX(1)
      N=NINT+1
      DO 1 I=2,N
    1 XX(I)=XX(I-1)+(XMAX-XMIN)/NINT
      NM1=NUMN-1
      GOTO(11,12,13)NM1
C
C     LINEAR ELEMENTS
C
   11 NNODE=N
      DO 2 I=1,NNODE
    2 X(I)=XX(I)
      DO 3 J=1,NELEM
      NODE(J,1)=J
    3 NODE(J,2)=J+1
      GOTO 20
C
C     QUADRATIC ELEMENTS
C
   12 NNODE=2*N-1
      DO 5 I=2,N
      J=2*I-1
    5 X(J)=XX(I)
      DO 6 I=2,NNODE,2
    6 X(I)=X(I-1)+(X(I+1)-X(I-1))/2.
      DO 7 J=1,NELEM
      K=2*J-1
      NODE(J,1)=K
      NODE(J,2)=K+1
    7 NODE(J,3)=K+2
      GOTO 20
C
C     CUBIC ELEMENTS
C
   13 NNODE=3*N-2
      DO 8 I=2,N
      J=3*I-2
    8 X(J)=XX(I)
      DO 9 I=2,NNODE,3
    9 X(I)=X(I-1)+(X(I+2)-X(I-1))/3.
      DO 10 I=3,NNODE,3
   10 X(I)=X(I-2)+2.*(X(I+1)-X(I-2))/3.
      DO 14 J=1,NELEM
      K=3*J-2
      NODE(J,1)=K
      NODE(J,2)=K+1
      NODE(J,3)=K+2
   14 NODE(J,4)=K+3
C
   20 CONTINUE
C
      READ(5,*)KBC1
      IF(KBC1.GE.1)THEN
      DO 15 KK=1,KBC1
   15 READ(5,*)NET(KK),TEMP(KK)
      ENDIF
C
      READ(5,*)KBC2
      IF(KBC2.GE.1)THEN
      DO 16 KK=1,KBC2
   16 READ(5,*)NEQ(KK),Q(KK)
      ENDIF
C
      READ(5,*)KBC3
      IF(KBC3.GE.1)THEN
      DO 17 KK=1,KBC3
   17 READ(5,*)NEH(KK),H(KK),TINF(KK)
      ENDIF
C
      IF(NVEL.EQ.1)THEN
      READ(5,*)(I,U(I),J=1,NNODE)
      ELSE
      DO 18 I=1,NNODE
   18 U(I)=UX
      ENDIF
C
      DO 55 J=1,NELEM
   55 DX(J)=DXX
      RETURN
      END
C -----------------------------------------------------------------
      SUBROUTINE OUTPUT
C    *****************************************************************
C    *                                                               *
C    *   SUBROUTINE OUTPUT: WRITE OUT MESH DATA TO FILE              *
C    *                                                               *
C    *****************************************************************
      COMMON/A/ X(51),QQ(51),Q(51),U(51),DX(51)
      COMMON/B/ NELEM,NNODE,NUMN,NINT,NTYPE,NSTOP,KPRNT,NVEL
      COMMON/C/ DT,AF,KBC1,KBC2,KBC3,TO,RHO,CP
      COMMON/D/ NET(2),NEQ(2),NEH(2),TEMP(2),TINF(2),H(2)
      COMMON/E/ XMIN,XMAX
      COMMON/F/ NODE(51,4)
      COMMON/G/ TITLE(18)
      CHARACTER*4 STOP(3),WORD
      DATA STOP/'DIRC','FLUX','CONV'/
      WRITE(*,100)
  100 FORMAT('',80('*'))
C
      WRITE(15,'(5X,18A4)')TITLE
      WRITE(*,'(5X,18A4)')TITLE
C
C     WRITE(15,1004)XMIN,XMAX
      WRITE(*,1004)XMIN,XMAX
      MTYPE=0
      NUMDIM=1
      IAXI=0
C
C     WRITE THE NUMBER OF SYSTEM NODES AND ELEMENTS
C
      WRITE(15,1005)MTYPE,NUMDIM,NNODE,NELEM,NUMN,NSTOP,KPRNT,NVEL
      WRITE(*,1005)NNODE,NELEM,NUMN,NSTOP,KPRNT,NVEL
      WRITE(15,1006)NTYPE,DT,AF,TO,RHO,CP,IAXI
      WRITE(*,1006)NTYPE,DT,AF,TO,RHO,CP
C
C     WRITE NODE NUMBERS AND NODAL COORDINATES
C
      DO 104 I=1,NNODE
      WRITE(15,1007)I,X(I)
      WRITE(*,1007)I,X(I)
  104 CONTINUE
C
      DO 105 J=1,NELEM
      WRITE(15,1010)J,QQ(J),DX(J),(NODE(J,MM),MM=1,NUMN)
      WRITE(*,1010)J,QQ(J),DX(J),(NODE(J,MM),MM=1,NUMN)
  105 CONTINUE
C
C     CYCLE FOR EACH NODE HAVING SPECIFIED VALUE
C
      DO 110 I=1,KBC1
      IF(NET(I).EQ.1)THEN
      WRITE(15,1015)NODE(1,1),TEMP(I)
      WRITE(*,1015)NODE(1,1),TEMP(I)
      ELSE
      IF(NET(I).EQ.2)THEN
      WRITE(15,1015)NODE(NELEM,NUMN),TEMP(I)
      WRITE(*,1015)NODE(NELEM,NUMN),TEMP(I)
      ENDIF
      ENDIF
  110 CONTINUE
      WRITE(15,'(6X,A4)')STOP(1)
      WRITE(*,'(6X,A4)')STOP(1)
C
C  WRITE NODE NUMBER THAT DEFINES BOUNDARY
C  WHERE FLUX (GRADIENT) IS SPECIFIED
C
      DO 120 I=1,KBC2
      IF(NEQ(I).EQ.1)THEN
      NEL=1
      WRITE(15,1020)Q(I),NEL,NEL
      WRITE(*,1020)Q(I),NEL,NEL
      ELSE
      IF(NEQ(I).EQ.2)THEN
      NEL=2
      WRITE(15,1020)Q(I),NELEM,NEL
      WRITE(*,1020)Q(I),NELEM,NEL
      ENDIF
      ENDIF
  120 CONTINUE
      WRITE(15,'(6X,A4)')STOP(2)
      WRITE(*,'(6X,A4)')STOP(2)
C
      DO 130 I=1,KBC3
      IF(NEH(I).EQ.1)THEN
      NEL=1
      WRITE(15,1025)H(I),TINF(I),NEL,NEL
      WRITE(*,1025)H(I),TINF(I),NEL,NEL
      ELSE
      IF(NEH(I).EQ.2)THEN
      NEL=2
      WRITE(15,1025)H(I),TINF(I),NELEM,NEL
      WRITE(*,1025)H(I),TINF(I),NELEM,NEL
      ENDIF
      ENDIF
  130 CONTINUE
      WRITE(15,'(6X,A4)')STOP(3)
      WRITE(*,'(6X,A4)')STOP(3)
C
C     FORMAT STATEMENTS
C
 1004 FORMAT(2X,2F8.3)
 1005 FORMAT(8I5)
 1006 FORMAT(I5,5F8.4)
 1007 FORMAT(4X,I6,F8.4)
 1010 FORMAT(I4,1X,F7.3,2X,F7.3,2X,4I4)
 1015 FORMAT(10X,I5,5X,F10.5)
 1020 FORMAT(10X,F10.5,2I5)
 1025 FORMAT(10X,2F10.5,2I5)
      RETURN
      END
C -----------------------------------------------------------------
      SUBROUTINE BANDW
C    *****************************************************************
C    *                                                               *
C    *   SUBROUTINE BANDW: BANDWIDTH FOR 1-D MESH IS EQUAL TO NUMBER *
C    *                     OF LOCAL NODES WITHIN THE ELEMENT         *
C    *                                                               *
C    *****************************************************************
      COMMON/B/ NELEM,NNODE,NUMN,NINT,NTYPE,NSTOP,KPRNT,NVEL
      NBANDW=NUMN
      WRITE(*,100)NBANDW
  100 FORMAT(2X,'THE BANDWIDTH IS ',I2,' WIDE')
      RETURN
      END
C -----------------------------------------------------------------
      SUBROUTINE CLEAR(A,N)
C    *****************************************************************
C    *                                                               *
C    *   SUBROUTINE CLEAR: CLEAR ARRAY                               *
C    *                                                               *
C    *****************************************************************
      DIMENSION A(N)
      DO 1 I=1,N
    1 A(I)=0.0
      RETURN
      END
C -----------------------------------------------------------------
      SUBROUTINE INFILE
C    *****************************************************************
C    *                                                               *
C    *   SUBROUTINE INFILE: INPUT FILE NAMES AND SET UP UNIT NUMBERS *
C    *                      FOR READ AND WRITE                       *
C    *                                                               *
C    *****************************************************************
      COMMON/I/ NELEM,NNODE,NUMN,NINT,NTYPE,NSTOP,DT,AF,KBC1,KBC2
      CHARACTER INFIL*12,OUTFIL*12
      DATA INFIL/'            '/,OUTFIL/'            '/
      WRITE(*,105)
  105 FORMAT(1X,'IF YOU WANT HARDCOPY OUTPUT, PRESS CTL-PRTSC BEFORE ENT
     1ERING INPUT FILE NAME')
      WRITE(*,'(/1X,A\)')' ENTER INPUT FILE NAME (e.g.,INPUTFIL.DAT): '
      READ(*,'(BN,A)')INFIL
      WRITE(*,'(/1X,A\)')' ENTER OUTPUT FILE NAME (e.g.,OUTFILE.DAT): '
      READ(*,'(BN,A)')OUTFIL
      OPEN(5,FILE=INFIL)
C
C     OPEN UNIT 15 FOR OUTPUT
C
      OPEN(15,FILE=OUTFIL,STATUS='NEW')
      WRITE(*,104)
  104 FORMAT(/,1X,'PLEASE WAIT, GRID IS DEVELOPING',/)
      RETURN
      END
