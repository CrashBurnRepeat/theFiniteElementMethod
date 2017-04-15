C  *********************************************************************
C  **                                                                 **
C  **                          MESH-2D.FOR                            **
C  **                     A 2-D MESH PROGRAM WHICH                     **
C  **                     ACCOMPANIES THE TEXTBOOK                    **
C  **                                                                 **
C  **                "THE FINITE ELEMENT METHOD: BASIC                **
C  **                    CONCEPTS AND APPLICATIONS"                   **
C  **                                BY                               **
C  **              DARRELL W. PEPPER AND JUAN C. HEINRICH             **
C  **                  HEMISPHERE PUBLISHING COMPANY                  **
C  **                                                                 **
C  **              GENERATES MESH AND BOUNDARY CONDITIONS             **
C  **               FOR 3, 4, 6, OR 8 NODED 2-D ELEMENTS              **
C  **                                                                 **
C  **    (NOTE: THIS PROGRAM DIMENSIONED TO A MAXIMUM OF 501 NODES)   **
C  **                                                                 **
C  **                                  NOV. 2, 1991 D. W. PEPPER      **
C  *********************************************************************
C                             SIDE=1
C                        ******************
C                        *                *
C                        *                *
C                  SIDE=4*  ELEMENT SIDES *SIDE=2
C                        *                *
C                        *                *
C                        ******************
C                             SIDE=3
C
C                      DATA SET UP FOR MESH-2D
C----------------------------------------------------------------------
C             NUMN   NTYPE   NSTOP   DT     AF     IAXI
C EXAMPLE     4      1       100     0.2    1.0    1
C----------------------------------------------------------------------
C             XMIN   XMAX    YMIN    YMAX
C EXAMPLE     0.0    1.00    0.0     1.00
C----------------------------------------------------------------------
C             REGION #   # OF NODES   TYPE     0,1     (TYPE=3,4,6,8)
C EXAMPLE     1          8            4        1
C----------------------------------------------------------------------
C (EA. X,Y)   REGION #   X     Y
C EXAMPLE     1          0     0
C----------------------------------------------------------------------
C (REG SIDE)  REGION #  DIRICHLET(-1)  NEUMANN(-2) OTHER
C (REGION #)  1         -1             -2             1
C----------------------------------------------------------------------
C             REGION #   # OF ROWS   # OF COLUMNS    ELEMENT NODE #S
C EXAMPLE     1          3           3              1 4 6 7 8 5 3 2
C----------------------------------------------------------------------
C             DXX        DYY         QQ
C EXAMPLE     15         15          0
C----------------------------------------------------------------------
C (EA. DIRIC) REGION #   SIDE        VALUE
C EXAMPLE     1          2           0
C----------------------------------------------------------------------
C             REGION #   SIDE        VALUE
C EXAMPLE     1          4           100
C----------------------------------------------------------------------
      CALL INFILE
      CALL GRID
      CALL BANDW
      CALL OUTPUT
      STOP
      END
C -----------------------------------------------------------------
      SUBROUTINE GRID
C    *****************************************************************
C    *                                                               *
C    *   SUBROUTINE GRID:   GENERATES MESH, BOUNDARY CONDITIONS,     *
C    *                      INITIAL CONDITIONS, AND PROPERTY DATA    *
C    *                      BASED ON LIMITED REGIONAL INPUT.         *
C    *                                                               *
C    *****************************************************************
C
C ARGUMENT LIST:
C                   LNODS(I,J)   -   CONNECTIVITY MATRIX
C                   KNODS(L,K)   -   CONNECTIVITY MATRIX FOR FLUX
C                   JNODS(L,K)   -   CONNECTIVITY MATRIX FOR CONVECTION
C                   XCORD(I)     -   X CO-ORDINATE NODAL ARRAY
C                   YCORD(I)     -   Y CO-ORDINATE NODAL ARRAY
C                    VALU(I)     -   NODAL ARRAY OF DEPENDENT VALUES
C                    LNODE       -   INTEGER NUMBER NODES PER ELEMENT
C                    IELEM       -   SPECIFIC ELEMEMT OF INTEREST
C                     NB         -   NUMBER OF DOMAIN NODE POINTS
C                     NEL        -   NUMBER OF DOMAIN ELEMENTS
C                     KBC1       -   NUMBER OF DIRICHLET BOUNDARY POINTS
C                     KBC2       -   NUMBER OF CHAUCHY BOUNDARY SURFACES
C                     KBC3       -   NUMBER OF CONVECTIVE SURFACES
C                      Z         -   ZETA LOCATION FOR INTERPOLATION
C                      E         -    ETA LOCATION FOR INTERPOLATION
C                      I         -   I=ELEMENTS
C                      J         -   J=DOMAIN NODE
C                      K         -   K=BOUNDARY NODE
C                      L         -   L=INDEX OF CHAUCHY SURFACES
C
      COMMON/A/ XP(501),YP(501),XRG(9),YRG(9),NDN(8)
      COMMON/B/ NN(20,20),YC(20,20),XC(20,20),NNRB(20,4,20),JT(20,4)
      COMMON/C/ NE(501),NR(9),ICOMP(4,4),SN(8)
      COMMON/D/ QQ(501),DX(501),DY(501)
      COMMON/E/ XCORD(501),YCORD(501)
      COMMON/F/ PHI(501),H(501),TINF(501)
      COMMON/G/ FIXED(501),NBA(501)
      COMMON/H/ LNODS(501,8),KNODS(501,3),JNODS(501,3)
      COMMON/I/ LNODE,NPOIN,NB,NEL,NELEM,NNODE,NUMN,JFLAG,IFLAG
      COMMON/J/ KBC1,KBC2,KBC3
      COMMON/P/ NSTOP,NTYPE,DT,AF,NVEL,KPRNT,IAXI
      COMMON/S/ ISI(501),LEM(501),ISIH(501),LME(501)
      COMMON/PIC/XMIN,XMAX,YMIN,YMAX
      COMMON/NEW/NORDER(501)
C
C      INITIALIZE COUNTERS.
C
      NB=0
      KBC1=0
      KBC2=0
      KBC3=0
      NEL=0
      NELOLD=1
      DO 1 I=1,501
      NORDER(I)=0
      DO 1 J=1,8
    1 LNODS(I,J)=0
      DO 2 I=1,501
      DO 2 J=1,3
      JNODS(I,J)=0
    2 KNODS(I,J)=0
      DO 3 I=1,20
      DO 3 J=1,20
    3 NN(I,J)=0
C
C READ INITIAL DATA
C
      READ(5,*)NUMN,NSTOP,NTYPE,KPRNT,NVEL,DT,AF,IAXI
      WRITE(*,94)NUMN,NSTOP,NTYPE,KPRNT,NVEL,DT,AF,IAXI
   94 FORMAT(1H1,5I5,2F10.5,I5)
      READ(5,*)XMIN,XMAX,YMIN,YMAX
      WRITE(*,96)XMIN,XMAX,YMIN,YMAX
   96 FORMAT(10X,'XMIN=',F6.2,2X,'XMAX=',F6.2,2X,'YMIN=',F6.2,2X,
     1'YMAX=',F6.2)
C
C READ #REGIONS, #INPUT PTS, ELEMENT SHAPE, PRINT FLAG
C INPUT POINTS ARE THE REGIONS CORNER AND MID POINTS
C
C      LNODE     ELEMENT GEOMETRY
C        3         3 NODE TRIANGLE
C        4         4 NODE QUADRILATERAL
C        6         6 NODE TRIANGLE
C        8         8 NODE QUADRILATERIAL
C
C IFLAG = 0 ; NO ECHO PRINTING
C IFLAG = 1 ; PRINTS COORDINATES, CONNECTIVITY MATRIX, REGIONAL INPUT
C IFLAG = 2 ; PRINTS INITIAL REGIONAL NODE NUMBERS, X,Y VALUES AND
C             CONNECTIVITY MATRIX FOR OUTPUT FILE
C
      READ(5,*) INRG,INBP,IFLAG
      LNODE=NUMN
C
C      SET INDICATOR OF NUMBER OF NODES PER ELEMENT SIDE.
C
      KNODE=2
      IF(LNODE.EQ.6.OR.LNODE.EQ.8)KNODE=3
      IF(LNODE.EQ.3.OR.LNODE.EQ.4.OR.LNODE.EQ.6.OR.LNODE.EQ.8) GOTO 15
      WRITE(*,12)
   12 FORMAT(10X,'ELEMENT SHAPE NOT PROPERLY DEFINED')
      STOP
C
C READ X AND Y COORDINATES OF INPUT POINTS
C
   15 DO 16 K=1,INBP
   16 READ(5,*)I,XP(I),YP(I)
C
C READ CONNECTIVITY MATRIX - COMMON SIDES OF ADJACENT REGIONS
C
      DO 30 I=1,INRG
   30 READ(5,*) NRG,(JT(NRG,J),J=1,4)
C
C CHECK PRINT FLAG FOR INPUT ECHO
C
      IF(IFLAG.EQ.0) GOTO 100
C
C     ECHO INPUT
C
      WRITE(*,50)
   50 FORMAT(//,1X,18HGLOBAL COORDINATES,//,1X,6HNUMBER,10X,
     1 7HX-COORD,10X,7HY-COORD)
      WRITE(*,60) (I,XP(I),YP(I),I=1,INBP)
   60 FORMAT(5X,I4,10X,F7.2,10X,F7.2)
      WRITE(*,70)
   70 FORMAT(///,1X,17HCONNECTIVITY DATA,/,1X,6HREGION,5X,1H1,
     1 5X,1H2,5X,1H3,5X,1H4)
      DO 80 I=1,INRG
   80 WRITE(*,90) I,(JT(I,J),J=1,4)
   90 FORMAT(3X,I4,2X,4(4X,I4))
C
C LOOP ON REGIONS TO GENERATE THE ELEMENTS
C
  100 DO 1000 KK=1,INRG
C
C READ REGION#, NO. ROWS OF NODES TO BE IN THE REGION, NO. COLUMNS,
C REGION BOUNDARY POINTS, STARTING IN LOWER LEFT CORNER PROGRESSING
C COUNTER CLOCKWISE
C
      READ(5,*) NRG,NROWS,NCOL,(NDN(I),I=1,8)
      READ(5,*) DXX,DYY,HGEN
C
C REGIONS MUST BE READ IN NUMERICAL ORDER
C
      IF(NRG.EQ.KK) GOTO 108
C
C PRINT MESSAGE AND TERMINATE PROGRAM IF NOT IN ORDER
C
      WRITE(*,105)
  105 FORMAT(10X,'PLACE REGIONS IN NUMERICAL ORDER - PROGRAM STOPPED')
      STOP
C
C ECHO INPUT IF DESIRED
C
  108 IF(IFLAG.NE.0)WRITE(*,110) NRG,NROWS,NCOL,(NDN(I),I=1,8)
  110 FORMAT(//,1X,6HREGION,I3,2X,4HROWS,I4,2X,7HCOLUMNS,I3,/,
     1 1X,14HBOUNDARY NODES,8I5)
C
C GENERATION OF THE ELEMENTAL NODE COORDINATES
C
      DO 120 I=1,8
      II = NDN(I)
      XRG(I) = XP(II)
  120 YRG(I) = YP(II)
      XRG(9) = XRG(1)
      YRG(9) = YRG(1)
C
C IF 6 NODED TRIANGLES OR 8 NODED QUADS ARE DESIRED CHANGE GRID
C TO INCLUDE MID POINTS.
C
      IF(LNODE.LT.5) GOTO 130
      NROWS = 2*NROWS - 1
      NCOL = 2*NCOL - 1
C
C CALCULATE GRID INCREMENTS; DELTA ETA AND DELTA SI
C
  130 DETA = 2.0 / (FLOAT(NROWS) - 1.0)
      DSI = 2.0 / (FLOAT(NCOL) - 1.0)
C
C STEP THRU LOCAL CO-ORDINATES ( ETA & SI ) AND USE SHAPE FUNCTIONS
C (SN(1) -SN(8) ) TO CALCULATE GLOBAL COORDINATES ( XC & YC ).
C START IN UPPER LEFT CORNER AND PROGRESS LEFT TO RIGHT AND TOP TO BOTTOM
C
      DO 140 I=1,NROWS
      ETA = 1.0 - DETA*(FLOAT(I) - 1.0 )
      DO 140 J=1,NCOL
      SI = -1.0 + DSI*(FLOAT(J) - 1.0)
C
C SHAPE FUNCTIONS
C
      SN(1) = -0.25 * (1.0-SI) * (1.0-ETA) * (SI+ETA+1.0)
      SN(2) = 0.50 * (1.0-SI*SI) * (1.0-ETA)
      SN(3) = 0.25 * (1.0+SI) * (1.0-ETA) * (SI-ETA-1.0)
      SN(4) = 0.50 * (1.0+SI) * (1.0-ETA*ETA)
      SN(5) = 0.25 * (1.0+SI) * (1.0+ETA) * (SI+ETA-1.0)
      SN(6) = 0.50 * (1.0-SI*SI) * (1.0+ETA)
      SN(7) = 0.25 * (1.0-SI) * (1.0+ETA) * (ETA-SI-1.0)
      SN(8) = 0.50 * (1.0-SI) * (1.0-ETA*ETA)
C
C INITIALIZE XC & YC TO ZERO PREVIOUS SUMMATION
C
      XC(I,J) = 0.0
      YC(I,J) = 0.0
C
C SUM SHAPE FUNCTIONS
C
      DO 140 K=1,8
      XC(I,J) = XC(I,J) + XRG(K)*SN(K)
  140 YC(I,J) = YC(I,J) + YRG(K)*SN(K)
C
C GENERATION OF THE REGION NODE NUMBERS
C
      KN1 = 1
      KS1 = 1
      KN2 = NROWS
      KS2 = NCOL
C
C CHECK EACH SIDE OF THE REGION FOR OUTSIDE BOUNDARY OR ANOTHER REGION
C
      DO 200 I=1,4
C
C CHECK CONNECTIVITY MATRIX
C
      NRT = JT(NRG,I)
C
C IF SIDE IS OUTSIDE BOUNDARY OR THIS REGION IS LOWER NUMBERED THAN
C ADJACENT REGION THIS REGION HAS PRIORITY IN NUMBERING THE BOUNDARY
C
      IF(NRT.LE.0.OR.NRT.GT.NRG) GOTO 200
C
C IF SIDE HAS BEEN NUMBERED PREVIOUSLY; NUMBER ACCORDING TO PREVIOUS WORK
C DETERMINE WHICH SIDE OF ADJACENT REGION IS COMMON
C
      DO 150 J=1,4
  150 IF(JT(NRT,J).EQ.NRG) NRTS = J
      K = NCOL
      IF(I.EQ.2.OR.I.EQ.4) K = NROWS
      JL = 1
      JK = ICOMP(I,NRTS)
      IF(JK.EQ.-1) JL = K
      DO 190 J=1,K
      GOTO(160,165,170,175),I
  160 NN(1,J) = NNRB(NRT,NRTS,JL)
      KN1 = 2
      GOTO 190
  165 NN(J,NCOL) = NNRB(NRT,NRTS,JL)
      KS2 = NCOL - 1
      GOTO 190
  170 NN(NROWS,J) = NNRB(NRT,NRTS,JL)
      KN2 = NROWS - 1
      GOTO 190
  175 NN(J,1) = NNRB(NRT,NRTS,JL)
      KS1 = 2
  190 JL = JL + JK
  200 CONTINUE
C
C IF ELEMENT = 8 NODED QUAD BRANCH TO DIFFERENT PART OF PROGRAM
C
      IF(LNODE.EQ.8) GOTO 800
C
C IF THERE IS ONLY ONE ROW OR ONE COLUMN ALL NODES HAVE BEEN NUMBERED,
C SKIP NUMBERING REMAINING NODES
C
      IF(KN1.GT.KN2) GOTO 270
      IF(KS1.GT.KS2) GOTO 270
C
C NUMBER NEW NODES; RELATING REGIONAL 2-D PLACEMENT WITH GLOBAL NUMBERING
C
      DO 210 I=KN1,KN2
      DO 210 J=KS1,KS2
C
C     FOR MULTIPLY CONNECTED REGIONS
C
C     IF(NRG.EQ.3)NBOW1=1
C     IF(NRG.EQ.5.AND.I.EQ.NROWS.AND.J.EQ.1)NN(I,J)=NNRB(3,2,NBOW1)
C     IF(NRG.EQ.5.AND.I.EQ.NROWS.AND.J.EQ.1)GO TO 210
C
C     IF(NRG.EQ.9)NBOW2=1
C     IF(NRG.EQ.14.AND.I.EQ.NROWS.AND.J.EQ.1)NN(I,J)=NNRB(9,2,NBOW2)
C     IF(NRG.EQ.14.AND.I.EQ.NROWS.AND.J.EQ.1)GO TO 210
C
C     END ALTERATION
C
      NB = NB + 1
      NN(I,J) = NB
  210 CONTINUE
C
C STORE BOUNDARY NODE NUMBERS
C
      DO 220 I=1,NCOL
      NNRB(NRG,3,I) = NN(NROWS,I)
  220 NNRB(NRG,1,I) = NN(1,I)
      DO 230 I=1,NROWS
      NNRB(NRG,2,I) = NN(I,NCOL)
  230 NNRB(NRG,4,I) = NN(I,1)
C
C     FLAGGED OUTPUT
C
      IF(IFLAG.EQ.0) GOTO 270
      WRITE(*,240)
  240 FORMAT(//,1X,19HREGION NODE NUMBERS)
      DO 250 I=1,NROWS
  250 WRITE(*,260) (NN(I,J),J=1,NCOL)
  260 FORMAT(16I5)
C
C NUMBER NODES CONSECUTIVELY; K - REGIONAL NO., NN - GLOBAL NO.
C
  270 K = 1
      DO 280 I=1,NROWS
      DO 280 J=1,NCOL
      NE(K) = NN(I,J)
      NEK = NE(K)
      NORDER(NEK)=NEK
      XCORD(NEK) = XC(I,J)
      YCORD(NEK) = YC(I,J)
  280 K = K + 1
      IF(LNODE.EQ.6) GOTO 600
C
C INITIALIZATION FOR 3 NODE TRIANGLES AND 4 NODE QUAD
C
      L = NROWS - 1
      DO 499 I=1,L
      DO 499 J=2,NCOL
      NR(1) = NCOL*I + J - 1
      NR(2) = NCOL*I + J
      NR(3) = NCOL*(I-1) + J
      NR(4) = NCOL*(I-1) + J -1
C
C     IF 4 NODE QUAD BRANCH
C
      IF(LNODE.EQ.4) GOTO 400
C
C     3 NODE TRIANGLE
C
C DETERMINE DIAGONAL LENGTHS
C
      DIAG1=(XC(I,J)-XC(I+1,J-1))**2+(YC(I,J)-YC(I+1,J-1))**2
      DIAG2=(XC(I+1,J)-XC(I,J-1))**2+(YC(I+1,J)-YC(I,J-1))**2
C
C DIVIDE QUADS INTO TRIANGLES
C
      DO 399 IJ=1,2
      NEL = NEL + 1
      DX(NEL)=DXX
      DY(NEL)=DYY
      QQ(NEL)=HGEN
C
C DIVIDE WITH SHORTEST DIAGONAL
C
      IF((DIAG1/DIAG2).GT.1.04) GOTO 350
      J1 = NR(1)
      J2 = NR(IJ+1)
      J3 = NR(IJ+2)
      LNODS(NEL,1) = NE(J1)
      LNODS(NEL,2) = NE(J2)
      LNODS(NEL,3) = NE(J3)
      GOTO 399
  350 J1 = NR(IJ)
      J2 = NR(IJ+1)
      J3 = NR(4)
      LNODS(NEL,1) = NE(J1)
      LNODS(NEL,2) = NE(J2)
      LNODS(NEL,3) = NE(J3)
  399 CONTINUE
      GOTO 499
C
C     4 NODE QUADS
C
  400 NEL = NEL + 1
      DX(NEL)=DXX
      DY(NEL)=DYY
      QQ(NEL)=HGEN
      J1 = NR(1)
      J2 = NR(2)
      J3 = NR(3)
      J4 = NR(4)
      LNODS(NEL,1) = NE(J1)
      LNODS(NEL,2) = NE(J2)
      LNODS(NEL,3) = NE(J3)
      LNODS(NEL,4) = NE(J4)
  499 CONTINUE
      GOTO 900
C
C     6 NODE TRIANGLES
C
  600 L = NROWS - 2
      DO 699 I=1,L,2
      DO 699 J=3,NCOL,2
C
C CALCULATE DIAGONAL LENTHS; SHORTER DIAGONAL WILL DIVIDE
C
      DIAG1=(XC(I,J)-XC(I+2,J-2))**2 +(YC(I,J)-YC(I+2,J-2))**2
      DIAG2=(XC(I+2,J)-XC(I,J-2))**2 +(YC(I+2,J)-YC(I,J-2))**2
C
C NUMBER NODES ACCORDING TO WHICH DIAGONAL WILL DIVIDE
C
      IF((DIAG1/DIAG2).GT.1.04) GOTO 610
      NR(1) = NCOL*(I+1) + J - 2
      NR(2) = NCOL*(I+1) + J - 1
      NR(3) = NCOL*(I+1) + J
      NR(4) = NCOL*I + J
      NR(5) = NCOL*I + J - 1
      NR(6) = NCOL*(I-1) + J
      NR(7) = NCOL*(I-1) + J - 1
      NR(8) = NCOL*I + J - 2
      NR(9) = NCOL*(I-1) + J - 2
      GOTO 620
  610 NR(1) = NCOL*(I-1) + J - 2
      NR(2) = NCOL*I     + J - 2
      NR(3) = NCOL*(I+1) + J - 2
      NR(4) = NCOL*(I+1) + J - 1
      NR(5) = NCOL* I    + J - 1
      NR(6) = NCOL*(I+1) + J
      NR(7) = NCOL* I    + J
      NR(8) = NCOL*(I-1) + J - 1
      NR(9) = NCOL*(I-1) + J
  620 DO 640 IJ=1,4,3
      NEL = NEL + 1
      DX(NEL)=DXX
      DY(NEL)=DYY
      QQ(NEL)=HGEN
      J1 = NR(1)
      J2 = NR(1+IJ)
      J3 = NR(2+IJ)
      J4 = NR(3+IJ)
      J5 = NR(5+IJ)
      J6 = NR(4+IJ)
      LNODS(NEL,1) = NE(J1)
      LNODS(NEL,2) = NE(J2)
      LNODS(NEL,3) = NE(J3)
      LNODS(NEL,4) = NE(J4)
      LNODS(NEL,5) = NE(J5)
  640 LNODS(NEL,6) = NE(J6)
  699 CONTINUE
      GOTO 900
C
C     8 NODED QUADS
C
C NUMBER NODES CONSECUTIVELY; GLOBAL# = NN(REGION ROW, REGION COL)
C SKIP NUMBERING IF ONLY ONE ROW OR ONE COLUMN PRESENT
C
  800 IF(KN1.GT.KN2) GOTO 816
      IF(KS1.GT.KS2) GOTO 816
      DO 814 I=1,KN2,2
      DO 812 INC=1,2
      I2 = I + INC - 1
      IF(I2.LT.KN1) GOTO 812
      IF(I2.GT.KN2) GOTO 812
      DO 810 J=1,KS2,INC
      IF(J.LT.KS1) GOTO 810
      NB = NB + 1
      NN(I2,J) = NB
  810 CONTINUE
  812 CONTINUE
  814 CONTINUE
  816 CONTINUE
C
C STORE BOUNDARY NODES
C
      DO 820 I=1,NCOL
      NNRB(NRG,3,I) = NN(NROWS,I)
  820 NNRB(NRG,1,I) = NN(1,I)
      DO 830 I=1,NROWS
      NNRB(NRG,2,I) = NN(I,NCOL)
  830 NNRB(NRG,4,I) = NN(I,1)
C
C     FLAGGED OUTPUT
C
      IF(IFLAG.EQ.0) GOTO 870
      WRITE(*,840)
  840 FORMAT(//,1X,19HREGION NODE NUMBERS)
      DO 850 I=1,NROWS
  850 WRITE(*,860) (NN(I,J),J=1,NCOL)
  860 FORMAT(16I5)
C
C NUMBER NODES CONSECUTIVELY; K=REGION#, NE(K)=GLOBAL#
C
  870 K = 1
      DO 882 I=1,NROWS,2
      DO 882 INC=1,2
      I2 = I + INC - 1
      IF(I2.GT.NROWS) GOTO 881
      DO 880 J=1,NCOL,INC
      NE(K) = NN(I2,J)
      NEK = NE(K)
      NORDER(NEK)=NEK
      XCORD(NEK) = XC(I2,J)
      YCORD(NEK) = YC(I2,J)
      K = K + 1
  880 CONTINUE
  881 CONTINUE
  882 CONTINUE
      L = ( NROWS - 1 ) / 2
      NCOL2= ( NCOL + 1 ) / 2
      DO 890 I=1,L
      DO 890 J=3,NCOL,2
      NEL = NEL + 1
      DX(NEL)=DXX
      DY(NEL)=DYY
      QQ(NEL)=HGEN
      J1 = (NCOL + NCOL2)*I + J - 2
      J2 = (NCOL + NCOL2)*I + J - 1
      J3 = (NCOL + NCOL2)*I + J
      J4 = NCOL*I + NCOL2*(I-1) + (J+1)/2
      J5 = (NCOL + NCOL2)*(I-1) + J
      J6 = (NCOL + NCOL2)*(I-1) + J - 1
      J7 = (NCOL + NCOL2)*(I-1) + J - 2
      J8 = NCOL*I + NCOL2*(I-1) + (J-1)/2
      LNODS(NEL,1) = NE(J1)
      LNODS(NEL,2) = NE(J2)
      LNODS(NEL,3) = NE(J3)
      LNODS(NEL,4) = NE(J4)
      LNODS(NEL,5) = NE(J5)
      LNODS(NEL,6) = NE(J6)
      LNODS(NEL,7) = NE(J7)
      LNODS(NEL,8) = NE(J8)
  890 CONTINUE
C
C CHECK EACH SIDE FOR BOUNDARY CONDITIONS
C
  900 DO 990 I=1,4
      J = JT(NRG,I)
      IF(J.GT.0) GOTO 990
      IF(J.EQ.-1) GOTO 910
      IF(J.EQ.-2) GOTO 960
      IF(J.EQ.-3) GOTO 1001
C
C     DIRICHLET BOUNDARY CONDITIONS
C
  910 READ(5,*) IRG,ISIDE,TEMP
      IF(IRG.NE.NRG.OR.ISIDE.NE.I) GOTO 2010
      K = NCOL
      IF(I.EQ.2.OR.I.EQ.4) K = NROWS
      DO 950 L=1,K
      KBC1 = KBC1 + 1
      GOTO(920,925,930,935),I
  920 NBA(KBC1) = NN(1,L)
      GOTO 940
  925 NBA(KBC1) = NN(L,NCOL)
      GOTO 940
  930 NBA(KBC1) = NN(NROWS,L)
      GOTO 940
  935 NBA(KBC1) = NN(L,1)
  940 FIXED(KBC1) = TEMP
C
C      OMIT ANY NODE POINTS PREVIOUSLY SPECIFIED.
C
      KM1=KBC1-1
      IF(KM1.EQ.0)GO TO 950
      DO 955 KKBC1=1,KM1
  955 IF(NBA(KBC1).EQ.NBA(KKBC1))KBC1=KM1
  950 CONTINUE
      GO TO 990
C
C     NEUMANN BOUNDARY CONDITIONS
C
  960 READ(5,*)IRG,ISIDE,PH
      IF(IRG.NE.NRG.OR.ISIDE.NE.I)GO TO 2010
      IF(LNODE.EQ.6.OR.LNODE.EQ.8)GO TO 980
      K=NCOL-1
      IF(I.EQ.2.OR.I.EQ.4)K=NROWS-1
      KSTEP=1
      GO TO 985
  980 K=NCOL-2
      IF(I.EQ.2.OR.I.EQ.4)K=NROWS-2
      KSTEP=2
  985 DO 986 IK=1,K,KSTEP
      KBC2=KBC2+1
      IF(LNODE.EQ.3.OR.LNODE.EQ.6)THEN
      IF(ISIDE.EQ.1)JSIDE=2
      IF(ISIDE.EQ.2)JSIDE=2
      IF(ISIDE.EQ.3)JSIDE=3
      IF(ISIDE.EQ.4)JSIDE=1
      ELSE IF(LNODE.EQ.4.OR.LNODE.EQ.8)THEN
      IF(ISIDE.EQ.1)JSIDE=2
      IF(ISIDE.EQ.2)JSIDE=3
      IF(ISIDE.EQ.3)JSIDE=4
      IF(ISIDE.EQ.4)JSIDE=1
      ENDIF
      NS1=0
      NS2=0
      NS3=0
      ISI(KBC2)=JSIDE
      PHI(KBC2)=PH
      KNODS(KBC2,1)=NNRB(NRG,I,IK)
      KNODS(KBC2,2)=NNRB(NRG,I,IK+1)
      IF(KNODE.EQ.3)KNODS(KBC2,3)=NNRB(NRG,I,IK+2)
      NS1=NNRB(NRG,I,IK)
      NS2=NNRB(NRG,I,IK+1)
      IF(KNODE.EQ.3)NS3=NNRB(NRG,I,IK+2)
      DO 987 IEL=NELOLD,NEL
      KE1=0
      KE2=0
      KE3=0
      DO 987 ILLN=1,LNODE
      ILL=LNODS(IEL,ILLN)
      IF(ILL.EQ.NS1)KE1=1
      IF(ILL.EQ.NS2)KE2=1
      IF(ILL.EQ.NS3)KE3=1
      KEALL=KE1+KE2+KE3
      IF(KEALL.GE.2)GOTO 988
  987 CONTINUE
  988 LEM(KBC2)=IEL
  986 CONTINUE
      GOTO 990
C
C     CONVECTION BOUNDARY CONDITIONS
C
 1001 READ(5,*)IRG,LSIDE,HH,TINFF
      IF(IRG.NE.NRG.OR.LSIDE.NE.I)GO TO 2010
      IF(LNODE.EQ.6.OR.LNODE.EQ.8)GO TO 1002
      K=NCOL-1
      IF(I.EQ.2.OR.I.EQ.4)K=NROWS-1
      KSTEP=1
      GO TO 1003
 1002 K=NCOL-2
      IF(I.EQ.2.OR.I.EQ.4)K=NROWS-2
      KSTEP=2
 1003 DO 1004 IK=1,K,KSTEP
      KBC3=KBC3+1
      IF(LNODE.EQ.3.OR.LNODE.EQ.6)THEN
      IF(LSIDE.EQ.1)JSIDE=2
      IF(LSIDE.EQ.2)JSIDE=2
      IF(LSIDE.EQ.3)JSIDE=3
      IF(LSIDE.EQ.4)JSIDE=1
      ELSE IF(LNODE.EQ.4.OR.LNODE.EQ.8)THEN
      IF(LSIDE.EQ.1)JSIDE=2
      IF(LSIDE.EQ.2)JSIDE=3
      IF(LSIDE.EQ.3)JSIDE=4
      IF(LSIDE.EQ.4)JSIDE=1
      ENDIF
      NS1=0
      NS2=0
      NS3=0
      ISIH(KBC3)=JSIDE
      H(KBC3)=HH
      TINF(KBC3)=TINFF
      JNODS(KBC3,1)=NNRB(NRG,I,IK)
      JNODS(KBC3,2)=NNRB(NRG,I,IK+1)
      IF(KNODE.EQ.3)JNODS(KBC3,3)=NNRB(NRG,I,IK+2)
      NS1=NNRB(NRG,I,IK)
      NS2=NNRB(NRG,I,IK+1)
      IF(KNODE.EQ.3)NS3=NNRB(NRG,I,IK+2)
      DO 1005 IEL=NELOLD,NEL
      KE1=0
      KE2=0
      KE3=0
      DO 1005 ILLN=1,LNODE
      ILL=LNODS(IEL,ILLN)
      IF(ILL.EQ.NS1)KE1=1
      IF(ILL.EQ.NS2)KE2=1
      IF(ILL.EQ.NS3)KE3=1
      KEALL=KE1+KE2+KE3
      IF(KEALL.GE.2)GOTO 1006
 1005 CONTINUE
 1006 LME(KBC3)=IEL
 1004 CONTINUE
C
  990 CONTINUE
      NELOLD=NEL+1
 1000 CONTINUE
C
C     OPTIMIZE MESH NUMBERING BY PASSING HORIZONTAL PLANES
C
      NSTEP=500
      IF(JFLAG.EQ.1)THEN
      WRITE(*,111)
  111 FORMAT(/,10X,'MESH BEING OPTIMIZED')
      CALL OPTNOD(NSTEP)
      ENDIF
C
C     FLAGGED OUTPUT
C
      IF(IFLAG.NE.1) GOTO 2000
      WRITE(*,1010)
 1010 FORMAT(//,'**** ELEMENTAL DATA ****')
      WRITE(*,1100)
 1100 FORMAT(///,1X,4HNODE,11X,7HX-COORD,10X,7HY-COORD,/)
      DO 1110 I=1,NB
      J=NORDER(I)
 1110 WRITE(*,1120) J,XCORD(I),YCORD(I)
 1120 FORMAT(1X,I4,10X,F9.4,8X,F9.4)
      WRITE(*,1130)
 1130 FORMAT(///,1X,10HTOPOGRAPHY,/,1X,10H----------,/,1X,
     1 7HELEMENT,15X,5HNODES,//)
      DO 1140 I=1,NEL
 1140 WRITE(*,1150) I,(NORDER(LNODS(I,J)),J=1,LNODE)
 1150 FORMAT(3X,I4,7X,8(I4))
C     WRITE(*,7000)
C7000 FORMAT (//,25X,'MATERIAL PROPERTIES ARE',//,
C    15X,'ELEMENT    THERMAL       HEAT',/,
C    25X,'NUMBER   DIFFUSIVITY   GENERATION',//)
C     DO 7080 JELEM=1,NEL
C7080 WRITE(*,7085) JELEM,QQ(JELEM),DX(JELEM),DY(JELEM)
C7085 FORMAT(5X,I5,5X,F10.5,2X,F10.5,F10.5,3X,F10.5)
      WRITE(*,1160)KBC1
 1160 FORMAT(//,1X,I6,'  DIRICHLET BOUNDARY CONDITIONS',//,1X,'NODE',
     1 5X,'VALUE',//)
      DO 1170 I=1,KBC1
 1170 WRITE(*,1180) NORDER(NBA(I)),FIXED(I)
 1180 FORMAT(1X,I4,5X,F9.4)
      WRITE(*,1190)KBC2
 1190 FORMAT(///,1X,I6,' NEUMANN BOUNDARY CONDITIONS',//,1X,
     *'FLUX          NODE NUMBERS',//)
      DO 1191 I=1,KBC2
 1191 WRITE(*,1192)PHI(I),(NORDER(KNODS(I,J)),J=1,KNODE)
 1192 FORMAT(F10.5,5X,8I6)
      WRITE(*,1193)KBC3
 1193 FORMAT(///,1X,I6,' CONVECTION BOUNDARY CONDITIONS',//,1X,
     *'H    TINF     NODE NUMBERS',//)
      DO 1195 I=1,KBC3
 1195 WRITE(*,1196)H(I),TINF(I),(NORDER(JNODS(I,J)),J=1,KNODE)
 1196 FORMAT(2F10.5,5X,8I6)
 2000 RETURN
 2010 WRITE(*,2020)
 2020 FORMAT(10X,'INCOMPATIBLE BOUNDARY CONDITION INPUT')
      STOP
      END
C -----------------------------------------------------------------
      SUBROUTINE OUTPUT
      COMMON/VH/NODES4(4,2),NODES8(4,3),NODES3(3,2),NODES6(3,3)
      COMMON/LH/NODE(501,8),KNODE(501),Q(501)
      COMMON/NEW/NORDER(501)
      COMMON/LL/X(501),Y(501)
      COMMON/D/ QQ(501),DX(501),DY(501)
      COMMON/E/ XCORD(501),YCORD(501)
      COMMON/F/ PHI(501),H(501),TINF(501)
      COMMON/G/ FIXED(501),NBA(501)
      COMMON/H/ LNODS(501,8),KNODS(501,3),JNODS(501,3)
      COMMON/I/ LNODE,NPOIN,NB,NEL,NELEM,NNODE,NUMN,JFLAG,IFLAG
      COMMON/J/ KBC1,KBC2,KBC3
      COMMON/P/ NSTOP,NTYPE,DT,AF,NVEL,KPRNT,IAXI
      COMMON/S/ ISI(501),LEM(501),ISIH(501),LME(501)
      COMMON/VE/VX(501),VY(501)
      COMMON/PIC/XMIN,XMAX,YMIN,YMAX
      CHARACTER*4 STOP(4)
      DATA STOP/'DIRC','FLUX','CONV','VELC'/
C
      NNODE=NB
      NELEM=NEL
      NUM=2
      IF(NUMN.EQ.6.OR.NUMN.EQ.8)NUM=3
      DO 33 I=1,NNODE
      NODPT=NORDER(I)
      X(NODPT)=XCORD(I)
   33 Y(NODPT)=YCORD(I)
      DO 34 J=1,NELEM
      DO 34 L=1,NUMN
   34 NODE(J,L)=NORDER(LNODS(J,L))
      DO 35 I=1,NNODE
   35 Q(I)=PHI(I)
      RHO=1.0
      CP=1.0
C
C     WRITE XMIN, XMAX, YMIN, YMAX
C
      WRITE(*,1004)XMIN,XMAX,YMIN,YMAX
C
C     WRITE THE NUMBER OF SYSTEM NODES AND ELEMENTS
C
      MTYPE=0
      NUMDIM=2
      WRITE(15,1005)MTYPE,NUMDIM,NNODE,NELEM,NUMN,NSTOP,KPRNT,NVEL
      TO=0.0
      WRITE(15,1006)NTYPE,DT,AF,TO,RHO,CP,IAXI
C
C     WRITE NODE NUMBERS AND NODAL COORDINATES
C
      WRITE(15,1007)(I,X(I),Y(I),I=1,NNODE)
C
C     WRITE SYSTEM TOPOLOGY(ELEMENT NO. AND NODE NUMBERS
C     IN COUNTER-CLOCKWISE FASHION STARTING AT ANY NODE)
C     AND INPUT SOURCE TERM FOR EACH ELEMENT
C
      DO 105 J=1,NELEM
      WRITE(15,1010)J,QQ(J),DX(J),DY(J),(NODE(J,MM),MM=1,NUMN)
  105 CONTINUE
C
C     CYCLE FOR EACH NODE HAVING SPECIFIED VALUE
C
      DO 110 I=1,KBC1
      WRITE(15,1015)NORDER(NBA(I)),FIXED(I)
  110 CONTINUE
      WRITE(15,'(6X,A4)')STOP(1)
C
C  WRITE ELEMENT SIDE THAT DEFINES BOUNDARY
C  SEGMENT WHERE FLUX (GRADIENT) IS SPECIFIED
C
      DO 130 I=1,KBC2
      WRITE(15,1025)Q(I),LEM(I),ISI(I)
  130 CONTINUE
      WRITE(15,'(6X,A4)')STOP(2)
C
C  WRITE ELEMENT SIDE THAT DEFINES BOUNDARY SEGMENT
C  WHERE CONVECTIVE SURFACE IS SPECIFIED
C
      DO 140 I=1,KBC3
      WRITE(15,1026)H(I),TINF(I),LME(I),ISIH(I)
  140 CONTINUE
      WRITE(15,'(6X,A4)')STOP(3)
C
C  PRINT ALL DATA THAT HAS BEEN READ IN
C
      WRITE(*,1035)NNODE,NELEM
      WRITE(*,1040)
      DO 150 I=1,NNODE
      WRITE(*,1045)I,X(I),Y(I)
  150 CONTINUE
      WRITE(*,1050)
      DO 155 I=1,NELEM
      WRITE(*,1055)I,QQ(I),DX(I),DY(I),(NODE(I,MM),MM=1,NUMN)
  155 CONTINUE
      WRITE(*,1060)
  555 DO 160 I=1,KBC1
      WRITE(*,1065)I,NORDER(NBA(I)),FIXED(I)
  160 CONTINUE
      WRITE(*,1080)
      DO 170 I=1,KBC2
      ISI1=ISI(I)
      L=LEM(I)
      IF(NUMN.EQ.3.OR.NUMN.EQ.6)THEN
      IF(NUM.EQ.2)WRITE(*,1085)I,Q(I),(NODE(L,NODES3(ISI1,K)),K=1,NUM)
      IF(NUM.EQ.3)WRITE(*,1085)I,Q(I),(NODE(L,NODES6(ISI1,K)),K=1,NUM)
      ELSE
      IF(NUM.EQ.2)WRITE(*,1085)I,Q(I),(NODE(L,NODES4(ISI1,K)),K=1,NUM)
      IF(NUM.EQ.3)WRITE(*,1085)I,Q(I),(NODE(L,NODES8(ISI1,K)),K=1,NUM)
      ENDIF
  170 CONTINUE
      WRITE(*,1081)
      WRITE(*,1083)
      DO 180 I=1,KBC3
      ISI1=ISIH(I)
      L=LME(I)
      IF(NUMN.EQ.3.OR.NUMN.EQ.6)THEN
      IF(NUM.EQ.2)WRITE(*,1086)I,H(I),TINF(I),(NODE(L,NODES3(ISI1,K)),
     1K=1,NUM)
      IF(NUM.EQ.3)WRITE(*,1086)I,H(I),TINF(I),(NODE(L,NODES6(ISI1,K)),
     1K=1,NUM)
      ELSE
      IF(NUM.EQ.2)WRITE(*,1086)I,H(I),TINF(I),(NODE(L,NODES4(ISI1,K)),
     1K=1,NUM)
      IF(NUM.EQ.3)WRITE(*,1086)I,H(I),TINF(I),(NODE(L,NODES8(ISI1,K)),
     1K=1,NUM)
      ENDIF
  180 CONTINUE
C
C  IF NVEL=1 READ SPECIFIED NODE NUMBERS AND NODAL VELOCITIES
C
      IF(NVEL.EQ.0)GO TO 28
      DO 26 J=1,NNODE
   26 READ(5,*)I,VX(I),VY(I)
      WRITE(*,4041)
      DO 27 I=1,NNODE
      WRITE(*,1018)I,VX(I),VY(I)
   27 WRITE(15,1018)I,VX(I),VY(I)
   28 CONTINUE
      WRITE(15,'(6X,A4)')STOP(4)
C
C     FORMAT STATEMENTS
C
 1004 FORMAT(//,5X,4F8.3)
 1005 FORMAT(8I5)
 1006 FORMAT(I4,2X,5F8.4,2X,I4)
 1007 FORMAT(5X,I5,5X,2F8.3)
 1010 FORMAT(I5,1X,3F8.3,2X,8I5)
 1015 FORMAT(10X,I5,5X,F10.5)
 1018 FORMAT(10X,I5,5X,2F10.5)
 1025 FORMAT(10X,F10.5,2I5)
 1026 FORMAT(10X,2F10.5,2I5)
 1035 FORMAT('0',10X,'NO. OF NODES=',I4,5X,'NO.OF ELEMENTS=',I4)
 1040 FORMAT('0',2X,'SUMMARY OF NODAL COORDINATES')
 1041 FORMAT('0',7X,'I',9X,'X',11X,'Y')
 1045 FORMAT(5X,I4,2(4X,F8.4))
 1050 FORMAT('0',3X,'ELEMENT',3X,'SOURCE',4X,'DX',5X,'DY',9X,'NODE NUMBE
     1RS')
 1055 FORMAT(5X,I4,3X,3F8.3,8(1X,I4))
 1060 FORMAT('0',7X,'NODES WHERE CNEW IS SPECIFIED')
 1065 FORMAT(2X,2(4X,I4),5X,F8.3)
 1080 FORMAT('0',15X,'NODES WHERE FLUX IS SPECIFIED')
 1082 FORMAT('0',8X,'I',6X,'FLUX',5X,'NODE NUMBERS')
 1081 FORMAT('0',15X,'NODES WHERE CONF IS SPECIFIED')
 1083 FORMAT('0',10X,'I',7X,'H',7X,'TINF',7X,'NODE NUMBERS')
 1085 FORMAT(6X,I4,3X,F8.3,2X,3(1X,I4))
 1086 FORMAT(8X,I4,2X,2(2X,F8.3),2X,3(1X,I4))
 4041 FORMAT('0',7X,'I',12X,'VX',12X,'VY')
      RETURN
      END
C -----------------------------------------------------------------
      BLOCK DATA
      COMMON/VH/NODES4(4,2),NODES8(4,3),NODES3(3,2),NODES6(3,3)
      COMMON/XY/X4(4),Y4(4),X8(8),Y8(8)
      COMMON/C/NE(501),NR(9),ICOMP(4,4),SN(8)
      DATA ICOMP/-1,-1,1,1,-1,-1,1,1,1,1,-1,-1,1,1,-1,-1/
      DATA NODES3/3,2,1,1,3,2/
      DATA NODES6/5,3,1,6,4,2,1,5,3/
      DATA NODES4/4,3,2,1,1,4,3,2/
      DATA NODES8/7,5,3,1,8,6,4,2,1,7,5,3/
      DATA X4/-1.0,1.0,1.0,-1.0/
      DATA Y4/-1.0,-1.0,1.0,1.0/
      DATA X8/-1.0,0.0,1.0,1.0,1.0,0.0,-1.0,-1.0/
      DATA Y8/-1.0,-1.0,-1.0,0.0,1.0,1.0,1.0,0.0/
      END
C -----------------------------------------------------------------
      SUBROUTINE BANDW
      COMMON/LH/NODE(501,8),KNODE(501),Q(501)
      COMMON/I/ LNODE,NPOIN,NB,NEL,NELEM,NNODE,NUMN,JFLAG,IFLAG
      COMMON/H/ LNODS(501,8),KNODS(501,3),JNODS(501,3)
      COMMON/NEW/NORDER(501)
      DIMENSION NOD(8)
C
      DO 34 K=1,NEL
      DO 34 I=1,NUMN
   34 NODE(K,I)=NORDER(LNODS(K,I))
C
C     CALCULATE BANDWIDTH FOR EACH ELEMENT
C
      NBW=0
      DO 1 K=1,NEL
      DO 2 I=1,NUMN
    2 NOD(I)=NODE(K,I)
      L=NUMN-1
      DO 3 I=1,L
      J=I+1
      DO 3 KK=J,NUMN
      NW=IABS(NOD(I)-NOD(KK))
      IF(NW.EQ.0)WRITE(*,100)K
  100 FORMAT(10X,'ELEMENT',I4,' HAS TWO IDENTICAL NODE NUMBERS')
      IF(NW.LE.NBW)GO TO 3
      INBW=K
      NBW=NW
    3 CONTINUE
    1 CONTINUE
      NBW=NBW+1
      WRITE(*,101) NBW,INBW
  101 FORMAT(//,5X,'BANDWIDTH IS',I4,' IN ELEMENT',I4)
      RETURN
      END
C -----------------------------------------------------------------
      SUBROUTINE OPTNOD(NSTEP)
C
C    *****************************************************************
C    *                                                               *
C    *   SUBROUTINE OPTNOD: OPTIMIZES NODE NUMBERING BY PASSING A    *
C    *                      SERIES OF HORIZONTAL LINES THROUGH THE   *
C    *                      MESH.                                    *
C    *                                                               *
C    *****************************************************************
C
      COMMON/I/ LNODE,NPOIN,NB,NEL,NELEM,NNODE,NUMN,JFLAG,IFLAG
      COMMON/E/ XCORD(501),YCORD(501)
      COMMON/PIC/XMIN,XMAX,YMIN,YMAX
      COMMON/NEW/NORDER(501)
      DIMENSION NODX(501),MORDER(501)
C
C     INITIALIZE COUNTERS.
C
      DELT=0.002
      YMAX1=YMAX+DELT
      YMIN1=YMIN-DELT
      NLAST=0
C
C     LOOP OVER THE SERIES OF HORIZONTAL LINES.
C
      STEP=(YMAX1-YMIN1)/FLOAT(NSTEP)
      DO 1000 I=1,NSTEP
      YDN=STEP*FLOAT(I-1)+YMIN1
      YUP=YDN+STEP
C
C     LOOP OVER ALL NODES AND LOCATE THOSE WITHIN GIVEN REGION.
C
      NUM=0
      DO 100 IG=1,NB
      YPT=YCORD(IG)
C     IF(I.EQ.1)PRINT *,IG,YPT
      IF(YPT.LT.YDN .OR. YPT.GE.YUP) GO TO 100
      NUM=NUM+1
      NODX(NUM)=IG
  100 CONTINUE
      IF(NUM.EQ.0) GO TO 1000
C
C     SORT THE NUM POINTS INTO ASCENDING ORDER.
C
      NUM1=NUM-1
      DO 300 NSTRT=1,NUM1
      XMIN1=XCORD(NODX(NSTRT))
      IORDER=NSTRT
      DO 200 J=NSTRT,NUM
      JNOD=NODX(J)
      XPT=XCORD(JNOD)
      IF(XPT.GT.XMIN1) GO TO 200
      XMIN1=XPT
      IORDER=J
  200 CONTINUE
      NTEMP=NODX(IORDER)
      NODX(IORDER)=NODX(NSTRT)
      NODX(NSTRT)=NTEMP
  300 CONTINUE
C
C     PLACE THIS SWEEP OF NODES INTO LARGE ARRAY.
C
      DO 400 J=1,NUM
      NORDER(J+NLAST)=NODX(J)
  400 CONTINUE
      NLAST=NLAST+NUM
C     PRINT *,NUM,NLAST,I,YDN,YUP
 1000 CONTINUE
      DO 301 I=1,NB
      DO 201 J=1,NB
      IF(I.NE.NORDER(J))GO TO 201
      MORDER(I)=J
  201 CONTINUE
  301 CONTINUE
      DO 302 I=1,NB
      NORDER(I)=MORDER(I)
  302 CONTINUE
C
C     PRINT OUT THE ARRAY FOR DEBUGGING
C
      IF(IFLAG.EQ.2)THEN
      WRITE(*,1102)
      DO 1100 I=1,NLAST
      WRITE(*,1101) I,NORDER(I)
 1100 CONTINUE
      ENDIF
 1101 FORMAT(5X,2I10)
 1102 FORMAT(//,5X,'NODE RENUMBERING DEBUG CHECK',//)
      RETURN
      END
C -----------------------------------------------------------------
      SUBROUTINE INFILE
      COMMON/I/ LNODE,NPOIN,NB,NEL,NELEM,NNODE,NUMN,JFLAG,IFLAG
      CHARACTER INFIL*12,OUTFIL*12,AJ*1
      DATA INFIL/'            '/,OUTFIL/'            '/,AJ/' '/
      WRITE(*,105)
  105 FORMAT(1X,'IF YOU WANT HARDCOPY OUTPUT, PRESS CTL-PRTSC BEFORE ENT
     1ERING INPUT FILE NAME')
      WRITE(*,'(/1X,A\)')' ENTER INPUT MESH FILE NAME: '
      READ(*,'(BN,A)')INFIL
      WRITE(*,'(/1X,A\)')' ENTER OUTPUT FILE NAME: '
      READ(*,'(BN,A)')OUTFIL
      WRITE(*,'(/1X,A\)')'DO YOU WANT TO OPTIMIZE THE MESH (Y/N)?'
      READ(*,'(BN,A)')AJ
      IF(AJ.EQ.'N'.OR.AJ.EQ.'n')JFLAG=0
      IF(AJ.EQ.'Y'.OR.AJ.EQ.'y')JFLAG=1
      OPEN(5,FILE=INFIL)
C
C     OPEN UNIT 15 FOR OUTPUT
C
      OPEN(15,FILE=OUTFIL,STATUS='NEW')
      WRITE(15,'(A)')OUTFIL
      WRITE(*,104)
  104 FORMAT(/,1X,'PLEASE WAIT, GRID IS DEVELOPING',/)
      RETURN
      END
