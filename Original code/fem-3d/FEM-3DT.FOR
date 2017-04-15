C  *********************************************************************
C  **                                                                 **
C  **                            FEM-3DT.FOR                          **
C  **                   A PROGRAM FOR SOLVING THE 3-D                 **
C  **                          TIME-DEPENDENT                         **
C  **                  TEMPERATURE CONDUCTION EQUATION                **
C  **                                                                 **
C  **                        WRITTEN BY D.W.PEPPER                    **
C  **                                                                 **
C  **       THIS PROGRAM USES A GAUSS-SEIDEL ITERATION SOLVER         **
C  **       OR GAUSSIAN ELIMINATION FOR TRIANGLES OR TETRAHEDRONS     **
C  **                                                                 **
C  **       RADIATION PROBLEMS ARE TIME DEPENDENT (AND SENSITIVE TO   **
C  **       TIME STEP ERRORS DUE TO THE DIFFERENCE OF VALUES TO THE   **
C  **       FOURTH POWER - BE CAREFUL)                                **
C  **                                                                 **
C  **       MAXIMUM NUMBER OF NODES AND/OR ELEMENTS = 81              **
C  **                                                                 **
C  *********************************************************************
      COMMON/VA/COLD(81),CNEW(81),KNODE(81),FIXED(81),E(81),AB(81)
      COMMON/VD/AF,AFM,DT,NNQS,NUMN,NUM,NNHC,NNST,NSTOP,KPRNT,IFLAG
      COMMON/VK/NNODE,NELEM,NTYPE,RHOCP,NTIME,TIME,NNRA,TO,TAMB
      CALL INIT
      DO 4 L=1,NNODE
    4 CNEW(L)=COLD(L)
      KOUNT=0
      NTIME=0
      TIME=DT
      CALL ASSEMB
      CALL BNDCON
      DO 2 NSTEP=1,NSTOP
      CALL MATRIX
      IF(KOUNT.EQ.KPRNT)THEN
        CALL PRINT
        KOUNT=0
      ENDIF
      TIME=TIME+DT
      KOUNT=KOUNT+1
      NTIME=NTIME+1
      CALL RESID
      DO 3 L=1,NNODE
    3 COLD(L)=CNEW(L)
    2 CONTINUE
C
    1 WRITE(*,10)
   10 FORMAT('0',1X,'SOLUTION IS FINISHED')
      STOP
      END
C ----------------------------------------------------------------------
      SUBROUTINE INIT
      COMMON/VA/COLD(81),CNEW(81),KNODE(81),FIXED(81),E(81),AB(81)
      COMMON/VC/Q(81),F(81),NTS(81),NQS(81),DX(81),DY(81),DZ(81)
      COMMON/VZ/X(81),Y(81),Z(81),QQ(81),AREA(81),VOLUME(81)
      COMMON/VD/AF,AFM,DT,NNQS,NUMN,NUM,NNHC,NNST,NSTOP,KPRNT,IFLAG
      COMMON/VM/H(81),TINF(81),TH(81),FORM(81)
      COMMON/VK/NNODE,NELEM,NTYPE,RHOCP,NTIME,TIME,NNRA,TO,TAMB
      COMMON/VL/NODE(81,4),NODEF(81,3),NODEC(81,3),NOR1(81,3),NOR2(81,3)
      COMMON/TI/TITLE(18)
      CHARACTER*4 STOP(4),WORD
      DATA STOP/'DIRC','FLUX','CONV','RADI'/
C
      CALL SETPC
C
      READ(5,100)TITLE
  100 FORMAT(18A4)
      WRITE(*,101)TITLE
  101 FORMAT('0',5X,18A4)
C
      READ(5,*)NNODE,NELEM,NUMN,NSTOP,KPRNT
      READ(5,*)NTYPE,DT,AF,RHOCP,TO,TAMB
      IF(NUMN.EQ.3)THEN
        NUM=2
         ELSE
        NUM=3
      ENDIF
      DO 33 I=1,NNODE
      CNEW(I)=0.
      COLD(I)=TO
      NTS(I)=0
      NQS(I)=0
      X(I)=0.
      Y(I)=0.
      Z(I)=0.
      QQ(I)=0.
      Q(I)=0.
      F(I)=0.
      FIXED(I)=0.
      KNODE(I)=0
      H(I)=0.
      TINF(I)=0.
      DX(I)=0.
      DY(I)=0.
      DZ(I)=0.
      DO 34 K=1,4
   34 NODE(I,K)=0
   33 CONTINUE
C
      DO 104 J=1,NNODE
      READ(5,*)I,X(I),Y(I),Z(I)
  104 CONTINUE
C
      DO 105 I=1,NELEM
      READ(5,*)J,QQ(J),DX(J),DY(J),DZ(J),TH(J),(NODE(J,MM),MM=1,NUMN)
  105 CONTINUE
C
      DO 110 I=1,NNODE
      READ(5,1015)WORD,NT,CNT
      IF (WORD.EQ.STOP(1))GOTO 120
      NTS(I)=NT
      COLD(NT)=CNT
      KNODE(NT)=1
      FIXED(NT)=CNT
  110 CONTINUE
  120 NNST=I-1
C
      DO 130 I=1,NELEM
      READ(5,1016)WORD,Q(I),(NODEF(I,J),J=1,NUM)
      IF(WORD.EQ.STOP(2))GO TO 135
  130 CONTINUE
  135 NNQS=I-1
C
      DO 140 I=1,NELEM
      READ(5,1017)WORD,H(I),TINF(I),(NODEC(I,J),J=1,NUM)
      IF(WORD.EQ.STOP(3))GO TO 145
  140 CONTINUE
  145 NNHC=I-1
C
      DO 147 I=1,NELEM
      READ(5,1018)WORD,E(I),AB(I),FORM(I),(NOR1(I,J),J=1,NUM)
     *            ,(NOR2(I,J),J=1,NUM)
      IF(WORD.EQ.STOP(4))GO TO 148
  147 CONTINUE
  148 NNRA=I-1
C
      IF(NTYPE.EQ.1)WRITE(*,14)NTYPE
      IF(NTYPE.EQ.2)WRITE(*,15)NTYPE
      AFM=AF-1.0
      IF(NTYPE.EQ.1)GO TO 1
      IF(AF.EQ.1.)WRITE(*,10)DT
      IF(AF.EQ.0.5)WRITE(*,11)DT
    1 WRITE(*,12)NSTOP,KPRNT,TO,TAMB
      WRITE(*,1035)NNODE,NELEM
      WRITE(*,1040)
      WRITE(*,1041)
      DO 150 I=1,NNODE
      WRITE(*,1045)I,X(I),Y(I),Z(I)
  150 CONTINUE
      WRITE(*,1050)
      DO 155 I=1,NELEM
      WRITE(*,1055)I,QQ(I),DX(I),DY(I),DZ(I),(NODE(I,MM),MM=1,NUMN)
  155 CONTINUE
      WRITE(*,1060)
      WRITE(*,1061)
      DO 160 I=1,NNST
      WRITE(*,1065)I,NTS(I),COLD(NTS(I))
  160 CONTINUE
      WRITE(*,1080)
      WRITE(*,1082)
      DO 170 I=1,NNQS
      WRITE(*,1085)I,Q(I),(NODEF(I,J),J=1,NUM)
  170 CONTINUE
      WRITE(*,1081)
      WRITE(*,1083)
      DO 180 I=1,NNHC
      WRITE(*,1086)I,H(I),TINF(I),(NODEC(I,J),J=1,NUM)
  180 CONTINUE
      WRITE(*,1087)
      WRITE(*,1088)
      DO 183 I=1,NNRA
      IF(NUM.EQ.2)THEN
      WRITE(*,1089)I,E(I),AB(I),FORM(I),(NOR1(I,J),J=1,NUM),
     *             (NOR2(I,J),J=1,NUM)
      ELSE
      WRITE(*,1090)I,E(I),AB(I),FORM(I),(NOR1(I,J),J=1,NUM),
     *             (NOR2(I,J),J=1,NUM)
      ENDIF
  183 CONTINUE
C
      WRITE(*,191)
  191 FORMAT('0',5X,'PLEASE WAIT, SOLUTION IS PROCEEDING')
C
   10 FORMAT('0',10X,'FULLY IMPLICIT ',2X,'DT=',F7.4)
   11 FORMAT('0',10X,'CRANK NICOLSON ',2X,'DT=',F7.4)
   12 FORMAT('0',5X,'NSTOP=',I5,2X,'KPRNT=',I4,2X,'TO=',F8.3,2X,
     *'T AMBIENT=',F8.3)
   14 FORMAT('0',10X,'NTYPE=',I2,2X,'STEADY STATE CNEW')
   15 FORMAT('0',10X,'NTYPE=',I2,2X,'TIME DEPENDENT CNEW')
 1015 FORMAT(6X,A4,I5,5X,F10.5)
 1016 FORMAT(6X,A4,F10.5,3I5)
 1017 FORMAT(6X,A4,2F10.5,3I5)
 1018 FORMAT(6X,A4,3F10.5,6I3)
 1035 FORMAT('0',10X,'NO. OF NODES=',I4,5X,'NO.OF ELEMENTS=',I4)
 1040 FORMAT('0',8X,'NODAL COORDINATE VALUES')
 1041 FORMAT('0',7X,'I',9X,'X',11X,'Y',11X,'Z')
 1045 FORMAT(5X,I4,3(4X,F8.4))
 1050 FORMAT('0',1X,'ELEMENT',3X,'SOURCE',4X,'DX',6X,'DY',6X,'DZ',6X,'NO
     $DE NUMBERS')
 1055 FORMAT(2X,I4,3X,4F8.3,3X,4(1X,I3))
 1060 FORMAT('0',7X,'NODES WHERE CNEW IS SPECIFIED')
 1061 FORMAT('0',10X,'I',5X,'NODE',6X,'CNEW')
 1065 FORMAT(4X,2(4X,I4),5X,F8.3)
 1080 FORMAT('0',7X,'NODES WHERE FLUX IS SPECIFIED')
 1081 FORMAT('0',7X,'NODES WHERE CONV IS SPECIFIED')
 1082 FORMAT('0',10X,'I',6X,'FLUX',5X,'NODE NUMBERS')
 1083 FORMAT('0',10X,'I',7X,'H',7X,'TINF',7X,'NODE NUMBERS')
 1085 FORMAT(8X,I4,2X,F8.3,4X,3I5)
 1086 FORMAT(8X,I4,2X,2(F8.3,2X),3X,3I5)
 1087 FORMAT('0',7X,'NODES WHERE RADIATION IS SPECIFIED')
 1088 FORMAT('0',10X,'I',7X,'E',7X,'AB',7X,'FORM',7X,'ELE1 NODES
     *',5X,'ELE2 NODES')
 1089 FORMAT(8X,I4,2X,2F8.3,2X,F8.3,5X,2I5,5X,2I5)
 1090 FORMAT(8X,I4,2X,2F8.3,2X,F8.3,5X,3I5,2X,3I5)
      RETURN
      END
C ----------------------------------------------------------------------
      SUBROUTINE ASSEMB
      COMMON/VA/COLD(81),CNEW(81),KNODE(81),FIXED(81),E(81),AB(81)
      COMMON/VV/A(81,81),P(81,81),R(81,81),B(81)
      COMMON/VZ/X(81),Y(81),Z(81),QQ(81),AREA(81),VOLUME(81)
      COMMON/VK/NNODE,NELEM,NTYPE,RHOCP,NTIME,TIME,NNRA,TO,TAMB
      COMMON/VL/NODE(81,4),NODEF(81,3),NODEC(81,3),NOR1(81,3),NOR2(81,3)
      COMMON/VC/Q(81),F(81),NTS(81),NQS(81),DX(81),DY(81),DZ(81)
      COMMON/VD/AF,AFM,DT,NNQS,NUMN,NUM,NNHC,NNST,NSTOP,KPRNT,IFLAG
      COMMON/VM/H(81),TINF(81),TH(81),FORM(81)
C
      DO 1 I=1,NNODE
      DO 1 J=1,NNODE
      P(I,J)=0.0
      R(I,J)=0.0
    1 A(I,J)=0.0
C
      DO 2 K=1,NELEM
C
C     CONDUCTIVITY AND MASS MATRICES
C
      I = NODE(K,1)
      J = NODE(K,2)
      N = NODE(K,3)
      XI = X(I)
      XJ = X(J)
      XK = X(N)
      YI = Y(I)
      YJ = Y(J)
      YK = Y(N)
      IF (NUMN .EQ. 4) GO TO 99
      AREA(K) = ABS(XJ*YK - XK*YJ - XI*(YK-YJ) + YI*(XK-XJ))/2.
C
C     INTERNAL HEAT GENERATION
C
      FACTH =  QQ(K) * AREA(K) * TH(K) /3.
      F(I) =  F(I)+FACTH
      F(J) =  F(J)+FACTH
      F(N) =  F(N)+FACTH
C
      BI = YJ-YK
      BJ = YK-YI
      BK = YI-YJ
      CI = XK-XJ
      CJ = XI-XK
      CK = XJ-XI
      CX = DX(K)*TH(K)/(4.*AREA(K))
      CY = DY(K)*TH(K)/(4.*AREA(K))
      A(I,I)=A(I,I)+ BI*CX*BI + CI*CY*CI
      A(I,J)=A(I,J)+ BJ*CX*BI + CJ*CY*CI
      A(I,N)=A(I,N)+ BK*CX*BI + CK*CY*CI
      A(J,J)=A(J,J)+ BJ*CX*BJ + CJ*CY*CJ
      A(J,N)=A(J,N)+ BK*CX*BJ + CK*CY*CJ
      A(N,N)=A(N,N)+ BK*CX*BK + CK*CY*CK
      A(J,I)=A(I,J)
      A(N,I)=A(I,N)
      A(N,J)=A(J,N)
C
C     MASS MATRIX - 2-D TRIANGLE
C
      IF(NTYPE.EQ.1)GOTO 2
      FACTOR = RHOCP*AREA(K)*TH(K) / 12.
      FACTO2 = FACTOR + FACTOR
      P(I,I) = P(I,I)+FACTO2
      P(I,J) = P(I,J)+FACTOR
      P(I,N) = P(I,N)+FACTOR
      P(J,J) = P(J,J)+FACTO2
      P(J,N) = P(J,N)+FACTOR
      P(N,N) = P(N,N)+FACTO2
      P(J,I) = P(I,J)
      P(N,I) = P(I,N)
      P(N,J) = P(J,N)
      GOTO 2
C
C     CONDUCTIVITY MATRIX - 3-D TETRAHEDRON
C
   99 M=NODE(K,4)
      XL = X(M)
      YL = Y(M)
      ZI = Z(I)
      ZJ = Z(J)
      ZK = Z(N)
      ZL = Z(M)
      DET11 = XJ*(YK*ZL-YL*ZK) - YJ*(XK*ZL-XL*ZK) + ZJ*(XK*YL-XL*YK)
      DET21 = -YK*ZL + YL*ZK + YJ*(ZL - ZK) - ZJ*(YL - YK)
      DET22=-(-YL*ZI + YI*ZL + YK*(ZI-ZL) - ZK*(YI-YL))
      DET23 = -YI*ZJ + YJ*ZI + YL*(ZJ - ZI) - ZL*(YJ - YI)
      DET24=-(-YJ*ZK + YK*ZJ + YI*(ZK - ZJ) - ZI*(YK - YJ))
      DET31=-(XJ*(ZL - ZK) -(XK*ZL - XL*ZK) + ZJ*(XK - XL))
      DET32=  XK*(ZI - ZL) - (XL*ZI - XI*ZL) + ZK*(XL - XI)
      DET33=-(XL*(ZJ - ZI) - (XI*ZJ - XJ*ZI) + ZL*(XI - XJ))
      DET34=  XI*(ZK - ZJ) - (XJ*ZK - XK*ZJ) +ZI*(XJ - XK)
      DET41 = -XJ*(YK-YL) + YJ*(XK - XL) - XK*YL + XL*YK
      DET42=-(-XK*(YL - YI) + YK*(XL - XI) - XL*YI + XI*YL)
      DET43 = -XL*(YI -YJ) + YL*(XI - XJ) - XI*YJ + XJ*YI
      DET44=-(-XI*(YJ - YK) + YI*(XJ - XK) - XJ*YK + XK*YJ )
C
      VOLUME(K)=ABS((DET11+DET21*XI+DET31*YI+DET41*ZI)/6.)
C
C     INTERNAL HEAT GENERATION
C
      FACTH =  QQ(K) * VOLUME(K)*.25
      F(I) =  F(I)+FACTH
      F(J) =  F(J)+FACTH
      F(N) =  F(N)+FACTH
      F(M) =  F(M)+FACTH
C
      CW = 36.*VOLUME(K)
      CX = DX(K)/CW
      CY = DY(K)/CW
      CZ = DZ(K)/CW
      A(I,I)=A(I,I)+ DET21*CX*DET21 + DET31*CY*DET31
     *     + DET41*CZ*DET41
      A(I,J)=A(I,J)+ DET22*CX*DET21 + DET32*CY*DET31
     *     + DET42*CZ*DET41
      A(I,N)=A(I,N)+ DET23*CX*DET21 + DET33*CY*DET31
     *     + DET43*CZ*DET41
      A(I,M)=A(I,M)+ DET21*CX*DET24 + DET31*CY*DET34
     *     + DET41*CZ*DET44
      A(J,J)=A(J,J)+ DET22*CX*DET22 + DET32*CY*DET32
     *     + DET42*CZ*DET42
      A(J,N)=A(J,N)+ DET23*CX*DET22 + DET33*CY*DET32
     *     + DET43*CZ*DET42
      A(J,M)=A(J,M)+ DET22*CX*DET24 + DET32*CY*DET34
     *     + DET42*CZ*DET44
      A(N,N)=A(N,N)+ DET23*CX*DET23 + DET33*CY*DET33
     *     + DET43*CZ*DET43
      A(N,M)=A(N,M)+ DET23*CX*DET24 + DET33*CY*DET34
     *     + DET43*CZ*DET44
      A(M,M)=A(M,M)+ DET24*CX*DET24 + DET34*CY*DET34
     *     + DET44*CZ*DET44
      A(J,I)=A(I,J)
      A(N,I)=A(I,N)
      A(M,I)=A(I,M)
      A(N,J)=A(J,N)
      A(M,J)=A(J,M)
      A(M,N)=A(N,M)
C
C     MASS MATRIX 3-D TETRAHEDRON
C
      IF(NTYPE.EQ.1)GOTO 2
      FACTOR = RHOCP*VOLUME(K)*.05
      FACTO2 = FACTOR + FACTOR
      P(I,I) = P(I,I) + FACTO2
      P(I,J) = P(I,J) + FACTOR
      P(I,N) = P(I,N) + FACTOR
      P(I,M) = P(I,M) + FACTOR
      P(J,J) = P(J,J) + FACTO2
      P(J,N) = P(J,N) + FACTOR
      P(J,M) = P(J,M) + FACTOR
      P(N,N) = P(N,N) + FACTO2
      P(N,M) = P(N,M) + FACTOR
      P(M,M) = P(M,M) + FACTO2
      P(J,I) = P(I,J)
      P(N,I) = P(I,N)
      P(N,J) = P(J,N)
      P(M,I) = P(I,M)
      P(M,J) = P(J,M)
      P(M,N) = P(N,M)
    2 CONTINUE
C
C     THESE WRITE STATEMENTS ALLOW YOU TO EXAMINE THE P AND A MATRICES
C     NOTE THAT THE GLOBAL MATRIX IS SET TO 12X12 IN THIS CASE (12F6.2)
C
C     write(*,100)
C     do 33 i=1,nnode
C  33 write(*,34)(p(i,j),j=1,nnode)
C  34 FORMat(2x,12f6.2)
C 100 FORMat(/,2x,'values of p')
C
C     write(*,103)
C     do 35 i=1,nnode
C  35 write(*,34)(a(i,j),j=1,nnode)
C 103 FORMat(/,2x,'values of a before convection')
C
C     CONVECTIVE HEAT FLUX INFLUENCE
C
      IF(NNHC.EQ.0)GOTO 6
      DO 4 K=1,NNHC
C
C     CONVECTION ACROSS BOUNDARY
C
      I = NODEC(K,1)
      J = NODEC(K,2)
      IF (NUMN .EQ. 4) GO TO 300
      H1 = SQRT((X(I)-X(J))**2+(Y(I)-Y(J))**2)*TH(K)*H(K)/6.
      H2 = 2. * H1
      A(I,I) = A(I,I)+H2
      A(I,J) = A(I,J)+H1
      A(J,J) = A(J,J)+H2
      A(J,I) = A(I,J)
      GOTO 4
C
C     CONVECTION OVER SURFACE
C
  300 N = NODEC(K,3)
      CALL SIDE(ASIDE,I,J,N)
      H1 = H(K) * ASIDE / 12.
      H2 = H1 * 2.
      A(I,I) = A(I,I)+H2
      A(I,J) = A(I,J)+H1
      A(I,N) = A(I,N)+H1
      A(J,J) = A(J,J)+H2
      A(J,N) = A(J,N)+H1
      A(N,N) = A(N,N)+H2
      A(J,I) = A(I,J)
      A(N,I) = A(I,N)
      A(N,J) = A(J,N)
    4 CONTINUE
C
C     THIS WRITE STATEMENT ALLOWS YOU TO EXAMINE THE A MATRIX AFTER THE
C     INFLUENCE OF CONVECTION.
C
C     write(*,101)
C     do 36 i=1,nnode
C  36 write(*,34)(a(i,j),j=1,nnode)
C 101 FORMat(/,2x,'values of a after convection')
C
C     CALCULATE LHS
C
    6 CONTINUE
      IF(IFLAG.EQ.2)THEN
      DO 5 L=1,NNODE
      DO 5 LL=1,NNODE
      R(L,LL)=R(L,LL)+AF*A(L,LL)+P(L,LL)/DT
    5 CONTINUE
      ENDIF
C
C     THIS WRITE STATEMENT ALLOWS YOU TO EXAMINE THE G GLOBAL MATRIX
C
C     write(*,102)
C     do 3 i=1,nnode
C   3 write(*,34)(r(i,j),j=1,nnode)
C 102 FORMat(/,2x,'values of r')
C
      RETURN
      END
C ----------------------------------------------------------------------
      SUBROUTINE MATRIX
      COMMON/VA/COLD(81),CNEW(81),KNODE(81),FIXED(81),E(81),AB(81)
      COMMON/VV/A(81,81),P(81,81),R(81,81),B(81)
      COMMON/VZ/X(81),Y(81),Z(81),QQ(81),AREA(81),VOLUME(81)
      COMMON/VK/NNODE,NELEM,NTYPE,RHOCP,NTIME,TIME,NNRA,TO,TAMB
      COMMON/VC/Q(81),F(81),NTS(81),NQS(81),DX(81),DY(81),DZ(81)
      COMMON/VD/AF,AFM,DT,NNQS,NUMN,NUM,NNHC,NNST,NSTOP,KPRNT,IFLAG
      COMMON/VM/H(81),TINF(81),TH(81),FORM(81)
C
      CALL CLEAR(B,NNODE)
C
C     CALCULATE RADIATION EFFECT
C
      IF(NNRA.NE.0)CALL RADIAN
C
C     CALCULATE LOAD VECTOR (RHS PLUS VALUES AT TIME LEVEL N)
C
    4 DO 2 L=1,NNODE
      B(L)=F(L)
      DO 2 LL=1,NNODE
      B(L)=B(L)+(AFM*A(L,LL)+P(L,LL)/DT)*COLD(LL)
    2 CONTINUE
C
C     THIS WRITE STATEMENT ALLOWS YOU TO EXAMINE THE B COLUMN VECTOR
C
C     write(*,102),ntime
C     do 3 i=1,nnode
C   3 write(*,34)b(i)
C 102 FORMat(/,2x,'values of b at n=',i4)
C  34 FORMat(2x,f8.4)
C
      IF(IFLAG.EQ.1)THEN
       CALL SEIDEL
        ELSE
          DO 3 K=1,NNST
         L=NTS(K)
    3   B(L)=0.0
       CALL GAUSSR(CNEW,NNODE)
      ENDIF
      RETURN
      END
C ----------------------------------------------------------------------
      SUBROUTINE BNDCON
      COMMON/VA/COLD(81),CNEW(81),KNODE(81),FIXED(81),E(81),AB(81)
      COMMON/VV/A(81,81),P(81,81),R(81,81),B(81)
      COMMON/VZ/X(81),Y(81),Z(81),QQ(81),AREA(81),VOLUME(81)
      COMMON/VK/NNODE,NELEM,NTYPE,RHOCP,NTIME,TIME,NNRA,TO,TAMB
      COMMON/VL/NODE(81,4),NODEF(81,3),NODEC(81,3),NOR1(81,3),NOR2(81,3)
      COMMON/VC/Q(81),F(81),NTS(81),NQS(81),DX(81),DY(81),DZ(81)
      COMMON/VD/AF,AFM,DT,NNQS,NUMN,NUM,NNHC,NNST,NSTOP,KPRNT,IFLAG
      COMMON/VM/H(81),TINF(81),TH(81),FORM(81)
C
      IF(NNQS.EQ.0)GOTO 200
      DO 1 K=1,NNQS
C
C     FLUX NORMAL TO SURFACE
C
      I = NODEF(K,1)
      J = NODEF(K,2)
      IF (NUMN.EQ.4) GO TO 100
      D= SQRT((X(I)-X(J))**2 + (Y(I)-Y(J))**2) * TH(K)
      FACTOR = Q(K) * D * .5
      F(I) = F(I) + FACTOR
      F(J) = F(J) + FACTOR
      GOTO 1
  100 N = NODEF(K,3)
      CALL SIDE(ASIDE,I,J,N)
      FACTOR = Q(K) * ASIDE/3.
      F(I) = F(I) + FACTOR
      F(J) = F(J) + FACTOR
      F(N) = F(N) + FACTOR
    1 CONTINUE
C
C     CONVECTION BOUNDARY CONDITIONS
C
  200 IF(NNHC.EQ.0)RETURN
      DO 2 K=1,NNHC
C
C     CONVECTION ACROSS BOUNDARY
C
      I = NODEC(K,1)
      J = NODEC(K,2)
      IF(NUMN.EQ.4)GOTO 300
      H1 = SQRT((X(I)-X(J))**2+(Y(I)-Y(J))**2)*TH(K)*H(K)/6.
      H4 =  H1*TINF(K)*3.
      F(I) = F(I) + H4
      F(J) = F(J) + H4
      GOTO 2
C
C     CONVECTION OVER SURFACE
C
  300 N = NODEC(K,3)
      CALL SIDE(ASIDE,I,J,N)
      H1 = H(K) * ASIDE / 12.
      H4 = H1 * 4. * TINF(K)
      F(I) = F(I) + H4
      F(J) = F(J) + H4
      F(N) = F(N) + H4
    2 CONTINUE
      RETURN
      END
C ----------------------------------------------------------------------
      SUBROUTINE SIDE(ASIDE,I,J,N)
      COMMON/VZ/X(81),Y(81),Z(81),QQ(81),AREA(81),VOLUME(81)
      XI = X(I)
      XJ = X(J)
      XK = X(N)
      YI = Y(I)
      YJ = Y(J)
      YK = Y(N)
      ZI = Z(I)
      ZJ = Z(J)
      ZK = Z(N)
      DET1=-(-YJ*ZK + YK*ZJ + YI*(ZK - ZJ) - ZI*(YK - YJ))
      DET2=  XI*(ZK - ZJ) - (XJ*ZK - XK*ZJ) +ZI*(XJ - XK)
      DET3=-(-XI*(YJ - YK) + YI*(XJ - XK) - XJ*YK + XK*YJ )
      ASIDE = SQRT(.25*(DET1**2+DET2**2+DET3**2))
      RETURN
      END
C ----------------------------------------------------------------------
      SUBROUTINE CLEAR(A,N)
      DIMENSION A(1)
      DO 1 I=1,N
    1 A(I)=0.0
      RETURN
      END
C ----------------------------------------------------------------------
      SUBROUTINE PRINT
      COMMON/VA/COLD(81),CNEW(81),KNODE(81),FIXED(81),E(81),AB(81)
      COMMON/VK/NNODE,NELEM,NTYPE,RHOCP,NTIME,TIME,NNRA,TO,TAMB
      IF(NTYPE.EQ.1)THEN
      WRITE(*,403)
      ELSE
      WRITE(*,402)NTIME,TIME
      ENDIF
      WRITE(*,400)
      DO 1 I=1,NNODE
    1 WRITE(*,401)I,CNEW(I)
  401 FORMAT(5X,I3,5X,4(F8.3))
  402 FORMAT('0',5X,'TIME STEPS=',I3,2X,'TIME=',F8.3)
  403 FORMAT('0',5X,'STEADY STATE')
  400 FORMAT('0',5X,'NODE',5X,'CNEW')
      RETURN
      END
C ----------------------------------------------------------------------
      SUBROUTINE RESID
      COMMON/VK/NNODE,NELEM,NTYPE,RHOCP,NTIME,TIME,NNRA,TO,TAMB
      COMMON/VD/AF,AFM,DT,NNQS,NUMN,NUM,NNHC,NNST,NSTOP,KPRNT,IFLAG
      COMMON/VA/COLD(81),CNEW(81),KNODE(81),FIXED(81),E(81),AB(81)
      INTEGER LOC,ERRPRN
      REAL MAXERR
      DATA ERRPRN/10/
C
      MAXERR = 0.0
      MAXRES=0
      ERR=0.0001
      DO 1 I=1,NNODE
      RR=ABS(CNEW(I)-COLD(I))
      MAXERR=MAX(RR,MAXERR)
      IF(RR.GT.ERR)THEN
          MAXRES = 1
          LOC    = I
      ENDIF
    1 CONTINUE
C
      IF(MOD(NTIME,ERRPRN).EQ.0)THEN
         WRITE(*,101)NTIME,MAXERR,LOC
      ENDIF
      IF(MAXRES.EQ.1)THEN
      IF(NTIME.EQ.NSTOP)THEN
        WRITE(*,102)NTIME,RR
         CALL PRINT
          CALL PRINTF
         STOP
        ENDIF
      RETURN
       ELSE
        WRITE(*,100)NTIME
         CALL PRINT
          CALL PRINTF
         STOP
        END IF
  100 FORMAT('0',2X,'PROGRAM HAS CONVERGED IN',I3,' STEPS')
  101 FORMAT('0',2X,'TIME=',I4,2X,'ERROR=',1PE12.5,2X,'NODE=',I4)
  102 FORMAT('0',2X,'PROGRAM HAS NOT CONVERGED IN ',I3,' STEPS',/,3X,'TH
     *E RESIDUAL VALUE IS ',F8.4)
      RETURN
      END
C ----------------------------------------------------------------------
      SUBROUTINE PRINTF
      COMMON/VK/NNODE,NELEM,NTYPE,RHOCP,NTIME,TIME,NNRA,TO,TAMB
      COMMON/VA/COLD(81),CNEW(81),KNODE(81),FIXED(81),E(81),AB(81)
      COMMON/VZ/X(81),Y(81),Z(81),QQ(81),AREA(81),VOLUME(81)
      COMMON/TI/TITLE(18)
      WRITE(8,'(18A)')TITLE
      WRITE(8,'(I4)')NTIME
      WRITE(8,'(2I4)')NNODE,NELEM
      DO 1 I=1,NNODE
    1 WRITE(8,100)I,X(I),Y(I),Z(I),CNEW(I)
  100 FORMAT(5X,I3,7(F8.3,1X))
      RETURN
      END
C ----------------------------------------------------------------------
      SUBROUTINE SETPC
      COMMON/VD/AF,AFM,DT,NNQS,NUMN,NUM,NNHC,NNST,NSTOP,KPRNT,IFLAG
      CHARACTER INFIL*12,OUTFIL*12
      DATA INFIL/'            '/,OUTFIL/'            '/
      WRITE(*,105)
  105 FORMAT(1X,'IF YOU WANT HARDCOPY OUTPUT, PRESS CTL-PRTSC NOW')
      WRITE(*,'(/1X,A\)')' ENTER NAME FOR INPUT FILE: '
      READ(*,'(BN,A)')INFIL
      OPEN(5,FILE=INFIL)
C
      WRITE(*,'(/1X,A\)')' ENTER NAME FOR PLOT OUTPUT FILE: '
      READ(*,'(BN,A)')OUTFIL
      OPEN(8,FILE=OUTFIL,STATUS='NEW')
      WRITE(8,'(A)')OUTFIL
C
      WRITE(*,'(/1X,A\)')' ENTER 1 FOR GAUSS-SEIDEL ITERATION OR 2 FOR G
     1AUSSIAN ELIMINATION: '
      READ(*,'(I1)')IFLAG
C
      RETURN
      END
C ----------------------------------------------------------------------
      SUBROUTINE RADIAN
      COMMON/VZ/X(81),Y(81),Z(81),QQ(81),AREA(81),VOLUME(81)
      COMMON/VA/COLD(81),CNEW(81),KNODE(81),FIXED(81),E(81),AB(81)
      COMMON/VD/AF,AFM,DT,NNQS,NUMN,NUM,NNHC,NNST,NSTOP,KPRNT,IFLAG
      COMMON/VK/NNODE,NELEM,NTYPE,RHOCP,NTIME,TIME,NNRA,TO,TAMB
      COMMON/VL/NODE(81,4),NODEF(81,3),NODEC(81,3),NOR1(81,3),NOR2(81,3)
      COMMON/VC/Q(81),F(81),NTS(81),NQS(81),DX(81),DY(81),DZ(81)
      COMMON/VM/H(81),TINF(81),TH(81),FORM(81)
C
C     RADIATION NORMAL TO BOUNDARY
C
C     UNITS MUST BE IN SI (METERS,WATTS,JOULES,AND K TEMPERATURE)
C
      SIG=5.669E-08
      DO 1 K=1,NNRA
      FNM=FORM(K)
      I = NOR1(K,1)
      J = NOR1(K,2)
      IF(NUMN.EQ.4)GOTO 100
      XL = SQRT((X(I)-X(J))**2 + (Y(I)-Y(J))**2) * TH(K)
      TAVEG = 0.5*(COLD(I) + COLD(J))
      TV1=TAVEG**4*E(K)
      IF(NOR2(K,1).EQ.0)THEN
      TV2=TAMB**4*AB(K)*FNM
      QR = (TV1-TV2)*XL*SIG*0.5
      F(I)=F(I)-QR
      F(J)=F(J)-QR
      ELSE
      L=NOR2(K,1)
      LL=NOR2(K,2)
      T1 = 0.5*(COLD(L) + COLD(LL))
      TV2=T1**4*AB(K)*FNM
      QR = (TV1-TV2)*XL*SIG*0.5
      F(I)=F(I)-QR
      F(J)=F(J)-QR
      ENDIF
      GOTO 1
C
  100 CONTINUE
      N= NOR1(K,3 )
      CALL SIDE(ASIDE,I,J,N)
      A1 = ASIDE
      TAVEG = (COLD(I) + COLD(J) + COLD(N)) / 3.
      TV1=TAVEG**4*E(K)
      IF(NOR2(K,1).NE.0)GOTO 220
      TV2=TAMB**4*AB(K)*FNM
      QR = (TV1-TV2)*A1*SIG/3.
      F(I)=F(I)-QR
      F(J)=F(J)-QR
      F(N)=F(N)-QR
      GOTO 1
  220 II = NOR2(K,1)
      JJ = NOR2(K,2)
      NN = NOR2(K,3)
      T1 = (COLD(II) + COLD(JJ) + COLD(NN)) / 3.
      TV2=T1**4*AB(K)*FNM
      QR = (TV1-TV2)*A1*SIG/3.
      F(I)=F(I)-QR
      F(J)=F(J)-QR
      F(N)=F(N)-QR
    1 CONTINUE
C
      RETURN
      END
C ----------------------------------------------------------------------
      SUBROUTINE SEIDEL
      COMMON/VA/COLD(81),CNEW(81),KNODE(81),FIXED(81),E(81),AB(81)
      COMMON/VV/A(81,81),P(81,81),R(81,81),B(81)
      COMMON/VK/NNODE,NELEM,NTYPE,RHOCP,NTIME,TIME,NNRA,TO,TAMB
      COMMON/VC/Q(81),F(81),NTS(81),NQS(81),DX(81),DY(81),DZ(81)
      COMMON/VD/AF,AFM,DT,NNQS,NUMN,NUM,NNHC,NNST,NSTOP,KPRNT,IFLAG
C
      RES=.0005
      IF(NTYPE.EQ.1)GOTO 500
C
  200 AMAX=0.0
      DO 400 L=1,NNODE
      IF(KNODE(L).EQ.1)GO TO 400
      OLDVAL=CNEW(L)
      SUM=0.0
      DO 300 JJ=1,NNODE
      IF(JJ.EQ.L)GO TO 300
      SUM=SUM+(AF*A(L,JJ)+P(L,JJ)/DT)*CNEW(JJ)
  300 CONTINUE
      S=AF*A(L,L)+P(L,L)/DT
      CNEW(L)=(-SUM+B(L))/S
      ERR=ABS(CNEW(L)-OLDVAL)
      IF(ERR.GT.AMAX)AMAX=ERR
  400 CONTINUE
      IF(AMAX.GT.RES)GO TO 200
      RETURN
C
  500 AMAX=0.0
      DO 850 L=1,NNODE
      IF(KNODE(L).EQ.1)GO TO 850
      OLDVAL=CNEW(L)
      SUM=0.0
      DO 900 JJ=1,NNODE
      IF(JJ.EQ.L)GO TO 900
      SUM=SUM+A(L,JJ)*CNEW(JJ)
  900 CONTINUE
      IF(A(L,L).EQ.0.)GO TO 850
      CNEW(L)=(-SUM+F(L))/A(L,L)
      ERR=ABS(CNEW(L)-OLDVAL)
      IF(ERR.GT.AMAX)AMAX=ERR
  850 CONTINUE
      IF(AMAX.GT.RES)GO TO 500
      RETURN
      END
C -----------------------------------------------------------------------
      SUBROUTINE GAUSSR(D,N)
      COMMON/VD/AF,AFM,DT,NNQS,NUMN,NUM,NNHC,NNST,NSTOP,KPRNT,IFLAG
      COMMON/VA/COLD(81),CNEW(81),KNODE(81),FIXED(81),E(81),AB(81)
      COMMON/VC/Q(81),F(81),NTS(81),NQS(81),DX(81),DY(81),DZ(81)
      COMMON/VV/A(81,81),P(81,81),R(81,81),B(81)
      DIMENSION S(81,81),RB(81),D(81)
C
      DO 1 I=1,N
      RB(I)=B(I)
      DO 1 J=1,N
    1 S(I,J)=R(I,J)
C
C     SETUP DIRICHLET VALUES
C
      DO 300 K=1,NNST
      I=NTS(K)
      DO 400 J=1,N
      IF(J.EQ.I)GOTO 400
      RB(J)=RB(J)-S(J,I)*FIXED(I)
      S(J,I)=0.
  400 CONTINUE
      DO 401 J=1,N
      IF(J.EQ.I)GOTO 401
      RB(I)=RB(I)-S(I,J)*FIXED(I)
      S(I,J)=0.
  401 CONTINUE
      S(I,I)=1.
      RB(I)=FIXED(I)
  300 CONTINUE
C
C     ELIMINATION ROUTINE
C
      DO 10 K=1,N
      K1=K+1
      RB(K)=RB(K)/S(K,K)
C
      IF(K.EQ.N)GOTO 100
      DO 20 J=K1,N
      IF(S(K,J).EQ.0)GOTO 20
      S(K,J)=S(K,J)/S(K,K)
      DO 30 I=K1,N
      S(I,J)=S(I,J)-S(I,K)*S(K,J)
   30 CONTINUE
      RB(J)=RB(J)-S(J,K)*RB(K)
   20 CONTINUE
   10 CONTINUE
C
C     BACKWARD SUBSTITUTION
C
  100 K1=K
      K=K-1
      IF(K.EQ.0)GOTO 200
      DO 40 J=K1,N
      RB(K)=RB(K)-S(K,J)*RB(J)
   40 CONTINUE
      GOTO 100
  200 CONTINUE
C
      DO 2 K=1,N
    2 D(K)=RB(K)
      RETURN
      END
