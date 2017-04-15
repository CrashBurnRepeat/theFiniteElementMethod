C  *********************************************************************
C  **                                                                 **
C  **                            FEM-3D.FOR                           **
C  **                   A PROGRAM FOR SOLVING THE 3-D                 **
C  **                          TIME-DEPENDENT                         **
C  **                     ADVECTION-DIFFUSION EQUATION                **
C  **                                                                 **
C  **                        WRITTEN BY D.W.PEPPER                    **
C  **                                                                 **
C  **       THIS PROGRAM USES A CHOLESKY SKYLINE MATRIX SOLVER        **
C  **       THE SHAPE FUNCTION IS THE 8 NODE LINEAR HEXAHEDRON        **
C  **                                                                 **
C  *********************************************************************
      COMMON/VA/COLD(51),CNEW(51)
      COMMON/VD/AF,AFM,DT,NNQS,NUMN,NUM,NGAUS,NNHC,NNST,NSTOP,KPRNT
      COMMON/VK/NNODE,NELEM,NTYPE,NTIME,TIME,NFRMAX
      CALL INIT
      CALL GAUSS
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
      COMMON/VA/COLD(51),CNEW(51)
      COMMON/VC/Q(51),F(51),NTS(51),NQS(51),DX(51),DY(51),DZ(51)
      COMMON/VB/A(6000),C(6000),G(6000)
      COMMON/VZ/X(51),Y(51),Z(51),QQ(51)
      COMMON/VEL/VX(51),VY(51),VZ(51)
      COMMON/VD/AF,AFM,DT,NNQS,NUMN,NUM,NGAUS,NNHC,NNST,NSTOP,KPRNT
      COMMON/VG/ISI(51),LEM(51),LME(51),ISIH(51)
      COMMON/VH/NODES4(6,4)
      COMMON/VM/H(51),TINF(51)
      COMMON/VK/NNODE,NELEM,NTYPE,NTIME,TIME,NFRMAX
      COMMON/VL/NODE(51,8)
      COMMON/TI/TITLE(18)
      COMMON/SK/NFR(51),JMIN(51),JMAX(51)
      CHARACTER*4 STOP(4),WORD
      DATA STOP/'DIRC','FLUX','CONV','VELC'/
C
      CALL SETPC
C
      READ(5,100)TITLE
  100 FORMAT(18A4)
      WRITE(*,101)TITLE
  101 FORMAT('0',5X,18A4)
C
      READ(5,*)NNODE,NELEM,NUMN,NSTOP,KPRNT,NVEL
      READ(5,*)NTYPE,DT,AF
        NGAUS=8
        NUM=4
      DO 33 I=1,NNODE
      CNEW(I)=0.
      COLD(I)=0.
      ISI(I)=0
      LEM(I)=0
      LME(I)=0
      ISIH(I)=0
      NTS(I)=0
      NQS(I)=0
      X(I)=0.
      Y(I)=0.
      Z(I)=0.
      QQ(I)=0.
      Q(I)=0.
      F(I)=0.
      H(I)=0.
      TINF(I)=0.
      DX(I)=0.
      DY(I)=0.
      DZ(I)=0.
      VX(I)=0.
      VY(I)=0.
      VZ(I)=0.
      DO 34 K=1,8
   34 NODE(I,K)=0
   33 CONTINUE
C
      DO 104 J=1,NNODE
      READ(5,*)I,X(I),Y(I),Z(I)
  104 CONTINUE
C
      DO 105 I=1,NELEM
      READ(5,*)J,QQ(J),DX(J),DY(J),DZ(J),(NODE(J,MM),MM=1,NUMN)
  105 CONTINUE
C
      DO 110 I=1,NNODE
      READ(5,1015)WORD,NT,CNT
      IF (WORD.EQ.STOP(1))GOTO 120
      NTS(I)=NT
      COLD(NTS(I))=CNT
  110 CONTINUE
  120 NNST=I-1
C
      DO 130 I=1,NNODE
      READ(5,1016)WORD,Q(I),LEM(I),ISI(I)
      IF(WORD.EQ.STOP(2))GO TO 135
  130 CONTINUE
  135 NNQS=I-1
C
      DO 140 I=1,NNODE
      READ(5,1017)WORD,H(I),TINF(I),LME(I),ISIH(I)
      IF(WORD.EQ.STOP(3))GO TO 145
  140 CONTINUE
  145 NNHC=I-1
C
      DO 146 J=1,NNODE
      READ(5,1018)WORD,I,VX(I),VY(I),VZ(I)
      IF(WORD.EQ.STOP(4))GO TO 147
  146 CONTINUE
  147 CONTINUE
C
      IF(NTYPE.EQ.1)WRITE(*,14)NTYPE
      IF(NTYPE.EQ.2)WRITE(*,15)NTYPE
      AFM=AF-1.0
      IF(NTYPE.EQ.1)GO TO 1
      IF(AF.EQ.1.)WRITE(*,10)DT
      IF(AF.EQ.0.5)WRITE(*,11)DT
    1 WRITE(*,12)NSTOP,KPRNT
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
      ISI1=ISI(I)
      WRITE(*,1085)I,Q(I),LEM(I),ISI1,
     *             (NODE(LEM(I),NODES4(ISI1,K)),K=1,NUM)
  170 CONTINUE
      WRITE(*,1081)
      WRITE(*,1083)
      DO 180 I=1,NNHC
      ISI1=ISIH(I)
      WRITE(*,1086)I,H(I),TINF(I),LME(I),ISI1,
     *            (NODE(LME(I),NODES4(ISI1,K)),K=1,NUM)
  180 CONTINUE
C
      IF(NVEL.EQ.0)GOTO 190
      WRITE(*,4040)
      WRITE(*,4041)
      DO 450 I=1,NNODE
  450 WRITE(*,1045)I,VX(I),VY(I),VZ(I)
C
C     CALCULATE BANDWIDTH FOR EACH ROW
C
  190 DO 620 I=1,NNODE
      JMIN(I)=NNODE
  620 JMAX(I)=1
      DO 630 K=1,NELEM
      DO 630 N=1,NUMN
      DO 630 M=1,NUMN
      I=NODE(K,N)
      J=NODE(K,M)
      IF(J.LT.JMIN(I))JMIN(I)=J
      IF(J.GT.JMAX(I))JMAX(I)=J
  630 CONTINUE
      NN1=NNODE-1
      DO 635 K=1,NN1
      KP1=K+1
      DO 635 I=KP1,NNODE
      IF(JMAX(I).LT.JMAX(K))JMAX(I)=JMAX(K)
  635 CONTINUE
C
C     TRANSFORM GLOBAL MATRIX TO A SINGLE INDEX
C
C     NFR(K) INDICATES THE POSITION OF FIRST NONZERO MEMBER IN
C     ROW K OF GLOBAL MATRIX
C
      NFR(1)=JMIN(1)
      DO 700 K=2,NNODE
  700 NFR(K)=NFR(K-1)+JMAX(K-1)-JMIN(K-1)+1
      NFRMAX=NFR(NNODE)+JMAX(NNODE)-JMIN(NNODE)
      WRITE(*,29)NFRMAX
   29 FORMAT(/,10X,'MAXIMUM 1-D ARRAY =',I6,/)
      IF(NFRMAX.GT.6000)THEN
        WRITE(*,30)
   30 FORMAT(/,10X,'SOLUTION STOPPED - MAXIMUM 1-D ARRAY EXCEEDS 6000')
         STOP
        ENDIF
C
      WRITE(*,191)
  191 FORMAT('0',5X,'PLEASE WAIT, SOLUTION IS PROCEEDING')
C
   10 FORMAT('0',10X,'FULLY IMPLICIT ',2X,'DT=',F7.4)
   11 FORMAT('0',10X,'CRANK NICOLSON ',2X,'DT=',F7.4)
   12 FORMAT('0',5X,'NSTOP=',I5,2X,'KPRNT=',I4)
   14 FORMAT('0',10X,'NTYPE=',I2,2X,'STEADY STATE CNEW')
   15 FORMAT('0',10X,'NTYPE=',I2,2X,'TIME DEPENDENT CNEW')
 1015 FORMAT(6X,A4,I5,5X,F10.5)
 1016 FORMAT(6X,A4,F10.5,2I5)
 1017 FORMAT(6X,A4,2F10.5,2I5)
 1018 FORMAT(6X,A4,I5,5X,3F10.5)
 1035 FORMAT('0',10X,'NO. OF NODES=',I4,5X,'NO.OF ELEMENTS=',I4)
 1040 FORMAT('0',8X,'NODAL COORDINATE VALUES')
 1041 FORMAT('0',7X,'I',9X,'X',11X,'Y',11X,'Z')
 1045 FORMAT(5X,I4,3(4X,F8.4))
 1050 FORMAT('0',1X,'ELEMENT',3X,'SOURCE',4X,'DX',6X,'DY',6X,'DZ',8X,'NO
     $DE NUMBERS')
 1055 FORMAT(2X,I4,3X,4F8.3,8(1X,I3))
 1060 FORMAT('0',7X,'NODES WHERE CNEW IS SPECIFIED')
 1061 FORMAT('0',8X,'I',5X,'NODE',6X,'CNEW')
 1065 FORMAT(2X,2(4X,I4),5X,F8.3)
 1080 FORMAT('0',15X,'NODES WHERE FLUX IS SPECIFIED')
 1081 FORMAT('0',15X,'NODES WHERE CONV IS SPECIFIED')
 1082 FORMAT('0',8X,'I',6X,'FLUX',5X,'ELE',2X,'SIDE',2X,'NODE NUMBERS')
 1083 FORMAT('0',10X,'I',7X,'H',7X,'TINF',7X,'ELE',2X,'SIDE',2X,'NODE NU
     *MBERS')
 1085 FORMAT(6X,I4,3X,F8.3,2X,2I5,4(1X,I4))
 1086 FORMAT(8X,I4,2X,2(2X,F8.3),2X,2I5,4(1X,I4))
 4040 FORMAT('0',5X,'SUMMARY OF SPECIFIED NODAL VELOCITIES ')
 4041 FORMAT('0',7X,'I',12X,'VX',12X,'VY',12X,'VZ')
      RETURN
      END
C ----------------------------------------------------------------------
      SUBROUTINE ASSEMB
      REAL NS,NX,NY,NZ,MASS
      COMMON/VA/COLD(51),CNEW(51)
      COMMON/VB/A(6000),C(6000),G(6000)
      COMMON/VZ/X(51),Y(51),Z(51),QQ(51)
      COMMON/VK/NNODE,NELEM,NTYPE,NTIME,TIME,NFRMAX
      COMMON/VL/NODE(51,8)
      COMMON/VC/Q(51),F(51),NTS(51),NQS(51),DX(51),DY(51),DZ(51)
      COMMON/VEL/VX(51),VY(51),VZ(51)
      COMMON/VD/AF,AFM,DT,NNQS,NUMN,NUM,NGAUS,NNHC,NNST,NSTOP,KPRNT
      COMMON/VE/NS(8),NX(8),NY(8),NZ(8),NXSI(8),NETA(8),NZTA(8)
      COMMON/VF/POSG1(8),POSG2(8),POSG3(8),W1(8),W2(8),W3(8)
      COMMON/VBC/POSB1(6,4),POSB2(6,4),WS1(6,4),WS2(6,4)
      COMMON/VG/ISI(51),LEM(51),LME(51),ISIH(51)
      COMMON/VM/H(51),TINF(51)
      COMMON/VH/NODES4(6,4)
      COMMON/SK/NFR(51),JMIN(51),JMAX(51)
      DIMENSION E(8),XS(4),YS(4),ZS(4)
C
      CALL CLEAR(G,NFRMAX)
      CALL CLEAR(A,NFRMAX)
      CALL CLEAR(C,NFRMAX)
C
C     CONSTRUCT AND STORE MASS AND STIFFNESS MATRICES
C
      DO 1 K=1,NELEM
      CALL NODSET(K,I,J,M,N,II,JJ,MM,NN)
      DO 1 IQ=1,NGAUS
      XSI=POSG1(IQ)
      ETA=POSG2(IQ)
      ZTA=POSG3(IQ)
      CALL SHAPE(K,XSI,ETA,ZTA,I,J,M,N,II,JJ,MM,NN,DET)
      DO 1 KK=1,NUMN
      L=NODE(K,KK)
C
C     SOURCE TERM
C
      F(L)=F(L)+NS(KK)*QQ(K)*DET*W1(IQ)*W2(IQ)*W3(IQ)
C
      DO 1 KKK=1,NUMN
      LL=NODE(K,KKK)
C
C     MASS MATRIX TERM
C
      MASS=NS(KKK)*NS(KK)
C
C     ADVECTION TERM
C
       VELX=NX(KKK)*NS(KK)*VX(L)
        VELY=NY(KKK)*NS(KK)*VY(L)
         VELZ=NZ(KKK)*NS(KK)*VZ(L)
          ADVEC=VELX+VELY+VELZ
C
C     DIFFUSION TERM
C
       DIFFX=NX(KKK)*NX(KK)*DX(K)
        DIFFY=NY(KKK)*NY(KK)*DY(K)
         DIFFZ=NZ(KKK)*NZ(KK)*DZ(K)
          DIFF=DIFFX+DIFFY+DIFFZ
C
      IJPT=NFR(L)+LL-JMIN(L)
      C(IJPT)=C(IJPT)+MASS*DET*W1(IQ)*W2(IQ)*W3(IQ)
    1 A(IJPT)=A(IJPT)+(DIFF+ADVEC)*DET*W1(IQ)*W2(IQ)*W3(IQ)
C
C     CONVECTIVE HEAT FLUX INFLUENCE
C
      IF(NNHC.EQ.0)GO TO 4
      DO 3 K=1,NNHC
      NEL=LME(K)
      ISI1=ISIH(K)
      CALL BCSET(NEL,ISI1,XS,YS,ZS)
      DO 3 IW=1,4
      CALL BCSIDE(IW,XS,YS,ZS,E,SIDE)
      DO 3 L1=1,NUM
         KN=NODES4(ISI1,L1)
      L=NODE(NEL,KN)
      DO 3 L2=1,NUM
         KK=NODES4(ISI1,L2)
      LL=NODE(NEL,KK)
      MASS=E(L1)*E(L2)
      IJPT=NFR(L)+LL-JMIN(L)
    3 A(IJPT)=A(IJPT)+MASS*H(K)*SIDE
    4 CONTINUE
C
C     FORM LHS
C
      DO 5 I=1,NFRMAX
    5 G(I)=AF*A(I)+C(I)/DT
C
      RETURN
      END
C ----------------------------------------------------------------------
      SUBROUTINE MATRIX
      COMMON/VB/A(6000),C(6000),G(6000)
      COMMON/VA/COLD(51),CNEW(51)
      COMMON/VZ/X(51),Y(51),Z(51),QQ(51)
      COMMON/VK/NNODE,NELEM,NTYPE,NTIME,TIME,NFRMAX
      COMMON/VC/Q(51),F(51),NTS(51),NQS(51),DX(51),DY(51),DZ(51)
      COMMON/VD/AF,AFM,DT,NNQS,NUMN,NUM,NGAUS,NNHC,NNST,NSTOP,KPRNT
      COMMON/SK/NFR(51),JMIN(51),JMAX(51)
      DIMENSION B(51),FIXED(51)
C
      CALL CLEAR(B,NNODE)
C
C     CONSTRUCT RHS
C
      DO 4 L=1,NNODE
      B(L)=F(L)
      JSTR=JMIN(L)
      JEND=JMAX(L)
      DO 4 LL=JSTR,JEND
      NGI=NFR(L)+LL-JSTR
      RHS=(AFM*A(NGI)+C(NGI)/DT)*COLD(LL)
    4 B(L)=B(L)+RHS
C
      DO 5 I=1,NNST
    5 FIXED(I)=COLD(NTS(I))
C
      CALL SKYLIN(B,CNEW,NTS,FIXED,NNST)
C
      RETURN
      END
C ----------------------------------------------------------------------
      BLOCK DATA
      COMMON/VBC/POSB1(6,4),POSB2(6,4),WS1(6,4),WS2(6,4)
      COMMON/VH/NODES4(6,4)
      DATA WS1/24*1.0/,WS2/24*1.0/
      DATA POSB1/
     *-0.5773502,0.5773502,0.5773502,-0.5773502,0.5773502,0.5773502,
     *0.5773502,-0.5773502,-0.5773502,0.5773502,-0.5773502,-0.5773502,
     *0.5773502,-0.5773502,-0.5773502,0.5773502,-0.5773502,-0.5773502,
     *-0.5773502,0.5773502,0.5773502,-0.5773502,0.5773502,0.5773502/
      DATA POSB2/
     *-0.5773502,-0.5773502,-0.5773502,-0.5773502,0.5773502,-0.5773502,
     *-0.5773502,-0.5773502,-0.5773502,-0.5773502,0.5773502,-0.5773502,
     *0.5773502,0.5773502,0.5773502,0.5773502,-0.5773502,0.5773502,
     *0.5773502,0.5773502,0.5773502,0.5773502,-0.5773502,0.5773502/
      DATA NODES4/1,6,5,5,8,5,
     *            2,7,6,8,4,1,
     *            3,3,7,4,3,2,
     *            4,2,8,1,7,6/
      END
C ----------------------------------------------------------------------
      SUBROUTINE BNDCON
      REAL NS,NX,NY,NZ
      COMMON/VK/NNODE,NELEM,NTYPE,NTIME,TIME,NFRMAX
      COMMON/VC/Q(51),F(51),NTS(51),NQS(51),DX(51),DY(51),DZ(51)
      COMMON/VL/NODE(51,8)
      COMMON/VD/AF,AFM,DT,NNQS,NUMN,NUM,NGAUS,NNHC,NNST,NSTOP,KPRNT
      COMMON/VE/NS(8),NX(8),NY(8),NZ(8),NXSI(8),NETA(8),NZTA(8)
      COMMON/VG/ISI(51),LEM(51),LME(51),ISIH(51)
      COMMON/VBC/POSB1(6,4),POSB2(6,4),WS1(6,4),WS2(6,4)
      COMMON/VM/H(51),TINF(51)
      COMMON/VH/NODES4(6,4)
      DIMENSION E(8),XS(4),YS(4),ZS(4)
C
C     IW=1,4 DENOTES FOUR GAUSS POINTS PER SIDE
C
   10 IF(NNQS.EQ.0)GO TO 35
      DO 3 K=1,NNQS
      NEL=LEM(K)
      ISI1=ISI(K)
      CALL BCSET(NEL,ISI1,XS,YS,ZS)
      DO 3 IW=1,4
      CALL BCSIDE(IW,XS,YS,ZS,E,SIDE)
      DO 3 L=1,NUM
         KK=NODES4(ISI1,L)
      LL=NODE(NEL,KK)
    3 F(LL)=F(LL)+E(L)*SIDE*Q(K)
C
   35 IF(NNHC.EQ.0)RETURN
      DO 9 K=1,NNHC
      NEL=LME(K)
      ISI1=ISIH(K)
      CALL BCSET(NEL,ISI1,XS,YS,ZS)
      DO 9 IW=1,4
      CALL BCSIDE(IW,XS,YS,ZS,E,SIDE)
      DO 9 L=1,NUM
         KK=NODES4(ISI1,L)
      LL=NODE(NEL,KK)
    9 F(LL)=F(LL)+E(L)*SIDE*H(K)*TINF(K)
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
      SUBROUTINE NODSET(K,I,J,M,N,II,JJ,MM,NN)
      COMMON/VL/NODE(51,8)
      COMMON/VD/AF,AFM,DT,NNQS,NUMN,NUM,NGAUS,NNHC,NNST,NSTOP,KPRNT
C
C     8-NODE LINEAR HEXAHEDRON
C
      I=NODE(K,1)
      J=NODE(K,2)
      M=NODE(K,3)
      N=NODE(K,4)
      II=NODE(K,5)
      JJ=NODE(K,6)
      MM=NODE(K,7)
      NN=NODE(K,8)
      RETURN
      END
C ----------------------------------------------------------------------
      SUBROUTINE BCSET(LEM,ISI,XS,YS,ZS)
      INTEGER  INDEX,ISI,K,LEM
      REAL XS(*), YS(*), ZS(*)
C
      COMMON/VZ/X(51),Y(51),Z(51),QQ(51)
      COMMON/VH/NODES4(6,4)
      COMMON/VL/NODE(51,8)
C
      DO 10 K = 1, 4
              INDEX = NODE(LEM,NODES4(ISI,K))
              XS(K) = X(INDEX)
              YS(K) = Y(INDEX)
              ZS(K) = Z(INDEX)
  10  CONTINUE
      RETURN
      END
C ----------------------------------------------------------------------
      SUBROUTINE BCSIDE(IW,XS,YS,ZS,NE,SIDE)
      INTEGER IW
      REAL XS(*), YS(*), ZS(*), NE(*)
      REAL SSIDE(4), TSIDE(4), WS1(4), WS2(4)
      REAL NSI(4), NTI(4)
      REAL POINT1, POINT2, XSI,ETA
C
      PARAMETER(POINT2 = 0.5773502, POINT1 = -POINT2 )
      DATA SSIDE/POINT1,POINT2,POINT2,POINT1/
      DATA tSIDE/POINT1,POINT1,POINT2,POINT2/
      DATA WS1,WS2/8*1.0/
      DATA D/0.25/
C
          XSI  = SSIDE(IW)
          ETA  = TSIDE(IW)
          NE(1)= D*(1.-XSI)*(1.-ETA)
          NE(2)= D*(1.+XSI)*(1.-ETA)
          NE(3)= D*(1.+XSI)*(1.+ETA)
          NE(4)= D*(1.-XSI)*(1.+ETA)
C
          NSI(1)=-D*(1.-ETA)
          NSI(2)= D*(1.-ETA)
          NSI(3)= D*(1.+ETA)
          NSI(4)=-D*(1.+ETA)
C
          NTI(1)=-D*(1.-XSI)
          NTI(2)=-D*(1.+XSI)
          NTI(3)= D*(1.+XSI)
          NTI(4)= D*(1.-XSI)
C
          XXSI=0.0
          XETA=0.0
          YXSI=0.0
          YETA=0.0
          ZXSI=0.0
          ZETA=0.0
          DO 1 L=1,4
              XXSI=XXSI+NSI(L)*XS(L)
              XETA=XETA+NTI(L)*XS(L)
              YXSI=YXSI+NSI(L)*YS(L)
              YETA=YETA+NTI(L)*YS(L)
              ZXSI=ZXSI+NSI(L)*ZS(L)
              ZETA=ZETA+NTI(L)*ZS(L)
    1     CONTINUE
C
          G11=XXSI*XXSI+YXSI*YXSI+ZXSI*ZXSI
          G22=XETA*XETA+YETA*YETA+ZETA*ZETA
          G12=XXSI*XETA+YXSI*YETA+ZXSI*ZETA
C
          DET=ABS(G11*G22-G12*G12)
C
          IF(DET.EQ.0.0)THEN
              WRITE(*,100)K
 100          FORMAT(2X,'THE DETERMINANT = 0.0 IN ELEMENT',I4)
              STOP
          ENDIF
C
C     NORMAL UNIT VECTORS
C
          XNORM = (YXSI*ZETA-YETA*ZXSI)/DET
          YNORM = (XETA*ZXSI-XXSI*ZETA)/DET
          ZNORM = (XXSI*YETA-YXSI*XETA)/DET
C
          SNORM = SQRT(DET)
C
          SIDE=SNORM*WS1(IW)*WS2(IW)
      RETURN
      END
C ----------------------------------------------------------------------
      SUBROUTINE SHAPE(K,XSI,ETA,ZTA,I,J,M,N,II,JJ,MM,NN,DET)
      REAL NS,NX,NY,NZ,NXSI,NETA,NZTA
      COMMON/VZ/X(51),Y(51),Z(51),QQ(51)
      COMMON/VD/AF,AFM,DT,NNQS,NUMN,NUM,NGAUS,NNHC,NNST,NSTOP,KPRNT
      COMMON/VE/NS(8),NX(8),NY(8),NZ(8),NXSI(8),NETA(8),NZTA(8)
C
C     8-NODE HEXAHEDRON ELEMENT
C
      D=0.125
C
      NS(1)=D*(1.-XSI)*(1.-ETA)*(1.+ZTA)
      NS(2)=D*(1.+XSI)*(1.-ETA)*(1.+ZTA)
      NS(3)=D*(1.+XSI)*(1.+ETA)*(1.+ZTA)
      NS(4)=D*(1.-XSI)*(1.+ETA)*(1.+ZTA)
      NS(5)=D*(1.-XSI)*(1.-ETA)*(1.-ZTA)
      NS(6)=D*(1.+XSI)*(1.-ETA)*(1.-ZTA)
      NS(7)=D*(1.+XSI)*(1.+ETA)*(1.-ZTA)
      NS(8)=D*(1.-XSI)*(1.+ETA)*(1.-ZTA)
C
      NXSI(1)=D*(ETA-1.)*(1.+ZTA)
      NXSI(2)=D*(1.-ETA)*(1.+ZTA)
      NXSI(3)=D*(1.+ETA)*(1.+ZTA)
      NXSI(4)=D*(-1.-ETA)*(1.+ZTA)
      NXSI(5)=D*(ETA-1.)*(1.-ZTA)
      NXSI(6)=D*(1.-ETA)*(1.-ZTA)
      NXSI(7)=D*(1.+ETA)*(1.-ZTA)
      NXSI(8)=D*(-1.-ETA)*(1.-ZTA)
C
      NETA(1)=D*(XSI-1.)*(1.+ZTA)
      NETA(2)=D*(-1.-XSI)*(1.+ZTA)
      NETA(3)=D*(1.+XSI)*(1.+ZTA)
      NETA(4)=D*(1.-XSI)*(1.+ZTA)
      NETA(5)=D*(XSI-1.)*(1.-ZTA)
      NETA(6)=D*(-1.-XSI)*(1.-ZTA)
      NETA(7)=D*(1.+XSI)*(1.-ZTA)
      NETA(8)=D*(1.-XSI)*(1.-ZTA)
C
      NZTA(1)=D*(1.-XSI)*(1.-ETA)
      NZTA(2)=D*(1.+XSI)*(1.-ETA)
      NZTA(3)=D*(1.+XSI)*(1.+ETA)
      NZTA(4)=D*(1.-XSI)*(1.+ETA)
      NZTA(5)=D*(XSI-1.)*(1.-ETA)
      NZTA(6)=D*(-1.-XSI)*(1.-ETA)
      NZTA(7)=D*(-1.-XSI)*(1.+ETA)
      NZTA(8)=D*(XSI-1.)*(1.+ETA)
C
      ASX=NXSI(1)*X(I)+NXSI(2)*X(J)+NXSI(3)*X(M)+NXSI(4)*X(N)+
     *    NXSI(5)*X(II)+NXSI(6)*X(JJ)+NXSI(7)*X(MM)+NXSI(8)*X(NN)
      ASY=NXSI(1)*Y(I)+NXSI(2)*Y(J)+NXSI(3)*Y(M)+NXSI(4)*Y(N)+
     *    NXSI(5)*Y(II)+NXSI(6)*Y(JJ)+NXSI(7)*Y(MM)+NXSI(8)*Y(NN)
      ASZ=NXSI(1)*Z(I)+NXSI(2)*Z(J)+NXSI(3)*Z(M)+NXSI(4)*Z(N)+
     *    NXSI(5)*Z(II)+NXSI(6)*Z(JJ)+NXSI(7)*Z(MM)+NXSI(8)*Z(NN)
C
      ATX=NETA(1)*X(I)+NETA(2)*X(J)+NETA(3)*X(M)+NETA(4)*X(N)+
     *    NETA(5)*X(II)+NETA(6)*X(JJ)+NETA(7)*X(MM)+NETA(8)*X(NN)
      ATY=NETA(1)*Y(I)+NETA(2)*Y(J)+NETA(3)*Y(M)+NETA(4)*Y(N)+
     *    NETA(5)*Y(II)+NETA(6)*Y(JJ)+NETA(7)*Y(MM)+NETA(8)*Y(NN)
      ATZ=NETA(1)*Z(I)+NETA(2)*Z(J)+NETA(3)*Z(M)+NETA(4)*Z(N)+
     *    NETA(5)*Z(II)+NETA(6)*Z(JJ)+NETA(7)*Z(MM)+NETA(8)*Z(NN)
C
      ARX=NZTA(1)*X(I)+NZTA(2)*X(J)+NZTA(3)*X(M)+NZTA(4)*X(N)+
     *    NZTA(5)*X(II)+NZTA(6)*X(JJ)+NZTA(7)*X(MM)+NZTA(8)*X(NN)
      ARY=NZTA(1)*Y(I)+NZTA(2)*Y(J)+NZTA(3)*Y(M)+NZTA(4)*Y(N)+
     *    NZTA(5)*Y(II)+NZTA(6)*Y(JJ)+NZTA(7)*Y(MM)+NZTA(8)*Y(NN)
      ARZ=NZTA(1)*Z(I)+NZTA(2)*Z(J)+NZTA(3)*Z(M)+NZTA(4)*Z(N)+
     *    NZTA(5)*Z(II)+NZTA(6)*Z(JJ)+NZTA(7)*Z(MM)+NZTA(8)*Z(NN)
C
      A11=ATY*ARZ-ATZ*ARY
      A12=ATZ*ARX-ATX*ARZ
      A13=ATX*ARY-ATY*ARX
C
      A21=ASZ*ARY-ASY*ARZ
      A22=ASX*ARZ-ASZ*ARX
      A23=ASY*ARX-ASX*ARY
C
      A31=ASY*ATZ-ASZ*ATY
      A32=ASZ*ATX-ASX*ATZ
      A33=ASX*ATY-ASY*ATX
C
      DET=A11*ASX+A12*ASY+A13*ASZ
C
      IF(DET.EQ.0.0)THEN
      WRITE(*,100)K
  100 FORMAT(2X,'THE DETERMINANT = 0.0 IN ELEMENT',I4)
      STOP
      ENDIF
C
      DO 8 L=1,NUMN
      NX(L)=(A11*NXSI(L)+A21*NETA(L)+A31*NZTA(L))/DET
      NY(L)=(A12*NXSI(L)+A22*NETA(L)+A32*NZTA(L))/DET
      NZ(L)=(A13*NXSI(L)+A23*NETA(L)+A33*NZTA(L))/DET
    8 CONTINUE
      RETURN
      END
C ----------------------------------------------------------------------
      SUBROUTINE GAUSS
      COMMON/VD/AF,AFM,DT,NNQS,NUMN,NUM,NGAUS,NNHC,NNST,NSTOP,KPRNT
      COMMON/VF/POSG1(8),POSG2(8),POSG3(8),W1(8),W2(8),W3(8)
C
C     8 NODE HEXAHEDRON
C
    4 W=1./SQRT(3.)
      POSG1(1)=-W
      POSG1(2)= W
      POSG1(3)= W
      POSG1(4)=-W
      POSG1(5)=-W
      POSG1(6)= W
      POSG1(7)= W
      POSG1(8)=-W
C
      POSG2(1)=-W
      POSG2(2)=-W
      POSG2(3)= W
      POSG2(4)= W
      POSG2(5)=-W
      POSG2(6)=-W
      POSG2(7)= W
      POSG2(8)= W
C
      POSG3(1)= W
      POSG3(2)= W
      POSG3(3)= W
      POSG3(4)= W
      POSG3(5)=-W
      POSG3(6)=-W
      POSG3(7)=-W
      POSG3(8)=-W
C
      W1(1)=1.0
      W1(2)=1.0
      W1(3)=1.0
      W1(4)=1.0
      W1(5)=1.0
      W1(6)=1.0
      W1(7)=1.0
      W1(8)=1.0
C
      W2(1)=1.0
      W2(2)=1.0
      W2(3)=1.0
      W2(4)=1.0
      W2(5)=1.0
      W2(6)=1.0
      W2(7)=1.0
      W2(8)=1.0
C
      W3(1)=1.0
      W3(2)=1.0
      W3(3)=1.0
      W3(4)=1.0
      W3(5)=1.0
      W3(6)=1.0
      W3(7)=1.0
      W3(8)=1.0
      RETURN
      END
C ----------------------------------------------------------------------
      SUBROUTINE PRINT
      COMMON/VA/COLD(51),CNEW(51)
      COMMON/VK/NNODE,NELEM,NTYPE,NTIME,TIME,NFRMAX
      COMMON/VEL/VX(51),VY(51),VZ(51)
      IF(NTYPE.EQ.1)THEN
      WRITE(*,403)
      ELSE
      WRITE(*,402)NTIME,TIME
      ENDIF
      WRITE(*,400)
      DO 1 I=1,NNODE
    1 WRITE(*,401)I,CNEW(I),VX(I),VY(I),VZ(I)
  401 FORMAT(5X,I3,5X,4(F8.3))
  402 FORMAT('0',5X,'TIME STEPS=',I3,2X,'TIME=',F8.3)
  403 FORMAT('0',5X,'STEADY STATE')
  400 FORMAT('0',5X,'NODE',5X,'CNEW',8X,'VX',6X,'VY',6X,'VZ')
      RETURN
      END
C ----------------------------------------------------------------------
      SUBROUTINE RESID
      COMMON/VK/NNODE,NELEM,NTYPE,NTIME,TIME,NFRMAX
      COMMON/VA/COLD(51),CNEW(51)
      INTEGER LOC,ERRPRN
      REAL MAXERR
      DATA ERRPRN/10/
C
      MAXERR = 0.0
      MAXRES=0
      ERR=0.0001
      DO 1 I=1,NNODE
      R=ABS(CNEW(I)-COLD(I))
      MAXERR=MAX(R,MAXERR)
      IF(R.GT.ERR)THEN
          MAXRES = 1
          LOC    = I
      ENDIF
    1 CONTINUE
C
      IF(MOD(NTIME,ERRPRN).EQ.0)THEN
         WRITE(*,101)NTIME,MAXERR,LOC
      ENDIF
      IF(MAXRES.EQ.1)THEN
      RETURN
      ELSE
      WRITE(*,100)NTIME
      CALL PRINT
      CALL PRINTF
      END IF
  100 FORMAT('0',2X,'PROGRAM HAS CONVERGED IN',I3,' STEPS')
  101 FORMAT('0',2X,'TIME=',I4,2X,'ERROR=',1PE12.5,2X,'NODE=',I4)
      STOP
      END
C ----------------------------------------------------------------------
      SUBROUTINE PRINTF
      COMMON/VK/NNODE,NELEM,NTYPE,NTIME,TIME,NFRMAX
      COMMON/VA/COLD(51),CNEW(51)
      COMMON/VZ/X(51),Y(51),Z(51),QQ(51)
      COMMON/VEL/VX(51),VY(51),VZ(51)
      COMMON/TI/TITLE(18)
      WRITE(8,'(18A)')TITLE
      WRITE(8,'(I4)')NTIME
      WRITE(8,'(2I4)')NNODE,NELEM
      DO 1 I=1,NNODE
    1 WRITE(8,100)I,X(I),Y(I),Z(I),CNEW(I),VX(I),VY(I),VZ(I)
  100 FORMAT(5X,I3,7(F8.3,1X))
      RETURN
      END
C ----------------------------------------------------------------------
      SUBROUTINE SETPC
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
      RETURN
      END
C ----------------------------------------------------------------------
      SUBROUTINE SKYLIN(B,VAR,NB,FIXED,NBOUND)
      COMMON/VB/A(6000),C(6000),G(6000)
      COMMON/VK/NNODE,NELEM,NTYPE,NTIME,TIME,NFRMAX
      COMMON/SK/NFR(51),JMIN(51),JMAX(51)
      DIMENSION AG(6000),B(51),VAR(51),NB(51),FIXED(51)
      NN1=NNODE-1
      DO 1 I=1,NFRMAX
    1 AG(I)=G(I)
C
C     CHECK FOR DIRICHLET VALUES
C
   14 IF(NBOUND.EQ.0)GO TO 700
      DO 550 K=1,NBOUND
      I=NB(K)
      NGA=NFR(I)+I-JMIN(I)
      AG(NGA)=1.0
      J1=JMIN(I)
      J2=JMAX(I)
      DO 550 J=J1,J2
      IF(J.EQ.I)GO TO 550
      NGA=NFR(I)+J-J1
      AG(NGA)=0.0
      IF(I.LT.JMIN(J))GO TO 550
      NGA=NFR(J)+I-JMIN(J)
      B(J)=B(J)-AG(NGA)*FIXED(K)
      AG(NGA)=0.0
  550 CONTINUE
      DO 600 K=1,NBOUND
      I=NB(K)
  600 B(I)=FIXED(K)
C
C     BACK CALCULATE FOR VARIABLE
C
  700 DO 50 K=1,NN1
      NGA=NFR(K)+K-JMIN(K)
      PIVOT=AG(NGA)
      JKMAX=JMAX(K)
      KP1=K+1
      DO 20 I=KP1,JKMAX
      IF(K.LT.JMIN(I))GO TO 20
      NGA=NFR(I)+K-JMIN(I)
      FACTR=AG(NGA)/PIVOT
      IF(ABS(FACTR).LT.1.0E-10)GO TO 20
      DO 10 J=K,JKMAX
      NGE=NFR(I)+J-JMIN(I)
      NGI=NFR(K)+J-JMIN(K)
      AG(NGE)=AG(NGE)-FACTR*AG(NGI)
   10 CONTINUE
      B(I)=B(I)-FACTR*B(K)
   20 CONTINUE
   50 CONTINUE
      DO 40 II=1,NN1
      I=NNODE+1-II
      IM1=I-1
      JSTR=JMIN(I)
      NGA=NFR(I)+I-JSTR
      VAR(I)=B(I)/AG(NGA)
      DO 30 J=JSTR,IM1
      NGI=NFR(J)+I-JMIN(J)
   30 B(J)=B(J)-VAR(I)*AG(NGI)
   40 CONTINUE
      NGA=NFR(1)+1-JMIN(1)
      VAR(1)=B(1)/AG(NGA)
      RETURN
      END
