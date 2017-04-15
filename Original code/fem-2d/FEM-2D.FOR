C  *********************************************************************
C  **                                                                 **
C  **                           FEM-2D.FOR                            **
C  **                     A 2-D FEM PROGRAM WHICH                     **
C  **                     ACCOMPANIES THE TEXTBOOK                    **
C  **                                                                 **
C  **                "THE FINITE ELEMENT METHOD: BASIC                **
C  **                    CONCEPTS AND APPLICATIONS"                   **
C  **                                BY                               **
C  **              DARRELL W. PEPPER AND JUAN C. HEINRICH             **
C  **                  HEMISPHERE PUBLISHING COMPANY                  **
C  **                                                                 **
C  **                        X-Y OR R-Z COORDINATES                   **
C  **                                                                 **
C  **               THIS PROGRAM USES CHOLESKY (SKYLINE)              **
C  **                 DECOMPOSITION TO SOLVE MATRICES                 **
C  **                                                                 **
C  **    (NOTE: THIS PROGRAM DIMENSIONED TO A MAXIMUM OF 201 NODES)   **
C  **                                                                 **
C  **                                 NOV. 30, 2005 D. W. PEPPER      **
C  *********************************************************************
      COMMON/VA/COLD(201),CNEW(201)
      COMMON/VD/AF,AFM,DT,NNQS,NUMN,NUM,NGAUS,NNHC,NNST,NSTOP,KPRNT,IAXI
      COMMON/VK/NNODE,NELEM,NTYPE,NTIME,TIME,NFRMAX
C
C     SET UP INITIAL DATA AND BOUNDARY CONDITIONS
C
      CALL INIT
      CALL GAUSS
      KOUNT=0
      NTIME=0
      TIME=DT
      IF(IAXI.EQ.0)THEN
       CALL MATSET
       CALL BNDCON
      ELSE
        CALL MATAXI
        CALL BNDAXI
      ENDIF
C
C     STEP THROUGH TIME
C
      DO 2 NSTEP=1,NSTOP
      CALL ASSEMB
      IF(KOUNT.EQ.KPRNT)THEN
        CALL PRINT
        KOUNT=0
      ENDIF
      TIME=TIME+DT
      KOUNT=KOUNT+1
      NTIME=NTIME+1
      CALL RESID
C
C     PUT NEW VALUES INTO OLD ARRAY FOR NEXT TIME STEP.
C
      DO 3 L=1,NNODE
    3 COLD(L)=CNEW(L)
C
    2 CONTINUE
      WRITE(*,10)
   10 FORMAT('0',1X,'SOLUTION IS FINISHED')
      STOP
      END
C ----------------------------------------------------------------------
      SUBROUTINE INIT
      COMMON/VP/P(4001),C(4001),R(4001)
      COMMON/VA/COLD(201),CNEW(201)
      COMMON/VB/X(201),Y(201),QQ(201)
      COMMON/VC/Q(201),F(201),NTS(201),NQS(201),DX(201),DY(201)
      COMMON/VEL/VX(201),VY(201)
      COMMON/VD/AF,AFM,DT,NNQS,NUMN,NUM,NGAUS,NNHC,NNST,NSTOP,KPRNT,IAXI
      COMMON/VG/ISI(201),LEM(201),LME(201),ISIH(201)
      COMMON/VH/NODES4(4,2),NODES8(4,3),NODES3(3,2),NODES6(3,3)
      COMMON/VM/H(201),TINF(201)
      COMMON/VK/NNODE,NELEM,NTYPE,NTIME,TIME,NFRMAX
      COMMON/VL/NODE(201,8)
      COMMON/SK/NFR(201),JMIN(201),JMAX(201)
      COMMON/VT/MTYPE,NUMDIM,TO,RHO,CP
      COMMON/TI/TITLE(18)
      CHARACTER*4 STOP(4),WORD
      DATA STOP/'DIRC','FLUX','CONV','VELC'/
C
      CALL SETPC
C
      READ(8,100)TITLE
  100 FORMAT(18A4)
      WRITE(*,101)TITLE
  101 FORMAT(5X,18A4)
C
      READ(8,*)MTYPE,NUMDIM,NNODE,NELEM,NUMN,NSTOP,KPRNT,NVEL
      READ(8,*)NTYPE,DT,AF,TO,RHO,CP,IAXI
      MTYPE=MTYPE+1
      IF(NUMN.EQ.3)NGAUS=3
      IF(NUMN.EQ.4)NGAUS=4
      IF(NUMN.EQ.6)NGAUS=7
      IF(NUMN.EQ.8)NGAUS=9
      NUM=2
      IF(NUMN.GE.6)NUM=3
      DO 33 I=1,NNODE
      CNEW(I)=0.0
      COLD(I)=TO
      ISI(I)=0
      LEM(I)=0
      LME(I)=0
      ISIH(I)=0
      NTS(I)=0
      NQS(I)=0
      X(I)=0.
      Y(I)=0.
      QQ(I)=0.
      Q(I)=0.
      F(I)=0.
      H(I)=0.
      TINF(I)=0.
      DX(I)=0.
      DY(I)=0.
      VX(I)=0.
      VY(I)=0.
      DO 34 K=1,NUMN
   34 NODE(I,K)=0
   33 CONTINUE
C
      READ(8,*)(I,X(I),Y(I),J=1,NNODE)
C
      DO 105 I=1,NELEM
      READ(8,*)J,QQ(J),DX(J),DY(J),(NODE(J,MM),MM=1,NUMN)
      DX(J)=DX(J)/(RHO*CP)
      DY(J)=DX(J)
      QQ(J)=QQ(J)/(RHO*CP)
  105 CONTINUE
C
      DO 110 I=1,NNODE
      READ(8,1015)WORD,NT,CNT
      IF (WORD.EQ.STOP(1))GOTO 120
      NTS(I)=NT
      COLD(NTS(I))=CNT
  110 CONTINUE
  120 NNST=I-1
C
      DO 130 I=1,NNODE
      READ(8,1016)WORD,Q(I),LEM(I),ISI(I)
      IF(WORD.EQ.STOP(2))GO TO 135
  130 CONTINUE
  135 NNQS=I-1
C
      DO 140 I=1,NNODE
      READ(8,1017)WORD,H(I),TINF(I),LME(I),ISIH(I)
      IF(WORD.EQ.STOP(3))GO TO 145
  140 CONTINUE
  145 NNHC=I-1
C
      DO 146 I=1,NNODE
      READ(8,1018)WORD,J,VXTEMP,VYTEMP
      IF(WORD.EQ.STOP(4))GO TO 147
	VX(J)=VXTEMP
	VY(J)=VYTEMP
  146 CONTINUE
  147 CONTINUE
C
      IF(NTYPE.EQ.1)WRITE(*,14)NTYPE
      IF(NTYPE.EQ.2)WRITE(*,15)NTYPE
      AFM=AF-1.0
      IF(NTYPE.EQ.1)GO TO 1
      IF(AF.EQ.1.)WRITE(*,10)DT
      IF(AF.EQ.0.5)WRITE(*,11)DT
    1 WRITE(*,12)NSTOP,KPRNT,TO,IAXI
      WRITE(*,1035)NNODE,NELEM
      WRITE(*,1040)
      WRITE(*,1041)
      DO 150 I=1,NNODE
      WRITE(*,1045)I,X(I),Y(I)
  150 CONTINUE
      WRITE(*,1050)
      DO 155 I=1,NELEM
      WRITE(*,1055)I,QQ(I),DX(I),DY(I),(NODE(I,MM),MM=1,NUMN)
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
      IF(NUMN.EQ.3.OR.NUMN.EQ.6)THEN
      IF(NUM.EQ.2)WRITE(*,1085)I,Q(I),(NODE(LEM(I),NODES3(ISI1,K)),K=1,N
     1UM)
      IF(NUM.EQ.3)WRITE(*,1085)I,Q(I),(NODE(LEM(I),NODES6(ISI1,K)),K=1,N
     1UM)
      ELSE
      IF(NUM.EQ.2)WRITE(*,1085)I,Q(I),(NODE(LEM(I),NODES4(ISI1,K)),K=1,N
     1UM)
      IF(NUM.EQ.3)WRITE(*,1085)I,Q(I),(NODE(LEM(I),NODES8(ISI1,K)),K=1,N
     1UM)
      ENDIF
  170 CONTINUE
      WRITE(*,1081)
      WRITE(*,1083)
      DO 180 I=1,NNHC
      ISI1=ISIH(I)
      IF(NUMN.EQ.3.OR.NUMN.EQ.6)THEN
      IF(NUM.EQ.2)WRITE(*,1086)I,H(I),TINF(I),(NODE(LME(I),NODES3(ISI1,K
     1)),K=1,NUM)
      IF(NUM.EQ.3)WRITE(*,1086)I,H(I),TINF(I),(NODE(LME(I),NODES6(ISI1,K
     1)),K=1,NUM)
      ELSE
      IF(NUM.EQ.2)WRITE(*,1086)I,H(I),TINF(I),(NODE(LME(I),NODES4(ISI1,K
     1)),K=1,NUM)
      IF(NUM.EQ.3)WRITE(*,1086)I,H(I),TINF(I),(NODE(LME(I),NODES8(ISI1,K
     1)),K=1,NUM)
      ENDIF
  180 CONTINUE
C
      IF(NVEL.EQ.0)GOTO 190
      WRITE(*,4040)
      WRITE(*,4041)
      DO 450 I=1,NNODE
  450 WRITE(*,1045)I,VX(I),VY(I)
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
      IF(NFRMAX.GT.4001)THEN
        WRITE(*,30)
         STOP
        ENDIF
C
      WRITE(*,191)
  191 FORMAT(/,5X,'PLEASE WAIT, SOLUTION IS PROCEEDING')
C
   10 FORMAT(/,10X,'FULLY IMPLICIT ',2X,'DT=',F7.4)
   11 FORMAT(/,10X,'CRANK NICOLSON ',2X,'DT=',F7.4)
   12 FORMAT(/,5X,'NSTOP=',I5,2X,'KPRNT=',I4,2X,'TO=',F7.2,2X,'AXI=',I4)
   14 FORMAT(/,10X,'NTYPE=',I2,2X,'STEADY STATE CNEW')
   15 FORMAT(/,10X,'NTYPE=',I2,2X,'TIME DEPENDENT CNEW')
   29 FORMAT(/,10X,'MAXIMUM 1-D ARRAY =',I10,/)
   30 FORMAT(/,10X,'SOLUTION STOPPED - MAXIMUM 1-D ARRAY EXCEEDS 4001')
 1015 FORMAT(6X,A4,I5,5X,F10.5)
 1016 FORMAT(6X,A4,F10.5,2I5)
 1017 FORMAT(6X,A4,2F10.5,2I5)
 1018 FORMAT(6X,A4,I5,5X,2F10.5)
 1035 FORMAT(/,10X,'NO. OF NODES=',I4,5X,'NO.OF ELEMENTS=',I4)
 1040 FORMAT(/,8X,'NODAL COORDINATE VALUES')
 1041 FORMAT(/,7X,'I',9X,'X',11X,'Y')
 1045 FORMAT(5X,I4,2(4X,F10.4))
 1050 FORMAT(/,3X,'ELEMENT',3X,'SOURCE',4X,'DX',6X,'DY',8X,'NODE NUMBE
     1RS')
 1055 FORMAT(5X,I4,3X,3F10.3,8(1X,I4))
 1060 FORMAT(/,7X,'NODES WHERE CNEW IS SPECIFIED')
 1061 FORMAT(/,8X,'I',5X,'NODE',6X,'CNEW')
 1065 FORMAT(2X,2(4X,I4),5X,F10.3)
 1080 FORMAT(/,15X,'NODES WHERE FLUX IS SPECIFIED')
 1081 FORMAT(/,15X,'NODES WHERE CONV IS SPECIFIED')
 1082 FORMAT(/,8X,'I',6X,'FLUX',5X,'NODE NUMBERS')
 1083 FORMAT(/,10X,'I',7X,'H',7X,'TINF',7X,'NODE NUMBERS')
 1085 FORMAT(6X,I4,3X,F10.3,2X,3(1X,I4))
 1086 FORMAT(8X,I4,2X,2(2X,F10.3),2X,3(1X,I4))
 4040 FORMAT(/,5X,'SUMMARY OF SPECIFIED NODAL VELOCITIES ')
 4041 FORMAT(/,7X,'I',12X,'VX',12X,'VY')
      RETURN
      END
C ----------------------------------------------------------------------
      SUBROUTINE MATSET
      REAL NS,NX,NY,MASS
      COMMON/VP/P(4001),C(4001),R(4001)
      COMMON/VA/COLD(201),CNEW(201)
      COMMON/VB/X(201),Y(201),QQ(201)
      COMMON/VC/Q(201),F(201),NTS(201),NQS(201),DX(201),DY(201)
      COMMON/VE/NS(8),NX(8),NY(8),NXSI(8),NETA(8)
      COMMON/VF/POSG1(9),POSG2(9),W1(9),W2(9)
      COMMON/VD/AF,AFM,DT,NNQS,NUMN,NUM,NGAUS,NNHC,NNST,NSTOP,KPRNT,IAXI
      COMMON/VG/ISI(201),LEM(201),LME(201),ISIH(201)
      COMMON/VH/NODES4(4,2),NODES8(4,3),NODES3(3,2),NODES6(3,3)
      COMMON/VEL/VX(201),VY(201)
      COMMON/VM/H(201),TINF(201)
      COMMON/VK/NNODE,NELEM,NTYPE,NTIME,TIME,NFRMAX
      COMMON/VL/NODE(201,8)
      COMMON/SK/NFR(201),JMIN(201),JMAX(201)
C
      CALL CLEAR(R,NFRMAX)
      CALL CLEAR(C,NFRMAX)
      CALL CLEAR(P,NFRMAX)
C
C     CONSTRUCT AND STORE MASS MATRIX P
C
      DO 1 K=1,NELEM
      CALL NODSET(K,I,II,J,JJ,M,MM,N,NN)
      DO 1 IQ=1,NGAUS
      XSI=POSG1(IQ)
      ETA=POSG2(IQ)
      CALL SHAPE(K,XSI,ETA,I,J,M,N,II,JJ,MM,NN,DET,YETA,XETA,XXSI,YXSI)
      DO 1 KK=1,NUMN
      L=NODE(K,KK)
      F(L)=F(L)+NS(KK)*QQ(K)*DET*W1(IQ)*W2(IQ)
      DO 1 KKK=1,NUMN
      LL=NODE(K,KKK)
      MASS=NS(KKK)*NS(KK)
      IJPT=NFR(L)+LL-JMIN(L)
    1 P(IJPT)=P(IJPT)+MASS*DET*W1(IQ)*W2(IQ)
C
C     CONSTRUCT AND STORE ADVECTION AND DIFFUSION MATRIX C
C
      DO 4 K=1,NELEM
      CALL NODSET(K,I,II,J,JJ,M,MM,N,NN)
      DO 4 IQ=1,NGAUS
      XSI=POSG1(IQ)
      ETA=POSG2(IQ)
      CALL SHAPE(K,XSI,ETA,I,J,M,N,II,JJ,MM,NN,DET,YETA,XETA,XXSI,YXSI)
      DO 4 KK=1,NUMN
      L=NODE(K,KK)
      DO 4 KKK=1,NUMN
      LL=NODE(K,KKK)
C
C     ADVECTION TERM
C
       VELX=VX(L)*NS(KK)*NX(KKK)
        VELY=VY(L)*NS(KK)*NY(KKK)
         ADVEC=VELX+VELY
C
C     DIFFUSION TERM
C
       DIFFX=DX(K)*NX(KKK)*NX(KK)
        DIFFY=DY(K)*NY(KKK)*NY(KK)
         DIFF=DIFFX+DIFFY
C
      IJPT=NFR(L)+LL-JMIN(L)
    4 C(IJPT)=C(IJPT)+(DIFF+ADVEC)*DET*W1(IQ)*W2(IQ)
C
C     CONVECTIVE HEAT FLUX INFLUENCE
C
      IF(NNHC.EQ.0)GO TO 6
      DO 5 K=1,NNHC
      NEL=LME(K)
      CALL NODSET(NEL,I,II,J,JJ,M,MM,N,NN)
      ISI1=ISIH(K)
      DO 5 IW=1,2
      CALL BCSIDE(NEL,ISI1,IW,SIDE,COSLX,COSLY,I,II,J,JJ,M,MM,N,NN)
      DO 5 JJ=1,NUM
      IF(NUMN.EQ.3.OR.NUMN.EQ.6)THEN
         IF(NUM.EQ.2)KN=NODES3(ISI1,JJ)
         IF(NUM.EQ.3)KN=NODES6(ISI1,JJ)
      ELSE
         IF(NUM.EQ.2)KN=NODES4(ISI1,JJ)
         IF(NUM.EQ.3)KN=NODES8(ISI1,JJ)
      ENDIF
      L=NODE(NEL,KN)
      DO 5 JJJ=1,NUM
      IF(NUMN.EQ.3.OR.NUMN.EQ.6)THEN
         IF(NUM.EQ.2)KK=NODES3(ISI1,JJJ)
         IF(NUM.EQ.3)KK=NODES6(ISI1,JJJ)
      ELSE
         IF(NUM.EQ.2)KK=NODES4(ISI1,JJJ)
         IF(NUM.EQ.3)KK=NODES8(ISI1,JJJ)
      ENDIF
      LL=NODE(NEL,KK)
      MASS=NS(KN)*NS(KK)
      IJPT=NFR(L)+LL-JMIN(L)
    5 C(IJPT)=C(IJPT)+MASS*H(K)*SIDE
    6 CONTINUE
C
C     FORM LHS
C
      DO 7 I=1,NFRMAX
    7 R(I)=AF*C(I)+P(I)/DT
C
      RETURN
      END
C ----------------------------------------------------------------------
      SUBROUTINE ASSEMB
      COMMON/VP/P(4001),C(4001),R(4001)
      COMMON/VA/COLD(201),CNEW(201)
      COMMON/VK/NNODE,NELEM,NTYPE,NTIME,TIME,NFRMAX
      COMMON/VC/Q(201),F(201),NTS(201),NQS(201),DX(201),DY(201)
      COMMON/VD/AF,AFM,DT,NNQS,NUMN,NUM,NGAUS,NNHC,NNST,NSTOP,KPRNT,IAXI
      COMMON/SK/NFR(201),JMIN(201),JMAX(201)
      DIMENSION B(201),FIXED(201)
C
      CALL CLEAR(B,NNODE)
C
C     CONSTRUCT RHS
C
      DO 4 L=1,NNODE
      B(L)=F(L)
      JSTR=JMIN(L)
      JEND=JMAX(L)
      DO 4 I1=JSTR,JEND
      NGI=NFR(L)+I1-JSTR
      RHS=(AFM*C(NGI)+P(NGI)/DT)*COLD(I1)
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
      SUBROUTINE MATAXI
      REAL NS,NX,NY,MASS
      COMMON/VP/P(4001),C(4001),R(4001)
      COMMON/VA/COLD(201),CNEW(201)
      COMMON/VB/X(201),Y(201),QQ(201)
      COMMON/VC/Q(201),F(201),NTS(201),NQS(201),DX(201),DY(201)
      COMMON/VE/NS(8),NX(8),NY(8),NXSI(8),NETA(8)
      COMMON/VF/POSG1(9),POSG2(9),W1(9),W2(9)
      COMMON/VD/AF,AFM,DT,NNQS,NUMN,NUM,NGAUS,NNHC,NNST,NSTOP,KPRNT,IAXI
      COMMON/VG/ISI(201),LEM(201),LME(201),ISIH(201)
      COMMON/VH/NODES4(4,2),NODES8(4,3),NODES3(3,2),NODES6(3,3)
      COMMON/VEL/VX(201),VY(201)
      COMMON/VM/H(201),TINF(201)
      COMMON/VK/NNODE,NELEM,NTYPE,NTIME,TIME,NFRMAX
      COMMON/VL/NODE(201,8)
      COMMON/SK/NFR(201),JMIN(201),JMAX(201)
C
      CALL CLEAR(R,NFRMAX)
      CALL CLEAR(C,NFRMAX)
      CALL CLEAR(P,NFRMAX)
C
C     CONSTRUCT AND STORE MASS MATRIX P
C
      DO 1 K=1,NELEM
      CALL NODSET(K,I,II,J,JJ,M,MM,N,NN)
      DO 1 IQ=1,NGAUS
      XSI=POSG1(IQ)
      ETA=POSG2(IQ)
      CALL SHAPE(K,XSI,ETA,I,J,M,N,II,JJ,MM,NN,DET,YETA,XETA,XXSI,YXSI)
      SUM=0.
      DO 2 KKKK=1,NUMN
      LR=NODE(K,KKKK)
    2 SUM=SUM+NS(KKKK)*Y(LR)
      DO 1 KK=1,NUMN
      L=NODE(K,KK)
      F(L)=F(L)+NS(KK)*QQ(K)*SUM*DET*W1(IQ)*W2(IQ)
      DO 1 KKK=1,NUMN
      LL=NODE(K,KKK)
      MASS=NS(KKK)*NS(KK)
      IJPT=NFR(L)+LL-JMIN(L)
    1 P(IJPT)=P(IJPT)+MASS*SUM*DET*W1(IQ)*W2(IQ)
C
C     CONSTRUCT AND STORE ADVECTION AND DIFFUSION MATRIX C
C
      DO 4 K=1,NELEM
      CALL NODSET(K,I,II,J,JJ,M,MM,N,NN)
      DO 4 IQ=1,NGAUS
      XSI=POSG1(IQ)
      ETA=POSG2(IQ)
      CALL SHAPE(K,XSI,ETA,I,J,M,N,II,JJ,MM,NN,DET,YETA,XETA,XXSI,YXSI)
      SUM=0.
      DO 9 KKKK=1,NUMN
      LR=NODE(K,KKKK)
    9 SUM=SUM+NS(KKKK)*Y(LR)
      DO 4 KK=1,NUMN
      L=NODE(K,KK)
      DO 4 KKK=1,NUMN
      LL=NODE(K,KKK)
C
C     ADVECTION TERM
C
       VELX=VX(L)*NS(KK)*NX(KKK)
        VELY=VY(L)*NS(KK)*NY(KKK)
         ADVEC=VELX+VELY
C
C     DIFFUSION TERM
C
       DIFFX=DX(K)*NX(KKK)*NX(KK)
        DIFFY=DY(K)*NY(KKK)*NY(KK)
         DIFF=DIFFX+DIFFY
C
      IJPT=NFR(L)+LL-JMIN(L)
    4 C(IJPT)=C(IJPT)+(DIFF+ADVEC)*SUM*DET*W1(IQ)*W2(IQ)
C
C     CONVECTIVE HEAT FLUX INFLUENCE
C
      IF(NNHC.EQ.0)GO TO 6
      DO 5 K=1,NNHC
      NEL=LME(K)
      CALL NODSET(NEL,I,II,J,JJ,M,MM,N,NN)
      ISI1=ISIH(K)
      DO 5 IW=1,2
      CALL BCSIDE(NEL,ISI1,IW,SIDE,COSLX,COSLY,I,II,J,JJ,M,MM,N,NN)
      SUM=0.
      DO 8 JR=1,NUM
      IF(NUMN.EQ.3.OR.NUMN.EQ.6)THEN
         IF(NUM.EQ.2)KR=NODES3(ISI1,JR)
         IF(NUM.EQ.3)KR=NODES6(ISI1,JR)
      ELSE
         IF(NUM.EQ.2)KR=NODES4(ISI1,JR)
         IF(NUM.EQ.3)KR=NODES8(ISI1,JR)
      ENDIF
      LR=NODE(NEL,KR)
    8 SUM=SUM+Y(LR)*NS(KR)
      DO 5 JJ=1,NUM
      IF(NUMN.EQ.3.OR.NUMN.EQ.6)THEN
         IF(NUM.EQ.2)KN=NODES3(ISI1,JJ)
         IF(NUM.EQ.3)KN=NODES6(ISI1,JJ)
      ELSE
         IF(NUM.EQ.2)KN=NODES4(ISI1,JJ)
         IF(NUM.EQ.3)KN=NODES8(ISI1,JJ)
      ENDIF
      L=NODE(NEL,KN)
      DO 5 JJJ=1,NUM
      IF(NUMN.EQ.3.OR.NUMN.EQ.6)THEN
         IF(NUM.EQ.2)KK=NODES3(ISI1,JJJ)
         IF(NUM.EQ.3)KK=NODES6(ISI1,JJJ)
      ELSE
         IF(NUM.EQ.2)KK=NODES4(ISI1,JJJ)
         IF(NUM.EQ.3)KK=NODES8(ISI1,JJJ)
      ENDIF
      LL=NODE(NEL,KK)
      MASS=NS(KN)*NS(KK)
      IJPT=NFR(L)+LL-JMIN(L)
    5 C(IJPT)=C(IJPT)+MASS*SUM*H(K)*SIDE
    6 CONTINUE
C
C     FORM LHS
C
      DO 7 I=1,NFRMAX
    7 R(I)=AF*C(I)+P(I)/DT
C
      RETURN
      END
C ----------------------------------------------------------------------
      BLOCK DATA
      COMMON/VH/NODES4(4,2),NODES8(4,3),NODES3(3,2),NODES6(3,3)
      COMMON/VI/POSG(2),POST(2),WQ(2),WT(2)
      COMMON/XY/X4(4),Y4(4),X8(8),Y8(8)
      DATA POSG/0.57735027,-0.57735027/
      DATA POST/0.7886751,0.1013249/
      DATA WQ/1.0,1.0/
      DATA WT/0.5,0.5/
      DATA NODES3/3,2,1,1,3,2/
      DATA NODES6/5,3,1,6,4,2,1,5,3/
      DATA NODES4/4,3,2,1,1,4,3,2/
      DATA NODES8/7,5,3,1,8,6,4,2,1,7,5,3/
      DATA X4/-1.0,1.0,1.0,-1.0/
      DATA Y4/-1.0,-1.0,1.0,1.0/
      DATA X8/-1.0,0.0,1.0,1.0,1.0,0.0,-1.0,-1.0/
      DATA Y8/-1.0,-1.0,-1.0,0.0,1.0,1.0,1.0,0.0/
      END
C ----------------------------------------------------------------------
      SUBROUTINE BNDCON
      REAL NS,NX,NY
      COMMON/VK/NNODE,NELEM,NTYPE,NTIME,TIME,NFRMAX
      COMMON/VC/Q(201),F(201),NTS(201),NQS(201),DX(201),DY(201)
      COMMON/VL/NODE(201,8)
      COMMON/VD/AF,AFM,DT,NNQS,NUMN,NUM,NGAUS,NNHC,NNST,NSTOP,KPRNT,IAXI
      COMMON/VE/NS(8),NX(8),NY(8),NXSI(8),NETA(8)
      COMMON/VG/ISI(201),LEM(201),LME(201),ISIH(201)
      COMMON/VM/H(201),TINF(201)
      COMMON/VH/NODES4(4,2),NODES8(4,3),NODES3(3,2),NODES6(3,3)
C
C     IW=1,2 DENOTES TWO GAUSS POINTS PER SIDE USING GAUSS QUADRATURE
C
   10 IF(NNQS.EQ.0)GO TO 35
      DO 3 K=1,NNQS
      NEL=LEM(K)
      CALL NODSET(NEL,I,II,J,JJ,M,MM,N,NN)
      ISI1=ISI(K)
      DO 3 IW=1,2
      CALL BCSIDE(NEL,ISI1,IW,SIDE,COSLX,COSLY,I,II,J,JJ,M,MM,N,NN)
      DO 3 L=1,NUM
      IF(NUMN.EQ.3.OR.NUMN.EQ.6)THEN
         IF(NUM.EQ.2)KK=NODES3(ISI1,L)
         IF(NUM.EQ.3)KK=NODES6(ISI1,L)
      ELSE
         IF(NUM.EQ.2)KK=NODES4(ISI1,L)
         IF(NUM.EQ.3)KK=NODES8(ISI1,L)
      ENDIF
      KKK=NODE(NEL,KK)
    3 F(KKK)=F(KKK)+NS(KK)*SIDE*Q(K)
C
   35 IF(NNHC.EQ.0)RETURN
      DO 9 K=1,NNHC
      NEL=LME(K)
      CALL NODSET(NEL,I,II,J,JJ,M,MM,N,NN)
      ISI1=ISIH(K)
      DO 9 IW=1,2
      CALL BCSIDE(NEL,ISI1,IW,SIDE,COSLX,COSLY,I,II,J,JJ,M,MM,N,NN)
      DO 9 L=1,NUM
      IF(NUMN.EQ.3.OR.NUMN.EQ.6)THEN
         IF(NUM.EQ.2)KK=NODES3(ISI1,L)
         IF(NUM.EQ.3)KK=NODES6(ISI1,L)
      ELSE
         IF(NUM.EQ.2)KK=NODES4(ISI1,L)
         IF(NUM.EQ.3)KK=NODES8(ISI1,L)
      ENDIF
      LL=NODE(NEL,KK)
    9 F(LL)=F(LL)+NS(KK)*SIDE*H(K)*TINF(K)
      RETURN
      END
C ----------------------------------------------------------------------
      SUBROUTINE BNDAXI
      REAL NS,NX,NY
      COMMON/VK/NNODE,NELEM,NTYPE,NTIME,TIME,NFRMAX
      COMMON/VC/Q(201),F(201),NTS(201),NQS(201),DX(201),DY(201)
      COMMON/VB/X(201),Y(201),QQ(201)
      COMMON/VL/NODE(201,8)
      COMMON/VD/AF,AFM,DT,NNQS,NUMN,NUM,NGAUS,NNHC,NNST,NSTOP,KPRNT,IAXI
      COMMON/VE/NS(8),NX(8),NY(8),NXSI(8),NETA(8)
      COMMON/VG/ISI(201),LEM(201),LME(201),ISIH(201)
      COMMON/VM/H(201),TINF(201)
      COMMON/VH/NODES4(4,2),NODES8(4,3),NODES3(3,2),NODES6(3,3)
C
C     IW=1,2 DENOTES TWO GAUSS POINTS PER SIDE USING GAUSS QUADRATURE
C
   10 IF(NNQS.EQ.0)GO TO 35
      DO 3 K=1,NNQS
      NEL=LEM(K)
      CALL NODSET(NEL,I,II,J,JJ,M,MM,N,NN)
      ISI1=ISI(K)
      DO 3 IW=1,2
      CALL BCSIDE(NEL,ISI1,IW,SIDE,COSLX,COSLY,I,II,J,JJ,M,MM,N,NN)
      SUM=0.
      DO 8 JR=1,NUM
      IF(NUMN.EQ.3.OR.NUMN.EQ.6)THEN
         IF(NUM.EQ.2)KR=NODES3(ISI1,JR)
         IF(NUM.EQ.3)KR=NODES6(ISI1,JR)
      ELSE
         IF(NUM.EQ.2)KR=NODES4(ISI1,JR)
         IF(NUM.EQ.3)KR=NODES8(ISI1,JR)
      ENDIF
      LR=NODE(NEL,KR)
    8 SUM=SUM+Y(LR)*NS(KR)
      DO 3 L=1,NUM
      IF(NUMN.EQ.3.OR.NUMN.EQ.6)THEN
         IF(NUM.EQ.2)KK=NODES3(ISI1,L)
         IF(NUM.EQ.3)KK=NODES6(ISI1,L)
      ELSE
         IF(NUM.EQ.2)KK=NODES4(ISI1,L)
         IF(NUM.EQ.3)KK=NODES8(ISI1,L)
      ENDIF
      KKK=NODE(NEL,KK)
    3 F(KKK)=F(KKK)+NS(KK)*SIDE*Q(K)*SUM
C
   35 IF(NNHC.EQ.0)RETURN
      DO 9 K=1,NNHC
      NEL=LME(K)
      CALL NODSET(NEL,I,II,J,JJ,M,MM,N,NN)
      ISI1=ISIH(K)
      DO 9 IW=1,2
      CALL BCSIDE(NEL,ISI1,IW,SIDE,COSLX,COSLY,I,II,J,JJ,M,MM,N,NN)
      SUM=0.
      DO 7 JR=1,NUM
      IF(NUMN.EQ.3.OR.NUMN.EQ.6)THEN
         IF(NUM.EQ.2)KR=NODES3(ISI1,JR)
         IF(NUM.EQ.3)KR=NODES6(ISI1,JR)
      ELSE
         IF(NUM.EQ.2)KR=NODES4(ISI1,JR)
         IF(NUM.EQ.3)KR=NODES8(ISI1,JR)
      ENDIF
      LR=NODE(NEL,KR)
    7 SUM=SUM+Y(LR)*NS(KR)
      DO 9 L=1,NUM
      IF(NUMN.EQ.3.OR.NUMN.EQ.6)THEN
         IF(NUM.EQ.2)KK=NODES3(ISI1,L)
         IF(NUM.EQ.3)KK=NODES6(ISI1,L)
      ELSE
         IF(NUM.EQ.2)KK=NODES4(ISI1,L)
         IF(NUM.EQ.3)KK=NODES8(ISI1,L)
      ENDIF
      LL=NODE(NEL,KK)
    9 F(LL)=F(LL)+NS(KK)*SIDE*H(K)*TINF(K)*SUM
      RETURN
      END
C ----------------------------------------------------------------------
      SUBROUTINE CLEAR(A,N)
      DIMENSION A(N)
      DO 1 I=1,N
    1 A(I)=0.
      RETURN
      END
C ----------------------------------------------------------------------
      SUBROUTINE BCSIDE(K,ISI1,IW,SIDE,COSLX,COSLY,I,II,J,JJ,M,MM,N,NN)
      COMMON/VD/AF,AFM,DT,NNQS,NUMN,NUM,NGAUS,NNHC,NNST,NSTOP,KPRNT,IAXI
      COMMON/VI/POSG(2),POST(2),WQ(2),WT(2)
      IF(NUMN.EQ.3.OR.NUMN.EQ.6)GOTO 1
      IF(NUMN.EQ.4.OR.NUMN.EQ.8)GOTO 9
C
C     3 OR 6 NODE TRIANGLES
C
    1 GO TO (2,3,4),ISI1
C
C     SIDE 1
C
    2 XSI= 0.
      ETA=POST(IW)
      CALL SHAPE(K,XSI,ETA,I,J,M,N,II,JJ,MM,NN,DET,YETA,XETA,XXSI,YXSI)
      R=-YETA
      S=XETA
      GO TO 49
C
C     SIDE 2
C
    3 XSI=POST(IW)
      ETA= 1.-XSI
      CALL SHAPE(K,XSI,ETA,I,J,M,N,II,JJ,MM,NN,DET,YETA,XETA,XXSI,YXSI)
      R=-YXSI
      S=XXSI
      GO TO 49
C
C     SIDE 3
C
    4 XSI=POST(IW)
      ETA= 0.
      CALL SHAPE(K,XSI,ETA,I,J,M,N,II,JJ,MM,NN,DET,YETA,XETA,XXSI,YXSI)
      R=YXSI
      S=-XXSI
   49 SIDEL=SQRT(R*R+S*S)
      COSLX=R/SIDEL
      COSLY=S/SIDEL
      SIDE=SIDEL*WT(IW)
      RETURN
C
C     4 OR 8 NODE QUADRILATERALS
C
    9 GO TO (10,20,30,40),ISI1
C
C     SIDE 1
C
   10 XSI=-1.
      ETA=POSG(IW)
      CALL SHAPE(K,XSI,ETA,I,J,M,N,II,JJ,MM,NN,DET,YETA,XETA,XXSI,YXSI)
      R=-YETA
      S=XETA
      GO TO 50
C
C     SIDE 2
C
   20 XSI=POSG(IW)
      ETA=1.
      CALL SHAPE(K,XSI,ETA,I,J,M,N,II,JJ,MM,NN,DET,YETA,XETA,XXSI,YXSI)
      R=-YXSI
      S=XXSI
      GO TO 50
C
C     SIDE 3
C
   30 XSI=1.
      ETA=POSG(IW)
      CALL SHAPE(K,XSI,ETA,I,J,M,N,II,JJ,MM,NN,DET,YETA,XETA,XXSI,YXSI)
      R=YETA
      S=-XETA
      GO TO 50
C
C     SIDE 4
C
   40 XSI=POSG(IW)
      ETA=-1.
      CALL SHAPE(K,XSI,ETA,I,J,M,N,II,JJ,MM,NN,DET,YETA,XETA,XXSI,YXSI)
      R=YXSI
      S=-XXSI
C
   50 SIDEL=SQRT(R*R+S*S)
      COSLX=R/SIDEL
      COSLY=S/SIDEL
      SIDE=SIDEL*WQ(IW)
      RETURN
      END
C ----------------------------------------------------------------------
      SUBROUTINE NODSET(K,I,II,J,JJ,M,MM,N,NN)
      COMMON/VL/NODE(201,8)
      COMMON/VD/AF,AFM,DT,NNQS,NUMN,NUM,NGAUS,NNHC,NNST,NSTOP,KPRNT,IAXI
      IF(NUMN.EQ.3)GO TO 1
      IF(NUMN.EQ.4)GO TO 2
      IF(NUMN.EQ.6)GO TO 3
      IF(NUMN.EQ.8)GO TO 4
C
C     3-NODE LINEAR TRIANGULAR ELEMENT
C
    1 I=NODE(K,1)
      J=NODE(K,2)
      M=NODE(K,3)
      RETURN
C
C     4-NODE LINEAR QUADRILATERAL ELEMENT
C
    2 I=NODE(K,1)
      J=NODE(K,2)
      M=NODE(K,3)
      N=NODE(K,4)
      RETURN
C
C     6-NODE QUADRATIC TRIANGULAR ELEMENT
C
    3 I=NODE(K,1)
      II=NODE(K,2)
      J=NODE(K,3)
      JJ=NODE(K,4)
      M=NODE(K,5)
      MM=NODE(K,6)
      RETURN
C
C     8-NODE QUADRATIC QUADRILATERAL ELEMENT
C
    4 I=NODE(K,1)
      II=NODE(K,2)
      J=NODE(K,3)
      JJ=NODE(K,4)
      M=NODE(K,5)
      MM=NODE(K,6)
      N=NODE(K,7)
      NN=NODE(K,8)
      RETURN
      END
C ----------------------------------------------------------------------
      SUBROUTINE SHAPE(K,XSI,ETA,I,J,M,N,II,JJ,MM,NN,DET,YETA,XETA,
     1XXSI,YXSI)
      REAL NS,NX,NY,NXSI,NETA
      COMMON/VB/X(201),Y(201),QQ(201)
      COMMON/VL/NODE(201,8)
      COMMON/VD/AF,AFM,DT,NNQS,NUMN,NUM,NGAUS,NNHC,NNST,NSTOP,KPRNT,IAXI
      COMMON/VE/NS(8),NX(8),NY(8),NXSI(8),NETA(8)
      COMMON/XY/X4(4),Y4(4),X8(8),Y8(8)
      D=0.25
      E=0.5
      IF(NUMN.EQ.3)GO TO 10
      IF(NUMN.EQ.4)GO TO 11
      IF(NUMN.EQ.6)GO TO 12
      IF(NUMN.EQ.8)GO TO 13
C
C     3-NODE LINEAR TRIANGULAR ELEMENT
C
   10 NS(1)=1-XSI-ETA
      NS(2)=XSI
      NS(3)=ETA
C
      NXSI(1)=-1
      NXSI(2)= 1
      NXSI(3)= 0
C
      NETA(1)=-1
      NETA(2)= 0
      NETA(3)= 1
C
      XXSI=NXSI(1)*X(I)+NXSI(2)*X(J)+NXSI(3)*X(M)
      XETA=NETA(1)*X(I)+NETA(2)*X(J)+NETA(3)*X(M)
      YXSI=NXSI(1)*Y(I)+NXSI(2)*Y(J)+NXSI(3)*Y(M)
      YETA=NETA(1)*Y(I)+NETA(2)*Y(J)+NETA(3)*Y(M)
      GO TO 7
C
C     4-NODE LINEAR QUADRILATERAL ELEMENT
C
   11 DO 2 L=1,NUMN
    2 NS(L)=D*(1+XSI*X4(L))*(1+ETA*Y4(L))
C
      NXSI(1)=-D*(1-ETA)
      NXSI(2)=D*(1-ETA)
      NXSI(3)=D*(1+ETA)
      NXSI(4)=-D*(1+ETA)
C
      NETA(1)=-D*(1-XSI)
      NETA(2)=-D*(1+XSI)
      NETA(3)=D*(1+XSI)
      NETA(4)=D*(1-XSI)
C
      XXSI=NXSI(1)*X(I)+NXSI(2)*X(J)+NXSI(3)*X(M)+NXSI(4)*X(N)
      XETA=NETA(1)*X(I)+NETA(2)*X(J)+NETA(3)*X(M)+NETA(4)*X(N)
      YXSI=NXSI(1)*Y(I)+NXSI(2)*Y(J)+NXSI(3)*Y(M)+NXSI(4)*Y(N)
      YETA=NETA(1)*Y(I)+NETA(2)*Y(J)+NETA(3)*Y(M)+NETA(4)*Y(N)
      GO TO 7
C
C     6-NODE QUADRATIC TRIANGULAR ELEMENT
C
   12 NS(1)=(1.-XSI-ETA)*(1.-2.*XSI-2.*ETA)
      NS(2)=4.*XSI*(1.-XSI-ETA)
      NS(3)=XSI*(2.*XSI-1.)
      NS(4)=4.*XSI*ETA
      NS(5)=ETA*(2.*ETA-1.)
      NS(6)=4.*ETA*(1.-XSI-ETA)
C
      NXSI(1)=4.*XSI+4.*ETA-3.
      NXSI(2)=4.-8.*XSI-4.*ETA
      NXSI(3)=4.*XSI-1.
      NXSI(4)=4.*ETA
      NXSI(5)=0.
      NXSI(6)=-4.*ETA
C
      NETA(1)=4.*ETA+4.*XSI-3.
      NETA(2)=-4.*XSI
      NETA(3)=0.
      NETA(4)=4.*XSI
      NETA(5)=4.*ETA-1.
      NETA(6)=4.-4.*XSI-8.*ETA
C
      XXSI=NXSI(1)*X(I)+NXSI(2)*X(II)+NXSI(3)*X(J)+NXSI(4)*X(JJ)+NXSI(5)
     1*X(M)+NXSI(6)*X(MM)
      XETA=NETA(1)*X(I)+NETA(2)*X(II)+NETA(3)*X(J)+NETA(4)*X(JJ)+NETA(5)
     1*X(M)+NETA(6)*X(MM)
      YXSI=NXSI(1)*Y(I)+NXSI(2)*Y(II)+NXSI(3)*Y(J)+NXSI(4)*Y(JJ)+NXSI(5)
     1*Y(M)+NXSI(6)*Y(MM)
      YETA=NETA(1)*Y(I)+NETA(2)*Y(II)+NETA(3)*Y(J)+NETA(4)*Y(JJ)+NETA(5)
     1*Y(M)+NETA(6)*Y(MM)
      GOTO 7
C
C     8-NODE QUADRATIC QUADRILATERAL ELEMENT
C
   13 DO 5 L=1,NUMN,2
    5 NS(L)=D*(1.+XSI*X8(L))*(1.+ETA*Y8(L))*(X8(L)*XSI+Y8(L)*ETA-1.)
      DO 6 L=2,NUMN,2
    6 NS(L)=E*X8(L)*X8(L)*(1.+X8(L)*XSI)*(1.-ETA*ETA)+E*Y8(L)*Y8(L)*(1.+
     1Y8(L)*ETA)*(1.-XSI*XSI)
C
      NXSI(1)=D*(ETA+2*XSI-2*XSI*ETA-ETA*ETA)
      NXSI(2)=-XSI+XSI*ETA
      NXSI(3)=D*(-ETA+2*XSI-2*XSI*ETA+ETA*ETA)
      NXSI(4)=E*(1.-ETA*ETA)
      NXSI(5)=D*(ETA+2*XSI+2*XSI*ETA+ETA*ETA)
      NXSI(6)=-XSI-XSI*ETA
      NXSI(7)=D*(-ETA+2*XSI+2*XSI*ETA-ETA*ETA)
      NXSI(8)=E*(-1.+ETA*ETA)
C
      NETA(1)=D*(XSI+2*ETA-2*XSI*ETA-XSI*XSI)
      NETA(2)=E*(-1.+XSI*XSI)
      NETA(3)=D*(-XSI+2*ETA-XSI*XSI+2*XSI*ETA)
      NETA(4)=-ETA-XSI*ETA
      NETA(5)=D*(XSI+2*ETA+XSI*XSI+2*XSI*ETA)
      NETA(6)=E*(1.-XSI*XSI)
      NETA(7)=D*(-XSI+2*ETA+XSI*XSI-2*XSI*ETA)
      NETA(8)=-ETA+XSI*ETA
C
      XXSI=NXSI(1)*X(I)+NXSI(2)*X(II)+NXSI(3)*X(J)+NXSI(4)*X(JJ)+NXSI(5)
     1*X(M)+NXSI(6)*X(MM)+NXSI(7)*X(N)+NXSI(8)*X(NN)
      XETA=NETA(1)*X(I)+NETA(2)*X(II)+NETA(3)*X(J)+NETA(4)*X(JJ)+NETA(5)
     1*X(M)+NETA(6)*X(MM)+NETA(7)*X(N)+NETA(8)*X(NN)
      YXSI=NXSI(1)*Y(I)+NXSI(2)*Y(II)+NXSI(3)*Y(J)+NXSI(4)*Y(JJ)+NXSI(5)
     1*Y(M)+NXSI(6)*Y(MM)+NXSI(7)*Y(N)+NXSI(8)*Y(NN)
      YETA=NETA(1)*Y(I)+NETA(2)*Y(II)+NETA(3)*Y(J)+NETA(4)*Y(JJ)+NETA(5)
     1*Y(M)+NETA(6)*Y(MM)+NETA(7)*Y(N)+NETA(8)*Y(NN)
C
    7 CONTINUE
      DET=XXSI*YETA-YXSI*XETA
      IF(DET.EQ.0.0)THEN
      WRITE(*,100)K
  100 FORMAT(2X,'THE DETERMINANT = 0.0 IN ELEMENT',I4)
      STOP
      ENDIF
      XSIX=YETA/DET
      XSIY=-XETA/DET
      ETAX=-YXSI/DET
      ETAY=XXSI/DET
      DO 8 KKK=1,NUMN
      NX(KKK)=NXSI(KKK)*XSIX+NETA(KKK)*ETAX
      NY(KKK)=NXSI(KKK)*XSIY+NETA(KKK)*ETAY
    8 CONTINUE
      RETURN
      END
C ----------------------------------------------------------------------
      SUBROUTINE GAUSS
      COMMON/VD/AF,AFM,DT,NNQS,NUMN,NUM,NGAUS,NNHC,NNST,NSTOP,KPRNT,IAXI
      COMMON/VF/POSG1(9),POSG2(9),W1(9),W2(9)
      IF(NUMN.EQ.3)GOTO 1
      IF(NUMN.EQ.4)GOTO 2
      IF(NUMN.EQ.6)GOTO 3
      IF(NUMN.EQ.8)GOTO 4
C
C     3 NODE LINEAR TRIANGLE (3 GAUSS POINTS)
C
    1 W=0.5
      POSG1(1)=W
      POSG1(2)=W
      POSG1(3)=0.
C
      POSG2(1)=0.
      POSG2(2)=W
      POSG2(3)=W
C
      W1(1)=1.0/6.
      W1(2)=1.0/6.
      W1(3)=1.0/6.
C
      W2(1)=1.0
      W2(2)=1.0
      W2(3)=1.0
      RETURN
C
C     4 NODE LINEAR QUADRILATERAL (4 GAUSS POINTS)
C
    2 W=SQRT(1.0/3.0)
      POSG1(1)=-W
      POSG1(2)=W
      POSG1(3)=W
      POSG1(4)=-W
C
      POSG2(1)=-W
      POSG2(2)=-W
      POSG2(3)=W
      POSG2(4)=W
C
      W1(1)=1.0
      W1(2)=1.0
      W1(3)=1.0
      W1(4)=1.0
C
      W2(1)=1.0
      W2(2)=1.0
      W2(3)=1.0
      W2(4)=1.0
      RETURN
C
C     6 NODE QUADRATIC TRIANGLE (7 GAUSS POINTS)
C
    3 POSG1(1)=0.
      POSG1(2)=1./2.
      POSG1(3)=1.
      POSG1(4)=1./2.
      POSG1(5)=0.
      POSG1(6)=0.
      POSG1(7)=1./3.
C
      POSG2(1)=0.
      POSG2(2)=0.
      POSG2(3)=0.
      POSG2(4)=1./2.
      POSG2(5)=1.
      POSG2(6)=1./2.
      POSG2(7)=1./3.
C
      W1(1)=1./40.
      W1(2)=1./15.
      W1(3)=1./40.
      W1(4)=1./15.
      W1(5)=1./40.
      W1(6)=1./15.
      W1(7)=9./40.
C
      W2(1)=1.0
      W2(2)=1.0
      W2(3)=1.0
      W2(4)=1.0
      W2(5)=1.0
      W2(6)=1.0
      W2(7)=1.0
      RETURN
C
C     8 NODE QUADRATIC QUADRILATERAL (9 GAUSS POINTS)
C
    4 W=SQRT(0.60)
      POSG1(1)=-W
      POSG1(2)=0.
      POSG1(3)=W
      POSG1(4)=W
      POSG1(5)=W
      POSG1(6)=0.
      POSG1(7)=-W
      POSG1(8)=-W
      POSG1(9)=0.
C
      POSG2(1)=-W
      POSG2(2)=-W
      POSG2(3)=-W
      POSG2(4)=0.
      POSG2(5)=W
      POSG2(6)=W
      POSG2(7)=W
      POSG2(8)=0.
      POSG2(9)=0.
C
      W1(1)=5.0/9.0
      W1(2)=8.0/9.0
      W1(3)=5.0/9.0
      W1(4)=5.0/9.0
      W1(5)=5.0/9.0
      W1(6)=8.0/9.0
      W1(7)=5.0/9.0
      W1(8)=5.0/9.0
      W1(9)=8.0/9.0
C
      W2(1)=5.0/9.0
      W2(2)=5.0/9.0
      W2(3)=5.0/9.0
      W2(4)=8.0/9.0
      W2(5)=5.0/9.0
      W2(6)=5.0/9.0
      W2(7)=5.0/9.0
      W2(8)=8.0/9.0
      W2(9)=8.0/9.0
      RETURN
      END
C ----------------------------------------------------------------------
      SUBROUTINE PRINT
      COMMON/VA/COLD(201),CNEW(201)
      COMMON/VK/NNODE,NELEM,NTYPE,NTIME,TIME,NFRMAX
      COMMON/VEL/VX(201),VY(201)
      IF(NTYPE.EQ.1)THEN
      WRITE(*,403)
      ELSE
      WRITE(*,402)NTIME,TIME
      ENDIF
      WRITE(*,400)
      DO 1 I=1,NNODE
    1 WRITE(*,100)I,CNEW(I),VX(I),VY(I)
  100 FORMAT(5X,I3,5X,3(F10.3))
  402 FORMAT(/,5X,'TIME STEPS= ',I6,2X,'TIME= ',F10.3)
  403 FORMAT(/,5X,'STEADY STATE')
  400 FORMAT(/,5X,'NODE',5X,'CNEW',8X,'VX',6X,'VY')
      RETURN
      END
C ----------------------------------------------------------------------
      SUBROUTINE RESID
      COMMON/VD/AF,AFM,DT,NNQS,NUMN,NUM,NGAUS,NNHC,NNST,NSTOP,KPRNT,IAXI
      COMMON/VK/NNODE,NELEM,NTYPE,NTIME,TIME,NFRMAX
      COMMON/VA/COLD(201),CNEW(201)
      INTEGER ERRPRN
      DATA ERRPRN/5/
C
      LOCE=0
      RESMAX=0.
      ERR=0.0001
C
      DO 1 I=1,NNODE
      R=ABS(CNEW(I)-COLD(I))
      IF(R.GT.RESMAX)THEN
         RESMAX=R
        LOCE=I
       ENDIF
    1 CONTINUE
C
      IF(RESMAX.GT.ERR)THEN
      IF(NTIME.GE.NSTOP)THEN
       WRITE(*,101)NTIME,RESMAX
        CALL PRINT
         CALL PRINTF
        STOP
       END IF
      IF(MOD(NTIME,ERRPRN).EQ.0)THEN
       WRITE(*,102)NTIME,RESMAX,LOCE
      ENDIF
      RETURN
       ELSE
        WRITE(*,100)NTIME
         CALL PRINT
          CALL PRINTF
         STOP
        END IF
  100 FORMAT(/,2X,'PROGRAM HAS CONVERGED IN',I6,' STEPS')
  102 FORMAT(/,2X,'TIME STEP = ',I6,' RESIDUAL = ',1PE10.3,2X,' NODE = '
     1 ,I4)
  101 FORMAT(/,2X,'PROGRAM HAS NOT CONVERGED IN ',I6,' STEPS',/,3X,'TH
     1E RESIDUAL VALUE IS ',1PE10.3)
      END
C ----------------------------------------------------------------------
      SUBROUTINE PRINTF
      COMMON/VK/NNODE,NELEM,NTYPE,NTIME,TIME,NFRMAX
      COMMON/VD/AF,AFM,DT,NNQS,NUMN,NUM,NGAUS,NNHC,NNST,NSTOP,KPRNT,IAXI
      COMMON/VC/Q(201),F(201),NTS(201),NQS(201),DX(201),DY(201)
      COMMON/VA/COLD(201),CNEW(201)
      COMMON/VB/X(201),Y(201),QQ(201)
      COMMON/VEL/VX(201),VY(201)
      COMMON/VL/NODE(201,8)
      COMMON/VT/MTYPE,NUMDIM,TO,RHO,CP
      NDUMMY=0
      WRITE(9,'(8I6)')MTYPE,NUMDIM,NNODE,NELEM,NUMN,NSTOP,KPRNT,NDUMMY
      WRITE(9,'(I4,2X,5(F8.4,1X),2X,I5)')NTYPE,DT,AF,TO,RHO,CP,IAXI
      DO 1 I=1,NNODE
    1 WRITE(9,101)I,X(I),Y(I)
      DO 2 I=1,NELEM
    2 WRITE(9,100)I,QQ(I),DX(I),DY(I),(NODE(I,J),J=1,NUMN)
      WRITE(9,103)NTIME
      DO 3 I=1,NNODE
    3 WRITE(9,102)I,CNEW(I),VX(I),VY(I)
      CALL TECPLT
  100 FORMAT(5X,I4,2X,3(F10.3,1X),8I4)
  101 FORMAT(5X,I4,2X,2(F10.3,1X))
  102 FORMAT(5X,I4,2X,3(F10.3,1X))
  103 FORMAT(/,2X,'NO. OF TIME STEPS = ',I5)
      RETURN
      END
C
	SUBROUTINE TECPLT
      COMMON/VK/NNODE,NELEM,NTYPE,NTIME,TIME,NFRMAX
      COMMON/VD/AF,AFM,DT,NNQS,NUMN,NUM,NGAUS,NNHC,NNST,NSTOP,KPRNT,IAXI
      COMMON/VC/Q(201),F(201),NTS(201),NQS(201),DX(201),DY(201)
      COMMON/VA/COLD(201),CNEW(201)
      COMMON/VB/X(201),Y(201),QQ(201)
      COMMON/VEL/VX(201),VY(201)
      COMMON/VL/NODE(201,8)
C
C     SAVE SOLUTION IN TECPLOT FILE FORMAT
C
      OPEN(10,FILE='tecplt.plt',STATUS='UNKNOWN')
        WRITE(10,*)'TITLE = "2-D FEM RESULTS"'
        WRITE(10,*)'VARIABLES = X,Y,CNEW,VX,VY'
        WRITE(10,25)NNODE,NELEM
        DO 30,I=1,NNODE
          WRITE(10,35)X(I),Y(I),CNEW(I),VX(I),VY(I)
   30   CONTINUE
        DO 40,I=1,NELEM
          WRITE(10,45) NODE(I,1),NODE(I,2),NODE(I,3),NODE(I,4)
   40   CONTINUE
      CLOSE(10)
   25 FORMAT('ZONE I=',I5,',J=',I5,',F=FEPOINT')
   35 FORMAT(1X,5(E10.3,1X))
   45 FORMAT(4I10)
      RETURN
      END
C ----------------------------------------------------------------------
      SUBROUTINE SETPC
      CHARACTER INFIL*12,OUTFIL*12
      DATA INFIL/'            '/,OUTFIL/'            '/
      WRITE(*,'(/1X,A\)')' ENTER NAME FOR INPUT FILE: '
      READ(*,'(BN,A)')INFIL
      OPEN(8,FILE=INFIL)
C
      WRITE(*,'(/1X,A\)')' ENTER NAME FOR OUTPUT FILE: '
      READ(*,'(BN,A)')OUTFIL
      OPEN(9,FILE=OUTFIL,STATUS='NEW')
      WRITE(9,'(A)')OUTFIL
      RETURN
      END
C ----------------------------------------------------------------------
      SUBROUTINE SKYLIN(B,VAR,NB,FIXED,NBOUND)
      COMMON/VP/P(4001),C(4001),R(4001)
      COMMON/VK/NNODE,NELEM,NTYPE,NTIME,TIME,NFRMAX
      COMMON/SK/NFR(201),JMIN(201),JMAX(201)
      DIMENSION AG(4001),B(201),VAR(201),NB(201),FIXED(201)
      NN1=NNODE-1
      DO 1 I=1,NFRMAX
    1 AG(I)=R(I)
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
