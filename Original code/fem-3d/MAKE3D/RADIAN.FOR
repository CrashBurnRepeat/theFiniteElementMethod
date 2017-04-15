      SUBROUTINE RADIAN
      COMMON/VZ/X(81),Y(81),Z(81),QQ(81),AREA(81),VOLUME(81)
      COMMON/VA/COLD(81),CNEW(81),KNODE(81),FIXED(81),E(81),AB(81)
      COMMON/VD/AF,AFM,DT,NNQS,NUMN,NUM,NNHC,NNST,NSTOP,KPRNT,IFLAG
      COMMON/VK/NNODE,NELEM,NTYPE,RHOCP,NTIME,TIME,NNRA,TO,TAMB
      COMMON/VL/NODE(81,4),NODEF(81,3),NODEC(81,3),NOR1(81,3),NOR2(81,3)
      COMMON/VC/Q(81),F(81),NTS(81),NQS(81),DX(81),DY(81),DZ(81)
      COMMON/VM/H(81),TINF(81),TH(81),FORM(81),HR(81)
C
C     RADIATION NORMAL TO BOUNDARY
C
C     UNITS MUST BE IN SI (METERS,WATTS,JOULES,AND K TEMPERATURE)
C
      SIG=5.669E-08
C
      DO 10 K=1,NNRA
      IF(HR(K).NE.0.0)THEN
      RFLAG=1.0
      GOTO 13
      ELSE
      RFLAG=0.0
      ENDIF
   10 CONTINUE
C
C     USES E, A, AND FNM
C
   13 IF(RFLAG.EQ.0.0)THEN
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
      ELSE
C
C     USES HR(K) FOR RADIATION FLUX COEFFICIENT
C
   11 DO 12 K=1,NNRA
      I = NOR1(K,1)
      J = NOR1(K,2)
      IF(NUMN.EQ.4)GOTO 101
      XL = SQRT((X(I)-X(J))**2 + (Y(I)-Y(J))**2) * TH(K)
      TAVEG = 0.5*(COLD(I) + COLD(J))
      TV1=TAVEG
      IF(NOR2(K,1).EQ.0)THEN
      TV2=TAMB
      QR = (TV1-TV2)*XL*HR(K)*0.5
      F(I)=F(I)-QR
      F(J)=F(J)-QR
      ELSE
      L=NOR2(K,1)
      LL=NOR2(K,2)
      T1 = 0.5*(COLD(L) + COLD(LL))
      TV2=T1
      QR = (TV1-TV2)*XL*HR(K)*0.5
      F(I)=F(I)-QR
      F(J)=F(J)-QR
      ENDIF
      GOTO 12
C
  101 CONTINUE
      N= NOR1(K,3 )
      CALL SIDE(ASIDE,I,J,N)
      A1 = ASIDE
      TAVEG = (COLD(I) + COLD(J) + COLD(N)) / 3.
      TV1=TAVEG
      IF(NOR2(K,1).NE.0)GOTO 221
      TV2=TAMB
      QR = (TV1-TV2)*HR(K)*A1/3.
      F(I)=F(I)-QR
      F(J)=F(J)-QR
      F(N)=F(N)-QR
      GOTO 12
  221 II = NOR2(K,1)
      JJ = NOR2(K,2)
      NN = NOR2(K,3)
      T1 = (COLD(II) + COLD(JJ) + COLD(NN)) / 3.
      TV2=T1
      QR = (TV1-TV2)*A1*HR(K)/3.
      F(I)=F(I)-QR
      F(J)=F(J)-QR
      F(N)=F(N)-QR
   12 CONTINUE
C
      ENDIF
C
      RETURN
      END
