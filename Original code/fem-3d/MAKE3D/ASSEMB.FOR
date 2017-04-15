      SUBROUTINE ASSEMB
      COMMON/VA/COLD(81),CNEW(81),KNODE(81),FIXED(81),E(81),AB(81)
      COMMON/VV/A(81,81),P(81,81),R(81,81),B(81)
      COMMON/VZ/X(81),Y(81),Z(81),QQ(81),AREA(81),VOLUME(81)
      COMMON/VK/NNODE,NELEM,NTYPE,RHOCP,NTIME,TIME,NNRA,TO,TAMB
      COMMON/VL/NODE(81,4),NODEF(81,3),NODEC(81,3),NOR1(81,3),NOR2(81,3)
      COMMON/VC/Q(81),F(81),NTS(81),NQS(81),DX(81),DY(81),DZ(81)
      COMMON/VD/AF,AFM,DT,NNQS,NUMN,NUM,NNHC,NNST,NSTOP,KPRNT,IFLAG
      COMMON/VM/H(81),TINF(81),TH(81),FORM(81),HR(81)
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
C  34 format(2x,12f6.2)
C 100 format(/,2x,'values of p')
C
C     write(*,103)
C     do 35 i=1,nnode
C  35 write(*,34)(a(i,j),j=1,nnode)
C 103 format(/,2x,'values of a before convection')
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
C 101 format(/,2x,'values of a after convection')
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
C 102 format(/,2x,'values of r')
C
      RETURN
      END
