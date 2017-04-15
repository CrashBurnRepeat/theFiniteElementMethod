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