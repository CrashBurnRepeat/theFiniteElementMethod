      SUBROUTINE BNDCON
      COMMON/VA/COLD(81),CNEW(81),KNODE(81),FIXED(81),E(81),AB(81)
      COMMON/VV/A(81,81),P(81,81),R(81,81),B(81)
      COMMON/VZ/X(81),Y(81),Z(81),QQ(81),AREA(81),VOLUME(81)
      COMMON/VK/NNODE,NELEM,NTYPE,RHOCP,NTIME,TIME,NNRA,TO,TAMB
      COMMON/VL/NODE(81,4),NODEF(81,3),NODEC(81,3),NOR1(81,3),NOR2(81,3)
      COMMON/VC/Q(81),F(81),NTS(81),NQS(81),DX(81),DY(81),DZ(81)
      COMMON/VD/AF,AFM,DT,NNQS,NUMN,NUM,NNHC,NNST,NSTOP,KPRNT,IFLAG
      COMMON/VM/H(81),TINF(81),TH(81),FORM(81),HR(81)
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
