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
