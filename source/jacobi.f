C################jacobi#############################################
      SUBROUTINE JACOBI(NM,N,A,W,EIVR,IERR)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(NM,N),W(N),EIVR(NM,N)
C...................................................................
C     RSJC -- use jacobi method to emulate eispack RS
C     THIS ROUTINE USES A VARIABLE THRESHOLD JACOBI METHOD
C     IT GIVES VERY GOOD EIGENVALUES AND EIGENVECTORS
C     THE ROUTINE IS MUCH FASTER THAN THE OLD HDIAG ROUTINE WRITTEN
C     AT M.I.T. THAT USES THE JACOBI METHOD BUT NOT THE VARIABLE
C     THRESHOLD TECHNIQUE THAT IS APPLIED HERE
C...................................................................
C-----parameters---------------
      SQ2INV = 1.D0/SQRT(2.D0)
      T1     = 1.D-14
      T2     = 1.D-14
      IEGEN  = 0
C------------------------------
      IERR   = 0
      AVGF   = DBLE(N*(N-1))*0.55
      IF(N.LT.1) THEN
        GOTO 599
      ELSE IF(N.LT.2) THEN
        EIVR(1,1)=1.0
        GOTO 599
      END IF
      IF(IEGEN.EQ.0) THEN
        DO 101 J=1,N
          DO 100 I=1,N
100         EIVR(I,J)=0.0
101       EIVR(J,J)=1.0
      END IF
C...................................................................
C           FIND THE ABSOLUTELY LARGEST ELEMENT OF A
C...................................................................
      ATOP=0.0
      DO 111 I=1,N
        DO 111 J=1,N
  111     IF(ATOP.LT. ABS(A(I,J)))ATOP= ABS(A(I,J))
      IF(ATOP.LE.0.00) THEN
        IERR = 1
        GOTO 599
      END IF
C...............................................
C     CALCULATE THE STOPPING CRITERION -- DSTOP.
C...............................................
      D = 0.0
      DO 114 JJ=2,N
        DO 114 II=2,JJ
114       D = D + A(II-1,JJ)**2
      DSTOP=T1*D
      IF(DSQRT(D/ATOP).LT.T1) GOTO 599
C....................................
C     CALCULATE THE THRESHOLD, THRSH.
C....................................
      THRSH= SQRT(D/AVGF)
C...................
C     START A SWEEP.
C...................
115   IFLAG=0
      DO 130 JCOL=2,N
        JCOL1=JCOL-1
        DO 130 IROW=1,JCOL1
          AIJ=A(IROW,JCOL)
C.................................................
C     COMPARE THE OFF-DIAGONAL ELEMENT WITH THRSH.
C.................................................
          IF( ABS(AIJ).LE.THRSH) GOTO 130
          AII=A(IROW,IROW)
          AJJ=A(JCOL,JCOL)
          S=AJJ-AII
C...................................................................
C     THE CHOSEN ROTATION IS LESS THAN THE ROUNDING.  DO NOT ROTATE.
C...................................................................
          IF (ABS(AIJ).LT.T2*ABS(S)) GOTO 130
          IFLAG=1
          IF(T2*ABS(AIJ).GE.ABS(S))THEN
C....................................................................
C     ROTATION IS VERY CLOSE TO 45 DEGREES, SIN AND COS = 1/(ROOT 2).
C....................................................................
            S = SQ2INV
            C = S
          ELSE
C..............................................
C     ROTATION IS NOT VERY CLOSE TO 45 DEGREES.
C..............................................
            T = AIJ/S
            U = 0.25D0/ SQRT(0.25D0+T*T)
            C =  SQRT(0.5D0+U)
            S = 2.D0*T*U/C
          END IF
C........................................
C     CALCULATE NEW ELEMENTS OF MATRIX A.
C........................................
          DO 121 I=1,IROW
            T         = A(I,IROW)
            U         = A(I,JCOL)
            A(I,IROW) = C*T-S*U
121         A(I,JCOL) = S*T+C*U
          I2 = IROW+2
          IF (I2.LE.JCOL) THEN
            DO 122 I=I2,JCOL
              T           = A(I-1,JCOL)
              U           = A(IROW,I-1)
              A(I-1,JCOL) = S*U+C*T
122           A(IROW,I-1) = C*U-S*T
          END IF
123       A(JCOL,JCOL) = S*AIJ+C*AJJ
          A(IROW,IROW) = C*A(IROW,IROW)-S*(C*AIJ-S*AJJ)
          DO 124 J=JCOL,N
            T         = A(IROW,J)
            U         = A(JCOL,J)
            A(IROW,J) = C*T-S*U
124         A(JCOL,J) = S*T+C*U
C................................................................
C     ROTATION COMPLETED. SEE IF EIGENVECTORS ARE WANTED BY USER.
C................................................................
            IF(IEGEN.EQ.0) THEN
              DO 125 I=1,N
                T=EIVR(I,IROW)
                EIVR(I,IROW)=C*T-EIVR(I,JCOL)*S
125             EIVR(I,JCOL)=S*T+EIVR(I,JCOL)*C
            END IF
C.....................................................
C     CALCULATE THE NEW NORM D AND COMPARE WITH DSTOP.
C.....................................................
            S=AIJ
            D=D-S*S
            IF(D.LT.DSTOP) THEN
C............................................................
C     RECALCULATE DSTOP AND THRSH TO DISGARD ROUNDING ERRORS.
C............................................................
              D=0.
              DO 128 JJ=2,N
                DO 128 II=2,JJ
128               D=D+A(II-1,JJ)**2
              DSTOP=T1*D
            END IF
            THRSH= SQRT(D/AVGF)
130   CONTINUE
      IF(IFLAG.NE.0) GOTO 115
C............................
C     FILL EIGENVALUE VECTOR.
C............................
  599 DO 600 I=1,N
  600   W(I) = A(I,I)
      RETURN
      END
