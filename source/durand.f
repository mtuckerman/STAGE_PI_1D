c
      SUBROUTINE DURAND(SEED,N,X)
      REAL*8 X(*)
      REAL*8 SEED
      REAL*8 S,T,XM,XMP,SCALEU,SCALED,ONEMEP,ONEPEP
      DATA XM/16807.D0/
      DATA XMP/16807.00000 78263 69255 78117 37060 54687 5D0/
      DATA SCALEU/2147483648.D0/
      DATA SCALED /.00000 00004 65661 28730 77392 57812 5D0/
      DATA ONEPEP/1.00000 00004 65661 28730 77392 57812 5D0/
      DATA ONEMEP/ .99999 99995 34338 71269 22607 42187 5D0/
      DATA X123/64 41594.D0/
      REAL*8 S1,T1
      REAL*8 A,XR,XADD,XSUB
      data xr/.5d0/
      data xadd/4 50359 96273 70497.d0/
      data xsub/4 50359 96273 70497.d0/
      IF(N.LT.0)GO TO 100
      IF(SEED.LT.1.D0)GO TO 100
      IF(SEED.GT.2147483646.D0)GO TO 100
      T=SEED*SCALED
      DO 1 J=1,N
      S=INT(T*XMP)
      T=T*XM
      T=T-S*ONEMEP
    1 X(J)=T*ONEPEP
      SEED=T*SCALEU
      RETURN
  100 CONTINUE
c      IF(N.LT.0)CALL EMON(1,2)
c      IF(SEED.LT.1.D0)CALL EMON(73,1)
c      IF(SEED.GT.2147483646.D0)CALL EMON(73,1)
c      CALL EMON(99)
      RETURN
      END