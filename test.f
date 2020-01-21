C----|--1---------2---------3---------4---------5---------6---------7-|
     
      PROGRAM TEST
C
      IMPLICIT NONE
C
      REAL X,Y
      INTEGER I
    
C
      DO I=1,10
        X=1.0*I
        CALL KWADRAAT(X,Y)
        WRITE(*,1000)X,Y
      END DO
 
 1000 FORMAT("BLABLABLA",F7.2,F7.2)

      END
