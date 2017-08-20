 
      SUBROUTINE mgf2eff(TJUL,X,Y,Z,LONG,DXT,DYT,DZT)

      IMPLICIT NONE


      REAL*8        TJUL
      REAL*8        X,Y,Z
      INTEGER*4     NEWCMB
      REAL*8        P(3,3)
      REAL*8        Q(3,3)
      REAL*8        A(3,3)
      INTEGER*4     IR
      INTEGER*2     I,J,K
      REAL*8        DXT, DYT, DZT
      REAL*8        LONG
      REAL*8        DTRA4(3)

Cf2py intent(in) TJUL
Cf2py intent(in) X
Cf2py intent(in) Y
Cf2py intent(in) Z
Cf2py intent(in) LONG

Cf2py intent(out) DXT
Cf2py intent(out) DYT
Cf2py intent(out) DZT

      TJUL=15198.9166667D0 - (1.D0/48.D0) !11789.5D0!11550.0D0 !0833D0
      NEWCMB=1

      CALL COOT

C  CONVERT EPOCH TO 2000.0 (and allow for possible FDS day offset)      
!       D = DAY - 0.5D0 - d_fd_day_offset                    ! mjhw
!            TJUL = TJUL - 0.5D0
C  DISTINGUISH DAY IN MJD1950 FROM MJD2000 BY SIGN IN XLSUNA:           
!           TJUL = TJUL-18262.D0                                
      CALL WOBBLE( TJUL, NEWCMB, P, IR )
C     See if it was successful - check flag 
C     (IR = RETURN CODE. 0=OK, 1=TOO EARLY, 2=TOO LATE)

      DXT = X                                                          
      DYT = Y                                                          
      DZT = Z

      if (LONG<0) THEN
          LONG=(360+LONG)*3.14159265359D0/180.
      END IF


C Use nominal longitude of the spacecraft                        

        A(1,1) = COS( LONG )                                           
        A(1,2) = SIN( LONG )                                           
        A(1,3) = 0.0D0                                                  
        A(2,1) = - SIN( LONG )                                         
        A(2,2) = COS( LONG )
        A(2,3) = 0.0D0                                                  
        A(3,1) = 0.0D0                                                  
        A(3,2) = 0.0D0                                                  
        A(3,3) = 1.0D0                                                  

                                                                      
        DO 17 I = 1, 3                                                 
          DO 16 J = 1, 3                                              
              Q(I,J) = 0.0D0                                            
              DO 15 K = 1, 3                                           
                Q(I,J) = Q(I,J) + A(I,K) * P(K,J)                     
15            CONTINUE                                                 
16         CONTINUE                                                    
17      CONTINUE                                                       

C Transform the coordinates of the spacecraft            
C from the mean geocentric frame into           
C the Earth Fixed Frame (EFF)                               
C                                                                
        DTRA4(1) = DXT                                 
        DTRA4(2) = DYT                                  
        DTRA4(3) = DZT                                 
C                                                   
        DXT =                                                               
     *    Q(1,1)*DTRA4(1) + Q(1,2)*DTRA4(2) + Q(1,3)*DTRA4(3)         
        DYT =                                                               
     *    Q(2,1)*DTRA4(1) + Q(2,2)*DTRA4(2) + Q(2,3)*DTRA4(3)         
        DZT =                                                               
     *    Q(3,1)*DTRA4(1) + Q(3,2)*DTRA4(2) + Q(3,3)*DTRA4(3)         

      END