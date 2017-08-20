C                                                                       
C...................................................................... 
C                                                                       
C SUBROUTINE  320230D1  WOBBLE   VERSION 02    890712                   
C OAD              FORTRAN                                              
C                                                                       
C TRANSFORMATION TO EARTH FIXED FRAME OF REFERENCE                      
C                                                                       
C CALLING SEQUENCE -                                                    
C CALL WOBBLE(TJUL,NEWCMB,P,IR)                                         
C                                                                       
C EXPLANATION -                                                         
CP  GIVES TRANSFORMATION MATRIX P() FROM INERTIAL MEAN SYSTEM OF        
CP  DATE TO TRUE SYSTEM OF DATE OR EARTH-FIXED ROTATING SYSTEM          
CP  BY ROTATION (OPTIONAL) AND NUTATION.                                
CP  USES AN ANALYTICAL MODEL WITH 6.D-7 RAD ACCURACY.                   
C                                                                       
CI  TJUL = MODIFIED JULIAN DAY                                          
CI  NEWCMB =-1 => TO TRUE SYSTEM OF DATE                                
CI  NEWCMB = 0 => INTERMEDIATE TRANSFORMATION TO EARTH-FIXED ROTA-      
CI   TING SYSTEM BUT WITHOUT NEWCOMB'S ANGLE. (NOT TRUE-OF-DATE)        
CI  NEWCMB = 1 => TO EARTH-FIXED ROTATING SYSTEM                        
C                                                                       
CO  P(3,3) = ROTATION MATRIX SUCH THAT                                  
CO  THE TRANSFORMATION IS    R(TRANSF) = P * R(INERT)                   
CO  IR = RETURN CODE. 0=OK, 1=TOO EARLY, 2=TOO LATE.                    
C                                                                       
C SPECIAL FEATURES -                                                    
C DOCUMENTATION ADDED BY R.A.BOWEN TO UNCHANGED OAD SOURCE              
C                                                                       
C CHANGES -                                                             
C     01 NEW                                                            
C     02  890712  NEW IDENTIFICATION NUMBER 023                         
C                                                                       
C...................................................................... 
C                                                                       
      SUBROUTINE WOBBLE(TJUL,NEWCMB,P,IR)                               
CP  GIVES TRANSFORMATION MATRIX P() FROM INERTIAL MEAN SYSTEM OF        
CP  DATE TO TRUE SYSTEM OF DATE OR EARTH-FIXED ROTATING SYSTEM          
CP  BY ROTATION (OPTIONAL) AND NUTATION.                                
CP  USES AN ANALYTICAL MODEL WITH 6.D-7 RAD ACCURACY.                   
C                                                                       
CI  TJUL = MODIFIED JULIAN DAY                                          
CI  NEWCMB =-1 => TO TRUE SYSTEM OF DATE                                
CI  NEWCMB = 0 => INTERMEDIATE TRANSFORMATION TO EARTH-FIXED ROTA-      
CI   TING SYSTEM BUT WITHOUT NEWCOMB'S ANGLE. (NOT TRUE-OF-DATE)        
CI  NEWCMB = 1 => TO EARTH-FIXED ROTATING SYSTEM                        
C                                                                       
CO  P(3,3) = ROTATION MATRIX SUCH THAT                                  
CO  THE TRANSFORMATION IS    R(TRANSF) = P * R(INERT)                   
CO  IR = RETURN CODE. 0=OK, 1=TOO EARLY, 2=TOO LATE.                    
C                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      INTEGER IR,NEWCMB
      DIMENSION P(9)                                                    
      COMMON /CBASIC/ PI,TWOPI,PIH,DEG,RAD                              
      COMMON/CDYNAE/ STD50,OMT50,OMQ50,OMROT,STD50R,OMT50R,OMQ50R       
!       PRINT*,"PI= "
!       PRINT*, PI
!       PRINT*,"STD50= "
!       PRINT*, STD50
      IR=1                                                              
      IF(TJUL.LT.0.D0) RETURN                                           
      IR=2                                                              
      IF(TJUL.GE.1.D5) RETURN                                           
      IR=0                                                              
C     NUTATION CORRECTION                                               
      C1=-0.0529539D0*TJUL                                              
      C2=(12.1128D0+C1)*RAD                                             
      C22=2.D0*C2                                                       
      SOM=DSIN(C2)                                                      
      C3=0.9856473D0*TJUL                                               
      C42 = 2.D0*(280.0812+C3)*RAD                                      
      SOM2=DSIN(C42)                                                    
      C5=13.176396D0*TJUL                                               
      C62 = 2.D0*(64.3824D0+C5)*RAD                                     
      SOM4= DSIN(C22) - DSIN(C62)                                       
      DNU=-33.3D-06*SOM+0.4D-06*SOM4-2.5D-06*SOM2                       
      DEPS=44.7D-06*DCOS(C2) + 2.7D-06*DCOS(C42)+                       
     & 0.4D-06*(DCOS(C62)-DCOS(C22))                                    
C     I-    ROTATION OF -DEPS ABOUT THE X-AXIS                          
      CDE=DCOS(DEPS)                                                    
      SDE=DSIN(DEPS)                                                    
C                                                                       
C     II-   ROTATION OF  DNU  ABOUT THE Y-AXIS                          
      CNU=DCOS(DNU)                                                     
      SNU=DSIN(DNU)                                                     
C                                                                       
      IF(NEWCMB) 10,40,20                                               
C                                                                       
C  ROTATION OF -DMU => ST  ABOUT THE Z-AXIS                             
  10  ST = 76.7D-6*SOM - 0.9D-6*SOM4 + 5.7D-6*SOM2                      
      GOTO 30                                                           
C                                                                       
C     NEWCOMB'S FORMULA                                                 
C     REF: ASTRONOMICAL EPHEMERIS 1976, P. 531.                         
  20  ST = DMOD(STD50R+(OMT50R+OMQ50R*TJUL)*TJUL,TWOPI)                 
C                                                                       
  30  SST=DSIN(ST)                                                      
      CST=DCOS(ST)                                                      
      P(1) = CNU*CST-SNU*SDE*SST                                        
      P(2) = -CNU*SST-SNU*SDE*CST                                       
      P(4) = SST*CDE                                                    
      P(5) = CST*CDE                                                    
      P(7) = -SNU*CST-CNU*SDE*SST                                       
      P(8) = SNU*SST-CNU*SDE*CST                                        
      GOTO 50                                                           
  40  P(1) = CNU                                                        
      P(2) = -SNU*SDE                                                   
      P(4) = 0.D0                                                       
      P(5) = CDE                                                        
      P(7) = -SNU                                                       
      P(8) = -CNU*SDE                                                   
  50  P(3) = SNU*CDE                                                    
      P(6) = SDE                                                        
      P(9) = CNU*CDE                                                    
      RETURN                                                            
      END                                                               
