C                                                                       
C...................................................................... 
C                                                                       
C SUBROUTINE  32023021  COOT     VERSION 02    890712                   
C OAD              FORTRAN                                              
C                                                                       
C SET VARIABLES FOR OAD ROUTINES                                        
C                                                                       
C CALLING SEQUENCE -                                                    
C CALL COOT                                                             
C                                                                       
C EXPLANATION -                                                         
C  PREPARATION OF SOME IMPORTANT CONSTANTS FOR EARTH SATELLITE DYNAMICS 
C  IF NOT ESPECIALLY MENTIONED, ALL VALUES ARE EXPRESSED IN KM, SEC, RAD
C  UPDATED IN JULY 1984.                                                
C                                                                       
C SPECIAL FEATURES -                                                    
C DOCUMENTATION ADDED BY R.A.BOWEN TO UNCHANGED OAD SOURCE              
C                                                                       
C CHANGES -                                                             
C     01 NEW                                                            
C     02 890712  NEW IDENTIFICATION NUMBER 023                          
C                                                                       
C...................................................................... 
C                                                                       
      SUBROUTINE COOT                                                   
C  PREPARATION OF SOME IMPORTANT CONSTANTS FOR EARTH SATELLITE DYNAMICS 
C  IF NOT ESPECIALLY MENTIONED, ALL VALUES ARE EXPRESSED IN KM, SEC, RAD
C  UPDATED IN JULY 1984.                                                
C
      IMPLICIT NONE      ! MTP
      REAL*8 STD50,OMT50,OMQ50,OMROT,STD50R,OMT50R,OMQ50R ! MTP      
      REAL*8 PI,TWOPI,PIH,DEG,RAD                        ! MTP
      REAL*8 RE,FLAT,FLATSQ,CLIGHT                            ! MTP  
      REAL*8 XMU,WMU,WMUIN,XMUM,XMUS                            ! MTP
      REAL*8 SOLPR,ECL,CECL,SECL,XLSUNA,XLSUNB,ECLECC,ECLOM     ! MTP
      REAL*8 PRACCY,PRACCZ                                      ! MTP
      REAL*8 TODEG,AN,BN2,CN     ! MTP
C
      COMMON /CDYNAE/ STD50,OMT50,OMQ50,OMROT,STD50R,OMT50R,OMQ50R      
      COMMON /CBASIC/ PI,TWOPI,PIH,DEG,RAD                              
      COMMON /CGEOME/ RE,FLAT,FLATSQ,CLIGHT                             
      COMMON /CPOEMS/ XMU,WMU,WMUIN,XMUM,XMUS                           
      COMMON /CSOLAR/ SOLPR,ECL,CECL,SECL,XLSUNA,XLSUNB,ECLECC,ECLOM    
      COMMON /CORIOL/ PRACCY,PRACCZ                                     
C                                                                       
      PIH = DATAN2(1.D0,0.D0)                                           
      PI = PIH + PIH                                                    
      TWOPI = PI+PI                                                     
      DEG = 180.D0/PI                                                   
      RAD = PI/180.D0                                                   
C                                                                       
C  RE = MEAN RADIUS OF EQUATOR;   FLAT = FLATTENING COEFFICIENT         
C  REF: GODDARD EARTH MODEL 6                                           
      RE = 6378.144 D0                                                  
      FLAT = 1.D0/298.257 D0                                            
      FLATSQ = FLAT*(2.D0-FLAT)                                         
C  CLIGHT = SPEED OF LIGHT. INTERNATIONALLY ADOPTED VALUE.              
      CLIGHT = 299792.458 D0                                            
C                                                                       
C  STD50, OMT50, OMQ50 GIVE THE ANGLE BETWEEN GREENWICH MERIDIAN        
C  AND MEAN EQUINOX (DEG) FROM  D = MODIFIED JULIAN DATE (1950)         
C  ACCORDING TO  W = STD50 + OMT50*D + OMQ50*D*D                        
C  STD50R, OMT50R, OMQ50R  ARE THE SAME ANGLES EXPRESSED IN RADIANS     
C  REF: THE ASTRONOMICAL EPHEMERIS 1976, P. 531. NEWCOMBS FORMULA:      
C  W = AN + BN*(O.5 + D/36525) + CN*(0.5 + D/35625)**2                  
      TODEG = 36.D1/864.D2                                              
C  AN = 6H + 38MIN + 45.836 SEC                                         
      AN = TODEG*((36.D1 + 38.D0)*6.D1 + 45.836D0)                      
C  BN = (100 TURNS + 184.542 SEC)/JULIAN CENTURY                        
      BN2 = TODEG*184.542 D0                                            
      CN = TODEG*0.0929D0                                               
      STD50 = AN + (BN2 + CN*0.5D0)*0.5D0                               
      OMT50 = 36.D1 + (36.D3 + BN2 + CN)/36525.D0                       
      OMQ50 = CN/36525.D0**2                                            
      STD50R = STD50*RAD                                                
      OMT50R = OMT50*RAD                                                
      OMQ50R = OMQ50*RAD                                                
C  OMROT = MEAN VELOCITY OF ROTATION OF THE EARTH (RAD/SEC)             
      OMROT = OMT50R/864. D2                                            
C                                                                       
C  XMU = CENTRAL EARTH POTENTIAL. REF: SMITHSONIAN STANDARD EARTH III   
C        AND GODDARD EARTH MODELS 5 AND 6                               
      XMU = 398601.3D0                                                  
      WMU = DSQRT(XMU)                                                  
      WMUIN=1.D0/WMU                                                    
C  XMUM, XMUS = MOON, SUN POTENTIAL. REF: JPL DEV. EPHEMERIDES NO. 19   
      XMUM = XMU/81.301 D0                                              
      XMUS = 132715.0D6                                                 
C                                                                       
C  SOLPR=  RADIATION PRESSURE AT MEAN EARTH DISTANCE FROM THE SUN       
C  IN KG*KM/(M*SEC)**2  REF: NASA TM-X 64627                            
      SOLPR = 4.51 D-9                                                  
C  ECL, CECL, SECL= INCLINATION OF THE ECLIPTIC AND ITS COS AND SIN .   
C  XLSUNA, XLSUNB, ECLECC, ECLOM GIVE WITH A PRECISION OF 0.01 DEG      
C  A VALUE FOR THE LONGITUDE OF THE SUN (LS(RAD)) IN THE ECLIPTIC       
C  FROM: LM=XLSUNA+XLSUNB*DAY; LS=LM+ECLECC*DSIN(LM-ECLOM)              
C  REF: THE ASTRONOMICAL EPHEMERIS 1978, CALCULATED FOR JAN 1985.       
      ECL = RAD*23.441242 D0                                            
      CECL = DCOS(ECL)                                                  
      SECL = DSIN(ECL)                                                  
      XLSUNA = RAD*280.081182 D0                                        
      XLSUNB = RAD*0.9856473495 D0                                      
      ECLECC = 0.33430838D-01                                           
      ECLOM = RAD*282.682456 D0                                         
C                                                                       
C  ROTATION *2 FROM PRECESSION OF MEAN SYSTEM OF DATE: CORIOLIS         
C  ACCELERATION = (-PRACCY*VY-PRACCZ*VZ,PRACCY*VX,PRACCZ*VX)            
      PRACCY = 0.14167 D-10                                             
      PRACCZ = 0.6158 D-11                                              
      RETURN                                                            
      END                                                               
