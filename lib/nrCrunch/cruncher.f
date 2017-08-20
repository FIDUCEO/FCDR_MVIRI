C FILE: cruncher.f
      SUBROUTINE radiance(N,C,CF,CS,R)
C
C     CALCULATE Reflectance for N x N pixel
C     AUTHOR: Frank Ruethrich EUMETSAT
C     COMPILE WITH:
c     f2py -c -m cruncher cruncher.f
C
      INTEGER N
      REAL*4 CF
      REAL*4 CS
      INTEGER*2 C(N,N)
      REAL*8 R(N,N)
Cf2py intent(in) N
Cf2py intent(in) C
Cf2py intent(in) CF
Cf2py intent(in) CS
Cf2py intent(out) R
Cf2py integer intent(hide),depend(C) :: N=shape(C,0)

      DO 1 I=1,N
        DO 2 J=1,N
         R(I,J)=CF*(C(I,J)-CS)
2       ENDDO
1     ENDDO
      END

      SUBROUTINE reflectance(N,C,CF,CS,SZA,MSK,IRR,SUNDIST,R)
C
C     CALCULATE Reflectance for N x N pixel
C     AUTHOR: Frank Ruethrich EUMETSAT
C     COMPILE WITH:
c     f2py -c -m cruncher cruncher.f
C
      INTEGER    N,I,J
      INTEGER*2  C(N,N)
      INTEGER    MSK(N,N)

      REAL*4     SZA(N,N)
      REAL*4     CF,CS,IRR
      REAL*8     SUNDIST

      REAL*8     PI,TORAD
      REAL*8     Z

      REAL*8     R(N,N)
Cf2py intent(in) N

Cf2py intent(in) C
Cf2py intent(in) CF
Cf2py intent(in) CS
Cf2py intent(in) SZA
Cf2py intent(in) MSK
Cf2py intent(in) IRR
Cf2py intent(in) SUNDIST

Cf2py intent(out) R
Cf2py integer intent(hide),depend(C) :: N=shape(C,0)

      PI    = 3.14159265359D0
      TORAD = PI/180D0
      DO 3 I=1,N
        DO 4 J=1,N
         IF (MSK(I,J).EQ.1) THEN
           R(I,J)= CF*(C(I,J)-CS)
           Z     = SZA(I,J)*TORAD
           R(i,J)=  (PI*R(I,J)*(SUNDIST*SUNDIST)) / (IRR*DCOS(Z))
         ELSE
           R(i,J)=-9.
         ENDIF
4       ENDDO
3     ENDDO
      END

      SUBROUTINE EvalGeoDLat(flLatitude,dbEarth_Polar_radius_meter,
     &dbEarth_radius_meter,dbGeoLat)
      REAL*8 flLatitude
      REAL*8 dbEarth_Polar_radius_meter
      REAL*8 dbEarth_radius_meter
      REAL*8 dbGeoLat
      REAL*8 dbRadEq
      REAL*8 dbLatRad
      REAL*8 dbDiff
      REAL*8 dbTmpout
      REAL*8 PI,PI2
      REAL*8 DOUBLE_EPSILON

      PI               =3.14159265359D0
      PI2              =PI*0.5
      TORAD            =PI/180.D0
      TODEG            =180.D0/PI
      DOUBLE_EPSILON = 1.0E-10

      dbRadEq = (dbEarth_radius_meter -dbEarth_Polar_radius_meter )
     &/dbEarth_radius_meter;
      dbLatRad = flLatitude*TORAD
      dbDiff = DABS((PI*0.5) - DABS(dbLatRad))

      if (dbDiff < DOUBLE_EPSILON) THEN
          dbTmpout = dbLatRad;
      else
        dbRadEq = (1. - dbRadEq)*(1. - dbRadEq);
        dbTmpout = atan(tan(dbLatRad)*dbRadEq);
        dbTmpout = dbTmpout*TODEG;
      ENDIF
      
      dbGeoLat = dbTmpout;
      END

      SUBROUTINE PixRadius(flLatitude,dbEarth_Polar_radius_meter
     &,dbEarth_radius_meter,dbRadius)

      REAL*8 dbEarth_Polar_radius_meter,dbEarth_radius_meter
      REAL*8 flLatitude,dbRadius
      REAL*8 dbTmpout
      REAL*8 dbExc
      REAL*8 dbExc2
      REAL*8 dbExc4
      REAL*8 dbSp1
      REAL*8 dbSp12
      REAL*8 dbSp14
      REAL*8 dbLatRad
      REAL*8 TORAD,PI
      
      PI               =3.14159265359D0
      TORAD=PI/180.D0

      dbExc = dsqrt((dbEarth_radius_meter * dbEarth_radius_meter) - 
     &(dbEarth_Polar_radius_meter * dbEarth_Polar_radius_meter))/
     &(dbEarth_radius_meter)
      dbExc2 = dbExc*dbExc
      dbExc4 = dbExc2*dbExc2
      
      dbLatRad = flLatitude*TORAD
      dbSp1 = dsin(dbLatRad)
      dbSp12 = dbSp1*dbSp1
      dbSp14 = dbSp12*dbSp12
      
      dbTmpout = dbEarth_radius_meter * (1. -(dbExc2*dbSp12)/2. + 
     &(dbExc4*dbSp12)/2. - (5.*dbSp14*dbExc4)/8.)
      
      dbRadius = dbTmpout
      END

      SUBROUTINE vza(N,M,SSLAT,SSLON,MSK,LAT,LON,S,Z)
C
C     CALCULATE VZA(=S) and VAZ(=Z) for N x N pixel
C     AUTHOR: Frank Ruethrich EUMETSAT
C     COMPILE WITH:
c     f2py -c -m cruncher cruncher.f
C
      implicit none
      INTEGER    N,M
      REAL*8     SSLAT,SSLON
      REAL*8     LAT(N,N),LON(N,N)
      INTEGER    MSK(M,M)

      INTEGER    I,II,J,JJ
      REAL*8     dbEarth_Polar_radius_meter,dbEarth_radius_meter
      REAL*8     dbSat_Eight_meter
      REAL*8     dbPixLonRad,dbPixLatRad
      REAL*8     dbGeoSatLonRad,dbGeoDSatLat,dbGeoDLat
      REAL*8     TODEG,TORAD            
      REAL*8     PI,PI2
      REAL*8     dbTmpVal,dbTmpValAzi,dbSaph1,dbSaph2,dbC1Sat
      REAL*8     dbRadEq
      REAL*8     DOUBLE_EPSILON,dbGeoDLatRad
!       REAL*8     dbExc,dbExc2,dbExc4,dbSp1,dbSp12,dbSp14,
      REAL*8     dbRadius
      REAL*8     dbGeoDSatLatRad
      REAL*8     dbXSat,dbYSat,dbZSat,dbPixHeight
      REAL*8     dbXPix,dbYPix,dbZPix
      REAL*8     dbDro,dbXGeo,dbYGeo,dbZGeo
      REAL*8     dbTGeo,dbTSat,dbTSatGeo,dbGeoSatDAngle
      REAL*8     dbArg2,flZenith
      INTEGER    iLabSat
      REAL*8  iFac,dbV90X,dbV90Y,dbV90Z,dbV90,dbY0X,dbY0Y,dbY0Z,dbY0
      REAL*8  dbSSX,dbSSY,dbSSZ,dbSS
      REAL*8  dbProdX,dbProdY,dbProdZ,dbProd
      REAL*8  dbC2Sat
      REAL*8  dbTmpAr(3)

      REAL*8     S(N,N),Z(N,N)
Cf2py intent(in) N
Cf2py intent(in) M

Cf2py intent(in) SSPLAT
Cf2py intent(in) SSPLON
Cf2py intent(in) MSK
Cf2py intent(in) LAT
Cf2py intent(in) LON

Cf2py intent(out) S
Cf2py intent(out) Z

Cf2py INTEGER intent(hide),depend(MSK) :: M=shape(MSK,0)
Cf2py INTEGER intent(hide),depend(LAT) :: N=shape(LAT,0)
      PI               =3.14159265359D0
      PI2              =PI*0.5
      TORAD            =PI/180.D0
      TODEG            =180.D0/PI
      dbEarth_Polar_radius_meter =  6356.755 * 1000.
      dbEarth_radius_meter       =  6378.140 * 1000.
      dbSat_Eight_meter          =  42164.0  * 1000.

      dbGeoSatLonRad   = SSLON*TORAD

      dbTmpVal    = 0.0D0
      dbTmpValAzi = 0.0D0
      dbSaph1 = 0.0D0
      dbSaph2 = 0.0D0
      dbC1Sat = 0.0D0
      dbRadEq = (dbEarth_radius_meter -dbEarth_Polar_radius_meter )
     &/dbEarth_radius_meter
      DOUBLE_EPSILON = 1.0E-10

!     Transform satellite latitude from geographical to geocentric
      CALL EvalGeoDLat(SSLAT,dbEarth_Polar_radius_meter,
     &dbEarth_radius_meter,dbGeoDSatLat)
      dbGeoDSatLatRad=dbGeoDSatLat*TORAD
!     Terrestrial coordinates (meters) of satellite
      dbXSat = dbSat_Eight_meter * dcos(dbGeoDSatLatRad)*
     &dcos(dbGeoSatLonRad)
      dbYSat = dbSat_Eight_meter * dcos(dbGeoDSatLatRad)*
     &dsin(dbGeoSatLonRad)
      dbZSat = dbSat_Eight_meter * dsin(dbGeoDSatLatRad)

      DO 7 I=1,N
        DO 8 J=1,N

         IF (N>M) THEN
           II=INT(I/2+1)
           JJ=INT(J/2+1)
         ELSE
           II=I
           JJ=J
         END IF
         IF (MSK(II,JJ).EQ.1) THEN
           dbPixLatRad = LAT(I,J)*TORAD
           dbPixLonRad = LON(I,J)*TORAD
!          Calculate geocentric latitude of observation point
           CALL EvalGeoDLat(LAT(I,J),
     &dbEarth_Polar_radius_meter,
     &dbEarth_radius_meter,
     &dbGeoDLat)
           dbGeoDLatRad = dbGeoDLat*TORAD
!          Earth radius (in meters) at geographical latitude of observation point
           CALL PixRadius(LAT(I,J),dbEarth_Polar_radius_meter
     &,dbEarth_radius_meter,dbRadius)
           dbPixHeight = dbRadius
!          Terrestrial coordinates (meters) of observation point
           dbXPix = dbPixHeight * dcos(dbGeoDLatRad)*dcos(dbPixLonRad)
           dbYPix = dbPixHeight * dcos(dbGeoDLatRad)*dsin(dbPixLonRad)
           dbZPix = dbPixHeight * dsin(dbGeoDLatRad)
!          Components of a vector orthogonal to the
!          surface of the earth geoid
           if (DABS(LON(I,J) - SSLON) == 90.) THEN 
             dbDro = 0.
           END IF
           if (DABS(LON(I,J) - SSLON) == 0.)  THEN
              dbDro = sqrt(dbXPix*dbXPix + dbYPix*dbYPix)
           END IF
           if (DABS(LON(I,J) - SSLON) .GT. 0. .AND.
     &DABS(LON(I,J) - SSLON) .LT. 90.) THEN
              dbDro = dbZPix/dtan(dbPixLatRad)
           END IF
           dbXGeo = dbDro * dcos(dbPixLonRad)
           dbYGeo = dbDro * dsin(dbPixLonRad)
           dbZGeo = dbZPix
!          Vector norms
           dbTGeo = dsqrt((dbXGeo*dbXGeo) + (dbYGeo*dbYGeo) + 
     &(dbZGeo*dbZGeo))
           dbTSat = dsqrt((dbXSat*dbXSat) + (dbYSat*dbYSat) + 
     &(dbZSat*dbZSat))
           dbTSatGeo = dsqrt((dbXSat-dbXPix)*(dbXSat-dbXPix) + 
     &(dbYSat-dbYPix)*(dbYSat-dbYPix) + 
     &(dbZSat-dbZPix)*(dbZSat-dbZPix))
!          Geocentric angle between directions of satellite 
!          and observation point
           dbGeoSatDAngle = (dbXSat*dbXPix + dbYSat*dbYPix +
     &dbZSat*dbZPix)/(dbTSat*dbPixHeight);
           if (DABS(1.- DABS(dbGeoSatDAngle)) .LT. DOUBLE_EPSILON) THEN
             if (dbGeoSatDAngle > 0.) THEN 
               dbGeoSatDAngle = 1.
             END IF
             if (dbGeoSatDAngle < 0.) THEN 
               dbGeoSatDAngle = -1.
             END IF
           END IF
           dbGeoSatDAngle = dacos(dbGeoSatDAngle)*TODEG
!          Check if geocentric directions of satellite and 
!          observer are coincident
           iLabSat = 0
           if (dbGeoSatDAngle > 179.9 .OR.dbGeoSatDAngle <= 0.) THEN
             iLabSat = 1
           end if
!          Satellite zenith angle
           dbArg2 = (dbXPix*(dbXSat-dbXPix) + dbYPix*(dbYSat-dbYPix)+
     &dbZPix*(dbZSat-dbZPix))/(dbPixHeight*dbTSatGeo)
           if (DABS(1. - DABS(dbArg2)) < DOUBLE_EPSILON) THEN
              if (dbArg2 > 0.) THEN 
                dbArg2 = 1.
              end if
              if (dbArg2 < 0.) THEN 
                dbArg2 = -1.
              end if
           end if
           flZenith=dacos(dbArg2)*TODEG
           S(I,J) = flZenith
           Z(I,J) = 0.

!          Satellite Azimuth angle
   
!          Define vector v0 in the meridional plane, perpendicular 
!          to vertical of obervation point and pointing to north pole 
!          [see that vertical direction is not necessarily geocentric direction]


           IF (DABS(90. - DABS(dbGeoDLat)) < DOUBLE_EPSILON) THEN
            iFac = -1
            if (dbGeoDLat .LT. 0.) THEN
              iFac = -iFac
            end if
            dbV90X = iFac*dbPixHeight*cos(dbPixLonRad)
            dbV90Y = iFac*dbPixHeight*sin(dbPixLonRad)
            dbV90Z = 0.
           else
            dbV90X = -dbXGeo*dbZGeo
            dbV90Y = -dbYGeo*dbZGeo
            dbV90Z = (dbXGeo*dbXGeo) + (dbYGeo*dbYGeo)
           end if
           dbV90 = dsqrt((dbV90X*dbV90X)+(dbV90Y*dbV90Y)+
     &(dbV90Z*dbV90Z))
   
!           Define vector dbY0 orthogonal to the vectors dbV90 and flGeo
!           and forming with them a 3-axis right system (dbY0 points 
!           eastwards)
            dbY0X = dbV90Y*dbZGeo - dbYGeo*dbV90Z
            dbY0Y = dbV90Z*dbXGeo - dbV90X*dbZGeo
            dbY0Z = dbV90X*dbYGeo - dbXGeo*dbV90Y
            dbY0 = dsqrt((dbY0X*dbY0X)+(dbY0Y*dbY0Y)+(dbY0Z*dbY0Z))
!           Define the vector dbSS pointing to the satellite at observation 
!           point
            dbSSX = dbXSat - dbXPix
            dbSSY = dbYSat - dbYPix
            dbSSZ = dbZSat - dbZPix
            dbSS = dsqrt((dbSSX*dbSSX)+(dbSSY*dbSSY)+(dbSSZ*dbSSZ))

!           Calculate vector on the plane tangent to the observation point
!           and pointing to the satellite = flGeo x dbSS x flGeo
            if (iLabSat == 0) THEN
                dbTmpAr(0+1) = dbYGeo*dbSSZ - dbZGeo*dbSSY
                dbTmpAr(1+1) = dbZGeo*dbSSX - dbXGeo*dbSSZ
                dbTmpAr(2+1) = dbXGeo*dbSSY - dbYGeo*dbSSX
                
                dbProdX = dbTmpAr(1+1)*dbZGeo - dbTmpAr(2+1)*dbYGeo
                dbProdY = dbTmpAr(2+1)*dbXGeo - dbTmpAr(0+1)*dbZGeo
                dbProdZ = dbTmpAr(0+1)*dbYGeo - dbTmpAr(1+1)*dbXGeo
                dbProd = dsqrt((dbProdX*dbProdX) + (dbProdY*dbProdY)
     & + (dbProdZ*dbProdZ))

                if (dbProd > 0.) THEN
!               Calculate angles wrt vectors dbY0 and dbV90
!               If angle wrt to y0 exceeds 90 deg => satellite 
!               azimuth > 180 deg 
                dbC1Sat = (dbV90X*dbProdX + dbV90Y*dbProdY + 
     &dbV90Z*dbProdZ)/(dbV90*dbProd)
                dbC2Sat = (dbY0X*dbProdX + dbY0Y*dbProdY + 
     &dbY0Z*dbProdZ)/(dbY0*dbProd)
                if (DABS(1. - DABS(dbC1Sat)) < DOUBLE_EPSILON) THEN
                  if (dbC1Sat > 0.) THEN 
                    dbC1Sat = 1. 
                  ENDIF
                  if (dbC1Sat < 0.) THEN 
                    dbC1Sat = -1. 
                  ENDIF
                ENDIF
                if (DABS(1. - DABS(dbC2Sat)) < DOUBLE_EPSILON) THEN
                  if (dbC2Sat > 0.) THEN 
                    dbC2Sat = 1.
                  ENDIF
                  if (dbC2Sat < 0.) THEN 
                    dbC2Sat = -1. 
                  ENDIF
                ENDIF
              ELSE  ! dbProd 
!             Special case of exterior product = 0 <> dbSS and flGeo have same
!             direction It occurs if the satellite is at the zenith (or nadir).
!             Its azimuth is undetermined. For geosynchronous satellites we 
!             adopt the same solution as for Sun
!             Minor modifications will be needed for non-geosynchronous platforms 
                if (dbZPix >= dbZSat) THEN
                    if (S(I,J) <= 90.) THEN
                      dbC1Sat = -1.
                      dbC2Sat = 0.
                    else
                      dbC1Sat = +1.
                      dbC2Sat = 0.
                    END IF
                END IF
                if (dbZPix < dbZSat) THEN
                    if (S(I,J) <= 90.) THEN
                      dbC1Sat = +1.
                      dbC2Sat = 0.
                    else
                      dbC1Sat = -1.
                      dbC2Sat = 0.
                    END IF
                END IF
              END IF ! dbProd
              dbSaph1 = (dacos(dbC1Sat))*TODEG
              dbSaph2 = (dacos(dbC2Sat))*TODEG
              if (dbSaph2 <= 90.) THEN
                dbTmpValAzi = dbSaph1
              end if
              if (dbSaph2 > 90.)  THEN
                dbTmpValAzi = 360. - dbSaph1
              end if
              Z(I,J) = dbTmpValAzi
            else ! iLabSat
!        Observation point and satellite directions are coincident
!        These are special cases of the satellite azimuth
!        The solutions used here are applicable to geosynchronous 
!        satellites only
!        Minor modifications will be needed for non-geosynchronous platforms
              if (dbZPix >= dbZSat .AND. dbGeoSatDAngle < 90.) THEN
                dbTmpValAzi = 180.
              end if
              if (dbZPix >= dbZSat .AND. dbGeoSatDAngle > 90.) THEN
                dbTmpValAzi = 0.
              end if
              if (dbZPix <  dbZSat .AND. dbGeoSatDAngle < 90.) THEN
                dbTmpValAzi = 0.
              end if
              if (dbZPix <  dbZSat .AND. dbGeoSatDAngle > 90.) THEN
                dbTmpValAzi = 180.
              end if
              Z(I,J) =dbTmpValAzi
            end if

         ELSE
           S(I,J) = -999.
           Z(I,J) = -999.
         END IF


8       ENDDO
7     ENDDO

      END


      SUBROUTINE sza(N,M,T,DOY,MSK,LAT,LON,S,Z)
C
C     CALCULATE SZA for N x N pixel
C     Frank Ruethrich EUMETSAT
C     COMPILE WITH:
c     f2py -c -m cruncher cruncher.f
c     DEBUG with:
c     rm *.so
c     f2py -m cruncher -h --overwrite-signature cruncher.pyf cruncher.f
c     f2py -c --debug --build-dir build cruncher.pyf cruncher.f
c     gdb python #then press run timeseries_cruncher.py
c     OR:
c     f2py --debug-capi -c -m cruncher cruncher.f

C
      INTEGER    N,M,DOY
      REAL*8     T(M,M),LAT(N,N),LON(N,N)
      INTEGER    MSK(M,M)

      INTEGER    I,II,J,JJ
      REAL*8     PI               
      REAL*8     TORAD            
      REAL*8     TODEG            
      REAL*8     flNumDayInTheYear
      REAL*8     MSA_EQTIME1      
      REAL*8     MSA_EQTIME2      
      REAL*8     MSA_EQTIME3      
      REAL*8     MSA_EQTIME4      
      REAL*8     MSA_EQTIME5      
      REAL*8     MSA_EQTIME6      
      REAL*8     MSA_DECL1        
      REAL*8     MSA_DECL2        
      REAL*8     MSA_DECL3        
      REAL*8     MSA_DECL4        
      REAL*8     MSA_DECL5        
      REAL*8     MSA_DECL6        
      REAL*8     MSA_DECL7        
      REAL*8     flTimeZone       

      REAL*8     flRadLat 
      REAL*8     flRadLon 
      REAL*8     flGamma,flEquTime,flDecli,flTimeOffset
      REAL*8     flTrueSolarTime,flHa,flHaRad,flCosZen
      REAL*8     flTmpZenRad,flTmpZen,flCosAzi,flTmpAzi

      REAL*8     S(N,N)
      REAL*8     Z(N,N)
! Cf2py intent(in) N
! Cf2py intent(in) M

Cf2py intent(in) T
Cf2py intent(in) DOY
Cf2py intent(in) MSK
Cf2py intent(in) LAT
Cf2py intent(in) LON

Cf2py intent(out) S
Cf2py intent(out) Z
Cf2py integer intent(hide),depend(T) :: M=shape(T,0)
Cf2py integer intent(hide),depend(LAT) :: N=shape(LAT,0)


      PI               =3.14159265359D0
      TORAD            =PI/180.D0
      TODEG            =180.D0/PI
      flNumDayInTheYear=365.D0
      MSA_EQTIME1      =229.18D0
      MSA_EQTIME2      =0.000075
      MSA_EQTIME3      =0.001868
      MSA_EQTIME4      =0.032077
      MSA_EQTIME5      =0.014615
      MSA_EQTIME6      =0.040849
      MSA_DECL1        =0.006918
      MSA_DECL2        =0.399912
      MSA_DECL3        =0.070257
      MSA_DECL4        =0.006758
      MSA_DECL5        =0.000907
      MSA_DECL6        =0.002697
      MSA_DECL7        =0.00148
      flTimeZone       =0.


      DO 5 I=1,N
        DO 6 J=1,N
         IF (N>M) THEN
           II=INT(I/2+1)
           JJ=INT(J/2+1)
         ELSE
           II=I
           JJ=J
         END IF
         IF (MSK(II,JJ).EQ.1) THEN
           flRadLat = LAT(I,J)*TORAD
           flRadLon = LON(I,J)*TORAD
!          Evaluate the fractional year in radians
           flGamma =  2*PI*((DOY)+((T(II,JJ))/24.))/flNumDayInTheYear !checked
!          Evaluate the Equation of time in minutes
           flEquTime = MSA_EQTIME1*(MSA_EQTIME2+MSA_EQTIME3*          !checked
     &dcos(flGamma)-MSA_EQTIME4*dsin(flGamma)-
     &MSA_EQTIME5*dcos(2.*flGamma)-MSA_EQTIME6*
     &dsin(2.*flGamma))      
!          Evaluate the solar declination angle in radians */
           flDecli = MSA_DECL1-MSA_DECL2*dcos(flGamma)+MSA_DECL3*     !checked
     &dsin(flGamma)- MSA_DECL4*dcos(2*flGamma)+MSA_DECL5
     &*dsin(2*flGamma)-MSA_DECL6*dcos(3*flGamma)+
     &MSA_DECL7*dsin(3*flGamma)
!          Time offset in minutes equivalent to 
!          here was an error
           flTimeOffset = flEquTime-(4.D0*(-1.D0*LON(I,J)))+
     &(60.D0*flTimeZone)!fixed and checked
!          True solar time in minutes */
           flTrueSolarTime = (T(II,JJ)*60.)+flTimeOffset
!          solar hour angle in degrees and in radians 
           flHa = (flTrueSolarTime/4.)-180.
           flHaRad = TORAD*flHa
!          Evaluate the Solar local Coordinates
           flCosZen = (dsin(flRadLat)*dsin(flDecli)+ dcos(flRadLat)*
     &dcos(flDecli)*dcos(flHaRad))
           flTmpZenRad = dacos(flCosZen)  
           flTmpZen = TODEG*(flTmpZenRad)      
           flCosAzi = -((dsin(flRadLat)*dcos(flTmpZenRad)-
     &dsin(flDecli))/(dcos(flRadLat)*dsin(flTmpZenRad)))
           flTmpAzi =  360. - TODEG*(dacos(flCosAzi))
!        Correct for SZA < 180
           IF (flTrueSolarTime < 720.) THEN   
              flTmpAzi = 360. - flTmpAzi       
           END IF
           if (LAT(I,J).eq.0.) THEN
!             print*, LON(I,J)," ",flTmpZen," ",flHa," ",flTimeOffset,
!      &" ", flEquTime ," ",flDecli*TODEG
           end if
           S(I,J) = flTmpZen    !checked
           Z(I,J) = flTmpAzi    !checked
!            if (LAT(I,J) > 28.53 .and. LAT(I,J) < 28.57) THEN !for lybia 4
!              if (LON(I,J) > 23.37 .and. LON(I,J) < 23.41) THEN
!                print*, LON(I,J)," ",LAT(I,J)," ",flTmpZen,
!      &" ",flTmpAzi
!              end if
!            end if
         ELSE
           S(I,J) = -999.
           Z(I,J) = -999.
         END IF
6       ENDDO
5     ENDDO
      END


      subroutine refgeo (N,line, pixel, lat, lon, visible)
C This subroutine converts digital to geographical co-ordinates.

C Input parameters:
C line - line number, measured from southern end of frame
C pixel - pixel number, measured from eastern end of frame
C Line and pixel values are real numbers to enable sub-pixel
C accuracy. Integer values correspond to the middle of the pixel,
C e.g.(500, 800) would correspond to the middle of the pixel with
C corners (499.5, 799.5), (499.5, 800.5), (500.5, 799.5), (500.5, 800.5).

C Output parameters
C lat - latitude of this pixel (degrees North from Equator)
C lon - longitude of this pixel (degrees East from Greenwich)
C visible - flag set to TRUE if pixel is on visible disc,
C         - flag set to FALSE if pixel is in space.

      implicit none
      real*4 lat, lon, line, pixel, cenlat,N
      real*4 aline, asamp, tanal, tanas, det, k
      real*4 altitude, req, rpol, oblate, pi, deg_to_rad, rad_to_deg
      real*4 step, x, y, z, a, b, c, p, q, r
      logical*1 visible

C Set up constants.
C altitude = distance from earth centre to satellite
C req = Equatorial earth radius
C rpol = Polar earth radius
C oblate = earth oblateness
C deg_to_rad and rad_to_deg are conversion factors


      altitude = 42164.0
      req = 6378.140
      rpol = 6356.755
      oblate = 1.0 / 298.257
      pi = 3.141592653
      deg_to_rad = pi / 180.0
      rad_to_deg = 180.0 / pi

C Step is the radiometer step as seen by the spacecraft,
C in degrees. The image represents an 18° x 18° field
C of view divided up on an equi-angular basis. For this
C program an IR channel of 2500 x 2500 is assumed but
C in the real code the size of each channel must be accounted C for.

      step = 18.0 / N
C Convert line/pixel values to angular offsets from centre point
c FRue: added +0.5 below to correct grid
      asamp = - (pixel - (N/2+0.5)) * step
      aline = (line - (N/2+0.5)) * step
      asamp = asamp * deg_to_rad
      aline = aline * deg_to_rad
C Calculate tangents of angles
      tanal = tan(aline)
      tanas = tan(asamp)
C Calculate components of an arbitrary vector from the spacecraft
C in the viewing direction.
      p = -1.
      q = tanas
      r = tanal * sqrt (1.+ q*q)
C The location of the point on the earth can be identified by
C solving a quadratic equation for the intersection between
C the earth's surface and the viewing line from the spacecraft.
C If this equation has no real roots then there is no intersection;
C otherwise the required root is the one nearer to the spacecraft
C (on the visible side of the earth).
      a = q*q + (r*req/rpol)**2 + p*p
      b = 2. * altitude * p
      c = altitude*altitude - req*req
C Calculate determinant. If it is negative (no real roots to
C quadratic equation) there is no intersection between the
C line of sight and the disc and so the pixel does not correspond
C to visible data.
      det = b*b - 4 * a * c
      if (det .le. 0.) then
        visible = .false.
        lat = 0.
        lon = 0.
        goto 999
      else
        visible = .true.
      endif
      k = (- b - sqrt(det)) / (2. * a)
      x = altitude + k * p
      y = k * q
      z = k * r
      lon = atan (y/x)
      cenlat = atan (z * cos(lon) / x)
C This is the geocentric latitude. Convert it to the geodetic
C (or geographic) latitude before returning it to the calling program
      lat = atan(tan(cenlat)/((1.0-oblate)**2))
C Convert from radians to degrees
      lat = lat * rad_to_deg
      lon = lon * rad_to_deg
!       if (lon.lt.0.) then
!         lon=360-lon*-1
!       end if
!       print*, lat
999   return
      end


      SUBROUTINE latlon(N,LAT,LON)
C
C     CALCULATE Lat/Lon for N x N pixel
C     AUTHOR: Frank Ruethrich EUMETSAT
C     COMPILE WITH:
c     f2py -c -m cruncher cruncher.f
C
      INTEGER N
      logical*1 visible

      INTEGER I,J

      REAL*4 P,L

      REAL*4 la,lo,LAT(N,N),LON(N,N)

      

Cf2py intent(in) N
Cf2py intent(out) LAT
Cf2py intent(out) LON
      call refgeo(5000.,2500,2500,la,lo,visible)
!       print*,la
!       print*,lo
      call refgeo(5000.,2500.5,2500.5,la,lo,visible)
!       print*,la
!       print*,lo
      DO 1 I=1,N
        DO 2 J=1,N
          P=I
          L=J
         call refgeo(float(N),P,L,la,lo,visible)
!          print*, la
!          print*,""
         IF (visible .EQV. .TRUE.) THEN
          LAT(I,J)=la
          LON(I,J)=lo
         ELSE
          LAT(I,J)=-999.
          LON(I,J)=-999.
         END IF
2       ENDDO
1     ENDDO
      END




      subroutine georef (rlat, rlong, N, NOMSSP, line, pixel, visible)

C      This subroutine converts pixel position from geographical
C      (lat / long) co-ordinates to digital (line / pixel) co-ordinates.

C  Input parameters:
C       rlat    - latitude of pixel (North is +ve, South is ve)
C       rlong   - longitude of pixel (East is +ve, West is -ve)
C      Note that these are standard geographic co-ordinates as would
C      be found in an atlas.

C  Output parameters

C       line    - line number, measured from southern end of frame
C       pixel   - pixel number, measured from eastern end of frame
C       visible - flag set to TRUE if pixel is on visible disc,
C               - flag set to FALSE if pixel is in space.
C       (c) EUMETSAT 1997

      implicit none
      real*4 rlat, rlong, lat, long, geolat, NOMSSP
      real*4 altitude, req, rpol, oblate, pi, deg_to_rad, rad_to_deg
      real*4 x, y, z, rtheta, aline, asamp, dotprod
      integer*4 line, pixel, nlines, nsamps,N
      logical*1 visible
C       Set up constants.
C       altitude = distance from earth centre to satellite
C       req = Equatorial earth radius
C       rpol = Polar earth radius
C       oblate = earth oblateness
C       deg_to_rad and rad_to_deg are conversion factors
      altitude = 42164.0
      req = 6378.140
      rpol = 6356.755
      oblate = 1.0 / 298.257
      pi = 3.141592653
      deg_to_rad = pi / 180.0
      rad_to_deg = 180.0 / pi
C       Convert inputs to radians

      geolat = rlat * deg_to_rad 
      long = (rlong - NOMSSP) * deg_to_rad
C      Convert geodetic latitudes (as input) to geocentric latitudes
C      for use within the algorithm

      lat = atan ( ((1.0 - oblate) **2) * tan (geolat) )

C      Calculate rtheta. This is the distance from the earth centre to
C      a point on the surface at latitude 'lat'.

      rtheta = (req*rpol)/sqrt(rpol**2*cos(lat)**2 + req**2*sin(lat)**2)

C      Calculate Cartesian co-ordinates of target point. This is
C      basic geometry. The co-ordinate system is geocentric with
C      the x-axis towards the spacecraft, the y-axis to the East
C       and the x-axis towards the N pole. 

      x       = rtheta * cos(lat) * cos(long) 
      y       = rtheta * cos(lat) * sin(long) 
      z       = rtheta * sin(lat) 

C      Check for invisibility. This is done using the basic geometric
C      theorem that the dot product of two vectors A and B is equal
C      to
C              |A||B| cos (theta)
C
C      where theta is the angle between them. In this case, the test
C      is simple. The horizon is defined as the locus of points where
C      the local normal is perpendicular to the spacecraft sightline
C      vector. All visible points have (theta) less than 90°
C      and all invisible points have (theta) greater than 90°.
C      The test therefore reduces to whether the sign of the dot
C      product is +ve or -ve; if it is -ve the point is invisible.
C      The vector from the point to the spacecraft has components
C      Rs-x, -y, -z where Rs is the distance from the origin to the
C      satellite. The vector for the normal has components
C      x  y  z(Re/Rp)^2

      dotprod = (altitude-x)*x - y*y - z*z*((req/rpol)**2) 
      if (dotprod .le. 0.) then
          visible = .false. 
          line = 0
          pixel = 0 
          goto 999
      else
          visible = .true. 
      endif
C      In this co-ordinate system the spacecraft (S) is at position
C      (altitude,0,0), the earth centre (O) at (0,0,0) and the point (P)
C      at (x,y,z). Two additional points need to be defined, so that the
C      angles from the reference planes to the target point (i.e. the
C      position of the point in the sensor FOV) can be extracted.
C      These points are defined by dropping lines perpendicularly from P
C      onto the equatorial plane and the Greenwich meridian plane.
C      Their co-ordinates are defined as:
C
C      O' = (x, y, 0)      and      O'' = (x, 0, z).
C
C With these points, right-angled triangles can be defined SO'P 
C and SO''P which can be used directly to determine the angular 
C co-ordinates (aline, asamp) of P in the FOV.

      asamp = atan (y / (altitude - x))
      aline = atan (z / sqrt (y**2 + (altitude - x)**2))

C      Convert back to degrees

      asamp = asamp * rad_to_deg 
      aline = aline * rad_to_deg

C  Calculate line, pixel. Note that since pixels are measured from
C the right of the image, and the angular conversion was measured in 
C the x (east) direction, a sign correction has to be included for
C pixels. The image represents an 18° x 18° field of view 
C divided up on an equi-angular basis.

      nlines = N 
      nsamps = N
      asamp = asamp / (18.0 / float(nsamps)) 
      aline = aline / (18.0 / float(nlines))
C     FRue: added +0.5 below to correct grid to center
      if (asamp .ge. 0.0) then
        pixel = nsamps / 2 + 0.5  - int (asamp) 
      else
        pixel = nsamps / 2 + 0.5 + 1 - int (asamp) 
      endif

      if (aline .ge. 0.0) then
        line = nlines / 2 + 1 + int (aline) 
      else
        line = nlines / 2 + int (aline) 
      endif

!            print*, line
!            print*, pixel

999   continue 
      return 
      end

      SUBROUTINE linepixel(N,LAT,LON,NOMSSP,LINE,PIXEL)
C
C     CALCULATE pixel/line of a pair of lat/lon coorinates
C     Basically a wrapper for georef()
C     AUTHOR: Frank Ruethrich EUMETSAT
C     COMPILE WITH:
c     f2py -c -m cruncher cruncher.f
C
      logical*1 visible

      REAL*4 LAT,LON,NOMSSP
      INTEGER*4 PIXEL,LINE
      

Cf2py intent(in) LAT
Cf2py intent(in) LON
Cf2py intent(in) NOMSSP
Cf2py intent(out) LINE
Cf2py intent(out) PIXEL

      call georef(LAT,LON,N,NOMSSP,LINE,PIXEL,visible)
!           print*, LINE
!           print*, PIXEL
      IF (visible .NEQV. .TRUE.) THEN
        LINE=-999.
        PIXEL=-999.
      ENDIF

      END


C END FILE cruncher.f