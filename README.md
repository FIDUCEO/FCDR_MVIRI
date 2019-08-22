# Meteosat Fundamental Climate Data Record  
  
Copyright (C) EUMETSAT, 2019  
This code was developed for the EC project “Fidelity and Uncertainty in   
Climate Data Records from Earth Observations (FIDUCEO)”. 
Grant Agreement: 638822  
  
Version 0.1.1 Reviewed and approved by Frank Rüthrich, EUMETSAT 19-08-2019  
  
 * This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option)
any later version.  
 * This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
more details.  
  
A copy of the GNU General Public License should have been supplied along
with this program; if not, see http://www.gnu.org/licenses/  


# FCDR generation Software for fiduceo  

## REQUIREMENTS  
#Strand 7!  
module purge  
source ~/.climproc/modules/00-generic_setup_with_menu.module  
### ENVIRONMENT VARIABLES  
export PREF="..."  
export RICalPyROOT=${PREF}/git/RICalPy  
#not on DSNNAS to be easily editable from X2go:  
export RICalPyCODE="..."  
#to install from git:  
export gitdir=${PREF}"/git/"  
export P27git=${PREF}"/git/python2.7"  
#to install from tarballs:  
tarballs=${PREF}"/tarballs/"  
export CONFIGPP=${P27git}/mviri_lib/mviri/config  

export PGDATA=${PREF}/data  

export PATH=$PATH:${PREF}/bin  
export PATH=$PATH:${PREF}/lib  

export PYTHONPATH=${PREF}/lib64/python2.7/site-packages/  
export PYTHONPATH=${PYTHONPATH}:${PREF}/lib/python2.7/site-packages/  

export PYTHON_DIR=${PREF}/pythons/RICalPy  

export LD_LIBRARY_PATH=${PREF}/lib/  
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$PYTHON_DIR/include/python2.7/
### STRUCTURES  
mkdir -p ${PREF}/lib/python2.7/site-packages/  
mkdir -p ${PREF}/lib64/python2.7/site-packages/  
mkdir -p $P27git  
mkdir -p $tarballs  
mkdir -p $P27git  
mkdir -p $PREF/DB  

ln -s $RICalPyCODE $RICalPyROOT #link from my home; may not be necessary once development is ready  

### get a decent Python 2.7 virtual evironment
recommended: anaconda2/4.3.1
### easy accessible libs
easy_install psycopg2  
pip install PyAstronomy
<!--pip install pyspectral-->

### Tool for MVIRI extraction of calibration coefficients  
#### Installation:  
cd $gitdir  
  git clone git@gitlab.eumetsat.int:frarue/sscc_calib_tool.git  
  cd $gitdir/sscc_calib_tool  
    python setup.py build  
    python setup.py install --prefix=$PREF  
    cd ..  
  python -c "import calib"  
### Fast tools for MVIRI numbercrunching written in fortran  
#### Installation:  
cd $P27git  
  git clone git@gitlab.eumetsat.int:USC_Climate/mviri_cruncher.git  
  cd $P27git/mviri_cruncher  
    python setup.py build  
    python setup.py install --prefix=${PREF}  
    cd ..  
  python -c "import cruncher; print dir(cruncher)"  
  python -c "import cruncher; print cruncher.vza.__doc__"  
### Tools for conversion of MFG orbit coordinates written in fortran  
#### Installation:  
cd $P27git  
  git clone git@gitlab.eumetsat.int:USC_Climate/mfg_orbit.git  
  cd $P27git/mfg_orbit  
    python setup.py build  
    python setup.py install --prefix=${PREF}  
    cd ..  
  python -c "import orbit2eff; print dir(orbit2eff)"  
  python -c "import orbit2eff; print orbit2eff.mgf2eff.__doc__"  
### Library for reading mviri L1.0 and 1.5 data  
#### Installation:  
cd $P27git
  git clone git@gitlab.eumetsat.int:USC_Climate/mviri_lib.git  
  cd $P27git/mviri_lib  
    python setup.py build  
    python setup.py install --prefix=${PREF}  
    cd ..  
  python -c "from mviri.mviri_l15 import read_rect2lp as rd; print dir(rd)"  
### FIDUCEO FCDR writer  
#### Test which version
python -c "from fiduceo.fcdr.writer.fcdr_writer import FCDRWriter; print FCDRWriter._version"  
#### Installation of working but outdated version:  
cd ${RICalPyROOT}/git/FCDRTools  
  python setup.py build  
  python setup.py install --prefix=${PREF}  
#### Installation of most recent version:  
#while developing:  
cd ${gitdir}/  
  git clone https://github.com/FIDUCEO/FCDRTools.git  
  cd FCDRTools  
    python setup.py build  
    python setup.py install --prefix=${PREF}  
#once development is ready:  
cd $P27git  
  git clone https://github.com/FIDUCEO/FCDRTools.git  
  cd FCDRTools  
    python setup.py build  
    python setup.py install --prefix=${PREF}  
### allantools  
cd $P27git  
 git clone https://github.com/aewallin/allantools  
 cd allantools/  
   python setup.py build  
   python setup.py install --prefix=${PREF}  
 cd..
 python -c "import allantools; print allantools.__version__"
### cython 
cd $P27git  
 git clone https://github.com/cython/cython.git  
 cd cython  
   python setup.py build  
   python setup.py install --prefix=${PREF}  
 cd ..  
 python -c "import cython; print cython.__version__"  
### netCDF4  
cd $P27git  
  git clone https://github.com/Unidata/netcdf4-python.git  
  cd netcdf4-python/  
    python setup.py build  
    python setup.py install --prefix=${PREF}   
  cd ..  
  python -c "import netCDF4; print netCDF4.__version__"  
### openCV  
<!--pip install cmake  
cd $P27git  
  git clone https://github.com/opencv/opencv.git  
  cd opencv/  
    mkdir build  
    cd build  
      cmake -D CMAKE_BUILD_TYPE=RELEASE -D CMAKE_INSTALL_PREFIX=${PREF} ..  
      make  
      make install-->  
  pip install opencv-python  
  python -c "import cv2; print cv2.__version__"  
  
### postgresql  
#cd $gitdir  
  #git clone https://github.com/postgres/postgres.git  
  #cd postgres  
    #./configure --prefix=$PREF --without-readline  
    #make  
    #make install  
cd $tarballs  
  wget https://ftp.postgresql.org/pub/source/v10.0/postgresql-10.0.tar.gz  
  tar -xzf postgresql-10.0.tar.gz   
  cd $tarballs/postgresql-10.0  
    ./configure --prefix=$PREF --with-python PYTHON=${PYTHON_DIR}/bin/python --without-readline  
    make j=12  
    make install  
    #best to start on cdrlnxv01  
    initdb -D ${PGDATA}  
    psql -h cdrlnxv01 -l#shows available db's on server  
    #to query content:  
    #psql -h cdrlnxv01 ricalpy_3_1  
    #SELECT * FROM logtable WHERE 1=0;  
    #COPY ( SELECT * FROM logtable WHERE idn>7200611010000) TO '/tcenas/home/frankr/forAlessio/MET7_IODC_db.csv' WITH CSV DELIMITER ',';  
    ##SELECT input_l10,input_l15,status_easy FROM logtable WHERE idn > 5200500000000 AND idn < 5200600000000;  
    #COPY ( SELECT * FROM logtable WHERE input_l10 = 1 AND input_l15 = -1) TO '/tcenas/home/frankr/QC/rect2lp_missing.csv' WITH CSV DELIMITER ',';  
    #COPY ( SELECT * FROM logtable WHERE idn>7000000000000 AND idn<8000000000000 AND status_easy=0 AND status_full=0 AND NOT problem_easy LIKE 'missing') TO '/tcenas/home/frankr/QC/fcdr_pending.csv' WITH CSV DELIMITER ',';  
    #COPY ( SELECT * FROM logtable WHERE idn>5000000000000 AND idn<6000000000000 AND status_easy=0 and input_l10=1 and  input_l15=1) TO '/tcenas/home/frankr/QC/fcdr_pendingbutpossible.csv' WITH CSV DELIMITER ',';  
    #DELETE FROM logtable WHERE idn>=6000000000000 and idn<7000000000000 AND status_easy=0;  
    #DELETE FROM logtable WHERE idn>7201612312359;  
    #DELETE FROM logtable WHERE idn>=2000000000000 and idn<3000000000000;  
    #DELETE FROM logtable WHERE idn>=3000000000000 and idn<4000000000000;  
    #DELETE FROM logtable WHERE idn>=2198212000000 and idn<2198301240000;  
    #DELETE FROM logtable WHERE idn>=2198710280100 and idn<21987 10280130;  
    #SELECT idn, problem_easy from logtable WHERE idn>=6000000000000 and idn<7000000000000 and status_easy=0 and input_l10=1 and  input_l15=1 and not problem_easy LIKE '%no SRF%';  
    #DELETE FROM logtable WHERE idn>=2000000000000 and idn<3000000000000 and (status_full <= 0 or status_easy<= 0);  
#while developing  
cd ${RICalPyCODE}  

## Configure postgres server to accept remote connections  
find ${PGDATA} -name "postgresql.conf"  
pg_ctl -D ${PGDATA} stop  
gedit /DSNNAS/Repro_Temp/users/frankr/CDR/data/postgresql.conf  
#in this file change  
#listen_addresses = 'localhost'  
#to  
#listen_addresses = '*'  
#and increase number of allowed connections to 600  
gedit /DSNNAS/Repro_Temp/users/frankr/CDR/data/pg_hba.conf  
#in this add to the end:  
#host    all             all              0.0.0.0/0                       trust  
#host    all             all              ::/0                            trust  
pg_ctl -D ${PGDATA} start  
netstat -na|grep 5432  

### Initialize database  
source ${gitdir}/RICalPy/environment.sh  
cd ${RICalPyROOT}/tools  
python generate_loggingDB.py  
#follow instructions  

### Initialize correlation file  
#copy over  
#e.g.:  
cp ${RICalPyROOT}/config/effect_correlation_MET7_3.1.csv ${RICalPyROOT}/config/effect_correlation_MET2_3.1.csv  
#or:  
#run /tcenas/home/frankr/GSCC/pbsify/pbsify.py  
#follow instructions to create MC result  
#make sure all old calibration and monte-carlo files are deleted  
#and then:  
cd ${RICalPyROOT}/tools  
python cf_errors.py  
#follow instructions  
#then copy files from $CALIB_STORE to the RICalPy/config folder:
cp ${CALIB_STORE}/M2*_R3.1* ${gitdir}/RICalPy/config/

### Create calibration files for all gains  
#This ensures a flawless processing of the software.  
#Although RICalPy would generate missing calibration files on the fly.  
#-->refer to sscc_calib_tool README.md
#basically the files will be generated in $CALIB_STORE and may have 
#to be copied over to RICalPy/config
#e.g. with:
cp ${gitdir}/sscc_calib_tool/config/M2_*_R3.1* ${gitdir}/RICalPy/config/
### Initialize LUT + static file  
cd ${RICalPyROOT}/tools  
python generate_static_and_LUT.py  
#follow instructions  

## RUNNING  

### on cdrlnxv01:  
#to allow logging  
pg_ctl -D ${PGDATA} start  
### on any machine:
source ${gitdir}/RICalPy/environment.sh  
cd ${RICalPyROOT}/src  
python ricals_main.py -s 2005180 -e 2005180 -m "MET7" -t 1200 -r 2.7 -v 1801 -a /DSNNAS/Repro_Temp/mviri_dev/FCDR/ -c "preliminary files only for illustration"  
### on cdrlnxv01:
#to stop logging db
pg_ctl -D ${PGDATA} stop  

## PBS  
source ${gitdir}/RICalPy/environment.sh  
cd ${RICalPyROOT}/pbsify
python pbsify.py
#follow instructions
