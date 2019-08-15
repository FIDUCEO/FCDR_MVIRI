#!/bin/sh
#module purge  
source ~/.climproc/modules/00-generic_setup_with_strand7.module
### ENVIRONMENT VARIABLES  
export PREF="/DSNNAS/Repro_Temp/users/frankr/CDR"  
echo $PREF
export RICalPyROOT=${PREF}/git/RICalPy  
#not on DSNNAS to be easily editable from X2go:  
export RICalPyCODE="/tcenas/home/frankr/git/RICalPy"  
#to install from git:  
export gitdir=${PREF}"/git/"  
export P27git=${PREF}"/git/python2.7"  
#to install from tarballs:  
tarballs=${PREF}"/tarballs/"  

export CALIB_STORE=${RICalPyROOT}
export CONFIGPP=${P27git}/mviri_lib/mviri/config  
export PGDATA=${PREF}/data  
export PATH=$PATH:${PREF}/bin  
export PATH=$PATH:${PREF}/lib  
export PYTHONPATH=${PREF}/lib64/python2.7/site-packages/  
export PYTHONPATH=${PYTHONPATH}:${PREF}/lib/python2.7/site-packages/  
export PYTHON_DIR=${PREF}/pythons/RICalPy  
export LD_LIBRARY_PATH=${PREF}/lib/  
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$PYTHON_DIR/include/python2.7/
source /DSNNAS/Repro_Temp/users/frankr/anaconda2/etc/profile.d/conda.sh
#source /DSNNAS/Repro_Temp/users/mgrant/anaconda2/etc/profile.d/conda.sh
conda activate $PYTHON_DIR
# module use -a /DSNNAS/Repro_Temp/users/mgrant/modulefiles
# module load anaconda2/4.3.1
# source activate $PYTHON_DIR



cd ${RICalPyROOT}/src
