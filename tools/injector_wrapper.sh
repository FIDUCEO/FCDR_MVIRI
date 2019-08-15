#!/bin/sh
#$1
#$2
source /tcenas/home/frankr/git/RICalPy/environment.sh
cd /tcenas/home/frankr/git/RICalPy/tools
echo python FCDR_injector.py 3.1 $1 $2
python FCDR_injector.py 3.1 $1 $2
