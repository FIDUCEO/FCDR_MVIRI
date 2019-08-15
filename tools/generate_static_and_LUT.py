import sys, getopt
import numpy as np
import netCDF4 as nc
import time
import timeit
import datetime as dt
import csv
import os
import glob
import shutil
#from sympy import mpmath as sm
import matplotlib.pyplot as plt
import functools

from functools import partial

#import psycopg2

import cruncher as cr #grid definition as described in mfg user manual

sys.path.insert(0, '../src/')
import ricals_netcdf as wr
import ricals_effects as ef
import ricals_tools as to
#import ricals_calib as ca
import calib as ca
import ricals_irrad as ir

from mviri.mviri_l10 import read_imag2tg as r0
from mviri.mviri_l15 import read_rect2lp as rd

from fiduceo.fcdr.writer.fcdr_writer import FCDRWriter
from fiduceo.common.writer.default_data import DefaultData

import calib  
debug=0


satellite= raw_input("Which satellite (e.g. MET5)?")
nomlon   = float(raw_input("at which nominal longitude?"))
rel      = raw_input("Create files for which file release? (e.g. 2.5) ")
Software = raw_input("Create files for which software release? (e.g. 2.2) ")
comment  = raw_input("Any comment to store in the file (e.g. first release - use with caution!)?")
a        = raw_input("Destination folder (e.g.: /DSNNAS/Repro/mviri/level1/MFG15_FCDR_V1/)")

#a        = "/DSNNAS/Repro/mviri/level1/MFG15_FCDR_V1/"
#a        = "/DSNNAS/Repro_Temp/mviri_dev/FCDR/"

srf_vers = "1000"
versionfolder    = ""#"MFG15_FCDR_V"+rel+"_SRF"+srf_vers
staticdir        =  a + versionfolder+'/aux/'
lutdir           =  a + versionfolder+'/aux/'

calli=calib.calib(satellite,rel)

N=5000

LAT,LON   = cr.latlon(N)
LON       = np.add(LON,nomlon)
LAT2,LON2 = cr.latlon(N/2)
LON2     = np.add(LON2,nomlon)
m1             = np.ones((2500,2500))
m1[(LAT2<-90)|(LAT2>90)] = 0
m2             = np.ones((5000,5000))
m2[ (LAT<-90)|(LAT>90) ]  = 0

s_Year  = calli.launchdate()[:4]
s_Month = calli.launchdate()[4:6]
s_Day   = calli.launchdate()[6:8]
s_Time  = "1200"

print "   ...first have to prepare static fcdr"
staticname=staticdir+wr.ncName(s_Year,s_Month,s_Day,s_Time,"STATIC",nomlon,Software,rel,satellite)

varnames=["latitude_vis","longitude_vis","latitude_ir_wv","longitude_ir_wv"]
vals=[LAT,LON,LAT2,LON2]
wr.save_static(staticname,nomlon,varnames,vals,m1,comment,rel,debug=debug)


print "      ...now have to prepare LUT file"
lutname=lutdir+wr.ncName(s_Year,s_Month,s_Day,s_Time,"LUT",nomlon,Software,rel,satellite)

varnames=['s_sza_latitude_vis','s_sza_longitude_vis','s_sza_time']
vals=[LAT,LON,0]
wr.save_static(lutname,nomlon,varnames,vals,m1,comment,rel,debug=debug)

