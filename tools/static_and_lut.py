# -*- coding utf-8 -*- 
"""
this file contains the main program
"""

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

import ricals_netcdf as wr
import ricals_effects as ef
import ricals_tools as to
import ricals_calib as ca
import ricals_irrad as ir
import read_rect2lp as rd
import read_imag2tg as r0

from fiduceo.fcdr.writer.fcdr_writer import FCDRWriter
from fiduceo.fcdr.writer.default_data import DefaultData

sys.path.insert(0, os.getcwd()[:-4]+'/lib/nrCrunch/')

#import cruncher_f as cr #grid definition as implemented in fortran routines in mfg user manual
#global grid
#grid="mfg_imp"

import cruncher as cr #grid definition as described in mfg user manual
global grid
grid="mfg_def"

global Software
global version
global parallel #flag to switch on/off parallel processing of multiple years
global N
global early
global inflated
global tstsite

#default globals
Software = "2.2"
N        = 5000
early    = ['M2','M3','P2']
parallel = "n"
save     = 1
inflated = 1
#tstsite   = [17.96,-18.84] #North atlantic; lat,lon
#tstsite   = [38.21,-27.42] #azores; lat,lon
#tstsite   = [23.8,-0.40]   #alg1; lat,lon
tstsite   = [28.55,23.39]  #lyb4; lat,lon
#tstsite = [0.,0.]
rel="2.5"
srf_vers="0954"
a="/DSNNAS/Repro_Temp/mviri_dev/FCDR/"
versionfolder    = "MFG15_FCDR_V"+rel+"_SRF"+srf_vers
s_Year="XXXX"
s_Month="XX"
s_Day="XX"
s_Time="XXXX"
comment="pre-beta release!! Use only for illustration; acqtime fixed; SZA fixed; solar spectrum from rayference"

staticdir  =  a + versionfolder+'/aux/'
lutdir     =  a + versionfolder+'/aux/'

for satellite,nomlon in zip(["MET7","MET5","MET5"],[57.0,0.0,63.0]):
  #lat, lon, VZA, VAZ and static+++++++++
  print "create static fcdr"
  staticname=staticdir+wr.ncName(s_Year,s_Month,s_Day,s_Time,"STATIC",nomlon,Software,rel,satellite)
  LAT,LON   = cr.latlon(N)
  LON       = np.add(LON,nomlon)
  LAT2,LON2 = cr.latlon(N/2)
  LON2      = np.add(LON2,nomlon)
  m1                       =  np.ones((2500,2500))
  m1[(LAT2<-90)|(LAT2>90)] =  0
  m2                       =  np.ones((5000,5000))
  m2[ (LAT<-90)|(LAT>90) ] =  0
  if not os.path.exists(staticname):
    print "      ...writing: "+staticname
    if not os.path.exists(staticdir):
      os.makedirs(staticdir)
    if float(rel)<2.0:
      varnames=["latitude_vis","longitude_vis","latitude_ir_wv","longitude_ir_wv"]
    else:
      varnames=["latitude_vis","longitude_vis","latitude_ir_wv","longitude_ir_wv","viewing_geometry"]
    vals=[LAT,LON,LAT2,LON2,0]
    wr.save_static(staticname,nomlon,varnames,vals,m1,comment,rel)
  #sensitivities and LUT+++++++++
  print "create SZA sensitivities LUT"
  lutname=lutdir+wr.ncName(s_Year,s_Month,s_Day,s_Time,"LUT",nomlon,Software,rel,satellite)
  if (not os.path.exists(lutname)): #create LUT file only if t is not set
    if not os.path.exists(staticdir):
      os.makedirs(lutdir)
    varnames=['s_sza_latitude_vis','s_sza_longitude_vis','s_sza_time']
    vals=[LAT,LON,0]
    wr.save_static(lutname,nomlon,varnames,vals,m1,comment,rel)
exit()