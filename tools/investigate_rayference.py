import sys, getopt
import numpy as np
import netCDF4 as nc
import time
import timeit
import datetime as dt
from ephem import *
import csv
import os
import glob
import shutil
#from sympy import mpmath as sm
import matplotlib.pyplot as plt
import functools
#import psycopg2
import cruncher as cr #grid definition as described in mfg user manual
sys.path.insert(0, '../src/')
import ricals_netcdf as wr
import ricals_effects as ef
import ricals_tools as to
import ricals_calib as ca
import ricals_irrad as ir
from mviri.mviri_l10 import read_imag2tg as r0
from mviri.mviri_l15 import read_rect2lp as rd
from fiduceo.fcdr.writer.fcdr_writer import FCDRWriter
from fiduceo.fcdr.writer.default_data import DefaultData
  
#def per_sr(R,year,doy,dist_to_source_m,rad_of_source_m=6.96342*10**8):# R:[W m^-2 or W m^-2 micron^-1]
  #planeangle_rad=2*np.sin(rad_of_source_m/dist_to_source_m)#[radians] sin alpha=gegenkathete/hypothenuse
  #solidangle_sr=2*np.pi * (1 - np.cos(planeangle_rad/2)) # [steradians]http://smart-unit-converter.com/angle-converter.php
  #temp=np.divide(R,solidangle_sr)#[W m^-2 sr^-1]
  #return temp

a="/DSNNAS/Repro_Temp/mviri_dev/FCDR/"
rel="2.6"
srf_vers="1000"
versionfolder    = "MFG15_FCDR_V"+rel+"_SRF"+srf_vers
runfolder=a+versionfolder
satellite="MET7"
year=2006
DOY=30
timestring="200601301200"
logs=logs=to.logger(rel)
debug=1
srf=to.makesrf(runfolder,satellite,year,DOY,srf_vers,timestring,logs,debug)
srf.srfonfly(0)
#sundist,sdu,irr,uirr,srf=ir.irrad(int(timestring[:4]),DOY,srffile=srf.expected)
ssfile="../config/solirr_yves.txt"
ssf=np.loadtxt(ssfile, skiprows=1, delimiter=';')
srffile=srf.expected
srf=np.loadtxt(srffile, skiprows=2)
mini = np.amin(srf[:,0])
maxi = np.amax(srf[:,0])
binN = len(srf[:,0])
step = step=((maxi-mini)/(binN-1))
bint = np.linspace(mini-(step/2), maxi+(step/2), binN+1)
bins = np.linspace(mini, maxi, binN)
import pandas as pd
df1=pd.DataFrame(ssf[:,0:2],columns=['wl','rad'])
df1=df1.loc[(df1.wl>=mini) & (df1.wl<=maxi)]
ssi = ir.interpolate_ssi(srf,df1,bint,binN)#[W m^-2 micron^-1]
AU=149597870700     #[m]
sundist,sdu= ir.sunearth(year,DOY)#[AU,AU]
sunradius=6.96342*10**8 #[m]
ssi_per_sr=ir.per_sr(ssi,year,DOY,sundist*AU,sunradius)

def rayference(sza,ssi):
  return np.divide(np.multiply(ssi,np.cos(np.radians(sza))),np.pi)

sza=60.08
ssi_sza=rayference(sza,ssi)
plt.plot(srf[:,0],srf[:,1]*2000,label="SRF [1 micron^-1; scaled by 2000]")
plt.plot(ssf[:,0],ssf[:,1],label="original [W m^-2 micron^-1]")
plt.plot(srf[:,0],ssi,label="interpolated [W m^-2 micron^-1]")
plt.plot(srf[:,0],ssi_sza,label="at SZA of "+str(sza)+" [W m^-2 micron^-1]")
plt.plot(srf[:,0],ssi_per_sr*10**-5,label="per sr [W m^-2 micron^-1 sr^-1; scaled by 10^-5]")
plt.legend()
plt.show()
f  = ssi*srf[:,1]
df = ssi*srf[:,3:]*ssi[:, None]
f_ray  = ssi_sza*srf[:,1]
df_ray=ssi_sza*srf[:,3:]*ssi_sza[:, None]

plt.imshow(srf[:,3:])
plt.title("COV MATR")
plt.colorbar()
plt.show()
plt.imshow(df)
plt.title("R*COV MATR*R.T")
plt.colorbar()
plt.show()
plt.imshow(df_ray)
plt.title("((cos(sza)*R)/pi)*COV MATR*((cos(sza)*R)/pi).T")
plt.colorbar()
plt.show()
plt.imshow(df-df_ray)
plt.title("diff")
plt.colorbar()
plt.show()

F  = np.trapz(f,x=srf[:,0])
F_ray  = np.trapz(f_ray,x=srf[:,0])
F_test=F
meteosatdist=36000000.
earthradius=12742000.
F_test_sr=ir.per_sr(F_test,year,DOY,meteosatdist,earthradius)

dF=np.sqrt(to.int_2D(df,srf[:,0],srf[:,0]))
dF_ray=np.sqrt(to.int_2D(df_ray,srf[:,0],srf[:,0]))

F_per_sr=ir.per_sr(F,year,DOY)
dF_per_sr=ir.per_sr(dF,year,DOY)

F_ray_per_sr=ir.per_sr(F_ray,year,DOY)
dF_ray_per_sr=ir.per_sr(dF_ray,year,DOY)

