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
Software = "2.4"
N        = 5000
early    = ['M2','M3','P2']
parallel = "n"
save     = 0
inflated = 1

timing=[]

LAT,LON   = cr.latlon(N)
LAT2,LON2 = cr.latlon(N/2)
m1             = np.ones((N/2,N/2))
m1[ (LAT2<-90)|(LAT2>90) ]  = 0
m2             = np.ones((N,N))
m2[ (LAT<-90)|(LAT>90) ]  = 0
t1=np.ones((N/2,N/2))*9.5

timing.append(timeit.default_timer())
SZA,SAZ=cr.sza(t1,122,np.array(m1),LAT,LON)
timing.append(timeit.default_timer())

R=10
SZA_tp=SZA[::R,::R]
LAT_tp=LAT[::R,::R]
LON_tp=LON[::R,::R]

N_tp=N/R
X_tp=np.arange(N_tp)
Y_tp=np.arange(N_tp)
X=np.arange(N)
Y=np.arange(N)
  
from scipy.interpolate import griddata
from scipy.interpolate import RectBivariateSpline

points=np.column_stack((np.ndarray.flatten(LON_tp),np.ndarray.flatten(LAT_tp)))

timing.append(timeit.default_timer())

SZA_cubic=griddata(points, np.ndarray.flatten(SZA_tp),(LON,LAT),method='cubic')
timing.append(timeit.default_timer())

SZA_lin=griddata(points, np.ndarray.flatten(SZA_tp),(LON,LAT),method='linear')
timing.append(timeit.default_timer())

spl=RectBivariateSpline(X_tp,Y_tp, SZA_tp)
xpoints, ypoints = np.meshgrid(X,Y)
SZA_spli=spl.ev(xpoints.ravel(), ypoints.ravel()).reshape(xpoints.shape).T
timing.append(timeit.default_timer())

print "calculation took   : "+str(timing[1]-timing[0])
print "cubic interpol took: "+str(timing[3]-timing[2])
print "linear interpol took: "+str(timing[4]-timing[3])
print "bivspline interpol took: "+str(timing[5]-timing[4])

SZA[m2==0]=np.nan
SZA_cubic[m2==0]=np.nan
SZA_lin[m2==0]=np.nan
SZA_spli[m2==0]=np.nan

plt.imshow(SZA,clim=(0,120))
plt.title("SZA")
plt.colorbar()
plt.show()

#absolute
Err_cub=SZA-SZA_cubic
Err_lin=SZA-SZA_lin
Err_spli=SZA-SZA_spli

plt.imshow(Err_cub,clim=(-0.05,0.05))
plt.title("Abs Err SZA cubic interpol from every "+str(R))
plt.colorbar()
plt.show()
plt.imshow(Err_lin,clim=(-0.05,0.05))
plt.title("Abs Err SZA lin interpol from every "+str(R))
plt.colorbar()
plt.show()
plt.imshow(Err_spli,clim=(-0.05,0.05))
plt.title("Abs Err SZA biv interpol from every "+str(R))
plt.colorbar()
plt.show()
#relative
Err_rel_cub=(SZA-SZA_cubic)/SZA
Err_rel_lin=(SZA-SZA_lin)/SZA
plt.imshow((SZA-SZA_cubic)/SZA,clim=(-0.01,0.01))
plt.title("Rel Err SZA cubic interpol from every "+str(R))
plt.colorbar()
plt.show()
plt.imshow((SZA-SZA_lin)/SZA,clim=(-0.01,0.01))
plt.title("Rel Err SZA lin interpol from every "+str(R))
plt.colorbar()
plt.show()
plt.imshow((SZA-SZA_spli)/SZA,clim=(-0.01,0.01))
plt.title("Rel Err SZA biv interpol from every "+str(R))
plt.colorbar()
plt.show()
