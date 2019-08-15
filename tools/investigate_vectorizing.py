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

from functools import partial

timing=[]
corr=np.array([[ 1.  ,  0.99,  0.9 ,  0.  ,  0.  ,  0.  ],\
               [ 0.99,  1.  ,  0.9 ,  0.  ,  0.  ,  0.  ],\
               [ 0.9 ,  0.9 ,  1.  ,  0.  ,  0.  ,  0.  ],\
               [ 0.  ,  0.  ,  0.  ,  1.  ,  0.  ,  0.  ],\
               [ 0.  ,  0.  ,  0.  ,  0.  ,  1.  ,  0.  ],\
               [ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  1.  ]])

#generate template array
N=10
E=6
temp=np.zeros((N,N,E))
#setup 3D arrays for uncertainties and sensitivities
collectu=np.copy(temp)
collects=np.copy(temp)
for efN,targ in zip([0,1,2,3,4,5],["E","a0","a1","Z","SZA","Cs"]):
  uval=np.random.rand(N,N)*(0.1*float(efN)+0.000001)
  collectu[:,:,efN]=uval
  sensitmp=np.random.rand(N,N)*(0.01*float(efN)+0.000001)
  collects[:,:,efN]=sensitmp

#watermark them
FV=-9
collectu[0,0,:]=FV
collectu[0,N-1,:]=FV
collectu[N-1,0,:]=FV
collectu[N-1,N-1,:]=FV
collectu[0,0,:]=FV
collectu[0,N-1,:]=FV
collectu[N-1,0,:]=FV
collectu[N-1,N-1,:]=FV
collects[0,0,:]=FV
collects[0,N-1,:]=FV
collects[N-1,0,:]=FV
collects[N-1,N-1,:]=FV
#broadcast multiply arrays
#https://stackoverflow.com/questions/32171917/copy-2d-array-into-3rd-dimension-n-times-python
cov=collectu[:,:,:,None]*corr[None,None,:,:]*collectu[:,:,None,:]
imp=collects[:,:,:,None]*cov*collects[:,:,None,:]
#convolute covariances
uncMS_syst=np.sqrt(imp.sum(axis=(2,3)))

#version that is easier to debug:
uncMS_syst2=np.zeros((N,N))
uncMS_syst3=np.zeros((N,N))
cov2=np.zeros((N,N,E,E))
cov3=np.zeros((N,N,E,E))
imp2=np.zeros((N,N,E,E))
imp3=np.zeros((N,N,E,E))
for x in range(N):
  print "line "+str(x)
  for y in range(N):
    uncMS_syst2[x,y],cov2[x,y,:,:],imp2[x,y,:,:]=ef.law_of_error_trivial(corr,collectu[x,y,:],collects[x,y,:])
    uncMS_syst3[x,y],cov3[x,y,:,:],imp3[x,y,:,:]=ef.law_of_error_normal(corr,collectu[x,y,:],collects[x,y,:])


#compare result:
print("both methods differ?")
for e in [0,1,2,3,4,5]:
  print ["E","a0","a1","Z","SZA","Cs"][e]
  print collectu[:,:,e].mean()
  plt.imshow(collectu[:,:,e])
  plt.show()

print("uncMS_syst")
plt.imshow(uncMS_syst)
plt.show()
print("uncMS_syst2")
plt.imshow(uncMS_syst2)
plt.show()
print("uncMS_syst3")
plt.imshow(uncMS_syst3)
plt.show()
print("uncMS_syst diff")
plt.imshow(uncMS_syst-uncMS_syst3)
plt.show()