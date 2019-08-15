#GSICS naming convention:
#W_JP-JMA-MSC,VIS+IR+SRF,Himawari8+AHI_C_RJTD.nc
#example:
#W_XX-EUMETSAT-DarmstadtVISIRSRFMET7MVIRI_C_EUMG.nc
#see: http://gsics.atmos.umd.edu/pub/Development/Srf4Giro/W_XX-EUMETSAT-DarmstadtVISIRSRFMET7MVIRI_C_EUMG.nc
import glob
from netCDF4 import Dataset
import pandas as pd
import numpy as np

import sys

sys.path.insert(0, '../src/')
import ricals_tools as to

ref="/tcenas/home/frankr/tmp/W_XX-EUMETSAT-DarmstadtVISIRSRFMET7MVIRI_C_EUMG.nc"
fr = Dataset(ref, mode='r')

srf=np.vstack((fr.variables["wavelength"][:,:].T,fr.variables["srf"][:,:].T)).T
df1=pd.DataFrame(srf,columns=['wl_old','wl_wv','wl_ir','resp_old','resp_wv','resp_ir'])
oril=len(df1)
for i in range(20):
  print len(df1)
  df1.loc[len(df1)] = [(df1.wl_old[len(df1)-1]+0.01),-9999.00,-9999.00,-9999.00,-9999.00,-9999.00]


def get_srf(srfnamesource="/DSNNAS/Repro/mviri/level1/MFG15_FCDR_V1/tmp/recovered_srf_v1000/srf_MET7_2005123_2005124_1000-RC00009_S10EE_10.dat",se="",headers=56):
  #read srf
  if se=="":
    allrows=np.loadtxt(srfnamesource,skiprows=headers)
  else:
    allrows=np.loadtxt(srfnamesource,skiprows=headers,delimiter=se)
  return allrows
def merge(df1,df2,col="resp_new",idl="wl_new",idr="wl_old"): #df1=spectrum df2=srf
    step=df2[idr][1]-df2[idr][0]
    print step
    tmp=(np.histogram(df1[idl]+(step/2.), df2[idr], weights=df1[col])[0] /
                np.histogram(df1[idl]+(step/2.), df2[idr])[0])
    tmpdf=pd.DataFrame({idr:df2[idr][:-1],col:tmp})
    d2=pd.merge(df2,tmpdf,on=idr)
    return d2

for f in glob.glob("/DSNNAS/Repro/mviri/level1/MFG15_FCDR_V1/tmp/FCDR_MVIRISRF.git/trunk/srf/srf_MET7_*_*_1801-RC00004_*_10/*MET7*.dat"):
  srf=get_srf(f)
  df2=pd.DataFrame(srf[:,0:2],columns=['wl_new','resp_new'])
  
  df3=merge(df2,df1)
  
  #import matplotlib.pyplot as plt
  #plt.plot(df2.wl_new,df2.resp_new)
  #plt.plot(df3.wl_old,df3.resp_new)
  #plt.show()
  #exit()
  
  fileout="/tcenas/home/frankr/tmp/W_XX-EUMETSAT-Darmstadt,VIS+IR+SRF,"+f[-46:-42]+"+MVIRI_C_EUMG_"+f[-41:-34]+".nc"
  out = Dataset(fileout, mode='w')

  for attname in fr.ncattrs():
    setattr(out,attname,getattr(fr,attname))

  for dimname,dim in fr.dimensions.iteritems():
    if not len(dim)==oril:
      out.createDimension(dimname,len(dim))
    else:
      out.createDimension(dimname,len(df3))
  out.close()

  for varname,ncvar in fr.variables.iteritems():
    print varname
    try:
      fv=getattr(ncvar,"_FillValue")
      fvflag=True
    except AttributeError:
      print "no FillValue"
      fvflag=False
    if fvflag:
      out = Dataset(fileout, mode='a')
      var = out.createVariable(varname,ncvar.dtype,ncvar.dimensions,fill_value=fv)
    else:
      out = Dataset(fileout, mode='a')
      print ncvar.dtype,ncvar.dimensions
      var = out.createVariable(varname,ncvar.dtype,ncvar.dimensions)

    for attname in ncvar.ncattrs():
      print "   "+attname
      if not "FillValue" in attname:
        setattr(var,attname,getattr(ncvar,attname))
    
    if varname=="wavelength":
      var[:,0]=df3.wl_old.values
      var[:,1]=df3.wl_wv.values
      var[:,2]=df3.wl_ir.values
    elif varname=="srf":
      var[:,0]=df3.resp_new.values
      var[:,1]=df3.resp_wv.values
      var[:,2]=df3.resp_ir.values
    elif varname=="wavenumber":
      var[:,0]=10000000 / (df3.wl_old.values *1000)
      var[:,1]=10000000 / (df3.wl_wv.values *1000)
      var[:,2]=10000000 / (df3.wl_ir.values *1000)
    else:
      var[:] = ncvar[:]
    out.close()
