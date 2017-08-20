# Routine to:
# 1. read a constant solar spectrum (SS) with default taken from eumetsat inhouse standard (EUM/MTG/DEF/10/0611)
# 2. read a spectral response function (SRF)
# 3. interpolate SS to SRF bins 
# 4. convolute both
#
#Author: Frank Ruethrich
import numpy as np
from ephem import *
from datetime import timedelta
from datetime import date
import pandas as pd

import ricals_tools as to
import ricals_effects as ef

#debugging possibilities
global debug
debug=0
if debug==1:
  import matplotlib.pyplot as plt

def band_irrad(sundist,sdu,srffile="../sensor_data/pre_launch/spectral_response_MET_7_VIS.dat",ssfile="../config/solirr_yves.txt",d=0):
  """
  The default spectrum is the one described in EUM/MTG/DEF/10/0611
  """
  AU=149597870700     #[m]
  sunradius=695.700*1000 #[m]
  sunarea=4*np.pi*sunradius**2
  
  # 1. read solar spectrum (SS) [W m^-2 micron^-1]
  ssf=np.loadtxt(ssfile, skiprows=1, delimiter=';') #http://penandpants.com/2012/03/09/reading-text-tables-with-python/
  # 2.b read a spectral response function (SRF)
  try:
    srf=np.loadtxt(srffile, skiprows=2)
  except ValueError:
    srf=np.loadtxt(srffile, skiprows=2, delimiter=';')

  mini = np.amin(srf[:,0])
  maxi = np.amax(srf[:,0])
  binN = len(srf[:,0])
  step = step=((maxi-mini)/(binN-1))

  # 3. interpolate SS to SRF bins [W m^-2 micron^-1]
    #interpolate between selected data points
  #ssi=np.interp(srf[:,0],ss[:,0],ss[:,1]) #http://docs.scipy.org/doc/numpy-1.10.1/reference/generated/numpy.interp.html
    #interpolate in groups and means
  bint = np.linspace(mini-(step/2), maxi+(step/2), binN+1)
  bins = np.linspace(mini, maxi, binN)
  df1=pd.DataFrame(ssf[:,0:2],columns=['wl','rad'])
  df1=df1.loc[(df1.wl>=mini) & (df1.wl<=maxi)]
  ssi = interpolate_ssi(srf,df1,bint,binN)
  
  # 4. convolute both 
  #weighting
  f  = ssi*srf[:,1]#*np.median(np.diff(srf[:,0]))
  #df = ssi*srf[:,2]#*np.median(np.diff(srf[:,0]))
  df = ssi*srf[:,3:]*ssi[:, None]#includes covariances

  if debug==1:
    plt.imshow(srf[:,3:])
    plt.title("COV MATR")
    plt.colorbar()
    plt.show()
    plt.imshow(df)
    plt.title("R*COV MATR*R.T")
    plt.colorbar()
    plt.show()

  #integrating [W m^-2]
  F  = np.trapz(f,x=srf[:,0])#np.sum(f)
  #dF = np.trapz(df,x=srf[:,0])#np.sum(df)
  #dF=np.sqrt(to.int_2D(srf[:,3:],srf[:,0],srf[:,0]))#includes covariances
  dF=np.sqrt(to.int_2D(df,srf[:,0],srf[:,0]))#includes covariances

  #debugging
  if debug==1:
    #get sscc spectrum
    ssccspec="../config/solirr_sscc.txt"
    ssf_sscc=np.loadtxt(ssccspec, skiprows=1, delimiter=';')
    df1_sscc=pd.DataFrame(ssf_sscc[:,0:2],columns=['wl','rad'])
    df1_sscc=df1_sscc.loc[(df1_sscc.wl>=mini) & (df1_sscc.wl<=maxi)]
    ssi_sscc = interpolate_ssi(srf,df1_sscc,bint,binN)
    f_sscc  = ssi_sscc*srf[:,1]
    F_sscc  = np.trapz(f_sscc,x=srf[:,0])
    df_sscc = ssi*srf[:,3:]*ssi_sscc[:, None]
    dF_sscc=np.sqrt(to.int_2D(df_sscc,srf[:,0],srf[:,0]))
    #plot spectra
    plt.plot(df1.wl,df1.rad,"-",label="full sol spec (in use)")
    plt.plot(df1_sscc.wl,df1_sscc.rad,"-",label="full sol spec (SSCC)")
    plt.plot(srf[:,0],ssi,"-",color="red",label="sol spec at "+str(len(ssi)))
    plt.legend()
    plt.title("Solspec in use [1011]: F="+str(F)+" dF="+str(dF)+" Solspec SSCC ["+str(len(df1_sscc.wl))+"]: F="+str(F_sscc)+" dF="+str(dF_sscc))
    plt.show()
    #plot differences
    yvals=ssi-ssi_sscc
    plt.fill_between(srf[:,0], 0,yvals,where=yvals > 0,facecolor='green')
    plt.fill_between(srf[:,0], 0,yvals,where=yvals < 0,facecolor='red')
    plt.plot(srf[:,0],ssi-ssi_sscc,"-",color="black",label="sol spec diff at "+str(len(ssi)))
    plt.legend()
    plt.title("Solspec difference: in use - SSCC")
    plt.show()
    #plot for different count values
    counts=range(10,250)
    uncert=[]
    a1  = 0.0124
    a0  = 0.98
    YSL = 9
    SZA = 30
    Cs  = 5
    for C in counts:
      ref =ef.reflectance(C,a1*YSL+a0,Cs,SZA,1,F,1)
      uref=abs(F-F_sscc)*abs(ef.sensi_p("E",C,Cs,YSL,a0,a1,SZA,F,1))
      uncert.append((uref/ref)*100)
    plt.plot(counts,np.array(uncert),".")
    plt.show()
    #plot for different SZA
    C=100
    SZAs=range(0,90)
    uncert=[]
    for SZA in SZAs:
      ref =ef.reflectance(C,a1*YSL+a0,Cs,SZA,1,F,1)
      uref=abs(F-F_sscc)*abs(ef.sensi_p("E",C,Cs,YSL,a0,a1,SZA,F,1))
      uncert.append((uref/ref)*100)
    plt.plot(SZAs,np.array(uncert),".")
    plt.show()
    
    print "Irradiance at "+str(len(ssi))+" bins = "+str(F)+" Wm^-2"
    print "Irradiance uncertainty at "+str(len(ssi))+" bins = "+str(dF)+" Wm^-2"
    
  if d==2:#consider sun distance
    irr_at_sun=F*AU**2/sunarea #[Wm^-2sr^-1] #https://www.google.de/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&ved=0ahUKEwjV3trwy-LRAhWKJcAKHZ12AaQQFggmMAA&url=http%3A%2F%2Fwww.ghinstruments.com%2Fwp-content%2Fuploads%2F2015%2F03%2FP23_RADIANCE-TO-IRRADIANCE-CONVERSION.ppt&usg=AFQjCNFzDeGH5WYCeKizTMj5VfigqoNAYQ&cad=rja
    dirr_at_sun=dF*AU**2/sunarea
    F=irr_at_sun*sunarea/(sundist*AU)**2 #[W m^-2]
    dF=dirr_at_sun*sunarea/(sundist*AU)**2
    
  return (F,dF,srf)


def sunearth(year,doy):
    #eccentricity=0.01672
    #sidereal=360/365.256363
    #F = (1-eccentricity*np.cos( np.deg2rad(sidereal*(doy-4)) )) #approximation
    d=date(year,1,1)+timedelta(doy-1)
    s = Sun(d)
    s.compute(d)
    F=s.earth_distance
    dF=0
    return F, dF
  
def irrad(year,doy,srffile="../sensor_data/pre_launch/spectral_response_MET_7_VIS.dat"):
    #estimate some stuff
    sundist,sdu= sunearth(year,doy)
    (irr,irru,srf) = band_irrad(sundist,sdu,srffile)
    irru=np.sqrt(np.square(irr*0.001)+np.square(irru))
    return sundist,sdu,irr,irru,srf
  
def interpolate_ssi(srf,df1,bint,binN):
  if len(df1)<binN:
    bint=srf[:,0]
    ss=np.vstack((bint,np.interp(bint,df1.wl,df1.rad))).T
    df2=pd.DataFrame(ss[:,0:2],columns=['wl','rad'])
    ssi=np.array(df2.rad)
  else:
    bin_means = (np.histogram(df1.wl, bint.T, weights=df1.rad)[0] /
                np.histogram(df1.wl, bint.T)[0])
    mask = np.isnan(bin_means)
    bin_means[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), bin_means[~mask])
    ss=np.vstack((bins,bin_means)).T
    df2=pd.DataFrame(ss,columns=['wl','rad'])
    #groups  = df1.groupby(np.digitize(df1.wl, bint))
    #df2=pd.DataFrame(groups.mean())
    ssi=np.array(df2.rad)
  return ssi