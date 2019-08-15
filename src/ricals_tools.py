"""
this file contains tools for plotting etc
"""
import sys 
import datetime
import numpy as np
import netCDF4 as nc
import time
import timeit
import csv
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors
from matplotlib.ticker import LogFormatter
import math
from PyAstronomy import pyasl
from numba import jit

import functools
from functools import partial

import orbit2eff as ob


global lim
lim=[0,2500]
global debug
debug=0
#----------------------------------------------------------------------------
#numba functions
#----------------------------------------------------------------------------
#@jit
#def flagger(ARR,N):
  #ARRB=np.zeros((N,N),dtype=np.uint8)
  #for m in range(np.shape(ARR)[1]):
    #for n in range(np.shape(ARR)[2]):
      #string="".join(ARR[:,m,n])
      #ARRB[m,n]=int(string,2)
  #return ARRB#map(functools.partial(bitmod,X=BIT), ARR)

#@jit
#def sumit(vals,axis):
  #return np.sum(vals,axis=axis)

def law_of_error_normal(corr,uncert):
  #generate covariance matrix from uncertainties and correlations
  imp_mat=uncert*corr*uncert[:, None]
  #sum
  comb=imp_mat.sum(axis=(0,1))
  #root
  u_comb=np.sqrt(comb)
  return u_comb#,cov_mat,imp_mat
  
@jit
def broadcast(u,corr,N,E,out):
  for x in range(N):
    print x
    for y in range(N):
      out[x,y]=law_of_error_normal(corr,u[x,y,:])
  return

#----------------------------------------------------------------------------
#Print a 2D or 3D array
#----------------------------------------------------------------------------
def printimage(pix,tit,mini=0,maxi=255,save=0,cmap=plt.cm.gray,direction="c",m2="n"):
    if m2=="n":
      if  isinstance(pix,list):
        print "3d"
        print np.shape(pix[0])
        m2=np.ones_like(pix[0])
        for p in range(len(pix)):
          tmp=np.array(pix[p])
          tmp[m2==0]=1
          pix[p]=tmp
      else:
        print "2d"
        m2=np.ones_like(pix)
        pix[m2==0]=1
    print np.shape(pix)
    try: 
      dims=pix.ndim
      fig, ax = plt.subplots(figsize=(18,18))
      im=ax.imshow(np.fliplr(pix),origin='lower', clim=(mini, maxi), cmap=cmap, interpolation='nearest', extent=[lim[0]+np.shape(pix)[0],lim[0],lim[0],lim[0]+np.shape(pix)[0]])
      plt.colorbar(im,fraction=0.046, pad=0.04)
      plt.title(tit,fontsize=10)
      if save==0:
        plt.show()
      else:
        try:
          plt.savefig("../test/"+tit+".png", bbox_inches='tight')
          plt.close()
        except:
          plt.close()
    except AttributeError:
      name="_".join(tit)
      subplts=len(pix)
      if direction=="h":
        subpltsX=len(pix)
        subpltsY=1
      elif direction=="v":
        subpltsX=1
        subpltsY=len(pix)
      elif direction=="c":
        subpltsX=int(np.sqrt(len(pix))+0.5)
        subpltsY=subpltsX
      #fig, ax = plt.subplots(subplts,sharey=True,figsize=(6,6))
      fig, ax = plt.subplots(subpltsY,subpltsX,figsize=(36,36)) 
      print maxi
      if isinstance(maxi,(list, tuple, np.ndarray)):
        for idx,mi,ma in zip(range(subplts),mini,maxi):
          idxY=int(idx/subpltsX)
          idxX=idx-(idxY*subpltsX)
          print idxY,idxX
          print mi,ma,idx
          if ma<1E+20:
            im=ax[idxY,idxX].imshow(np.fliplr(pix[idx]),origin='lower', \
              clim=(mi, ma), cmap=cmap[idx], interpolation='nearest', \
                extent=[lim[0]+np.shape(pix[idx])[0],lim[0],lim[0],lim[0]+np.shape(pix[idx])[0]])
          else:
            print "plotting logscale 1"
            im=ax[idxY,idxX].imshow(np.fliplr(pix[idx]),origin='lower', \
              clim=(mi, ma), cmap=cmap[idx], interpolation='nearest', \
              extent=[lim[0]+np.shape(pix[idx])[0],lim[0],lim[0],lim[0]+np.shape(pix[idx])[0]], \
              norm=matplotlib.colors.LogNorm() )
          ax[idxY,idxX].set_title(tit[idx])
          if ma<1E+20:
            fig.colorbar(im,ax=ax[idxY,idxX],fraction=0.046, pad=0.04)
          else:
            formatter = LogFormatter(10, labelOnlyBase=False)
            fig.colorbar(im,ax=ax[idxY,idxX],fraction=0.046, pad=0.04, format=formatter)
            #cb=fig.colorbar(im,ax=ax[idxY,idxX],fraction=0.046, pad=0.04)#,ticks=[0,0.5*scale,1*scale,10*scale],ticklabels=["0","0.5","1","10"])
            #cb.set_ticks=([0,0.5*scale,1*scale,10*scale])
            #cb.set_ticklabels=(["0","0.5","1","10"])
      else:
        for idx in range(subplts):
          idxY=int(idx/subpltsX)
          idxX=idx-(idxY*subpltsX)
          print idxY,idxX
          if maxi<1E+20:
            try:
              im=ax[idxY,idxX].imshow(np.fliplr(pix[idx]),origin='lower',\
                clim=(mini, maxi), cmap=cmap, interpolation='nearest', \
                  extent=[lim[0]+np.shape(pix[idx])[0],lim[0],lim[0],lim[0]+np.shape(pix[idx])[0]])
            except IndexError:
              im=ax[idx].imshow(np.fliplr(pix[idx]),origin='lower',\
                clim=(mini, maxi), cmap=cmap, interpolation='nearest', \
                  extent=[lim[0]+np.shape(pix[idx])[0],lim[0],lim[0],lim[0]+np.shape(pix[idx])[0]])
          else:
            print "plotting logscale 2"
            im=ax[idxY,idxX].imshow(np.fliplr(pix[idx]),origin='lower',\
              clim=(mini, maxi), cmap=cmap, interpolation='nearest', \
              extent=[lim[0]+np.shape(pix[idx])[0],lim[0],lim[0],lim[0]+np.shape(pix[idx])[0]], \
              norm=matplotlib.colors.LogNorm() )
          try:
            ax[idxY,idxX].set_title(tit[idx])
          except IndexError:
            ax[idx].set_title(tit[idx])
        if maxi<1E+20:
          fig.colorbar(im, ax=ax.ravel().tolist())
        else:
          formatter = LogFormatter(10, labelOnlyBase=False) 
          fig.colorbar(im, ax=ax.ravel().tolist(),ticks=[mini+0.0001,1,10], format=formatter)
      if save==0:
        plt.show()
      else:
        try:
          plt.savefig("../test/"+name+".png", bbox_inches='tight')
          plt.close()
        except:
          plt.close()

#----------------------------------------------------------------------------
#Print 4 matrices
#----------------------------------------------------------------------------
def print4mat(pixs,tits,m2,save=0):
  subpltsX=2
  subpltsY=2
  subplts=subpltsX*subpltsY
  centered=['n','n','n','y']
  cmaps=['Greys_r','Blues','Reds','RdBu_r']
  #initialize figure
  fig, ax = plt.subplots(subpltsX,subpltsY,figsize=(18,18)) 
  for idx,pix,tit in zip(range(subplts),pixs,tits):
    #get absolute spread and set min/max
    pix[m2==0]=np.nan
    low=abs(np.nanpercentile(np.array(pix),5))
    up=abs(np.nanpercentile(np.array(pix),95))
    print low,up
    maxi=np.amax([low,up])
    mini=0
    if centered[idx]=='y':
      mini=maxi*-1
    print mini,maxi
    stepwidth=(mini-maxi)/11
    #define cmap
    cmap=plt.get_cmap(cmaps[idx], 11)
    #define axes object
    idxY=int(idx/subpltsX)
    idxX=idx-(idxY*subpltsX)
    #plot
    ax[idxY,idxX].set_axis_bgcolor('grey')
    im=ax[idxY,idxX].matshow(np.flipud(np.fliplr(pix)),\
                                          cmap=cmap, \
                                          vmin = mini-stepwidth/2, \
                                          vmax = maxi+stepwidth/2)
    #colorbar
    ticks=np.around(np.linspace(mini,maxi,11),decimals=4)
    ticklabels=ticks.astype(str)
    ticklabels[0]="<="+ticklabels[0]
    ticklabels[len(ticklabels)-1]=">="+np.around(ticks[len(ticklabels)-1],decimals=2).astype(str)
    print ticklabels
    cax = plt.colorbar(im,ax=ax[idxY,idxX],extend='both',fraction=0.046, pad=0.04)#,extent= [mini-stepwidth/2,maxi+stepwidth/2,mini,maxi])
    #cax = mpl.colorbar.ColorbarBase(ax=ax[idxY,idxX], 
                                #cmap=cmap,
                                #norm=mpl.colors.normalize(mini,maxi),
                                ## to use 'extend', you must
                                ## specify two extra boundaries:
                                #boundaries=mini-stepwidth/2 + ticks + maxi+stepwidth/2,
                                #extend='both',
                                #ticks=ticks)
    #cax.tick_params(labelsize=12)
    #title
    ax[idxY,idxX].set_title(tit,fontsize= 25)
    #append filename
  name="_".join(tits)
  if save==0:
    plt.show()
  else:
    try:
      plt.savefig("../test/"+name+".png", bbox_inches='tight')
      plt.close()
    except:
      plt.close()


#----------------------------------------------------------------------------
#Print a 2D array and quivers
#----------------------------------------------------------------------------
def printquiver(pix,px,py,valuesx,valuesy,tit,mini=0,maxi=255,save=0):
    #display image
    print np.shape(pix)
    fig, ax = plt.subplots(figsize=(18,18))
    im=ax.imshow(np.fliplr(pix),origin='lower', clim=(mini, maxi), cmap=plt.cm.gray, interpolation='nearest', extent=[np.shape(pix)[1],1,1,np.shape(pix)[0]])
    plt.colorbar(im,fraction=0.046, pad=0.04)
    Q = plt.quiver(px,py,valuesx,valuesy)
    plt.title(tit,fontsize=10)
    if save == 0:
      plt.show()
    else:
      plt.savefig("../test/"+tit+".png", bbox_inches='tight')
      plt.close()

#----------------------------------------------------------------------------
#Print a 2D array georeferenced with contours
#----------------------------------------------------------------------------
def printcontours(ssp_lat,ssp_lon,lats,lons,cont_arr,pix,m2,tit,maxi=0,steps=20,save=0):
    try:
      pix=np.array(pix).astype(float)
      cont_arr=np.array(cont_arr).astype(float)
      pix[m2==0]=np.nan
      cont_arr[m2==0]=np.nan
      from mpl_toolkits.basemap import Basemap, cm
      import matplotlib.pyplot as plt
      #map = Basemap(projection='ortho',lat_0=ssp_lat,lon_0=ssp_lon,resolution='l')
      map = Basemap(projection='geos',lon_0=ssp_lon,resolution='l')
      map.drawcoastlines(linewidth=0.25)
      map.drawcountries(linewidth=0.25)
      x, y = map(lons, lats)
      if not maxi==0:
        clevs = np.linspace(np.nanmin(pix),maxi,steps)
      else:
        clevs = np.linspace(np.nanmin(pix),np.nanmax(pix),steps)
      clevs=np.round(clevs,decimals=3)
      #print clevs
      cs = map.contourf(x,y,pix,clevs,cmap=plt.cm.Reds)
      cbar = map.colorbar(cs,location='bottom',pad="5%")
      cbar.ax.set_xticklabels(clevs,rotation=90)
      cs2 = map.contour(x,y,cont_arr,colors='k',linewidths=1)
      plt.clabel(cs2, inline=1, fontsize=10)
      plt.title(tit)
      if save == 0:
        plt.show()
      else:
        plt.savefig("../test/"+tit+".png", bbox_inches='tight')
        plt.close()
    except ImportError:
      print "ERROR: basemap not installed!!"
      exit()
#----------------------------------------------------------------------------
#Print a 2D array georeferenced with contours
#----------------------------------------------------------------------------
def printdisk(ssp_lon,lats,lons,pix,tit,maxi=0,steps=20,save=0):
    try:
      pix=np.array(pix).astype(float)
      pix[lats<-999]=np.nan
      from mpl_toolkits.basemap import Basemap, cm
      import matplotlib.pyplot as plt
      #map = Basemap(projection='ortho',lat_0=ssp_lat,lon_0=ssp_lon,resolution='l')
      map = Basemap(projection='geos',lon_0=ssp_lon,resolution='l')
      map.drawcoastlines(linewidth=0.25)
      map.drawcountries(linewidth=0.25)
      x, y = map(lons, lats)
      if not maxi==0:
        clevs = np.linspace(np.nanmin(pix),maxi,steps)
      else:
        clevs = np.linspace(np.nanmin(pix),np.nanmax(pix),steps)
      clevs=np.round(clevs,decimals=3)
      #print clevs
      cs = map.contourf(x,y,pix,clevs,cmap=plt.cm.Reds)
      cbar = map.colorbar(cs,location='bottom',pad="5%")
      cbar.ax.set_xticklabels(clevs,rotation=90)
      plt.title(tit)
      if save == 0:
        plt.show()
      else:
        plt.savefig("../test/"+tit+".png", bbox_inches='tight')
        plt.close()
    except ImportError:
      print "ERROR: basemap not installed!!"
      exit()
#----------------------------------------------------------------------------
#Print a 2D array georeferenced with coastlines
#----------------------------------------------------------------------------
def printrecti(pix,LON,LAT,subX,subY,tit,save=0):
  from mpl_toolkits.basemap import Basemap, cm
  LAT=LAT[subY[0]:subY[1],subX[0]:subX[1]]
  LON=LON[subY[0]:subY[1],subX[0]:subX[1]]
  pix=pix[subY[0]:subY[1],subX[0]:subX[1]]
  fig = plt.figure(figsize=(18,18))
  ax = fig.add_axes([0.05,0.05,0.9,0.9])
  #m = Basemap(projection='ortho',lat_0=0,lon_0=0,resolution='h')
  m = Basemap(llcrnrlon=np.amin(LON[LON>-99]),llcrnrlat=np.amin(LAT[LAT>-99]),urcrnrlon=np.amax(LON),urcrnrlat=np.amax(LAT),projection='mill',resolution='h')
  im1 = m.pcolormesh(LON,LAT,pix,shading='flat',cmap=plt.cm.gray,latlon=True)
  m.drawcoastlines(color='r')
  m.drawparallels(np.arange(-90.,99.,30.),color='w')
  m.drawmeridians(np.arange(-180.,180.,60.),color='w')
  ax.set_title(tit)
  cb = m.colorbar(im1,"bottom", size="5%", pad="2%")
  if save==0:
    plt.show()
  else:
    plt.savefig("../test/"+tit+".png", bbox_inches='tight')
    plt.close()
#----------------------------------------------------------------------------
#Julian Date conversion
#----------------------------------------------------------------------------

def julian_day2(now):
    """
    1. Get current values for year, month, and day
    2. Same for time and make it a day fraction
    3. Calculate the julian day number via   https://en.wikipedia.org/wiki/Julian_day
    4. Add the day fraction to the julian day number

    """
    year = now.year
    month = now.month
    day = now.day
    day_fraction = now.hour + now.minute / 60.0 + now.second / 3600.0 / 24.0

    # The value 'march_on' will be 1 for January and February, and 0 for other months.
    march_on = math.floor((14 - month) / 12)
    year = year + 4800 - march_on
    # And 'month' will be 0 for March and 11 for February. 0 - 11 months
    month = month + 12 * march_on - 3

    y_quarter = math.floor(year / 4)
    jdn = day + math.floor((month * 153 + 2) / 5) + 365 * year + y_quarter

    julian = year < 1582 or year == (1582 and month < 10) or (month == 10 and day < 15)
    if julian:
        reform = 32083 # might need adjusting so needs a test
    else:
        reform = math.floor(year / 100) + math.floor(year / 400) + 32030.1875 # fudged this

    return jdn - reform + day_fraction
  
def julian_day(dt):
  return pyasl.jdcnv(dt)
  #return pyasl.juldate(dt)
#----------------------------------------------------------------------------
#DOY to Date conversion
#----------------------------------------------------------------------------
def doy2date(y,doy):
  dt=datetime.datetime(y, 1, 1) + datetime.timedelta(doy - 1)
  d=int(dt.strftime("%d"))
  m=int(dt.strftime("%m"))
  return (d,m)

#----------------------------------------------------------------------------
#Date to DOY conversion
#----------------------------------------------------------------------------
def date2doy(y,m,d):
  dt=datetime.datetime(y, m, d)
  doy=int(dt.strftime("%j"))
  return (doy)

#----------------------------------------------------------------------------
#Time to seconds conversion
#----------------------------------------------------------------------------
def time_to_seconds_since(year,m,d,t1):
  #FIXME: eventually the first slot is problematic because doy is already 
  #new day, while line is acquired in old day
  t   = datetime.datetime(year,m, d, 0, 0)
  t1s = np.multiply(t1,60*60)#convert t1 to seconds since midnight
  td=(t-datetime.datetime(1970,1,1))
  td=(td.microseconds + (td.seconds + td.days * 24 * 3600) * 10**6) / 10**6
  return (t1s,td)

#----------------------------------------------------------------------------
#Satname conversions
#----------------------------------------------------------------------------
def satnames(sat):
  satdict={"M2":"MET2","F2":"MET2","P2":"MET3","M3":"MET3","M4":"MET4","M5":"MET5","M6":"MET6","M7":"MET7"}
  return satdict[sat]

def satnamesR(sat):
  satdict={"MET2":"M2","MET3":"P2","MET4":"M4","MET5":"M5","MET6":"M6","MET7":"M7"}
  return satdict[sat]

def satnamesIRMon(sat):
  satdict={"MET2":"mviri_m2","MET3":"mviri_m3","MET4":"mviri_m4","MET5":"mviri_m5","MET6":"mviri_m6","MET7":"mviri_m7"}
  return satdict[sat]

def satnamesIRRef(sat):
  satdict={"MET2":"*_*","MET3":"*_*","MET4":"*_*","MET5":"*_*","MET6":"*_*","MET7":"*_*"}
  return satdict[sat]
#----------------------------------------------------------------------------
#Sat specific BT conversion
#----------------------------------------------------------------------------
def BTstuff(sat):
  bt_A_IR={"MET2":9.0107671,"MET3":8.9999579,"MET4":9.0228739,  "MET5":9.0211808, "MET6":9.0112547, "MET7":8.9819224}
  bt_A_WV={"MET2":10.528172,"MET3":10.517181,"MET4":10.654340,  "MET5":10.659564, "MET6":10.664768, "MET7":10.596521}
  bt_B_IR={"MET2":-1264.2322,"MET3":-1259.9131,"MET4":-1269.4460,"MET5":-1269.3131,"MET6":-1264.4411,"MET7":-1252.6850}
  bt_B_WV={"MET2":-2177.3002,"MET3":-2164.6560,"MET4":-2253.0677,"MET5":-2263.6804,"MET6":-2261.2885,"MET7":-2221.2987}
  return (bt_A_IR[sat],bt_A_WV[sat],bt_B_IR[sat],bt_B_WV[sat])
#----------------------------------------------------------------------------
#Get IR-Calfile filename
#----------------------------------------------------------------------------
def findcal(s_Sat,basepath):
  #calfile = MET7_calibration_coefficients.nc
  #basepath = /DSNNAS/Repro/aux/ricals/calibration_coefficients/
  wildcard = basepath+s_Sat+'_calibration_coefficients.nc'
  calfile  = glob.glob(wildcard)
  try:
    if len(calfile)<=1:
      return calfile[0]
    else:
      print "WARNING: too many Calibration Files!!!"
      print "...choosing arbitrarily one"
      return calfile[0]
  except IndexError:
      print "searched for: "+wildcard
      print "\nERROR: No Calibration File!!!\n" 
      return -1

#class makesrf():
  #def __init__(self,runfolder,sat,year,DOY,proc_vers,timestring,logs,debug=0):
    #self.home      = "/DSNNAS/Repro/mviri/level1/MFG15_FCDR_V1/tmp/"
    #self.sat       = sat
    #self.year      = year
    #self.DOY       = DOY
    #self.proc_vers = proc_vers
    #self.timestring= timestring
    #self.prelaunch = "../sensor_data/pre_launch/"
    #self.expected  = runfolder+"/aux/sensor_data/"+sat+timestring+"_srf.dat"
    #self.expected_ir=self.prelaunch+"spectral_response_"+self.sat[:3]+"_"+self.sat[3:]+"_IR.dat"
    #self.expected_wv=self.prelaunch+"spectral_response_"+self.sat[:3]+"_"+self.sat[3:]+"_WV.dat"
    #self.debug=debug
    #self.logs=logs
    #if not os.path.exists(runfolder+"/aux/sensor_data/"):
              #os.makedirs(runfolder+"/aux/sensor_data/")
    #self.features={#"MET7":["S10EE","0954","02",34],
                    ##"MET7":["S10EE","0956","02",54],
                    ##"MET7":["S10EE","0956","03",30],
                    ##"MET7":["S10EE","0956","10",54],
                    #"MET7":["S10EE","1000","10",55],
                    #"MET6":["S10GL","1000","10",55],
                    ##"MET5":["S10EL","0954","02",34],
                    #"MET5":["S10EL","1000","10",55],
                    ##"MET2":["S10GL","0955","01",34]
                    #"MET4":["S10EL","1000","10",55],
                    #"MET3":["S10EE","1801","10",55],
                    #"MET2":["S10EE","1000","11",56]
                    #}[self.sat]
    #self.bestsrf   =self.bestsrfname()
    #if self.debug!=0:
      #print (self.bestsrf)
  #def bestsrfname(self):
    #"""
    #function to reconstruct the filename of the closest matching reconstructed srf.
    #names should be like:
    #srf_MET7_1998168_1998169_0946_S10EE_01.dat
    #"""
    ##proc_vers  = "0954"#"0953"#"0951"#"0946"
    #if self.proc_vers in ["pre-launch","pre_launch","pre-flight","pre_flight"]:
      #return self.prelaunch+"spectral_response_"+self.sat[:3]+"_"+self.sat[3:]+"_VIS.dat"
    #else:
      #folder = self.home+"/recovered_srf_v"+self.proc_vers+"/"
      #model = self.features[0]
      #pover =self.features[1]
      #jobid = self.features[2]
      #try:
        #wildcard1=folder + "srf_"+self.sat+"_"
        #wildcard2=str(self.year) +"*"+"_"+pover+"*_"+model+"_"+jobid+".dat"
        #wildcard=wildcard1+wildcard2
        ##if self.debug==1:
        #print wildcard
        #lst=glob.glob(wildcard)
        #if self.debug==1:
          #print lst
        #diff=[]
        #for f in lst:
          #diff.append(int(f[len(wildcard1):len(wildcard1)+7])-int(str(self.year)+str(self.DOY)))
        #srfnamesource=lst[np.argmin(np.abs(diff))]
        ##print srfnamesource
      #except NameError:
        #print "nothing found for "+ wildcard
        #srfnamesource=-1
      #except ValueError:
        #print "ERROR: no recovered SRF"
        #srfnamesource=-1
      #return srfnamesource
  #def srfonfly(self,binN=0):
    #"""
    #function to tweak srf on the fly in terms of bins and uncertainty representation
    #"""
    #srfnamesource=self.bestsrf
    #srfnametarget=self.expected
    #if str(srfnamesource)=="-1":
      #return
    #if "pre_launch" in srfnamesource:
      #print "WARNING: using pre-launch SRF (did you want that????)"
      #headers=2
      #cov="n"
      #err=self.logs.update(satellite,curr_datetime,"problem_easy","pre-launch SRF")
    #else:
      #headers=self.features[3]
      #cov="y"
    ##read
    #allrows=np.loadtxt(srfnamesource,skiprows=headers)
    #if np.shape(allrows)[1]<=3:#apparently no cov was stored
      #print "WARNING: no covariance information"
      #err=self.logs.update(satellite,curr_datetime,"problem_easy","no cov in SRF")
      #covart=np.power(np.diag(allrows[:,2]),2)
      #allrows=np.hstack((allrows,covart))
      #if self.debug>=2:
        #print np.shape(allrows)
        #plt.imshow(covart)
        #plt.show()
    ##allrows[:,2]=np.diag(allrows[:,3:])#testing
    #if cov=="n":
      #df=pd.DataFrame(allrows[:,0:3],columns=['wl','res','err'])
    #else:
      #colnr=np.arange(-2,np.shape(allrows)[1]-2)
      #colnr=colnr.astype(str)
      #colnr[0:3]=['wl','res','err']
      #df=pd.DataFrame(allrows,columns=colnr)
    #if binN!=0:
      ##interpolate
      #xth=int((len(df)/binN)+0.5)
      #binN=len(df)/xth
      #if binN!=len(df)/xth:
        #print "WARNING: updated SRF resolution to be multiple of SRF N: "+str(binN)
      #if cov=="n":
        #allrows2 = np.array(df)[::xth,:][:-1,:]
      #else:
        ##according to FastOpt, reduction of the covariance via taking
        ##every Xth value is more correct than averaging the bins
        #allrowsA = np.array(df)[::xth,:3]
        #allrowsB = np.array(df)[::xth,3::xth]
        #allrows2 = np.column_stack((allrowsA,allrowsB))
    #else:
      #binN = len(df.wl)
      #allrows2 = np.array(df)
    #if self.debug>=2:
      #try:
        #print "orig: "+srfnamesource
        #print "plotting...."
        #print "integral of original SSR: " + str( np.trapz(allrows[:,1],x=allrows[:,0]))
        #print "integral of sub-sampled SSR: " + str( np.trapz(allrows2[:,1],x=allrows2[:,0]))
        #plt.plot(allrows[:,0],allrows[:,1],"-",label="original")
        #plt.plot(allrows2[:,0],allrows2[:,1],"-",label="resized")
        #plt.legend()
        #plt.show()
        #plt.plot(allrows[:,0],allrows[:,2],"-",label="FOpt uncert")
        #plt.plot(allrows[:,0],np.diag(allrows[:,3:]),"-",label="Diagonal")
        #plt.legend()
        #plt.show()
        #plt.imshow(allrows[:,3:])
        #plt.title("original")
        #plt.colorbar()
        #plt.show()
        #plt.imshow(allrows2[:,3:])
        #plt.title("sub sample")
        #plt.colorbar()
        #plt.show()
        #print np.shape(allrows2)
      #except:
        #print "WARNING: plotting failed; matplotlib installed?"
    ##write
    #print "creating "+srfnametarget
    #mini=np.amin(allrows2[0,:])
    #maxi=np.amax(allrows2[0,:])
    #step=((maxi-mini)/(binN-1))
    #head="  DataRelease = 1.0\n"+"  "+str(binN)+"    "+str(step)
    #np.savetxt(srfnametarget,allrows2, fmt= '%5.5f',delimiter='   ',header=head,newline="\n",comments='')	
    #return

class makesrf():
  def __init__(self,runfolder,sat,year,DOY,proc_vers,timestring="YYYYMMDD",logs=None,mode="nominal",member=0,cov="y",debug=0):
    #self.home      = "/DSNNAS/Repro/mviri/level1/MFG15_FCDR_V1/tmp/"
    self.home      = "/DSNNAS/Repro/mviri/level1/MFG15_FCDR_V1/tmp/FCDR_MVIRISRF.git/trunk/srf/"
    self.cov       = cov
    self.sat       = sat
    self.year      = year
    self.DOY       = DOY
    self.proc_vers = proc_vers
    self.timestring= timestring
    self.mode=mode
    if "nominal"in self.mode:
      self.prelaunch = "../sensor_data/pre_launch/"
      self.expected  = runfolder+"/aux/sensor_data/"+sat+timestring+"_srf.dat"
      if not os.path.exists(runfolder+"/aux/sensor_data/"):
          os.makedirs(runfolder+"/aux/sensor_data/")
    elif "sscc" in self.mode :
      self.prelaunch = "../sensor_data/"
      self.expected  = runfolder+"spectral_response_"+sat[:3]+"_"+sat[3:]+"_VIS.dat"
    self.expected_ir=self.prelaunch+"spectral_response_"+self.sat[:3]+"_"+self.sat[3:]+"_IR.dat"
    self.expected_wv=self.prelaunch+"spectral_response_"+self.sat[:3]+"_"+self.sat[3:]+"_WV.dat"
    self.debug=debug
    self.logs=logs
    self.features={ "MET7":["S10EE","1801","10",55],
                    "MET6":["S10EL","1801","10",55],
                    "MET5":["S10EL","1801","10",55],
                    "MET4":["S10EL","1801","10",55],
                    "MET3":["S10EE","1801","10",55],
                    "MET2":["S10EL","1801","10",55]
                    }[self.sat]
    self.bestsrf=self.bestsrfname()
    if self.bestsrf==-1:
      self.year=self.year-1
      self.DOY=365
      self.bestsrf=self.bestsrfname()
    if self.bestsrf==-1:
      self.year=self.year+2
      self.DOY=001
      self.bestsrf=self.bestsrfname()
  def bestsrfname(self):
    """
    function to reconstruct the filename of the closest matching reconstructed srf.
    names should be like:
    srf_MET7_1998168_1998169_0946_S10EE_01.dat
    """
    #proc_vers  = "0954"#"0953"#"0951"#"0946"
    if self.proc_vers in ["pre-launch","pre_launch","pre-flight","pre_flight"]:
      return self.prelaunch+"spectral_response_"+self.sat[:3]+"_"+self.sat[3:]+"_VIS.dat"
    else:
      model  = self.features[0]
      jobid  = self.features[2]
      wildc=self.home+"srf_"+self.sat+"_*_*_"+self.proc_vers+"-*_"+model+"_"+jobid+"/"
      print wildc
      folder = glob.glob(wildc)[0]
      try:
        wildcard1=folder + "srf_"+self.sat+"_"
        wildcard2=str(self.year) +"*"+"_"+self.proc_vers+"*_"+model+"_"+jobid+".dat"
        wildcard=wildcard1+wildcard2
        print wildcard
        if self.debug==1:
          print wildcard
        lst=glob.glob(wildcard)
        if self.debug==1:
          print lst
        diff=[]
        for f in lst:
          diff.append(int(f[len(wildcard1):len(wildcard1)+7])-int(str(self.year)+str(self.DOY)))
        srfnamesource=lst[np.argmin(np.abs(diff))]
        print srfnamesource
      except NameError:
        print "nothing found for "+ wildcard
        srfnamesource=-1
      except ValueError:
        print "ERROR: no recovered SRF"
        srfnamesource=-1
    return srfnamesource
  def srfonfly(self,binN=0):
    """
    function to tweak srf on the fly in terms of bins and uncertainty representation
    """
    srfnamesource=self.bestsrf
    srfnametarget=self.expected
    if str(srfnamesource)=="-1":
      return
    if "pre_launch" in srfnamesource:
      print "WARNING: using pre-launch SRF (did you want that????)"
      headers=2
      self.cov="n"
      if self.logs:
        err=self.logs.update(satellite,curr_datetime,"problem_easy","pre-launch SRF")
    else:
      headers=self.features[3]
      self.cov="y"
    #read
    allrows=np.loadtxt(srfnamesource,skiprows=headers)
    if np.shape(allrows)[1]<=3:#apparently no cov was stored
      print "WARNING: no covariance information"
      if self.logs:
        err=self.logs.update(satellite,curr_datetime,"problem_easy","no cov in SRF")
      covart=np.power(np.diag(allrows[:,2]),2)
      allrows=np.hstack((allrows,covart))
      if self.debug>=2:
        print np.shape(allrows)
        plt.imshow(covart)
        plt.show()
    #allrows[:,2]=np.diag(allrows[:,3:])#testing
    if self.cov=="n":
      df=pd.DataFrame(allrows[:,0:3],columns=['wl','res','err'])
    else:
      colnr=np.arange(-2,np.shape(allrows)[1]-2)
      colnr=colnr.astype(str)
      colnr[0:3]=['wl','res','err']
      df=pd.DataFrame(allrows,columns=colnr)
    if binN!=0:
      #interpolate
      xth=int((len(df)/binN)+0.5)
      binN=len(df)/xth
      if binN!=len(df)/xth:
        print "WARNING: updated SRF resolution to be multiple of SRF N: "+str(binN)
      if self.cov=="n":
        allrows2 = np.array(df)[::xth,:][:-1,:]
      else:
        #according to FastOpt, reduction of the covariance via taking
        #every Xth value is more correct than averaging the bins
        allrowsA = np.array(df)[::xth,:3]
        allrowsB = np.array(df)[::xth,3::xth]
        allrows2 = np.column_stack((allrowsA,allrowsB))
    else:
      binN = len(df.wl)
      allrows2 = np.array(df)
    if self.debug>=2:
      try:
        print "orig: "+srfnamesource
        print "plotting...."
        print "integral of original SSR: " + str( np.trapz(allrows[:,1],x=allrows[:,0]))
        print "integral of sub-sampled SSR: " + str( np.trapz(allrows2[:,1],x=allrows2[:,0]))
        plt.plot(allrows[:,0],allrows[:,1],"-",label="original")
        plt.plot(allrows2[:,0],allrows2[:,1],"-",label="resized")
        plt.legend()
        plt.show()
        plt.plot(allrows[:,0],allrows[:,2],"-",label="FOpt uncert")
        plt.plot(allrows[:,0],np.diag(allrows[:,3:]),"-",label="Diagonal")
        plt.legend()
        plt.show()
        plt.imshow(allrows[:,3:])
        plt.title("original")
        plt.colorbar()
        plt.show()
        plt.imshow(allrows2[:,3:])
        plt.title("sub sample")
        plt.colorbar()
        plt.show()
        print np.shape(allrows2)
      except:
        print "WARNING: plotting failed; matplotlib installed?"
    #write
    print "creating "+srfnametarget
    binN=len(allrows2[:,0])
    mini=np.amin(allrows2[:,0])
    maxi=np.amax(allrows2[:,0])
    step=((maxi-mini)/(binN-1))
    head="  DataRelease = 1.0\n"+"  "+str(binN)+"    "+str(step)
    np.savetxt(srfnametarget,allrows2, fmt= '%5.5f',delimiter='   ',header=head,newline="\n",comments='')       
    return
  def store(self,member,freqcol,weightcol,ucol=None,binN=0):
    """
    function to generate member of srf ensemble on the fly
    """
    srfnametarget=self.expected
    print member
    
    print "SRF will be stored to: "
    print srfnametarget
    if not ucol:
      ucol=np.ones(len(weightcol))*0.001
    allrows=np.hstack((freqcol[:,None],weightcol[:,None],ucol[:,None]))
    if self.cov=="n":
      df=pd.DataFrame(allrows[:,0:3],columns=['wl','res','err'])
    else:
      colnr=np.arange(-2,np.shape(allrows)[1]-2)
      colnr=colnr.astype(str)
      colnr[0:3]=['wl','res','err']
      df=pd.DataFrame(allrows,columns=colnr)
    if binN!=0:
      #interpolate
      xth=int((len(df)/binN)+0.5)
      binN=len(df)/xth
      if binN!=len(df)/xth:
        print "WARNING: updated SRF resolution to be multiple of SRF N: "+str(binN)
      if self.cov=="n":
        allrows2 = np.array(df)[::xth,:][:-1,:].astype(float)
      else:
        #according to FastOpt, reduction of the covariance via taking
        #every Xth value is more correct than averaging the bins
        allrowsA = np.array(df)[::xth,:3]
        allrowsB = np.array(df)[::xth,3::xth]
        allrows2 = np.column_stack((allrowsA,allrowsB)).astype(float)
    else:
      binN = len(df.wl)
      allrows2 = np.array(df).astype(float)
    if self.debug>=2:
      try:
        print "orig member: "+str(member)+self.timestring
        print "plotting...."
        print "integral of original SSR: " + str( np.trapz(allrows[:,1],x=allrows[:,0]))
        print "integral of sub-sampled SSR: " + str( np.trapz(allrows2[:,1],x=allrows2[:,0]))
        plt.plot(allrows[:,0],allrows[:,1],"-",label="original")
        plt.plot(allrows2[:,0],allrows2[:,1],"-",label="resized")
        plt.legend()
        plt.show()
        print np.shape(allrows2)
      except AttributeError:
        print "WARNING: plotting failed; matplotlib installed?"
        raise
    #write
    print "creating "+srfnametarget
    binN=len(allrows2[:,0])
    mini=np.amin(allrows2[:,0])
    maxi=np.amax(allrows2[:,0])
    step=((maxi-mini)/(binN-1)).astype(float)
    head="  DataRelease = 1.0\n"+"  "+str(binN)+"    "+str(step)
    print step
    np.savetxt(srfnametarget,allrows2, fmt= '%5.5f',delimiter='   ',header=head,newline="\n",comments='')       
    return
  def randomvectors(self,picklepath,N):
    '''
    largely adopted from Elisa Pinat, Rayference
    '''
    error_matrix = (self.bestsrf)
    datContent = [i.strip().split() for i in open(error_matrix).readlines()[self.features[3]:]]
    srf_matrix = {}
    srf = []
    for i in datContent:
      h = i[0].split('   ')
      info = (float(h[0]), float(i[1]), float(i[2]))
      srf.append(info)
      srf_matrix[float(h[0])] = []
      list_errs = []
      for j in i[3:]:
          list_errs.append(float(j))
      srf_matrix[float(h[0])] = list_errs  
    wavelength_ssr = []
    nominal_ssr = []
    u_ssr = []
    for i in srf:
        wavelength_ssr.append(i[0])
        nominal_ssr.append(i[1])
        u_ssr.append(i[2])
    wavelength_ssr = np.array(wavelength_ssr)
    nominal_ssr = np.array(nominal_ssr)
    u_ssr = np.array(u_ssr)
    Nsrf=len(wavelength_ssr)
    
    pickle.dump(wavelength_ssr, open(picklepath+'/'+self.sat+'_wavelength_ssr.p', 'wb'))

    # EIGENVALUES AND EIGENVECTOR DECOMPOSITION
    try:
      od = collections.OrderedDict(sorted(srf_matrix.items()))
    except AttributeError:
      od = OrderedDict(sorted(srf_matrix.items()))
    
    full_array = []
    for i, j in od.items():
        full_array.extend(j)

    final_matrix = np.array(full_array).reshape(Nsrf, Nsrf)
    
    pickle.dump(final_matrix, open(picklepath+'/'+self.sat+'_final_matrix.p', 'wb'))
    
    # # D = eigenvalues
    # # B = eigenvectors

    D,B = linalg.eig(final_matrix)
    #flip sign
    for z in range(5):
      if np.sum(B[:,z])<=0:
        B[:,z]=np.multiply(B[:,z],-1)
    
    pickle.dump(D, open(picklepath+'/'+self.sat+'_eigen_values.p', 'wb'))
    pickle.dump(B, open(picklepath+'/'+self.sat+'_eigen_vectors.p', 'wb'))

    # # Generate random vectors
    random_vectors = {}
    for i in range(N):
        r_v = np.empty(shape=(len(wavelength_ssr),1))
        r_v = np.random.normal(0, 1, r_v.shape)
        random_vectors[i] = r_v
    pickle.dump(random_vectors, open(picklepath+'/'+self.sat+'_random_vectors.p', 'wb'))
  def montemake(self,member,picklepath):
    '''
    largely adopted from Elisa Pinat, Rayference
    '''
    error_matrix = (self.bestsrf)
    print error_matrix
    datContent = [i.strip().split() for i in open(error_matrix).readlines()[self.features[3]:]]
    srf_matrix = {}
    srf = []
    for i in datContent:
      h = i[0].split('   ')
      info = (float(h[0]), float(i[1]), float(i[2]))
      srf.append(info)
      srf_matrix[float(h[0])] = []
      list_errs = []
      for j in i[3:]:
          list_errs.append(float(j))
      srf_matrix[float(h[0])] = list_errs  
    wavelength_ssr = []
    nominal_ssr = []
    u_ssr = []
    for i in srf:
        wavelength_ssr.append(i[0])
        nominal_ssr.append(i[1])
        u_ssr.append(i[2])
    wavelength_ssr = np.array(wavelength_ssr)
    nominal_ssr = np.array(nominal_ssr)
    u_ssr = np.array(u_ssr)
    Nsrf=len(wavelength_ssr)
    # EIGENVALUES AND EIGENVECTOR DECOMPOSITION
    try:
      od = collections.OrderedDict(sorted(srf_matrix.items()))
    except AttributeError:
      od = OrderedDict(sorted(srf_matrix.items()))
    
    full_array = []
    for i, j in od.items():
        full_array.extend(j)

    final_matrix = np.array(full_array).reshape(Nsrf, Nsrf)

    # # D = eigenvalues
    # # B = eigenvectors
    D=pickle.load(open(picklepath+'/'+self.sat+'_eigen_values.p', 'rb'))
    B=pickle.load(open(picklepath+'/'+self.sat+'_eigen_vectors.p', 'rb'))

    #D,B = linalg.eig(final_matrix)
    ##flip sign
    #for z in range(5):
    #  if np.sum(B[:,z])<=0:
    #    B[:,z]=np.multiply(B[:,z],-1)
        
    A = np.dot(B,np.dot(np.diag(D), linalg.inv(B)))

    #fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(16,7))
    #axes[0].grid(linestyle='dashed')
    ##axes[0].plot(D[0:5], 'o')
    #x = np.array(range(1,6))
    #y = np.real(D[0:5])
    #axes[0].scatter(x,y)
    ## axes[0]set(xticks=[1,2,3,4,5]), xticklabels=names)
    ## # axes[0].set_xticks(
    ## # names =[1,2,3,4,5])
    ##axes[0].semilogy(True)
    #axes[0].set_yscale('log')
    #axes[0].set_ylim(1e-3, 1)
    #plt.xticks(fontsize = 15) 
    #plt.yticks(fontsize = 15) 
    #axes[0].set_xlabel('Eigenvalues ranking', fontsize =22)
    #axes[0].set_ylabel('Eigenvalues', fontsize =22)
    ##plt.figure()
    #axes[1].plot(wavelength_ssr, B[:,0])
    #axes[1].plot(wavelength_ssr, B[:,1])
    #axes[1].plot(wavelength_ssr, B[:,2])
    #axes[1].plot(wavelength_ssr, B[:,3])
    #axes[1].plot(wavelength_ssr, B[:,4])
    #axes[1].set_xlabel('$\lambda$ ($\mu$m)', fontsize =22)
    #axes[1].set_ylabel('Eigenvector amplitude', fontsize =22)
    #plt.grid(linestyle='dashed')
    #plt.savefig('../tmp/eigenvalues/'+str(member)+'_'+self.timestring+'eigenvalues_eigenvectors.png', dpi=100)
    #plt.close()

    # # NOTATION of the DOCUMENT

    D_sqr = np.sqrt(np.diag(D))
    D_sqr.shape
    B.shape

    # # PERTURBATE THE SPECTRAL RESPONSE
    random_vectors=pickle.load(open(picklepath+'/'+self.sat+'_random_vectors.p', 'rb'))
    BD = np.dot(B,D_sqr)
    perturbed_ssr = {}
    temp_ssr = nominal_ssr.reshape((Nsrf, 1))
    perturbed_ssr = np.add(temp_ssr,np.dot(BD,random_vectors[member]))
    perturbed_ssr[np.real_if_close(perturbed_ssr)<0] = 0.0 # setting the spectral response function to zero when negative
    normalized_to_one = perturbed_ssr/perturbed_ssr.max()
    perturbed_ssr = normalized_to_one
    
    
    #fig, ax = plt.subplots(figsize=(15,9))
    #ax.plot(wavelength_ssr, np.diagonal(final_matrix), lw=2, alpha=0.5)
    #plt.xticks(fontsize = 15) 
    #plt.yticks(fontsize = 15) 
    #plt.xlabel('$\lambda$', fontsize =22)
    #plt.ylabel('diagonal of final matrix', fontsize =22)
    #plt.grid(linestyle='dashed')
    #plt.savefig('../tmp/final_matrix/'+str(member)+'_'+self.timestring+'final_matrix.png', dpi=100)
    #plt.close()
    
    #fig, ax = plt.subplots(figsize=(15,9))
    #ax.plot(wavelength_ssr, np.real(perturbed_ssr), lw=2, alpha=0.5)
    #ax.plot(wavelength_ssr, nominal_ssr,ls='dashed', lw=2)
    #ax.fill_between(wavelength_ssr, nominal_ssr+u_ssr, nominal_ssr-u_ssr,alpha=0.2)    
    #plt.xticks(fontsize = 15) 
    #plt.yticks(fontsize = 15) 
    #plt.xlabel('$\lambda$', fontsize =22)
    #plt.ylabel('Perturbed spectral response functions', fontsize =22)
    #plt.grid(linestyle='dashed')
    #plt.savefig('../tmp/perturbed_ssr/'+str(member)+'_'+self.timestring+'perturbed_ssr.png', dpi=100)
    #plt.close()
    #plt.show()
    
    self.store(member,wavelength_ssr,np.squeeze(perturbed_ssr),binN=200)


def calcSSP(head,head0,nomlon,early,logs,curr_datetime):
  imagtrace= "".join(getattr(head0[0],'szOriginalFormat'))
  sat="".join(getattr(head[0], 'szSatCode'))
  satellite=satnames(sat)
  orbx=getattr(head[1], 'dsOrbitCoordinatesFixedEarthImageStart')[0]
  orby=getattr(head[1], 'dsOrbitCoordinatesFixedEarthImageStart')[1]
  orbz=getattr(head[1], 'dsOrbitCoordinatesFixedEarthImageStart')[2]
  orbxE=getattr(head[1], 'dsOrbitCoordinatesFixedEarthImageEnd')[0]
  orbyE=getattr(head[1], 'dsOrbitCoordinatesFixedEarthImageEnd')[1]
  orbzE=getattr(head[1], 'dsOrbitCoordinatesFixedEarthImageEnd')[2]
  
  if not sat in early:
    #try simple calculation
    try:
      sslat = np.rad2deg(np.arctan(orbz/np.sqrt(orbx**2+orby**2)))
      sslon = nomlon+np.rad2deg(np.arctan(orby/orbx))
      sslatE = np.rad2deg(np.arctan(orbzE/np.sqrt(orbxE**2+orbyE**2)))
      sslonE = nomlon+np.rad2deg(np.arctan(orbyE/orbxE))
    except ZeroDivisionError:
      err=logs.append(satellite,curr_datetime,"problem_l15","ssp calc fails; orbit from suspicious recovery")
      sslat=900
      sslon=900
    if debug==1:
      print sslat,sslon,sslatE,sslonE
      print 'dsOrbitCoordinatesFixedEarthImageStart:',getattr(head[1],'dsOrbitCoordinatesFixedEarthImageStart')
      print 'dsOrbitCoordinatesFixedEarthImageEnd:',getattr(head[1],'dsOrbitCoordinatesFixedEarthImageEnd')
      print 'flEarthCentreInPixels:',getattr(head[1],'flEarthCentreInPixels')
    #check if totally off
    thr=8.
    if (abs(sslat-0.)>thr or abs(sslon-nomlon)>thr or abs(sslatE-0.)>thr or abs(sslonE-nomlon)>thr):
        err=logs.append(satellite,curr_datetime,"problem_l15","ssp suspicious; orbit from suspicious recovery")
        orbx=getattr(head[1], 'dsOrbitCoordinatesMgfImageStart')[0]
        orby=getattr(head[1], 'dsOrbitCoordinatesMgfImageStart')[1]
        orbz=getattr(head[1], 'dsOrbitCoordinatesMgfImageStart')[2]
        tjul=((getattr(head[1], 'iJulianSlotNum')[0]/48.)-(1./48.))+0.000694444
        orbxE=getattr(head[1], 'dsOrbitCoordinatesMgfImageEnd')[0]
        orbyE=getattr(head[1], 'dsOrbitCoordinatesMgfImageEnd')[1]
        orbzE=getattr(head[1], 'dsOrbitCoordinatesMgfImageEnd')[2]
        tjulE=((getattr(head[1], 'iJulianSlotNum')[0]/48.)-0.003472222)-0.000694444
        orbx,orby,orbz    = ob.mgf2eff(tjul, orbx, orby, orbz, nomlon)
        orbxE,orbyE,orbzE = ob.mgf2eff(tjulE,orbxE,orbyE,orbzE,nomlon)
        sslat = np.rad2deg(np.arctan(orbz/np.sqrt(orbx**2+orby**2)))
        sslon = nomlon+np.rad2deg(np.arctan(orby/orbx))
        sslatE = np.rad2deg(np.arctan(orbzE/np.sqrt(orbxE**2+orbyE**2)))
        sslonE = nomlon+np.rad2deg(np.arctan(orbyE/orbxE))
        if debug==1:
          print sslat,sslon,sslatE,sslonE
  else:
    #for early satellites convert coordinates first
    orbx=getattr(head[1], 'dsOrbitCoordinatesMgfImageStart')[0]
    orby=getattr(head[1], 'dsOrbitCoordinatesMgfImageStart')[1]
    orbz=getattr(head[1], 'dsOrbitCoordinatesMgfImageStart')[2]
    orbxE=getattr(head[1], 'dsOrbitCoordinatesMgfImageEnd')[0]
    orbyE=getattr(head[1], 'dsOrbitCoordinatesMgfImageEnd')[1]
    orbzE=getattr(head[1], 'dsOrbitCoordinatesMgfImageEnd')[2]
    if imagtrace == "IMAG2TG":#imag2tg ->"Start" tjul referring to southern horizon
      err=logs.append(satellite,curr_datetime,"problem_l15","ssp suspicious; orbit from suspicious recovery")
      print "WARNING: Coordinate recovery for MET2/MET3 not yet 100% accurate"
      tjul=((getattr(head[1], 'iJulianSlotNum')[0]/48.)-(1./48.))+0.000694444
      tjulE=((getattr(head[1], 'iJulianSlotNum')[0]/48.)-0.003472222)-0.000694444
    else:#imag1tg ->"Start" tjul still referring to first line 
      err=logs.append(satellite,curr_datetime,"problem_l15","ssp suspicious; orbit from good recovery")
      print "WARNING: Coordinate recovery for MET2/MET3 not yet 100% accurate"
      tjul=((getattr(head[1], 'iJulianSlotNum')[0]/48.)-(1./48.))
      tjulE=((getattr(head[1], 'iJulianSlotNum')[0]/48.)-0.003472222)
    orbx,orby,orbz    = ob.mgf2eff(tjul, orbx, orby, orbz, nomlon)
    orbxE,orbyE,orbzE = ob.mgf2eff(tjulE,orbxE,orbyE,orbzE,nomlon)
    try:
      sslat = np.rad2deg(np.arctan(orbz/np.sqrt(orbx**2+orby**2)))
      sslon = nomlon+np.rad2deg(np.arctan(orby/orbx))
      sslatE = np.rad2deg(np.arctan(orbzE/np.sqrt(orbxE**2+orbyE**2)))
      sslonE = nomlon+np.rad2deg(np.arctan(orbyE/orbxE))
    except ZeroDivisionError:
      err=logs.update(satellite,curr_datetime,"problem_l15","No MGF coordinates in header?")
      err=logs.update(satellite,curr_datetime,"status_easy",-1)
      err=logs.update(satellite,curr_datetime,"status_full",-1)
      err=logs.update(satellite,curr_datetime,"input_L1",0)
      raise
  return sslat, sslon, sslatE, sslonE

#----------------------------------------------------------------------------
#Integrate
#----------------------------------------------------------------------------
def int_2D(arr,xlin,ylin):
    """
    Integrate over an 2D array (e.g. covariance matrix)
    """
    Ny   = len(ylin)
    # do a 1-D integral over every row
    #-----------------------------------
    I = np.zeros( Ny )
    for i in range(Ny):
        I[i] = np.trapz( arr[i,:], xlin )
    # then an integral over the result
    #-----------------------------------    
    F = np.trapz( I, ylin )
    return F
  
#----------------------------------------------------------------------------
#Memory usage of the current process in kilobytes
#----------------------------------------------------------------------------
def memory_usage():
    """Memory usage of the current process in kilobytes."""
    status = None
    result = {'peak': 0, 'rss': 0}
    try:
        # This will only work on systems with a /proc file system
        # (like Linux).
        status = open('/proc/self/status')
        for line in status:
            parts = line.split()
            key = parts[0][2:-1].lower()
            if key in result:
                result[key] = int(parts[1])
    finally:
        if status is not None:
            status.close()
    return result

#----------------------------------------------------------------------------
"""
simple 3D interpolation of Z value at point I between point A and B
input:
  X - array of 3 X coordinates (X point A, X point I   , X point B)
  Y - array of 3 Y coordinates (Y point A, Y point I   , Y point B)
  Z - array of 2 Z coordinates (Z point A, NAN         , Z point B)
output:
  Z - array of 3 Z coordinates (Z point A, Interpolated, Z point B)
"""
#----------------------------------------------------------------------------
def interpol(X,Y,Z):
  print "interpol("+str(X)+","+str(Y)+","+str(Z)+")"
  V=Z[2]-Z[0]
  A=abs(Y[2]-Y[0])
  B=abs(X[2]-X[0])
  C=np.sqrt(A*A+B*B)
  VC=V/C
  P=abs(Y[1]-Y[0])
  Q=abs(X[1]-X[0])
  R=np.sqrt(P*P+Q*Q)
  Z[1]=R*VC+Z[0]
  return Z
#----------------------------------------------------------------------------
"""
simple 3D interpolation of Z value at point I using X and Y grids
input:
  X  - grid of X coordinates 
  Y  - grid of Y coordinates
  Z  - grid of Z coordinates
  XI - X coordinate at point I
  YI - Y coordinate at point I
output:
  ZI - Interpolated value at I
"""
#----------------------------------------------------------------------------
def interpol_cubic(X,Y,Z,XI,YI):
  from scipy.interpolate import griddata
  grid_x, grid_y = np.mgrid[np.amin(X):np.amax(X):100j, np.amin(Y):np.amax(Y):100j]
  points=np.column_stack((np.ndarray.flatten(X),np.ndarray.flatten(Y)))
  Z=griddata(points, np.ndarray.flatten(Z), (grid_x, grid_y), method='cubic')
  ilat = (np.abs(grid_y[0,:]-YI)).argmin()
  ilon = (np.abs(grid_x[:,0]-XI)).argmin()
  return Z[ilon,ilat]
#----------------------------------------------------------------------------
"""
simple 3D interpolation of Z value for a grid using X and Y grids
input:
  X  - grid of X coordinates 
  Y  - grid of Y coordinates
  Z  - grid of Z coordinates
output:
  ZI - Interpolated values
"""
#----------------------------------------------------------------------------
def interpol_grid(X,Y,Z):
  from scipy.interpolate import griddata
  grid_x, grid_y = np.mgrid[0:1:2500j, 0:1:2500j]
  points=np.column_stack((np.ndarray.flatten(X),np.ndarray.flatten(Y)))
  ZI=griddata(points, np.ndarray.flatten(Z), (grid_x, grid_y), method='linear')
  return ZI
#----------------------------------------------------------------------------
"""
manipulate a bit of an integer at location X to Y
  I  - Input integer
  X  - location of bit (MSB 0 --> from left to right)
  Y  - value for bit
  N  - nr of bits
output:
  MI - the modified integer value
"""
#----------------------------------------------------------------------------
def bitmod(I,X,Y=1,N=8):
  SI=list('%0*d' % (N, int(bin(I)[2:])))
  SI[X]=str(Y)
  MI=int("".join(SI),2)
  return MI

#----------------------------------------------------------------------------
"""
get range for a 2D fulldisk image for plotting
  I  - Input array
  C  - if true the MI/MA will be centered around zero
output:
  MI - the useful minimum value
  MA - the useful maximum value
"""
#----------------------------------------------------------------------------
def getrange(I,C=False):
  #infer X and Y dimensions
  X=np.shape(I)[0]
  Y=np.shape(I)[1]
  #discard outer D%
  D=10
  pX=(X/100)*D
  pY=(Y/100)*D
  I=I[pX:X-pX,pY:Y-pY]
  #compute percentiles
  MI=np.nanpercentile(I,10.0)
  MA=np.nanpercentile(I,90.0)
  if C:
    A=np.amax([abs(MI),abs(MA)])
    MI=A*-1
    MA=A
  return MI,MA
  
#----------------------------------------------------------------------------
from sqlalchemy import *
from sqlalchemy.ext.declarative import declarative_base
import psycopg2
class logger():
  def __init__(self,rel):
    self.dbserv='cdrlnxv01'#'localhost'
    self.table_name  = 'logtable'
    self.rel=rel
    self.dbname='ricalpy_{v1}'.format(v1=self.rel.replace('.','_'))
    self.user=os.environ.get('USER')
    err=self.connecting()
    Base = declarative_base()
    class ToDo(Base):
      __tablename__     = self.table_name
      idn               = Column(BigInteger,unique=True,primary_key=True)
      rect2lp_filename  = Column(String(255))
      platform          = Column(String(4))
      filename_easy     = Column(String(255))
      filename_full     = Column(String(255))
      scan_ts           = Column(DateTime)
      ricalpy_ts        = Column(DateTime)
      input_l10         = Column(Integer) #shall be: 1=available -1=notavailable
      input_l15         = Column(Integer) #shall be: 1=available  0=anomaly -1=notavailable
      status_easy       = Column(Integer) #shall be: 0=pending -1=failed 1=successful
      status_full       = Column(Integer) #shall be: 0=pending -1=failed 1=successful
      problem_easy      = Column(String(1000))
      problem_full      = Column(String(1000))
      problem_l15       = Column(String(1000))#comment about L15 file
  def connecting(self):
    conn = psycopg2.connect(dbname=self.dbname,user=self.user, host=self.dbserv, password='')
    conn.autocommit = True
    self.c=conn.cursor()
    return 0
  def createDB(self):
    try:
      self.c.execute('CREATE DATABASE {db};'.format(db=self.dbname))
    except psycopg2.ProgrammingError:
      print "DB already existing"
  def createTable(self):
    db_string='postgresql://{n}@{h}:5432/ricalpy_{v1}'.format(n=self.user,h=self.dbserv,v1=self.rel.replace('.','_'))
    engine = create_engine(db_string)
    self.Base.metadata.create_all(bind=engine)
    return 0
  def add(self,sat,current_datetime,
          L15_file="",easy_file="",full_file="",
          ricalpy_ts=datetime.datetime.utcnow(),
          L10_error=1,L15_error=1,easy_error=0,full_error=0,
          prob_easy="",prob_full="",prob_L15=""):
    idn=int(sat[3:]+current_datetime.strftime("%Y%m%d%H%M"))
    if L15_file=="":
      L15_file=str(idn)
    columnstring="({idx},{a},{b},{c},{d},{e},{f},{g},{h},{i},{j},{k},{l},{m})".format(
                    idx="idn",
                    a="rect2lp_filename",
                    b="platform",
                    c="filename_easy",
                    d="filename_full",
                    e="scan_ts",
                    f="ricalpy_ts",
                    g="input_l10",
                    h="input_l15",
                    i="status_easy",
                    j="status_full",
                    k="problem_easy",
                    l="problem_full",
                    m="problem_l15")
    valuestring="({idx},'{a}','{b}','{c}','{d}','{e}','{f}',{g},{h},{i},{j},'{k}','{l}','{m}')".format(
                    idx=idn,
                    a=L15_file,
                    b=sat,
                    c=easy_file,
                    d=full_file,
                    e=current_datetime,
                    f=ricalpy_ts,
                    g=L10_error,
                    h=L15_error,
                    i=easy_error,
                    j=full_error,
                    k=prob_easy,
                    l=prob_full,
                    m=prob_L15)
    try:
        sql_string="INSERT INTO {tn} {col} VALUES {val}".format(tn=self.table_name, col=columnstring, val=valuestring)
        self.c.execute(sql_string)
    except psycopg2.IntegrityError:
        #print 'ID already exists in PRIMARY KEY column {idx}'.format(idx=idn)
        return 1#ID already exists in PRIMARY KEY column {idx}'.format(idx=idn)
    return 0
  def query(self,sat,current_datetime,cols="*"):
    idn=int(sat[3:]+current_datetime.strftime("%Y%m%d%H%M"))
    sql_string="SELECT {col} FROM {tn} WHERE idn = {idx}".format(col=cols,tn=self.table_name, idx=idn)
    self.c.execute(sql_string)
    row = self.c.fetchone()
    colnames=self.c.description
    return colnames,row
  def update(self,sat,current_datetime,field,value):
    idn=int(sat[3:]+current_datetime.strftime("%Y%m%d%H%M"))
    if 'datetime.datetime' in str(type(value)) or 'str' in str(type(value)):
      sql_string="UPDATE {tn} SET {fi} = '{va}' WHERE idn = {idx}".format(tn=self.table_name,fi=field,va=value,idx=idn)
    else:
      sql_string="UPDATE {tn} SET {fi} = {va} WHERE idn = {idx}".format(tn=self.table_name,fi=field,va=value,idx=idn)
    self.c.execute(sql_string)
    return 0
  def append(self,sat,current_datetime,field,value):
    nam,cont=self.query(sat,current_datetime,field)
    try:
      cont=cont[0]+" "+value
    except TypeError:
      cont=value
    try:
      err=self.update(sat,current_datetime,field,cont)
    except psycopg2.DataError:
      cont="cleared; "+value
      err=self.update(sat,current_datetime,field,cont)
    return err
  def verify(self,sat,current_datetime):
    idn=int(sat[3:]+current_datetime.strftime("%Y%m%d%H%M"))
    sql_string="SELECT input_l10,input_l15,status_easy,status_full FROM {tn} WHERE idn = {idx}".format(tn=self.table_name, idx=idn)
    self.c.execute(sql_string)
    out=self.c.fetchone()
    if not out:
      self.add(sat,current_datetime,ricalpy_ts=datetime.datetime.utcnow())
      return False,out #false means no product for verification yet there!
    else:
      if out[0]==1 and out[1]==1 and (out[2]<1 or out[3]<1):
        return False,out #false means no product for verification yet there!
      elif out[0]==1 and out[1]==1 and out[2]>=1 and out[3]>=1:
        return True,out #True means product for verification there!
      #elif (out[0]!=1 or out[1]!=1) and (out[2]<1 or out[3]<1):
      elif (out[0]!=1 or out[1]==-1) and (out[2]<1 or out[3]<1):#this one does not exclude anomalies
        return True,out #true means no product there because no input
  def query2pandas(self):
    sql_string="SELECT * FROM logtable WHERE input_l10 = 1 AND input_l15 = 1"
    conn = psycopg2.connect(dbname=self.dbname,user=self.user, host=self.dbserv, password='')
    df=pd.read_sql(sql_string,conn)
    return df
  def query4injection(self):
    sql_string="SELECT * FROM logtable WHERE status_easy = 1"
    conn = psycopg2.connect(dbname=self.dbname,user=self.user, host=self.dbserv, password='')
    df=pd.read_sql(sql_string,conn)
    return df
  def check_easy(self,filename):
    '''
    checker for easy fcdr file
    '''
    error=1
    prob=""
    #first check filesize
    st = os.stat(filename)
    if (st.st_size/1000000)<30:
      error=-1
      prob="file too small"
    return (error,prob)
  def check_full(self,filename):
    '''
    checker for full fcdr file
    '''
    error=1
    prob=""
    #first check filesize
    st = os.stat(filename)
    if (st.st_size/1000000)<30:
                  error=-1
                  prob="file too small"
    return (error,prob)
  def check_rect2lp(self,filename,comp):
    '''
    checker for rect2lp file
    '''
    error=1
    prob=""
    #first check filesize
    st = os.stat(filename)
    if not (st.st_size/1000000)==38:
            error=-1
            prob="filesize wrong"
    #read file
    with open(filename, "rb") as fi:
            fi, head = rd.header(fi, 1, 3,t=0)
            nomlon = np.rad2deg(getattr(head[1],'dsNominalSubSatelliteLongitude')[0])
            lininfo, pixvis, pixir, pixwv = rd.image_3(fi)
            fi, trail = rd.header(fi, 2504, 2504, t=0)
    #get comparefile
    solartimedelta = nomlon/15 #sun moves 15 degrees lon per hour
    slotdelta      = int(solartimedelta*2+0.5)
    compareslot    = int(getattr(head[0],'iSlotNum')[0]) + slotdelta
    if compareslot>=48:
            compareslot=compareslot-48
    try:
            #H1=np.histogram(comp[compareslot,0,:,:], bins=256, range=[0, 256],density=True)
            #H1 = cv2.calcHist([cv2.equalizeHist(comp[compareslot,0,:,:].astype(np.uint8))],[0],None,[256],[0,256])
            H1 = cv2.calcHist([comp[compareslot,0,:,:].astype(np.uint8)],[0],None,[256/2],[0,256/2])
            H1=cv2.normalize(H1,None)
    except IndexError:
            error = -1
            prob="slot in metadata wrong: "+getattr(head[0],'iSlotNum')[0]
            return (error,prob)
    #H2=np.histogram(np.array(pixvis), bins=256, range=[0, 256],density=True)
    #H2 = cv2.calcHist([cv2.equalizeHist(np.array(pixvis).astype(np.uint8))],       [0],None,[256],[0,256])
    H2 = cv2.calcHist([np.array(pixvis).astype(np.uint8)],       [0],None,[256/2],[0,256/2])
    H2 = cv2.normalize(H2,None)
    #chi2=chi2_distance(H1[0],H2[0])
    #correlation test
    test=cv2.compareHist(H1, H2, cv2.HISTCMP_CORREL)
    print test
    if "METEOSAT2" in filename or "METEOSAT3" in filename:
            thr=0.4#for MET2 and 3 the difference might be quite large
    else:
            thr=0.6#relatively arbitrary threshold
    if test<thr:
            error = 0
            prob="histogram anomalous with R2="+str(test)
            return (error,prob)
    #chi2 test  
    test=cv2.compareHist(H1, H2, cv2.HISTCMP_CHISQR)
    print test
    if "METEOSAT2" in filename or "METEOSAT3" in filename:
            thr=300#for MET2 and 3 the difference might be quite large
    else:
            thr=200#relatively arbitrary threshold
    if test>thr:
            error = 0
            prob="histogram anomalous with Chi2="+str(test)
            return (error,prob)     
    #intersect test
    test=cv2.compareHist(H1, H2, cv2.HISTCMP_INTERSECT)
    print test
    if "METEOSAT2" in filename or "METEOSAT3" in filename:
            thr=0.25#for MET2 and 3 the difference might be quite large
    else:
            thr=0.5#relatively arbitrary threshold
    if test<thr:
            error = 0
            prob="histogram anomalous with intersect="+str(test)
            return (error,prob)
    test=cv2.compareHist(H1, H2, cv2.HISTCMP_BHATTACHARYYA)
    print test
    if "METEOSAT2" in filename or "METEOSAT3" in filename:
            thr=1#for MET2 and 3 the difference might be quite large
    else:
            thr=0.5#relatively arbitrary threshold
    if test>thr:
            error = 0
            prob="histogram anomalous with Bhattacharya="+str(test)
    return (error,prob)

