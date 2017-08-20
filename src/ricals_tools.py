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

try:
  sys.path.insert(0, os.getcwd()[:-4]+'/lib/orbit2eff/')
  import orbit2eff as ob
except ImportError:
  print "WARNING: ricals_tools without orbit correction functionality"

global lim
lim=[0,2500]
global debug
debug=0

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
  satdict={"M2":"MET2","P2":"MET3","M4":"MET4","M5":"MET5","M6":"MET6","M7":"MET7"}
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
  #FIXME: not complete for other sats then MET7
  unit_conv_IR={"MET2":0,"MET3":0,"MET4":0,"MET5":0,"MET6":0,"MET7":7.55978}
  unit_conv_WV={"MET2":0,"MET3":0,"MET4":0,"MET5":0,"MET6":0,"MET7":3.90293}
  bt_A_IR={"MET2":0,"MET3":0,"MET4":0,"MET5":0,"MET6":0,"MET7":6.9752}
  bt_A_WV={"MET2":0,"MET3":0,"MET4":0,"MET5":0,"MET6":0,"MET7":9.2477}
  bt_B_IR={"MET2":0,"MET3":0,"MET4":0,"MET5":0,"MET6":0,"MET7":-1255.4469}
  bt_B_WV={"MET2":0,"MET3":0,"MET4":0,"MET5":0,"MET6":0,"MET7":-2233.4882}
  return (unit_conv_IR[sat],unit_conv_WV[sat],bt_A_IR[sat],bt_A_WV[sat],bt_B_IR[sat],bt_B_WV[sat])
#----------------------------------------------------------------------------
#Get IR-Calfile filename
#----------------------------------------------------------------------------
def findcal(s_Sat,s_datestring,s_Time,basepath):
  #calfile = "/tcc1/proj/eraclim/ricals/mviri_m7_hirs_noaa14_2003_cal_values.nc"
  #basepath = "/tcc1/proj/eraclim/ricals/"
  satMon   = satnamesIRMon(s_Sat)
  satRef   = satnamesIRRef(s_Sat)
  wildcard = basepath+satMon+'_'+satRef+'_'+s_datestring[:4]+'_cal_values.nc'
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
      exit()
      #return -1
  print calfile
  exit()
#----------------------------------------------------------------------------
#Get SRF filename with closest date
#----------------------------------------------------------------------------
def bestsrfname(home,sat,year,DOY,proc_vers):
  """
  function to reconstruct the filename of the closest matching reconstructed srf.
  names should be like:
  srf_MET7_1998168_1998169_0946_S10EE_01.dat
  """
  
  #proc_vers  = "0954"#"0953"#"0951"#"0946"
  folder     = home+"/recovered_srf_v"+proc_vers+"/"
  model      = "S10EE"
  
  if proc_vers  == "0953":
    jobid      = "01"
  if proc_vers  == "0954":
    jobid      = "02"
    
  try:
    wildcard1=folder + "srf_"+sat+"_"
    wildcard2=str(year) +"*"+"_"+proc_vers+"*_"+model+"_"+jobid+".dat"
    wildcard=wildcard1+wildcard2
    if debug==1:
      print wildcard
    lst=glob.glob(wildcard)
    if debug==1:
      print lst
    diff=[]
    for f in lst:
      diff.append(int(f[len(wildcard1):len(wildcard1)+7])-int(str(year)+str(DOY)))
    #print sorted(diff)
    srfnamesource=lst[np.argmin(np.abs(diff))]
    #print srfnamesource
  except NameError:
    print "nothing found for "+ wildcard
    srfnamesource=-1
  except ValueError:
    srfnamesource=-1
  
  return srfnamesource

#----------------------------------------------------------------------------
#Modify srf on the fly
#----------------------------------------------------------------------------
#def srfonfly(srfnamesource,srfnametarget,binN=101):
  #"""
  #function to tweak srf on the fly in terms of bins and uncertainty representation

  #"""
  #cov="y"
  #headers=34
  #mini=0.3
  #maxi=1.3
  #step=((maxi-mini)/(binN-1))
  ##print step
  #maxi=maxi
  ##read
  #allrows=np.loadtxt(srfnamesource,skiprows=headers)
  #if cov=="n":
    #df=pd.DataFrame(allrows[:,0:3],columns=['wl','res','err'])
  #else:
    #colnr=np.arange(-2,np.shape(allrows)[1]-2)
    #colnr=colnr.astype(str)
    #colnr[0:3]=['wl','res','err']
    #print np.shape(colnr)
    #print np.shape(allrows)
    #df=pd.DataFrame(allrows,colnr)
  ##interpolate
  #if binN!=0:
    ##bins    = np.linspace(df.wl.min(), df.wl.max(), binN)
    #bint    = np.linspace(mini-(step/2), maxi+(step/2), binN+1) #for averaging it needs one step more
    #bins    = np.linspace(mini, maxi, binN)
    ##print bins
    #bins    = np.around(bins,2)
    #groups  = df.groupby(np.digitize(df.wl, bint))
    #df=pd.DataFrame(groups.mean())
    #allrows2 = np.array(df)[1:,:]
    #allrows2[:,0]=bins
  #else:
    #binN = len(df.wl)
    #allrows2 = np.array(df)
  #if debug==1:
    #try:
      #plt.plot(allrows[:,0],allrows[:,1],"-")
      #plt.plot(allrows2[:,0],allrows2[:,1],"-")
      #plt.show()
    #except:
      #print "WARNING: plotting failed; matplotlib installed?"
  ##write
  #head="  DataRelease = 1.0\n"+"  "+str(binN)+"    "+str(step)
  #np.savetxt(srfnametarget,allrows2, fmt= '%5.4f',delimiter='   ',header=head,newline="\n",comments='')
  ##with open(srfnametarget, 'w') as srftarget:
    ##srftarget.write("  DataRelease = 1.0\n")
    ##srftarget.write("  "+str(binN)+"    "+str((df.wl.max()-df.wl.min())/ binN)+"\n")
    ##for line in allrows:
      ##srftarget.write(str(line[0]))#note: there has to be an empty line below the last one!
def srfonfly(srfnamesource,srfnametarget,binN=101,debug=0):
  """
  function to tweak srf on the fly in terms of bins and uncertainty representation
  """
  cov="y"
  if "pre_launch" in srfnamesource:
    print "WARNING: using pre-launch SRF (did you want that????)"
    headers=2
  else:
    headers=34
  binN=0
  mini=0.3
  maxi=1.3
  #read
  allrows=np.loadtxt(srfnamesource,skiprows=headers)
  if np.shape(allrows)[1]<=3:#apparently no cov was stored
    print "WARNING: no covariance information"
    covart=np.power(np.diag(allrows[:,2]),2)
    allrows=np.hstack((allrows,covart))
    if debug>=2:
      print np.shape(allrows)
      plt.imshow(covart)
      plt.show()
  #allrows[:,2]=np.diag(allrows[:,3:])#testing
  if cov=="n":
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
    if cov=="n":
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

  if debug>=2:
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
  mini=np.amin(allrows2[0,:])
  maxi=np.amax(allrows2[0,:])
  step=((maxi-mini)/(binN-1))
  head="  DataRelease = 1.0\n"+"  "+str(binN)+"    "+str(step)
  np.savetxt(srfnametarget,allrows2, fmt= '%5.5f',delimiter='   ',header=head,newline="\n",comments='')
  #with open(srfnametarget, 'w') as srftarget:
    #srftarget.write("  DataRelease = 1.0\n")
    #srftarget.write("  "+str(binN)+"    "+str((df.wl.max()-df.wl.min())/ binN)+"\n")
    #for line in allrows:
      #srftarget.write(str(line[0]))#note: there has to be an empty line below the last one!



def calcSSP(head,head0,nomlon,early):
  imagtrace= "".join(getattr(head0[0],'szOriginalFormat'))
  sat="".join(getattr(head[0], 'szSatCode'))
  orbx=getattr(head[1], 'dsOrbitCoordinatesFixedEarthImageStart')[0]
  orby=getattr(head[1], 'dsOrbitCoordinatesFixedEarthImageStart')[1]
  orbz=getattr(head[1], 'dsOrbitCoordinatesFixedEarthImageStart')[2]
  orbxE=getattr(head[1], 'dsOrbitCoordinatesFixedEarthImageEnd')[0]
  orbyE=getattr(head[1], 'dsOrbitCoordinatesFixedEarthImageEnd')[1]
  orbzE=getattr(head[1], 'dsOrbitCoordinatesFixedEarthImageEnd')[2]
  
  if not sat in early:
    sslat = np.rad2deg(np.arctan(orbz/orbx))
    sslon = nomlon+np.rad2deg(np.arctan(orby/orbx))
    sslatE = np.rad2deg(np.arctan(orbzE/orbxE))
    sslonE = nomlon+np.rad2deg(np.arctan(orbyE/orbxE))
    if debug==1:
      print sslat,sslon,sslatE,sslonE
      print 'dsOrbitCoordinatesFixedEarthImageStart:',getattr(head[1],'dsOrbitCoordinatesFixedEarthImageStart')
      print 'dsOrbitCoordinatesFixedEarthImageEnd:',getattr(head[1],'dsOrbitCoordinatesFixedEarthImageEnd')
      print 'flEarthCentreInPixels:',getattr(head[1],'flEarthCentreInPixels')
    thr=8.
    if (abs(sslat-0.)>thr or abs(sslon-nomlon)>thr or abs(sslatE-0.)>thr or abs(sslonE-nomlon)>thr):
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
        sslat = np.rad2deg(np.arctan(orbz/orbx))
        sslon = nomlon+np.rad2deg(np.arctan(orby/orbx))
        sslatE = np.rad2deg(np.arctan(orbzE/orbxE))
        sslonE = nomlon+np.rad2deg(np.arctan(orbyE/orbxE))
        if debug==1:
          print sslat,sslon,sslatE,sslonE
  else:
    if imagtrace == "IMAG2TG":
      print "WARNING: Coordinate recovery for MET2/MET3 not yet implemented"
      sslat = np.rad2deg(np.arctan(orbz/orbx))
      sslon = nomlon+np.rad2deg(np.arctan(orby/orbz))
      sslatE = np.rad2deg(np.arctan(orbzE/orbxE))
      sslonE = nomlon+np.rad2deg(np.arctan(orbyE/orbzE))
    else:
      print "WARNING: Coordinate recovery for MET2/MET3 not yet implemented"
      sslat = np.rad2deg(np.arctan(orbz/orbx))
      sslon = nomlon+np.rad2deg(np.arctan(orby/orbz))
      sslatE = np.rad2deg(np.arctan(orbzE/orbxE))
      sslonE = nomlon+np.rad2deg(np.arctan(orbyE/orbzE))
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
