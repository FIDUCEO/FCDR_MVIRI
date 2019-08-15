#preparations:
#1. run sscc for 15 srf ensembles
#2. copy each ensemble to individual path:
#M="00"
#mkdir /tcenas/proj/eraclim/fiduceo_tmp/SSCC_Archive/results/MET7_MC$M
#cp -r /tcenas/proj/eraclim/fiduceo_tmp/SSCC_Archive/results/MET7_4.01/*_*_*_m$M /tcenas/proj/eraclim/fiduceo_tmp/SSCC_Archive/results/MET7_MC$M/
#or:
#sat=6
#rsync -av --ignore-existing /tcenas/proj/eraclim/fiduceo_tmp/SSCC_Archive/results/MET${sat}_4.01/*_*_*_m$M /tcenas/proj/eraclim/fiduceo_tmp/SSCC_Archive/results/MET${sat}_MC$M/
#or
#sat=2
#for M in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14; do rsync -av --ignore-existing /tcenas/proj/eraclim/fiduceo_tmp/SSCC_Archive/results/MET${sat}_4.01/*_*_*_m$M /tcenas/proj/eraclim/fiduceo_tmp/SSCC_Archive/results/MET${sat}_MC$M/; done
#DEL=/tcenas/proj/eraclim/fiduceo_tmp/SSCC_Archive/results/
#for M in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 14; do find ${DEL}MET${sat}_MC$M/ -empty -type d -delete; done
#DEL=/tcenas/proj/eraclim/fiduceo_tmp/SSCC_Archive/results/MET2_4.01/
#find ${DEL} -empty -type d -delete
#3. set environment
#source /tcenas/home/frankr/git/RICalPy/environment.sh
##cd /tcenas/home/frankr/git/sscc_calib_tool
#export CALIB_MC=/DSNNAS/Repro_Temp/proj/FIDUCEO/MonteCarloSSCC/
#export CALIB_STORE=/tcenas/home/frankr/git/sscc_calib_tool/
#rm ${CALIB_MC}/config/M${sat}*_RMC*
#python /DSNNAS/Repro_Temp/users/frankr/CDR/git/RICalPy/tools/cf_errors.py

def interpolate_ssi(srf,df1,bint,binN):
  if len(df1)<binN:
    bint=srf[:,0]
    ss=np.vstack((bint,np.interp(bint,df1.wl,df1.rad))).T
    df2=pd.DataFrame(ss[:,0:2],columns=['wl','rad'])
    ssi=np.array(df2.rad)
  else:
    bint=srf[:,0]
    step=((np.amax(bint)-np.amin(bint))/(len(bint)-1))
    bin_borders=np.hstack((np.subtract(bint,step/2),np.amax(bint)+(step/2)))
    bin_means = (np.histogram(df1.wl, bin_borders.T, weights=df1.rad)[0] /
                np.histogram(df1.wl, bin_borders.T)[0])
    mask = np.isnan(bin_means)
    bin_means[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), bin_means[~mask])
    #print np.shape(bint),np.shape(bin_means)
    ss=np.vstack((bint,bin_means)).T
    df2=pd.DataFrame(ss,columns=['wl','rad'])
    #groups  = df1.groupby(np.digitize(df1.wl, bint))
    #df2=pd.DataFrame(groups.mean())
    ssi=np.array(df2.rad)
  return ssi
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
def band_irrad(sundist,sdu,srffile="../sensor_data/pre_launch/spectral_response_MET_7_VIS.dat",ssfile="../config/solirr_yves.txt",d=0):
  """
  The default spectrum is the one described in EUM/MTG/DEF/10/0611
  """
  AU=149597870700     #[m]
  sunradius=6.96342*10**8 #[m]
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
  bint = np.linspace(mini-(step/2), maxi+(step/2), binN+1)
  bins = np.linspace(mini, maxi, binN)
  df1=pd.DataFrame(ssf[:,0:2],columns=['wl','rad'])
  df1=df1.loc[(df1.wl>=mini) & (df1.wl<=maxi)]
  ssi = interpolate_ssi(srf,df1,bint,binN)
  # 4. convolute both 
  #weighting
  f  = ssi*srf[:,1]#*np.median(np.diff(srf[:,0]))
  #df = ssi*srf[:,2]#*np.median(np.diff(srf[:,0]))
  #print np.shape(ssi),np.shape(srf[:,3:]),np.shape(ssi[:, None])
  try:
    df = ssi*srf[:,3:]*ssi[:, None]#includes covariances
  except ValueError:
    df = np.diag(srf[:,2])
  F  = np.trapz(f,x=srf[:,0])#np.sum(f)
  dF=np.sqrt(int_2D(df,srf[:,0],srf[:,0]))#includes covariances
  return F,dF,srf

import numpy as np
import pandas as pd
import calib
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['axes.linewidth'] = 2

datestring=raw_input("datestring (should have a sscc run in all ensembles and in nominal e.g.: 199806101200): ")
sat=raw_input("sat (e.g.: MET7): ")
release=raw_input("release (e.g.: 3.1): ")
outfile="/tcenas/home/frankr/git/RICalPy/config/effect_correlation_"+sat+"_"+str(release)+".csv"

print "load nominal as ref"
#load nominal as reference
y=int(datestring[:4])
m=int(datestring[4:6])
d=int(datestring[6:8])
nom=calib.calib(sat,release)
doy1=nom.date2doy(y,m,d)
doy2=doy1+4
nom.datestring=datestring
nom.prepcf()
#load nominal irradiance
ssf="/tcenas/home/frankr/git/RICalPy/config/solirr_yves.txt"
srff=nom.ssccfolder()+nom.sat+"_4.01/"+nom.sat+"_"+str(y)+"_"+str(doy1).zfill(3)+"_"+str(doy2).zfill(3)+"/spectral_response_"+nom.sat[:3]+"_"+nom.sat[3:]+"_VIS.dat"
nomI=band_irrad(1,0,srffile=srff,ssfile=ssf)

print "load members"
#load members
cf_obj=[]
a0_err=[]
a1_err=[]
a2_err=[]
cf_err=[]
ensI=[]
I_err=[]
N=15


for member in range(N):
  release="MC"+str(member).zfill(2)
  #calib coefs
  co=calib.calib(sat,release)  
  co.datestring=datestring
  co.prepcf()
  cf_obj.append(co)
  a0_err.append(co.a0-nom.a0)
  a1_err.append(co.a1-nom.a1)
  a2_err.append(co.a2-nom.a2)
  cf_err.append(co.cf-nom.cf)
  #Irradiance
  co=cf_obj[member]
  srff=co.ssccfolder()+co.sat+"_MC"+str(member).zfill(2)+"/"+co.sat+"_"+str(y)+"_"+str(doy1).zfill(3)+"_m"+str(member).zfill(2)+"/spectral_response_"+nom.sat[:3]+"_"+nom.sat[3:]+"_VIS.dat"
  I=band_irrad(1,0,srffile=srff,ssfile=ssf)
  ensI.append(I)
  I_err.append(I[0]-nomI[0])
EFs=[I_err,a0_err,a1_err,a2_err,cf_err]
EFnames=["E","a0","a1","a2","Z"]
#correlation
r=np.corrcoef(np.array(EFs))
print EFnames
print r

#save correlation matrix
Es=["E","a0","a1","a2","Z","SZA","Cs"]
q=np.ones(len(Es))
qd=np.diag(q)
qd[:len(EFnames),:len(EFnames)]=r

np.savetxt(outfile,qd,delimiter=";")

plt.imshow(qd)
plt.xticks(range(0,len(Es)), Es)
plt.yticks(range(0,len(Es)), Es)
cbar=plt.colorbar()
cbar.set_label('correlation', rotation=270)
plt.show()

for ci in cf_obj:
  plt.plot(ci.a0)
plt.show()

#plot results
fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(16,16))
axes[0,0].grid(linestyle='dashed')
x = I_err
y = a0_err
axes[0,0].scatter(x,y)
plt.xticks(fontsize = 10) 
plt.yticks(fontsize = 10) 
axes[0,0].set_xlabel('$\delta$ Irradiance [Wm^-2]', fontsize =22)
axes[0,0].set_ylabel('$\delta$ a0 [Wm^-2sr^-1/DC]', fontsize =22)
plt.grid(linestyle='dashed')
axes[0,0].text(-0.1, 1.0, "A", transform=axes[0,0].transAxes, size=20, weight='bold')

axes[0,1].grid(linestyle='dashed')
x = I_err
y = a1_err
axes[0,1].scatter(x,y)
plt.xticks(fontsize = 10) 
plt.yticks(fontsize = 10) 
axes[0,1].set_xlabel('$\delta$ Irradiance [Wm^-2]', fontsize =22)
axes[0,1].set_ylabel('$\delta$ a1 [Wm^-2sr^-1/DCyr^-1]', fontsize =22)
plt.grid(linestyle='dashed')
axes[0,1].text(-0.1, 1.0, "B", transform=axes[0,1].transAxes, size=20, weight='bold')

axes[0,2].grid(linestyle='dashed')
x = I_err
y = a2_err
axes[0,2].scatter(x,y)
plt.xticks(fontsize = 10) 
plt.yticks(fontsize = 10) 
axes[0,2].set_xlabel('$\delta$ Irradiance [Wm^-2]', fontsize =22)
axes[0,2].set_ylabel('$\delta$ a2 [Wm^-2sr^-1/DCyr^-1]', fontsize =22)
plt.grid(linestyle='dashed')
axes[0,2].text(-0.1, 1.0, "C", transform=axes[0,2].transAxes, size=20, weight='bold')

axes[1,0].grid(linestyle='dashed')
x = a1_err
y = a0_err
axes[1,0].scatter(x,y)
plt.xticks(fontsize = 10) 
plt.yticks(fontsize = 10) 
axes[1,0].set_ylabel('$\delta$ a0 [Wm^-2sr^-1/DC]', fontsize =22)
axes[1,0].set_xlabel('$\delta$ a1 [Wm^-2sr^-1/DCyr^-1]', fontsize =22)
plt.grid(linestyle='dashed')
axes[1,0].text(-0.1, 1.0, "D", transform=axes[1,0].transAxes, size=20, weight='bold')

axes[1,1].grid(linestyle='dashed')
x = a2_err
y = a0_err
axes[1,1].scatter(x,y)
plt.xticks(fontsize = 10) 
plt.yticks(fontsize = 10) 
axes[1,1].set_ylabel('$\delta$ a0 [Wm^-2sr^-1/DC]', fontsize =22)
axes[1,1].set_xlabel('$\delta$ a2 [Wm^-2sr^-1/DCyr^-1]', fontsize =22)
plt.grid(linestyle='dashed')
axes[1,1].text(-0.1, 1.0, "E", transform=axes[1,1].transAxes, size=20, weight='bold')

axes[1,2].grid(linestyle='dashed')
nomsrf=nomI[2]
axes[1,2].plot(nomsrf[:,0], nomsrf[:,1])
axes[1,2].fill_between(nomsrf[:,0],nomsrf[:,1]-nomsrf[:,2] , nomsrf[:,1]+nomsrf[:,2], alpha=0.3)
for I in ensI:
  srf=I[2]
  axes[1,2].plot(srf[:,0], srf[:,1])
plt.xticks(fontsize = 10) 
plt.yticks(fontsize = 10) 
axes[1,2].set_xlabel('$\lambda$ ($\mu$m)', fontsize =22)
axes[1,2].set_ylabel('SRFs at t=0', fontsize =22)
plt.grid(linestyle='dashed')
axes[1,2].text(-0.1, 1.0, "F", transform=axes[1,2].transAxes, size=20, weight='bold')

plt.show()

