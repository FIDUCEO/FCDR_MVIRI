"""
this file contains functions for effect calculations
"""
import sys 
import numpy as np
import netCDF4 as nc
import time
import timeit
import csv
import os
import glob
import allantools

import pickle

import ricals_tools as to

import cruncher as cr
import matplotlib.pyplot as plt

def law_of_error_trivial(corr,uncert,sensi):
  #generate covariance matrix from uncertainties and correlations
  cov_mat=np.copy(corr)
  imp_mat=np.copy(corr)
  for i in range(corr.shape[0]):
    for j in range(corr.shape[1]):
      cov_mat[i,j]=uncert[i]*corr[i,j]*uncert[j]
      imp_mat[i,j]=sensi[i]*cov_mat[i,j]*sensi[j]
  #integrate
  #comb=to.int_2D(imp_mat,range(np.shape(imp_mat)[0]),range(np.shape(imp_mat)[1]))
  #sum
  comb=imp_mat.sum(axis=(0,1))
  #root
  u_comb=np.sqrt(comb)
  return u_comb,cov_mat,imp_mat

def law_of_error_normal(corr,uncert,sensi):
  #generate covariance matrix from uncertainties and correlations
  cov_mat=uncert*corr*uncert[:, None]
  #multiply with the sensitivity coefficients
  imp_mat=sensi*cov_mat*sensi[:, None]
  #integrate
  #comb=to.int_2D(imp_mat,range(np.shape(imp_mat)[0]),range(np.shape(imp_mat)[1]))
  #sum
  comb=imp_mat.sum(axis=(0,1))
  #root
  u_comb=np.sqrt(comb)
  return u_comb,cov_mat,imp_mat

def rad_noise(head,trail):
    #calculate radiometric noise
    meanS=getattr(head[2], 'flSpaceCornersVISS')
    meanN=getattr(head[2], 'flSpaceCornersVISN')
    stdS =getattr(head[2], 'flStdSpaceCornersVISN')
    stdN =getattr(head[2], 'flStdSpaceCornersVISS')
    a=np.divide(np.sum([np.power(stdS[k],2)+np.power(stdS[k],2) for k in range(4)]),8)
    k0=np.divide(np.mean(meanS)+np.mean(meanN),2)
    u_rad=np.sqrt(a+np.power(np.divide(np.mean(meanS)-np.mean(meanN),2),2))

    return u_rad,k0

def comp_allan(vals):
    mn=np.mean(vals)
    sd=np.std(vals)
    ad=allantools.adev(np.array(vals), rate=1,taus=2)[1]
    return mn,sd,ad
  
def space_eval(values,corner_means,corner_stds,logs,satellite,curr_datetime,debug=0):
  #dump
  #pickle.dump(values, open('./tmp/values.p', 'wb'))
  #pickle.dump(corner_means, open('./tmp/corner_means.p', 'wb'))
  #pickle.dump(corner_stds, open('./tmp/corner_stds.p', 'wb'))
  #values=pickle.load( open('./tmp/values.p', 'rb'))
  #corner_means=pickle.load( open('./tmp/corner_means.p', 'rb'))
  #corner_stds=pickle.load( open('./tmp/corner_stds.p', 'rb'))
  #VIS1
  mn1,sd1,ad1 = comp_allan(values[0,:])
  #VIS2
  mn2,sd2,ad2 = comp_allan(values[1,:])
  usable=[True,True]
  if mn1<0:#if only detector 2 was on
    usable[0]=False
  if mn2<0:#if only detector 1 was on
    usable[1]=False
  #all
  mn,sd,ad    = comp_allan(np.ndarray.flatten(values[values>=0]))
  #check corners
  rigor=4
  uthr1=mn1+(rigor*sd1)
  uthr2=mn2+(rigor*sd2)
  lthr1=mn1-(rigor*sd1)
  lthr2=mn2-(rigor*sd2)
  #first sensor
  flags=[]
  if usable[0]:
    uthr=uthr1
    lthr=lthr1
    for corner in [0,1,2,3]:
      cm=corner_means[corner]
      if cm<(lthr) or cm>(uthr):
        flags.append(1)
      else:
        flags.append(0)
  else:
      flags=[1,1,1,1]
  if usable[1]:
    #second sensor
    uthr=uthr2
    lthr=lthr2
    for corner in [4,5,6,7]:
      cm=corner_means[corner]
      if cm<(lthr) or cm>(uthr):
        flags.append(1)
      else:
        flags.append(0)
  else:
    flags=flags+[1,1,1,1]
  #evaluate
  if 1 in flags[:4]:
    err=logs.append(satellite,curr_datetime,"problem_l15","space corner suspicious; flagged: "+str(flags))
    Cs_flag=1
    X=[0,1,2,3,4,5,6,7]
    if debug>0:
      plt.plot(X,[uthr1,uthr1,uthr1,uthr1,uthr2,uthr2,uthr2,uthr2],"--",color="red",label="2sd")
      plt.plot(X,[lthr1,lthr1,lthr1,lthr1,lthr2,lthr2,lthr2,lthr2],"--",color="red",label="2sd")
      plt.plot(X,[mn1,mn1,mn1,mn1,mn2,mn2,mn2,mn2],"--",color="green",label="detector mean")
      plt.errorbar(X,corner_means,yerr=corner_stds,fmt=".",label="corner means")
      plt.show()
    #redo VIS1
    if usable[0]:
      vals=values[0,:]
      mn1,sd1,ad1 = comp_allan(vals[(vals>int(lthr1-1.5))&(vals<int(uthr1+1.5))])
  elif 1 in flags[4:]:
    err=logs.append(satellite,curr_datetime,"problem_l15","space corner suspicious; flagged: "+str(flags))
    Cs_flag=1
    #redo VIS2
    if usable[1]:
      vals=values[1,:]
      mn2,sd2,ad2 = comp_allan(vals[(vals>int(lthr2-1.5))&(vals<int(uthr2+1.5))])
  else:
    Cs_flag=0
  #combine
  #u_all->measure of the noise level at earth counts
  if not False in usable:
    a=np.divide(np.power(ad1,2)+np.power(ad2,2),2)
    u_all=np.sqrt(a+np.power(np.divide(mn1-mn2,2),2))
    #k0->dark signal
    k0=np.divide(mn1+mn2,2)
  else:
    if usable[0]:
      u_all=ad1
      k0=mn1
    elif usable[1]:
      u_all=ad2
      k0=mn2
    elif not True in usable:
      err=logs.append(satellite,curr_datetime,"problem_easy","no valid spacecorner (masked?); ")
      err=logs.append(satellite,curr_datetime,"problem_full","no valid spacecorner (masked?); ")
      u_all=4
      k0=0
  #u_Cs -> appropriateness of space corner counts as dark signal for entire image
    #std between means of 4 corners of detector 1
  tmp1=np.std(corner_means[:4])
    #std between means of 4 corners of detector 2
  tmp2=np.std(corner_means[4:])
    #std between means of detector 1 and 2
  tmpD=np.std([np.mean(corner_means[:4]),np.mean(corner_means[4:])])
    #combine (variance of each detector with variance between detectors)
  if not False in usable:
    u_Cs=np.sqrt(tmp1**2+tmp1**2+tmpD**2)
  else:
    if usable[0]:
      u_Cs=tmp1
    elif usable[1]:
      u_Cs=tmp2
    elif not True in usable:
      u_Cs=2
  return u_all,k0,u_Cs,Cs_flag

def space_header_eval(sc_meanS_head,sc_meanN_head,sc_sdevS_head,sc_sdevN_head,logs,satellite,curr_datetime,debug=0):
  corner_means=np.hstack([sc_meanS_head,sc_meanN_head])
  corner_means=np.hstack([sc_meanS_head,sc_meanN_head])
  #VIS1
  mn1=np.mean(sc_meanS_head)
  sd1=np.sqrt(np.sum(np.power(sc_sdevS_head,2))/2)
  #VIS2
  mn2=np.mean(sc_meanN_head)
  sd2=np.sqrt(np.sum(np.power(sc_sdevN_head,2))/2)
  #all
  mn    = np.mean([mn1,mn2])
  sd    = np.sqrt(0.5*(sd1**2+sd2**2)+((mn1-mn2)/2)**2)
  usable=[True,True]
  #Assumption: if mean is 0. detector was off (even though normally it is then -1)
  if mn1<=0:#if only detector 2 was on
    usable[0]=False
    mn=mn2
    sd=sd2
  if mn2<=0:#if only detector 1 was on
    usable[1]=False
    mn=mn1
    sd=sd1
  #check corners
  rigor=4#how many standard deviations for thresholds
  uthr1=mn1+(rigor*sd1)
  uthr2=mn2+(rigor*sd2)
  lthr1=mn1-(rigor*sd1)
  lthr2=mn2-(rigor*sd2)
  #first sensor
  flags=[]
  if usable[0]:#if detector 1 was on check if all corners in range
    uthr=uthr1
    lthr=lthr1
    for corner in [0,1,2,3]:
      cm=corner_means[corner]
      if cm<(lthr) or cm>(uthr):
        flags.append(1)
      else:
        flags.append(0)
  else:        #if detector 1 was off set all flags to 1
    flags=[1,1,1,1]
  #second sensor
  if usable[1]:#if detector 1 was on check if all corners in range
    uthr=uthr2
    lthr=lthr2
    for corner in [4,5,6,7]:
      cm=corner_means[corner]
      if cm<(lthr) or cm>(uthr):
        flags.append(1)
      else:
        flags.append(0)
  else:#if detector 1 was off set all flags to 1
    flags=flags+[1,1,1,1]
  #evaluate
  if 1 in flags[:4]:#if any corner was flagged: something is strange -> mitigate
    err=logs.append(satellite,curr_datetime,"problem_l15","header space corner suspicious; flagged: "+str(flags))
    Cs_flag=1
    X=[0,1,2,3,4,5,6,7]
    if debug>0:
      plt.plot(X,[uthr1,uthr1,uthr1,uthr1,uthr2,uthr2,uthr2,uthr2],"--",color="red",label="2sd")
      plt.plot(X,[lthr1,lthr1,lthr1,lthr1,lthr2,lthr2,lthr2,lthr2],"--",color="red",label="2sd")
      plt.plot(X,[mn1,mn1,mn1,mn1,mn2,mn2,mn2,mn2],"--",color="green",label="detector mean")
      plt.errorbar(X,corner_means,yerr=corner_stds,fmt=".",label="corner means")
      plt.show()
    #redo VIS1 with corners that are not flagged
    if usable[0]:
      mn1=np.mean(sc_meanS_head[flags[:4]==0])
      sd1=np.sqrt(np.sum(np.power(sc_sdevS_head[flags[:4]==0],2))/2)
  elif 1 in flags[4:]:#if any corner was flagged: something is strange -> mitigate
    err=logs.append(satellite,curr_datetime,"problem_l15","header space corner suspicious; flagged: "+str(flags))
    Cs_flag=1
    #redo VIS2 with corners that are not flagged
    if usable[1]:
      mn2=np.mean(sc_meanN_head[flags[4:]==0])
      sd2=np.sqrt(np.sum(np.power(sc_sdevN_head[flags[:4]==0],2))/2)
  else:
    Cs_flag=0   #if no corner was flagged: evertything fine
  #combine
  #u_all->measure of the noise level at earth counts
  if not False in usable: #if both detectors were switched on
    a=np.divide(np.power(sd1,2)+np.power(sd2,2),2)
    u_all=np.sqrt(a+np.power(np.divide(mn1-mn2,2),2))
    #k0->dark signal
    k0=np.divide(mn1+mn2,2)
  else:#if only one detector was switched on
    if usable[0]:
      u_all=sd1
      k0=mn1
    elif usable[1]:
      u_all=sd2
      k0=mn2
    elif not True in usable:#if no detector was switched on
      err=logs.append(satellite,curr_datetime,"problem_easy","no valid header spacecorner; ")
      err=logs.append(satellite,curr_datetime,"problem_full","no valid header spacecorner; ")
      u_all=4
      k0=0
  #u_Cs -> appropriateness of space corner counts as dark signal for entire image
    #std between means of 4 corners of detector 1
  tmp1=np.std(corner_means[:4])
    #std between means of 4 corners of detector 2
  tmp2=np.std(corner_means[4:])
    #std between means of detector 1 and 2
  tmpD=np.std([np.mean(corner_means[:4]),np.mean(corner_means[4:])])
    #combine (variance of each detector with variance between detectors)
  if not False in usable:
    u_Cs=np.sqrt(tmp1**2+tmp1**2+tmpD**2)
  else:
    if usable[0]:
      u_Cs=tmp1
    elif usable[1]:
      u_Cs=tmp2
    elif not True in usable:
      u_Cs=2
  return u_all,k0,u_Cs,Cs_flag

def shot_noise(image,cf,Cs,srf,m2,debug):
  """
  calculates photonic noise for a given radiance value
  """
  debug=1
  if debug>=1:
    import ricals_tools as to
    
  radiance=cr.radiance(image,cf,Cs)
  if debug>=1:
    to.printimage(radiance,"radiance 4 shot noise [W/m^2/sr]",mini=np.amin(radiance[m2==1]),maxi=np.amax(radiance[m2==1]))
  #filterintegral
  srfint=np.trapz(srf[:,1],x=srf[:,0])#microns
  srfint=srfint*1E-6#meters
  if debug>=1:
    print "integral of srf ",srfint
    print "integral of srf ",np.trapz(np.ones(len(srf[:,0])),x=srf[:,0])
    print "min/max",np.amin(srf[:,0]),"/",np.amax(srf[:,0])
  #radiance conversion to spectral
  radiance=radiance/srfint
  if debug>=1:
    to.printimage(radiance,"spectral radiance 4 shot noise [W/m^2/sr/m^-1]",mini=np.amin(radiance[m2==1]),maxi=np.amax(radiance[m2==1]))

  #speed of light
  c          = 299792458 #m/s
  #Planck's constant
  h          = 6.62607004E-34 #m^2 kg/s
  #entrance aperture of detector
  radiusA    = (653./2.)/10./100. #m as given on page 142 in MTP.88D.304
  A          = np.pi*np.power(radiusA,2)#m^2 primary mirror
  if debug>=1:
    print "Area of primary mirror im m2 ",A
  #integration time
  EartWidth  = 2.*8.7 #degree as given on page 102 in MTP.88D.304
  RPM        = 100.
  REVOLUTION = 60./RPM #seconds per revolution
  LINE       = (REVOLUTION/360.)*EartWidth #seconds per line (approx 30ms)
  if debug>=1:
    print "line time in seconds ",LINE
  T          = LINE/2500. #seconds per pixel (integration time)
  if debug>=1:
    print "integration time in seconds ",T
  #solid angle
  Area       = 5000.*5000.#m^2 Area of sensor on Earth FIXME: should be changing with lat/lon
  r          = 36000000.#m satellite altitude
  Om         = Area / np.power(r,2) #m^2 / m^2
  Omsr       = (Om / 4*np.pi)*12.57 #sr
  if debug>=1:
    print "Solid angle in sr ",Omsr
  #number of photons
  N          =  (A*Omsr*T) / (h*c)  * radiance # ((m^2*sr*s) / (m^2*kg/s*m/s))* W/m^2/sr
  if debug>=1:
    to.printimage(N,"Nr of photons",mini=np.amin(N[m2==1]),maxi=np.amax(N[m2==1]))
    print N[2500,2500]
  #noise in number of photons
  dN         = np.sqrt(N)
  if debug>=1:
    print dN[2500,2500]
    to.printimage(dN,"Noise in Nr of photons",mini=np.nanmin(dN[m2==1]),maxi=np.nanmax(dN[m2==1]))
  #noise in radiance
  PnoiseRad  = dN/( (A*Omsr*T) / (h*c) )
  #radiance conversion to effective
  PnoiseRad=PnoiseRad*srfint
  if debug>=1:
    to.printimage(PnoiseRad,"photon shot noise in radiance",mini=np.nanmin(PnoiseRad[m2==1]),maxi=np.nanmax(PnoiseRad[m2==1]))
  #noise in counts
  Pnoise     = (PnoiseRad/cf)
  Pnoisefract= Pnoise/image
  if debug>=1:
    to.printimage(Pnoise,"photon shot noise in counts [absolute]",mini=np.nanmin(Pnoise[m2==1]),maxi=np.nanmax(Pnoise[m2==1]))
    to.printimage(Pnoisefract,"photon shot noise in counts [fractional]",mini=np.nanmin(Pnoisefract[m2==1]),maxi=np.nanmax(Pnoisefract[m2==1]))

  return Pnoise

def digit_noise(head,inflated=1):
    #calculate digitalisation noise
    if inflated==1:
    #if the old data are inflated:
      early=['M2','M3','P2']
      sat=getattr(head[0], 'szSatCode')
      if not sat in early:
        b=8 #8bit encoding
        u_dig=0.288675#https://www.wolframalpha.com/input/?i=standard+deviation+of+uniform+distribution+with+min%3D-0.5+and+max%3D0.5
      else:
        b=6 #6bit encoding
        u_dig=1.1547#https://www.wolframalpha.com/input/?i=standard+deviation+of+uniform+distribution+with+min%3D-2+and+max%3D2
      #diff_dig=0.5*np.divide(np.power(2,8),np.power(2,b))
    else:
    #if that is not the case:
      u_dig=0.288675#https://www.wolframalpha.com/input/?i=standard+deviation+of+uniform+distribution+with+min%3D-0.5+and+max%3D0.5
    return u_dig

def timing(lininfo,slot,v=1):
  '''
  fastest way to calculate acquisiton time precisely
  v indicates whether lininfo is expected from 
    image() (v=0) or image_3() (v=1) reader.
  '''
  t=((slot/2.)-0.5)*60*60
  line=np.arange(2500)
  u_t=[]
  im=[]
  mask=[]
  maskvis=[]
  miss=0
  nofit=0
  timscale=1
  for l in np.arange(0,2500):
    if v==0:
      flag=getattr(lininfo[l], 'shTimefitFlag')[0]
      coef=getattr(lininfo[l], 'dsFitCoef')
      frect=getattr(lininfo[l], 'shFirstRectPix')[0]
      ffit=getattr(lininfo[l], 'shFirstFitPix')[0]
      lrect=getattr(lininfo[l], 'shLastRectPix')[0]
      lfit=getattr(lininfo[l], 'shLastFitPix')[0]
      ftim=getattr(lininfo[l], 'dsTimeFirstRectPix')[0]
      ltim=getattr(lininfo[l], 'dsTimeLastRectPix')[0]
      maxdev=getattr(lininfo[l], 'dsMaxDevFit')[0]
    else:
      flag=getattr(lininfo[l], 'shTimefitFlag')
      coef=getattr(lininfo[l], 'dsFitCoef')
      frect=getattr(lininfo[l], 'shFirstRectPix')
      ffit=getattr(lininfo[l], 'shFirstFitPix')
      lrect=getattr(lininfo[l], 'shLastRectPix')
      lfit=getattr(lininfo[l], 'shLastFitPix')
      ftim=getattr(lininfo[l], 'dsTimeFirstRectPix')
      ltim=getattr(lininfo[l], 'dsTimeLastRectPix')
      maxdev=getattr(lininfo[l], 'dsMaxDevFit')
    linetimes=np.zeros((2500))
    linemask=np.zeros((2500))
    linemask[ffit:lfit]=1
    linemaskvis=np.zeros((2500*2))
    linemaskvis[ffit*2:lfit*2]=1
    linemask=linemask[::-1]
    linemaskvis=linemaskvis[::-1]
    if flag==5:
      try:
        #linetimes[ffit:lfit]=timscale*(coef[0] + coef[1]*np.power(line,2) + coef[2]*np.power(line,3) + coef[3]*np.power(line,4) + coef[4]*np.power(line,5))[ffit:lfit]
        linetimes[ffit:lfit]=timscale*np.add(coef[0] ,np.multiply( coef[1],line) + np.multiply(coef[2],np.power(line,2)) + np.multiply(coef[3],np.power(line,3)) + np.multiply(coef[4],np.power(line,4)))[ffit:lfit]
        linetimes[frect]=ftim
        linetimes[lrect]=ltim
        linetimes[frect:ffit]=np.linspace(linetimes[frect],linetimes[ffit],len(linetimes[frect:ffit]))
        linetimes[lfit-1:lrect]=np.linspace(linetimes[lfit-1],linetimes[lrect],len(linetimes[lfit-1:lrect]))
        linetimes=linetimes[::-1]
      except ValueError:
        miss=miss+1
      except TypeError:
        miss=miss+1
    else:
      nofit=nofit+1
    #fill rest(but no confidence in uncertainties-->should be flagged)
    linetimes[linetimes==0]=t+(l*((25.*60.)/2500))

    #store
    u_t.append(maxdev)
    im.append(linetimes)
    mask.append(linemask)
    maskvis.append(linemaskvis)
    maskvis.append(linemaskvis)
  #import matplotlib.pyplot as plt
  #plt.plot(np.array(im)[:,1250])
  #plt.show()
  #plt.plot(np.array(im)[1250,:])
  #plt.show()
  #print nofit,miss
  return np.divide(im,60*60),np.array(mask),np.array(maskvis),np.array(u_t)

def dz_easy(x,dx,axis=0):
  """
  calculate partial derivatives from 2D Array
  """
  if axis == 1:
    #y = np.zeros(x.shape, dtype='double')
    #y[:-1] = (x[1:] - x[:-1]) / dx
    #y[-1]  = -x[-1] / dx
    y=np.gradient(x, 1,axis=0) / dx
  else:
    #y = np.empty_like(x)
    #y[:, :-1] = (x[:, 1:] - x[:, :-1]) / dx
    #y[:, -1] = -x[:, -1] / dx
    y=np.gradient(x, 1,axis=1) / dx
  return y

def dz(pix,axis,d=1):
    """
    calculate partial derivatives from 2D Array
    """
    der = np.empty_like(pix)
    der2 = np.empty_like(pix)  
    if axis==0:
      der[:-1]    = np.diff(pix, axis=axis) / d
      der[-1]     = -pix[-1] / d
      pix2=np.fliplr(pix)
      der2[:-1]    = np.diff(pix2, axis=axis) / d
      der2[-1]     = -pix2[-1] / d
      der2=np.fliplr(der2)
    elif axis==1:
      der[:, :-1] = np.diff(pix, axis=axis) / d
      der[:, -1]  = -pix[:, -1] / d
      pix2=np.flipud(pix)
      der2[:, :-1] = np.diff(pix2, axis=axis) / d
      der2[:, -1]  = -pix2[:, -1] / d
      der2=np.flipud(der2)
    der=np.fabs(der)
    der2=np.fabs(der2)
    out=np.mean( np.array([der,der2]),axis=0)
    del der
    del der2
    return outplt.close()

def reflectance(k,slope,intercept,sza,msk,irr,sundist):
    radiance =  np.multiply ( slope  ,  np.subtract(  k  , intercept  ) ) #[W m^-2 sr^-1]
    z        =  np.radians(sza) #[rad]
    result   =  np.divide(np.multiply(np.pi  , np.multiply(radiance , np.power(sundist,2))),np.multiply(irr,np.cos(z)))
    try:
      result[k<=0]=-9.
      result[msk==0]=-9.
    except TypeError:
      print "WARNING: scalar count value provided?"
    return result
  
def geolocation(trail,pix,timestamp,sat,debug=0,return_buffer=False):
      #FIXME: read landmark thresholds for current nom SSP
      suspicious=0
      early=["MET2","MET3"]
      if sat in early:
        ul=64
      else:
        ul=255
      #detail analysis
      #print "---evaluating geolocation uncertainty in detail"
      py=[]
      px=[]
      valuesy=[]
      valuesx=[]
      bufflm=[]
      bump=getattr(trail[0], 'flAbsLandmarkResults')
      #pickle.dump(bump, open('./tmp/bump.p', 'wb'))
      #exit()
      #bump=pickle.load( open('./tmp/bump.p', 'rb'))
      for lm in range(0,128):#int(NirCorr[0]+2)):
        buff=bump[lm*8:((lm+1)*8)-1]#get 8 bytes at current offset
        if buff[2] < 2500. and buff[3] < 2500.:
          #FIXME: apply landmark thresholds
          py.append(buff[0])#lmrk x and y are lower right corner of 43x43 pixel landmark area
          px.append(buff[1])
          valuesy.append(buff[2])
          valuesx.append(buff[3])
          bufflm.append(buff)
      #if small number of landmarks, result may be wrong
      if len(valuesx)<5:
        suspicious=1
      #if std of landmarkdeviations is high it is suspicious
      if np.std(valuesx)>1. or np.std(valuesy)>1.:
        suspicious=1
      #compute rmse
      rmsX_vis=np.sqrt(((np.array(valuesx)*2) ** 2).mean())
      rmsY_vis=np.sqrt(((np.array(valuesy)*2) ** 2).mean())
      rmsX_ir=np.sqrt((np.array(valuesx) ** 2).mean())
      rmsY_ir=np.sqrt((np.array(valuesy) ** 2).mean())
      #if rmse of landmarks high geoloc is suspicious
      if rmsX_vis>1.5 or rmsY_vis>1.5:
        suspicious=1
      if debug==1:
        print "RMS of IR landmarks in X direction:"+str(rmsX_ir)
        print "RMS of IR landmarks in Y direction:"+str(rmsY_ir)
      
      #plot
      if debug>=2:
        to.printquiver(pix,np.add(px,21.5),np.add(py,21.5),np.multiply(valuesx,-1),valuesy,str(timestamp),mini=0,maxi=ul,save=0)
      #print rmsX_vis,rmsY_vis,rmsX_ir,rmsY_ir
      if return_buffer:
        return np.array(bufflm),np.column_stack((py,px,valuesy,valuesx))
      else:
        return [rmsX_vis,rmsY_vis,rmsX_ir,rmsY_ir],suspicious#np.array(bufflm),np.vstack((py,px,valuesy,valuesx))

def sensi_SZA(t,dt,DOY,m1,LAT,dLAT,LON,dLON,fill_value=0,debug=0):
  save=1
  if debug>=1:
    print "     relevant delta: ",np.amax((dt,dLAT,dLON))
  t1=np.ones((2500,2500))*(t/(60*60))
  #calculations
  if debug>=1:
    print "     calculate ref"
  SZA_ref,SAZ_ref=cr.sza(t1,DOY,m1,LAT,LON)
  del(SAZ_ref)
  if debug>=1:
    print "     calculate comp"
  SZA_tmp,SAZ_tmp=cr.sza(t1+(dt/(60*60)),DOY,m1,LAT+dLAT,LON+dLON)
  del(SAZ_tmp)
  if debug>=1:
    print "     calculate delta"
  s_sza_x=(SZA_tmp-SZA_ref)/np.amax((dt,dLAT,dLON))
  s_sza_x[SZA_ref<0]=fill_value
  if debug>=1:
    to.printimage(SZA_ref,str(int(t/(60*60))*2).zfill(2)+"_SZA_ref",mini=0,maxi=90,save=save)
    to.printimage(SZA_tmp,str(int(t/(60*60))*2).zfill(2)+"_SZA_tmp",mini=0,maxi=90,save=save)
    to.printimage(s_sza_x,str(int(t/(60*60))*2).zfill(2)+"_s_sza_x",mini=np.amin(s_sza_x),maxi=np.amax(s_sza_x),save=save)
  del(SZA_tmp)
  del(SZA_ref)
  return s_sza_x
      
def sensi_p(target,Ce,Cs,YSL,co,sza,E,d,Z=0):
  a0=co.a0
  a1=co.a1
  a2=co.a2
  sza=np.radians(sza)
  #p=pi, 
  #d=distance_sun_earth, 
  #f=solar_irradiance_vis, 
  #z=solar_zenith_angle, 
  #c=count_vis, 
  #s= mean_count_space_vis, 
  #a=a0_vis, 
  #b=a1_vis, 
  #t= years_since_launch, 
  #r=0
  def forCe(d,Ce,Cs,a0,a1,a2,YSL,sza,E,Z):
    #(d^2*p*(g*t^2+b*t+r+a))/(f*cos(z))
    return np.divide(np.power(d,2)*np.pi*(a2*np.power(YSL,2)+a1*YSL+a0+Z) ,\
                              np.cos(sza)*E)
  def forCs(d,Ce,Cs,a0,a1,a2,YSL,sza,E,Z):
    #-(d^2*p*(g*t^2+b*t+r+a))/(f*cos(z))
    return -np.divide(np.power(d,2)*np.pi*(a2*np.power(YSL,2)+a1*YSL+a0+Z) ,\
                              np.cos(sza)*E)
  def forE(d,Ce,Cs,a0,a1,a2,YSL,sza,E,Z):
    #-(d^2*p*(c-s)*(g*t^2+b*t+r+a))/(cos(z)*f^2)
    return -np.divide(np.power(d,2)*np.pi*(Ce-Cs)*(a2*np.power(YSL,2)+a1*YSL+a0+Z) ,\
                              np.cos(sza)*np.power(E,2))
  def forSZA(d,Ce,Cs,a0,a1,a2,YSL,sza,E,Z):
    #(d^2*p*(c-s)*(g*t^2+b*t+r+a)*sin(z))/(f*cos(z)^2)
    return np.divide(np.divide(np.power(d,2)*np.pi*(Ce-Cs)*(a2*np.power(YSL,2)+a1*YSL+a0+Z)*np.sin(sza),\
                              np.power(np.cos(sza),2)*E),180)*np.pi
  def fora0(d,Ce,Cs,a0,a1,a2,YSL,sza,E,Z):
    #(d^2*p*(c-s))/(f*cos(z))
    return np.divide(np.power(d,2)*np.pi*(Ce-Cs),\
                              np.cos(sza)*E)
  def fora1(d,Ce,Cs,a0,a1,a2,YSL,sza,E,Z):
    #(d^2*p*(c-s)*t)/(f*cos(z))
    return np.divide(np.power(d,2)*np.pi*(Ce-Cs)*YSL,\
                              np.cos(sza)*E)
  def fora2(d,Ce,Cs,a0,a1,a2,YSL,sza,E,Z):
    #(d^2*p*(c-s)*t^2)/(f*cos(z))
    return np.divide(np.power(d,2)*np.pi*(Ce-Cs)*np.power(YSL,2),\
                              np.cos(sza)*E)
  def forZ(d,Ce,Cs,a0,a1,a2,YSL,sza,E,Z):
    #(d^2*p*(c-s))/(f*cos(z))
    return np.divide(np.power(d,2)*np.pi*(Ce-Cs),\
                              np.cos(sza)*E)
  options= {  
              "Ce"  :  forCe,
              "digit"  :  forCe,
              "Cs"  :  forCs,
              "E"   :  forE,
              "SZA" :  forSZA,
              "a0"  :  fora0,
              "a1"  :  fora1,
              "a2"  :  fora2,
              "Z"   :  forZ,
            }
  return options[target](d,Ce,Cs,a0,a1,a2,YSL,sza,E,Z)
