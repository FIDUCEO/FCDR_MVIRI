"""
this file contains functions for calibration information
note that some effects may be implemented as fortran libs in folder lib/
"""
from datetime import date,datetime,timedelta
import numpy as np

import sys 
#import datetime
import scipy as sp
import netCDF4 as nc
import time
import timeit
import csv
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from collections import namedtuple
##from statsmodels.stats.weightstats import DescrStatsW
##import statsmodels.api as sm
from scipy import stats


import ricals_tools as to

#globals:
global debug
global SSCCversion
global parallel #flag to switch on/off parallel processing of multiple years
global N
global early

#default globals:
debug=0
N=5000
early            = ['M2','P2']
SSCCversion          = "4.01"#"4.01_101_60_V2"#
parallel         = "n"#flag to switch on(y)/off(n) parallel processing


def launchdate(sat):
  l_yr={"M2":1981,"M3":1988,"P2":1988,"M4":1989,"M5":1991,"M6":1993,"M7":1997}
  l_mo={"M2":6   ,"M3":6   ,"P2":6   ,"M4":3   ,"M5":3   ,"M6":11  ,"M7":9   }
  l_dy={"M2":19  ,"M3":15  ,"P2":15  ,"M4":2   ,"M5":2   ,"M6":20  ,"M7":2   }
  return str(l_yr[sat])+str(l_mo[sat]).zfill(2)+str(l_dy[sat]).zfill(2)

def decommissiondate(sat):
  l_yr={"M2":1988,"M3":1995,"P2":1995,"M4":1994,"M5":2007,"M6":2011,"M7":2017}
  l_mo={"M2":8   ,"M3":5   ,"P2":5   ,"M4":12  ,"M5":4   ,"M6":4   ,"M7":5   }
  l_dy={"M2":11  ,"M3":31  ,"P2":31  ,"M4":31   ,"M5":16  ,"M6":13  ,"M7":2   }
  return str(l_yr[sat])+str(l_mo[sat]).zfill(2)+str(l_dy[sat]).zfill(2)

def daysincelaunch(sat,timestring):
  ld=launchdate(sat)
  l=date(int(ld[0:4]),int(ld[4:6]),int(ld[6:8]))
  n=date(int(timestring[0:4]),int(timestring[4:6]),int(timestring[6:8]))
  return (n-l).days

def gainperiods(sat,gain):
  s_yr={"M2":{0:1981,1:1987},"M3":{1:1988,0:1990},"M4":{4:1989},"M5":{5:1991},"M6":{5:1996},"M7":{6:1997}}
  s_mo={"M2":{0:6   ,1:5   },"M3":{1:6   ,0:1   },"M4":{4:3   },"M5":{5:11  },"M6":{5:10  },"M7":{6:9   }}
  s_dy={"M2":{0:19  ,1:12  },"M3":{1:15  ,0:25  },"M4":{4:2   },"M5":{5:26  },"M6":{5:21  },"M7":{6:2   }}
  s_sl={"M2":{0:0   ,1:800 },"M3":{1:0   ,0:830 },"M4":{4:0   },"M5":{5:1230},"M6":{5:0   },"M7":{6:0   }}
  e_yr={"M2":{0:1987,1:1988},"M3":{1:1990,0:1991},"M4":{4:1994},"M5":{5:2007},"M6":{5:2011},"M7":{6:2017}}
  e_mo={"M2":{0:5   ,1:8   },"M3":{1:1   ,0:7   },"M4":{4:12  },"M5":{5:4   },"M6":{5:4   },"M7":{6:5   }}
  e_dy={"M2":{0:12  ,1:11  },"M3":{1:25  ,0:31  },"M4":{4:31  },"M5":{5:16  },"M6":{5:13  },"M7":{6:2   }}
  e_sl={"M2":{0:730 ,1:630 },"M3":{1:800 ,0:2330},"M4":{4:2330},"M5":{5:2330},"M6":{5:2330},"M7":{6:2330  }}
  #return str(s_yr[sat][int(gain)])+str(s_mo[sat][int(gain)]).zfill(2)+str(s_dy[sat][int(gain)]).zfill(2)+str(s_sl[sat][int(gain)]).zfill(2), \
         #str(e_yr[sat][int(gain)])+str(e_mo[sat][int(gain)]).zfill(2)+str(e_dy[sat][int(gain)]).zfill(2)+str(e_sl[sat][int(gain)]).zfill(2)
  return str(s_yr[sat][int(gain)])+str(s_mo[sat][int(gain)]).zfill(2)+str(s_dy[sat][int(gain)]).zfill(2)+str(s_sl[sat][int(gain)]).zfill(4), \
         str(e_yr[sat][int(gain)])+str(e_mo[sat][int(gain)]).zfill(2)+str(e_dy[sat][int(gain)]).zfill(2)+str(e_sl[sat][int(gain)]).zfill(4)

def gains(sat,datetimestring):
  for gain in range(16):
    print gainperiods(sat,gain), datetimestring
    if (int(gainperiods(sat,gain)[0])<= int(datetimestring)) and (int(datetimestring)<int(gainperiods(sat,gain)[1])):
      return str(gain)

def ssccfolder(sat,rel):
  if float(rel)<2.0:
    ret={ "M2":"/tcenas/proj/eraclim/fiduceo_tmp/SSCC_Archive/results/current_best/",\
          "M3":"/tcenas/proj/eraclim/fiduceo_tmp/SSCC_Archive/results/current_best/",\
          "M4":"/tcenas/proj/eraclim/fiduceo_tmp/SSCC_Archive/results/current_best/",\
          "M5":"/tcenas/proj/eraclim/fiduceo_tmp/SSCC_Archive/results/current_best/",\
          "M6":"/tcenas/proj/eraclim/fiduceo_tmp/SSCC_Archive/results/current_best/",\
          "M7":"/tcenas/proj/eraclim/fiduceo_tmp/SSCC_Archive/results/current_best/"}
  else:
    ret={ "M2":"/tcenas/proj/eraclim/fiduceo_tmp/SSCC_Archive/results/current_best_dynamicSRF/",\
          "M3":"/tcenas/proj/eraclim/fiduceo_tmp/SSCC_Archive/results/current_best_dynamicSRF/",\
          "M4":"/tcenas/proj/eraclim/fiduceo_tmp/SSCC_Archive/results/current_best_dynamicSRF/",\
          "M5":"/tcenas/proj/eraclim/fiduceo_tmp/SSCC_Archive/results/current_best_dynamicSRF/",\
          "M6":"/tcenas/proj/eraclim/fiduceo_tmp/SSCC_Archive/results/current_best_dynamicSRF/",\
          "M7":"/tcenas/proj/eraclim/fiduceo_tmp/SSCC_Archive/results/current_best_dynamicSRF/"}
  return ret[sat]

def cf(t,datetimestring,sat,gain,rel):
  #check gain consistency
  gainv1=gain[4:6]
  gainv2=gain[6:8]
  try:
    if int(gainv1)==int(gainv2):
      gainv=gainv1
    else:
      print "ERROR: gain settings inconsistent"
      print gain
      exit()
  except ValueError:
    gainv=gains(sat,datetimestring)
  calfolder="../config/"
  calfile=calfolder+sat+"_vis-calfile_R"+rel+"_gain"+gainv+".csv"
  try:
    print "reading calfile-VIS: "+calfile +" at t="+str(t) 
    f=np.loadtxt(calfile,skiprows=0)
    a=f[:,0]
    cov=f[:,1:]
    slope=a[0]+(a[1]*(t/365.))
    slopeu=np.sqrt(np.diag(cov)[0]+(t*np.sqrt(np.diag(cov)[0]))**2)
  except:
    if rel<=1.0:
      print "ERROR: calfile not available and for Release 1.0 no creation is foreseen!"
      exit()
    else:
      print " ...calfile not available; creating calfile"
      #s = launchdate(sat)
      s,e=gainperiods(sat,gainv)
      s=s[:-4]
      e=e[:-4]
      #e  = decommissiondate(sat)
      i  = ssccfolder(sat,rel)
      o  = calfolder
      satellite=to.satnames(sat)
      argsca=([s,e,i,o,satellite,calfile,1])
      if debug>=1:
        print "launch"
      a,cov=launch(argsca)
      if debug>=1:
        print "slope"
      slope=a[0]+(a[1]*(t/365.))
      if debug>=1:
        print "slopeu"
      slopeu=np.sqrt(np.diag(cov)[0]+(t*np.sqrt(np.diag(cov)[0]))**2)
      if debug>=1:
        print "return"
  return slopeu,slope,a,cov

#----------------------------------------------------------------------------
#get cf from SSCC results
#----------------------------------------------------------------------------
def launch((s,e,i,o,satellite,calfile,pool)): #args=[su,eu,a,f,sat,ye]
  """
  this launches the acquisition; the purpose of seperating 
  it from main was to allow parallel processing if desired
  """
  print "**************************"
  print "Pool "+str(pool)+" startdate  ",s
  print "Pool "+str(pool)+" enddate    ",e
  print "Pool "+str(pool)+" input      ",i
  print "Pool "+str(pool)+" output     ",o
  print "Pool "+str(pool)+" sat        ",satellite
  print "**************************"

  path=i
  satdir=satellite+"_"+SSCCversion+"/"
  sat=to.satnamesR(satellite)

  #get folderlist:
  print "searching results in: \n"+path+satdir+"*"
  resultdirs=sorted(glob.glob(path+satdir+"*"))
  
  #setup variable collectors
  dts   =[] #for dates
  dsl   =[] #for days since launch
  wcf_c =[] #combined desert cf
  unc   =[] #cf uncertainty
  
  #loop over the sscc-runs+++++++++++++++++++++++++++++++++++++++++++++++
  first=1
  for resultdir in sorted(resultdirs):
    
    #check date++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    year      = int(resultdir[-12:-8])
    startdoy  = int(resultdir[-7:-4])
    res=to.doy2date(year,startdoy)
    currdate=str(year)+str(res[1]).zfill(2)+str(res[0]).zfill(2)
    if debug>=1:
      print "resultdir ="+resultdir
      print "year      ="+str(year)
      print "doy       ="+str(startdoy)
      print "currdate  ="+str(currdate)
      print "startdate ="+str(s)
      print "enddate   ="+str(e)
      if not (int(currdate) >= int(s)) & (int(currdate) <= int(e)):
        print "...not in daterange; skipping....\n"

    if (int(currdate) >= int(s)) & (int(currdate) <= int(e)):

      #get all result files+++++++++++++++++++++++++++++++++++++++++++++
      resultfiles=glob.glob(resultdir+"/Cal_Results_MVIRI_BAND_VIS_DES*")
      if len(resultfiles)>0:

        
        #main task+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        #read the results to a dataframe (one df=one sscc run)
        if debug>=1:
          print resultdir
        df=read_summary(resultdir)
        if debug>=2:
          print df
        if not df.Q[0]==-1:
          #and append to dataframe collector
          if first==1:
            all_df=df#in first loop set it up
            first=0
          else:
            all_df=all_df.append(df)
          #handle date:
          dts.append(datetime(int(currdate[:4]), int(currdate[4:6]), int(currdate[6:8])))
          dsl.append(daysincelaunch(sat,currdate)/365.)
        #main task-----------------------------------------------------------------

  
  #prepare dates
  dts = np.array(dts)
  dsl = np.array(dsl)

  
  #get desert and sea df
  dfD=all_df[all_df.target=="DESERT"]
  dfS=all_df[all_df.target=="SEA"]

  
  #mask
  idx  = np.isfinite(dsl.astype(int))
  idx2 = np.isfinite(dsl.astype(int)) & (dfS.Cf.values.astype(float)>0)
  if debug>=1:
    print dsl[idx]
    print dfD.Cf.values.astype(float)[idx]
    print dfD.Cf
    print dfD.Cf.values.astype(float)
    #print dir(dfD.Cf)
  
  #fit drift
  y,dy,a,cov=linfit_err(dsl,dfD.Cf.values.astype(float),idx,dfD.CfErrT.values.astype(float))
  yS,dyS,aS,covS=linfit_err(dsl,dfS.Cf.values.astype(float),idx2,dfS.CfErrT.values.astype(float))
  cf0=a[0]
  dr =a[1]
  ucf0=np.sqrt(np.diag(cov)[0])
  udr=np.sqrt(np.diag(cov)[1])
  cf0S=aS[0]
  drS =aS[1]
  ucf0S=np.sqrt(np.diag(covS)[0])
  udrS=np.sqrt(np.diag(covS)[1])
  
  #store
  out=np.hstack((np.array([a,]).T,cov))
  outS=np.hstack((np.array([aS,]).T,covS))
  print np.shape(out)
  print out
  np.savetxt(calfile,out)
  
  #choose what is x-axis
  #temp=dsl.astype(int)
  temp=dts
  textstr="DESERT: Cf0=%.2f+/-%.2f Drift=%.4f+/-%.5f per yr\n\
  SEA   : Cf0=%.2f+/-%.2f Drift=%.4f+/-%.5f per yr" \
  %(cf0,ucf0,dr,udr,cf0S,ucf0S,drS,udrS)
  #plot the ts of the single value
  f, axarr = plt.subplots(1)
  axarr.xaxis.set_major_formatter(mdates.DateFormatter('%d/%m/%Y'))
  axarr.xaxis.set_major_locator(mdates.YearLocator())
  axarr.errorbar(temp[idx],dfD.Cf.values.astype(float)[idx],dfD.CfErrT.values.astype(float)[idx],fmt='.',color='black',label="individual cf Desert")
  axarr.errorbar(temp[idx2],dfS.Cf.values.astype(float)[idx2],dfS.CfErrT.values.astype(float)[idx2],fmt='.',color='red',label="individual cf Sea")
  axarr.fill_between(temp[idx], y[idx]-dy[idx], y[idx]+dy[idx],facecolor='green', alpha=0.5)
  axarr.fill_between(temp[idx], yS[idx]-dyS[idx], yS[idx]+dyS[idx],facecolor='red', alpha=0.5)
  axarr.plot(temp[idx],y[idx],'-',color="green",label="linear fit cf DESERT")
  axarr.plot(temp[idx],yS[idx],'-',color="red",label="linear fit cf SEA")
  axarr.grid()
  axarr.set_ylim((0.5,1.5))
  axarr.set_xlabel("date")
  axarr.set_ylabel("cf (Wm-2sr-1/DC)")
  axarr.legend(loc=4)
  props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
  axarr.text(0.05, 0.95, textstr, transform=axarr.transAxes, fontsize=14,
        verticalalignment='top', bbox=props)
  if debug>=1:
    plt.show()
  else:
    plt.savefig(calfile+".png", bbox_inches='tight')
    plt.close()
  return a,cov
#----------------------------------------------------------------------------
#read SSCC summary files
#----------------------------------------------------------------------------
def read_summary(s_folder,radiometer="MVIRI",band="VIS"):
  mpeffile=glob.glob(s_folder+"/*MPEF_out_*")[0]
  summaryfile=glob.glob(s_folder+"/*Cal_Results_Summary_"+radiometer+"_BAND_"+band+"*")[0]
  #quality check in mpef file
  with open(mpeffile,'r') as f:
    for line in f:
      if "First Julian Day :" in line:
        line=line.strip()
        columns=line.split()
        (DOY,YYYY)=columns[4].split("-")
      if "Overall Quality Flag" in line:
        if "Successfull" in line:
          q=1
        elif "Failed" in line:
          q=0
  #read summary file
  if q==1:
    with open(summaryfile,'r') as f:
      j=0
      col_titles="Q YYYYJJJ target spare1 Cf CfErrM CfErrE CfErrR CfErrT CfStd CfWeight".split()
      tuples=[]
      for line in f:
        line=line.strip()
        columns=line.split()
        if "Target Type DESERT :" in line:
          recID       = 'r'+str(j)
          struc       = namedtuple(recID, col_titles, verbose=False)
          struc       = struc(*columns)
          struc       = struc._replace(YYYYJJJ=YYYY+"/"+DOY)
          struc       = struc._replace(Q=1)
          j=j+1
          #print struc
          tuples.append(struc)
          
        elif "Target Type SEA    :" in line:
          recID       = 'r'+str(j)
          struc       = namedtuple(recID, col_titles, verbose=False)
          struc       = struc(*columns)
          struc       = struc._replace(YYYYJJJ=YYYY+"/"+DOY)
          struc       = struc._replace(Q=1)
          j=j+1
          #print struc
          tuples.append(struc)
      
    df=pd.DataFrame(tuples,columns=struc._fields)
    #print df
    dates=pd.to_datetime(df.YYYYJJJ, format="%Y/%j")
    df['DateTime']=pd.Series(dates,index=df.index)
    return df
  else:
    d = {'Q' : pd.Series([-1])}
    return pd.DataFrame(d)


#----------------------------------------------------------------------------
#Regression as implemented in GUI
#----------------------------------------------------------------------------
def  SSCC_GUI_TL_LREGRESS2(fx, fy, fw):
#; OUTPUT:
#;       a       = estimate of a
#;       b       = estimate of b
#;       sa      = estimate of a error standard deviation
#;       sb      = estimate of b error standard deviation
#;       sab     = stda and stdb quadratic sum
#;       r       = x y correlation coefficient
#;       chi2    = residual chi square 
  #convert fw as in GUI
  thr=sp.stats.t.cdf((0.95+1.)/2,len(fw))
  fw=1./(fw**2*thr**2)
  #regression
  nval = len(fx)
  fSx    = np.sum(fw*fx)
  fSxx   = np.sum((fw*np.power(fx,2)))
  fSy    = np.sum(fw*fy)
  fSxy   = np.sum(fw*fx*fy)
  fS     = np.sum(fw)
  fdd    = fS*fSxx - np.power(fSx,2)
  fa     = (fSxx*fSy - fSx*fSxy) / fdd
  fb     = (fS*fSxy - fSx*fSy) / fdd
  fsig_y = np.sqrt(1./(fS-2.)*np.sum(np.power((fy-fa-fb*fx),2)))
  fsa    = np.sqrt(fSxx/fdd)
  fsb    = np.sqrt(fS/fdd)
  fsab   = np.sqrt(np.power(fsa,2)+np.power(fsb,2))      
  fr     = fSx/np.sqrt(fS*fSxx)       
  fchi2  = np.sum(np.power((fw*(fy-fa-fb*fx)),2)) 
  return fa, fb, fsa, fsb, fsab, fr, fchi2

#----------------------------------------------------------------------------
#wrapper for linear regression
#----------------------------------------------------------------------------
def linfit_err(x,wcf_c,idx,use_u_i):
  #perform linear fit
  #as in GUI++++++++++:
  
  result=SSCC_GUI_TL_LREGRESS2(x[idx],wcf_c[idx],use_u_i[idx])
  #as in sm+++++++++++:
  ##X = sm.add_constant(x[idx])
  ##model=sm.WLS(wcf_c[idx],X,weights=use_u_i[idx])
  ##model=sm.OLS(wcf_c[idx],X)
  ##result2=model.fit()
  #if debug>=1:
  #  print(result2.summary())
  #  print result2.cov_params()
  #print np.sqrt(result2.cov_params())
  cf_0=result[0]
  cf_S=result[1]
  ucf_0=result[2]
  ucf_S=result[3]
  #cf_0=result2.params[0]
  #cf_S=result2.params[1]
  #ucf_0=result2.bse[0]
  #ucf_S=result2.bse[1]
  if debug>=1:
    print('Drift           : ', cf_S)  
    print('cf at launch    : ', cf_0)
    print('drift uncert    : ', ucf_S)
    print('drift uncert(%) : ', ucf_S/(cf_S/100))
    print('cf0 uncert      : ', ucf_0) 
    print('cf0 uncert(%)   : ', ucf_0/(cf_0/100))
  y=cf_0+cf_S*x
  dy=np.sqrt(ucf_0**2+(x*ucf_S)**2)
  a=[cf_0,cf_S]
  cov=np.diag([ucf_0**2,ucf_S**2])
  #a  =result2.params
  #cov=result2.cov_params()
  #print np.sqrt(int_2D(result2.cov_params(),[1,2],[1,2])[0])
  return y,dy,a,cov

