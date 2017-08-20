# -*- coding utf-8 -*- 
"""
this file contains the main program
"""

import sys, getopt
import numpy as np
from scipy.interpolate import Rbf
from scipy.interpolate import griddata
import netCDF4 as nc
import time
import timeit
import csv
import os
import glob
import shutil
import pandas as pd
#from sympy import mpmath as sm
import matplotlib.pyplot as plt

sys.path.insert(0, os.getcwd()[:-5]+'src/')
sys.path.insert(0, os.getcwd()[:-5]+'lib/nrCrunch/')

#import ricals_netcdf as wr
import ricals_effects as ef
import ricals_tools as to
#import ricals_calib as ca
#import ricals_irrad as ir
import read_rect2lp as rd
import read_imag2tg as r0

global badfiles
badfiles=[]
badfiles.append("/DSNNAS/Repro/mviri/level1/HR-MFG15/data/MET7/2000/08/METEOSAT7-MVIRI-MTP15-NA-NA-20000819053000.000000000Z")
badfiles.append("/DSNNAS/Repro/mviri/level1/HR-MFG15/data/MET7/2000/06/METEOSAT7-MVIRI-MTP15-NA-NA-20000601040000.000000000Z")

global parallel #flag to switch on/off parallel processing of multiple years
global N
global early
global inflated

#default globals
N        = 5000
early    = ['M2','M3','P2']
parallel = "n"
save     = 1
inflated = 1


def main(argv):
  #defaults
  global plt
  global debug 
  debug = 0

  f         = '/DSNNAS/Repro/mviri/level1/HR-MFG15/data/'
  a         = '../test/'
  t         = ''
  sorter   = 'line'
  irCALpath = '/tcc1/proj/eraclim/ricals/'
  
  try:
    opts, args = getopt.getopt(argv, "s:e:m:i:a:t:f:d:c:h",["startdate=","enddate=",\
                                                          "meteosat=",\
                                                          "irCAL=","archive=","time=", \
                                                          "filebuffer=","debug=", \
                                                          "sorter="])
  except getopt.GetoptError:
    print "ERROR"
    print "\n correct usage:"
    print '\nricals_main.py -s <startdate as yyyymmdd or yyyydoy> '
    print '                 -e <enddate as yyyymmdd or yyyydoy> '
    print '                 -m <meteosat such as MET7> '
    print '                 (-a <path to archive folder> )'
    print '                 (-t <specific time as HHMM> )'
    print '                 (-f <path to filebuffer> )'
    print '                 (-d <debug mode>)' 
    print '                 (-c <sorter>)'

    sys.exit(2)
  if len(opts)<4:
    opts=[["-h",""]]
  
  for opt, arg in opts:
    if opt == '-h':
        print '\nprogram to automatically run RICalPy on MVIRI data for a period between date S and date E'
        print "\nusage:"
        print '\nricals_main.py -s <startdate as yyyymmdd or yyyydoy> '
        print '                 -e <enddate as yyyymmdd or yyyydoy> '
        print '                 -m <meteosat such as MET7> '
        print '                 (-a <path to archive folder> )'
        print '                 (-t <specific time as HHMM> )'
        print '                 (-f <path to filebuffer> )'
        print '                 (-d <debug mode>)' 
        print '                 (-c <sorter>)'
        sys.exit()
    elif opt in ("-s", "--startdate"):
        s = str(arg)
    elif opt in ("-e", "--enddate"):
        e = str(arg)
    elif opt in ("-m", "--meteosat"):
        sat = str(arg)
    elif opt in ("-i", "--irCAL"):
        irCALpath = str(arg)
    elif opt in ("-a", "--archive"):
        a = str(arg)
    elif opt in ("-t", "--time"):
        t = str(arg)
    elif opt in ("-f", "--filebuffer"):
        f = str(arg)
    elif opt in ("-d", "--debug"):
        debug = int(str(arg))
    elif opt in ("-c", "--sorter"):
        sorter = arg

  
  if len(s)<8: #If the date is given as doy, the string of the date will be smaller 8
    res=to.doy2date(int(s[:4]),int(s[4:]))
    s=s[:4]+str(res[1]).zfill(2)+str(res[0]).zfill(2)
    
  if len(e)<8: #If the date is given as doy, the string of the date will be smaller 8
    res=to.doy2date(int(e[:4]),int(e[4:]))
    e=e[:4]+str(res[1]).zfill(2)+str(res[0]).zfill(2)

  args=([s,e,a,t,f,sat,1,sorter,irCALpath])
  timing=[]
  timing.append(timeit.default_timer()) #0
  launch(args) #standard launch without parallel threads
  timing.append(timeit.default_timer()) #1
  print "this call took "+str(timing[1]-timing[0])+" seconds"
  
def launch((s,e,a,t,fstr,satellite,pool,sorter,irCALpath)): #args=[su,eu,a,f,sat,ye]

  print "**************************"
  print "Pool "+str(pool)+" startdate ",s
  print "Pool "+str(pool)+" enddate   ",e
  print "Pool "+str(pool)+" archive   ",a
  print "Pool "+str(pool)+" filebuffer",fstr
  print "Pool "+str(pool)+" sat       ",satellite
  print "**************************"

#===============================================================================================
# Part 1 :Setup parameters
#-----------------------------------------------------------------------------------------------  
  startyear             = int(s[:4])
  endyear               = int(e[:4])
  startdoy              = to.date2doy(startyear,int(s[4:6]),int(s[6:8]))
  enddoy                = to.date2doy(  endyear,int(e[4:6]),int(e[6:8]))
  Nd                    = ((((endyear)-startyear)*365)+enddoy-startdoy)+1
  sensor                = "MVIRI"
#===============================================================================================
# Part 2 :Setup landmarks
#-----------------------------------------------------------------------------------------------  
  lmrkfiles=[]
  lmrkfiles.append("Landmarks/LMRKLIST.E000")
  #lmrkfiles.append("Landmarks/LMRKLIST.E057")
  #lmrkfiles.append("Landmarks/LMRKLIST.E063")

  lmrklist=pd.concat([pd.read_csv(lmrkfile, header=None,delim_whitespace=True,usecols=[0,1,2]) for lmrkfile in lmrkfiles])
  lmrklist.columns = ["Name", "line", "pixel"]
  notsorter="pixel"
  if sorter=="pixel":
    notsorter="line"
  lmrklist=lmrklist.sort(sorter)
  #lmrklist=lmrklist.drop_duplicates(cols=["Name"])
  index     = np.arange(0,Nd*48)
  resultsetx = pd.DataFrame(index=index, columns=lmrklist["Name"])
  resultsety = pd.DataFrame(index=index, columns=lmrklist["Name"])
  
  lmrktables=[]
  for lmrkfile in lmrkfiles:
    tmp=pd.read_csv(lmrkfile, header=None,delim_whitespace=True)
    tmp.columns = ["Name", "line", "pixel"]
    lmrktables.append(tmp.sort(sorter))
  
  print lmrktables
#===============================================================================================
# Part 3 :Start looping
#-----------------------------------------------------------------------------------------------  
  if not os.path.exists("../test/resultsety.txt"):
    idx=0
    for year in range(startyear,endyear+1):
      
      if year == startyear: # first year starts somewhere
        sd=startdoy
      else:                # the other years start on first day
        sd=1
        
      if year < endyear:   # consider the year to always have 365 days
        ed=365
      else:                # except for the last year
        ed=enddoy
        
      s_Year      = str(year)
        
      for doy in range(sd,ed+1):
        s_Doy            = str(doy).zfill(3)
        i_Day,i_Month    = to.doy2date(year,doy)
        s_Month          = str(i_Month).zfill(2)
        s_Day            = str(i_Day).zfill(2)
        s_datestring = s_Year+s_Month+s_Day
        #check if in iodc service
        iodc=0
        if "MET5" in satellite and int(s_Year)==1998 and int(s_Doy)>182:
          iodc=1
        elif "MET5" in satellite and int(s_Year)>1998:
          iodc=1
        if "MET7" in satellite and int(s_Year)==2006 and int(s_Doy)>305:
          iodc=1
        elif "MET7" in satellite and int(s_Year)>2006:
          iodc=1
        print "iodc:       ",iodc
        print "debug mode: ",str(debug)+ " \n"

        filepath   =  fstr+'/'+satellite+'/'+s_Year+'/'+s_Month+'/'
        if debug>=2:
          print "searching "+filepath
        if os.path.exists(filepath):
          if t=='':
            wildcard=filepath+"*"+s_datestring+"*"
          else:
            wildcard=filepath+"*"+s_datestring+t+"*"
          fod      = sorted(glob.glob(wildcard))
          if debug>=2:
            print "wildcard is: "+wildcard
            
        else:
          fod      = -1
          print "WARNING: not existing "+filepath
        
        if len(fod)==0:
          print "\nERROR: file not existing for wildcard: "+wildcard+"\n"
        

        for fname in fod:
          if not fname in badfiles:
            print fname
            exists=False
            if debug>=1:
              print "infile: "+fname
            i_Tindex =  fname.index(s_datestring)+len(s_datestring)
            s_Time   =  fname[i_Tindex:i_Tindex+4]
            
            filename = fname #s_datestring+s_Time
            
            if debug==1:
              print "time   : "+s_Time
            
    #===============================================================================================
    # Part 2 :Get Data
    #-----------------------------------------------------------------------------------------------  
            
            print "Reading data..."
            #read rect2lp++++++++++
            try:
              with open(fname, "rb") as f:
                f, head = rd.header(f, 1, 3,t=0)
                #f, lininfo, pixvis, pixir, pixwv = rd.image(f,head,4,2503, t=0)
                f, trail = rd.header(f, 2504, 2504, t=0)
              
              sat="".join(getattr(head[0], 'szSatCode'))
              tloc=fname.index(s_datestring)
              timestring=fname[tloc:tloc+12]
              nomlon = np.rad2deg(getattr(head[1],'dsNominalSubSatelliteLongitude')[0])
              lonidx={0:0,57:1,63:2}[int(nomlon)]
              lm_def=ef.geolocation(trail,np.zeros((2500,2500)),timestring,sat,debug,return_buffer=True)
              for line in lm_def[1]:
                #get name
                sset     = lmrktables[lonidx].loc[lmrktables[lonidx]['line'] == line[0]]
                try:
                  name     = (sset.loc[sset['pixel'] == line[1]]['Name']).values[0]
                  #print name
                  #if name=="XAAFUM":
                    #print idx,line[3]
                  #put in resultset
                  #print resultsetx[idx]
                  resultsetx.set_value(idx, name, line[3])
                  resultsety.set_value(idx, name, line[2])
                except IndexError:
                  print "WARNING: No entry for: " + str(line)
              #increment slot index
              idx=idx+1
            except:
              idx=idx+1
          else:
            idx=idx+1
    
    resultsetx.to_csv("../test/resultsetx.txt", mode='w')
    resultsety.to_csv("../test/resultsety.txt", mode='w')
  else:
    print "just reading ../test/resultsetx.txt and ../test/resultsety.txt"
    resultsetx=pd.read_csv("../test/resultsetx.txt")
    resultsety=pd.read_csv("../test/resultsety.txt")
  
  #explore covariances
  if not os.path.exists('../test/'+sorter+'/COV_X.txt'):
    COV_X=np.zeros((len(lmrklist["Name"]),len(lmrklist["Name"])))
    COV_Y=np.zeros((len(lmrklist["Name"]),len(lmrklist["Name"])))
    DIST=np.zeros((len(lmrklist["Name"]),len(lmrklist["Name"])))
    DIST2=np.zeros((len(lmrklist["Name"]),len(lmrklist["Name"])))
    print np.shape(COV_X)
    i=0
    for n in lmrklist["Name"]:
      j=0
      for m in lmrklist["Name"]:
        loc1=lmrklist.loc[lmrklist['Name'] == n][sorter].values[0]
        loc2=lmrklist.loc[lmrklist['Name'] == m][sorter].values[0]
        DIST[i,j]=abs(loc1-loc2)
        loc1=lmrklist.loc[lmrklist['Name'] == n][notsorter].values[0]
        loc2=lmrklist.loc[lmrklist['Name'] == m][notsorter].values[0]
        DIST2[i,j]=abs(loc1-loc2)
        a=resultsetx[n]
        b=resultsetx[m]
        a=pd.concat([a, b], axis=1).dropna(axis=0,thresh=2).as_matrix()
        c=resultsety[n]
        d=resultsety[m]
        c=pd.concat([c, d], axis=1).dropna(axis=0,thresh=2).as_matrix()
        try:
          covax =np.corrcoef(a.astype(float), rowvar=False)#
          covay =np.corrcoef(c.astype(float), rowvar=False)
          try:
            COV_X[i,j]=covax[0][1]
            COV_Y[i,j]=covay[0][1]
          except IndexError:
            COV_X[i,j]=0
            COV_Y[i,j]=0
        except TypeError:
          COV_X[i,j]=0
          COV_Y[i,j]=0
        j=j+1
      i=i+1
    DISTANCE=np.sqrt(np.power(DIST,2)+np.power(DIST2,2))
    np.savetxt('../test/'+sorter+'/DISTANCE.txt',DISTANCE , delimiter=';')
    np.savetxt('../test/'+sorter+'/DIST_'+sorter+'.txt',DIST , delimiter=';')
    np.savetxt('../test/'+sorter+'/DIST_'+notsorter+'.txt',DIST2 , delimiter=';')
    np.savetxt('../test/'+sorter+'/COV_X.txt',COV_X , delimiter=';')
    np.savetxt('../test/'+sorter+'/COV_Y.txt',COV_Y , delimiter=';')
  else:
    print "just reading covariances"
    DISTANCE= np.genfromtxt('../test/'+sorter+'/DISTANCE.txt', delimiter=';')
    DIST    = np.genfromtxt('../test/'+sorter+'/DIST_'+sorter+'.txt', delimiter=';')
    DIST2   = np.genfromtxt('../test/'+sorter+'/DIST_'+notsorter+'.txt', delimiter=';')
    COV_X   = np.genfromtxt('../test/'+sorter+'/COV_X.txt', delimiter=';')
    COV_Y   = np.genfromtxt('../test/'+sorter+'/COV_Y.txt', delimiter=';')
  

  grid_x, grid_y = np.mgrid[0:2500:2500j, 0:2500:2500j]
  points=np.column_stack((DIST2[~np.isnan(COV_X)],DIST[~np.isnan(COV_X)]))
  grid_z0 = griddata(points, COV_X[~np.isnan(COV_X)], (grid_x.astype(int), grid_y.astype(int)), method='linear',rescale=False)
  plt.imshow(grid_z0.T, extent=(0,2500,0,2500), origin='lower')
  plt.colorbar()
  plt.xlabel(notsorter)
  plt.ylabel(sorter)
  plt.plot(points[:,0], points[:,1], 'k.', ms=1)
  plt.title("Correlation of error in longitude-direction for landmark pairs, plotted as a function of their distance in line and pixel direction")
  #plt.scatter(points[:,0], points[:,1],c=np.power(COV_X[~np.isnan(COV_X)],1))
  plt.show()
  points=np.column_stack((DIST2[~np.isnan(COV_Y)],DIST[~np.isnan(COV_Y)]))
  grid_z0 = griddata(points, COV_Y[~np.isnan(COV_Y)], (grid_x.astype(int), grid_y.astype(int)), method='linear',rescale=False)
  plt.imshow(grid_z0.T, extent=(0,2500,0,2500), origin='lower')
  plt.colorbar()
  plt.xlabel(notsorter)
  plt.ylabel(sorter)
  plt.plot(points[:,0], points[:,1], 'k.', ms=1)
  plt.title("Correlation of error in latitude-direction for landmark pairs, plotted as a function of their distance in line and pixel direction")
  #plt.scatter(points[:,0], points[:,1],c=np.power(COV_Y[~np.isnan(COV_Y)],1))
  plt.show()
  
  #plt.scatter(DIST[COV_Y!=0],DIST2[COV_Y!=0],c=np.power(COV_Y[COV_Y!=0],1))
  #plt.title("COV_X vs. LMRK distance")
  #plt.xlabel(sorter)
  #plt.ylabel(notsorter)
  #plt.colorbar()
  #plt.show()
  
  #plt.scatter(DIST[COV_Y!=0],DIST2[COV_Y!=0],c=np.power(COV_Y[COV_Y!=0],1))
  #plt.title("COV_Y vs. LMRK distance")
  #plt.xlabel(sorter)
  #plt.ylabel(notsorter)
  #plt.colorbar()
  #plt.show()
  
  plt.plot(np.ndarray.flatten(DIST[COV_X!=0]),np.power(np.ndarray.flatten(COV_X[COV_X!=0]),1),'.')
  plt.title("COV_X vs. inter-LMRK distance")
  plt.xlabel("distance")
  plt.ylabel("Correlation")
  plt.show()
  
  plt.plot(np.ndarray.flatten(DISTANCE[COV_Y!=0]),np.power(np.ndarray.flatten(COV_Y[COV_Y!=0]),1),'.')
  plt.title("COV_Y vs. inter-LMRK distance")
  plt.xlabel("distance")
  plt.ylabel("Correlation")
  plt.show()
  print np.nanmax(DIST)
  plt.imshow(DIST)
  plt.title("DISTANCES in direction "+sorter)
  plt.colorbar()
  plt.show()
  print np.nanmax(DISTANCE)
  plt.imshow(DISTANCE)
  plt.title("DISTANCES both directions")
  plt.colorbar()
  plt.show()
  print np.nanmax(COV_X)
  plt.imshow(COV_X)
  plt.title("COV MATR X")
  plt.colorbar()
  plt.show()
  print np.nanmax(COV_Y)
  plt.imshow(COV_Y)
  plt.title("COV MATR Y")
  plt.colorbar()
  plt.show()
  
  #explore autocorrelation
  for nm in lmrklist["Name"]:
    print nm
    #print resultsetx[nm].dropna()
    #print resultsetx.index.values[:]
    #autocorrelation
    acx=autocorr(resultsetx[nm].values[:])
    acy=autocorr(resultsety[nm].values[:])
    #plot
    try:
      if len(resultsetx[nm].dropna())>1:
        plt.plot(resultsetx.index.values[:],resultsetx[nm].values[:],'.')
        plt.title(nm+" autocorrelation X: "+str(acx))
        plt.xlabel("slots")
        plt.ylabel("Landmark deviation [pixel]")
        plt.savefig("../test/LMRK_"+nm+"X.png", bbox_inches='tight')
        plt.close()
    except TypeError:
      print "WARNING: no results for this target in X"
    try:
      if len(resultsety[nm].dropna())>1:
        plt.plot(resultsety.index.values[:],resultsety[nm].values[:],'.')
        plt.title(nm+" autocorrelation Y: "+str(acy))
        plt.xlabel("slots")
        plt.ylabel("Landmark deviation [pixel]")
        plt.savefig("../test/LMRK_"+nm+"Y.png", bbox_inches='tight')
        plt.close()
    except TypeError:
      print "WARNING: no results for this target "
  

def autocorr(x, t=1):
    x=np.array([x[0:len(x)-t],x[t:len(x)]])
    a=pd.DataFrame({'x':x[:,0], 'y':x[:,1]})
    #x=a.dropna(axis=0,thresh=2).as_matrix()
    x=a.as_matrix()
    if len(x)>5:
      cov=np.corrcoef(x.astype(float), rowvar=False)
      print cov[0][1]
      return cov[0][1]
    else:
      return 0

if __name__ == "__main__":
   main(sys.argv[1:])