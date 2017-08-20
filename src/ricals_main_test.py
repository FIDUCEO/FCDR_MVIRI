# -*- coding utf-8 -*- 
"""
this file contains the main program
"""

import sys, getopt
import numpy as np
import netCDF4 as nc
import time
import timeit
import csv
import os
import glob
import shutil
#from sympy import mpmath as sm
import matplotlib.pyplot as plt

import ricals_netcdf as wr
import ricals_effects as ef
import ricals_tools as to
import ricals_calib as ca
import ricals_irrad as ir
import read_rect2lp as rd
import read_imag2tg as r0

from fiduceo.fcdr.writer.fcdr_writer import FCDRWriter


sys.path.insert(0, os.getcwd()[:-4]+'/lib/nrCrunch/')
import cruncher as cr

global Software
global version
global parallel #flag to switch on/off parallel processing of multiple years
global N
global early
global inflated

#default globals
Software = "2.1"
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
  comment   = 'created without comment'
  irCALpath = '/tcc1/proj/eraclim/ricals/'
  
  try:
    opts, args = getopt.getopt(argv, "s:e:m:r:v:i:a:t:f:d:c:h",["startdate=","enddate=",\
                                                          "meteosat=","release=","srfVersion=",\
                                                          "irCAL=","archive=","time=", \
                                                          "filebuffer=","debug=", \
                                                          "comment="])
  except getopt.GetoptError:
    print "ERROR"
    print "\n correct usage:"
    print '\nricals_main.py -s <startdate as yyyymmdd or yyyydoy> '
    print '                 -e <enddate as yyyymmdd or yyyydoy> '
    print '                 -m <meteosat such as MET7> '
    print '                 -r <release (=format version; Software: '+Software+'>)'
    print '                 -v <srfversion (=srf version; e.g. 0954>)'
    print '                 (-i <path to ir-calibration folder> )'
    print '                 (-a <path to archive folder> )'
    print '                 (-t <specific time as HHMM> )'
    print '                 (-f <path to filebuffer> )'
    print '                 (-d <debug mode>)' 
    print '                 (-c <comment>)'

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
        print '                 -r <release (=format version; Software: '+Software+'>)'
        print '                 -v <srfversion (=srf version; e.g. 0954>)'
        print '                     (-i <path to ir-calibration folder> )'
        print '                 (-a <path to archive folder> )'
        print '                 (-t <specific time as HHMM> )'
        print '                 (-f <path to filebuffer> )'
        print '                 (-d <debug mode>)' 
        print '                 (-c <comment>)'
        sys.exit()
    elif opt in ("-s", "--startdate"):
        s = str(arg)
    elif opt in ("-e", "--enddate"):
        e = str(arg)
    elif opt in ("-m", "--meteosat"):
        sat = str(arg)
    elif opt in ("-r", "--release"):
        rel = arg
    elif opt in ("-v", "--srfVersion"):
        srf_vers = arg
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
    elif opt in ("-c", "--comment"):
        comment = arg
        
        if debug==1:
          print "importing matplotlib..."
          try:
            import matplotlib.pyplot as plt
          except:
            print "WARNING: matplotlib not installed"
  
  if float(rel)<2.0:
    try:
      if not srf_vers in ["pre-launch","pre_launch","pre-flight","pre_flight"]:
        print "ERROR: Release "+rel+" is not meant to contain reconstructed SRF"
        return
    except:
      srf_vers="pre-launch"
  
  if len(s)<8: #If the date is given as doy, the string of the date will be smaller 8
    res=to.doy2date(int(s[:4]),int(s[4:]))
    s=s[:4]+str(res[1]).zfill(2)+str(res[0]).zfill(2)
    
  if len(e)<8: #If the date is given as doy, the string of the date will be smaller 8
    res=to.doy2date(int(e[:4]),int(e[4:]))
    e=e[:4]+str(res[1]).zfill(2)+str(res[0]).zfill(2)

  if parallel=="n":
    args=([s,e,a,t,f,sat,1,comment,rel,srf_vers,irCALpath])
    timing=[]
    timing.append(timeit.default_timer()) #0
    launch(args) #standard launch without parallel threads
    timing.append(timeit.default_timer()) #1
    print "this call took "+str(timing[1]-timing[0])+" seconds"
  else:
    nr_years=int(e[:4])-int(s[:4])+1
    args=[]
    if nr_years <=1: #parallel threads only if multiple years
      print "************************************************************************"
      print "WARNING: too little years for parallel processing; doing all in one pool"
      print "************************************************************************"
      args=([s,e,a,t,f,sat,1,comment,rel,srf_vers,irCALpath])
      timing=[]
      timing.append(timeit.default_timer()) #0
      launch(args)
      timing.append(timeit.default_timer()) #1
      print "this call took "+str(timing[1]-timing[0])+" seconds"
    else:
      print "************************************************************************"
      print "             parallel processing in "+str(nr_years)+" pools"
      print "************************************************************************"
      for ye in range(1,nr_years+1):#prepare arguments for the runs
        if ye == 1:
          su =s
          end=to.doy2date(int(s[:4]),365)
          eu =s[:4]+str(end[1]).zfill(2)+str(end[0]).zfill(2)
        elif ye==nr_years:
          eu =e
          su =e[:4]+"0101"
        else:
          yu=str(int(s[:4])+(ye-1)).zfill(4)
          su=yu+"0101"
          end=to.doy2date(int(yu),365)
          eu=yu+str(end[1]).zfill(2)+str(end[0]).zfill(2)

        args.append((su,eu,a,t,f,sat,ye,comment,rel,srf_vers,irCALpath))
        print args
    
      pool = multiprocessing.Pool(processes=nr_years)
      pool.map(launch,args)



def launch((s,e,a,t,f,satellite,pool,comment,rel,srf_vers,irCALpath)): #args=[su,eu,a,f,sat,ye]

  print "**************************"
  print "Pool "+str(pool)+" startdate ",s
  print "Pool "+str(pool)+" enddate   ",e
  print "Pool "+str(pool)+" archive   ",a
  print "Pool "+str(pool)+" filebuffer",f
  print "Pool "+str(pool)+" sat       ",satellite
  print "**************************"

#===============================================================================================
# Part 1 :Setup parameters
#-----------------------------------------------------------------------------------------------  
  startyear             = int(s[:4])
  endyear               = int(e[:4])
  startdoy              = to.date2doy(startyear,int(s[4:6]),int(s[6:8]))
  enddoy                = to.date2doy(endyear,int(e[4:6]),int(e[6:8]))
  
  sensor           = "MVIRI"
  versionfolder    = "MFG15_FCDR_V"+rel+"_SRF"+srf_vers

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

      filepath   =  f+'/'+satellite+'/'+s_Year+'/'+s_Month+'/'
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
      
      staticdir  =  a + versionfolder+'/aux/'
      lutdir     =  a + versionfolder+'/aux/'
      outdir     =  a + versionfolder+'/data/'+satellite+'/'+s_Year+'/'+s_Month+'/'
      if not os.path.exists(outdir):
        os.makedirs(outdir)

      for fname in fod:
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
        with open(fname, "rb") as f:
            f, head = rd.header(f, 1, 3,t=0)
            nomlon = np.rad2deg(getattr(head[1],'dsNominalSubSatelliteLongitude')[0])
            f, lininfo, pixvis, pixir, pixwv = rd.image(f,head,4,2503, t=0)
            f, trail = rd.header(f, 2504, 2504, t=0)

        t1,m1,m2,ut=ef.timing(lininfo)
 
        #LAT,LON   = cr.latlon(N)
        #LON       = np.add(LON,nomlon)
        #DOY=getattr(head[0],'iYday')[0]
        #SZA,SAZ=cr.sza(t1,DOY,m1,LAT,LON)
        #sslat, sslon, sslatE, sslonE = to.calcSSP(head,head,nomlon,early)
        #sslatM=(sslat+sslatE)/2
        #sslonM=(sslon+sslonE)/2
        #cont_arr=SZA
        #cont_arr[cont_arr<0]=np.nan
        #cont_arr[cont_arr>100]=np.nan
        #to.printcontours(sslatM,sslonM,LAT,LON,cont_arr,pixvis,m2,"test (absolute)",save=0)

        tits=["BRF","random uncertainty","non-random uncertainty","diff non-random - random"]
        pixs  =[pixvis,pixvis,pixvis,pixvis]
        to.print4mat(pixs,tits,m2,save=0)
        
if __name__ == "__main__":
   main(sys.argv[1:])
