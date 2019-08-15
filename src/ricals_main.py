# -*- coding utf-8 -*- 
"""
this file contains the main program
"""

import sys, getopt
import numpy as np
import netCDF4 as nc
import time
import timeit
import datetime as dt
import csv
import os
import glob
import shutil
#from sympy import mpmath as sm
import matplotlib.pyplot as plt
import functools

from functools import partial

import cruncher as cr #grid definition as described in mfg user manual

#import psycopg2
#from guppy import hpy
import ricals_netcdf as wr
import ricals_effects as ef
import ricals_tools as to
#import ricals_calib as ca
import calib as ca
import ricals_irrad as ir

from mviri.mviri_l10 import read_imag2tg as r0
from mviri.mviri_l15 import read_rect2lp as rd

from fiduceo.fcdr.writer.fcdr_writer import FCDRWriter
from fiduceo.common.writer.default_data import DefaultData

#warnings.filterwarnings("ignore",category =RuntimeWarning)

global grid
grid="mfg_def"

global Software
global version
global parallel #flag to switch on/off parallel processing of multiple years
global N
global early
global inflated
global tstsite

#default globals
Software = "2.6"
N        = 5000
early    = ['M2','M3','P2']
parallel = "n"
save     = 0
inflated = 1
#tstsite   = [17.96,-18.84] #North atlantic; lat,lon
#tstsite   = [38.21,-27.42] #azores; lat,lon
#tstsite   = [23.8,-0.40]   #alg1; lat,lon
tstsite   = [28.55,23.39]  #lyb4; lat,lon
#tstsite = [0.,0.]

def main(argv):
  #defaults
  global plt
  global debug 
  debug = 0

  f         = '/DSNNAS/Repro/mviri/level1/HR-MFG15/data/'
  a         = '../test/'
  t         = ''
  comment   = 'created without comment'
  irCALpath = '/tcenas/home/vijuj/work/satcal/mike2check/'#'/DSNNAS/Repro/aux/ricals/calibration_coefficients/'
  # launch
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
        sys.exit(1)
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
        debug = float(str(arg))
    elif opt in ("-c", "--comment"):
        comment = arg
        if debug==1:
          print "importing matplotlib..."
          try:
            import matplotlib.pyplot as plt
          except:
            print "WARNING: matplotlib not installed"
  
  #logging DB init
  logs=to.logger(rel)
  err=logs.createDB()
  if err==0: #in case DB was newly created
    err=logs.create_Table()

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
    args=([s,e,a,t,f,sat,1,comment,rel,srf_vers,irCALpath,logs])
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
      args=([s,e,a,t,f,sat,1,comment,rel,srf_vers,irCALpath,logs])
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

        args.append((su,eu,a,t,f,sat,ye,comment,rel,srf_vers,irCALpath,logs))
        print args
    
      pool = multiprocessing.Pool(processes=nr_years)
      pool.map(launch,args)


#@profile
def launch((s,e,a,t,f,satellite,pool,comment,rel,srf_vers,irCALpath,logs)): #args=[su,eu,a,f,sat,ye]
  #hp = hpy()
  #print "Heap at the beginning of the function\n", hp.heap()
  print "**************************"
  print "Pool "+str(pool)+" startdate ",s
  print "Pool "+str(pool)+" enddate   ",e
  print "Pool "+str(pool)+" archive   ",a
  print "Pool "+str(pool)+" filebuffer",f
  print "Pool "+str(pool)+" sat       ",satellite
  print "Pool "+str(pool)+" t         ",t
  print "**************************"

#===============================================================================================
# Part 1 :Setup parameters
#-----------------------------------------------------------------------------------------------  
  
  startyear             = int(s[:4])
  endyear               = int(e[:4])
  startdoy              = to.date2doy(startyear,int(s[4:6]),int(s[6:8]))
  enddoy                = to.date2doy(endyear,int(e[4:6]),int(e[6:8]))
  
  sensor           = "MVIRI"
  versionfolder    = ""#"MFG15_FCDR_V"+rel+"_SRF"+srf_vers

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
      s_datestring     = s_Year+s_Month+s_Day
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
      if t=='':
        s_datetimestring=s_datestring+"0000"
      else:
        s_datetimestring=s_datestring+t
      if os.path.exists(filepath):
        wildcard=filepath+"*"+s_datestring+t+"*"
        fod      = sorted(glob.glob(wildcard))
        if debug>=2:
          print "wildcard is: "+wildcard
      else:
        fod      = -1
        print "ERROR: not existing path: "+filepath
        logs.add(satellite,dt.datetime.strptime(s_datetimestring,"%Y%m%d%H%M"),L15_error=-1,prob_L15="ERROR: not existing path: "+filepath)
        exit()
      try:
        if len(fod)==0:
          print "\nERROR: not any existing file for wildcard: "+wildcard+"\n"
          logs.add(satellite,dt.datetime.strptime(s_datetimestring,"%Y%m%d%H%M"),L15_error=-1,prob_L15="ERROR: not any existing file for wildcard: "+wildcard)
      except TypeError:
        print "\nERROR: not any existing file for wildcard: "+wildcard+"\n"
        logs.add(satellite,dt.datetime.strptime(s_datetimestring,"%Y%m%d%H%M"),L15_error=-1,prob_L15="ERROR: not any existing file for wildcard: "+wildcard)
      staticdir  =  a + versionfolder+'/aux/'
      lutdir     =  a + versionfolder+'/aux/'
      outdir     =  a + versionfolder+'/data/'+satellite+'/'+s_Year+'/'+s_Month+'/'
      if not os.path.exists(outdir):
        try:
          os.makedirs(outdir)
        except OSError:
          print "done in the meantime"

#===============================================================================================
# Part 2 :Check Data
#-----------------------------------------------------------------------------------------------  
      for fname in fod:
        
        #get current date
        if debug>=1:
          print "infile: "+fname
        i_Tindex =  fname.index(s_datestring)+len(s_datestring)
        s_Time   =  fname[i_Tindex:i_Tindex+4]
        
        filename = fname.replace('//','/') #s_datestring+s_Time
        
        curr_datetime=dt.datetime.strptime(s_datestring+s_Time,"%Y%m%d%H%M")
        
        if debug==1:
          print "time   : "+s_Time
        
        #check if file has to be processed:
        exists,info=logs.verify(satellite,curr_datetime)
        if not info:
          logs.add(satellite,curr_datetime,L15_file=fname)
          #err=logs.update(satellite,curr_datetime,"status_easy",-1)
          #err=logs.update(satellite,curr_datetime,"status_full",-1)
#===============================================================================================
# Part 2b :Get Data
#-----------------------------------------------------------------------------------------------  
        if exists==True:
          print exists,info
          print fname
          print "EXIT: data already there or missing because input is missing -> no processing done now"
        if exists==False:
          err=logs.update(satellite,curr_datetime,"ricalpy_ts",curr_datetime)
          err=logs.update(satellite,curr_datetime,"status_easy",0)
          err=logs.update(satellite,curr_datetime,"status_full",0)
          d_pixel_mask=np.zeros((8,N,N),dtype=np.uint8)
          print "Reading data..."
          timing = []
          timing.append(timeit.default_timer())#0
          #fetch calfile+++++++++
          calfile = to.findcal(satellite,irCALpath)
          if debug>=0.1:
            print "calfile: "+calfile
          try:
                        calfile_content=wr.read_calfile(calfile,satellite,doy,s_datestring,s_Time,logs,debug)
          except IOError:
                        calfile_content=-1
                        #report to db
                        err=logs.update(satellite,curr_datetime,"status_easy",-1)
                        err=logs.update(satellite,curr_datetime,"problem_easy","no IR calfile")
                        err=logs.update(satellite,curr_datetime,"status_full",-1)
                        err=logs.update(satellite,curr_datetime,"problem_full","no IR calfile")
          #read rect2lp++++++++++
          timing.append(timeit.default_timer())#1
          with open(fname, "rb") as fi:
              fi, head = rd.header(fi, 1, 3,t=0)
              nomlon = round(np.rad2deg(getattr(head[1],'dsNominalSubSatelliteLongitude')[0]),2)
              lininfo, pixvis, pixir, pixwv = rd.image_3(fi)
              fi, trail = rd.header(fi, 2504, 2504, t=0)
              if not info:
                err=logs.update(satellite,curr_datetime,"input_l15",1)
              else:
                err=logs.update(satellite,curr_datetime,"rect2lp_filename",fname)
          #convert lybia coordinates to pix/line for debugging
          if debug>=1:
            tstsitepixlin=cr.linepixel(5000,tstsite[0],tstsite[1],nomlon)
            tstsitepixIRlin=cr.linepixel(2500,tstsite[0],tstsite[1],nomlon)
            print "\nTestsite LINE/PIXEL: " + str( tstsitepixlin )+"\n"
            tstsitepixlin=list(tstsitepixlin)
            tstsitepixIRlin=list(tstsitepixIRlin)
            tstsitepixlin[0]=tstsitepixlin[0]-1
            tstsitepixlin[1]=tstsitepixlin[1]-1
            tstsitepixIRlin[0]=tstsitepixIRlin[0]-1
            tstsitepixIRlin[1]=tstsitepixIRlin[1]-1
            tstsitepixlin=tuple(tstsitepixlin)
            tstsitepixIRlin=tuple(tstsitepixIRlin)
            print "\nTestsite LINE/PIXEL Idx: " + str( tstsitepixlin )+"\n"

          if not exists:
            #some preps++++++++++++
            sat="".join(getattr(head[0],'szSatCode'))
            satellite=to.satnames(sat)
            tloc=fname.index(s_datestring)
            timestring=fname[tloc:tloc+12]
            year=int(timestring[:4])
            DOY=getattr(head[0],'iYday')[0]
            slot=int(getattr(head[0],'iSlotNum')[0])
            
            #get imag2tg name++++++
            timing.append(timeit.default_timer())#2
            iname=r0.l15_to_L10_name(fname,s_datestring)
            
            #read imag2tg++++++++++
            with open(iname, "rb") as fi:
              try:
                fi,head0   = r0.header(fi)
              except:
                err=logs.update(satellite,curr_datetime,"input_l10",-1)
                exists=True#fcdr does not exist, but because L10 is also missing/corrupt)
              try:
                lininfo0  = r0.telem_mmap2(fi)
                if not info:
                  err=logs.update(satellite,curr_datetime,"input_l10",1)
              except ValueError:
                err=logs.update(satellite,curr_datetime,"input_l10",-1)
                exists=True #fcdr does not exist, but because L10 is also missing/corrupt)
          if not exists:
            #extract telemetry+++++
            timing.append(timeit.default_timer())#3
            try:
              TMerr,telem,telem_descr=r0.decode_telem(lininfo0,head0,mode="rm")
            except IndexError:
              TMerr=-1
            if TMerr!=0:
                err=logs.update(satellite,curr_datetime,"input_l10",0)
                err=logs.update(satellite,curr_datetime,"status_easy",-1)
                err=logs.append(satellite,curr_datetime,"problem_easy","no TM in L1.0 input file")
                err=logs.update(satellite,curr_datetime,"status_full",-1)
                err=logs.append(satellite,curr_datetime,"problem_full","no TM in L1.0")
            #read space corners++++
            timing.append(timeit.default_timer())#4
            with open(iname, "rb") as fi:
              space,corner_means,corner_stds = r0.space_corner_mmap2(fi)
            #prepare SRF+++++++++++
            runfolder=a+versionfolder
            srf=to.makesrf(runfolder,satellite,year,DOY,proc_vers=srf_vers,timestring=timestring,logs=logs,debug=debug)
            if srf.bestsrf==-1:
              err=logs.update(satellite,curr_datetime,"status_easy",-1)
              err=logs.append(satellite,curr_datetime,"problem_easy","no SRF")
              err=logs.update(satellite,curr_datetime,"status_full",-1)
              err=logs.append(satellite,curr_datetime,"problem_full","no SRF")
              if satellite=="MET7" and int(year)>2014:
                #srf=to.makesrf(runfolder,satellite,2014,365,srf_vers,timestring,logs,debug)
                srf=to.makesrf(runfolder,satellite,2014,365,proc_vers=srf_vers,timestring=timestring,logs=logs,debug=debug)
                print "ERROR: no SRF but mitigate using from 2014"
              else:
                print "fatal ERROR: no SRF "
                exit(1)
            srf.srfonfly(0)#file "srf.expect" will be read later
            #reading done++++++++++++++++++++++++++++++++++++++++++++++++++++
            timing.append(timeit.default_timer())#5
            
            #check timing++++++++++
            print "reading took       "+ str((timing[5]-timing[0]))
            print "reading ir-calfile "+ str((timing[1]-timing[0]))
            print "reading L1         "+ str((timing[2]-timing[1]))
            print "reading L0         "+ str((timing[3]-timing[2]))
            print "decoding telem     "+ str((timing[4]-timing[3]))
            print "reading space      "+ str((timing[5]-timing[4]))
    #===============================================================================================
    # Part 3 :calculate contents
    #-----------------------------------------------------------------------------------------------  
            print "Calculating data..."
            timing=[]
            
            #acquisition time++++++
            print "...calc acq time"
            timing.append(timeit.default_timer())#0
            t1,f1,f2,ut=ef.timing(lininfo,slot,v=1)
            if debug>0.1:
              plt.plot(t1[1250,:])
              plt.plot(t1[:,1250])
              plt.show()
            secsince1970,timeoffset=to.time_to_seconds_since(year,i_Month,i_Day,t1)
            if debug>=1:
              print 'timeoffset: ',timeoffset
              print "secsince1970: ",secsince1970[1250,1250]
              print 'time_unc shape: ', np.shape(ut)
              days=(secsince1970[1250,1250]+timeoffset)/86400.
              print dt.datetime(1970, 1, 1,0,0) + dt.timedelta(days)
            
            #space corner evaluation++++++++
            print "...do space corner evaluation"
            timing.append(timeit.default_timer())#1
            #spacecorner from header
            sc_meanS_head=getattr(head[2],'flSpaceCornersVISS')
            sc_meanN_head=getattr(head[2],'flSpaceCornersVISN')
            sc_sdevS_head=getattr(head[2],'flStdSpaceCornersVISS')
            sc_sdevN_head=getattr(head[2],'flStdSpaceCornersVISN')
            adevh,Csh,u_Csh,spaceflagh=ef.space_header_eval(sc_meanS_head,sc_meanN_head,sc_sdevS_head,sc_sdevN_head,logs,satellite,curr_datetime,debug=debug)
            #spacecorner from corners
            adev,Cs,u_Cs,spaceflag=ef.space_eval(space,corner_means,corner_stds,logs,satellite,curr_datetime,debug=debug)
            #compare
            if len(space[space>0])<10000:#too little valid variability; probably space masked in L1.0
              (adev,Cs,u_Cs,spaceflag)=(adevh,Csh,u_Csh,spaceflagh)
            else:#at least some valid space corner pixel
              if u_Cs==0.:#this is suspicious, better use headers
                (adev,Cs,u_Cs,spaceflag)=(adevh,Csh,u_Csh,spaceflagh)
                err=logs.append(satellite,curr_datetime,"problem_easy","no valid spacecorner (masked?)->fallback on value from header; ")
                err=logs.append(satellite,curr_datetime,"problem_full","no valid spacecorner (masked?)->fallback on value from header; ")
              if Csh==0. and Csh>0.:#this is suspicious, better use headers
                (adev,Cs,u_Cs,spaceflag)=(adevh,Csh,u_Csh,spaceflagh)
                err=logs.append(satellite,curr_datetime,"problem_easy","no valid spacecorner (masked?)->fallback on value from header; ")
                err=logs.append(satellite,curr_datetime,"problem_full","no valid spacecorner (masked?)->fallback on value from header; ")
            if debug>=1:
              print "adev="+str(adev)
              print np.mean(space,axis=1)
              print np.std(space,axis=1)
              gr=np.ones(np.shape(space))
              gr[0,:]=gr[0,:]*np.std(space,axis=1)[0]
              gr[1,:]=gr[1,:]*np.std(space,axis=1)[1]
              print np.sqrt(np.sum(gr**2,axis=1))*(1./320000.)
              rz=raw_input("continue?")
            
            #digitalization noise
            digidev=ef.digit_noise(head,inflated)
            
            #SSP+++++++++++++++++++
            print "...calc SSP"
            timing.append(timeit.default_timer())#2
            if debug==1:
              print "".join(getattr(head[0],'szOriginalFormat'))
              print "".join(getattr(head0[0],'szOriginalFormat'))
            #nomlon = np.rad2deg(getattr(head[1],'dsNominalSubSatelliteLongitude')[0])
            sslat, sslon, sslatE, sslonE = to.calcSSP(head,head0,nomlon,early,logs,curr_datetime)
            sslatM                       = (sslat+sslatE)/2
            sslonM                       = (sslon+sslonE)/2
            if debug==1:
              print ' actual SubSatellite Latitude Start: '+str(sslat)
              print ' actual SubSatellite Longitude Start: '+str(sslon)
              print ' actual SubSatellite Latitude End: '+str(sslatE)
              print ' actual SubSatellite Longitude End: '+str(sslonE)
              print ' actual SubSatellite Latitude Mean: '+str(sslatM)
              print ' actual SubSatellite Longitude Mean: '+str(sslonM)
              print ' nominal SubSatellite Longitude: '+str(nomlon)

            #lat, lon, VZA, VAZ and static+++++++++
            print "...read lat/lon from static fcdr"
            timing.append(timeit.default_timer())#3
            staticname=staticdir+wr.ncName(s_Year,s_Month,s_Day,s_Time,"STATIC",nomlon,Software,rel,satellite)
            if not os.path.exists(staticname) and t=='':
              print "   ...first have to prepare static fcdr"
              print "consider to run: ../tools/generate_static_and_lut.py"
            try:
              with nc.Dataset(staticname,"r") as static:
                print "   ...reading: "+staticname
                #slot=int(getattr(head[0],'iSlotNum')[0])
                LAT=static.variables['latitude_vis'][:]
                LON=static.variables['longitude_vis'][:]
                LAT2=static.variables['latitude_ir_wv'][:]
                LON2=static.variables['longitude_ir_wv'][:]
                if debug>=2:
                  to.printimage(LAT,"lat",mini=-90,maxi=90)
                  to.printimage(LON,"lon",mini=-90,maxi=90)
            except:
              print "WARNING: could not read static FCDR file; calculating on the fly"
              print "consider to run: ../tools/generate_static_and_lut.py"
              LAT,LON = cr.latlon(N)
              LAT2,LON2 = cr.latlon(N/2)
              LON     = np.add(LON,nomlon)
              LON2     = np.add(LON2,nomlon)
            #set final mask for visible area
            m1             = np.ones((2500,2500))
            m1[(LAT2<-90)|(LAT2>90)] = 0
            m2             = np.ones((5000,5000))
            m2[ (LAT<-90)|(LAT>90) ]  = 0
            #VZA on tie-point grid
            VZA10,VAZ10 = cr.vza(sslat,sslon,m1[::10,::10],LAT[::10,::10],LON[::10,::10])
            #eventually test static resolution
            if debug==0.5:
              wr.TSTgeomLUT(staticname,nomlon,m1,LAT,LON)
              
            #sensitivities and LUT+++++++++
            gain="".join(getattr(head[2],'szAllChannelGains'))
            if float(rel) < 2.0:
              print "Release "+rel+" --> only calculate cf"
              co=ca.calib(satellite,rel)
              co.datestring=timestring
              try:
                co.prepcf()#[W m^-2 micron^-1]
                dsl=co.daysincelaunch()
                ysl=dsl/365.
                cf=co.cf #[W m^-2 micron^-1]
              except IndexError:
                err=logs.update(satellite,curr_datetime,"input_l15",0)
                err=logs.update(satellite,curr_datetime,"status_easy",-1)
                err=logs.append(satellite,curr_datetime,"problem_easy","before launch - no CF")
                err=logs.update(satellite,curr_datetime,"status_full",-1)
                err=logs.append(satellite,curr_datetime,"problem_full","before launch - no CF")
                raise IndexError
            else:
              print "...read SZA sensitivities from LUT"
              timing.append(timeit.default_timer())#4
              lutname=lutdir+wr.ncName(s_Year,s_Month,s_Day,s_Time,"LUT",nomlon,Software,rel,satellite)
              if (not os.path.exists(lutname) and t==''): #create LUT file only if t is not set
                print "      ...LUT file not there!!"
                print "      run: ../tools/generate_static_and_lut.py"
              try:
                with nc.Dataset(lutname,"r") as lut:
                  print "     ...reading: "+lutname
                  #slot=int(getattr(head[0],'iSlotNum')[0])
                  #s_sza_t1=lut.variables['s_sza_time'][i_Month-1,slot-1,:,:]
                  s_sza_LAT=lut.variables['s_sza_latitude_vis'][i_Month-1,slot-1,:,:]#FIXME: flipping should be done also for sensitivities
                  s_sza_LON=lut.variables['s_sza_longitude_vis'][i_Month-1,slot-1,:,:]#FIXME: flipping should be done also for sensitivities
              except:
                print "WARNING: could not read SZA sensitivities; calculating on the fly"
                print "consider to run: ../tools/generate_static_and_lut.py"
                t=(slot/2.)*(60*60)
                s_sza_LAT=ef.sensi_SZA(t,0,DOY,m1,LAT,0.05,LON,0,debug=debug)
                s_sza_LON=ef.sensi_SZA(t,0,DOY,m1,LAT,0,LON,0.05,debug=debug)
                
              if debug >=1:
                print "\ncentral LAT: \n" + \
                  str( (LAT)[2500,2500] )+' '+str( (LAT)[2500,2499] )+"\n"+ \
                  str( (LAT)[2499,2500] )+' '+str( (LAT)[2499,2499] )+"\n"
                print "\ncentral LON: \n" + \
                  str( (LON)[2500,2500] )+' '+str( (LON)[2500,2499] )+"\n"+ \
                  str( (LON)[2499,2500] )+' '+str( (LON)[2499,2499] )+"\n"
                print "\nTestsite LAT: " + str( LAT[tstsitepixlin] )+"\n"
                print "\nTestsite LON: " + str( LON[tstsitepixlin])+"\n"
                print "\nTestsite s_sza_LAT: " + str( s_sza_LAT[tstsitepixlin] )+"\n"
                print "\nTestsite s_sza_LON: " + str( s_sza_LON[tstsitepixlin] )+"\n"
                auxname="/DSNNAS/Repro/aux/stamp/ancillary/MFG_"+str(int(nomlon)).zfill(3)+"_LatLon.nc"
                print auxname
                with nc.Dataset(auxname,"r") as aux:
                  #to be flipped because in that file origin is left
                  LAT_old=np.fliplr(aux.variables['lat'][:,:])
                  LON_old=np.fliplr(aux.variables['lon'][:,:])
                print "\nold central LAT: \n" + \
                  str( (LAT_old)[1250,1250] )+' '+str( (LAT_old)[1250,1249] )+"\n"+ \
                  str( (LAT_old)[1249,1250] )+' '+str( (LAT_old)[1249,1249] )+"\n"
                print "\nold central LON: \n" + \
                  str( (LON_old)[1250,1250] )+' '+str( (LON_old)[1250,1249] )+"\n"+ \
                  str( (LON_old)[1249,1250] )+' '+str( (LON_old)[1249,1249] )+"\n"
                print "\nIR central LAT: \n" + \
                  str( (LAT2)[1250,1250] )+' '+str( (LAT2)[1250,1249] )+"\n"+ \
                  str( (LAT2)[1249,1250] )+' '+str( (LAT2)[1249,1249] )+"\n"
                print "\nIR central LON: \n" + \
                  str( (LON2)[1250,1250] )+' '+str( (LON2)[1250,1249] )+"\n"+ \
                  str( (LON2)[1249,1250] )+' '+str( (LON2)[1249,1249] )+"\n"
                loc="Kap_Verde"
                suvX=np.array([3400,3800])
                suvY=np.array([3000,3400])
                to.printrecti(np.array(pixvis),LON,LAT,suvX,suvY,timestring+"_"+grid+"_"+loc,save=save)
                #im=plt.imshow(LAT_old-LAT2, clim=(-0.1, 0.1))
                #plt.colorbar(im,fraction=0.046, pad=0.04)
                #plt.show()
                #im=plt.imshow(LON_old-LON2, clim=(-0.1, 0.1))
                #plt.colorbar(im,fraction=0.046, pad=0.04)
                #plt.show()

              print "...calc angles"
              #angles (SZA,SAZ)++++++
              timing.append(timeit.default_timer())#5
              SZA,SAZ=cr.sza(t1,DOY,m1,LAT,LON)
              SAZ10=SAZ[::10,::10]
              del SAZ
              #angles (VZA,VAZ)++++++
              timing.append(timeit.default_timer())#6
              if debug>=2:
                VZA,VAZ=cr.vza(sslat,sslon,m1,LAT,LON)#25.4.2017: now calculated in static FCDR
                to.printimage(VZA,"vza",mini=np.amin(VZA[VZA>-999.]),maxi=np.amax(VZA))
                VZAE,VAZE=cr.vza(sslatE,sslonE,m1,LAT,LON)
                VZD=VZA-VZAE
                VZD[VZA==-999.]=0
                VZD[VZAE==-999.]=0
                to.printimage(VZD,"vza-diff",mini=np.amin(VZD),maxi=np.amax(VZD))
              
            #sun stuff+++++++++++++
            print "...calc sun stuff"
            timing.append(timeit.default_timer())#7
            sundist,sdu,irr,uirr,srf_vis=ir.irrad(int(timestring[:4]),DOY,srffile=srf.expected)#[W m^-2]
            #read IR SRFs++++++++++
            srf_ir=np.loadtxt(srf.expected_ir, skiprows=4)
            srf_wv=np.loadtxt(srf.expected_wv, skiprows=4)
            #convert wavelengths
            srf_wv[:,0]=(1/srf_wv[:,0])*10000
            srf_ir[:,0]=(1/srf_ir[:,0])*10000
            #swap
            srf_wv[:,:]=srf_wv[::-1,:]
            srf_ir[:,:]=srf_ir[::-1,:]
            srfx=[srf_vis,srf_wv,srf_ir]
            if float(rel) >= 2.0:
              #reflectance+++++++++++
              print "...calc reflectance"
              timing.append(timeit.default_timer())#8
              #set calib object
              co=ca.calib(satellite,rel)
              co.datestring=timestring
              try:
                co.prepcf()#[W m^-2 sr^-1/DC]
                print "cf: ",co.cf
                print "gain: ",co.gain
                dsl=co.daysincelaunch()
                ysl=dsl/365.
                cf=co.cf #[W m^-2 sr^-1/DC]
              except IndexError:
                print "logging"
                err=logs.update(satellite,curr_datetime,"status_easy",-1)
                err=logs.append(satellite,curr_datetime,"problem_easy","before launch - no CF")
                err=logs.update(satellite,curr_datetime,"status_full",-1)
                err=logs.append(satellite,curr_datetime,"problem_full","before launch - no CF")
                raise IndexError
              except UnboundLocalError:
                print "logging"
                err=logs.update(satellite,curr_datetime,"status_easy",-1)
                err=logs.append(satellite,curr_datetime,"problem_easy","cf issue...strange gain?")
                err=logs.update(satellite,curr_datetime,"status_full",-1)
                err=logs.append(satellite,curr_datetime,"problem_full","cf issue...strange gain?")
                raise UnboundLocalError
              timing.append(timeit.default_timer())#9
              ref=ef.reflectance(pixvis,cf,Cs,SZA,m2,irr,sundist) #uses numpy
              timing.append(timeit.default_timer())#10
              #ref2=cr.reflectance(pixvis,cf[1],Cs,SZA,m2,irr,sundist)#uses fortran
              # [W m^-2 sr^-1]   /    [W m^-2]*
              #check timing++++++++++
              timing.append(timeit.default_timer())#11
              # contents done+++++++++++++++++++++++++++++++++++++++++++++++++++

              print "acqtime calculation:     " + str(timing[1]-timing[0])
              print "space corner evaluation: " + str(timing[2]-timing[1])
              print "sslat/sslon calculation: " + str(timing[3]-timing[2])
              print "static handling:         " + str(timing[4]-timing[3])
              print "LUT handling:            " + str(timing[5]-timing[4])
              print "sza/saz calculation:     " + str(timing[6]-timing[5])
              print "vza/vaz calculation:     " + str(timing[7]-timing[6])
              print "sun stuff calculation:   " + str(timing[8]-timing[7])
              print "reflectance preparation: " + str(timing[9]-timing[8])
              print "reflectance1 calculation:" + str(timing[10]-timing[9])
              #print "reflectance2 calculation:" + str(timing[11]-timing[10])
              if debug>0:
                to.printimage(ref,"BRF_"+timestring,mini=0,maxi=1,save=save)
                tmp=s_sza_LAT
                print np.median(tmp[m2==1]),np.std(tmp[m2==1]), np.amin(tmp[m2==1]),np.amax(tmp[m2==1])
                to.printimage(abs(s_sza_LAT),"s_sza_LAT",mini=np.median(abs(s_sza_LAT)[m2==1])-(np.std(abs(s_sza_LAT)[m2==1])),maxi=np.median(abs(s_sza_LAT)[m2==1])+(np.std(abs(s_sza_LAT)[m2==1])),save=save)
                to.printimage(abs(s_sza_LON),"s_sza_LON",mini=np.median(abs(s_sza_LON)[m2==1])-(np.std(abs(s_sza_LON)[m2==1])),maxi=np.median(abs(s_sza_LON)[m2==1])+(np.std(abs(s_sza_LON)[m2==1])),save=save)
                to.printimage([pixvis,pixir,pixwv],["VIS","IR","WV"],save=save,direction="h")
              if debug>=2:
                print "further debugging"
                to.printimage(np.array(pixvis),"pixvis")
                radi=cr.radiance(pixvis,cf[1],Cs)
                to.printimage(radi,"radiance",mini=np.amin(radi),maxi=np.amax(radi))
                #to.printimage(radi*np.pi,"radiance x pi",mini=np.amin(radi*np.pi),maxi=np.amax(radi*np.pi))
                print "1"
                to.printimage(t1,"t1",mini=np.amin(t1[t1>0]),maxi=np.amax(t1))
                print "2"
                plt.plot(ut)
                plt.show()
                plt.plot(t1[1250,:][t1[1250,:]>0])
                plt.show()
                #arr2=np.ones((2500,2500))
                #arr=t1[:,1250] #central column
                #arr=np.multiply(arr2,arr).T
                #arr[t1<=0]=t1[t1<=0]
                #tdiff=(arr*3600)-(t1*3600)
                #to.printimage(tdiff,"error of time storing [seconds]",mini=np.amin(tdiff),maxi=np.amax(tdiff))
                to.printimage(LAT,"lat",mini=-90,maxi=90)
                to.printimage(LON,"lon",mini=-90,maxi=90)
                plt.plot(LON[2500,:])#plot East-West
                plt.show()
                print "sza"
                to.printimage(SZA,"sza",mini=np.amin(SZA[SZA>-999.]),maxi=np.amax(SZA))
                plt.plot(SZA[2500,:])
                plt.show()
                print "saz"
                to.printimage(SAZ10,"saz",mini=np.amin(SAZ10[SAZ10>-999.]),maxi=np.amax(SAZ10))
                to.printimage(VAZ10,"vaz",mini=np.amin(VZA10[VZA10>-999.]),maxi=np.amax(VAZ10))
                #to.printimage(np.cos(np.radians(SZA)),"cossza",mini=0,maxi=2)
                #to.printimage(np.cos(np.radians(SZA))*irr,"cossza x irr",mini=0,maxi=irr)

                
            
    #===============================================================================================
    # Part 5 :calculate remaining effects
    #-----------------------------------------------------------------------------------------------  
            print "Calculating remaining effects..."
            
            ##photon shot noise
            #commented out because not significant
            #Pnoise = ef.shot_noise(pixvis,cf[1],Cs,srf[:,0:2],m2,debug=1)
            
            if float(rel) >= 2.0:
              
              timing=[]
              #get lat/lon uncertainty++
              print "...calc lat/lon uncertainty"
              timing.append(timeit.default_timer())#0
              lm_def,geoflag=ef.geolocation(trail,pixir,timestring,sat,debug)
              #print lm_def
              #lm_def: 0=x; 1=y;2=x_ir; 3=y_ir
              ulats=np.sqrt(np.power(ef.dz_easy(LAT,lm_def[0]),2)+np.power(ef.dz_easy(LAT,lm_def[1],axis=1),2))
              ulons=np.sqrt(np.power(ef.dz_easy(LON,lm_def[0]),2)+np.power(ef.dz_easy(LON,lm_def[1],axis=1),2))
              if debug>=1:
                print "\nTestsite u_LAT: " + str( ulats[tstsitepixlin] )+"\n"
                print "\nTestsite u_LON: " + str( ulons[tstsitepixlin] )+"\n"
              if debug>=2:
                tmp=ef.dz_easy(LAT,lm_def[0])
                to.printimage(tmp,"ulats part1",mini=np.amin(tmp[m2==1]),maxi=np.amax(tmp[m2==1]))
                tmp=ef.dz_easy(LAT,lm_def[1],axis=1)
                to.printimage(tmp,"ulats part2",mini=np.amin(tmp[m2==1]),maxi=np.amax(tmp[m2==1]))
              print "...calc sza uncertainty"
              timing.append(timeit.default_timer())#1
              #from lat
              uSZA_lat=s_sza_LAT*ulats
              uSZA_lat[m2==0]=0
              #from lon
              uSZA_lon=s_sza_LON*ulons
              uSZA_lon[m2==0]=0
              #from time
                #neglectable
              #combined
              uSZA=np.sqrt(np.power(uSZA_lat,2)+np.power(uSZA_lon,2))
              del uSZA_lat
              del uSZA_lon
              if debug>=1:
                print "\nTestsite u_SZA: " + str( uSZA[tstsitepixlin] )
                print "Testsite CE: " + str( np.array(pixvis)[tstsitepixlin] )
                print "Testsite CS: " + str( Cs )
                print "Testsite a0: " + str( co.a0 )
                print "Testsite a1: " + str( co.a1 )
                print "Testsite a2: " + str( co.a2 )
                print "Testsite ysl: " + str( ysl )
                print "Testsite SZA: " + str( SZA[tstsitepixlin] )
                print "Testsite UTC: " + str( t1[tstsitepixIRlin] )
                print "Testsite E: " + str( irr)
                print "sundist: " + str( sundist)
                print "Testsite BRF: " + str( ref[tstsitepixlin] )+"\n"
              
              uSZA[m2==0]=0
              print "...calc sensitivities"
              timing.append(timeit.default_timer())#2
              
              if debug>=1:
                for targ in ["Ce","Cs","E","SZA","a0","a1","a2","Z"]:
                  sensitmp=ef.sensi_p(targ,pixvis,Cs,ysl,co,SZA,irr,sundist)
                  sensitmp[sensitmp>1]=1
                  sensitmp[sensitmp<-1]=-1
                  sensitmp[m2==0]=0
                  print "\nTestsite sensi "+targ+": " + str( sensitmp[tstsitepixlin] )+"\n"
                  to.printimage(sensitmp,"sensi"+targ,mini=np.median(sensitmp)-(np.std(sensitmp)),maxi=np.median(sensitmp)+(np.std(sensitmp)),save=save)

              print "...combine random"
              timing.append(timeit.default_timer())#3
              uncMS_rand=np.zeros( (N,N) )
              for targ,uval in zip(["Ce","digit"],[adev,digidev]):
                sensitmp=ef.sensi_p(targ,pixvis,Cs,ysl,co,SZA,irr,sundist)
                unctmp=np.multiply(uval,sensitmp)
                uncMS_rand=uncMS_rand+np.power(unctmp,2)
                if debug>=1:
                  print "\nTestsite uncert "+targ+": " + str( unctmp[tstsitepixlin] )+"\n"
                  if debug>=1.5:
                    cont_arr=SZA
                    cont_arr[cont_arr<0]=np.nan
                    cont_arr[cont_arr>100]=np.nan
                    fill_arr=sensitmp
                    fill_arr[fill_arr>1]=1
                    fill_arr[fill_arr<-1]=-1
                    #to.printcontours(sslatM,sslonM,LAT,LON,cont_arr,fill_arr,m2,"SZA and s_BRF to "+targ+" (absolute)",save=save)
                    #to.printcontours(sslatM,sslonM,LAT,LON,cont_arr,fill_arr/ref,m2,"SZA and s_BRF to "+targ+" (fractional)",save=save)
                    fill_arr=abs(unctmp)
                    fill_arr[fill_arr>1]=1
                    to.printcontours(sslatM,sslonM,LAT,LON,cont_arr,fill_arr,m2,"SZA and u_BRF from "+targ+" (absolute)",maxi=np.nanpercentile(fill_arr,90),save=save)
                    fill_arr=abs(fill_arr/ref)
                    fill_arr[fill_arr>1]=1
                    to.printcontours(sslatM,sslonM,LAT,LON,cont_arr,fill_arr,m2,"SZA and u_BRF from "+targ+" (fractional)",maxi=np.nanpercentile(fill_arr,90),save=save)
                    #plot uncertainty against SZA
                    if debug>=2:
                      plt.plot(np.ndarray.flatten(SZA[m2==1]),np.ndarray.flatten(abs(unctmp[m2==1])),".")
                      plt.title("u_BRF_from_"+targ+" (absolute)")
                      plt.show()
                      plt.plot(np.ndarray.flatten(SZA[m2==1]),np.ndarray.flatten(abs(unctmp[m2==1])/ref[m2==1]),".")
                      plt.title("u_BRF_from_"+targ+" (fractional)")
                      plt.show()
              del sensitmp
              del unctmp
              uncBRF_RMS_rand=np.sqrt(uncMS_rand)
              del uncMS_rand
              uncBRF_RMS_rand[m2==0]=-1
              uncBRF_RMS_rand[uncBRF_RMS_rand>1]=1
              if debug>=1:
                print "\nTestsite random uncert: " + str( uncBRF_RMS_rand[tstsitepixlin] )+"\n"
              if debug>=2:
                to.printimage(uncBRF_RMS_rand,"random uncertainty",mini=np.median(uncBRF_RMS_rand[m2==1])-(np.std(uncBRF_RMS_rand[m2==1])),maxi=np.median(uncBRF_RMS_rand[m2==1])+(np.std(uncBRF_RMS_rand[m2==1])),save=save)

              
              print "...combine structured"
              timing.append(timeit.default_timer())#4
              #corr=np.array([[ 1.  ,  0.99,  0.9 ,  0.  ,  0.  ,  0.  ],\
                             #[ 0.99,  1.  ,  0.9 ,  0.  ,  0.  ,  0.  ],\
                             #[ 0.9 ,  0.9 ,  1.  ,  0.  ,  0.  ,  0.  ],\
                             #[ 0.  ,  0.  ,  0.  ,  1.  ,  0.  ,  0.  ],\
                             #[ 0.  ,  0.  ,  0.  ,  0.  ,  1.  ,  0.  ],\
                             #[ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ,  1.  ]])
              try:
                corrfile="../config/effect_correlation_"+satellite+"_"+str(rel)+".csv"
                corr=np.loadtxt(corrfile,delimiter=";")
              except:
                err=logs.update(satellite,curr_datetime,"status_easy",-1)
                err=logs.update(satellite,curr_datetime,"problem_easy","no correlation file")
                err=logs.update(satellite,curr_datetime,"status_full",-1)
                err=logs.update(satellite,curr_datetime,"problem_full","no correlation file")
                print("fatal ERROR: Error with effect correlations file: "+corrfile)
                raise
                exit(1)
              Es=["E","a0","a1","a2","Z","SZA","Cs"]
              E=len(Es)
              ua=np.diagonal(co.acov)
              try:
                collectu=np.zeros( (N,N,E) )
              except MemoryError:
                err=logs.update(satellite,curr_datetime,"status_easy",-1)
                err=logs.update(satellite,curr_datetime,"problem_easy","MemoryError during collectu def")
                err=logs.update(satellite,curr_datetime,"status_full",-1)
                err=logs.update(satellite,curr_datetime,"problem_full","MemoryError during collectu def")
                raise
              for efN,targ,uval in zip([0,1,2,3,4,5,6],Es,[uirr,np.sqrt(ua[0]),np.sqrt(ua[1]),np.sqrt(ua[2]),co.Z,uSZA,u_Cs]):
                #get sensitivity of current effect
                sensitmp=ef.sensi_p(targ,pixvis,Cs,ysl,co,SZA,irr,sundist)
                #collect uncertainties in 3D array
                collectu[:,:,efN]=np.multiply(sensitmp,uval)
                if debug>=0.5:
                  to.printimage(abs(uval*sensitmp)/ref,"u_BRF_from_"+targ+" (fractional)",mini=0,maxi=0.04,save=save)
              uSZA10=uSZA[::10,::10]
              del uSZA
              del sensitmp
              del uval
              print "   ...apply corr in fortran"
              try:
                uncMS_syst=cr.urrur(collectu,corr)
              except MemoryError:
                err=logs.update(satellite,curr_datetime,"status_easy",-1)
                err=logs.update(satellite,curr_datetime,"problem_easy","MemoryError during cov broadcasting")
                err=logs.update(satellite,curr_datetime,"status_full",-1)
                err=logs.update(satellite,curr_datetime,"problem_full","MemoryError during cov broadcasting")
                raise
              del collectu
              uncBRF_RMS_syst=np.sqrt(uncMS_syst)
              del uncMS_syst
              uncBRF_RMS_syst[m2==0]=-1
              uncBRF_RMS_syst[uncBRF_RMS_syst>1]=1
              
              if debug>=0.5:
                print "plot uncertainties"
                #to.printimage([ref,uncBRF_RMS_rand,uncBRF_RMS_syst],["BRF","random uncertainty","non-random uncertainty"],mini=[0,-0,-0],maxi=[1,0.05,0.05],save=0)
                #plottitels=["BRF","random uncertainty","non-random uncertainty","ratio non-random / random"]
                plottitels=["BRF","Independent uncertainty [absolute]","Structured uncertainty [absolute]","Diff (Structured - Independent)"]
                plotdata  =[ref,uncBRF_RMS_rand,uncBRF_RMS_syst,uncBRF_RMS_syst-uncBRF_RMS_rand]
                to.print4mat(plotdata,plottitels,f2,save=save)
              if debug>=0.1:
                if debug>=1.:
                  print "\nTestsite non-random uncert: " + str( uncBRF_RMS_syst[tstsitepixlin] )+"\n"
                #plottitels=["BRF","fractional random uncertainty","fractional non-random uncertainty","ratio non-random / random"]
                plottitels=["BRF","Independent uncertainty[%]","Structured uncertainty[%]","Diff (Structured - Independent)"]
                uncBRF_RMS_rand_rel=(uncBRF_RMS_rand/ref)*100
                uncBRF_RMS_syst_rel=(uncBRF_RMS_syst/ref)*100
                plotdata  =[ref,uncBRF_RMS_rand_rel,uncBRF_RMS_syst_rel,uncBRF_RMS_syst_rel-uncBRF_RMS_rand_rel]
                to.print4mat(plotdata,plottitels,f2,save=save)
              if debug>=2:
                mi= np.median(uncBRF_RMS_syst[m2==1])-np.std(uncBRF_RMS_syst[m2==1])
                ma= np.median(uncBRF_RMS_syst[m2==1])+np.std(uncBRF_RMS_syst[m2==1])
                #to.printimage(uncBRF_RMS_syst,"uncBRF_RMS_syst",mini=np.median(uncBRF_RMS_syst[m2==1])-(np.std(uncBRF_RMS_syst[m2==1])),maxi=np.median(uncBRF_RMS_syst[m2==1])+(np.std(uncBRF_RMS_syst[m2==1])),save=save)
                to.printimage(uncBRF_RMS_syst,"non-random uncertainty",mini=mi,maxi=ma,save=save)
                #to.printimage(uncBRF_RMS_syst-uncBRF_RMS_rand,"difference non-random minus random uncertainty",mini=mi,maxi=ma,save=save)

              timing.append(timeit.default_timer())#5
              
              print "lat/lon uncert calculation:     " + str(timing[1]-timing[0])
              print "sza     uncert calculation:     " + str(timing[2]-timing[1])
              print "sensitivity    calculation:     " + str(timing[3]-timing[2])
              print "random         calculation:     " + str(timing[4]-timing[3])
              print "systematic     calculation:     " + str(timing[5]-timing[4])
              print "combined       calculation:     commented out"# + str(timing[6]-timing[5])
              
              if debug>=1:
                print np.median(ulats[m2==1]),np.std(ulats[m2==1]), np.amin(ulats[m2==1]),np.amax(ulats[m2==1])
                to.printimage(ulats,"u_lat",mini=np.median(ulats[m2==1])-(0.5*np.std(ulats[m2==1])),maxi=np.median(ulats[m2==1])+(0.5*np.std(ulats[m2==1])),save=save)#np.amax(ulats))
                to.printimage(ulons,"u_lon",mini=np.median(ulons[m2==1])-(0.5*np.std(ulons[m2==1])),maxi=np.median(ulons[m2==1])+(0.5*np.std(ulons[m2==1])),save=save)
                to.printimage(abs(uSZA_lat),"uSZA_lat",mini=np.median(abs(uSZA_lat[m2==1]))-(0.5*np.std(abs(uSZA_lat[m2==1]))),maxi=np.median(abs(uSZA_lat[m2==1]))+(0.5*np.std(abs(uSZA_lat[m2==1]))),save=save)
                to.printimage(abs(uSZA_lon),"uSZA_lon",mini=np.median(abs(uSZA_lon[m2==1]))-(0.5*np.std(abs(uSZA_lon[m2==1]))),maxi=np.median(abs(uSZA_lon[m2==1]))+(0.5*np.std(abs(uSZA_lon[m2==1]))),save=save)
                to.printimage(uSZA,"u_SZA",mini=np.median(uSZA[m2==1])-(0.5*np.std(uSZA[m2==1])),maxi=np.median(uSZA[m2==1])+(0.5*np.std(uSZA[m2==1])),save=save)
              print "...masking"
              #pixel_mask flags++++++
                #q_pixel_mask = quality_pixel_bitmask
                  #automatically set
                #d_pixel_mask = data_quality_bitmask
                  #uncertainty_suspicious
              print "mask SZA"
              d_pixel_mask[0,:,:][SZA>90]=1
              print "mask large u"
              d_pixel_mask[1,:,:][(m2==1)&(pixvis<=Cs)]=1
              print "mask space"
              if spaceflag!=0:
                  d_pixel_mask[2,:,:][m2<2]=1
              print "mask not on earth"
              d_pixel_mask[3,:,:][m2==0]=1
              print "mask suspect_time"
              d_pixel_mask[4,:,:][f2==0]=1
              print "mask suspect_geo"
              if geoflag!=0:
                d_pixel_mask[5,:,:][m2<2]=1
              print "merge"
              ARRB=np.zeros((N,N),dtype=np.uint8)
              for bit in range(8):
                ARRB=ARRB+d_pixel_mask[bit,:,:]*2**(7-bit)
              d_pixel_mask=ARRB
              timing.append(timeit.default_timer())#6
              SZA10=SZA[::10,::10]
              del SZA

    #===============================================================================================
    # Part 6 :write netcdf
    #-----------------------------------------------------------------------------------------------  
            u_srf  =  srf_vis[:, 3:]
            if float(rel) >= 2.0:
              print "Writing easyFCDR data..."

              #filename++++++++++++++
              outeasy   =  outdir+wr.ncName(s_Year,s_Month,s_Day,s_Time,"EASY",nomlon,Software,rel,satellite)
              
              if debug>0:
                print "outfile:"+outeasy
              if os.path.exists(outeasy):
                os.remove(outeasy)
              #doi+++++++++++++++++++
              #in case the metadata has the wrong exact 57 ssp it is there two times
              try:
                doistr={0.0:"10.15770/EUM_SEC_CLM_0009",\
                  57.0:"10.15770/EUM_SEC_CLM_0012",\
                  57.5:"10.15770/EUM_SEC_CLM_0012",\
                  63.0:"10.15770/EUM_SEC_CLM_0013"}[nomlon]
              except KeyError:
                doistr="no doi because of unexpected Sub Satellite Longitude"
              # preparing netcdf file++
              writer = FCDRWriter()
              dataset = writer.createTemplateEasy("MVIRI", 12835,corr_dx=200, corr_dy=200)
              dataset.attrs["institution"] = "EUMETSAT"
              dataset.attrs["title"] = "MVIRI Easy FCDR"
              dataset.attrs["source"] = "Produced from UMARF RECT2LP and IMAG2TG data with MVIRI FCDR code RICalPy, version "+str(Software)
              dataset.attrs["history"] = "Created: " + time.ctime(time.time())+"; added doi"+"; updated IR/WV"
              dataset.attrs["references"] = "in preparation"
              dataset.attrs["comment"] = comment
              dataset.attrs["authors"] = "Frank Ruethrich,Viju John,Rob Roebeling and Joerg Schulz"
              dataset.attrs["email"] = "frank.ruethrich@eumetsat.int"
              dataset.attrs["satellite"] = satellite
              dataset.attrs["url"] = "http://www.fiduceo.eu"
              dataset.attrs["fcdr_software_version"] = str(Software)
              dataset.attrs["data_version"] = str(rel)
              dataset.attrs["RECT2LP_file_name"] = os.path.basename(filename)
              dataset.attrs["channels"]          = "vis, ir, wv"
              dataset.attrs["description"]       = "Meteosat First Generation Rectified (Level 1.5) Image"
              dataset.attrs["doi"]               = doistr
              #preparing relative uncertainties
              uncBRF_RMS_rand_rel=abs(uncBRF_RMS_rand/ref)
              uncBRF_RMS_syst_rel=abs(uncBRF_RMS_syst/ref)
              del uncBRF_RMS_rand
              del uncBRF_RMS_syst
              # populating netcdf file++
              #flag
              dataset.variables["data_quality_bitmask"].data[:,:]                   = (d_pixel_mask)
              #time
              secsince1970[m1==0]=dataset.variables["time"].attrs['_FillValue']
              dataset.variables["time"].data[:,:]                   = (secsince1970)
              dataset.variables["time"].attrs['add_offset']         = timeoffset
              #reflectance
              ref[m2==0]=np.nan#DefaultData.get_default_fill_value(ref.dtype)
              dataset.variables["toa_bidirectional_reflectance_vis"].data[:,:] =(ref)
              del ref
              dataset.variables["toa_bidirectional_reflectance_vis"].attrs["rms_landmarks_x_vis"] = lm_def[0]
              dataset.variables["toa_bidirectional_reflectance_vis"].attrs["rms_landmarks_y_vis"] = lm_def[1]
              #srf
              #dataset.variables["spectral_response_function_vis"].data[:]=srf[:,1]
              #dataset.variables["spectral_response_function_vis"].attrs['version']  = srf_vers
              #dataset.variables["spectral_response_function_vis"].attrs['source']  = os.path.basename(srfsource)
              #dataset.variables["spectral_response_function_vis"].attrs['valid(YYYYDDD)']  = os.path.basename(srfsource)[9:9+7]
              dataset.variables["SRF_frequencies"].data[0,:len(srf_vis[:,0])]= srf_vis[:,0]
              dataset.variables["SRF_frequencies"].data[1,:len(srf_wv[:,0])] = srf_wv[:,0]
              dataset.variables["SRF_frequencies"].data[2,:len(srf_ir[:,0])] = srf_ir[:,0]
              dataset.variables["SRF_weights"].data[0,:len(srf_vis[:,1])]    = srf_vis[:,1]
              dataset.variables["SRF_weights"].data[1,:len(srf_wv[:,1])]     = srf_wv[:,1]
              dataset.variables["SRF_weights"].data[2,:len(srf_ir[:,1])]     = srf_ir[:,1]
              dataset.variables["covariance_spectral_response_function_vis"].data[:np.shape(u_srf)[0],:np.shape(u_srf)[1]] = u_srf
              #independent uncertainties
              uncBRF_RMS_rand_rel[m2==0]  = np.nan#DefaultData.get_default_fill_value(uncBRF_RMS_rand_rel.dtype)
              dataset.variables["u_independent_toa_bidirectional_reflectance"].data[:,:] = (uncBRF_RMS_rand_rel)
              del uncBRF_RMS_rand_rel
              #systematic uncertainties
              uncBRF_RMS_syst_rel[m2==0] = np.nan#DefaultData.get_default_fill_value(uncBRF_RMS_syst_rel.dtype)
              dataset.variables["u_structured_toa_bidirectional_reflectance"].data[:,:]  = (uncBRF_RMS_syst_rel)
              del uncBRF_RMS_syst_rel
              ##other 2D
              pixwv[m1==0]  = DefaultData.get_default_fill_value(pixwv.dtype)
              pixir[m1==0]  = DefaultData.get_default_fill_value(pixir.dtype)
              dataset.variables["count_wv"].data[:,:]               = (np.array(pixwv).astype(int))
              dataset.variables["count_ir"].data[:,:]               = (np.array(pixir).astype(int))
              #lmrks are calculated on IR image; so store them here as attributes
              dataset.variables["count_ir"].attrs["rms_landmarks_x_ir"] = lm_def[2]
              dataset.variables["count_ir"].attrs["rms_landmarks_y_ir"] = lm_def[3]
              dataset.variables["count_wv"].attrs["rms_landmarks_x_wv"] = lm_def[2]
              dataset.variables["count_wv"].attrs["rms_landmarks_y_wv"] = lm_def[3]       
              #tie-point grids
              SAZ10[f2[::10,::10]==0]=np.nan#DefaultData.get_default_fill_value(SAZ.dtype)
              dataset.variables["solar_azimuth_angle"].data[:,:]    = (SAZ10)
              SZA10[f2[::10,::10]==0]=np.nan#DefaultData.get_default_fill_value(SZA.dtype)
              dataset.variables["solar_zenith_angle"].data[:,:]     = (SZA10)
              VAZ10[f2[::10,::10]==0]=np.nan#DefaultData.get_default_fill_value(SAZ.dtype)
              dataset.variables["satellite_azimuth_angle"].data[:,:]    = (VAZ10)
              VZA10[f2[::10,::10]==0]=np.nan#DefaultData.get_default_fill_value(SZA.dtype)
              dataset.variables["satellite_zenith_angle"].data[:,:]     = (VZA10)
              #scalars
              dataset.variables["years_since_launch"].data   = ysl
              dataset.variables["solar_irradiance_vis"].data   = irr
              dataset.variables["u_solar_irradiance_vis"].data = uirr
              dataset.variables["distance_sun_earth"].data     = sundist
              dataset.variables["sub_satellite_longitude_start"].data = sslat
              dataset.variables["sub_satellite_latitude_start"].data  = sslon
              dataset.variables["sub_satellite_longitude_end"].data   = sslatE
              dataset.variables["sub_satellite_latitude_end"].data    = sslonE
              try:
                wr.write_calfile2nc( dataset, calfile_content[0], calfile_content[1],\
                            calfile_content[2], calfile_content[3])
              except TypeError:
                print "WARNING: something wrong with IR calfile"
                err=logs.update(satellite,curr_datetime,"status_easy",-1)
                err=logs.update(satellite,curr_datetime,"problem_easy","something wrong with IR calfile")
              dataset.variables["bt_a_ir"].data = to.BTstuff(satellite)[0]
              dataset.variables["bt_a_wv"].data = to.BTstuff(satellite)[1]
              dataset.variables["bt_b_ir"].data = to.BTstuff(satellite)[2]
              dataset.variables["bt_b_wv"].data = to.BTstuff(satellite)[3]
              #writing netcdf file+++++
              try:
                writer.write(dataset, outeasy)
              except RuntimeError:
                err=logs.update(satellite,curr_datetime,"status_easy",-1)
                err=logs.update(satellite,curr_datetime,"problem_easy","RuntimeError: NetCDF: HDF error")
                print("fatal ERROR : RuntimeError: NetCDF: HDF error - easyFCDR datetime: "+curr_datetime)
                raise
                exit(1)
              #update in db
              nam,status=logs.query(satellite,curr_datetime,"status_easy")
              print status[0]
              if status[0]==0:
                err=logs.update(satellite,curr_datetime,"status_easy",1)

            print "Writing regular/full FCDR data..."
            #full fcdr+++++++++++++++++
            
            #filename++++++++++++++++++
            if float(rel) < 2.0:
              outnc   =  outdir+wr.ncName(s_Year,s_Month,s_Day,s_Time,"",nomlon,Software,rel,satellite)
            else:
              outnc   =  outdir+wr.ncName(s_Year,s_Month,s_Day,s_Time,"FULL",nomlon,Software,rel,satellite)
            if debug==1:
              print "outfile:"+outnc
              
            #prepare metadata++++++
            head1 = head[0]._asdict()
            head2 = head[1]._asdict()
            head3 = head[2]._asdict()
            trail = trail[0]._asdict()
            
            #prepare file
            if float(rel) < 2.0:
              print "dropped support for Release 1!!"
              exit(1)
              
            if debug>0:
              print "outfile:"+outnc
              
            if os.path.exists(outnc):
              os.remove(outnc)
            
            writer = FCDRWriter()
            mfggrp = writer.createTemplateFull("MVIRI", 12835)
            
            # adding global attributes
            mfggrp.attrs["institution"] = "EUMETSAT"
            mfggrp.attrs["title"] = "MVIRI Full FCDR"
            mfggrp.attrs["source"] = "Produced from UMARF RECT2LP and IMAG2TG data with MVIRI FCDR code RICalPy, version "+str(Software)
            mfggrp.attrs["history"] = "Created: " + time.ctime(time.time())+"; added doi"+"; updated IR/WV"
            mfggrp.attrs["references"] = "in preparation"
            mfggrp.attrs["comment"] = comment
            mfggrp.attrs["authors"] = "Frank Ruethrich, Viju John, Rob Roebeling and Joerg Schulz"
            mfggrp.attrs["email"] = "frank.ruethrich@eumetsat.int"
            mfggrp.attrs["satellite"] = satellite
            mfggrp.attrs["url"] = "http://www.fiduceo.eu"
            mfggrp.attrs["fcdr_software_version"] = str(Software)
            mfggrp.attrs["data_version"] = str(rel)
            mfggrp.attrs["RECT2LP_file_name"] = os.path.basename(filename)
            mfggrp.attrs["channels"]          = "vis, ir, wv"
            mfggrp.attrs["description"]       = "Meteosat First Generation Rectified (Level 1.5) Image"
            mfggrp.attrs["doi"]               = doistr
            #adding srf
            mfggrp.variables["SRF_frequencies"].data[0,:len(srfx[0][:,0])] = srfx[0][:,0]
            mfggrp.variables["SRF_frequencies"].data[1,:len(srfx[1][:,0])] = srfx[1][:,0]
            mfggrp.variables["SRF_frequencies"].data[2,:len(srfx[2][:,0])] = srfx[2][:,0]
            mfggrp.variables["SRF_weights"].data[0,:len(srfx[0][:,1])]     = srfx[0][:,1]
            mfggrp.variables["SRF_weights"].data[1,:len(srfx[1][:,1])]     = srfx[1][:,1]
            mfggrp.variables["SRF_weights"].data[2,:len(srfx[2][:,1])]     = srfx[2][:,1]
            mfggrp.variables["covariance_spectral_response_function_vis"].data[:np.shape(u_srf)[0],:np.shape(u_srf)[1]] = u_srf
              
            #mask
            mfggrp.variables["data_quality_bitmask"].data[:,:] = (d_pixel_mask)
            
            #time
            #secsince1970[m1==0]=mfggrp.variables["time"].attrs['_FillValue']
            mfggrp.variables["time"].data[:,:]                   = (secsince1970)
            mfggrp.variables["time"].attrs['add_offset']         = timeoffset
              
            #images
            pixvis[m2==0]  = mfggrp.variables["count_vis"].attrs['_FillValue']#DefaultData.get_default_fill_value(pixvis.astype(int).dtype)
            #pixwv[m1==0]  = DefaultData.get_default_fill_value(pixwv.dtype)
            #pixir[m1==0]  = DefaultData.get_default_fill_value(pixir.dtype)
            mfggrp.variables["count_wv"].data[:,:]               = (np.array(pixwv).astype(int))
            mfggrp.variables["count_ir"].data[:,:]               = (np.array(pixir).astype(int))
            mfggrp.variables["count_vis"].data[:,:]              = (np.array(pixvis).astype(int))
            
            #tie-point grid
            ulats[f2==0]=np.nan
            mfggrp.variables["u_latitude"].data[:,:]    = (ulats[::10,::10])
            ulons[f2==0]=np.nan
            mfggrp.variables["u_longitude"].data[:,:]   = (ulons[::10,::10])
            SAZ10[f2[::10,::10]==0]=np.nan#DefaultData.get_default_fill_value(SAZ.dtype)
            mfggrp.variables["solar_azimuth_angle"].data[:,:]    = (SAZ10)
            SZA10[f2[::10,::10]==0]=np.nan
            mfggrp.variables["solar_zenith_angle"].data[:,:]     = (SZA10)
            uSZA10[f2[::10,::10]==0]=np.nan
            mfggrp.variables["u_solar_zenith_angle"].data[:,:]     = (uSZA10)
            VAZ10[f2[::10,::10]==0]=np.nan#DefaultData.get_default_fill_value(SAZ.dtype)
            mfggrp.variables["satellite_azimuth_angle"].data[:,:]    = (VAZ10)
            VZA10[f2[::10,::10]==0]=np.nan#DefaultData.get_default_fill_value(SZA.dtype)
            mfggrp.variables["satellite_zenith_angle"].data[:,:]     = (VZA10)
            
            #small 2D variables
            mfggrp.variables["covariance_a_vis"].data   = co.acov
            mfggrp.variables["effect_correlation_matrix"].data     = corr
            
            # 1D arrays
            mfggrp.variables["u_time"].data[:]   = ut
            
            #scalars
            mfggrp.variables["years_since_launch"].data   = ysl
            mfggrp.variables["a0_vis"].data   = co.a0
            mfggrp.variables["a1_vis"].data   = co.a1
            mfggrp.variables["a2_vis"].data   = co.a2
            mfggrp.variables["u_a0_vis"].data   = np.sqrt(co.acov[0,0])
            mfggrp.variables["u_a1_vis"].data   = np.sqrt(co.acov[1,1])
            mfggrp.variables["u_a2_vis"].data   = np.sqrt(co.acov[2,2])
            mfggrp.variables["u_zero_vis"].data   = co.Z
            print adev
            try:
              mfggrp.variables["allan_deviation_counts_space_vis"].data   = adev[0]
              mfggrp.variables["u_electronics_counts_vis"].data   = adev[0]
            except IndexError:
              mfggrp.variables["allan_deviation_counts_space_vis"].data   = float(adev)
              mfggrp.variables["u_electronics_counts_vis"].data   = float(adev)
            except TypeError:
              mfggrp.variables["allan_deviation_counts_space_vis"].data   = float(adev)
              mfggrp.variables["u_electronics_counts_vis"].data   = float(adev)
            mfggrp.variables["solar_irradiance_vis"].data   = irr
            mfggrp.variables["u_solar_irradiance_vis"].data = uirr
            mfggrp.variables["distance_sun_earth"].data     = sundist
            print Cs
            mfggrp.variables["mean_count_space_vis"].data   = float(Cs)
            mfggrp.variables["u_digitization_counts_vis"].data   = digidev
            print u_Cs
            mfggrp.variables["u_mean_counts_space_vis"].data   = float(u_Cs)
            #ir-calfile
            try:
              wr.write_calfile2nc( mfggrp, calfile_content[0], calfile_content[1],\
                                      calfile_content[2], calfile_content[3])
            except TypeError:
              print "WARNING: something wrong with IR calfile"
              err=logs.update(satellite,curr_datetime,"status_full",-1)
              err=logs.update(satellite,curr_datetime,"problem_full","something wrong with IR calfile")
            mfggrp.variables["bt_a_ir"].data = to.BTstuff(satellite)[0]
            mfggrp.variables["bt_a_wv"].data = to.BTstuff(satellite)[1]
            mfggrp.variables["bt_b_ir"].data = to.BTstuff(satellite)[2]
            mfggrp.variables["bt_b_wv"].data = to.BTstuff(satellite)[3]
            
            #historic data
            wr.write_manager( mfggrp, filename, head1, head2, head3, lininfo, trail, telem, telem_descr, \
                              satellite,logs)
            
            #write
            try:
              err=writer.write( mfggrp, outnc )
            except RuntimeError:
              err=logs.update(satellite,curr_datetime,"status_full",-1)
              err=logs.update(satellite,curr_datetime,"problem_full","RuntimeError: NetCDF: HDF error")
              print("fatal ERROR: RuntimeError: NetCDF: HDF error - fullFCDR datetime: "+curr_datetime)
              raise
              exit(1)
            nam,status=logs.query(satellite,curr_datetime,"status_full")
            print status[0]
            if status[0]==0:
              print "success"
              err=logs.update(satellite,curr_datetime,"status_full",1)
            timing.append(timeit.default_timer())#7
            
            print "masking took: "+ str(timing[6]-timing[5])
            print "writing took: "+str(timing[7]-timing[6])
            
if __name__ == "__main__":
   main(sys.argv[1:])
