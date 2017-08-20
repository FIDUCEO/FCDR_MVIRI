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

import ricals_netcdf as wr
import ricals_effects as ef
import ricals_tools as to
import ricals_calib as ca
import ricals_irrad as ir
import read_rect2lp as rd
import read_imag2tg as r0

from fiduceo.fcdr.writer.fcdr_writer import FCDRWriter
from fiduceo.fcdr.writer.default_data import DefaultData

sys.path.insert(0, os.getcwd()[:-4]+'/lib/nrCrunch/')

#import cruncher_f as cr #grid definition as implemented in fortran routines in mfg user manual
#global grid
#grid="mfg_imp"

import cruncher as cr #grid definition as described in mfg user manual
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
Software = "2.2"
N        = 5000
early    = ['M2','M3','P2']
parallel = "n"
save     = 1
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
        debug = float(str(arg))
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
        print "ERROR: not existing path: "+filepath
        exit()
      try:
        if len(fod)==0:
          print "\nERROR: not existing file for wildcard: "+wildcard+"\n"
          exit()
      except TypeError:
        print "\nERROR: not existing file for wildcard: "+wildcard+"\n"
        exit()
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
        
        pixel_mask=np.zeros((N,N)).astype(np.uint8)
        
#===============================================================================================
# Part 2 :Get Data
#-----------------------------------------------------------------------------------------------  
        print "Reading data..."
        timing = []
        timing.append(timeit.default_timer())#0
        #fetch calfile+++++++++
        calfile = to.findcal(satellite,s_datestring,s_Time,irCALpath)
        if debug==1:
          print "calfile: "+calfile
        calfile_content=wr.read_calfile(calfile,satellite,doy,s_datestring,s_Time,debug)
        #read rect2lp++++++++++
        timing.append(timeit.default_timer())#1
        with open(fname, "rb") as fi:
            fi, head = rd.header(fi, 1, 3,t=0)
            nomlon = np.rad2deg(getattr(head[1],'dsNominalSubSatelliteLongitude')[0])
            #testname=outdir+'FCDR_'+sensor+'_'+satellite+'_'+str(nomlon).zfill(4)+'_'+s_Year+s_Month+s_Day+s_Time+'_full'+'.nc'
            if float(rel) < 2.0:
              testname   =  outdir+wr.ncName(s_Year,s_Month,s_Day,s_Time,"",nomlon,Software,rel,satellite)
            else:
              testname   =  outdir+wr.ncName(s_Year,s_Month,s_Day,s_Time,"FULL",nomlon,Software,rel,satellite)
            if os.path.exists(testname):
              exists=True
              print "\n\n**WARNING: outfile already existing...skipping\n\n"
            if not exists:
              fi, lininfo, pixvis, pixir, pixwv = rd.image(fi,head,4,2503, t=0)
              fi, trail = rd.header(fi, 2504, 2504, t=0)
              
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
          iname=r0.l15_to_L10_name(fname,s_datestring,satellite)
          
          #read imag2tg++++++++++
          with open(iname, "rb") as fi:
            fi,head0   = r0.header(fi)
            lininfo0  = r0.telem_mmap2(fi)
          
          #extract telemetry+++++
          timing.append(timeit.default_timer())#3
          telem,telem_descr=r0.decode_telem(lininfo0,head0,mode="rm")
          
          #read space corners++++
          timing.append(timeit.default_timer())#4
          with open(iname, "rb") as fi:
            space = r0.space_corner_mmap2(fi)
          
          #prepare SRF+++++++++++
          
          srfname=a+versionfolder+"/aux/sensor_data/"+sat+timestring+"_srf.dat"
          if srf_vers in ["pre-launch","pre_launch","pre-flight","pre_flight"]:
            srfsource="../sensor_data/pre_launch/spectral_response_"+satellite[:3]+"_"+satellite[3:]+"_VIS.dat"
          else:
            srfsource=to.bestsrfname("/DSNNAS/Repro/mviri/level1/MFG15_FCDR_V1/tmp/",satellite,year,DOY,srf_vers)
            if srfsource==-1:
              print "ERROR: no recovered SRF (testloc=1)"
              return
          if not os.path.exists(a+versionfolder+"/aux/sensor_data/"):
            os.makedirs(a+versionfolder+"/aux/sensor_data/")
          to.srfonfly(srfsource,srfname,0,debug=debug)#file "srfname" will be read later
          
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
          t1,f1,f2,ut=ef.timing(lininfo,slot)
          secsince1970,timeoffset=to.time_to_seconds_since(year,i_Month,i_Day,t1)
          if debug>=1:
            print 'timeoffset: ',timeoffset
            print "secsince1970: ",secsince1970[1250,1250]
            print 'time_unc shape: ', np.shape(ut)
            days=(secsince1970[1250,1250]+timeoffset)/86400.
            print dt.datetime(1970, 1, 1,0,0) + dt.timedelta(days)
          
          #put f2 as flag into pixel_mask
          pixel_mask[f2==0]=map(functools.partial(to.bitmod,X=1), pixel_mask[f2==0])
          
          #allandeviation++++++++
          print "...calc allan dev"
          timing.append(timeit.default_timer())#1
          adev,Cs=ef.allan_noise(space)
          if debug>=1:
            print "adev="+str(adev)
          
          #digitalization noise
          digidev=ef.digit_noise(head,inflated)
          
          #SSP+++++++++++++++++++
          print "...calc SSP"
          timing.append(timeit.default_timer())#2
          if debug==1:
            print "".join(getattr(head[0],'szOriginalFormat'))
            print "".join(getattr(head0[0],'szOriginalFormat'))
          #nomlon = np.rad2deg(getattr(head[1],'dsNominalSubSatelliteLongitude')[0])
          sslat, sslon, sslatE, sslonE = to.calcSSP(head,head0,nomlon,early)
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
            LAT,LON   = cr.latlon(N)
            LON       = np.add(LON,nomlon)
            LAT2,LON2 = cr.latlon(N/2)
            LON2      = np.add(LON2,nomlon)
            m1                       =  np.ones((2500,2500))
            m1[(LAT2<-90)|(LAT2>90)] =  0
            m2                       =  np.ones((5000,5000))
            m2[ (LAT<-90)|(LAT>90) ] =  0
            print "      ...writing: "+staticname
            if not os.path.exists(staticdir):
              os.makedirs(staticdir)
            if float(rel)<2.0:
              varnames=["latitude_vis","longitude_vis","latitude_ir_wv","longitude_ir_wv"]
            else:
              varnames=["latitude_vis","longitude_vis","latitude_ir_wv","longitude_ir_wv","viewing_geometry"]
            vals=[LAT,LON,LAT2,LON2,0]
            wr.save_static(staticname,nomlon,varnames,vals,m1,comment,rel,debug=debug)
            del(LAT)
            del(LON)
            del(LAT2)
            del(LON2)
          try:
            with nc.Dataset(staticname,"r") as static:
              print "   ...reading: "+staticname
              #slot=int(getattr(head[0],'iSlotNum')[0])
              LAT=np.fliplr(static.variables['latitude_vis'][:])#FIXME: flipping should be done also for sensitivities
              LON=np.fliplr(static.variables['longitude_vis'][:])#FIXME: flipping should be done also for sensitivities
              LAT2=np.fliplr(static.variables['latitude_ir_wv'][:])
              if debug>=1:
                LON2=np.fliplr(static.variables['longitude_ir_wv'][:])
              #VZA not needed for processing
              if debug>=2:
                to.printimage(LAT,"lat",mini=-90,maxi=90)
                to.printimage(LON,"lon",mini=-90,maxi=90)
          except:
            print "WARNING: could not read static FCDR file; calculating on the fly"
            LAT,LON = cr.latlon(N)
            LAT2,LON2 = cr.latlon(N/2)
            LON     = np.add(LON,nomlon)
            LON2     = np.add(LON2,nomlon)
            #VZA,VAZ = cr.vza(sslat,sslon,m1,LAT,LON)#VZA not needed for processin
          #set final mask for visible area
          m1             = np.ones((2500,2500))
          m1[(LAT2<-90)|(LAT2>90)] = 0
          m2             = np.ones((5000,5000))
          m2[ (LAT<-90)|(LAT>90) ]  = 0
          #eventually test static resolution
          if debug==0.5:
            wr.TSTgeomLUT(staticname,nomlon,m1,LAT,LON)
          #sensitivities and LUT+++++++++
          gain="".join(getattr(head[2],'szAllChannelGains'))
          if float(rel) < 2.0:
            print "Release "+rel+" --> only calculate cf"
            dsl=ca.daysincelaunch(sat,timestring)
            ysl=dsl/365.
            cf=ca.cf(ysl,timestring,sat,gain,rel) #[W m^-2 micron^-1] #content cf: 0=cf-uncert(t);1=cf(t);2=a;3=cov;
          else:
            print "...read SZA sensitivities from LUT"
            timing.append(timeit.default_timer())#4
            lutname=lutdir+wr.ncName(s_Year,s_Month,s_Day,s_Time,"LUT",nomlon,Software,rel,satellite)
            if (not os.path.exists(lutname) and t==''): #create LUT file only if t is not set
              print "      ...first have to prepare LUT file"
              print "      ...writing: "+lutname
              if not os.path.exists(staticdir):
                os.makedirs(lutdir)
              varnames=['s_sza_latitude_vis','s_sza_longitude_vis','s_sza_time']
              vals=[LAT,LON,0]
              wr.save_static(lutname,nomlon,varnames,vals,m1,comment,rel,debug=debug)
            try:
              with nc.Dataset(lutname,"r") as lut:
                print "     ...reading: "+lutname
                #slot=int(getattr(head[0],'iSlotNum')[0])
                #s_sza_t1=lut.variables['s_sza_time'][i_Month-1,slot-1,:,:]
                s_sza_LAT=lut.variables['s_sza_latitude_vis'][i_Month-1,slot-1,:,:]#FIXME: flipping should be done also for sensitivities
                s_sza_LON=lut.variables['s_sza_longitude_vis'][i_Month-1,slot-1,:,:]#FIXME: flipping should be done also for sensitivities
            except:
              print "WARNING: could not read SZA sensitivities; calculating on the fly"
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
          sundist,sdu,irr,uirr,srf=ir.irrad(int(timestring[:4]),DOY,srffile=srfname)#[W m^-2]
          

          if float(rel) >= 2.0:
            #reflectance+++++++++++
            print "...calc reflectance"
            timing.append(timeit.default_timer())#8
            dsl=ca.daysincelaunch(sat,timestring)
            ysl=dsl/365.
            cf=ca.cf(ysl,timestring,sat,gain,rel) #[W m^-2 micron^-1] #content cf: 0=cf-uncert(t);1=cf(t);2=a;3=cov;
            timing.append(timeit.default_timer())#9
            ref=ef.reflectance(pixvis,cf[1],Cs,SZA,m2,irr,sundist) #uses numpy
            timing.append(timeit.default_timer())#10
            #ref2=cr.reflectance(pixvis,cf[1],Cs,SZA,m2,irr,sundist)#uses fortran
            # [W m^-2 micron^-1]   /    [W m^-2]*
            #check timing++++++++++
            timing.append(timeit.default_timer())#11
            # contents done+++++++++++++++++++++++++++++++++++++++++++++++++++

            print "acqtime calculation:     " + str(timing[1]-timing[0])
            print "allandev calculation:    " + str(timing[2]-timing[1])
            print "sslat/sslon calculation: " + str(timing[3]-timing[2])
            print "static handling:         " + str(timing[4]-timing[3])
            print "LUT handling:            " + str(timing[5]-timing[4])
            print "sza/saz calculation:     " + str(timing[6]-timing[5])
            print "vza/vaz calculation:     " + str(timing[7]-timing[6])
            print "sun stuff calculation:   " + str(timing[8]-timing[7])
            print "reflectance preparation: " + str(timing[9]-timing[8])
            print "reflectance1 calculation:" + str(timing[10]-timing[9])
            #print "reflectance2 calculation:" + str(timing[11]-timing[10])
            if debug>=1:
              to.printimage(ref,"BRF_"+timestring,mini=0,maxi=1,save=save)
              tmp=s_sza_LAT
              print np.median(tmp[m2==1]),np.std(tmp[m2==1]), np.amin(tmp[m2==1]),np.amax(tmp[m2==1])
              to.printimage(abs(s_sza_LAT),"s_sza_LAT",mini=np.median(abs(s_sza_LAT)[m2==1])-(np.std(abs(s_sza_LAT)[m2==1])),maxi=np.median(abs(s_sza_LAT)[m2==1])+(np.std(abs(s_sza_LAT)[m2==1])),save=save)
              to.printimage(abs(s_sza_LON),"s_sza_LON",mini=np.median(abs(s_sza_LON)[m2==1])-(np.std(abs(s_sza_LON)[m2==1])),maxi=np.median(abs(s_sza_LON)[m2==1])+(np.std(abs(s_sza_LON)[m2==1])),save=save)
              to.printimage([pixvis,pixir,pixwv],["VIS","IR","WV"],save=save,direction="h")
            if debug>=2:
              print "further debugging"
              to.printimage(np.array(pixvis),"pixvis")
              to.printimage(t1,"t1",mini=np.amin(t1[t1>0]),maxi=np.amax(t1))
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
              to.printimage(SAZ,"saz",mini=np.amin(SAZ[SAZ>-999.]),maxi=np.amax(SAZ))
              to.printimage(VAZ,"vaz",mini=np.amin(VZA[VZA>-999.]),maxi=np.amax(VAZ))
              #to.printimage(np.cos(np.radians(SZA)),"cossza",mini=0,maxi=2)
              #to.printimage(np.cos(np.radians(SZA))*irr,"cossza x irr",mini=0,maxi=irr)
              radi=cr.radiance(pixvis,cf[1],Cs)
              to.printimage(radi,"radiance",mini=np.amin(radi),maxi=np.amax(radi))
              #to.printimage(radi*np.pi,"radiance x pi",mini=np.amin(radi*np.pi),maxi=np.amax(radi*np.pi))
              
          
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
            lm_def=ef.geolocation(trail,pixir,timestring,sat,debug)
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
            if debug>=1:
              print "\nTestsite u_SZA: " + str( uSZA[tstsitepixlin] )
              print "Testsite CE: " + str( np.array(pixvis)[tstsitepixlin] )
              print "Testsite CS: " + str( Cs )
              print "Testsite a0: " + str( cf[2][0] )
              print "Testsite a1: " + str( cf[2][1] )
              print "Testsite ysl: " + str( ysl )
              print "Testsite SZA: " + str( SZA[tstsitepixlin] )
              print "Testsite SAZ: " + str( SAZ[tstsitepixlin] )
              print "Testsite UTC: " + str( t1[tstsitepixIRlin] )
              print "Testsite E: " + str( irr)
              print "sundist: " + str( sundist)
              print "Testsite BRF: " + str( ref[tstsitepixlin] )+"\n"
            
            uSZA[m2==0]=0
            print "...calc sensitivities"
            timing.append(timeit.default_timer())#2
            
            if debug>=1:
              for targ in ["Ce","Cs","E","SZA","a0","a1","Z"]:
                sensitmp=ef.sensi_p(targ,pixvis,Cs,ysl,cf[2][0],cf[2][1],SZA,irr,sundist)
                sensitmp[sensitmp>1]=1
                sensitmp[sensitmp<-1]=-1
                sensitmp[m2==0]=0
                print "\nTestsite sensi "+targ+": " + str( sensitmp[tstsitepixlin] )+"\n"
                to.printimage(sensitmp,"sensi"+targ,mini=np.median(sensitmp)-(np.std(sensitmp)),maxi=np.median(sensitmp)+(np.std(sensitmp)),save=save)

            print "...combine random"
            timing.append(timeit.default_timer())#3
            uncMS_rand=np.zeros( (N,N) )
            for targ,uval in zip(["Ce","Cs"],[adev,adev]):
              sensitmp=ef.sensi_p(targ,pixvis,Cs,ysl,cf[2][0],cf[2][1],SZA,irr,sundist)
              unctmp=uval*sensitmp
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
                  to.printcontours(sslatM,sslonM,LAT,LON,cont_arr,fill_arr,m2,"SZA and s_BRF to "+targ+" (absolute)",save=save)
                  to.printcontours(sslatM,sslonM,LAT,LON,cont_arr,fill_arr/ref,m2,"SZA and s_BRF to "+targ+" (fractional)",save=save)
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
            uncBRF_RMS_rand=np.sqrt(uncMS_rand)
            uncBRF_RMS_rand[m2==0]=-1
            uncBRF_RMS_rand[uncBRF_RMS_rand>1]=1
            if debug>=1:
              print "\nTestsite random uncert: " + str( uncBRF_RMS_rand[tstsitepixlin] )+"\n"
            if debug>=2:
              to.printimage(uncBRF_RMS_rand,"random uncertainty",mini=np.median(uncBRF_RMS_rand[m2==1])-(np.std(uncBRF_RMS_rand[m2==1])),maxi=np.median(uncBRF_RMS_rand[m2==1])+(np.std(uncBRF_RMS_rand[m2==1])),save=save)

            
            print "...combine structured"
            #E.W.:I think the aim is to decide if a particular error 
            #structure is "mostly random" or "mostly systematic" for 
            #scales over which a user is most interested, and then 
            #to combine (simply by adding in quadrature their effect 
            #on the main variable - so sensitivity coefficient times 
            #uncertainty) all the ones you think are "mostly 
            #systematic" and all the ones you think are "mostly random".
            timing.append(timeit.default_timer())#4
            uncMS_syst=np.zeros((N,N))
            for targ,uval in zip(["E","a0","a1","Z","SZA"],[uirr,np.sqrt(cf[3][0,0]),np.sqrt(cf[3][1,1]),0,uSZA]):
              sensitmp=ef.sensi_p(targ,pixvis,Cs,ysl,cf[2][0],cf[2][1],SZA,irr,sundist)
              unctmp=uval*sensitmp
              uncMS_syst=uncMS_syst+np.power(unctmp,2)
              if debug>=1:
                print "\nTestsite uncert "+targ+": " + str( unctmp[tstsitepixlin] )+"\n"
                if debug>=1.5:
                  cont_arr=SZA
                  cont_arr[cont_arr<0]=np.nan
                  cont_arr[cont_arr>100]=np.nan
                  fill_arr=sensitmp
                  fill_arr[fill_arr>1]=1
                  fill_arr[fill_arr<-1]=-1
                  to.printcontours(sslatM,sslonM,LAT,LON,cont_arr,fill_arr,m2,"SZA and s_BRF to "+targ+" (absolute)",save=save)
                  to.printcontours(sslatM,sslonM,LAT,LON,cont_arr,fill_arr/ref,m2,"SZA and s_BRF to "+targ+" (fractional)",save=save)
                  fill_arr=abs(unctmp)
                  fill_arr[fill_arr>1]=1
                  to.printcontours(sslatM,sslonM,LAT,LON,cont_arr,fill_arr,m2,"SZA and u_BRF from "+targ+" (absolute)",maxi=0.1,save=save)
                  fill_arr=abs(fill_arr/ref)
                  fill_arr[fill_arr>1]=1
                  to.printcontours(sslatM,sslonM,LAT,LON,cont_arr,fill_arr,m2,"SZA and u_BRF from "+targ+" (fractional)",maxi=0.1,save=save)
              if debug>=2:
                to.printimage(np.sqrt(uncMS_syst),"unc"+targ,mini=0,maxi=0.3,save=0)
            uncBRF_RMS_syst=np.sqrt(uncMS_syst)
            uncBRF_RMS_syst[m2==0]=-1
            uncBRF_RMS_syst[uncBRF_RMS_syst>1]=1
            if debug>=0.5:
             print "plot uncertainties"
             #to.printimage([ref,uncBRF_RMS_rand,uncBRF_RMS_syst],["BRF","random uncertainty","non-random uncertainty"],mini=[0,-0,-0],maxi=[1,0.05,0.05],save=0)
             #plottitels=["BRF","random uncertainty","non-random uncertainty","ratio non-random / random"]
             plottitels=["BRF","random uncertainty [absolute]","non-random uncertainty [absolute]","diff non-random - random [absolute]"]
             plotdata  =[ref,uncBRF_RMS_rand,uncBRF_RMS_syst,uncBRF_RMS_syst-uncBRF_RMS_rand]
             to.print4mat(plotdata,plottitels,f2,save=save)
             if debug>=1:
               print "\nTestsite non-random uncert: " + str( uncBRF_RMS_syst[tstsitepixlin] )+"\n"
             #plottitels=["BRF","fractional random uncertainty","fractional non-random uncertainty","ratio non-random / random"]
             plottitels=["BRF","Random Uncertainty[%]","Fractional Non-Random Uncertainty[%]","Diff (Non-Random - Random)"]
             uncBRF_RMS_rand_rel=(uncBRF_RMS_rand/ref)*100
             uncBRF_RMS_syst_rel=(uncBRF_RMS_syst/ref)*100
             plotdata  =[ref,uncBRF_RMS_rand_rel,uncBRF_RMS_syst_rel,uncBRF_RMS_syst_rel-uncBRF_RMS_rand_rel]
             to.print4mat(plotdata,plottitels,f2,save=save)
            if debug>=2:
              mi= np.median(uncBRF_RMS_syst[m2==1])-np.std(uncBRF_RMS_syst[m2==1])
              ma= np.median(uncBRF_RMS_syst[m2==1])+np.std(uncBRF_RMS_syst[m2==1])
              #to.printimage(uncBRF_RMS_syst,"uncBRF_RMS_syst",mini=np.median(uncBRF_RMS_syst[m2==1])-(np.std(uncBRF_RMS_syst[m2==1])),maxi=np.median(uncBRF_RMS_syst[m2==1])+(np.std(uncBRF_RMS_syst[m2==1])),save=save)
              to.printimage(uncBRF_RMS_syst,"non-random uncertainty",mini=mi,maxi=ma,save=save)
              to.printimage(uncBRF_RMS_syst-uncBRF_RMS_rand,"difference non-random minus random uncertainty",mini=mi,maxi=ma,save=save)

            #print "...calc combined BRF uncertainties"
            #timing.append(timeit.default_timer())#5
            #uncMS=np.zeros((N,N))
            #for targ,uval in zip(["Ce","Cs","E","SZA","a0","a1","Z"],[adev,adev,uirr,uSZA,np.sqrt(cf[3][0,0]),np.sqrt(cf[3][1,1]),0]):
              #unctmp=uval*ef.sensi_p(targ,pixvis,Cs,ysl,cf[2][0],cf[2][1],SZA,irr,sundist)
              #uncMS=uncMS+np.power(unctmp,2)
              #unctmp[unctmp>1]=1
              #unctmp[unctmp<-1]=-1
              #if debug>=1:
                #to.printimage(unctmp,"unc"+targ,mini=np.median(unctmp)-(np.std(unctmp)),maxi=np.median(unctmp)+(np.std(unctmp)),save=save)
                ##to.printimage(uncMS,"uncMS"+targ,mini=np.median(uncMS)-(np.std(uncMS)),maxi=np.median(uncMS)+(np.std(uncMS)))
            #uncBRF_RMS=np.sqrt(uncMS)
            #uncBRF_RMS[m2==0]=-1
            #uncBRF_RMS[uncBRF_RMS>1]=1
            #if debug>=1:
              #to.printimage(uncBRF_RMS,"uncBRF_RMS",mini=np.median(uncBRF_RMS[m2==1])-(np.std(uncBRF_RMS[m2==1])),maxi=np.median(uncBRF_RMS[m2==1])+(np.std(uncBRF_RMS[m2==1])),save=save)
            timing.append(timeit.default_timer())#6
            
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


  #===============================================================================================
  # Part 6 :write netcdf
  #-----------------------------------------------------------------------------------------------  
          u_srf  =  srf[:, 3:]
          if float(rel) >= 2.0:
            print "Writing easyFCDR data..."

            #filename++++++++++++++
            outeasy   =  outdir+wr.ncName(s_Year,s_Month,s_Day,s_Time,"EASY",nomlon,Software,rel,satellite)
            #if os.path.exists(outeasy):
              #os.remove(outeasy)
            
            if debug==1:
              print "outfile:"+outeasy
              
            # preparing netcdf file++
            writer = FCDRWriter()
            dataset = writer.createTemplateEasy("MVIRI", 12835)
            dataset.attrs["institution"] = "EUMETSAT"
            dataset.attrs["title"] = "MVIRI Easy FCDR"
            dataset.attrs["source"] = "Produced from UMARF RECT2LP and IMAG2TG data with MVIRI FCDR code RICalPy, version "+str(Software)
            dataset.attrs["history"] = "Created: " + time.ctime(time.time())
            dataset.attrs["references"] = "in preparation"
            dataset.attrs["comment"] = comment
            dataset.attrs["authors"] = "Frank Ruethrich, Viju John, Rob Roebeling and Joerg Schulz"
            dataset.attrs["email"] = "frank.ruethrich@eumetsat.int"
            dataset.attrs["satellite"] = satellite
            dataset.attrs["url"] = "http://www.fiduceo.eu"
            dataset.attrs["fcdr_software_version"] = str(Software)
            dataset.attrs["data_version"] = str(rel)
            dataset.attrs["RECT2LP_file_name"] = os.path.basename(filename)
            dataset.attrs["channels"]          = "vis, ir, wv"
            dataset.attrs["description"]       = "Meteosat First Generation Rectified (Level 1.5) Image"
            #preparing relative uncertainties
            uncBRF_RMS_rand_rel=abs(uncBRF_RMS_rand/ref)
            uncBRF_RMS_syst_rel=abs(uncBRF_RMS_syst/ref)
            # populating netcdf file++
            #flag
            dataset.variables["quality_pixel_bitmask"].data[:,:]                   = np.fliplr(pixel_mask)
            #time
            secsince1970[m1==0]=dataset.variables["time"].attrs['_FillValue']
            dataset.variables["time"].data[:,:]                   = np.fliplr(np.divide(secsince1970,1))#.astype(int))
            dataset.variables["time"].attrs['add_offset']         = timeoffset
            #angles
            #dataset.variables["satellite_azimuth_angle"].data[:,:]=np.fliplr(np.divide(VAZ,0.005493164).astype(int)) #now stored in static
            #dataset.variables["satellite_zenith_angle"].data[:,:] =np.fliplr(np.divide(VZA,0.005493248).astype(int)) #now stored in static
            SAZ=np.divide(SAZ,0.005493164)
            SAZ[m2==0]=dataset.variables["solar_azimuth_angle"].attrs['_FillValue']
            dataset.variables["solar_azimuth_angle"].data[:,:]    = np.fliplr(SAZ.astype(int))
            SZA=np.divide(SZA,0.005493164)
            SZA[m2==0]=dataset.variables["solar_zenith_angle"].attrs['_FillValue']
            dataset.variables["solar_zenith_angle"].data[:,:]     = np.fliplr(SZA.astype(int))
            #reflectance
            dtype=np.uint16
            ref=np.divide(ref,1.52588E-5)
            ref[m2==0]=DefaultData.get_default_fill_value(dtype)
            dataset.variables["toa_bidirectional_reflectance_vis"].data[:,:] =np.fliplr(ref.astype(dtype))
            dataset.variables["toa_bidirectional_reflectance_vis"].attrs["rms_landmarks_x_vis"] = lm_def[0]
            dataset.variables["toa_bidirectional_reflectance_vis"].attrs["rms_landmarks_y_vis"] = lm_def[1]
            #dataset.variables["count_vis"].data[:,:]             =np.fliplr(np.divide(pixvis,1).astype(int))
            #srf
            dataset.variables["spectral_response_function_vis"].data[:]=srf[:,1]#.astype(int))
            dataset.variables["spectral_response_function_vis"].attrs['version']  = srf_vers
            dataset.variables["spectral_response_function_vis"].attrs['source']  = os.path.basename(srfsource)
            dataset.variables["spectral_response_function_vis"].attrs['valid(YYYYDDD)']  = os.path.basename(srfsource)[9:9+7]
            dataset.variables["covariance_spectral_response_function_vis"].data[:,:]=u_srf#.astype(int))
            #uncertainties
            #random
            uncBRF_RMS_rand_rel_scale  = dataset.variables["u_random_toa_bidirectional_reflectance"].attrs['scale_factor']
            uncBRF_RMS_rand_rel        = np.divide(uncBRF_RMS_rand_rel,uncBRF_RMS_rand_rel_scale)
            uncBRF_RMS_rand_rel[m2==0] = dataset.variables["u_random_toa_bidirectional_reflectance"].attrs['_FillValue']
            dtype_u                    = dataset.variables["u_random_toa_bidirectional_reflectance"].dtype
            dataset.variables["u_random_toa_bidirectional_reflectance"].data[:,:]               = np.fliplr(uncBRF_RMS_rand_rel.astype(dtype_u))
            #non-random
            uncBRF_RMS_syst_rel_scale  = dataset.variables["u_non_random_toa_bidirectional_reflectance"].attrs['scale_factor']
            uncBRF_RMS_syst_rel        = np.divide(uncBRF_RMS_syst_rel,uncBRF_RMS_syst_rel_scale)
            uncBRF_RMS_syst_rel[m2==0] = dataset.variables["u_non_random_toa_bidirectional_reflectance"].attrs['_FillValue']
            dtype_u                    = dataset.variables["u_non_random_toa_bidirectional_reflectance"].dtype
            dataset.variables["u_non_random_toa_bidirectional_reflectance"].data[:,:]           = np.fliplr(uncBRF_RMS_syst_rel.astype(dtype_u))
            #u_latitude/u_longitude
            #ulats  = np.divide(ulats,dataset.variables["u_latitude"].attrs['scale_factor']).astype(int)
            #ulons  = np.divide(ulons,dataset.variables["u_longitude"].attrs['scale_factor']).astype(int)
            #ulats[m2==0]  = dataset.variables["u_latitude"].attrs['_FillValue']
            #ulons[m2==0]  = dataset.variables["u_longitude"].attrs['_FillValue']
            #dataset.variables["u_latitude" ].data[:,:]            = np.fliplr(ulats)
            #dataset.variables["u_longitude"].data[:,:]            = np.fliplr(ulons)
            #other 2D
            dataset.variables["count_wv"].data[:,:]               = np.fliplr(np.array(pixwv).astype(int))
            dataset.variables["count_ir"].data[:,:]               = np.fliplr(np.array(pixir).astype(int))
            #lmrks are calculated on IR image; so store them here as attributes
            dataset.variables["count_ir"].attrs["rms_landmarks_x_ir"] = lm_def[2]
            dataset.variables["count_ir"].attrs["rms_landmarks_y_ir"] = lm_def[3]
            dataset.variables["count_wv"].attrs["rms_landmarks_x_wv"] = lm_def[2]
            dataset.variables["count_wv"].attrs["rms_landmarks_y_wv"] = lm_def[3]
                        
            #scalars
            dataset.variables["solar_irradiance_vis"].data   = irr
            dataset.variables["u_solar_irradiance_vis"].data = uirr
            dataset.variables["distance_sun_earth"].data     = sundist
            
            dataset.variables["sub_satellite_longitude_start"].data = sslat
            dataset.variables["sub_satellite_latitude_start"].data  = sslon
            dataset.variables["sub_satellite_longitude_end"].data   = sslatE
            dataset.variables["sub_satellite_latitude_end"].data    = sslonE
            
              #day_number,a_ir(time),b_ir(time),sigma_a_ir(time),sigma_b_ir(time),q_ir(time),
              #chi_sq_ir(time),n_ir(time),rmsd_ir(time),corr_ir(time),min_mon_cnt3_ir(time),
              #max_mon_cnt3_ir(time),min_ref_rad_ir(time),max_ref_rad_ir(time),
              #num_per_bin_ir(time, nbin)a_wv(time)b_wv(time),sigma_a_wv(time),sigma_b_wv(time),
              #q_wv(time),chi_sq_wv(time),n_wv(time),rmsd_wv(time),corr_wv(time),
              #min_mon_cnt3_wv(time),max_mon_cnt3_wv(time),min_ref_rad_wv(time),
              #max_ref_rad_wv(time),num_per_bin_wv(time, nbin)
              #Of those are to be stored: [1,2,3,4,5,15,16,17,18,19]
            
            dataset.variables["a_ir"].data   = calfile_content[0][0]
            dataset.variables["b_ir"].data   = calfile_content[0][1]
            dataset.variables["u_a_ir"].data = calfile_content[0][2]
            dataset.variables["u_b_ir"].data = calfile_content[0][3]
            dataset.variables["q_ir"].data   = calfile_content[0][4]
            dataset.variables["a_wv"].data   = calfile_content[0][5]
            dataset.variables["b_wv"].data   = calfile_content[0][6]
            dataset.variables["u_a_wv"].data = calfile_content[0][7]
            dataset.variables["u_b_wv"].data = calfile_content[0][8]
            dataset.variables["q_wv"].data   = calfile_content[0][9]
            dataset.variables["unit_conversion_ir"].data = to.BTstuff(satellite)[0]
            dataset.variables["unit_conversion_wv"].data = to.BTstuff(satellite)[1]
            dataset.variables["bt_a_ir"].data = to.BTstuff(satellite)[2]
            dataset.variables["bt_a_wv"].data = to.BTstuff(satellite)[3]
            dataset.variables["bt_b_ir"].data = to.BTstuff(satellite)[4]
            dataset.variables["bt_b_wv"].data = to.BTstuff(satellite)[5]

            #writing netcdf file+++++
            writer.write(dataset, outeasy)
            
            
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
          #if float(rel) >= 2.0:


          ## content cf: 0=cf-uncert(t);1=cf(t);2=a;3=cov;

          if float(rel) < 2.0:
            images    =  ( pixvis, pixir, pixwv, t1) 
            uncert    =  ( cf[3][0,0],cf[3][1,1], digidev)
            refl_aux  =  ( adev, Cs, ysl, cf[2][0], cf[2][1])
            dimvals   =  []
          else:
            #images    =  ( pixvis, pixir, pixwv, ref, VZA, VAZ, SZA, SAZ, t1)
            images    =  ( pixvis, pixir, pixwv, SZA, SAZ, t1)
            uncert    =  ( uirr, u_srf, cf[3][0,0],cf[3][1,1],cf[3][0,1], ulats, ulons, ut, uSZA, digidev)
            #refl_aux  =  ( sundist, irr, srf[:,1:2].flatten(), adev, Cs, sslat, sslon, sslatE, sslonE, cf[2][0], cf[2][1] )
            refl_aux  =  ( sundist, irr, srf[:,1:2].flatten(), adev, Cs, ysl, cf[2][0], cf[2][1] )
            dimvals   =  ( srf[:,0].flatten(),0)
          # writing netcdf file++
          wr.write_manager( outnc, filename, head1, head2, head3, lininfo, trail, telem, telem_descr, \
                            images, uncert, refl_aux,calfile_content,dimvals,comment,rel,m2,pixel_mask,Software,satellite)


if __name__ == "__main__":
   main(sys.argv[1:])
