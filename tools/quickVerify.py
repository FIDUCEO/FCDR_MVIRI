
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import datetime as dt
#from mpl_toolkits.basemap import Basemap

def find_nearest(arraylat,arraylon,valuelat,valuelon):
    idx = np.argmin((np.abs(arraylat-valuelat)+np.abs(arraylon-valuelon)))
    idxX=idx/np.shape(arraylat)[0]
    idxY=idx%np.shape(arraylat)[1]
    print valuelat,valuelon
    print idxY,idxX
    print "LAT= "+str(arraylat[idxX,idxY])+" LON= "+str(arraylon[idxX,idxY])
    return idxX,idxY

startdate= "199801"#raw_input("enter startdate (YYYYMM)  : ")
enddate  = "200601"#raw_input("enter enddate (YYYYMM)    : ")
dset     = "EASY"#"raw_input("EASY or FULL              : ")
flag     = "y"#raw_input("test also contents (y/n)  : ")
if flag=="y":
  SSP      = "00.0"#raw_input("enter Sub Sat Point(XX.X) : ")
  targetVAR= "u_non_random_toa_bidirectional_reflectance"#raw_input("enter target variable     : ")
  ulLAT    = 28.7#float(raw_input("enter upper LAT     : "))
  lrLAT    = 28.4#float(raw_input("enter lower LAT     : "))
  lrLON    = 23.54#float(raw_input("enter bigger LON    : "))
  ulLON    = 23.24#float(raw_input("enter lower  LON    : "))
PBASE    = "/DSNNAS/Repro_Temp/mviri_dev/FCDR/MFG15_FCDR_V2.5_SRF0954/"
PATH     = PBASE+"data/MET7/"
PNAS     = "/DSNNAS/Repro/mviri/level1/HR-MFG15/data/MET7/"
yearS    = int(startdate[:4])
yearE    = int(enddate[:4])
monthS   = int(startdate[4:6])
monthE   = int(enddate[4:6])

ts=[]
dates=[]

if flag=="y":
  fstatic=PBASE+"aux/FIDUCEO_FCDR_L15_MVIRI_MET7-"+SSP+"_STATIC_v2.2_fv2.5.nc"
  static = Dataset(fstatic, mode='r')
  lons_vis = static.variables['longitude_vis'][:]
  lats_vis = static.variables['latitude_vis'][:]
  lons_ir  = static.variables['longitude_ir_wv'][:]
  lats_ir  = static.variables['latitude_ir_wv'][:]
  static.close()
  #get index
  if "_ir" in targetVAR or "_wv" in targetVAR or "time" in targetVAR:
    IDXulLAT = find_nearest(lats_ir,lons_ir,ulLAT,ulLON)[0]
    IDXulLON = find_nearest(lats_ir,lons_ir,ulLAT,ulLON)[1]
    IDXlrLAT = find_nearest(lats_ir,lons_ir,lrLAT,lrLON)[0]
    IDXlrLON = find_nearest(lons_ir,lons_ir,lrLAT,lrLON)[1]
    disp=lats_ir
  else:
    IDXulLAT = find_nearest(lats_vis,lons_vis,ulLAT,ulLON)[0]
    IDXulLON = find_nearest(lats_vis,lons_vis,ulLAT,ulLON)[1]
    IDXlrLAT = find_nearest(lats_vis,lons_vis,lrLAT,lrLON)[0]
    IDXlrLON = find_nearest(lats_vis,lons_vis,lrLAT,lrLON)[1]
    disp=lats_vis
  #to display selected area:
  disp[IDXulLAT-5:IDXulLAT+5,IDXulLON-5:IDXulLON+5]=180
  disp[IDXlrLAT-5:IDXlrLAT+5,IDXulLON-5:IDXulLON+5]=180
  disp[IDXulLAT-5:IDXulLAT+5,IDXlrLON-5:IDXlrLON+5]=-180
  disp[IDXlrLAT-5:IDXlrLAT+5,IDXlrLON-5:IDXlrLON+5]=-180
  plt.imshow(disp)
  plt.title("selected area: ")
  plt.colorbar()
  plt.show()
q=1
for year in range(yearS,yearE+1):
  PATHy=PATH+str(year)+"/"
  if year==yearS:
    firstmonth=monthS
  else:
    firstmonth=1
  if year==yearE:
    lastmonth=monthE
  else:
    lastmonth=1
  for month in range(firstmonth,lastmonth+1):
    PATHm = PATHy+str(month).zfill(2)+"/"
    try:
      Ndays=(dt.date(year, month+1, 1) - dt.date(year, month, 1)).days
    except ValueError:
      Ndays=(dt.date(year+1, 1, 1) - dt.date(year, month, 1)).days
    for day in range(1,Ndays+1):
      wildcard="*"+str(year)+str(month).zfill(2)+str(day).zfill(2)+"*_"+dset+"*.nc"
      #print wildcard
      filelist=sorted(glob.glob(PATHm+wildcard))
      
      if len(filelist)<48:
        print "\n"+str(48-len(filelist))+" files missing for "+str(day)+"/"+str(month)+"/"+str(year)
        wildcardNAS="*"+str(year)+str(month).zfill(2)+str(day).zfill(2)+"*"
        PATHn=PNAS+str(year)+"/"+str(month).zfill(2)+"/"
        filelistNAS=sorted(glob.glob(PATHn+wildcardNAS))
        print str(48-len(filelistNAS))+" files missing on NAS\n"
        if not len(filelist)==len(filelistNAS):
          print "WARNING: not all files from NAS processed!"
      if flag=="y":
        for f in filelist:
          print f
          tstrIDX=f.index(str(year)+str(month).zfill(2)+str(day).zfill(2))
          tstr   =f[tstrIDX+8:tstrIDX+12]
          hour   =int(tstr[:2])
          minute =int(tstr[2:4])
          #print tstr
          with Dataset(f, mode='r') as current:
            if q==1 and hour==12:
              disp=current.variables[targetVAR][:,:]
              #to display selected area:
              mi=np.amin(disp)
              mx=np.amax(disp)
              pc=(mx-mi)/100*30
              disp[IDXulLAT-5:IDXulLAT+5,IDXulLON-5:IDXulLON+5]=mx+pc
              disp[IDXlrLAT-5:IDXlrLAT+5,IDXulLON-5:IDXulLON+5]=mx+pc
              disp[IDXulLAT-5:IDXulLAT+5,IDXlrLON-5:IDXlrLON+5]=mi-pc
              disp[IDXlrLAT-5:IDXlrLAT+5,IDXlrLON-5:IDXlrLON+5]=mi-pc
              plt.imshow(disp)
              plt.title("selected area: ")
              plt.colorbar()
              plt.show()
              q=q+1
            if IDXulLAT<IDXlrLAT:
              val = current.variables[targetVAR][IDXulLAT:IDXlrLAT,IDXulLON:IDXlrLON]
              #try:
                #offset=current.variables[targetVAR].add_offset
                #val=val[:]+offset
              #except:
                #print "attr missing"
            else:
              val = current.variables[targetVAR][IDXlrLAT:IDXulLAT,IDXulLON:IDXlrLON]
              #try:
                #offset=current.variables[targetVAR].add_offset
                #val=val[:]+offset
              #except:
                #print "attr missing"
            if "time" in targetVAR:
              #print np.shape(val[:])
              i,j=np.shape(val[:])
              v=np.empty((i,j),dtype=object)
              days=val[int(i/2),int(j/2)]/86400.
              #print days
              v=dt.datetime(1970, 1, 1,0,0) + dt.timedelta(days)
              #print v
              val=v
               
            #print offset
            #print val[:]
            #print dir(val)
            #plt.imshow(val)
            #plt.colorbar()
            #plt.show()
            try:
              ts.append(np.mean(val))
            except TypeError:
              ts.append(val)
            dates.append(dt.datetime(year, month, day, hour, minute))
            #print f
            #print dt.datetime(year, month, day, hour, minute)
            #print val
            
if flag=="y":
 for i in range(len(dates)):
   print dates[i],ts[i]
 plt.plot(dates,ts,".")
 plt.show()
