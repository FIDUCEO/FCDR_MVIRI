 
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = (18,(18/3)*1.8)
import sys
sys.path.insert(0, '../src/')
import ricals_netcdf as wr
from datetime import datetime
#from mpl_toolkits.basemap import Basemap

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

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Software="2.6"
rel="3.1"
path="/DSNNAS/Repro/mviri/level1/MFG15_FCDR_V1/"

print "Software = ",Software
print "rel      = ",rel
default      =raw_input("Use default scene? (Y/n):")
defproj='geos'#'hammer'
if default=="n":
  satellite    =raw_input("Enter satellite (MET7):")
  nomlon       =raw_input("Enter nominal longtude (00.0):")
  s_Year       =raw_input("Enter year: (2000)")
  s_Month      =raw_input("Enter month: (03)")
  s_Day        =raw_input("Enter day: (15)")
  s_Time       =raw_input("Enter time:(1200)")
else:
  satellite    ="MET7"
  nomlon       ="00.0"
  s_Year       ="2000"
  s_Month      ="03"
  s_Day        ="15"
  s_Time       ="1100"

date = datetime(int(s_Year),int(s_Month),int(s_Day),int(s_Time[:2]),int(s_Time[2:]))
datadir=path+"data/"+satellite+'/'+s_Year+'/'+s_Month+'/'
auxdir=path+"aux/"
filename_easy=datadir+wr.ncName(s_Year,s_Month,s_Day,s_Time,"EASY",nomlon,Software,rel,satellite)
filename_static=auxdir+wr.ncName(s_Year,s_Month,s_Day,s_Time,"STATIC",nomlon,Software,rel,satellite)

print filename_easy
print filename_static

#read

fh_e = Dataset(filename_easy, mode='r')
fh_s = Dataset(filename_static, mode='r')

#print fh_e.variables

lons = fh_s.variables['longitude_vis'][:]
lats = fh_s.variables['latitude_vis'][:]

#msk = fh_e.variables['quality_pixel_bitmask'][:]
msk = fh_e.variables['data_quality_bitmask'][:]
data = fh_e.variables['toa_bidirectional_reflectance_vis'][:]
indep = fh_e.variables['u_independent_toa_bidirectional_reflectance'][:]
struct = fh_e.variables['u_structured_toa_bidirectional_reflectance'][:]
IR = fh_e.variables['count_ir'][:]
WV = fh_e.variables['count_wv'][:]
IR_a = fh_e.variables['a_ir'][:]
WV_a = fh_e.variables['a_wv'][:]
IR_b = fh_e.variables['b_ir'][:]
WV_b = fh_e.variables['b_wv'][:]

data=np.ma.masked_where((msk != 0)&(msk != 4)&(msk != 32)&(msk != 36), data)
indep=np.ma.masked_where((msk != 0)&(msk != 4)&(msk != 32)&(msk != 36), indep)
struct=np.ma.masked_where((msk != 0)&(msk != 4)&(msk != 32)&(msk != 36), struct)
IR=np.ma.masked_where((IR == 0)|(IR == 256)|(np.isnan(IR)), IR)*IR_b+IR_a
WV=np.ma.masked_where((WV == 0)|(WV == 256)|(np.isnan(WV)), WV)*WV_b+WV_a



fh_e.close()
fh_s.close()

#plot
dec=raw_input("plot BRF?(Y/n)")
if not dec=="n":
  comap=plt.cm.viridis
  comap.set_bad(color='white')
  fig = plt.figure(figsize=(10,10))
  ax = fig.add_axes([0.1,0.1,0.8,0.8])
  m = Basemap(projection=defproj,lon_0=float(nomlon))
  m.drawmapboundary(fill_color='white')
  m.drawcoastlines()#linewidth=0.25)
  im1 = m.pcolormesh(lons[::2,::2],lats[::2,::2],data[::2,::2],shading='gouraud',edgecolors='None',cmap=comap,vmin=0.01,vmax=1.0,latlon=True)
  #im1 = m.pcolor(lons,lats,data,shading='flat',cmap=comap,latlon=True)
  m.drawparallels(np.arange(-90.,99.,30.))
  m.drawmeridians(np.arange(-180.,180.,30.))
  cb = m.colorbar(im1,"bottom",extend='both', size="4%", pad="2%")
  cb.ax.tick_params(labelsize=15)
  cb.set_label("reflectance", fontsize=20)
  # add a title.
  ax.set_title('toa reflectance vis band \n%s'%date, fontsize=30)
  plt.show()
dec=raw_input("plot independent uncertainty?(Y/n)")
if not dec=="n":
  indep=indep*100.0
  comap=plt.cm.plasma
  comap.set_bad(color='white')
  fig = plt.figure(figsize=(10,10))
  ax = fig.add_axes([0.1,0.1,0.8,0.8])
  m = Basemap(projection=defproj,lon_0=float(nomlon))
  m.drawmapboundary(fill_color='white')
  m.drawcoastlines()#linewidth=0.25)
  im1 = m.pcolormesh(lons[::2,::2],lats[::2,::2],indep[::2,::2],shading='flat',cmap=comap,vmin=0.0,vmax=3.0,latlon=True)
  #im1 = m.pcolor(lons,lats,indep,shading='flat',cmap=comap,latlon=True)
  m.drawparallels(np.arange(-90.,99.,30.))
  m.drawmeridians(np.arange(-180.,180.,30.))
  cb = m.colorbar(im1,"bottom",extend='both', size="4%", pad="2%")
  cb.ax.tick_params(labelsize=15)
  cb.set_label("%", fontsize=20)
  # add a title.
  ax.set_title('independent uncertainty \n%s'%date, fontsize=30)
  plt.show()
dec=raw_input("plot structured uncertainty?(Y/n)")
if not dec=="n":
  struct=struct*100.0
  comap=plt.cm.plasma
  comap.set_bad(color='white')
  fig = plt.figure(figsize=(10,10))
  ax = fig.add_axes([0.1,0.1,0.8,0.8])
  m = Basemap(projection=defproj,lon_0=float(nomlon))
  m.drawmapboundary(fill_color='white')
  m.drawcoastlines()#linewidth=0.25)
  im1 = m.pcolormesh(lons[::2,::2],lats[::2,::2],struct[::2,::2],shading='flat',cmap=comap,vmin=0.0,vmax=1.0,latlon=True)
  #im1 = m.pcolor(lons,lats,struct,shading='flat',cmap=comap,latlon=True)
  m.drawparallels(np.arange(-90.,99.,30.))
  m.drawmeridians(np.arange(-180.,180.,30.))
  cb = m.colorbar(im1,"bottom",extend='both', size="4%", pad="2%")
  cb.ax.tick_params(labelsize=15)
  cb.set_label("%", fontsize=20)
  # add a title.
  ax.set_title('structured uncertainty \n%s'%date, fontsize=30)
  plt.show()
dec=raw_input("plot IR?(Y/n)")
if not dec=="n":
  comap=plt.cm.viridis
  comap.set_bad(color='white')
  fig = plt.figure(figsize=(10,10))
  ax = fig.add_axes([0.1,0.1,0.8,0.8])
  m = Basemap(projection=defproj,lon_0=float(nomlon))
  m.drawmapboundary(fill_color='white')
  m.drawcoastlines()#linewidth=0.25)
  im1 = m.pcolormesh(lons[::2,::2],lats[::2,::2],IR,shading='flat',cmap=comap,latlon=True)
  #im1 = m.pcolor(lons,lats,IR,shading='flat',cmap=comap,latlon=True)
  m.drawparallels(np.arange(-90.,99.,30.))
  m.drawmeridians(np.arange(-180.,180.,30.))
  cb = m.colorbar(im1,"bottom",extend='both', size="4%", pad="2%")
  cb.ax.tick_params(labelsize=15)
  cb.set_label("mW/m2/sr/cm-1 ", fontsize=20)
  # add a title.
  ax.set_title('toa radiance infrared band \n%s'%date, fontsize=30)
  plt.show()
dec=raw_input("plot WV?(Y/n)")
if not dec=="n":
  comap=plt.cm.viridis
  comap.set_bad(color='white')
  fig = plt.figure(figsize=(10,10))
  ax = fig.add_axes([0.1,0.1,0.8,0.8])
  m = Basemap(projection=defproj,lon_0=float(nomlon))
  m.drawmapboundary(fill_color='white')
  m.drawcoastlines()#linewidth=0.25)
  im1 = m.pcolormesh(lons[::2,::2],lats[::2,::2],WV,shading='flat',cmap=comap,latlon=True)
  #im1 = m.pcolor(lons,lats,WV,shading='flat',cmap=comap,latlon=True)
  m.drawparallels(np.arange(-90.,99.,30.))
  m.drawmeridians(np.arange(-180.,180.,30.))
  cb = m.colorbar(im1,"bottom",extend='both', size="4%", pad="2%")
  cb.ax.tick_params(labelsize=15)
  cb.set_label("mW/m2/sr/cm-1", fontsize=20)
  # add a title.
  ax.set_title('toa radiance water vapor absorption band \n%s'%date, fontsize=30)
  plt.show()
