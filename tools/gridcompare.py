import matplotlib.pyplot as plt
import cruncher as cr
import msg_cruncher as msg_cr
import numpy as np
from numba import jit
from netCDF4 import Dataset

@jit
def indexsmallestdist(y,x,y1,x1,lat2,lon2,lat1_max,lon1_max):
  '''
  function that returns index of best match 
  of current seviri (y,x) in mviri line/pixel 
  coordinates and the distance in degree or meters
  inputs:
  y        = scalar
  x        = scalar
  lat2     = latitude grid of seviri
  lon2     = longitude grid of seviri
  lat1_max = latitude grid of mviri (corners)
  lon1_max = longitude grid of mviri (corners)
  '''
  #min dist of each seviri pixel center to anymviri pixel corner
  yc=3712/2-1
  xc=3712/2-1
  #searchradius around central seviri pixel as function of distance from center
  dc=np.sqrt((xc-x)**2+(yc-y)**2)
  r=int(dc*0.25+5)
  #distance in lat
  dlat=lat2[y,x]-lat1_max[y1-r:y1+r,x1-r:x1+r]
  #distance in lon
  dlon=lon2[y,x]-lon1_max[y1-r:y1+r,x1-r:x1+r]
  ##actual offset in degree
  #distmat = np.sqrt((dlat)**2+(dlon)**2)
  #actual offset in km
  distmat = np.sqrt((dlat*110.54)**2+(dlon*111.320*np.cos(np.radians(lat2[y,x])))**2)
  #prepare 2D index
  ya,xa   = np.unravel_index(np.nanargmin(distmat),distmat.shape)
  #return indext of best match in mviri line/pixel coordinates and distance
  return y1-r+ya , x1-r+xa , np.nanmin(distmat),dc

def smallestdistanceloop(N,lat2,lon2,lat1_max,lon1_max):
  '''
  this function loops over a grid and applies 
  the terms in function indexsmallestdist on 
  each cell
  inputs:
  N = gridsize (3712 for Seviri)
  '''
  bestx=np.zeros((N,N))
  besty=np.zeros((N,N))
  bestd=np.zeros((N,N))
  distcenter=np.zeros((N,N))
  bestx[:]=np.nan
  besty[:]=np.nan
  bestd[:]=np.nan
  distcenter=np.nan
  off=0
  for y in range(off,N-off):
    print y
    for x in range(off,N-off):
      try:
        #most likely mviri match (sort of roughly shifts the grid):
        y1,x1=cr.linepixel(2500,lat2[y,x],lon2[y,x],0)
        #check
        besty[y,x],bestx[y,x],bestd[y,x],distcenter[y,x]=indexsmallestdist(y,x,y1,x1,lat2,lon2,lat1_max,lon1_max)
      except ValueError:
        besty[y,x]=np.nan
        bestx[y,x]=np.nan
        bestd[y,x]=np.nan
        distcenter[y,x]=np.nan
  return bestx,besty,bestd,distcenter

#define grids
lat,lon=cr.latlon(1249)
lat1,lon1=cr.latlon(2500)
lat2,lon2=msg_cr.latlon(3712)

#handle image corners
lat1[lat1<-90]=np.nan
lat2[lat2<-90]=np.nan
lon1[lon1<-90]=np.nan
lon2[lon2<-90]=np.nan
lat[lat<-90]=np.nan
lon[lon<-90]=np.nan

#gradients 
lat1grad=np.gradient(lat1,axis=0)
lat2grad=np.gradient(lat2,axis=0)
lon1grad=np.gradient(lon1,axis=1)
lon2grad=np.gradient(lon2,axis=1)
#gradients in meters
ym1grad=np.gradient(lat1,axis=0)*110.54
ym2grad=np.gradient(lat2,axis=0)*110.54
xm1grad=np.gradient(lon1,axis=1)*111.320*np.cos(np.radians(lat1))
xm2grad=np.gradient(lon2,axis=1)*111.320*np.cos(np.radians(lat2))

#pixel corners
lat1_min=lat1-(lat1grad/2)
lat1_max=lat1+(lat1grad/2)
lat2_min=lat2-(lat2grad/2)
lat2_max=lat2+(lat2grad/2)
lon1_min=lon1-(lon1grad/2)
lon1_max=lon1+(lon1grad/2)
lon2_min=lon2-(lon2grad/2)
lon2_max=lon2+(lon2grad/2)

##assumption: 9 seviri pixel are always contained in 4 mviri pixel
##get pixel heights
#lat1siz=lat1_max[:,:]-lat1_min[:,:]
#lat2siz=lat2_max[:,:]-lat2_min[:,:]
##get pixel widths
#lon1siz=lon1_max[:,:]-lon1_min[:,:]
#lon2siz=lon2_max[:,:]-lon2_min[:,:]
##get pixel areas
#siz1=lat1siz*lon1siz
#siz2=lat2siz*lon2siz
##plot first test: are 4 mviri always as big as 9 seviri?
#plt.plot(lat1[:,2500/2],siz1[:,(2500)/2]*4)
#plt.plot(lat2[:,3712/2],siz2[:,(3712)/2]*9)
#plt.show()
#plt.plot(lon1[2500/2,:],siz1[(2500)/2,:]*4)
#plt.plot(lon2[3712/2,:],siz2[(3712)/2,:]*9)
#plt.show()
##assumption -> TRUE

#where do corners of MVIRI grid coincide with centers of SEVIRI grid?
N=3712
x,y,d,dc=smallestdistanceloop(3712,lat2,lon2,lat1_max,lon1_max)

#dlat=lat2-lat1_max[y.astype(np.int16),x.astype(np.int16)]
#dlon=lon2-lon1_max[y.astype(np.int16),x.astype(np.int16)]
#dist=np.sqrt(dlat**2+dlon**2)

#plt.imshow(d,vmax=1)
#plt.colorbar()
#plt.show()

#save as netcdf
filename    ="/DSNNAS/Repro_Temp/users/frankr/seviri2mviri_full.nc"
try:
  fo          = Dataset(filename,'w',format='NETCDF4')
except IOError:
  import os
  os.remove(filename)
  fo          = Dataset(filename,'w',format='NETCDF4')
fo.pixel    = np.arange(N)
fo.line     = np.arange(N)
line        = fo.createDimension('line', N)
pixel       = fo.createDimension('pixel', N)
dist        = fo.createVariable('minimum_distance',np.float64,('line','pixel'),fill_value=np.nan)
x_idx       = fo.createVariable('mviri_pixel',np.float64,('line','pixel'),fill_value=np.nan)
y_idx       = fo.createVariable('mviri_line' ,np.float64,('line','pixel'),fill_value=np.nan)

dist[:,:]   = d
dist.units  = "kilometers"
dist.longname = "distance of seviri pixel center to closest mviri corner"
dist.valid_min=0
dist.valid_max=3
x_idx[:,:]  = x
x_idx.units = "pixel"
x_idx.longname = "pixel number starting at 0 of mviri pixel the NW-corner of which matches with seviri pixel center"
x_idx.valid_min=0
x_idx.valid_max=N
y_idx[:,:]  = y
y_idx.units = "lines"
y_idx.longname = "line number starting at 0 of mviri pixel the NW-corner of which matches with seviri pixel center"
y_idx.valid_min=0
y_idx.valid_max=N
distcenter[:,:]   = dc
dist.units  = "pixel"
dist.longname = "distance of current seviri pixel to nadir pixel"
dist.valid_min=0
dist.valid_max=3712/2
fo.close()

exit()
