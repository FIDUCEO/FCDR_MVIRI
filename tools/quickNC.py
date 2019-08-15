 
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
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

filename=raw_input("enter filename (full path):")

fh = Dataset(filename, mode='r')
print fh.variables
try:
  lons = fh.variables['lon'][:]
  lats = fh.variables['lat'][:]
except:
  print ""
Dims=raw_input("No of Dimensions: ")
m="y"
while m=="y":
  Target=raw_input("Target Var name: ")

  if int(Dims) ==2:
    val = fh.variables[Target][:]
    mi,ma=getrange(val)
    plt.imshow(val,vmin=mi,vmax=ma)
    plt.colorbar()
    plt.show()

  if int(Dims) ==3:
    fst=raw_input("Plot which slide on first axis? ")
    val = fh.variables[Target][int(fst),:,:]
    plt.imshow(val)
    plt.colorbar()
    plt.show()
    
  if int(Dims) ==4:
    fst=raw_input("Plot which slide on first axis? ")
    scn=raw_input("Plot which slide on second axis? ")
    val = fh.variables[Target][int(fst),int(scn),:,:]
    print val [2500,2500]
    plt.imshow(val)
    plt.colorbar()
    plt.show()
    
  m=raw_input("plot more? (y/n)")
fh.close()