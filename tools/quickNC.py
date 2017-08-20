 
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap

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
    plt.imshow(val)
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