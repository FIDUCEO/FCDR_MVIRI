 
import cruncher as cr
import numpy as np
import sys
import matplotlib.pyplot as plt

sys.path.insert(0, '../../src/')
sys.path.insert(0, '../../lib/')
import ricals_tools as to

global lim

print cr.latlon.__doc__
print cr.sza.__doc__
print cr.radiance.__doc__
print cr.reflectance.__doc__

N=5000
lim=[0,N]
DOY=76

#lat/lon using fortran (very fast)
LAT,LON=cr.latlon(N)
LAT2,LON2=cr.latlon(N/2)

plt.plot(LAT[:,2410])
plt.show()
plt.plot(LON[2410,:])
plt.show()


to.printimage(LAT,"LAT",mini=-90,maxi=90)
to.printimage(LON,"LON",mini=-90,maxi=90)

T = 9.82

C = np.random.randint(0,255,(N,N))

#sza using fortran
m1                       =  np.ones((2500,2500),dtype=int)
m1[(LAT2<-90)|(LAT2>90)] =  0
print np.shape(m1)

SZA=cr.sza(T,DOY,m1,LAT,LON)
to.printimage(SZA,"SZA",mini=0,maxi=160)
print SZA[3803,3448]
print LAT[3803,3448]
print LON[3803,3448]
#reflectance with fortran (slightly slower)
rf=cr.reflectance(C,0.91,4.5,SZA,100,1)

#reflectance in python (faster)

rf2=ef.reflectance(C,0.91,4.5,1,100,SZA)

#radiance with fortran (fast)
r=cr.radiance(C,0.91,4.5)



