conversion of orbit coordinates 
from 
  mean geocentric frame of date 
to 
  earth fixed frame
AUTHOR: Frank Ruethrich EUMETSAT
COMPILE WITH:
f2py -c -m orbit2eff orbit2eff.f COOT.FOR WOBBLE.FOR

in python use as:
import orbit2eff
orbit2eff.mgf2eff.__doc__