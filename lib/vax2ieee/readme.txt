 
Example how the MVIRI VIS rect2lp reader can be used in python.

reads the header and returns RMS of landmarks

For compiling:

ssh -X tclxs18
cd ~/src_c/vax2ieee/SWIG
#compile
rm *_wrap.c*
rm *.so
rm *.o
swig -python py_vax2ieee.i
g++ -fPIC -c STAMP_PDL_Vax2Ieee.c py_vax2ieee_wrap.c -I/usr/include/python
g++ -fPIC STAMP_PDL_Vax2Ieee.c convert_vax_data.c -c
g++ -shared convert_vax_data.o STAMP_PDL_Vax2Ieee.o py_vax2ieee_wrap.o -o _vax2ieee.so
#test
python test_swig.py


