import subprocess
import multiprocessing
import sys

sat=raw_input("which satellite? (e.g.:"+'"MET7"'+"):")
comment='"pre-beta release!! Use only for illustration; acqtime fixed; SZA fixed; solar spectrum from rayference"'
release="2.7"
srfvers="1000"
s_s=sys.argv[1]
s_e=sys.argv[2]
s=int(s_s[4:])
e=int(s_e[4:])
sy=int(s_s[:4])
ey=int(s_e[:4])
for year in range(sy,ey+1):
  s_y=str(year)
  sj=1
  ej=366
  if year==sy:
    sj=s
  elif year==ey:
    ej=e
  for doy in range(sj,ej):
    for t in range(0,24):
      for m in ["00","30"]:
        #t=12
        #m ="00"
        time=str(t).zfill(2)+m
        #execute
        sd=s_y+str(doy).zfill(3)
        ed=s_y+str(doy).zfill(3)
        #sat='"MET7"'
        cmd="python ricals_main.py -s "+sd+" -e "+ed+" -m "+sat+" -t "+time+" -r "+release+" -v "+srfvers+" -a /DSNNAS/Repro_Temp/mviri_dev/FCDR/ -c "+comment
        print cmd
        try:
          dt=subprocess.call(cmd,shell=True)
        except OSError:
          print "ERROR"
