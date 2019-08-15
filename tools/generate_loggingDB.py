import os,sys

from sqlalchemy import *
from sqlalchemy.ext.declarative import declarative_base

import psycopg2

import datetime

global badfiles

badfiles=[
  '/DSNNAS/Repro/mviri/level1/HR-MFG15/data/MET7/2000/06/METEOSAT7-MVIRI-MTP15-NA-NA-20000601040000.000000000Z',
  '/DSNNAS/Repro/mviri/level1/HR-MFG15/data/MET7/2000/06/METEOSAT7-MVIRI-MTP15-NA-NA-20000609123000.000000000Z',
  '/DSNNAS/Repro/mviri/level1/HR-MFG15/data/MET7/2000/08/METEOSAT7-MVIRI-MTP15-NA-NA-20000819053000.000000000Z',
  '/DSNNAS/Repro/mviri/level1/HR-MFG15/data/MET5/2000/05/METEOSAT5-MVIRI-MTP15-NA-NA-20000501180000.000000000Z',
  '/DSNNAS/Repro/mviri/level1/HR-MFG15/data/MET5/2000/06/METEOSAT5-MVIRI-MTP15-NA-NA-20000604160000.000000000Z',
  '/DSNNAS/Repro/mviri/level1/HR-MFG15/data/MET5/2000/08/METEOSAT5-MVIRI-MTP15-NA-NA-20000815180000.000000000Z'
  ]

Base = declarative_base()

import numpy as np

from mviri.mviri_l15 import read_rect2lp as rd

import cv2

sys.path.insert(0, '../src/')
import ricals_netcdf as wr

import calib

class ToDo(Base):
  __tablename__     = 'logtable'
  idn               = Column(BigInteger,unique=True,primary_key=True)
  rect2lp_filename  = Column(String(255))
  platform          = Column(String(4))
  filename_easy     = Column(String(255))
  filename_full     = Column(String(255))
  scan_ts           = Column(DateTime)
  ricalpy_ts        = Column(DateTime)
  input_l10         = Column(Integer) #shall be: 1=available -1=notavailable
  input_l15         = Column(Integer) #shall be: 1=available  0=anomaly -1=notavailable
  status_easy       = Column(Integer) #shall be: 0=pending -1=failed 1=successful
  status_full       = Column(Integer) #shall be: 0=pending -1=failed 1=successful
  problem_easy      = Column(String(1000))
  problem_full      = Column(String(1000))
  problem_l15       = Column(String(1000))
def satnamesR(sat):
  satdict={"MET2":"M2","MET3":"P2","MET4":"M4","MET5":"M5","MET6":"M6","MET7":"M7"}
  return satdict[sat]

def satnamesL(sat):
  satdict={"MET2":"METEOSAT2","MET3":"METEOSAT3","MET4":"METEOSAT4","MET5":"METEOSAT5","MET6":"METEOSAT6","MET7":"METEOSAT7"}
  return satdict[sat]


def check_easy(filename):
  '''
  checker for easy fcdr file
  '''
  error=1
  prob=""
  #first check filesize
  st = os.stat(filename)
  if (st.st_size/1000000)<10:
    error=-1
    prob="file too small"
  return (error,prob)

def check_full(filename):
  '''
  checker for full fcdr file
  '''
  error=1
  prob=""
  #first check filesize
  st = os.stat(filename)
  if (st.st_size/1000000)<5:
		error=-1
		prob="file too small"
  return (error,prob)

def check_rect2lp(filename,comp,anomalycheck=0):
  '''
  checker for rect2lp file
  '''
  error=1
  prob=""
  #first check filesize
  st = os.stat(filename)
  if not (st.st_size/1000000)==38:
    error=-1
    prob="filesize wrong"
  if anomalycheck != 0:
    #read file
    print "..read.."
    with open(filename, "rb") as fi:
      fi, head = rd.header(fi, 1, 3,t=0)
      nomlon = np.rad2deg(getattr(head[1],'dsNominalSubSatelliteLongitude')[0])
      lininfo, pixvis, pixir, pixwv = rd.image_3(fi)
      fi, trail = rd.header(fi, 2504, 2504, t=0)
    #get comparefile
    solartimedelta = nomlon/15 #sun moves 15 degrees lon per hour
    slotdelta      = int(solartimedelta*2+0.5)
    compareslot    = int(getattr(head[0],'iSlotNum')[0]) + slotdelta
    if compareslot>=48:
      compareslot=compareslot-48
    try:
      #H1=np.histogram(comp[compareslot,0,:,:], bins=256, range=[0, 256],density=True)
      #H1 = cv2.calcHist([cv2.equalizeHist(comp[compareslot,0,:,:].astype(np.uint8))],[0],None,[256],[0,256])
      H1 = cv2.calcHist([comp[compareslot,0,:,:].astype(np.uint8)],[0],None,[256/2],[0,256/2])
      H1=cv2.normalize(H1,None)
    except IndexError:
      error = -1
      prob="slot in metadata wrong: "+getattr(head[0],'iSlotNum')[0]
      return (error,prob)
    #H2=np.histogram(np.array(pixvis), bins=256, range=[0, 256],density=True)
    #H2 = cv2.calcHist([cv2.equalizeHist(np.array(pixvis).astype(np.uint8))],       [0],None,[256],[0,256])
    H2 = cv2.calcHist([np.array(pixvis).astype(np.uint8)],       [0],None,[256/2],[0,256/2])
    H2 = cv2.normalize(H2,None)
    #chi2=chi2_distance(H1[0],H2[0])
    #correlation test
    test=cv2.compareHist(H1, H2, cv2.HISTCMP_CORREL)
    print test
    if "METEOSAT2" in filename or "METEOSAT3" in filename:
      thr=0.4#for MET2 and 3 the difference might be quite large
    else:
      thr=0.6#relatively arbitrary threshold
    if test<thr:
      error = 0
      prob="histogram anomalous with R2="+str(test)
      return (error,prob)
    #chi2 test  
    test=cv2.compareHist(H1, H2, cv2.HISTCMP_CHISQR)
    print test
    if "METEOSAT2" in filename or "METEOSAT3" in filename:
      thr=300#for MET2 and 3 the difference might be quite large
    else:
      thr=200#relatively arbitrary threshold
    if test>thr:
      error = 0
      prob="histogram anomalous with Chi2="+str(test)
      return (error,prob)	
    #intersect test
    test=cv2.compareHist(H1, H2, cv2.HISTCMP_INTERSECT)
    print test
    if "METEOSAT2" in filename or "METEOSAT3" in filename:
      thr=0.25#for MET2 and 3 the difference might be quite large
    else:
      thr=0.5#relatively arbitrary threshold
    if test<thr:
      error = 0
      prob="histogram anomalous with intersect="+str(test)
      return (error,prob)
    test=cv2.compareHist(H1, H2, cv2.HISTCMP_BHATTACHARYYA)
    print test
    if "METEOSAT2" in filename or "METEOSAT3" in filename:
      thr=1#for MET2 and 3 the difference might be quite large
    else:
      thr=0.5#relatively arbitrary threshold
    if test>thr:
      error = 0
      prob="histogram anomalous with Bhattacharya="+str(test)
  return (error,prob)

def chi2_distance(histA, histB, eps = 1e-10):
	# compute the chi-squared distance
	d = 0.5 * np.sum([((a - b) ** 2) / (a + b + eps)
		for (a, b) in zip(histA, histB)])
	# return the chi-squared distance
	return d

def readcomparefiles():
	'''
	reads a number of files that can act as 
	a sample to find anomalies in other files
	'''
	print "get samples for histogram comparison"
	comp=np.empty((48,3,5000,5000))
	for slot in range(0,48):
		print slot
		comparetime_dec= slot/2
		comparetime_str= str(int(comparetime_dec)).zfill(2)+str((comparetime_dec-int(comparetime_dec))*60).zfill(2)
		filename="/DSNNAS/Repro/mviri/level1/HR-MFG15/data/MET7/2005/12/METEOSAT7-MVIRI-MTP15-NA-NA-20051206{sl}00.000000000Z".format(sl=comparetime_str)
		with open(filename, "rb") as fi:
			lininfo, pixvis, pixir, pixwv = rd.image_3(fi)
		print "store vis"
		comp[slot,0,:,:]        = np.array(pixvis)
		print "store ir"
		comp[slot,1,:2500,:2500]= np.array(pixir)
		print "store wv"
		comp[slot,2,:2500,:2500]= np.array(pixwv)
	return comp

def enter(sat,FCDR_pth,rel,anomalycheck=0):
  tablename="logtable"
  sat_short=satnamesR(sat)
  sat_long =satnamesL(sat)
  co       = calib.calib(sat,rel)
  s_str    = co.launchdate()
  e_str    = co.fcdr_enddate()
  sy,sm,sd = (int(s_str[:4]),int(s_str[4:6]),int(s_str[6:]))
  ey,em,ed = (int(e_str[:4]),int(e_str[4:6]),int(e_str[6:]))
  start_datetime   = datetime.datetime(year=sy,month=sm,day=sd)
  current_datetime = start_datetime
  end_datetime     = datetime.datetime(year=ey,month=em,day=ed)
  delta_time       = datetime.timedelta(minutes=30)
  host='10.5.194.201'
  conn = psycopg2.connect(dbname='ricalpy_{v1}'.format(v1=rel.replace('.','_')),user=os.environ.get('USER'), host=host, password='')
  conn.autocommit = True
  c=conn.cursor()
  if anomalycheck !=0:
    comp=readcomparefiles()
  else:
    comp=None
  while current_datetime<=end_datetime:
    idn=int(sat[3:]+current_datetime.strftime("%Y%m%d%H%M"))
    outdir=FCDR_pth+"data/"+sat+"/"+current_datetime.strftime("%Y/%m/")
    print outdir
    c.execute("SELECT * from {tn} WHERE {idf} = {fn}".format(tn=tablename,idf='idn',fn=idn))
    row = c.fetchone()
    if row:
      print "Already done {id}".format(id=idn)
    else:
      L15_longpath="{p}{s}/{y}/{m}/".format(p=L15_path,s=sat,y=str(current_datetime.year),m=str(current_datetime.month).zfill(2))
      L15_file="{lp}{ls}-MVIRI-MTP15-NA-NA-{ds}00.000000000Z".format(lp=L15_longpath,ls=sat_long,ds=current_datetime.strftime("%Y%m%d%H%M"))
      try:
        L10_file=mviri_l10.l15_to_L10_name(L15_file,current_datetime.strftime("%Y%m%d"))
        L10_error=1
      except IOError:
        L10_error=-1
      print "L15"
      if not L15_file in badfiles:
        try:
          with open(L15_file, "rb") as fi:
            print "..check.. "+L15_file
            L15_error,prob_L15=check_rect2lp(L15_file,comp,anomalycheck)
            fi, head  = mviri_l15.header(fi, 1, 3,t=0)
            #L15_error=1#L15_file is there
        except IOError:
          L15_error=-1 #L15_file is not there
          prob_L15 = "missing"
      else:
        L15_error=-1
        prob_L15 = "missing"
      try:
        nomlon      = np.rad2deg(getattr(head[1],'dsNominalSubSatelliteLongitude')[0])
        easy_file   = outdir+wr.ncName(str(current_datetime.year),str(current_datetime.month).zfill(2),str(current_datetime.day).zfill(2),current_datetime.strftime("%H%M"),"EASY",nomlon,Software,rel,sat)
        full_file   = outdir+wr.ncName(str(current_datetime.year),str(current_datetime.month).zfill(2),str(current_datetime.day).zfill(2),current_datetime.strftime("%H%M"),"FULL",nomlon,Software,rel,sat)
      except NameError:
        nomlon= np.nan
        easy_file=""
        full_file=""
      try:
        ricalpy_ts=datetime.datetime.fromtimestamp(os.path.getmtime(easy_file))
        easy_error=1#easy is there
        prob_easy=""
        if os.path.exists(full_file):#check also full
          full_error=1
          prob_full=""
        else:
          full_error=-1
          prob_full="missing"
      except OSError:
        easy_error=-1#easy is not there but eventually full
        prob_easy="missing"
        try:
          ricalpy_ts=datetime.datetime.fromtimestamp(os.path.getmtime(full_file))
          full_error=1
          prob_full=""
        except OSError:#both missing
          easy_error=0
          full_error=0
          prob_easy="missing"
          prob_full="missing"
          ricalpy_ts=datetime.datetime(1970,1,1)
      print "content"
      if full_error==1 and easy_error==1:#if both files there, check content
        easy_error,prob_easy=check_easy(easy_file)
        full_error,prob_full=check_full(full_file)
      columnstring="({idx},{a},{b},{c},{d},{e},{f},{g},{h},{i},{j},{k},{l},{m})".format(
                      idx="idn",
                      a="rect2lp_filename",
                      b="platform",
                      c="filename_easy",
                      d="filename_full",
                      e="scan_ts",
                      f="ricalpy_ts",
                      g="input_l10",
                      h="input_l15",
                      i="status_easy",
                      j="status_full",
                      k="problem_easy",
                      l="problem_full",
                      m="problem_l15")
      valuestring="({idx},'{a}','{b}','{c}','{d}','{e}','{f}',{g},{h},{i},{j},'{k}','{l}','{m}')".format(
                      idx=idn,
                      a=L15_file,
                      b=sat,
                      c=easy_file,
                      d=full_file,
                      e=current_datetime,
                      f=ricalpy_ts,
                      g=L10_error,
                      h=L15_error,
                      i=easy_error,
                      j=full_error,
                      k=prob_easy,
                      l=prob_full,
                      m=prob_L15)
      print "write db"
      try:
        sql_string="INSERT INTO {tn} {col} VALUES {val}".format(tn=tablename, col=columnstring, val=valuestring)
        print sql_string
        c.execute(sql_string)
      except psycopg2.IntegrityError:
        print('ID already exists in PRIMARY KEY column {idx}'.format(idx=idn))
    current_datetime = current_datetime+delta_time
  conn.close()
  return

if __name__ == "__main__":
  rel      = raw_input("Create DB for which file release? (e.g. 2.5) ")
  Software = raw_input("Create DB for which software release? (e.g. 2.2) ")
  tablename="logtable"
  
  host='10.5.194.201'#'localhost'
  
  con = psycopg2.connect(dbname='postgres', user=os.environ.get('USER'), host=host, password='')
  con.autocommit = True
  cur = con.cursor()
  a="N"
  
  #res=cur.execute("select exists(select * from information_schema.tables where table_name=%s)", (tablename,))
  try:
    print "creating db"
    cur.execute('CREATE DATABASE ricalpy_{v1};'.format(v1=rel.replace('.','_')))
    res=False
  except psycopg2.ProgrammingError:
    res=True
  
  if res==False:
    a='y'
  else:
    a=raw_input("db exists; reset? (y/N) ") 
    if a=='y':
      print "re-creating db"
      cur.execute('DROP DATABASE ricalpy_{v1};'.format(v1=rel.replace('.','_')))
      cur.execute('CREATE DATABASE ricalpy_{v1};'.format(v1=rel.replace('.','_')))
  if a=='y':
    db_string='postgresql://{n}@{h}:5432/ricalpy_{v1}'.format(n=os.environ.get('USER'),h=host,v1=rel.replace('.','_'))
    engine = create_engine(db_string)
    print "creating table"
    Base.metadata.create_all(bind=engine)
  
  a=raw_input("Also fill with existing FCDR content? (y/N) ")
  anomalycheck = int(raw_input("include anomaly-check? (0=n/1=y)"))
    
  if a=="y":
    
    FCDR_pth=raw_input("Where is FCDR basedir? (e.g./DSNNAS/Repro/mviri/level1/MFG15_FCDR_V1/) ")  
    
    from mviri.mviri_l10 import read_imag2tg as mviri_l10
    from mviri.mviri_l15 import read_rect2lp as mviri_l15
    import numpy as np
    
    import multiprocessing
    import psutil
    
    L15_path="/DSNNAS/Repro/mviri/level1/HR-MFG15/data/"
  
    sats=["MET7","MET6","MET5","MET4","MET3","MET2"]
    #sats=["MET5","MET4","MET3","MET2"]
    #proc = 6 #psutil.cpu_count()
    #pool = multiprocessing.Pool(processes=proc)
    #pool.map(enter,sats)
    #multiprocessing somehow fails so:
    for sat in sats:
      enter(sat,FCDR_pth,rel,anomalycheck)
    exit()



    
