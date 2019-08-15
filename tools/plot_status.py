import pandas as pd
import numpy as np
import os,sys
import matplotlib.pylab as plt
import time

cont="y"

sys.path.insert(0, '../src/')
import ricals_tools as to

while cont=="y":
  rel=raw_input("which release (e.g. 3.0)?")
  sat=raw_input("which sat (e.g. MET5 or all)?")

  logs=to.logger(rel)

  df=logs.query2pandas()

  #print df.columns

  df=df.sort_values('scan_ts')
  if not sat=="all":
    df=df[df.platform==sat]
  print "\ncurrent timestamp: ",time.time()
  #print "\nsuccess difference:",np.sum(df.status_easy)
  print "easy all entries:  ",df.shape[0]
  print "easy updated:      ",df[df.status_easy==2].shape[0]
  print "easy successful:   ",df[df.status_easy==1].shape[0]
  print "easy pending:      ",df[df.status_easy==0].shape[0]
  print "easy errors:       ",df[df.status_easy==-1].shape[0]

  plt.plot(df.scan_ts,df.status_easy,'+',label="easy")
  plt.plot(df.scan_ts,df.status_full,'x',label="full")
  plt.yticks([-1,0,1,2], ["error","pending","successful","format updated"])
  plt.legend()
  plt.show()

  q=raw_input("\nprint error descriptions? (y/N)")
  if (q=="y") or (q=="Y"):
    dft=df[df.status_easy==-1]
    try:
      errors=dft.groupby(dft.problem_easy.tolist(),as_index=False).size()
      #errors=errors.drop_duplicates
      print errors
    except ValueError:
      print "\nNo ERRORS happened :)!\n"

  q=raw_input("\nprint pending descriptions? (y/N)")
  if (q=="y") or (q=="Y"):
    dft=df[df.status_easy==0]
    try:
      errors=dft.groupby(dft.problem_easy.tolist(),as_index=False).size()
      #errors=errors.drop_duplicates
      print errors
    except ValueError:
      print "\nNo ERRORS happened :)!\n"

  q=raw_input("\nprint successful logs? (y/N)")
  if (q=="y") or (q=="Y"):
    dft=df[df.status_easy==1]
    errors=dft.groupby(dft.problem_easy.tolist(),as_index=False).size()
    #errors=errors.drop_duplicates
    print errors
  
  q=raw_input("\nprint format update logs? (y/N)")
  if (q=="y") or (q=="Y"):
    dft=df[df.status_easy==2]
    errors=dft.groupby(dft.problem_easy.tolist(),as_index=False).size()
    #errors=errors.drop_duplicates
    print errors
  
  cont=raw_input("\nContinue? (y/n)")
