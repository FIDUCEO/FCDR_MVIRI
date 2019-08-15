#source /tcenas/home/frankr/git/RICalPy/environment.sh
#cd /tcenas/home/frankr/git/RICalPy/tools
import traceback
import pandas as pd
import numpy as np
import os,sys
import matplotlib.pylab as plt
import time
import glob
import datetime

import netCDF4 as nc
import xarray as xr
import datetime

from fiduceo.fcdr.writer.fcdr_writer import FCDRWriter
from fiduceo.common.writer.default_data import DefaultData
from fiduceo.common.writer.templates.templateutil import TemplateUtil as tu
from xarray import Variable

import shutil
import os

sys.path.insert(0, '../src/')
import ricals_tools as to

from mviri.mviri_l10 import read_imag2tg as r0
from mviri.mviri_l15 import read_rect2lp as rd


global early
global rlppath
global itgpath
global sqlcollect
global commandcollect

early    = ['M2','M3','P2']
rlppath="/DSNNAS/Repro/mviri/level1/HR-MFG15/data/"
itgpath="/DSNNAS/Repro/mviri/level0/IMAG2TG/data/"

sqlcollect="/DSNNAS/Repro_Temp/users/frankr/SATsqlcollect.txt"
commandcollect="/DSNNAS/Repro_Temp/users/frankr/SATcommandcollect.txt"


def fixfile(filename):
  '''
  problem: netCDF4/HDF5 problem indicates a corrupted file.
  solution: recreate the file or put on blacklist for deletion
  inputs:
  filename
  '''
  src=nc.Dataset(filename,"r")
  
  dst=nc.Dataset(filename+"b",'w')#,diskless=True,persist=True)
  
  dst.setncatts(src.__dict__)
  # copy dimensions
  for name, dimension in src.dimensions.items():
      dst.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))
  
  # copy all file data 
  ef=[]
  for name, variable in src.variables.items():
          x = dst.createVariable(name, variable.datatype, variable.dimensions,zlib=True,complevel=5)
          # copy variable attributes all at once via dictionary
          dst[name].setncatts(src[name].__dict__)
          try:
            dst[name][:] = src[name][:]
            ef.append(0)
          except IndexError:
            ef.append(1)
  
  src.close()
  dst.close()
  
  os.remove(filename)
  shutil.move(filename+"b",filename)
  if 1 in ef:
    return 1
  else:
    return 0
  
def fixHDF5limit(filename,var,to_add):
  '''
  problem: netCDF4/HDF5 only allows limited number of edits of attributes.
  solution: recreate the file to set log to u_zero_vis
  inputs:
  filename
  var: string with attribute name
  var: string that is to be added to attribute
  '''
  src=nc.Dataset(filename,"r")

  dst=nc.Dataset(filename+"b",'w')#,diskless=True,persist=True)
  dst.setncatts(src.__dict__)
  # copy dimensions
  for name, dimension in src.dimensions.items():
      dst.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))
  
  # copy all file data 
  for name, variable in src.variables.items():
          x = dst.createVariable(name, variable.datatype, variable.dimensions,zlib=True,complevel=5)
          # copy variable attributes all at once via dictionary
          dst[name].setncatts(src[name].__dict__)
          try:
            dst[name][:] = src[name][:]
          except IndexError:
            dst[name] = src[name]
  
  src.close()
  varval=dst.getncattr(var)
  dst.setncattr(var, varval+to_add)
  dst.close()
  os.remove(filename)
  shutil.move(filename+"b",filename)

def modDOI(filename,branch,base,sat,format_con,curr_datetime,logs):
  datetime_now=datetime.datetime.now()
  #check if attribute exists using netcdf
  try:
    ds = nc.Dataset(filename, "a", format="NETCDF4")
    ef=0
  except IOError,i:
      print i
      traceback.print_exc()
      print "...fix..."
      try:
        ef=fixfile(filename)
      except IOError:#if fixing failed, ef will be set to 1
        ef=1
        print "...fix failed"
      except IndexError:#if fixing failed, ef will be set to 1
        ef=1
        print "...fix failed"
      if ef==0:
        print "...try again"
        ds = nc.Dataset(filename, "a", format="NETCDF4")
  if ef==0:#if fixing failed, ef will be set to 1
    try:
      currdoistr=ds.getncattr("doi")
      format_con=format_con+1 #count as successful
      err=0
    except AttributeError:
      err=1
    #if not exists create
    if err==1 or "no doi" in currdoistr:
      SSP=float(base[base.index(sat)+5:base.index(sat)+9])
      try:
        doistr={0.0 :"10.15770/EUM_SEC_CLM_0009",\
                57.5:"10.15770/EUM_SEC_CLM_0012",\
                57.0:"10.15770/EUM_SEC_CLM_0012",\
                63.0:"10.15770/EUM_SEC_CLM_0013"}[SSP]
      except KeyError:
        doistr="no doi for this sub satellite point"
        if branch=="EASY":
          err=logs.append(sat,curr_datetime,"problem_easy",";no doi for SSP")
        if branch=="FULL":
          err=logs.append(sat,curr_datetime,"problem_full",";no doi for SSP")
      #using netcdf
      ds.setncattr("doi",doistr)
      history=ds.getncattr("history")
      #annotate this in logging db and file:
      ds.setncattr("history",history+";added doi: "+datetime_now.strftime("%Y/%m/%d"))
      try:
        ds.close()
      except RuntimeError,i:
        print i
        traceback.print_exc()
      format_con=format_con+1
      if branch=="EASY":
        err=logs.append(sat,curr_datetime,"problem_easy",";added doi")
      if branch=="FULL":
        err=logs.append(sat,curr_datetime,"problem_full",";added doi")
  return format_con,ef

def modAUTHORS(filename,branch,base,sat,format_con,curr_datetime,logs):
  datetime_now=datetime.datetime.now()
  ds = nc.Dataset(filename, "a", format="NETCDF4")
  history=ds.getncattr("history")
  try:
    if "updated authorship" in history:
      format_con=format_con+1
    else:
      refs="Product User Guide Document reference: EUM/USC/DOC/17/906121\n"+\
            "Methods:\n"+\
            "Algorithm Theorethical Basis Document Reference: EUM/OPS/DOC/18/990143\n"+\
            "Ruethrich, F.; John, V.O.; Roebeling, R.A.; Quast, R.; Govaerts, Y.; Woolliams, E.R.; Schulz, J. Climate Data Records from Meteosat First Generation Part III: Recalibration and Uncertainty Tracing of the Visible Channel on Meteosat-27 Using Reconstructed, Spectrally Changing Response Functions. Remote Sens. 2019, 11, 1165.\n"+\
            "Quast, R.; Giering, R.; Govaerts, Y.; Ruethrich, F.; Roebeling, R. Climate Data Records from Meteosat First Generation Part II: Retrieval of the In-Flight Visible Spectral Response. Remote Sens. 2019, 11, 480. \n"+\
            "Govaerts, Y.M.; Ruethrich, F.; John, V.O.; Quast, R. Climate Data Records from Meteosat First Generation Part I: Simulation of Accurate Top-of-Atmosphere Spectral Radiance over Pseudo-Invariant Calibration Sites for the Retrieval of the In-Flight Visible Spectral Response. Remote Sens. 2018, 10, 1959. \n"+\
            "John, V.O.; Tabata, T.; Ruethrich, F.; Roebeling, R.; Hewison, T.; Stoeckli, R.; Schulz, J. On the Methods for Recalibrating Geostationary Longwave Channels Using Polar Orbiting Infrared Sounders. Remote Sens. 2019, 11, 1171. \n"+\
            "Original Data:\n"+\
            "Level 1.5 Format and Metadata Document Reference: EUM/TSS/ICD/14/778737\n"+\
            "Level 1.0 Format and Metadata Document Reference: EUM/OPS-MTP/MAN/16/854401\n"+\
            "Technical Note Orbit Coordinates Document Reference: EUM/OPS/DOC/18/1000912\n"+\
            "Ruethrich, F.; John, V.O.; Roebeling, R.; Wagner, S.; Viticchie, B.; Hewison, T.; Govaerts, Y.; Quast, R.; Giering, R.; Schulz, J. A Fundamental Climate Data Record that accounts for Meteosat First Generation Visible Band Spectral Response Issues. In Proceedings of the 2016 EUMETSAT Meteorological Satellite Conference, Darmstadt, Germany, 2630 September 2016."
      authorshipdict={  "authors" :"EUMETSAT",\
                        "email"   :"ops@eumetsat.int",\
                        "url"     :"www.eumetsat.int for the full dataset and www.fiduceo.eu for the VIS channel",\
                        "references":refs}
      errcount=0
      for item in authorshipdict:
        try:
          insertstr=authorshipdict[item]
          ds.setncattr(item,insertstr)
        except:
          errorcount=errcount+1
      if errcount<len(authorshipdict):#at least one of authorship items written
        #annotate this in logging db and file:
        ds.setncattr("history",history+";updated authorship : "+datetime_now.strftime("%Y/%m/%d"))
      else:
        ds.close()
        raise RuntimeError
      try:
        ds.close()
      except RuntimeError,i:
          print i
          traceback.print_exc()
          raise
      format_con=format_con+1
      if branch=="EASY":
        err=logs.append(sat,curr_datetime,"problem_easy",";updated authorship")
      if branch=="FULL":
        err=logs.append(sat,curr_datetime,"problem_full",";updated authorship")
  except RuntimeError:
    print("RuntimeError; probably: NetCDF: HDF error when trying to close dataset [modAUTHOR]")
    if branch=="EASY":
      err=logs.append(sat,curr_datetime,"problem_easy",";corrupted: could not update authorship")
    if branch=="FULL":
      err=logs.append(sat,curr_datetime,"problem_full",";corrupted: could not update authorship")
  return format_con

def modLICENSE(filename,branch,base,sat,format_con,curr_datetime,logs):
  datetime_now=datetime.datetime.now()
  ds = nc.Dataset(filename, "a", format="NETCDF4")
  history=ds.getncattr("history")
  try:
    if "updated license" in history:
      format_con=format_con+1
    else:
      licensevis="Content in this file that is related to the visible channel is released for use under CC-BY licence (https://creativecommons.org/licenses/by/4.0/) and was developed in the EC FIDUCEO project \"Fidelity and Uncertainty in Climate Data Records from Earth Observations\". Grant Agreement: 638822."
      licenseir ="Content in this file that is related to the infrared and water vapour channel is released for use according to the EUMETSAT data policy. Access to this product is granted to all users without charge and without conditions on use if a licence agreement has been signed. For the full EUMETSAT data policy, please refer to the Product User Guide and the corresponding EUMETSAT webpage: https://www.eumetsat.int/website/home/AboutUs/WhoWeAre/LegalFramework/DataPolicy/index.html"
      license=licensevis+"\n"+licenseir
      licensedict={"licence" :license}
      errcount=0
      for item in licensedict:
        try:
          insertstr=licensedict[item]
          ds.setncattr(item,insertstr)
        except:
          errorcount=errcount+1
      if errcount<len(licensedict):#at least one of authorship items written
        #annotate this in logging db and file:
        ds.setncattr("history",history+";updated license: "+datetime_now.strftime("%Y/%m/%d"))
      else:
        ds.close()
        raise RuntimeError
      try:
        ds.close()
      except RuntimeError,i:
          print i
          traceback.print_exc()
          raise
      format_con=format_con+1
      if branch=="EASY":
        err=logs.append(sat,curr_datetime,"problem_easy",";updated license")
      if branch=="FULL":
        err=logs.append(sat,curr_datetime,"problem_full",";updated license")
  except RuntimeError:
    print("RuntimeError; probably: NetCDF: HDF error when trying to close dataset [modLICENSE]")
    if branch=="EASY":
      err=logs.append(sat,curr_datetime,"problem_easy",";corrupted: could not update license")
    if branch=="FULL":
      err=logs.append(sat,curr_datetime,"problem_full",";corrupted: could not update license")
  return format_con

def modSSP(filename,branch,base,sat,format_con,curr_datetime,logs):
  datetime_now=datetime.datetime.now()
  if branch=="EASY":
    ds = nc.Dataset(filename, "a", format="NETCDF4")
    history=ds.getncattr("history")
    try:
      if "corrected sub_satellite_latitude_* and sub_satellite_longitude_*" in history:
        format_con=format_con+1
      else:
        try:
          #tmplats=ds.variables["sub_satellite_longitude_start"][:]
          #tmplate=ds.variables["sub_satellite_longitude_end"][:]
          #tmplons=ds.variables["sub_satellite_latitude_start"][:]
          #tmplone=ds.variables["sub_satellite_latitude_end"][:]
          r2lp_datetime=curr_datetime+datetime.timedelta(minutes=30)
          rect2lpf=ds.getncattr("RECT2LP_file_name")
          rect2lp=rlppath+sat+"/"+str(r2lp_datetime.year)+"/"+str(r2lp_datetime.month).zfill(2)+"/"+rect2lpf
          imag2tg=r0.l15_to_L10_name(rect2lp,str(r2lp_datetime.year)+str(r2lp_datetime.month).zfill(2)+str(r2lp_datetime.day).zfill(2))
          with open(rect2lp, "rb") as fi:
            fi, head = rd.header(fi, 1, 3,t=0)
          with open(imag2tg, "rb") as fi:
            try:
              fi,head0   = r0.header(fi)
            except:
              head0='none'
              print "no L0 file"
          nomlon = round(np.rad2deg(getattr(head[1],'dsNominalSubSatelliteLongitude')[0]),2)
          tmplats, tmplons, tmplate, tmplone = to.calcSSP(head,head0,nomlon,early,logs,r2lp_datetime)
          
          ds.variables["sub_satellite_longitude_start"][:]=tmplons
          ds.variables["sub_satellite_longitude_end"][:]=tmplone
          ds.variables["sub_satellite_latitude_start"][:]=tmplats
          ds.variables["sub_satellite_latitude_end"][:]=tmplate
        except KeyError,i:
          print("ERROR: Var somehow not in FCDR file: sub_satellite_l* -->needs to be regenerated")
          print i
          with open( sqlcollect.replace("SAT",sat+branch) , 'a') as sql:
            with open( commandcollect.replace("SAT",sat+branch) , 'a') as com:
              sql.write("DELETE FROM logtable WHERE idn=%s%s;\n" %(sat[3:],r2lp_datetime.strftime("%Y%m%d%H%M")))
              com.write("/tcenas/home/frankr/git/RICalPy/pbsify/shellwrap.sh -s %s -e %s -m %s -r 3.1 -v 1801 -a /DSNNAS/Repro/mviri/level1/MFG15_FCDR_V1/ -c 'first release - use with caution'\n" %(curr_datetime.strftime("%Y%m%d%H%M"), r2lp_datetime.strftime("%Y%m%d%H%M"),sat))
        except IndexError,i:
          print i
          ds.variables["sub_satellite_longitude_start"][:]=ds.variables["sub_satellite_longitude_start"].getncattr("_FillValue")
          ds.variables["sub_satellite_longitude_end"][:]=ds.variables["sub_satellite_longitude_end"].getncattr("_FillValue")
          ds.variables["sub_satellite_latitude_start"][:]=ds.variables["sub_satellite_latitude_start"].getncattr("_FillValue")
          ds.variables["sub_satellite_latitude_end"][:]=ds.variables["sub_satellite_latitude_end"].getncattr("_FillValue")
        #update logger and history
        try:
          ds.setncattr("history", history+";corrected sub_satellite_latitude_* and sub_satellite_longitude_*: "+datetime_now.strftime("%Y/%m/%d"))
          ds.close()
        except RuntimeError,i:
          print i
          traceback.print_exc()
          to_add=";corrected sub_satellite_latitude_* and sub_satellite_longitude_*: "+datetime_now.strftime("%Y/%m/%d")
          fixHDF5limit(filename,"history",to_add)
        format_con=format_con+1
        if branch=="EASY":
          err=logs.append(sat,curr_datetime,"problem_easy",";corrected SSP")
        if branch=="FULL":
          err=logs.append(sat,curr_datetime,"problem_full",";corrected SSP")
    except RuntimeError,i:
      print i
      traceback.print_exc()
      #print ("RuntimeError; probably: NetCDF: HDF error when trying to close dataset [modSSP]")
      if branch=="EASY":
        err=logs.append(sat,curr_datetime,"problem_easy",";corrupted: could not correct SSP")
      if branch=="FULL":
        err=logs.append(sat,curr_datetime,"problem_full",";corrupted: could not correct SSP")
  return format_con

def modUNITS(filename,branch,base,sat,format_con,curr_datetime,logs):
  datetime_now=datetime.datetime.now()
  variabdict ={"count_ir":"count",\
                "count_wv":"count",\
                "solar_irradiance_vis":"W*m^-2",\
                "u_solar_irradiance_vis":"W*m^-2",\
                "a_ir":"mW*m^-2*sr^-1*cm^-1",\
                "b_ir":"mW*m^-2*sr^-1*cm^-1*count^-1",\
                "u_a_ir":"mW*m^-2*sr^-1*cm^-1",\
                "u_b_ir":"mW*m^-2*sr^-1*cm^-1*count^-1",\
                "a_wv":"mW*m^-2*sr^-1*cm^-1",\
                "b_wv":"mW*m^-2*sr^-1*cm^-1*count^-1",\
                "u_a_wv":"mW*m^-2*sr^-1*cm^-1",\
                "u_b_wv":"mW*m^-2*sr^-1*cm^-1*count^-1",\
                "count_vis":"count",\
                "a0_vis":"W*m^-2*sr^-1*count^-1",\
                "a1_vis":"W*m^-2*sr^-1*count^-1*year^-1",\
                "a2_vis":"W*m^-2*sr^-1*count^-1*year^-2",\
                "u_a0_vis":"W*m^-2*sr^-1*count^-1",\
                "u_a1_vis":"W*m^-2*sr^-1*count^-1*year^-1",\
                "u_a2_vis":"W*m^-2*sr^-1*count^-1*year^-2",\
                "u_zero_vis":"W*m^-2*sr^-1*count^-1",\
                "covariance_a_vis":"W*m^-2*sr^-1*count^-1"\
                }
  
  ds = nc.Dataset(filename, "a", format="NETCDF4")
  history=ds.getncattr("history")
  try:
    if "updated units" in history:
      format_con=format_con+1
    else:
      errcount=0
      for v in variabdict:
        try:
          ds.variables[v].units=variabdict[v]
        except KeyError:
          #var not in file; expected for easy
          errcount=errcount+1
      #update logger and history
      if errcount<len(variabdict):#at least one unit written
        ds.setncattr("history", history+";updated units: "+datetime_now.strftime("%Y/%m/%d"))
      else:
        raise RuntimeError
      try:
        ds.close()
      except RuntimeError,i:
          print i
          traceback.print_exc()
          raise
      format_con=format_con+1
      if branch=="EASY":
        err=logs.append(sat,curr_datetime,"problem_easy",";units updated")
      if branch=="FULL":
        err=logs.append(sat,curr_datetime,"problem_full",";units updated")
  except RuntimeError:
    print ("RuntimeError; probably: NetCDF: HDF error when trying to close dataset [modUNITS]")
    if branch=="EASY":
      err=logs.append(sat,curr_datetime,"problem_easy",";corrupted: could not update units")
    if branch=="FULL":
      err=logs.append(sat,curr_datetime,"problem_full",";corrupted: could not update units")
  return format_con

def modANC(filename,branch,base,sat,format_con,curr_datetime,logs):
  datetime_now=datetime.datetime.now()
  #for comment
  commentdict ={"count_ir":"convert to radiance with radiance=a_ir+count_ir*b_ir",\
                "count_wv":"convert to radiance with radiance=a_wv+count_wv*b_wv",\
                "count_vis":"convert to radiance with radiance=(count_vis-mean_count_space_vis)*(a0_vis+a1_vis*years_since_launch+a2_vis*years_since_launch^2)",\
                "time":      "acquisition time in IR/WV grid; can be used for VIS channel by linear interpolation in x-direction and by duplicating each line in y-direction",\
                "time_ir_wv":"acquisition time in IR/WV grid; can be used for VIS channel by linear interpolation in x-direction and by duplicating each line in y-direction",\
                "satellite_azimuth_angle":"tie-point grid contains every 10th entry of full VIS grid, starting at index [0,0]. We recommend cubic spline interpolation to reconstruct full grid.",\
                "satellite_zenith_angle":"tie-point grid contains every 10th entry of full VIS grid, starting at index [0,0]. We recommend cubic spline interpolation to reconstruct full grid.",\
                "solar_azimuth_angle":"tie-point grid contains every 10th entry of full VIS grid, starting at index [0,0]. We recommend cubic spline interpolation to reconstruct full grid.",\
                "solar_zenith_angle":"tie-point grid contains every 10th entry of full VIS grid, starting at index [0,0]. We recommend cubic spline interpolation to reconstruct full grid.",\
                }
  #for ancillary_variables
  relateddict ={"count_ir":"a_ir b_ir u_a_ir u_b_ir",\
                "count_wv":"a_wv b_wv u_a_wv u_b_wv",\
                "count_vis":"mean_count_space_vis a0_vis a1_vis a2_vis u_a0_vis u_a1_vis u_a2_vis covariance_a_vis years_since_launch",\
                "time":"",\
                "time_ir_wv":"",\
                "satellite_azimuth_angle":"satellite_zenith_angle",\
                "satellite_zenith_angle":"satellite_azimuth_angle",\
                "solar_azimuth_angle":"solar_zenith_angle",\
                "solar_zenith_angle":"solar_azimuth_angle",\
                }
  ds = nc.Dataset(filename, "a", format="NETCDF4")
  history=ds.getncattr("history")
  try:
    if "updated comments and ancillary" in history:
      format_con=format_con+1
    else:
      errcount=0
      for e in commentdict:
        try:
          ds.variables[e].comment=commentdict[e]
          ds.variables[e].ancillary_variables=relateddict[e]
        except KeyError:
          #var not in file; likely to happen for easy
          errcount=errcount+1
      #update logger and history
      if errcount<len(commentdict):#at least one expression written
        ds.setncattr("history", history+";updated comments and ancillary: "+datetime_now.strftime("%Y/%m/%d"))
      else:
        raise RuntimeError
      try:
        ds.close()
      except RuntimeError,i:
          print i
          traceback.print_exc()
          raise
      format_con=format_con+1
      if branch=="EASY":
        err=logs.append(sat,curr_datetime,"problem_easy",";comments and ancillary updated updated")
      if branch=="FULL":
        err=logs.append(sat,curr_datetime,"problem_full",";comments and ancillary updated updated")
  except RuntimeError:
    print ("RuntimeError; probably: NetCDF: HDF error when trying to close dataset [modANC]")
    if branch=="EASY":
      err=logs.append(sat,curr_datetime,"problem_easy",";corrupted: could not update comments")
    if branch=="FULL":
      err=logs.append(sat,curr_datetime,"problem_full",";corrupted: could not update comments")
  return format_con

def modIRWV(filename,branch,base,sat,format_con,curr_datetime,dc,time,logs):
  ef=0
  datetime_now=datetime.datetime.now()
  #get curr_idx:
  curr_idx=np.abs(time-curr_datetime).argmin(axis=-1)
  ds = nc.Dataset(filename, "a", format="NETCDF4")
  history=ds.getncattr("history")
  if "updated IR/WV" in history:
    format_con=format_con+1
  else:
    try:
      for irVar in ["a_ir","a_wv","b_ir","b_wv","u_a_ir","u_a_wv","u_b_ir","u_b_wv"]:
        try:
          val=dc.variables[irVar][curr_idx]
        except KeyError:
          try:
            val=ds.variables[irVar].getncattr("_FillValue")#uncertainties no longer in calib file from viju
          except KeyError:#var also not in FCDR file to get fillvalue
            ef=1
        try:
          ds.variables[irVar][:]=val
        except KeyError:
          print("ERROR: Var somehow not in FCDR file: "+irVar)
          ef=1
        except IndexError:
          ds.variables[irVar][:]=ds.variables[irVar].getncattr("_FillValue")
      #update logger and history
      ds.setncattr("history", history+";updated IR/WV: "+datetime_now.strftime("%Y/%m/%d"))
      try:
        ds.close()
      except RuntimeError,i:
          print i
          traceback.print_exc()
          raise
      format_con=format_con+1
      if branch=="EASY":
          err=logs.append(sat,curr_datetime,"problem_easy",";updated IR/WV")
      if branch=="FULL":
          err=logs.append(sat,curr_datetime,"problem_full",";updated IR/WV")
    except RuntimeError:
      print ("RuntimeError; probably: NetCDF: HDF error when trying to close dataset [modIRWV]")
      if branch=="EASY":
        err=logs.append(sat,curr_datetime,"problem_easy",";corrupted: could not update IR/WV")
      if branch=="FULL":
        err=logs.append(sat,curr_datetime,"problem_full",";corrupted: could not update IR/WV")
  return format_con,ef

def modNAMES(filename,branch,base,sat,format_con,curr_datetime,logs):
  datetime_now=datetime.datetime.now()
  variabdict ={"time":"time_ir_wv",\
                }
  ds = nc.Dataset(filename, "a", format="NETCDF4")
  history=ds.getncattr("history")
  try:
    if "updated names" in history:
      format_con=format_con+1
    else:
      errcount=0
      for v in variabdict:
        try:
          ds.renameVariable(v,variabdict[v])
        except KeyError:
          #var not in file; expected for easy
          errcount=errcount+1
      #update logger and history
      if errcount<len(variabdict):#at least one unit written
        ds.setncattr("history", history+";updated names: "+datetime_now.strftime("%Y/%m/%d"))
      else:
        raise RuntimeError
      try:
        ds.close()
      except RuntimeError,i:
          print i
          traceback.print_exc()
          raise
      format_con=format_con+1
      if branch=="EASY":
        err=logs.append(sat,curr_datetime,"problem_easy",";names updated")
      if branch=="FULL":
        err=logs.append(sat,curr_datetime,"problem_full",";names updated")
  except RuntimeError:
    print ("RuntimeError; probably: NetCDF: HDF error when trying to close dataset [modNAMES]")
    if branch=="EASY":
      err=logs.append(sat,curr_datetime,"problem_easy",";corrupted: could not update names")
    if branch=="FULL":
      err=logs.append(sat,curr_datetime,"problem_full",";corrupted: could not update names")
  return format_con

def modFLAGS(filename,branch,base,sat,format_con,curr_datetime,logs):
  datetime_now=datetime.datetime.now()
  variabdict ={"data_quality_bitmask":"",\
               "quality_pixel_bitmask":"",\
              }
  ds = nc.Dataset(filename, "a", format="NETCDF4")
  history=ds.getncattr("history")
  try:
    if "updated flag_masks" in history:
      format_con=format_con+1
    else:
      errcount=0
      for v in variabdict:
        try:
          masks    = ds.variables[v].flag_masks
          try:
            ds.variables[v].flag_masks=np.array(masks.split(",")).astype(int)
          except AttributeError: 
            ds.variables[v].flag_masks=masks.astype(int)
          #meanings = ds.variables[v].flag_meanings
          #ds.variables[v].flag_meanings=np.array(meanings.split(" ")).astype(str)
        except KeyError:
          #var not in file; expected for easy
          errcount = errcount+1
      #update logger and history
      if errcount<len(variabdict):#at least one unit written
        ds.setncattr("history", history+";updated flag_masks: "+datetime_now.strftime("%Y/%m/%d"))
      else:
        raise RuntimeError
      try:
        ds.close()
      except RuntimeError,i:
          print i
          traceback.print_exc()
          raise
      format_con=format_con+1
      if branch=="EASY":
        err=logs.append(sat,curr_datetime,"problem_easy",";flag_masks updated")
      if branch=="FULL":
        err=logs.append(sat,curr_datetime,"problem_full",";flag_masks updated")
  except RuntimeError:
    print ("RuntimeError; probably: NetCDF: HDF error when trying to close dataset [modFLAGS]")
    if branch=="EASY":
      err=logs.append(sat,curr_datetime,"problem_easy",";corrupted: could not update flag_masks")
    if branch=="FULL":
      err=logs.append(sat,curr_datetime,"problem_full",";corrupted: could not update flag_masks")
  return format_con

def modLONG(filename,branch,base,sat,format_con,curr_datetime,logs):
  datetime_now=datetime.datetime.now()
  variabdict ={"count_wv":"Water vapour image counts",\
               "count_vis":"Visible image counts",\
              }
  ds = nc.Dataset(filename, "a", format="NETCDF4")
  history=ds.getncattr("history")
  try:
    if "updated longnames" in history:
      format_con=format_con+1
    else:
      errcount=0
      for v in variabdict:
        try:
          ds.variables[v].long_name=variabdict[v]
        except KeyError:
          #var not in file; expected for easy
          errcount=errcount+1
      #update logger and history
      if errcount<len(variabdict):#at least one unit written
        ds.setncattr("history", history+";updated longnames: "+datetime_now.strftime("%Y/%m/%d"))
      else:
        raise RuntimeError
      try:
        ds.close()
      except RuntimeError,i:
          print i
          traceback.print_exc()
          raise
      format_con=format_con+1
      if branch=="EASY":
        err=logs.append(sat,curr_datetime,"problem_easy",";longnames updated")
      if branch=="FULL":
        err=logs.append(sat,curr_datetime,"problem_full",";longnames updated")
  except RuntimeError:
    print ("RuntimeError; probably: NetCDF: HDF error when trying to close dataset [modLONG]")
    if branch=="EASY":
      err=logs.append(sat,curr_datetime,"problem_easy",";corrupted: could not update longnames")
    if branch=="FULL":
      err=logs.append(sat,curr_datetime,"problem_full",";corrupted: could not update longnames")
  return format_con

def iterate(sat,rel,branch):
  '''
  this does the injection for one satellite for one branch (FULL/EASY)
  '''
  
  open( sqlcollect.replace("SAT",sat+branch) , 'w+')
  open( commandcollect.replace("SAT",sat+branch) , 'w+')
  
  datetime_now=datetime.datetime.now()
  
  #prepare IR calibratrion
  calfolder="/tcenas/home/vijuj/work/satcal/mike2check/"
  calfile=calfolder+sat+"_calibration_coefficients.nc"
  try:
    dc=nc.Dataset(calfile)
    time=dc.variables["time"][:]
    time=np.array([datetime.datetime(1970,1,1,0,0)+datetime.timedelta(seconds=t) for t in time])
  except:
    print "cannot open IR calibration file"
    exit(1)
  #open logger db
  logs=to.logger(rel)
  
  #get filenames from db
  #df=logs.query4injection()
  #df=df.sort_values('scan_ts')
  #if not sat=="all":
    #df=df[df.platform==sat]
  #filenames=df.filename_easy
  
  #get filenames from iterating through all dirs:
  print("get filenames; takes a bit...")
  filenames=glob.glob("/DSNNAS/Repro/mviri/level1/MFG15_FCDR_V1/data/"+sat+"/????/??/*"+branch+"*.nc")
  print("...got 'em")
  for filename in sorted(filenames):
    print(filename)
    #format-controller - used to count successfull updates of the file
    #for modifying the logger database afterwards
    exN=0#expected number of modifications
    format_con=0#counter for successful modifications
    #get datetime from filename
    base=os.path.basename(filename)
    try:
      datestring=base[base.index(sat)+10:base.index(sat)+22]
      shift=1
      while "_" in datestring:
        datestring=base[base.index(sat)+10+shift:base.index(sat)+22+shift] #happens sometimes for MET5
        shift=shift+1
      curr_datetime=datetime.datetime.strptime(datestring,"%Y%m%d%H%M")
    except ValueError,i:
      print i
    process=True
    if sat=="MET7":
      print "starting MET7 only in 2000!!!"
      if curr_datetime<datetime.datetime(2000,4,26,0):
        process=False
    if sat=="MET5":
      print "starting MET5 only in 1995!!!"
      if curr_datetime<datetime.datetime(1995,7,25,0):
        process=False
    if process:
      try:
        #0. test/blacklist file
        try:
          src=nc.Dataset(filename,"r")
          src.close()
          ef=0
        except IOError:
          ef=1
          r2lp_datetime=curr_datetime+datetime.timedelta(minutes=30)
          with open( sqlcollect.replace("SAT",sat+branch) , 'a') as sql:
            with open( commandcollect.replace("SAT",sat+branch) , 'a') as com:
              sql.write("DELETE FROM logtable WHERE idn=%s%s;\n" %(sat[3:],r2lp_datetime.strftime("%Y%m%d%H%M")))
              com.write("/tcenas/home/frankr/git/RICalPy/pbsify/shellwrap.sh -s %s -e %s -m %s -r 3.1 -v 1801 -a /DSNNAS/Repro/mviri/level1/MFG15_FCDR_V1/ -c 'first release - use with caution'\n" %(curr_datetime.strftime("%Y%m%d%H%M"), r2lp_datetime.strftime("%Y%m%d%H%M"),sat))

        #1. Modify Global Attributes

        #1.1 DOI----------------------------------------
        if sat=="MET7": print "1"
        exN=exN+1
        if ef ==0:#ef is set to 1 if file is on blacklist and has to be recreated
          format_con,ef = modDOI(filename,branch,base,sat,format_con,curr_datetime,logs)
        #1.2 authorship----------------------------------------
        if sat=="MET7": print "2"
        exN=exN+1
        if ef ==0:#ef is set to 1 if file is on blacklist and has to be recreated
          format_con=modAUTHORS(filename,branch,base,sat,format_con,curr_datetime,logs)
        
        #1.3 license----------------------------------------
        if sat=="MET7": print "3"
        exN=exN+1
        if ef ==0:#ef is set to 1 if file is on blacklist and has to be recreated
          format_con=modLICENSE(filename,branch,base,sat,format_con,curr_datetime,logs)

        #2. Modify Variables and Variable-attributes

        #2.1 SSPlat/lon swapping--------------------------
        if sat=="MET7": print "4"
        exN=exN+1
        if ef ==0:#ef is set to 1 if file is on blacklist and has to be recreated
          format_con=modSSP(filename,branch,base,sat,format_con,curr_datetime,logs)

        #2.2 units----------------------------
        if sat=="MET7": print "5"
        exN=exN+1
        if ef ==0:#ef is set to 1 if file is on blacklist and has to be recreated
          format_con=modUNITS(filename,branch,base,sat,format_con,curr_datetime,logs)

        #2.3 comments and ancillary variables----
        if sat=="MET7": print "6"
        exN=exN+1
        if ef ==0:#ef is set to 1 if file is on blacklist and has to be recreated
          format_con=modANC(filename,branch,base,sat,format_con,curr_datetime,logs)

        #2.4 IR calibratrion------------------------------
        if sat=="MET7": print "7"
        exN=exN+1
        if ef ==0:#ef is set to 1 if file is on blacklist and has to be recreated
          format_con,ef=modIRWV(filename,branch,base,sat,format_con,curr_datetime,dc,time,logs)
        
        #2.5 Variable renaming----------------------------
        if sat=="MET7": print "8"
        exN=exN+1
        if ef ==0:#ef is set to 1 if file is on blacklist and has to be recreated
          format_con=modNAMES(filename,branch,base,sat,format_con,curr_datetime,logs)
        
        #2.6 Flags string----------------------------
        if sat=="MET7": print "9"
        exN=exN+1
        if ef ==0:#ef is set to 1 if file is on blacklist and has to be recreated
          format_con=modFLAGS(filename,branch,base,sat,format_con,curr_datetime,logs)
        
        #2.7 long names----------------------------
        if sat=="MET7": print "10"
        exN=exN+1
        if ef ==0:#ef is set to 1 if file is on blacklist and has to be recreated
          format_con=modLONG(filename,branch,base,sat,format_con,curr_datetime,logs)
          
        if ef==1:
          print "...blacklisting file"
          r2lp_datetime=curr_datetime+datetime.timedelta(minutes=30)
          with open( sqlcollect.replace("SAT",sat+branch) , 'a') as sql:
            with open( commandcollect.replace("SAT",sat+branch) , 'a') as com:
              sql.write("DELETE FROM logtable WHERE idn=%s%s;\n" %(sat[3:],r2lp_datetime.strftime("%Y%m%d%H%M")))
              com.write("/tcenas/home/frankr/git/RICalPy/pbsify/shellwrap.sh -s %s -e %s -m %s -r 3.1 -v 1801 -a /DSNNAS/Repro/mviri/level1/MFG15_FCDR_V1/ -c 'first release - use with caution'\n" %(curr_datetime.strftime("%Y%m%d%H%M"), r2lp_datetime.strftime("%Y%m%d%H%M"),sat))
          
      except RuntimeError:
        print("some issue with the file")
        exN=1
        format_con=0
      #update db-entry for current filename
      if format_con>=exN:#update with actual number of expected mods
        if branch=="EASY":
          err=logs.update(sat,curr_datetime,"status_easy",2)
        if branch=="FULL":
          err=logs.update(sat,curr_datetime,"status_full",2)
      else:
        if branch=="EASY":
          err=logs.update(sat,curr_datetime,"status_easy",1)
        if branch=="FULL":
          err=logs.update(sat,curr_datetime,"status_full",1)

if "-h" in sys.argv[1:]:
  print ("usage: python FCDR_injector.py <release (e.g. 3.1)> <sat (e.g. MET7)> <branch (e.g. EASY)>\n"\
          "best to run with:\n"\
          "qsub injector_wrapper.sh -F MET5 -l walltime=96:00:00  -l mem=2000mb"\
          "or:\n"\
          "cd /tcenas/home/frankr/git/job_runners\n"+\
          "./qsub_2_level_batched_job.sh workspaceinjetion /tcenas/home/frankr/git/RICalPy/tools/injector_batch.sh --limits frankr=2 -- request_memory=2000mb")
  exit()

rel   =sys.argv[1]
sat   =sys.argv[2]
branch=sys.argv[3]

if sat=="all":
  for sat in["MET2","MET3","MET4","MET5","MET6","MET7"]:
    open( sqlcollect.replace("SAT",sat+branch) , 'w').close()
    open( commandcollect.replace("SAT",sat+branch) , 'w').close()
    iterate(sat,rel,branch)
else:
  iterate(sat,rel,branch)
exit()
