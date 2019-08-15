# -*- coding utf-8 -*-
"""
this file contains functions for writing the fcdr to netcdf
"""


import sys
import os
import numpy as np
import netCDF4 as nc
import xarray as xr
import datetime
import time
import timeit
import csv
import ricals_tools as to
import calendar
import time
import matplotlib.pyplot as plt

import ricals_effects as ef

from fiduceo.fcdr.writer.fcdr_writer import FCDRWriter
from fiduceo.common.writer.default_data import DefaultData
from fiduceo.common.writer.templates.templateutil import TemplateUtil as tu

import cruncher as cr

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#writing-related
#-----------------------------------------------------------------------
def ncName(s_Year,s_Month,s_Day,s_Time,dset,nomlon,ProcVersion,FormVersion,satellite="MET7",radiometer="MVIRI"):

  sepField        = "_"
  sepWord         = "-"
  
  if dset in ["EASY","easy","FULL","full"]:
    s_timeend       = s_Year+s_Month+s_Day+s_Time
    datetime_object = datetime.datetime.strptime(s_timeend, '%Y%m%d%H%M')
    td              = datetime.timedelta(minutes=30)
    datetime_object = datetime_object -td
    s_timestart     = datetime_object.strftime('%Y%m%d%H%M')
    fname           = "FIDUCEO"+sepField+"FCDR"+sepField+"L15"+sepField+radiometer+sepField \
                      +satellite+sepWord+str(nomlon).zfill(4)+sepField+s_timestart+sepField \
                      +s_timeend+sepField+dset+sepField+"v"+str(ProcVersion)+sepField \
                      +"fv"+str(FormVersion)+".nc"
    
  elif dset in ["STATIC","static","LUT","lut","aux","AUX"]:
    fname           = "FIDUCEO"+sepField+"FCDR"+sepField+"L15"+sepField+radiometer+sepField \
                      +satellite+sepWord+str(nomlon).zfill(4)+sepField \
                      +dset+sepField+"v"+str(ProcVersion)+sepField \
                      +"fv"+str(FormVersion)+".nc"
                    
  elif dset =="": #Release 1
    s_timeend       = s_Year+s_Month+s_Day+s_Time
    datetime_object = datetime.datetime.strptime(s_timeend, '%Y%m%d%H%M')
    td              = datetime.timedelta(minutes=30)
    datetime_object = datetime_object -td
    s_timestart     = datetime_object.strftime('%Y%m%d%H%M')
    fname           = "FCDR"+sepField+"L15"+sepField+radiometer+sepField \
                      +satellite+sepWord+str(nomlon).zfill(4)+sepField+s_timestart+sepField \
                      +s_timeend+sepField \
                      +"v"+str(ProcVersion)+sepField \
                      +"fv"+str(FormVersion)+".nc"
  else:
    print "ERROR: FCDR not defined for dset = "+dset
    exit()
    
  return fname

def write_manager( mfggrp,filename,head1,head2,head3,lininfo,trail,telem,telem_descr, \
                   satellite, logs, debug=0):
  '''
  This function organizes the writing of the complex fullFCDR files
  '''

  # Reading the excel file which gives information on whether attribute or variable
  store_info = np.genfromtxt('../config/vars_for_ricals.csv', names=True,delimiter=';',dtype=None)

  if debug>=1:
    print np.shape( store_info )

  namelin = lininfo[0]._asdict().keys()

  print 'Writing Header 1'
  write_head2nc(mfggrp, head1, 'header1_', store_info)

  print 'Writing Header 2'
  write_head2nc(mfggrp, head2, 'header2_', store_info)

  print 'Writing Header 3'
  write_head2nc(mfggrp, head3, 'header3_', store_info)

  print 'Writing Trailer'
  write_head2nc(mfggrp, trail, 'trailer_', store_info)

  print 'Writing Line Info'
  write_lineinfo2nc( mfggrp, lininfo, namelin )

  print 'Writing Telemetry'
  write_telem2nc( mfggrp, telem, telem_descr )
  
  return 0


def write_calfile2nc( mfggrp, data, names, unit, lname):
  """
  Fills the netCDF file with the content from Viju's IR/WV calibration files.
  Note: this has to be predefined in the netcdf template
  """
  for i,v,u,l in zip(names, data, unit, lname):
    dtype = 'f4'
    full_name  = str(i)
    long_name  = l
    scales     = 1
    size       = 1
    try:
      mfggrp.variables[full_name].data=v
    except KeyError:
      mfggrp.variables[full_name.replace("sigma","u")].data=v
  return

def write_head2nc( mfggrp, data, name, store_info ):
  """
  Writes and fills the netCDF file with the content of RECT2LP headers.
  Note: Does not have to be predefined in the netcdf template
  """
  for i,v in zip(data.keys(), data.values()):
        #print i
        if (not "Spare" in i) and (not "spare" in i):
            if i[:2] == "sh":
                full_name = name + i[2:]
                dtype = np.uint16
            if i[:2] == "fl":
                full_name = name + i[2:]
                dtype = np.float32
            if i[:2] == "ds":
                full_name = name + i[2:]
                dtype = np.float64
            if i[:1] == "i":
                full_name = name + i[1:]
                dtype = np.uint32
            if i[:2] == "sz":
                full_name = name + i[2:]
                dtype ='S1'
                v = ''.join(v)
            if i[:2] == "by":
                full_name = name + i[2:]
                dtype = np.ubyte

            j=store_info['name'].tolist().index(full_name)

            s_i=store_info[j]

            long_name  =  s_i[1]
            scales=float(s_i[10])
            
            # Have to add certain variables as global attributes. 
            # Also repeated variables in headers and trailers must kept only once
            
            #Some preparations
            if scales!=1:
                if v.isalnum():
                    try:
                      v=float(v)*scales
                    except:
                      v=np.nan
                    dtype=np.float32

            #STORING
            if s_i[8] == 1:#store current item? 1=yes
                if s_i[9] == 1: #how to store? 1=variable 0=global attribute
            #Store dimensions
                    fillval=DefaultData.get_default_fill_value(dtype)
                    if "DeformationMatrix"  in i:
                        v = np.array(v).reshape((26, 26))
                        default_array = DefaultData.create_default_array(26, 26, dtype)
                        variable = Variable(["y_"+"26", "x_"+"26"], default_array)
                    else:
                        default_array = DefaultData.create_default_vector(np.size(v), dtype)
                        variable = Variable(["y_"+str(np.size(v))], default_array)
                        
                    tu.add_fill_value(variable, fillval)
                    variable.attrs["long_name"] = long_name
                    mfggrp[full_name] = variable
                    #Store data
                    mfggrp.variables[full_name].data   = v
                elif s_i[9] == 0:
                    mfggrp.attrs[full_name]=v

def write_lineinfo2nc( mfggrp, data, name ):
  '''
  This function writes lineinfos to an open netcdf file
  Note: Does not have to be predefined in the netcdf template
  HISTORY:
  -originally designed to treat data as containing lininfo from image() reading function
  -XX.XX.2017 modified to handle xarray netcdf4 object
  -27.12.2017 modified to treat data as containing lininfo from image_3() reading function
  '''
  temp = 'lineinfo_'
  for v in name:
      linvar = []
      for i in data:
          linvar.append( getattr(i, v) )
      linvar=np.array(linvar)
      if (not "Spare" in v) and (not "spare" in v):
          if v[:2] == "sh":
              dtype=np.uint16
              fillval=DefaultData.get_default_fill_value(dtype)
              full_name=v[2:]
              default_array = DefaultData.create_default_vector(np.size(linvar), dtype)
              variable = Variable(["y_" + str(np.size(linvar))], default_array)
              values=linvar
          if v[:2] == "fl":
              dtype=np.float32
              fillval=DefaultData.get_default_fill_value(dtype)
              if "TimesArray" in v:
                  full_name=v[1:]
                  default_array = DefaultData.create_default_array(np.shape(linvar)[1],np.shape(linvar)[0], dtype)
                  variable = Variable(["y_" + str(np.shape(linvar)[0]),"x_"+str(np.shape(linvar)[1])], default_array)
                  values=linvar
              else:
                  full_name=v[2:]
                  default_array = DefaultData.create_default_vector(np.size(linvar), dtype)
                  variable = Variable(["y_" + str(np.size(linvar))], default_array)
                  values=linvar
              
          if v[:2] == "ds":
              dtype=np.float64
              fillval=DefaultData.get_default_fill_value(dtype)
              full_name=v[2:]
              if "FitCoef" in v:
                  default_array = DefaultData.create_default_array(np.shape(linvar)[1],np.shape(linvar)[0], dtype)
                  variable = Variable(["y_" + str(np.shape(linvar)[0]),"x_" + str(np.shape(linvar)[1])], default_array)
                  values=linvar
              else:
                  default_array = DefaultData.create_default_vector(np.size(linvar), dtype)
                  variable = Variable(["y_" + str(np.size(linvar))], default_array)
                  values=linvar
          if v[:1] == "i": 
              dtype=np.int32
              fillval=DefaultData.get_default_fill_value(dtype)
              full_name=v[1:]
              default_array = DefaultData.create_default_vector(np.size(linvar), dtype)
              variable = Variable(["y_"+ str(np.size(linvar))], default_array)
              values=linvar
          tu.add_fill_value(variable, fillval)
          mfggrp[temp + full_name] = variable
          #Store data
          mfggrp.variables[temp + full_name].data   = values

def write_telem2nc( mfggrp, data, descr ):
    '''
    This function writes extracted telemetry to an open netcdf file
    Note: Does not have to be predefined in the netcdf template
    '''
    temp = 'telem'
    names = data.dtype.names
    badnames = ['index', 'jday', 'HH', 'MM', 'SS', 'QA', 'missing', 'Frame']
    for name in names:
      if not name in badnames:
        try:
          telem_des     =  descr[name]
          valid_telem   =  data[name][(data["QA"]==0)&(data["index"]!=0)]
          mean_telem    =  np.mean(valid_telem)  
          median_telem  =  np.median(valid_telem)  
          std_telem     =  np.std(valid_telem)  
          num_telem     =  len(valid_telem)
          telem2wr = [mean_telem, median_telem, std_telem, num_telem]
        except:
          telem2wr = [np.nan,np.nan,np.nan,np.nan]
        dtype=np.float32
        fillval=DefaultData.get_default_fill_value(dtype)
        default_array = DefaultData.create_default_vector(np.size(telem2wr), dtype)
        variable = Variable(["x_4"], default_array)
        tu.add_fill_value(variable, fillval)
        variable.attrs["long_name"] = telem_des
        mfggrp[temp +"_"+ name] = variable
        #Store data
        mfggrp.variables[temp +"_"+ name].data   = telem2wr


def write_uncert(mfggrp,uncert,dimvals,rel,m2,debug=0):
  '''
  writing uncertainties
  Note: this has to be predefined in the netcdf template if Release>2
  HISTORY:
  5.1.2018: modified to match new writer in case of Release>2
  '''
  if float(rel) < 2.0:
    unc_names      =  ( 'u_a0_vis','u_a1_vis',"u_digitization_counts_vis")
    unc_scale_factor = ( 1.0, 1.0, 1.0)
    unc_add_offset   = ( 0.0, 0.0, 0.0)
  else:
    unc_names      =  ( 'u_solar_irradiance_vis', 'covariance_spectral_response_function_vis', 'u_a0_vis',  \
                        'u_a1_vis','covariance_a0_a1_vis','u_latitude', 'u_longitude', 'u_time', \
                        'u_solar_zenith_angle', \
                        "u_digitization_counts_vis" )
  for item in range(len(uncert)):
    if debug>=1:
      print unc_names[item]
    unc_inp  =  np.array( uncert[item] )
    if debug>=1:
     print np.amin( unc_inp ), np.amax( unc_inp )
    if float(rel) < 2.0:
      try:
        unc_data_type=mfggrp.variables[unc_names[item]].dtype
      except KeyError:
        print mfggrp.keys()
        exit(2)
    if float(rel) < 2.0:
      scale_factor = unc_scale_factor[item]
      add_offset   = unc_add_offset[item]    
      mfggrp.variables[unc_names[item]].attrs["scale_factor"] = scale_factor
      mfggrp.variables[unc_names[item]].attrs["add_offset"]   = add_offset
      unc_inp  =  ( unc_inp / scale_factor ) - add_offset
    #for ulat/lon also set fillvalues
    if( unc_names[item] == "u_latitude" or unc_names[item] == "u_longitude" ):
      if float(rel) < 2.0:
        unc_fillvalue        = mfggrp.variables[unc_names[item]].attrs['_FillValue']
        unc_inp[m2==0] = unc_fillvalue
      else:
        unc_fillvalue=np.nan#DefaultData.get_default_fill_value(unc_inp.dtype)
        unc_inp[m2==0] = unc_fillvalue
    if debug>=1:
      print np.shape(unc_inp)
    if float(rel) < 2.0:
      try:
        mfggrp.variables[unc_names[item]].data = np.fliplr(np.array( unc_inp )).astype(unc_data_type)
      except ValueError:
        mfggrp.variables[unc_names[item]].data = np.array( unc_inp ).astype(unc_data_type)
    else:
      mfggrp.variables[unc_names[item]].data = np.array( unc_inp )

def write_aux(mfggrp,refl_aux,dimvals,rel,debug=0):
  ### Writing auxiliary data needed for computation of reflectance
  
  if float(rel) < 2.0:
    
    refl_aux_names      =  ( 'allan_deviation_counts_space_vis', 'mean_count_space_vis','years_since_launch','a0_vis','a1_vis') #FIXME: add all other missing vars that are in template
    #refl_aux_long_names =  ( 'decimal years since launch of satellite',\
                             #'a0 - calibration coefficient at launch date of vis channel', \
                             #'a1 - time drift of calibration coefficient - vis channel' )
    #refl_aux_units      =  ( 'years','Wm-2Sr-1/DC','Wm-2Sr-1/DC/Y')
    
    refl_aux_scale_factor = ( 1.0, 1.0, 1.0, 1.0, 1.0)
    refl_aux_add_offset   = ( 0.0, 0.0, 0.0, 0.0, 0.0)
    
    #refl_aux_data_type    =  ( np.float32, np.float32, np.float32)

  else:
  
    refl_aux_names = ('distance_sun_earth', 'solar_irradiance_vis', \
                      'allan_deviation_counts_space_vis', 'mean_count_space_vis', \
                      #'sub_sat_latitude_start', \ #FIXME: have to be added to writer template
                      #'sub_sat_longitude_start', 'sub_sat_latitude_end', 'sub_sat_longitude_end', \
                      'years_since_launch', 'a0_vis', 'a1_vis' )

    #refl_aux_long_names = ( 'distance_sun_earth', 'solar_irradiance_vis', \
                            #'spectral_response_function_vis (normalised to max = 1)', \
                            #'count_vis_space_allan_deviation (for details see: http://allantools.readthedocs.io/en/latest/)', \
                            #'count_vis_space_mean', 'sub_sat_latitude_start', \
                            #'sub_sat_longitude_start', 'sub_sat_latitude_end', 'sub_sat_longitude_end', \
                            #'decimal years since launch of satellite',\
                            #'a0 - calibration coefficient at launch date of vis channel', \
                            #'a1 - time drift of calibration coefficient - vis channel' )
                              
    #refl_aux_units = ( 'au', 'Wm-2', 'micrometre', 'DC', 'DC', 'degree_north', 'degree_east', \
                      #'degree_north', 'degree_east', 'years' 'Wm-2Sr-1/DC', 'Wm-2Sr-1/DC/Year' )
    
    refl_aux_scale_factor = ( 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)#, 1.0, 1.0, 1.0 )
    refl_aux_add_offset   = ( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)#, 0.0, 0.0, 0.0 )

    #refl_aux_data_type    =  ( np.float32, np.float32, np.float32, np.float32, np.float32, np.float32, np.float32, np.float32, np.float32, np.float32, np.float32,np.float32 )
  
  #refl_aux_fillvalue   =  9.96920996839e+36

  for item in range(len(refl_aux)):
    if debug>=1:
      print refl_aux_names[item]
      print refl_aux[item]
      print dimvals[0]
   
    shape = np.size(refl_aux[item])
    dims = len( np.shape( refl_aux[item] ) )
    full_name=refl_aux_names[item]
    #long_name=refl_aux_long_names[item]
    #dtype=refl_aux_data_type[item]
    #fillval=DefaultData.get_default_fill_value(dtype)
    #auxvar=refl_aux[item]
    #units=refl_aux_units[item] 
    #test if it is in template
    if not full_name in mfggrp.variables:
      print "\nERROR: inconsistency between write_aux and template. Fix this!\n"
      #print mfggrp.variables
      print full_name
      exit()
    else:
      refl_aux_inp  =  np.array( refl_aux[item] )
      if debug>=1:
        print np.amin( refl_aux_inp ), np.amax( refl_aux_inp )
      
      try:#first try to get template's scaling
        scale_factor = mfggrp.variables[full_name].attrs["scale_factor"]
        add_offset   = mfggrp.variables[full_name].attrs["add_offset"]
      except KeyError:#if not defined get it from hardcoded dict
        scale_factor = refl_aux_scale_factor[item]
        add_offset   = refl_aux_add_offset[item]
        mfggrp.variables[full_name].attrs["scale_factor"] = scale_factor
        mfggrp.variables[full_name].attrs["add_offset"]   = add_offset
        
      refl_aux_inp  =  ( refl_aux_inp / scale_factor ) - add_offset

      mfggrp.variables[full_name].data   = refl_aux_inp



def save_static( filename,nomlon,varnames,vals,m1,comment,rel,debug=0):
  #print '    static filename: ', filename
  with nc.Dataset(filename, 'w',  format='NETCDF4') as ncset:
    pixel_VIS = ncset.createDimension('x_vis', 5000)
    line_VIS  = ncset.createDimension('y_vis', 5000)
    pixel     = ncset.createDimension('x_ir_wv', 2500)
    line      = ncset.createDimension('y_ir_wv', 2500)
    if debug==1:
      print "dims created"
    #fill dimensions
    #pixel_VIS=np.arange(4999,0,-1)
    #line_VIS =np.arange(0,4999,1)
    #pixel    =np.arange(2499,0,-1)
    #line     =np.arange(0,2499,1)

    if debug==1:
      print "dims filled"
    # Global Attributes 
    ncset.description = 'MVIRI Level 1.5 static fundamental climate data record'  
    ncset.history     = comment+' \nCreated: ' + time.ctime(time.time())  
    ncset.source      = 'EUMETSAT'
    ncset.license     = "This dataset is released for use under CC-BY licence (https://creativecommons.org/licenses/by/4.0/) and was developed in the EC " \
                        "FIDUCEO project \"Fidelity and Uncertainty in Climate Data Records from Earth " \
                        "Observations\". Grant Agreement: 638822."
    LAT=vals[0]
    LON=vals[1]
    if debug==1:
      print "Global Atts set"
    for varname,val in zip(varnames,vals):
      if debug>=1:
        print varname
      #if viewing geometry is to be written do this:+++++++++++++++++++++++++++
      if "viewing_geometry" in varname:
        ncset=geomLUT(ncset,nomlon,m1,LAT,LON,halfrangelat=5,halfrangelon=3,step=0.25,debug=debug)
      #if sensitivities are to be written do this:+++++++++++++++++++++++++++++
      elif "s_" in varname:
        ncset=sensiLUT(ncset,varname,m1,LAT,LON)
      #if lats/lons are written do this:+++++++++++++++++++++++++++++++++++++++
      else:
        print varname
        val=val
        val[val<(-180)]=nc.default_fillvals['i2']
        if np.shape(val)[0]>2500:
          if debug>=1:
            print "vis"
          arr  = ncset.createVariable(varname,'i2',dimensions=('y_vis','x_vis'),zlib=True)
        else:
          if debug>=1:
            print "ir"
          arr  = ncset.createVariable(varname,'i2',dimensions=('y_ir_wv','x_ir_wv'),zlib=True)
        # Variable specific Attributes
        if "latitude" in varname:
          arr.units = 'degree_north'
          arr.scale_factor=0.0027466658
          arr.standard_name="latitude"
        if "longitude" in varname:
          arr.units = 'degree_east'
          arr.scale_factor=0.0054933317
          arr.standard_name="longitude"
        # Variable common Attributes
        arr.fill_value=nc.default_fillvals['i2']
        #fill content
        if debug>=2:
          to.printimage(val,varname,mini=-90,maxi=90)
        arr[:]=val

#----------------------------------------------------------------------------
#Write VZA and VAZ for several ssp's as LUT
#----------------------------------------------------------------------------
def geomLUT(ncset,nomlon,m1,LAT,LON,halfrangelat=5,halfrangelon=3,step=0.25,debug=0):
  scale=0.0054933317
  latrange= np.arange(0-halfrangelat,0+halfrangelat+step,step)
  lonrange= np.arange(nomlon-halfrangelon,nomlon+halfrangelon+step,step)
  ncset.createDimension('sub_sat_latitude' , len(latrange))
  ncset.createDimension('sub_sat_longitude', len(lonrange))
  sslats   = ncset.createVariable('sub_sat_latitude','d',('sub_sat_latitude',))
  sslats[:]= latrange
  sslons   = ncset.createVariable('sub_sat_longitude','d',('sub_sat_longitude',))
  sslons[:]= lonrange
  arrVZA   = ncset.createVariable('satellite_zenith_angle' ,'u2',dimensions=('sub_sat_latitude','sub_sat_longitude','y_vis','x_vis'),zlib=True)
  arrVAZ   = ncset.createVariable('satellite_azimuth_angle','u2',dimensions=('sub_sat_latitude','sub_sat_longitude','y_vis','x_vis'),zlib=True)
  arrVZA.units = 'degree'
  arrVAZ.units = 'degree'
  arrVAZ.scale_factor=scale
  arrVZA.scale_factor=scale
  ilat    = 0
  for sslat in latrange:
    ilon    = 0
    for sslon in lonrange:
      VZA,VAZ               = cr.vza(sslat,sslon,m1,LAT,LON)
      arrVZA[ilat,ilon,:,:] = np.fliplr(VZA)#/scale).astype(np.uint16)
      arrVAZ[ilat,ilon,:,:] = np.fliplr(VAZ)#/scale).astype(np.uint16)
      if debug>=1:
        print sslat,sslon,VZA[2500,2500],VZA[2500,2500]/scale,(VZA[2500,2500]/scale).astype(np.uint16),arrVZA[ilat,ilon,2500,2500]
        if debug>=2:
          testarr=np.fliplr(VZA/scale).astype(np.uint16)
          to.printimage(testarr*scale,"testarr",mini=0,maxi=120)
          to.printimage(arrVZA[ilat,ilon,:,:],"ncarr",mini=0,maxi=120)
      ilon=ilon+1
    ilat=ilat+1
    if debug>=1:
      if ilat==len(latrange)/2 and ilon==len(lonrange)/2:
        plt.plot(arrVZA[:,ilon,2500,2500])
        plt.show()
  return ncset
#----------------------------------------------------------------------------
#test ssp resolution for VZA and VAZ LUT
#----------------------------------------------------------------------------
def TSTgeomLUT(staticname,nomlon,m1,LAT,LON,calc=0):
  from scipy.interpolate import RegularGridInterpolator
  #define pixel of interest
  x=2500
  y=2500
  #define resolution of lut
  halfrangelat=5
  halfrangelon=3
  step=0.25
  if calc==0:
    static=nc.Dataset(staticname,"r")
    lutssplats=static.variables['sub_sat_latitude'][:]
    lutssplons=static.variables['sub_sat_longitude'][:]
    static.close()
  else:
    lutssplats = np.arange(0-halfrangelat,0+halfrangelat+step,step)
    lutssplons = np.arange(nomlon-halfrangelon,nomlon+halfrangelon+step,step)
  #read example SSP TS
  exfile="../aux/MET7_example.txt"
  a=np.loadtxt(exfile)
  a=np.array(a)
  a[:,1][np.isnan(a[:,1])]=0
  a=a[~np.isnan(a).any(axis=1)]
  latrange=a[:,2]
  lonrange=a[:,3]
  #prepare full array
  arrVZA=np.zeros((len(latrange)))
  arrVAZ=np.zeros((len(latrange)))
  lutlatrange=np.zeros((len(latrange)))
  lutlonrange=np.zeros((len(latrange)))
  #prepare lut array
  lutVZA   = np.zeros((len(latrange)))
  lutVZA_A = np.zeros((4,4))
  lutVAZ   = np.zeros((len(latrange)))
  lutVAZ_A = np.zeros((4,4))
  lutlat   = np.zeros((4,4))
  lutlon   = np.zeros((4,4))
  print len(latrange)
  for i in range(len(latrange)):
    print i
    sslat  = latrange[i]
    sslon  = lonrange[i]
    #calculate VZA/VAZ
    VZA,VAZ   = cr.vza(sslat,sslon,m1,LAT,LON)
    arrVZA[i] = np.fliplr(VZA)[x,y]#/scale).astype(np.uint16)
    arrVAZ[i] = np.fliplr(VAZ)[x,y]#/scale).astype(np.uint16)
    #get VZA/VAZ from lut
    #first get closest match index
    ilat = (np.abs(lutssplats-sslat)).argmin()
    ilon = (np.abs(lutssplons-sslon)).argmin()
    sslatlut=lutssplats[ilat]
    sslonlut=lutssplons[ilon]
    #get distance to closest match
    distlat=sslatlut-sslat
    distlon=sslonlut-sslon
    #get borders
    if distlat>0:
      ilat1=ilat-2
      ilat2=ilat+1
    else:
      ilat2=ilat+2
      ilat1=ilat-1
    if distlon>0:
      ilon1=ilon-2
      ilon2=ilon+1
    else:
      ilon2=ilon+2
      ilon1=ilon-1
    sslatlut2=lutssplats[ilat1:ilat2+1]
    sslonlut2=lutssplons[ilon1:ilon2+1]
    #get value of closest match
    if calc==0:
      static=nc.Dataset(staticname,"r")
      lutVZA_A = (static.variables['satellite_zenith_angle'][ilat1:ilat2+1,ilon1:ilon2+1,x,y])
      lutVAZ_A = (static.variables['satellite_azimuth_angle'][ilat1:ilat2+1,ilon1:ilon2+1,x,y])
      static.close()
    for lat,k in zip(sslatlut2,range(len(sslatlut2))):
      for lon,j in zip(sslonlut2,range(len(sslatlut2))):
        if calc!=0:
          VZA,VAZ   = cr.vza(lat,lon,m1,LAT,LON)
          lutVZA_A[k,j] = np.fliplr(VZA)[x,y]
          lutVAZ_A[k,j] = np.fliplr(VAZ)[x,y]
        lutlat[k,j]=lat
        lutlon[k,j]=lon
    #interpolate 
    Z_VZA=to.interpol_cubic(lutlon,lutlat,lutVZA_A,sslon,sslat)
    Z_VAZ=to.interpol_cubic(lutlon,lutlat,lutVAZ_A,sslon,sslat)
    #store
    print np.shape(Z_VZA)
    lutVZA[i]=Z_VZA
    lutVAZ[i]=Z_VAZ
    #save also lut-ssp
    lutlatrange[i]=sslatlut
    lutlonrange[i]=sslonlut
  #plot
  plt.plot(latrange   ,lonrange)
  plt.plot(lutlatrange,lutlonrange)
  plt.show()
  plt.plot(arrVZA[:])
  plt.plot(lutVZA[:])
  plt.show()
  plt.plot(arrVZA[:]-lutVZA[:])
  plt.show()
  #write to txtfile
  outfile='../test/LUTresolutiontest_x'+str(x)+'_y'+str(y)+'resXY'+str(step)+'.txt'
  out=latrange
  for arr in [lonrange,lutlatrange,lutlonrange,arrVZA,lutVZA,arrVAZ,lutVAZ]:
    print np.shape(arr)
    out=np.column_stack((out,arr))
  np.savetxt(outfile,out,delimiter=';')
  return
#----------------------------------------------------------------------------
#Write solar geometries for several slots/doys
#eventually not needed because based on wrong acquisition time calculations
#----------------------------------------------------------------------------
def solarLUT(ncset,varname,m1,LAT,LON,debug=0):
  month    = ncset.createDimension('month', 12)
  slot     = ncset.createDimension('slot', 48)
  month    = np.arange(0,12,1)
  slot     = np.arange(0,48,1)
  arrSZA   = ncset.createVariable("solar_zenith_angle",'f4',dimensions=('month','slot','y_vis','x_vis'),zlib=True)
  arrSAZ   = ncset.createVariable("solar_azimuth_angle",'f4',dimensions=('month','slot','y_vis','x_vis'),zlib=True)
  arrSZA.units = 'degree'
  arrSAZ.units = 'degree'
  
  for mo in range (1,13):
    DOY=15+(mo*30) #DOY for every middle of month
    for slot in range (0,48):
      if debug>=1:
        print "slot ",slot
      t=(slot/2.)*(60*60)#time in seconds since midnight
      if "solar" in varname:
        t1=array([range(0,5000),]*5000).transpose()
        t1=t1*60/100+t
        SZA,SAZ=cr.sza(t1,DOY,m1,LAT,LON)

      if debug>=1:
        print np.shape(SZA)
        print np.shape(arrSZA[mo-1,slot,:,:])
      arrSZA[mo-1,slot,:,:]=SZA
      arrSAZ[mo-1,slot,:,:]=SAZ
      
  return ncset

#----------------------------------------------------------------------------
#Write sensitivities for several slots/doys
#----------------------------------------------------------------------------
def sensiLUT(ncset,varname,m1,LAT,LON,debug=0):
  try:
    month    = ncset.createDimension('month', 12)
    slot     = ncset.createDimension('slot', 48)
  except RuntimeError: 
    print "Error: NetCDF: String match to name in use"
    print "--> probably because dimension already defined. Ignoring..."
  month    = np.arange(0,12,1)
  slot     = np.arange(0,48,1)
  arr      = ncset.createVariable(varname,'f4',dimensions=('month','slot','y_vis','x_vis'),zlib=True)
  for mo in range (1,13):
    DOY=15+(mo*30) #DOY for every middle of month
    for slot in range (0,48):
      if debug>=1:
        print "slot ",slot
      t=(slot/2.)*(60*60)#time in seconds since midnight
      if "time" in varname:
        tmp=ef.sensi_SZA(t,120,DOY,m1,LAT,0,LON,0,debug=debug)
        arr.units = 'degree per second'
      elif "latitude" in varname:
        tmp=ef.sensi_SZA(t,0,DOY,m1,LAT,0.05,LON,0,debug=debug)
        arr.units = 'degree per degree'
      elif "longitude" in varname:
        tmp=ef.sensi_SZA(t,0,DOY,m1,LAT,0,LON,0.05,debug=debug)
        arr.units = 'degree per degree'
      if debug>=1:
        print np.shape(tmp)
        print np.shape(arr[mo-1,slot,:,:])
      arr[mo-1,slot,:,:]=tmp
  return ncset


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#reading-related
#-----------------------------------------------------------------------


def ncdump(nc_fid, verb=True):
    '''
    Ref: http://schubert.atmos.colostate.edu/~cslocum/netcdf_example.html
    ncdump outputs dimensions, variables and their attribute information.
    The information is similar to that of NCAR's ncdump utility.
    ncdump requires a valid instance of Dataset.

    Parameters
    ----------
    nc_fid : netCDF4.Dataset
        A netCDF4 dateset object
    verb : Boolean
        whether or not nc_attrs, nc_dims, and nc_vars are printed

    Returns
    -------
    nc_attrs : list
        A Python list of the NetCDF file global attributes
    nc_dims : list
        A Python list of the NetCDF file dimensions
    nc_vars : list
        A Python list of the NetCDF file variables
    '''
    def print_ncattr(key):
        """
        Prints the NetCDF file attributes for a given key

        Parameters
        ----------
        key : unicode
            a valid netCDF4.Dataset.variables key
        """
        try:
            print "\t\ttype:", repr(nc_fid.variables[key].dtype)
            for ncattr in nc_fid.variables[key].ncattrs():
                print '\t\t%s:' % ncattr,\
                      repr(nc_fid.variables[key].getncattr(ncattr))
        except KeyError:
            print "\t\tWARNING: %s does not contain variable attributes" % key

    #NetCDF global attributes
    nc_attrs = nc_fid.ncattrs()
    if verb:
        print "NetCDF Global Attributes:"
        for nc_attr in nc_attrs:
            print '\t%s:' % nc_attr, repr(nc_fid.getncattr(nc_attr))
    nc_dims = [dim for dim in nc_fid.dimensions]  # list of nc dimensions
    # Dimension shape information.
    if verb:
        print "NetCDF dimension information:"
        for dim in nc_dims:
            print "\tName:", dim 
            print "\t\tsize:", len(nc_fid.dimensions[dim])
            print_ncattr(dim)
    # Variable information.
    nc_vars = [var for var in nc_fid.variables]  # list of nc variables
    if verb:
        print "NetCDF variable information:"
        for var in nc_vars:
            if var not in nc_dims:
                print '\tName:', var
                print "\t\tdimensions:", nc_fid.variables[var].dimensions
                print "\t\tsize:", nc_fid.variables[var].size
                print_ncattr(var)
    nc_units = []
    nc_lnames= []
    for var in nc_vars:
      try:
        nc_units.append(nc_fid.variables[var].getncattr("units"))
      except AttributeError:
        nc_units.append(u"1")
      try:
        nc_lnames.append(nc_fid.variables[var].getncattr("long_name"))
      except AttributeError:
        nc_lnames.append(u"")
    return nc_attrs, nc_dims, nc_vars, nc_units, nc_lnames


def read_calfile(calfile,satellite,doy,s_datestring,s_Time,logs,debug=0):
  """
  calfile_content:
  [u'year', u'month', u'day', u'slot', u'julian_time', u'a_ir', u'b_ir', u'a_wv', u'b_wv']
  NOTE: slot  umber in calfile is wrong! use julian time (end of slot) for indexing.
  """
  d=datetime.datetime.strptime(s_datestring+s_Time, '%Y%m%d%H%M')
  jultime=to.julian_day(d)
  with nc.Dataset( calfile ) as cal_data:
    nc_attrs, nc_dims, nc_vars, nc_units, nc_lnames= ncdump(cal_data,verb=False)
    try:
      jultimes = cal_data.variables["julian_time"][:]
    except:
      err=logs.update(satellite,curr_datetime,"status_easy",-1)
      err=logs.update(satellite,curr_datetime,"problem_easy","strange thing with IR calfile")
      err=logs.update(satellite,curr_datetime,"status_full",-1)
      err=logs.update(satellite,curr_datetime,"problem_full","strange thing with IR calfile")
      raise
    ix= (np.abs(jultimes - jultime)).argmin()
    calfile_content = []
    calfile_names   = []
    calfile_units   = []
    calfile_lnames  = []
    for var,unit,lnm in zip(nc_vars,nc_units,nc_lnames):
        if var in 'a_ir b_ir a_wv b_wv':
          tmp=cal_data.variables[var][ix]#read
          calfile_content.append(tmp)
          calfile_names.append(var)
          calfile_units.append(unit)
          calfile_lnames.append(lnm)
    if debug>=2: 
      print calfile_content
    return (calfile_content,calfile_names,calfile_units,calfile_lnames)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#classes from fiduceo adapted to release 1.0
#-----------------------------------------------------------------------

from xarray import Variable

FULL_DIMENSION = 5000
IR_DIMENSION = 2500
SRF_SIZE = 101

TIME_FILL_VALUE = -32768

class FCDRWriterInternal:
    _version = "1.0.5"
    @classmethod
    def write(cls, ds, file, compression_level=None, overwrite=False):
        """
        Save a dataset to NetCDF file.
        :param ds: The dataset
        :param file: File path
        :param compression_level: the file compression level, 0 - 9, default is 5
        :param overwrite: set true to overwrite existing files
         """
        if os.path.isfile(file):
            if overwrite is True:
                os.remove(file)
            else:
                raise IOError("The file already exists: " + file)
        
        # set up compression parameter for ALL variables. Unfortunately, xarray does not allow
        # one set of compression params per file, only per variable. tb 2017-01-25
        if compression_level is None:
            compression_level = 5

        comp = dict(zlib=True, complevel=compression_level)
        encoding = dict((var, comp) for var in ds.data_vars)

        ds.to_netcdf(file, format='netCDF4', engine='netcdf4', encoding=encoding)
    @classmethod
    def createTemplateBasic(cls, sensorType, height):
        """
        Create a template dataset in EASY FCDR format for the sensor given as argument.
        :param sensorType: the sensor type to create the template for
        :param height the hheight in pixels of the data product
        :return the template dataset
        """
        dataset = xr.Dataset()
        cls._add_standard_global_attributes(dataset)

        template_factory = TemplateFactoryInternal()

        sensor_template = template_factory.get_sensor_template(sensorType)
        sensor_template.add_basic_fcdr_variables(dataset, height)

        return dataset
    @classmethod
    def _add_standard_global_attributes(cls, dataset):
        dataset.attrs["Conventions"] = "CF-1.6"
        dataset.attrs["licence"]     = "This dataset is released for use ..."
        dataset.attrs["writer_version"] = FCDRWriter._version
        # @todo tb/tb 2 the following dictionary entries have to be supplied by the data generators
        dataset.attrs["institution"] = None
        dataset.attrs["title"]       = None
        dataset.attrs["source"]      = None
        dataset.attrs["history"]     = None
        dataset.attrs["references"]  = None
        dataset.attrs["comment"]     = None
      
class TemplateFactoryInternal:
    def __init__(self):
        self.templates = dict([("MVIRI", MVIRIInternal)])

    def get_sensor_template(self, name):
        return self.templates[name]
      
class MVIRIInternal:
    @staticmethod
    def add_basic_fcdr_variables(dataset, height):
        # height is ignored - supplied just for interface compatibility tb 2017-02-05
        # time
        default_array = DefaultData.create_default_array(IR_DIMENSION, IR_DIMENSION, np.uint16)
        variable = Variable(["y_ir", "x_ir"], default_array)
        tu.add_fill_value(variable, DefaultData.get_default_fill_value(np.uint16))
        variable.attrs["standard_name"] = "time"
        variable.attrs["long_name"] = "Acquisition time of pixel"
        tu.add_units(variable, "seconds since 1970-01-01 00:00:00")
        tu.add_offset(variable, TIME_FILL_VALUE)
        dataset["time"] = variable

        dataset["years_since_launch"] = tu.create_scalar_float_variable(long_name="fractional year since launch of satellite", units="years")

        #dataset["satellite_azimuth_angle"] = MVIRIInternal._create_angle_variable_int(0.005493164, standard_name="sensor_azimuth_angle", unsigned=True)
        #dataset["satellite_zenith_angle"] = MVIRIInternal._create_angle_variable_int(0.005493248, standard_name="sensor_zenith_angle")
        #dataset["solar_azimuth_angle"] = MVIRIInternal._create_angle_variable_int(0.005493164, standard_name="solar_azimuth_angle", unsigned=True)
        #dataset["solar_zenith_angle"] = MVIRIInternal._create_angle_variable_int(0.005493248, standard_name="solar_zenith_angle")

        # count_ir
        default_array = DefaultData.create_default_array(IR_DIMENSION, IR_DIMENSION, np.uint8)
        variable = Variable(["y_ir", "x_ir"], default_array)
        tu.add_fill_value(variable, DefaultData.get_default_fill_value(np.uint8))
        variable.attrs["long_name"] = "Infrared Image Counts"
        tu.add_units(variable, "count")
        dataset["count_ir"] = variable

        # count_wv
        default_array = DefaultData.create_default_array(IR_DIMENSION, IR_DIMENSION, np.uint8)
        variable = Variable(["y_ir", "x_ir"], default_array)
        tu.add_fill_value(variable, DefaultData.get_default_fill_value(np.uint8))
        variable.attrs["long_name"] = "WV Image Counts"
        tu.add_units(variable, "count")
        dataset["count_wv"] = variable

        ## distance_sun_earth
        #dataset["distance_sun_earth"] = tu.create_scalar_float_variable(long_name="Sun-Earth distance", units="au")

        ## sol_eff_irr
        #dataset["solar_irradiance_vis"] = tu.create_scalar_float_variable(standard_name="solar_irradiance_vis", long_name="Solar effective Irradiance", units="W*m-2")

        ## u_sol_eff_irr
        #dataset["u_solar_irradiance_vis"] = tu.create_scalar_float_variable(standard_name="u_solar_irradiance_vis", long_name="Uncertainty of Solar effective Irradiance", units="W*m-2")

        ## srf
        #default_array = DefaultData.create_default_vector(SRF_SIZE, np.float32)
        #variable = Variable(["srf_size"], default_array)
        #tu.add_fill_value(variable, DefaultData.get_default_fill_value(np.float32))
        #variable.attrs["long_name"] = "Spectral Response Function"
        #dataset["spectral_response_function_vis"] = variable

        ## srf covariance
        #default_array = DefaultData.create_default_array(SRF_SIZE, SRF_SIZE, np.float32)
        #variable = Variable(["srf_size", "srf_size"], default_array)
        #tu.add_fill_value(variable, DefaultData.get_default_fill_value(np.float32))
        #variable.attrs["long_name"] = "Covariance of the Visible Band Spectral Response Function"
        #dataset["covariance_spectral_response_function_vis"] = variable

        #IR stuff
        dataset["a_ir"] = tu.create_scalar_float_variable(long_name="Calibration parameter a for IR Band", units="mWm^-2sr^-1cm^-1")
        dataset["b_ir"] = tu.create_scalar_float_variable(long_name="Calibration parameter b for IR Band", units="mWm^-2sr^-1cm^-1/DC")
        dataset["u_a_ir"] = tu.create_scalar_float_variable(long_name="Uncertainty of calibration parameter a for IR Band", units="mWm^-2sr^-1cm^-1")
        dataset["u_b_ir"] = tu.create_scalar_float_variable(long_name="Uncertainty of calibration parameter b for IR Band", units="mWm^-2sr^-1cm^-1/DC")
        dataset["a_wv"] = tu.create_scalar_float_variable(long_name="Calibration parameter a for WV Band", units="mWm^-2sr^-1cm^-1")
        dataset["b_wv"] = tu.create_scalar_float_variable(long_name="Calibration parameter b for WV Band", units="mWm^-2sr^-1cm^-1/DC")
        dataset["u_a_wv"] = tu.create_scalar_float_variable(long_name="Uncertainty of calibration parameter a for WV Band", units="mWm^-2sr^-1cm^-1")
        dataset["u_b_wv"] = tu.create_scalar_float_variable(long_name="Uncertainty of calibration parameter b for WV Band", units="mWm^-2sr^-1cm^-1/DC")
        dataset["q_ir"] = tu.create_scalar_float_variable(long_name="IR Band Calibration quality flag", units="1")
        dataset["q_wv"] = tu.create_scalar_float_variable(long_name="WV Band Calibration quality flag", units="1")
        dataset["unit_conversion_ir"] = tu.create_scalar_float_variable(long_name="IR Unit conversion factor", units="1")
        dataset["unit_conversion_wv"] = tu.create_scalar_float_variable(long_name="WV Unit conversion factor", units="1")
        dataset["bt_a_ir"] = tu.create_scalar_float_variable(long_name="IR Band BT conversion parameter A", units="1")
        dataset["bt_b_ir"] = tu.create_scalar_float_variable(long_name="IR Band BT conversion parameter B", units="1")
        dataset["bt_a_wv"] = tu.create_scalar_float_variable(long_name="WV Band BT conversion parameter A", units="1")
        dataset["bt_b_wv"] = tu.create_scalar_float_variable(long_name="WV Band BT conversion parameter B", units="1")

        #VIS stuff
        
        # count_vis
        default_array = DefaultData.create_default_array(FULL_DIMENSION, FULL_DIMENSION, np.uint8)
        variable = Variable(["y_vis", "x_vis"], default_array)
        tu.add_fill_value(variable, DefaultData.get_default_fill_value(np.uint8))
        variable.attrs["long_name"] = "Image counts"
        tu.add_units(variable, "count")
        dataset["count_vis"] = variable
        
        dataset["a0_vis"] = tu.create_scalar_float_variable("Calibration Coefficient at Launch", units="Wm^-2sr^-1/count")
        dataset["a1_vis"] = tu.create_scalar_float_variable("Time variation of a0", units="Wm^-2sr^-1/count day^-1 10^5")
        dataset["mean_count_space_vis"] = tu.create_scalar_float_variable("Space count combined for both VIS detectors", units="count")
        dataset["allan_deviation_counts_space_vis"] = tu.create_scalar_float_variable("Uncertainty of space count combined for both VIS detectors", units="count")
        
        # u_combined_counts_vis
        #dataset["u_combined_counts_vis"] = tu.create_scalar_float_variable("Uncertainty of counts combined for both VIS detectors", units="count")

        # u_a0_vis
        variable = tu.create_scalar_float_variable("Uncertainty in a0", units="Wm^-2sr^-1/count")
        dataset["u_a0_vis"] = variable

        # u_a1_vis
        variable = tu.create_scalar_float_variable("Uncertainty in a1", units="Wm^-2sr^-1/count day^-1 10^5")
        dataset["u_a1_vis"] = variable

        # covariance_a0_a1_vis
        dataset["covariance_a0_a1_vis"] = tu.create_scalar_float_variable(long_name="Covariance of calibration coefficients")
        #default_array = DefaultData.create_default_array(1, 1, np.float32)
        #variable = Variable(["calib_coeff_cov_size", "calib_coeff_cov_size"], default_array)
        #tu.add_fill_value(variable, DefaultData.get_default_fill_value(np.float32))
        #variable.attrs["long_name"] = "Covariance of calibration coefficients"
        #dataset["covariance_a0_a1_vis"] = variable

        #dataset["u_electronics_counts_vis"] = tu.create_scalar_float_variable("Uncertainty due to Electronics noise", units="count")
        dataset["u_digitization_counts_vis"] = tu.create_scalar_float_variable("Uncertainty due to digitization", units="count")

    @staticmethod
    def _create_angle_variable_int(scale_factor, standard_name=None, long_name=None, unsigned=False, fill_value=None):
        if unsigned is True:
            data_type = np.uint16
        else:
            data_type=np.int16

        default_array = DefaultData.create_default_array(FULL_DIMENSION, FULL_DIMENSION, data_type, fill_value=fill_value)
        variable = Variable(["y_vis", "x_vis"], default_array)

        if fill_value is not None:
            tu.add_fill_value(variable, fill_value)
        else:
            tu.add_fill_value(variable, DefaultData.get_default_fill_value(data_type))

        if standard_name is not None:
            variable.attrs["standard_name"] = standard_name

        if long_name is not None:
            variable.attrs["long_name"] = long_name

        tu.add_units(variable, "degree")
        tu.add_scale_factor(variable, scale_factor)
        return variable
