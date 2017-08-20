import os
import numpy as np
import struct
import sys
from collections import namedtuple
  
import _vax2ieee as vi

global recordsize
global rec_fmt

recordsize=15360 #records are described further down in structures()
rec_fmt=[['<8siiiiii8ciii12c8s2s6cih4c2c2chi2h5ii72cii3030c3030c3030c3030c4i4i4i4i4i2952c'],
	 ['<iiiiii676f676f5412c4s2s2c2f2f2f2f2f2f2f2f64c2f2f2f2f2f2f2f2f64cff84cddddd4c756c6d6d3f3f6d6dffi272ci4f4f6d6diiiiii2704c'],
	 ['<ihhhhh2shh316c316c316c316ch3h3h3h28s12b3hh6c3c6c3c5c5c6c3c6c3c5c5c3c5c5c3c5c15c8c4cdddd3d2difbbbbbbbbbbbb20c16c16b12hd3d3dffff4b4h4hff4fi4f4f4f4f4h4h88c800sh2cfiiiiiifff2f2f2i2f2f1024fff2f2f2i2f2f1024ff4f4fihi2c4f4f4f4f4f4f4f4f3d22c4102c']]


def vax(value,typ):
  #floats and doubles have to be treated specially 
  #because of the vax format
  if typ=='d':
    return vi.vax2d(value)  #for double: convert from vax
  elif typ=='f':
    return vi.vax2f(value)#for float: convert from vax
  else:
    return value#for rest: just take the value


def header(f):
  #open file
    records=range(1,1+1)#records=[1,2,3]
    

    for r in records:

      #handle records>3
      rq=r
      if r>3:
	rq=3
	
      #get information on record structure    
      names,dims=structures(rq)
      recID	='r'+str(r)
      struc 	= namedtuple(recID, names, verbose=False)
      
      #read and unpack record
      rec_bin = f.read(recordsize)
      record=struct.unpack(str(rec_fmt[rq-1][0]),rec_bin)
      
      #debug
      #print "\nrecord: "+str(r)+"\n"
      #print struct.calcsize(rec_fmt[rq-1][0])#print size of the structure
      #print len(record)
      #print len(names)
      
      #go through fields
      g=0
      lis=[]
      for i in range(0,len(names)):
	
	currname=names[i]
	currtype=currname[:1]
	
	#debug
	#print "name: "+str(names[i])
	#print "dims: "+str(dims[i])
	#print currtype
	
	#fill buffer to capture dimensional fields
	buf=[]
	for j in range(dims[i]):
	  buf.append(vax(record[g+j],currtype))
	#collect content in one list
	lis.append(buf)
	g=g+dims[i]
      
      #assign content to named tuple
      struc=struc(*lis)
      
      #debug
      print getattr(struc, 'szFilename')
      #print ', '.join(['{0}={1}'.format(k, getattr(struc, k)) for k in struc._fields])
      
def structures(z):
  import re
  #define names
  str_r1=['szFilename','iYear','iYday','iSlot_no','iData_type',\
	  'iIp_contents','iSpectral_content','szSpares1',\
	  'iRec_length','iNo_of_rec','iNo_of_rec_in_head',\
	  'szSpares2','szCtype','szSat_code','szCpares3',\
	  'iJulian_slot_no','shLast_line_of_image','Vis_channnel_in_use',\
	  'Ir_channnel_in_use','Wv_channnel_in_use','shScan_mode',\
	  'iSecs_past_midnight_for_shl','shVisible_south_north',\
	  'iFish_removel','iNo_lin_with_radpos_anomalies',\
	  'szSpares4','iFirst_line_used','iLast_line_used',\
	  'szVis_south_lost_line_table',\
	  'szVis_north_lost_line_table',\
	  'szIr_lost_line_table','szWv_lost_line_table',\
	  'lDrop_bit_test_control_codes','iDrop_bit_handling_stats_for_vis',\
	  'iDrop_bit_handling_stats_for_visn','iDrop_bit_handling_stats_for_ir',\
	  'iDrop_bit_handling_stats_for_wv','szSpare']
  dim_r1=[1,1,1,1,1,1,1,8,1,1,1,12,1,1,6,1,1,4,2,2,1,1,2,\
	  5,1,72,1,1,3030,3030,3030,3030,4,4,4,4,4,2952]

  #replace non-alphanumeric values
  i=0
  for stri in str_r1:
    str_r1[i]=re.sub(r'\W+', '',str_r1[i])
    i=i+1
  ##create tuples
  ##...with dimensions
  #d1		= namedtuple('d1', str_r1, verbose=False)
  #d1._make(dim_r1)
    
  #return
  if z==1:
    return (str_r1,dim_r1)

    
 #record1:
    #char szFilename[8];
    #int iYear;
    #int iYday;
    #int iSlot_no;
    #int iData_type,iIp_contents,iSpectral_content;
    #char szSpares1[8];
    #int iRec_length,iNo_of_rec,iNo_of_rec_in_head;
    #char szSpares2[12];
    #char szCtype[8],szSat_code[2];
    #char szCpares3[6];
    #int iJulian_slot_no;
    #short shLast_line_of_image;
    #char Vis_channnel_in_use[4];
    #char Ir_channnel_in_use[2];
    #char Wv_channnel_in_use[2];
    #short shScan_mode;
    #int iSecs_past_midnight_for_shl;
    #short shVisible_south_north[2];
    #int iFish_removel[5];
    #int iNo_lin_with_radpos_anomalies;
    #char szSpares4[72];
    #int iFirst_line_used,iLast_line_used;
    #char szVis_south_lost_line_table[3030];
    #char szVis_north_lost_line_table[3030];
    #char szIr_lost_line_table[3030];
    #char szWv_lost_line_table[3030];
    #int lDrop_bit_test_control_codes[4];
    #int iDrop_bit_handling_stats_for_vis[4];
    #int iDrop_bit_handling_stats_for_visn[4];  
    #int iDrop_bit_handling_stats_for_ir[4];
    #int iDrop_bit_handling_stats_for_wv[4];
    #char szSpare[2952];


#record2:
    #int iJulian_slot_no,iGeometric_quality,iNo_grid_points_in_def_matrix;
    #int iDeform_starts,iDeform_ends,iDeform_step;
    #float flX_deform_matrix[26][26],flY_deform_matrix[26][26];
    #char szSpares1[5412];
    #char szWefax_sat_id[4],szHr_sat_id[2];
    #char szSpares2[2];
    #float flVis1_reg_param_const[2];
    #float flVis2_reg_param_const[2];
    #float flVis3_reg_param_const[2];
    #float flVis4_reg_param_const[2];
    #float flIr1_reg_param_const[2];
    #float flIr2_reg_param_const[2];
    #float flWv1_reg_param_const[2];
    #float flWv2_reg_param_const[2];
    #char szSpares3[64];
    #float flVis1_reg_param_linear[2];
    #float flVis2_reg_param_linear[2];
    #float flVis3_reg_param_linear[2];
    #float flVis4_reg_param_linear[2];
    #float flIr1_reg_param_linear[2];
    #float flIr2_reg_param_linear[2];
    #float flWv1_reg_param_linear[2];
    #float flWv2_reg_param_linear[2];
    #char szSpares4[64];
    #float flSeasonal_amplitude,Seasonal_phase;
    #char szSpares5[84];
    #double dsEarth_equatorial_radius,dsEarth_polar_radius,dsDist_from_earth_centre;
    #double dsRadiometer_steppin_angle,dsNominal_sub_satellite_longitude;
    #char szCode_for_nominal_sub_satellite_point[4];
    #char szSpares6[756];
    #double dsOrbit_coordinates_in_mgf_at_image_start[6];
    #double dsOrbit_coordinates_in_mgf_at_image_end[6];
    #float flCart_comp_of_attitude_at_image_start[3];
    #float flCart_comp_of_attitude_at_image_end[3];
    #double dsOrbit_coordinates_in_fixed_earth_at_image_start[6];
    #double dsOrbit_coordinates_in_fixed_earth_at_image_end[6];
    #float flEarth_centre_in_pixels,flEarth_centre_in_lines;
    #int iProgram_error_level_stop_code;
    #char szSpares7[272];
    #int iLens_coefficient_type;
    #float flPixel_lens_correction_coefficient[4];
    #float flLine_lens_correction_coefficient[4];
    #double dsSun_coord_in_mgf_at_image_start[6];
    #double dsSun_coord_in_mgf_at_image_end[6];
    #int iNo_vis_south_horiz_pairs;
    #int iNo_vis_south_horiz_pairs_quality;
    #int iNo_vis_north_horiz_pairs;
    #int iNo_vis_north_horiz_pairs_quality;
    #int iNo_ir_horiz_pairs;
    #int iNo_ir_horiz_pairs_quality;
    #char szSpares8[2704];


    
#record3/2504
    #int iJulian_slot_no;
    #short shJulian_day,shSlot_no,shNominal_image_time;
    #short shSlot_type,shImage_quality_flag;
    #char szSat_id[2];
    #short shSub_sat_point_displacement_in_pixels;
    #short shSub_sat_point_displacement_in_lines;
    #char szVis_missing_line_table[316];
    #char szVisn_missing_line_table[316];
    #char szIr_missing_line_table[316];
    #char szWv_missing_line_table[316];
    #short shNo_missing_lines_replaced;
    #short shNo_black_lines_1[3],shNo_black_lines_2[3],shNo_black_lines_3[3];
    #char szWefax_annotation[28],szImage_quality_annotation[12];
    #short shSpectral_contents[3],shDeformation_used;
    #char szBb_ir_count_for_space_view[6],szStd_Bb_ir_count_for_space_view[3];
    #char szBb_wv_count_for_space_view[6],szStd_Bb_wv_count_for_space_view[3];
    #char szBb_cold_temp[5],szBb_warm_temp[5];
    #char szBb_ir_count_nominal[6],szStd_Bb_ir_count_nominal[3];
    #char szBb_wv_count_nominal[6],szStd_Bb_wv_count_nominal[3];
    #char szBb_calibration_timestamp[5];
    #char szMpef_abs_ir_cal[5],szMpef_ir_space_count[3],szMpef_ir_calibration_timestamp[5];
    #char szMpef_abs_wv_cal[5],szMpef_wv_space_count[3],szMpef_wv_calibration_timestamp[5];
    #char szSpares1[15];
    #char szAll_channel_gains[8],szSpares2[4];
    #double dsRight_ascension_att_south,dsDeclination_att_south;
    #double dsRight_ascension_att_north,dsDeclination_att_north;
    #double dsRefined_attitude_xyz[3];
    #double dsRight_ascension_declination_of_mean_att[2];
    #int iNo_slots_with_refined_attitude;
    #float flSpin_dur_minus_nominal_spin_dur;
    #char szEclipse_operation,szDecontamination,szManoeuvre,szView;
    #char szIr1_on,szIr2_on,szWv1_on,szWv2_on,szVis1_on,szVis2_on,szVis3_on,szVis4_on;
    #char szMpef_spares[20],szSpares4[16],szImage_status[16];
    #short shLine_pixel_orientation[12];
    #double dsSat_earcth_centre_distance;
    #double dsOrbit_offset_at_southern_horizon[3];
    #double dsOrbit_offset_at_northern_horizon[3];
    #float flMax_deformation_diff_x_inside_column;
    #float flMax_deformation_diff_y_inside_line;
    #float flMax_deformation_diff_x_inside_line;
    #float flMax_deformation_diff_y_inside_column;
    #char  szImage_conditions[4];
    #short shMin_count_in_histogram[4],shMax_count_in_histogram[4];
    #float flMean_vis_1_4,Mean_vis_2_3,flSnr_in_space_corners[4];
    #int iNo_lines_in_snr_calc;
    #float flSnr_eastern_part[4],flSnr_western_part[4];
    #float flMean_noise_count_eastern_part[4],flMean_noise_count_western_part[4];
    #short shMax_space_count_eastern_part[4],shMax_space_count_western_part[4];
    #char szSpares5[88],szSpares6[800];
    #short shNominal_sat_longitude;
    #char szSpares7[2];
    #float flNominal_sat_longitude_in_degrees;
    #int iCode_for_sub_sat_point,iNo_landmarks;
    #int iNo_vis_cloudfree_landmarks,iNo_ir_cloudfree_landmarks;
    #int   iGQASource;
    #int   iHeaderVersionNumber;
    #float flNo_vis_landmarks_with_corr;
    #float flNo_ir_landmarks_with_corr;
    #float flAbs_std_visible_landmarks[2];
    #float flAbs_max_visible_landmarks[2];
    #int iNo_landmarks_with_max_deviation[2];
    #float flAbs_std_ir_landmarks[2];
    #float flAbs_max_ir_landmarks[2];
    #float flAbs_landmark_results[128][8];
    #float flNo_vis_landmarks_with_curr;
    #float flNo_ir_landmarks_with_curr;
    #float flRel_std_visible_landmarks[2];
    #float flRel_max_visible_landmarks[2];
    #int iNo_landmarks_with_rel_max_deviation[2];
    #float flRel_std_ir_landmarks[2];
    #float flRel_max_ir_landmarks[2];
    #float flRel_landmark_results[128][8];
    #float flVis_equalisation_coef;
    #float flDif_direct_in_pixel[4];
    #float flDif_direct_in_line[4];
    #int iRectification_flag;
    #short shGQA_check_flag;
    #int iImage_processing_status;
    #char  szAlignPadding[2];
    #float flSpaceCornersVISS[4];
    #float flSpaceCornersVISN[4];
    #float flSpaceCornersIR[4];
    #float flSpaceCornersWV[4];
    #float flStdSpaceCornersVISS[4];
    #float flStdSpaceCornersVISN[4];
    #float flStdSpaceCornersIR[4];
    #float flStdSpaceCornersWV[4];
    #char szSpares9[4152]; 


filename="/DSNNAS/Repro/mviri/level1/HR-MFG15/data/MET7/2006/01/METEOSAT7-MVIRI-MTP15-NA-NA-20060101000000.000000000Z-1101580-1"
#filename="/DSNNAS/Repro/mviri/level1/HR-MFG15/data/MET7/2000/08/METEOSAT7-MVIRI-MTP15-NA-NA-20000819053000.000000000Z-1101882-2"      
with open(filename, "rb") as f:
  
  header(f)
    
