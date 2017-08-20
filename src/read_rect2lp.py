"""
this file contains functions reading the rect2lp file
"""

import sys,os
import struct
from collections import namedtuple
import numpy as np
import timeit
try:
  sys.path.insert(0, os.getcwd()[:-4]+'/lib/vax2ieee/')
  import _vax2ieee as vi
except ImportError:
  sys.path.insert(0, os.getcwd()[:-6]+'/lib/vax2ieee/')
  import _vax2ieee as vi

global recordsize

recordsize=15360 #records are described further down in structures()

def vax(value,typ,eness="le"):
  #floats and doubles have to be treated specially 
  #because of the vax format
  try:
    if typ=='d':
      return vi.vax2d(value)  #for double: convert from vax
    elif typ=='f':
      try:
	return vi.vax2f(value)#for float: convert from vax
      except OverflowError:
	print value
	return value
    else:
      return value#for rest: just take the value
  except TypeError:
    return -9999


def header(f,start=1,end=2504,eness="le",t=0):
#def header(f,start=1,end=3,eness="le",t=0):
    #open file
    records=range(start,end+1)#records=[1,2,3]
    f.seek((start-1)*recordsize, 0)
    recs=[]
    for r in records:
      
      #handle records>3
      rq=r
      if r==3 or r==2504:
	rq=3
      elif r>3 and r<2504:
	rq=4
	
      #get information on record structure    
      names,dims,rec_fmt=structures(rq)

      recID  = 'r'+str(r)
      struc  = namedtuple(recID, names, verbose=False)
      
      #read and unpack record
      if t==1:
	timing=[]
	timing.append(timeit.default_timer())
      rec_bin = f.read(recordsize)
      if t==1:
	timing.append(timeit.default_timer())

      record=struct.unpack(str(rec_fmt),rec_bin)
      #record=np.fromstring(rec_bin,dtype=str(rec_fmt[0]))
      if t==1:
	timing.append(timeit.default_timer())
	print "reading took: "+str(timing[1]-timing[0])+"seconds"
	print "unpacking took: "+str(timing[2]-timing[1])+"seconds"

      g=0
      lis=[]
      for i in range(0,len(names)):
	
	currname=names[i]
	currtype=currname[:1]
	
	
	#fill buffer to capture dimensional fields
	buf=[]
        
	for j in range(dims[i]):
	  buf.append(vax(record[g+j],currtype,eness))
	  
	#collect content in one list
	lis.append(buf)
	g=g+dims[i]

      
      #assign content to named tuple
      struc=struc(*lis)
      
      #store record
      recs.append(struc)
    
    return f,recs
    
def image(f,head,start=4,end=2503,eness="le",t=0):
  
    
    vistargets=['shVisSPixel','shVisNPixel']
    pixeldata=['shVisSPixel','shVisNPixel','shIRPixel','shWVPixel']
    irtargets =['shIRPixel']
    wvtargets =['shWVPixel']
    #open file
    records=range(start,end+1)#records=[1,2,3]
    
    recs=[]
    visim=[]
    irim=[]
    wvim=[]
    
    rq = 4
    names,dims,rec_fmt=structures(rq)
    
    for r in records:
      recID ='r'+str(r)
      struc = namedtuple(recID, names[:5]+names[9:], verbose=False)
      
      #read and unpack record
      if t==1:
	timing=[]
	timing.append(timeit.default_timer())
      rec_bin = f.read(recordsize)
      if t==1:
	timing.append(timeit.default_timer())
      record=struct.unpack(str(rec_fmt),rec_bin)
      
      if t==1:
	timing.append(timeit.default_timer())
	print "reading took: "+str(timing[1]-timing[0])+"seconds"
	print "unpacking took: "+str(timing[2]-timing[1])+"seconds"
      
      early=['M2','M3']
      sat=getattr(head[0], 'szSatCode')[0]
 
      #go through fields
      g=0
      lis=[]
      for i in range(0,len(names)): 
	currname=names[i]
	currtype=currname[:1]

	#fill buffer to capture dimensional fields
	buf=[]
	for j in range(dims[i]):
	  if sat in early:
	    if currname in pixeldata:
	      buf.append(int(format(vax(record[g+j],currtype,eness), 'b').zfill(8)[:6],2))#
	    else:
	      buf.append(vax(record[g+j],currtype,eness))
	  else:
	    buf.append(vax(record[g+j],currtype,eness))    
	  
	#store line headers and trailers in list...
	if not currname in pixeldata:
	  lis.append(buf)
	#...and store pixel values in image
	else:
	  if currname in vistargets:
	    visim.append(buf)
	  if currname in irtargets:
	    irim.append(buf)
	  if currname in wvtargets:
	    wvim.append(buf)
	g=g+dims[i]

      
      #assign content to named tuple
      struc=struc(*lis)

      #store record
      recs.append(struc)

    return f,recs,visim,irim,wvim
    
def structures(z):
  import re
  #define names
  all_r1=[['szOriginalFormat',                           8,         '8c'],   
          ['iYear',                                      1,         'i'],    
          ['iYday',                                      1,         'i'],    
          ['iSlotNum',                                   1,         'i'],    
          ['iDataType',                                  1,         'i'],    
	  ['iImageProcessingStatus',                     1,         'i'],    
          ['iSpectralContent',                           1,         'i'],    
          ['iNumScans',                                  1,         'i'],
          ['iNumScanValidRapid',                         1,         'i'],
	  ['iRecLength',                                 1,         'i'],    
          ['iNumRec',                                    1,         'i'],    
          ['iNumRecHead',                                1,         'i'],    
	  ['szSpares2',                                  12,        '12c'],  
          ['szFileType',                                 8,         '8c'],   # Do we really need this?
          ['szSatCode',                                  2,         '2c'],   
          ['szSpares3',                                  6,         '6c'],   
	  ['iJulianSlotNum',                             1,         'i'],   
          ['shLastLineImage',                            1,         'h'],    
          ['byVisChannnelInUse',                         4,         '4B'],   
	  ['byIrChannnelInUse',                          2,         '2B'],   
          ['byWvChannnelInUse',                          2,         '2B'],   
          ['shScanMode',                                 1,         '2h'],    
	  ['iSecsPastMidnightForSHL',                    1,         '5i'],    
          ['shVisibleSouthNorth',                        2,         'h'],   
	  ['iFishRemoval',                               5,         'i'],   
          ['iNumLineWithRadposAnomalies',                1,         'i'],    
          ['iRectInterpMethod',                          1,         'i'],    
	  ['iMixedVisModeFlag',                          1,         'i'],   
          ['ilEastWestCorrAppFlag',                      1,         'i'],   
          ['iLensRotationCorrApplFlag',                  1,         'i'],   
          ['iVariableGainCorrAppFlag',                   1,         'i'],   
          ['szSpares4',                                  12,        '12c'],  
          ['iNullDefMatFlag',                            1,         'i'],    
          ['iIPSdataVersion',                            1,         'i'],    
          ['iIPSalgoVersion',                            1,         'i'],    
          ['iIPSheadRec1Version',                        1,         'i'],    
          ['iWhereRectifiedIndicator',                   1,         'i'],    
          ['szSpares5',                                  20,        '20c'],  
          ['iFirstLineUsed',                             1,         'i'],    
          ['iLastLineUsed',                              1,         'i'],    
	  ['byVisSouthLostLineTable',                    3030,      '3030B'],
	  ['byVisNorthLostLineTable',                    3030,      '3030B'],
	  ['byIrLostLineTable',                          3030,      '3030B'],
          ['byWvLostLineTable',                          3030,      '3030B'],
	  ['iDropBitTestControlCodes',                   4,         '4i'],  
          ['iDropBitHandlingStatsViss',                  4,         '4i'],   
	  ['iDropBitHandlingStatsVisn',                  4,         '4i'],   
          ['iDropBitHandlingStatsIr',                    4,         '4i'],   
	  ['iDropBitHandlingStatsWv',                    4,         '4i'],   
          ['szSpare',                                    2952,      '2952c']]

  all_r1 = np.array(all_r1)
  str_r1 = all_r1[:, 0]
  dim_r1 = all_r1[:, 1].astype(int)
  fmt_r1 = '<' + ''.join(all_r1[:, 2])
  
  all_r2=[['iJulianSlotNum',                                    1,            'i'],      
          ['iGeometricQuality',                                 1,            'i'],      
          ['iNoGridPointsDefMat',                               1,            'i'],      
	  ['iDeformStarts',                                     1,            'i'],      
          ['iDeformEnds',                                       1,            'i'],      
          ['iDeformStep',                                       1,            'i'],      
	  ['flXDeformationMatrix',                            676,          '676f'], # These are always set to zero?
          ['flYDeformationMatrix',                            676,          '676f'], # These are always set to zero?  
	  ['szSpares1',                                      5412,         '5412c'],  
          ['szWefaxSatId',                                      4,            '4c'],     
          ['szHrSatId',                                         2,            '2c'],     
	  ['szSpares2',                                         2,            '2c'],     
          ['flVis1RegParamConst',                               2,            '2f'],     
	  ['flVis2RegParamConst',                               2,            '2f'],     
          ['flVis3RegParamConst',                               2,            '2f'],     
	  ['flVis4RegParamConst',                               2,            '2f'],     
          ['flIr1RegParamConst',                                2,            '2f'],     
	  ['flIr2RegParamConst',                                2,            '2f'],     
          ['flWv1RegParamConst',                                2,            '2f'],     
	  ['flWv2RegParamConst',                                2,            '2f'],     
          ['szSpares3',                                        64,           '64c'],    
	  ['flVis1RegParamLinear',                              2,            '2f'],     
          ['flVis2RegParamLinear',                              2,            '2f'],     
	  ['flVis3RegParamLinear',                              2,            '2f'],     
          ['flVis4RegParamLinear',                              2,            '2f'],     
	  ['flIr1RegParamLinear',                               2,            '2f'],     
          ['flIr2RegParamLinear',                               2,            '2f'],     
	  ['flWv1RegParamLinear',                               2,            '2f'],     
          ['flWv2RegParamLinear',                               2,            '2f'],     
	  ['szSpares4a',                                       64,           '64c'],    
          ['flSeasonaalAmplitude',                              1,            'f'],      
          ['flSeasonalPhase',                                   1,            'f'],      
          ['flDefMataFineTuneCorrYdir',                         1,            'f'],      
          ['flDefMatFineTuneCorrXdir',                          1,            'f'],      
          ['szSpares5',                                        76,           '76c'],    
          ['dsEarthEquatorialRadius',                           1,            'd'],      
          ['dsEarthPolarRadius',                                1,            'd'],      
	  ['dsDistFromEarthCentre',                             1,            'd'],      
	  ['dsRadiometerSteppinAngle',                          1,            'd'],      
          ['dsNominalSubSatelliteLongitude',                    1,            'd'],      
	  ['szCodeForNominalSubSatellitePoint',                 4,            '4c'],     
          ['szNullDefMatUsedFlag',                              4,            '4c'],     
          ['dsAlpha0',                                          1,            'd'],     
          ['dsSat_OffsetFromSSP',                               1,            'd'],      
          ['dsAbsPositionSat',                                  1,            'd'],      
          ['szSpares6',                                         8,            '8c'],     
          ['dsSICClockPeriod',                                  1,            'd'],      
          ['dsTimeDurationPixel',                               1,            'd'],      
          ['szSpare7',                                         16,           '16c'],    
          ['dsEquatorialRadEastAtmHeight',                      1,            'd'],      
          ['dsEquatorialRadWestAtmHeight',                      1,            'd'],      
          ['dsEquatorialRadSouthAtmHeight',                     1,            'd'],      
          ['dsEquatorialRadNorthAtmHeight',                     1,            'd'],      
          ['szSpares8',                                       260,          '260c'],   
          ['flRefinedSouthHorizRADPOS',                         1,            'f'],      
          ['flRefinedNorthHorizRADPOS',                         1,            'f'],      
          ['iRADPosLineIDdiff',                                 1,            'i'],      
          ['dsTimeRADPosRelationCoeff',                         2,            '2d'],     
          ['dsF0ScanningParam',                                 1,            'd'],      
          ['dsF1ScanningParam',                                 1,            'd'],      
          ['dsAttitude',                                        3,            '3d'],     
          ['dsPixelSamplingRateCoeff',                          3,            '3d'],     
          ['dsSpinDeviationCoeff',                              3,            '3d'],     
          ['dsRawImageSkewnessCoeff',                           3,            '3f'],     
          ['szSpare8',                                          4,            '4c'],     
          ['dsAttitudeAttDriftRate',                            4,            '4d'],     
          ['dsEpochAMVG3',                                      1,            'd'],      
          ['dsF2ScanningParam',                                 1,            'd'],      
          ['dsF3ScanningParam',                                 1,            'd'],      
          ['szSpares9',                                       208,          '208c'],   
          ['dsOrbitCoordinatesMgfImageStart',                   6,            '6d'],     
	  ['dsOrbitCoordinatesMgfImageEnd',                     6,            '6d'],     
	  ['flCartCompAttitudeImageStart',                      3,            '3f'],     
	  ['flCartCompAttitudeImageEnd',                        3,            '3f'],     
	  ['dsOrbitCoordinatesFixedEarthImageStart',            6,            '6d'],     
	  ['dsOrbitCoordinatesFixedEarthImageEnd',              6,            '6d'],     
	  ['flEarthCentreInPixels',                             1,            'f'],     
          ['flEarthCentreInLines',                              1,            'f'],      
	  ['iProgramErrorLevelStopCode',                        1,            'i'],      
          ['szSpares10',                                      272,          '272c'],   
          ['iLensCoeffType',                                    1,            'i'],      
          ['szSpares11',                                       32,           '32c'],    
          ['dsSunCoordMgfImage_start',                          6,            '6d'],    
          ['dsSunCoordMgfImage_end',                            6,            '6d'],     
	  ['iNoVisSouthHorizPairs',                             1,            'i'],     
          ['iNoVisSouthHorizPairsQuality',                      1,            'i'],      
	  ['iNoVisNorthHorizPairs',                             1,            'i'],      
          ['iNoVisNorthHorizPairsQuality',                      1,            'i'],      
	  ['iNoIrHorizPairs',                                   1,            'i'],      
          ['iNoIrHorizPairsQuality',                            1,            'i'],      
	  ['szSpares12',                                       44,           '44c'],    
          ['dsEastWestCenteringBiasCoeff',                    3,            '3d'],     
          ['dsEquatorialAtmosphericHeight',                    1,            'd'],     
          ['dsTHETA1ParamFromB0Wraparound',                     1,            'd'],      
          ['dsEastWestLensCorrSlotEnd',                         1,            'd'],      
          ['dsNorthSouthLensCorrSlotEnd',                       1,            'd'],      
	  ['byStableAttitudeAvailFlag',                         4,            '4B'],     
          ['byEeastWestBiasCoeffAvailThisSlotFlag',             4,            '4B'],     
          ['byCalibratedTH1G3AvailFlag',                        4,            '4B'],    
          ['byCalibratedT1HG3AvailFlag',                        4,            '4B'],     
          ['dsTHETA1ParamFromHousekeeping',                     1,            'd'],
          ['szSpares13',                                     2580,         '2580c']]    
  
  all_r2 = np.array(all_r2)
  str_r2 = all_r2[:, 0]
  dim_r2 = all_r2[:, 1].astype(int)
  fmt_r2 = '<' + ''.join(all_r2[:, 2])
  
  def get_num(x):
    return int(filter(str.isdigit), x)

  #  abc = []  
  # for s in all_r2[:, 2]:
  #  if len(s) > 1:
  #    abc.append(int(filter(str.isdigit, s)))
  #  else:
  #   abc.append(1)

  #Header 3
  all_r3=[['iJulianSlotNum',                                           1,        'i'],    
          ['shJulianDay',                                              1,        'h'],    
          ['shSlotNum',                                                1,        'h'],    
          ['shNominalImageTime',                                       1,        'h'],    
	  ['shSlotType',                                               1,        'h'],    
          ['shImageQualityFlag',                                       1,        'h'],    
          ['szSatId',                                                  2,        '2c'],   
	  ['shSubSatPointDisplacementInPixels',                        1,        'h'],    
          ['shSubSatPointDisplacementInLines',                         1,        'h'],    
	  ['byVissMissingLineTable',                                  316,      '316B'], 
          ['byVisn_missingLineTable',                                316,      '316B'], 
	  ['byIrMissingLineTable',                                   316,      '316B'], 
          ['byWvMissingLineTable',                                   316,      '316B'], 
	  ['shNumMissingLinesReplaced',                                 1,        'h3'],   
          ['shNumBlackLines1',                                         3,        'h3'],   
          ['shNumBlackLines2',                                         3,        'h3'],   
          ['shNumBlackLines3',                                         3,        'h2'],   
	  ['szWefaxAnnotation',                                      1,        '8s'],   
          ['szImageQualityAnnotation',                               1,        '12s'],  
	  ['shSpectralContent',                                      3,        '3h'],   
          ['shDeformationUsed',                                        1,        'h'],   
	  ['szBbIrCountForSpaceView',                                 6,        '6c'],   
          ['szStdBbIrCountForSpaceview',                              3,        '3c'],   
	  ['szBbWvCountForSpaceview',                                 6,        '6c'],   
          ['szStdBbWvCountForSpaceView',                              3,        '3c'],   
	  ['szBbColdTemp',                                            5,        '5c'],   
          ['szBbWarmTemp',                                            5,        '5c'],   
          ['szBbIrCountNominal',                                      6,        '6c'],   
          ['szStdBbIrCountNominal',                                   3,        '3c'],   
	  ['szBbWvCountNominal',                                      6,        '6c'],   
          ['szStdBbWvCountNominal',                                   3,        '3c'],   
          ['szBbCalibrationTimestamp',                                5,        '5c'],   
	  ['szMpefAabsIrCal',                                         5,        '5c'],   
          ['szMpefIrSpaceCount',                                      3,        '3c'],   
          ['szMpefIrCalibrationTimestamp',                            5,        '5c'],   
	  ['szMpefAbsWvCal',                                          5,        '5c'],   
          ['szMpefWvSpaceCount',                                      3,        '3c'],   
          ['szMpefWvCalibrationTimestamp',                            5,        '5c'],   
	  ['szSpares14',                                              15,       '15c'],  
          ['szAllChannelGains',                                       8,        '8c'],   
          ['szSpares15',                                               4,        '4c'],   
	  ['dsRightAscensionAattSouth',                                1,        'd'],    
          ['dsDeclinationAttSouth',                                    1,        'd'],    
          ['dsRightAscensionAttNorth',                                 1,        'd'],    
          ['dsDeclinationAttNorth',                                    1,        'd'],    
	  ['dsRefinedAttitudeXYZ',                                     3,        '3d'],   
          ['dsRightAscensionDeclinationMeanAtt',                      2,        '2d'],   
	  ['iNoSlotsWithRefinedAttitude',                              1,        'i'],    
          ['flSpinDurMinusNominalSpinDur',                             1,        'f'],    
	  ['byEclipseOperation',                                       1,        'B'],    
          ['byDecontamination',                                        1,        'B'],    
          ['byManoeuvre',                                              1,        'B'],    
          ['byView',                                                   1,        'B'],    
          ['byIr1On',                                                  1,        'B'],    
          ['byIr2On',                                                  1,        'B'],    
          ['byWv1On',                                                  1,        'B'],    
          ['byWv2On',                                                  1,        'B'],    
          ['byVis1On',                                                 1,        'B'],    
          ['byVis2On',                                                 1,        'B'],    
          ['byVis3On',                                                 1,        'B'],    
          ['byVis4On',                                                 1,        'c'],    
	  ['szsparesMpef16',                                          20,       '20c'],  
          ['szSpares17',                                              16,       '16c'],  
          ['byImageStatus',                                           16,       '16B'],  
          ['shLinePixelOrientation',                                  12,       '12h'],  
	  ['dsSatEarcthCentreDistance',                                1,        'd'],    
          ['dsOrbitOffsetSouthernHorizon',                             3,        '3d'],   
	  ['dsOrbitOffsetNorthernHorizon',                             3,        '3d'],   
          ['flMaxDeformationDiffXInsideColumn',                       1,        'f'],    
	  ['flMaxDeformationDiffYInsideLine',                         1,        'f'],    
          ['flMaxDeformationDiffXInsideLine',                         1,        'f'],    
	  ['flMaxDeformationDiffYInsideColumn',                       1,        'f'],    
          ['byImageConditions',                                       4,        '4B'],   
	  ['shMinCountValue',                                         4,        '4h'],   
          ['shMaxCountValue',                                         4,        '4h'],   
          ['flMeanVis1or4',                                           1,        'f'],    
          ['flMeanVis2or3',                                           1,        'f'],    
          ['flSnrSpaceCorners',                                        4,        '4f'],   
	  ['iNumLinesSnrCalc',                                         1,        'i'],    
          ['flSnrEastern',                                             4,        '4f'],   
          ['flSnrWestern',                                             4,        '4f'],   
          ['flMeanNoiseCountEastern',                                  4,        '4f'],   
          ['flMeanNoiseCountWestern',                                  4,        '4f'],   
	  ['shMaxSpaceCountEastern',                                   4,        '4h'],   
          ['shMaxSpaceCountWestern',                                   4,        '4h'],   
          ['szSpares18',                                              88,       '88c'],  
          ['szMessage',                                              800,      '800c'], 
	  ['shNominalSatLongitude',                                    1,        'h'],    
          ['szSpares19',                                               2,        '2c'],   
          ['flNominalSatLongitudeDeg',                                 1,        'f'],   
	  ['szCodeSSP',                                                4,        '4c'],   
          ['iNumLandmarks',                                            1,        'i'],    
          ['iNumVisCloudFreeLandmarks',                                1,        'i'],    
          ['iNumIrCloudfreeLandmarks',                                 1,        'i'],    
	  ['iGQASource',                                               1,        'i'],    
          ['iHeaderVersionNumber',                                     1,        'i'],    
          ['flNumVisLandmarksWithGoodCorr',                                1,        'f'],    
          ['flNumIrLandmarksWithGoodCorr',                                 1,        'f'],    
	  ['flAbsStdVisLandmarks',                                     2,        '2f'],   
          ['flAbsMaxVisLandmarks',                                     2,        '2f'],   
          ['iNumLandmarksWithMaxDeviation',                            2,        '2i'],   
	  ['flAbsStdIrLandmarks',                                      2,        '2f'],   
          ['flAbsMaxIrLandmarks',                                      2,        '2f'],   
	  ['flAbsLandmarkResults',                                 1024,     '1024f'],
          ['flNumVisLandmarksWithGoodCorr2',                            1,        'f'],    
	  ['flNumIrLandmarksWithGoodCorr2',                             1,        'f'],    
          ['flRelStdVisLandmarks',                                     2,        '2f'],   
	  ['flRelMaxVisLandmarks',                                     2,        '2f'],   
          ['iNumLandmarksWithRelMaxDeviation',                         2,        '2i'],   
	  ['flRelStdIrLandmarks',                                      2,        '2f'],   
          ['flRelMaxIrLandmarks',                                      2,        '2f'],   
          ['flRelLandmarkResults',                                  1024,     '1024f'],
	  ['flVisEqualisationCoef',                                    1,        'f'],    
          ['flDifInPixelDirection',                                    4,        '4f'],   
          ['flDifInLineDirection',                                     4,        '4f'],   
	  ['iRectificationFlag',                                       1,        'i'],   
          ['shGqaCheckFlag',                                           1,        'h'],    
          ['iImageProcessingStatus',                                   1,        'i'],    
	  ['szAlignPadding',                                           2,        '2c'],   
          ['flSpaceCornersVISS',                                       4,        '4f'],   
          ['flSpaceCornersVISN',                                       4,        '4f'],   
	  ['flSpaceCornersIR',                                         4,        '4f'],   
          ['flSpaceCornersWV',                                         4,        '4f'],   
          ['flStdSpaceCornersVISS',                                    4,        '4f'],   
          ['flStdSpaceCornersVISN',                                    4,        '4f'],  
	  ['flStdSpaceCornersIR',                                      4,        '4f'],   
          ['flStdSpaceCornersWV',                                      4,        '4f'],   
          ['dsQuadraticEqCoeff',                                       3,        '3d'],   
          ['szM6CCalibData',                                          26,       '26c'],  
          ['szSpares20',                                            4102,     '4102c']]

  all_r3 = np.array(all_r3)
  str_r3 = all_r3[:, 0]
  dim_r3 = all_r3[:, 1].astype(int)
  fmt_r3 = '<' + ''.join(all_r3[:, 2])
  
  # Record 4
  all_r4=[['iJulianSlotNum',                        1,        'i'],     
          ['shLineNr',                           1,        'h'],     
          ['shChannelIdVIS',                     1,        'h'],     
          ['shChannelIdIR',                      1,        'h'],     
          ['shChannelIdWV',                      1,        'h'],     
	  ['shVisSPixel',                        5000,     '5000B'], 
          ['shVisNPixel',                        5000,     '5000B'], 
          ['shIRPixel',                          2500,     '2500B'], 
          ['shWVPixel',                          2500,     '2500B'], 
	  ['shRectLineStatus',                   1,        'h'],     
          ['shLineStatus',                       1,        'h'],     
	  ['shTimefitFlag',                      1,        'h'],     
          ['shFirstRectPix',                     1,        'h'],     
          ['dsTimeFirstRectPix',                  1,        'd'],     
          ['shLastRectPix',                      1,        'h'],     
	  ['dsTimeLastRectPix',                   1,       'd'],     
          ['shFirstFitPix',                      1,        'h'],     
          ['shLastFitPix',                       1,        'h'],     
          ['dsMaxDevFit',                         1,        'd'],     
	  ['dsFitCoef',                           5,       '5d'],    
          ['shFirstTime',                        1,        'h'],     
          ['shInterval',                         1,        'h'],     
          ['flTimesArray',                       49,       '49f'],   
          ['szSpares21',                         70,       '70c']]   

  
  all_r4 = np.array(all_r4)
  str_r4 = []
  for item  in all_r4[:, 0]:
    str_r4.append(item)

  dim_r4 = all_r4[:, 1].astype(int)
  fmt_r4 = '<' + ''.join(all_r4[:, 2])
          
  # Unpacking descriptions can be found here:
  # https://docs.python.org/2/library/struct.html
          
  rec_fmt=[fmt_r1, fmt_r2, fmt_r3, fmt_r4]
  
  #use record specific data
  if z==1:
    str_r=str_r1
    dim_r=dim_r1
  elif z==2:
    str_r=str_r2
    dim_r=dim_r2
  elif z==3: #here automatically also 2504 will be handled
    str_r=str_r3
    dim_r=dim_r3
  elif z==4: #here automatically also 4-2503 will be handled
    str_r=str_r4
    dim_r=dim_r4

  #replace non-alphanumeric values
  i=0
  for stri in str_r:
    str_r[i]=re.sub(r'\W+', '',str_r[i])
    i=i+1
        
  #return
  return (str_r,dim_r,rec_fmt[z-1])

def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return np.isnan(y), lambda z: z.nonzero()[0]

def deform(head):
  n=int(getattr(head[1], 'iNo_grid_points_in_def_matrix')[0])
  print n
  n=26
  allX = np.reshape(np.array(getattr(head[1], 'flX_deform_matrix')),(n,n)) #seem all zero
  allY = np.reshape(np.array(getattr(head[1], 'flY_deform_matrix')),(n,n)) #seem all zero
  
  
  print allX.shape
  print allY.shape
  
  #interpol
  nans, x= nan_helper(allX)
  allX[nans]= np.interp(x(nans), x(~nans), allX[~nans])
  nans, x= nan_helper(allY)
  allY[nans]= np.interp(x(nans), x(~nans), allY[~nans])

  
  #make 2d
  dX=allX.reshape([n,n])
  dY=allY.reshape([n,n])
  
  
  return dX,dY
  
