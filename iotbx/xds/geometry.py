#!/usr/bin/env libtbx.python
#
# iotbx.xds.geometry.py
#
# Code to read experimental geometry from a number of XDS output files
# (XPARM, reflection files) returning geometry object tree in the existing
# coordinate frame. Use method structured_xds_geometry(filename)
#
# Aim of code: to save writing this over and over...
#
# Graeme Winter, Diamond Light Source, 2012/OCT/16

from __future__ import absolute_import, division, print_function
from six.moves import range

class bucket:
  pass

def structured_xds_geometry(filename):
  '''
  Structured view of the geometry stored in an XDS data or metadata file.
  '''

  # read the geometric information from the file

  if open(filename, 'r').read(26) == '!OUTPUT_FILE=INTEGRATE.HKL':
    _osc, _det, _ub, _pix, _beam = read_geometry_integrate_hkl(filename)[:5]
  elif open(filename, 'r').read(17) == '!FORMAT=XDS_ASCII':
    _osc, _det, _ub, _pix, _beam = read_geometry_integrate_hkl(filename)[:5]
  else:
    _osc, _det, _ub, _pix, _beam = read_geometry_xparm_xds(filename)[:5]

  # create some data structures for this
  oscillation = bucket()

  # FIXME I should remove the need from frame0 here - recalculate phi0
  oscillation.phi0 = _osc['phi0'] - (_osc['frame0'] - 1) * _osc['dphi']
  oscillation.dphi = _osc['dphi']
  oscillation.axis = _osc['axis']

  beam = bucket()
  beam.direction = _beam['direction'].normalize()
  beam.wavelength = _beam['wavelength']

  detector = bucket()
  detector.fast = _det['fast']
  detector.slow = _det['slow']
  detector.normal = _det['normal']
  detector.origin = _det['origin']
  detector.size = _pix['nx'], _pix['ny']
  detector.pixel_size = _pix['qx'], _pix['qy']

  sample = bucket()
  sample.ub = _ub['ub']
  sample.cell = (_ub['a'].length(), _ub['b'].length(), _ub['c'].length(),
           _ub['b'].angle(_ub['c'], deg = True),
           _ub['c'].angle(_ub['a'], deg = True),
           _ub['a'].angle(_ub['b'], deg = True))

  # gather everything together
  geometry = bucket()
  geometry.oscillation = oscillation
  geometry.beam = beam
  geometry.detector = detector
  geometry.sample = sample

  return geometry

def read_geometry_xparm_xds(filename):
  '''
  Read XDS XPARM file from indexing or global refinement and generate a model
  of the geometry. Returns object using cctbx matrices.
  '''

  from scitbx import matrix

  oscillation = { }
  detector = { }
  ub = { }
  pixel = { }
  beam = { }
  variance = { }

  # read
  xparm_values = open(filename, 'r').read().split()

  # validate
  assert(len(xparm_values) == 42)

  # understand
  xparm_values = [float(x) for x in xparm_values]
  oscillation['frame0'] = int(xparm_values[0])
  oscillation['phi0'] = xparm_values[1]
  oscillation['dphi'] = xparm_values[2]
  oscillation['axis'] = matrix.col(xparm_values[3:6])
  beam['wavelength'] = xparm_values[6]
  beam['direction'] = matrix.col(xparm_values[7:10])
  pixel['nx'] = int(xparm_values[10])
  pixel['ny'] = int(xparm_values[11])
  pixel['qx'] = xparm_values[12]
  pixel['qy'] = xparm_values[13]
  detector['distance'] = xparm_values[14]
  detector['orgx'] = xparm_values[15] * pixel['qx']
  detector['orgy'] = xparm_values[16] * pixel['qy']
  detector['fast'] = matrix.col(xparm_values[17:20])
  detector['slow'] = matrix.col(xparm_values[20:23])
  detector['normal'] = matrix.col(xparm_values[23:26])
  detector['origin'] = (detector['distance'] * detector['normal'] -
              detector['orgx'] * detector['fast'] -
              detector['orgy'] * detector['slow'])
  ub['a'] = matrix.col(xparm_values[33:36])
  ub['b'] = matrix.col(xparm_values[36:39])
  ub['c'] = matrix.col(xparm_values[39:42])
  ub['ub'] = matrix.sqr(ub['a'].elems + ub['b'].elems +
              ub['c'].elems).inverse()

  return oscillation, detector, ub, pixel, beam

def read_geometry_integrate_hkl(filename):
  '''
  Read XDS INTEGRATE.HKL or XDS_ASCII.HKL from integrate or correct steps,
  return model of geometry. Returns object using cctbx matrices.
  '''

  from scitbx import matrix

  oscillation = { }
  detector = { }
  ub = { }
  pixel = { }
  beam = { }

  for record in open(filename):
    if not record.startswith('!'):
      break

    tokens = record[1:].split()

    if tokens[0] == 'NX=':
      for j in range(0, 8, 2):
        pixel[tokens[j].replace('=', '').lower()] = float(
          tokens[j + 1])
    elif tokens[0] == 'STARTING_FRAME=':
      oscillation['frame0'] = int(tokens[1])
    elif tokens[0] == 'STARTING_ANGLE=':
      oscillation['phi0'] = float(tokens[1])
    elif tokens[0] == 'OSCILLATION_RANGE=':
      oscillation['dphi'] = float(tokens[1])
    elif tokens[0] == 'INCIDENT_BEAM_DIRECTION=':
      # TODO verify if matrix.col accepts iterator in both py2/3
      beam['direction'] = matrix.col( [float(t) for t in tokens[-3:]])
      beam['wavelength'] = 1.0 / beam['direction'].length()
    elif tokens[0] == 'ROTATION_AXIS=':
      oscillation['axis'] = matrix.col( [float(t) for t in tokens[-3:]])
    elif tokens[0] == 'DIRECTION_OF_DETECTOR_X-AXIS=':
      detector['fast'] = matrix.col([float(t) for t in tokens[-3:]])
    elif tokens[0] == 'DIRECTION_OF_DETECTOR_Y-AXIS=':
      detector['slow'] = matrix.col([float(t) for t in tokens[-3:]])
    elif tokens[0] == 'ORGX=':
      orgx, orgy = float(tokens[1]), float(tokens[3])
      detector['_org'] = orgx, orgy
    elif tokens[0] == 'DETECTOR_DISTANCE=':
      detector['distance'] = float(tokens[1])
    elif tokens[0] == 'UNIT_CELL_A-AXIS=':
      ub['a'] = matrix.col([float(t) for t in tokens[-3:]])
    elif tokens[0] == 'UNIT_CELL_B-AXIS=':
      ub['b'] = matrix.col([float(t) for t in tokens[-3:]])
    elif tokens[0] == 'UNIT_CELL_C-AXIS=':
      ub['c'] = matrix.col([float(t) for t in tokens[-3:]])

  # postprocess results
  orgx, orgy = detector['_org']
  detector['normal'] = detector['fast'].cross(detector['slow'])
  detector['orgx'], detector['orgy'] = (orgx * pixel['qx'],
                      orgy * pixel['qy'])
  detector['origin'] = (detector['distance'] * detector['normal'] -
              detector['orgx'] * detector['fast'] -
              detector['orgy'] * detector['slow'])
  ub['ub'] = matrix.sqr(ub['a'].elems + ub['b'].elems +
              ub['c'].elems).inverse()

  return oscillation, detector, ub, pixel, beam

# code here to allow test cases - static example texts

__integrate_text = '''!OUTPUT_FILE=INTEGRATE.HKL    DATE=11-Oct-2012
!Generated by INTEGRATE    (VERSION  March 15, 2012)
!PROFILE_FITTING= TRUE
!SPACE_GROUP_NUMBER=   75
!UNIT_CELL_CONSTANTS=    57.965    57.965   150.135  90.000  90.000  90.000
!NAME_TEMPLATE_OF_DATA_FRAMES=/Users/graeme/data/i03-setup-110711/x1/thau_2_????
!DETECTOR=PILATUS
!NX=  2463  NY=  2527    QX=  0.172000  QY=  0.172000
!STARTING_FRAME=       1
!STARTING_ANGLE=     0.000
!OSCILLATION_RANGE=  0.100000
!ROTATION_AXIS=  0.999999 -0.001342 -0.000955
!X-RAY_WAVELENGTH=  0.976300
!INCIDENT_BEAM_DIRECTION=  0.001787  0.003138  1.024269
!DIRECTION_OF_DETECTOR_X-AXIS=  1.000000  0.000000  0.000000
!DIRECTION_OF_DETECTOR_Y-AXIS=  0.000000  1.000000  0.000000
!ORGX=   1220.98  ORGY=   1263.78
!DETECTOR_DISTANCE=   268.288
!UNIT_CELL_A-AXIS=    41.592    10.019    39.111
!UNIT_CELL_B-AXIS=   -13.760    56.308     0.209
!UNIT_CELL_C-AXIS=   -98.311   -24.435   110.808
!BEAM_DIVERGENCE_E.S.D.=     0.027
!REFLECTING_RANGE_E.S.D.=     0.059
!MINPK=  75.00000
!CUT=    2.00
!VARIANCE_MODEL=  4.000E+00  1.000E-04
!NUMBER_OF_ITEMS_IN_EACH_DATA_RECORD=20
!H,K,L,IOBS,SIGMA,XCAL,YCAL,ZCAL,RLP,PEAK,CORR,MAXC,
!             XOBS,YOBS,ZOBS,ALF0,BET0,ALF1,BET1,PSI
!Items are separated by a blank and can be read in free-format
!END_OF_HEADER
 -48 -5 -9 1.240E+01 2.380E+01 8.8 2512.5 777.2 0.73498 99.99931 9 5 9.0 2511.7 777.3 60.34 0.20 134.15 48.13 -21.18
 -48 -4 -10 4.495E+01 2.502E+01 10.9 2508.8 760.9 0.73320 100.00000 21 5 12.3 2508.4 761.1 60.34 0.20 134.18 48.06 -19.17
 -48 -3 -11 3.779E+01 2.556E+01 12.5 2506.9 744.6 0.73225 100.00000 25 4 13.2 2506.3 744.7 60.34 0.20 134.19 48.02 -17.22
 -48 -2 -12 1.510E+02 3.357E+01 13.4 2506.3 728.2 0.73204 100.00000 51 9 13.1 2506.3 727.8 60.34 0.20 134.18 48.00 -15.32
 -48 -2 -11 3.612E+01 2.472E+01 4.4 2494.2 728.3 0.72505 82.11671 24 6 6.4 2493.7 728.7 60.34 0.20 134.68 47.97 -15.59
 -48 -1 -13 7.332E+00 2.338E+01 13.8 2507.4 711.8 0.73262 99.99820 11 6 14.8 2506.6 710.7 60.34 0.20 134.15 48.01 -13.46
 -48 -1 -12 6.146E+00 2.490E+01 4.8 2495.0 711.8 0.72553 87.16513 21 9 5.4 2493.8 712.3 60.34 0.20 134.65 47.97 -13.71
 -48 0 -14 6.579E+01 2.778E+01 13.5 2510.0 695.3 0.73400 99.99902 29 9 13.2 2508.7 694.7 60.34 0.20 134.09 48.05 -11.66
 -48 0 -13 -1.972E+00 2.098E+01 4.6 2497.5 695.2 0.72682 84.58707 -14 5 5.7 2494.4 695.0 60.34 0.20 134.59 48.00 -11.89
 '''

__correct_text = '''!FORMAT=XDS_ASCII    MERGE=FALSE    FRIEDEL'S_LAW=TRUE
!OUTPUT_FILE=XDS_ASCII.HKL        DATE=11-Oct-2012
!Generated by CORRECT   (VERSION  March 15, 2012)
!PROFILE_FITTING= TRUE
!NAME_TEMPLATE_OF_DATA_FRAMES=/Users/graeme/data/i03-setup-110711/x1/thau_2_????
!DATA_RANGE=       1    1800
!ROTATION_AXIS=  0.999999 -0.001347 -0.000967
!OSCILLATION_RANGE=  0.100000
!STARTING_ANGLE=     0.000
!STARTING_FRAME=       1
!INCLUDE_RESOLUTION_RANGE=    68.140     1.189
!SPACE_GROUP_NUMBER=   89
!UNIT_CELL_CONSTANTS=    57.967    57.967   150.137  90.000  90.000  90.000
!UNIT_CELL_A-AXIS=    41.594    10.018    39.112
!UNIT_CELL_B-AXIS=   -13.760    56.310     0.209
!UNIT_CELL_C-AXIS=   -98.311   -24.435   110.810
!REFLECTING_RANGE_E.S.D.=     0.059
!BEAM_DIVERGENCE_E.S.D.=     0.027
!X-RAY_WAVELENGTH=  0.976300
!INCIDENT_BEAM_DIRECTION=  0.001783  0.003139  1.024269
!FRACTION_OF_POLARIZATION=   0.990
!POLARIZATION_PLANE_NORMAL=  0.000000  1.000000  0.000000
!AIR=  0.001000
!SILICON=  3.670924
!SENSOR_THICKNESS=  0.320000
!DETECTOR=PILATUS
!OVERLOAD=   1048500
!DIRECTION_OF_DETECTOR_X-AXIS=   1.00000   0.00000   0.00000
!DIRECTION_OF_DETECTOR_Y-AXIS=   0.00000   1.00000   0.00000
!DETECTOR_DISTANCE=   268.292
!ORGX=   1220.99  ORGY=   1263.77
!NX=  2463  NY=  2527    QX=  0.172000  QY=  0.172000
!VARIANCE_MODEL=  1.761E+00  6.348E-04
!NUMBER_OF_ITEMS_IN_EACH_DATA_RECORD=12
!ITEM_H=1
!ITEM_K=2
!ITEM_L=3
!ITEM_IOBS=4
!ITEM_SIGMA(IOBS)=5
!ITEM_XD=6
!ITEM_YD=7
!ITEM_ZD=8
!ITEM_RLP=9
!ITEM_PEAK=10
!ITEM_CORR=11
!ITEM_PSI=12
!END_OF_HEADER
     0     0     3  3.094E+00  6.167E-01  1203.7  1245.6    780.2 0.01474 100  49  -27.01
     0     0    -4  1.658E+04  5.552E+02  1250.3  1299.2    762.9 0.01965 100  53   25.88
     0     0     4  1.609E+04  5.387E+02  1197.0  1238.0    782.6 0.01965 100  55  -27.18
     0     0    -5  1.144E+00  7.676E-01  1257.0  1306.9    760.4 0.02457 100  39   25.72
     0     0     5  4.891E-02  7.921E-01  1190.4  1230.3    785.1 0.02457 100  10  -27.34
     0     0    -6  6.769E-01  9.447E-01  1263.7  1314.5    758.0 0.02949 100  29   25.56
     0     0     6 -6.594E-01  8.800E-01  1183.7  1222.7    787.6 0.02949 100 -31  -27.50
     0     0    -7  9.411E-01  9.236E-01  1270.3  1322.2    755.5 0.03441 100  25   25.40
     0     0     7  2.855E-01  1.032E+00  1177.0  1215.0    790.0 0.03441 100  27  -27.66
     0     0    -8  2.716E+04 -9.096E+02  1277.0  1329.9    753.0 0.03932 100  62   25.24
     0     0     8  3.547E+04 -1.187E+03  1170.4  1207.3    792.5 0.03932 100  64  -27.82
     0     0    -9  3.287E+00  1.238E+00  1283.7  1337.6    750.5 0.04424 100  40   25.07
     0     0     9  2.630E+00  1.346E+00  1163.7  1199.7    795.0 0.04425 100  32  -27.98
'''

__xparm_text = '''     1        0.0000    0.1000  0.999999 -0.001347 -0.000967
       0.976300       0.001783       0.003139       1.024269
      2463      2527    0.172000    0.172000
     268.291748    1220.988037    1263.773804
       1.000000       0.000000       0.000000
       0.000000       1.000000       0.000000
       0.000000       0.000000       1.000000
    89     57.9672     57.9672    150.1366  90.000  90.000  90.000
      41.594490      10.018485      39.111748
     -13.759865      56.310059       0.209473
     -98.310692     -24.435324     110.810379
'''

def work():
  '''
  Test out the code above.
  '''

  import os
  import tempfile, math

  cell_ref = (57.97, 57.97, 150.14, 90.00, 90.00, 90.00)

  for text in __integrate_text, __correct_text, __xparm_text:
    fd, filename = tempfile.mkstemp()
    f = os.fdopen(fd, 'w')
    f.write(text)
    f.close()
    geometry = structured_xds_geometry(filename)
    cell = geometry.sample.cell
    for j in range(6):
      assert(math.fabs(cell[j] - cell_ref[j]) < 0.1)
    os.remove(filename)

  print('OK')

if __name__ == '__main__':
  work()
