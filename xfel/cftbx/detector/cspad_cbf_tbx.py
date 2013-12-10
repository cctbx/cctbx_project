from __future__ import division
#-*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# $Id: cspad_cbf_tbx.py
#

import pycbf, os
from scitbx import matrix

# need to define these here since it not defined in SLAC's metrology definitions
asic_dimension = (194,185)
asic_gap = 3

class cbf_wrapper(pycbf.cbf_handle_struct):
  """ Wrapper class that provids convience functions for working with cbflib"""

  def add_category(self, name, columns):
    """ Create a new category and populate it with column names """
    self.new_category(name)
    for column in columns:
      self.new_column(column)

  def add_row(self, data):
    """ Add a row to the current category.  If data contains more entries than
      there are columns in this category, then the remainder is truncated
      Use '.' for an empty value in a row. """
    self.new_row()
    self.rewind_column()
    for item in data:
      self.set_value(item)
      if item == '.':
        self.set_typeofvalue("null")
      try:
        self.next_column()
      except Exception:
        break

  def add_frame_shift(self, basis, axis_settings):
    """Add an axis representing a frame shift (a rotation axis with an offset)"""
    angle, axis = angle_and_axis(basis)

    if angle == 0:
      axis = (0,0,1)

    self.add_row([basis.axis_name,"rotation","detector",basis.depends_on,
                  str(axis[0]),str(axis[1]),str(axis[2]),
                  str(basis.translation[0]),
                  str(basis.translation[1]),
                  str(basis.translation[2]),
                  basis.equipment_component])

    axis_settings.append([basis.axis_name, "FRAME1", str(angle), "0"])

def angle_and_axis(basis):
  """Normalize a quarternion and return the angle and axis
  @param params metrology object"""
  q = matrix.col(basis.orientation).normalize()
  return q.unit_quaternion_as_axis_and_angle(deg=True)

def center(coords):
  """ Returns the average of a list of vectors
  @param coords List of vectors to return the center of
  """
  for c in coords:
    if 'avg' not in locals():
      avg = c
    else:
      avg += c
  return avg / len(coords)

class basis(object):
  """ Bucket for detector element information """
  def __init__(self, orientation = None, translation = None, panelgroup = None):
    if orientation is None or translation is None:
      assert orientation is None and translation is None

      d_mat = panelgroup.get_local_d_matrix()
      fast = matrix.col((d_mat[0],d_mat[3],d_mat[6])).normalize()
      slow = matrix.col((d_mat[1],d_mat[4],d_mat[7])).normalize()
      orig = matrix.col((d_mat[2],d_mat[5],d_mat[8]))

      v3 = fast.cross(slow).normalize()

      r3 = matrix.sqr((fast[0],slow[0],v3[0],
                       fast[1],slow[1],v3[1],
                       fast[2],slow[2],v3[2]))

      self.orientation = r3.r3_rotation_matrix_as_unit_quaternion()
      self.translation = orig

    else:
      assert panelgroup is None
      self.orientation = orientation
      self.translation = translation

  def as_homogenous_transformation(self):
    """ Returns this basis change as a 4x4 transformation matrix in homogenous coordinates"""
    r3 = self.orientation.normalize().unit_quaternion_as_r3_rotation_matrix()
    return matrix.sqr((r3[0],r3[1],r3[2],self.translation[0],
                       r3[3],r3[4],r3[5],self.translation[1],
                       r3[6],r3[7],r3[8],self.translation[2],
                       0,0,0,1))

def read_optical_metrology_from_flat_file(path, detector, pixel_size, asic_dimension, asic_gap, plot = False, old_style_diff_path = None):
  """ Read a flat optical metrology file from LCLS and apply some corrections.  Partly adapted from xfel.metrology.flatfile
  @param path to the file to read
  @param detector Choice of CxiDs1 and XppDs1.  Affect how the quadrants are laid out (rotated manually or in absolute coordinates, respectively)
  @param pixel_size Size of each pixel in mm
  @param asic_dimension: the size of each sensor in pixels
  @param asic_gap: the pixel gap between the two asics on each pixel
  @param plot: if True, will plot the read and corrected metrology
  @param old_style_diff_path: if set to an old-style calibration directory, will print out and plot some comparisons between the metrology in this
  file and the metrology in the given directory
  @return dictionary of tuples in the form of (detector id), (detector id,quadrant id), (detector id,quadrant id,sensor id), and
  (detector id,quadrant id,sensor id,asic_id), mapping the levels of the hierarchy to basis objects
  """
  assert detector in ['CxiDs1', 'XppDs1']

  from xfel.metrology.flatfile import parse_metrology
  quadrants = parse_metrology(path, detector, plot, old_style_diff_path=old_style_diff_path)

  # rotate sensors 6 and 7 180 degrees
  for q_id, sensor in quadrants.iteritems():
    six = sensor[6]
    svn = sensor[7]
    assert len(six) == 4 and len(svn) == 4
    sensor[6] = [six[2],six[3],six[0],six[1]]
    sensor[7] = [svn[2],svn[3],svn[0],svn[1]]
    quadrants[q_id] = sensor

  quadrants_trans = {}
  if detector == "CxiDs1":
    # apply transformations: bring to order (slow, fast) <=> (column,
     # row).  This takes care of quadrant rotations
    for (q, sensors) in quadrants.iteritems():
      quadrants_trans[q] = {}

      q_apa = q

      if q == 0:
        # Q0:
        #   x -> -slow
        #   y -> -fast
        for (s, vertices) in sensors.iteritems():
          quadrants_trans[q][s] = [matrix.col((-v[1]/1000, +v[0]/1000, v[2]/1000))
                                   for v in quadrants[q_apa][s]]
      elif q == 1:
        # Q1:
        #   x -> +fast
        #   y -> -slow
        for (s, vertices) in sensors.iteritems():
          quadrants_trans[q][s] = [matrix.col((+v[0]/1000, +v[1]/1000, v[2]/1000))
                                   for v in quadrants[q_apa][s]]
      elif q == 2:
        # Q2:
        #   x -> +slow
        #   y -> +fast
        for (s, vertices) in sensors.iteritems():
          quadrants_trans[q][s] = [matrix.col((+v[1]/1000, -v[0]/1000, v[2]/1000))
                                   for v in quadrants[q_apa][s]]
      elif q == 3:
        # Q3:
        #   x -> -fast
        #   y -> +slow
        for (s, vertices) in sensors.iteritems():
          quadrants_trans[q][s] = [matrix.col((-v[0]/1000, -v[1]/1000, v[2]/1000))
                                   for v in quadrants[q_apa][s]]
      else:
        # NOTREACHED
        raise RuntimeError(
          "Detector does not have exactly four quadrants")
  else:
    assert detector == "XppDs1"
    # The XPP CSPAD has fixed quadrants.  For 2013-01-24 measurement,
    # they are all defined in a common coordinate system.
    #
    #   x -> +fast
    #   y -> -slow
    #
    # Fix the origin to the center of mass of sensor 1 in the four
    # quadrants.
    o = matrix.col((0, 0, 0))
    N = 0
    for (q, sensors) in quadrants.iteritems():
      o += center(sensors[1])
      N += 1
    o /= N

    for (q, sensors) in quadrants.iteritems():
      quadrants_trans[q] = {}
      for (s, vertices) in sensors.iteritems():
        quadrants_trans[q][s] = [matrix.col(((-v[1] - (-o[1]))/1000, (+v[0] - (+o[0]))/1000, v[2]/1000))
                                     for v in quadrants[q][s]]


  if plot:
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')

    a = []; b = []; c = []; d = []; cents = []

    for q_id, q in quadrants_trans.iteritems():
      q_c = matrix.col((0,0,0))
      for s_id, s in q.iteritems():
        q_c += center(s)
      q_c /= len(q)
      cents.append(q_c)
      for s_id, s in q.iteritems():
        sensor = ((s[0][0], s[0][1]),
                  (s[1][0], s[1][1]),
                  (s[2][0], s[2][1]),
                  (s[3][0], s[3][1]))
        cents.append(center([matrix.col(v) for v in sensor]))
        ax.add_patch(Polygon(sensor, closed=True, color='green', fill=False, hatch='/'))

        a.append(s[0]); b.append(s[1]); c.append(s[2]); d.append(s[3])

    ax.set_xlim((-100, 100))
    ax.set_ylim((-100, 100))
    plt.scatter([v[0] for v in a], [v[1] for v in a], c = 'black')
    for i, v in enumerate(a):
        ax.annotate(i*2, (v[0],v[1]))
    plt.scatter([v[0] for v in b], [v[1] for v in b], c = 'yellow')
    #plt.scatter([v[0] for v in c], [v[1] for v in c], c = 'yellow')
    plt.scatter([v[0] for v in d], [v[1] for v in d], c = 'yellow')
    plt.scatter([v[0] for v in cents], [v[1] for v in cents], c = 'red')
    plt.show()

  null_ori = matrix.col((0,0,1)).axis_and_angle_as_unit_quaternion(0, deg=True)
  metro = { (0,): basis(null_ori, matrix.col((0,0,0))) }

  for q_id, q in quadrants_trans.iteritems():
    # calculate the center of the quadrant
    q_c = matrix.col((0,0,0))
    for s_id, s in q.iteritems():
      q_c += center(s)
    q_c /= len(q)

    metro[(0,q_id)] = basis(null_ori,q_c)

    for s_id, s in q.iteritems():
      sensorcenter_wrt_detector = center(s)
      sensorcenter_wrt_quadrant = sensorcenter_wrt_detector - q_c

      # change of basis from a sensor in the plane of the detector, lying on its side, to oriented
      # as it should be, relative to the frame of the detector
      v1 = (s[1] - s[0]).normalize() # +x
      v2 = (s[0] - s[3]).normalize() # -y
      v3 = (v1.cross(v2)).normalize()
      rotation = matrix.sqr((v1[0],v2[0],v3[0],
                             v1[1],v2[1],v3[1],
                             v1[2],v2[2],v3[2]))

      metro[(0,q_id,s_id)] = basis(rotation.r3_rotation_matrix_as_unit_quaternion(),sensorcenter_wrt_quadrant)

      w = pixel_size * (asic_dimension[0]/2 + asic_gap/2)
      metro[(0,q_id,s_id,0)] = basis(null_ori,matrix.col((-w,0,0)))
      metro[(0,q_id,s_id,1)] = basis(null_ori,matrix.col((+w,0,0)))

  return metro

def cbf_file_to_basis_dict(path):
  """ Maps a cbf file to a dictionary of tuples and basis objects, in the same form as the above
  read_optical_metrology_from_flat_file
  @param path cbf file path """
  from dxtbx.format.Registry import Registry
  reader = Registry.find(path)
  instance = reader(path)
  root = instance.get_detector().hierarchy()

  d = 0 # only allow one detector for now
  metro = {(d,):basis(panelgroup=root)}
  metro[(d,)].translation = matrix.col((0,0,0))

  for q, quad in enumerate(root):
    metro[(d,q)] = basis(panelgroup=quad)
    for s, sensor in enumerate(quad):
      metro[(d,q,s)] = basis(panelgroup=sensor)
      for a, asic in enumerate(sensor):
        # at the asic level, need to subtract off the d0 vector so this asic's basis does not include
        # the shift from the asic center to the asic corner.  Also need to flip the Y back around
        # to be consistent with how it was originally stored
        d_mat = asic.get_local_d_matrix()
        fast = matrix.col((d_mat[0],d_mat[3],d_mat[6])).normalize()
        slow = matrix.col((d_mat[1],d_mat[4],d_mat[7])).normalize()
        orig = matrix.col((d_mat[2],d_mat[5],d_mat[8]))

        v3 = fast.cross(slow).normalize()

        r3 = matrix.sqr((fast[0],slow[0],v3[0],
                         fast[1],slow[1],v3[1],
                         fast[2],slow[2],v3[2]))

        transform = matrix.sqr((1, 0, 0,
                                0,-1, 0,
                                0, 0,-1))

        pix_size = asic.get_pixel_size()
        img_size = asic.get_image_size()

        offset = matrix.col((-pix_size[1]*(img_size[1]-1)/2,
                             +pix_size[0]*(img_size[0]-1)/2,0))

        metro[(d,q,s,a)] = basis(orientation=(r3*transform).r3_rotation_matrix_as_unit_quaternion(),
                                 translation=orig-offset)

  return metro

def metro_phil_to_basis_dict(metro):
  """ Maps a phil object from xfel.cftbx.detector.metrology2phil to a dictionary of tuples and basis objects, in the same form
  as the above read_optical_metrology_from_flat_file
  @param metro unextracted phil scope object that contains translations and rotations for each asic """
  for o in metro.objects:
    if o.is_scope:
      #one of the subkeys of the root object will be the detector phil. it will be the only one not extracted.
      detector_phil = o.extract()
      break
  #metro = metro.extract() # not needed

  bd = {(detector_phil.serial,): basis(matrix.col(detector_phil.orientation),
                                       matrix.col(detector_phil.translation)*1000) }
  for p in detector_phil.panel:
    bd[(detector_phil.serial,p.serial)] = basis(matrix.col(p.orientation),
                                                matrix.col(p.translation)*1000)
    for s in p.sensor:
      bd[(detector_phil.serial,p.serial,s.serial)] = basis(matrix.col(s.orientation),
                                                           matrix.col(s.translation)*1000)
      for a in s.asic:
        bd[(detector_phil.serial,p.serial,s.serial,a.serial)] = basis(matrix.col(a.orientation),
                                                                      matrix.col(a.translation)*1000)

  return bd

def write_cspad_cbf(tiles, metro, metro_style, timestamp, destpath, wavelength, distance, verbose = True, header_only = False):
  assert metro_style in ['calibdir','flatfile','cbf']
  if metro_style == 'calibdir':
    metro = metro_phil_to_basis_dict(metro)

  # set up the metrology dictionary to include axis names, pixel sizes, and so forth
  try:
    from xfel.cxi.cspad_ana.cspad_tbx import dynamic_range as dr
    from xfel.cxi.cspad_ana.cspad_tbx import pixel_size as ps
    dynamic_range = dr
    pixel_size = ps
  except ImportError:
    dynamic_range = 2**14 - 1
    pixel_size = 110e-3

  dserial = None
  dname = None
  detector_axes_names = [] # save these for later
  for key in sorted(metro):
    basis = metro[key]
    if len(key) == 1:
      assert dserial is None # only one detector allowed for now
      dserial = key[0]

      dname = "AXIS_D%d"%dserial #XXX check if DS1 is here
      for a in ["_X","_Y","_Z","_R"]: detector_axes_names.append(dname+a)
      basis.equipment_component = "detector_arm"
      basis.depends_on = dname+"_X"
    elif len(key) == 2:
      detector_axes_names.append("FS_D%dQ%d"%key)
      basis.equipment_component = "detector_quadrant"
      basis.depends_on = dname+"_R"
    elif len(key) == 3:
      detector_axes_names.append("FS_D%dQ%dS%d"%(key))
      basis.equipment_component = "detector_sensor"
      basis.depends_on = "FS_D%dQ%d"%key[0:2]
    elif len(key) == 4:
      detector_axes_names.append("FS_D%dQ%dS%dA%d"%(key))
      basis.equipment_component = "detector_asic"
      basis.depends_on = "FS_D%dQ%dS%d"%key[0:3]
      basis.pixel_size = (pixel_size,pixel_size)
      if tiles is None:
        basis.dimension = asic_dimension
      else:
        basis.dimension = tuple(reversed(tiles[key].focus()))
      basis.saturation = dynamic_range
    else:
      assert False # shouldn't be reached as it would indicate more than four levels of hierarchy for this detector
    basis.axis_name = detector_axes_names[-1]


  # the data block is the root cbf node
  cbf=cbf_wrapper()
  cbf.new_datablock(os.path.splitext(os.path.basename(destpath))[0])

  # Each category listed here is preceded by the imageCIF description taken from here:
  # http://www.iucr.org/__data/iucr/cifdic_html/2/cif_img.dic/index.html

  """Data items in the DIFFRN category record details about the
   diffraction data and their measurement."""
  cbf.add_category("diffrn",["id"])
  cbf.add_row(["DS1"])

  """Data items in the DIFFRN_SOURCE category record details of
  the source of radiation used in the diffraction experiment."""
  cbf.add_category("diffrn_source", ["diffrn_id","source","type"])
  cbf.add_row(["DS1","xfel","LCLS CXI Endstation"])

  """Data items in the DIFFRN_RADIATION category describe
   the radiation used for measuring diffraction intensities,
   its collimation and monochromatization before the sample.

   Post-sample treatment of the beam is described by data
   items in the DIFFRN_DETECTOR category."""
  if not header_only:
    cbf.add_category("diffrn_radiation", ["diffrn_id","wavelength_id","probe"])
    cbf.add_row(["DS1","WAVELENGTH1","x-ray"])

  """ Data items in the DIFFRN_RADIATION_WAVELENGTH category describe
   the wavelength of the radiation used in measuring the diffraction
   intensities. Items may be looped to identify and assign weights
   to distinct wavelength components from a polychromatic beam."""
  if not header_only:
    cbf.add_category("diffrn_radiation_wavelength", ["id","wavelength","wt"])
    cbf.add_row(["WAVELENGTH1","%f"%wavelength,"1.0"])

  """Data items in the DIFFRN_DETECTOR category describe the
   detector used to measure the scattered radiation, including
   any analyser and post-sample collimation."""
  cbf.add_category("diffrn_detector", ["diffrn_id","id","type","details","number_of_axes"])
  cbf.add_row(["DS1","CSPAD_FRONT","CS PAD",".",str(len(detector_axes_names))])

  """Data items in the DIFFRN_DETECTOR_AXIS category associate
     axes with detectors."""
  # Note, does not include the fast and the slow axes
  cbf.add_category("diffrn_detector_axis",["detector_id","axis_id"])
  for name in detector_axes_names:
    cbf.add_row(["CSPAD_FRONT",name])

  tilestrs = []
  tilekeys = sorted([key for key in metro if len(key) == 4])

  # create a series of strings representing each asic.  Q here is for quadrant instead of panel.
  for tilekey in tilekeys:
    tilestrs.append("D%dQ%dS%dA%d"%tilekey)

  """Data items in the DIFFRN_DETECTOR_ELEMENT category record
   the details about spatial layout and other characteristics
   of each element of a detector which may have multiple elements."""
  cbf.add_category("diffrn_detector_element",["id","detector_id"])

  for tilename in tilestrs:
    cbf.add_row(["ELE_" + tilename, "CSPAD_FRONT"])

  """Data items in the DIFFRN_DATA_FRAME category record
   the details about each frame of data."""
  cbf.add_category("diffrn_data_frame",["id","detector_element_id","array_id","binary_id"])

  for i, tilename in enumerate(tilestrs):
    cbf.add_row(["FRAME1","ELE_"+tilename,"ARRAY_"+tilename,"%d"%(i+1)])

  """Data items in the DIFFRN_MEASUREMENT category record details
   about the device used to orient and/or position the crystal
   during data measurement and the manner in which the
   diffraction data were measured."""
  if not header_only:
    cbf.add_category("diffrn_measurement",["diffrn_id","id","number_of_axes","method","details"])
    cbf.add_row(["DS1","INJECTION","0","electrospray","crystals injected by electrospray"])

  """ Data items in the DIFFRN_SCAN category describe the parameters of one
     or more scans, relating axis positions to frames."""
  if not header_only:
    cbf.add_category("diffrn_scan",["id","frame_id_start","frame_id_end","frames"])
    cbf.add_row(["SCAN1","FRAME1","FRAME1","1"])

  """Data items in the DIFFRN_SCAN_FRAME category describe
   the relationships of particular frames to scans."""
  if not header_only:
    cbf.add_category("diffrn_scan_frame",["frame_id","frame_number","integration_time","scan_id","date"])
    cbf.add_row(["FRAME1","1","0.0","SCAN1",timestamp])

  """Data items in the AXIS category record the information required
     to describe the various goniometer, detector, source and other
     axes needed to specify a data collection.  The location of each
     axis is specified by two vectors: the axis itself, given as a unit
     vector, and an offset to the base of the unit vector.  These vectors
     are referenced to a right-handed laboratory coordinate system with
     its origin in the sample or specimen"""
  # More detail here: http://www.iucr.org/__data/iucr/cifdic_html/2/cif_img.dic/Caxis.html
  # Note we also use two new columns not in the latest imageCIF dictionary: rotation and rotation_axis.
  # We use them to specify an translation and a rotation in a single axis to describe a change in setting
  # when moving from one frame (say a quadrant) to another (say a sensor)
  cbf.add_category("axis",["id","type","equipment","depends_on","vector[1]","vector[2]","vector[3]",
                                                                "offset[1]","offset[2]","offset[3]",
                                                                "equipment_component"])

  # Keep a list of rows to add to the scan frame axis table
  axis_settings = []
  # keep a list of rows to add to the scan axis table
  axis_names = []

  # Create a series of axis describing frame shifts from each level of the detector to the next
  cbf.add_row( "AXIS_SOURCE  general     source   .        0  0  1 . . . .".split())                           ; axis_names.append("AXIS_SOURCE")
  cbf.add_row( "AXIS_GRAVITY general     gravity  .        0 -1  0 . . . .".split())                           ; axis_names.append("AXIS_GRAVITY")
  cbf.add_row(("%s_Z         translation detector .        0  0  1 . . . detector_arm"%(dname)).split())       ; axis_names.append("%s_Z"%dname)
  cbf.add_row(("%s_Y         translation detector %s_Z     0  1  0 . . . detector_arm"%(dname,dname)).split()) ; axis_names.append("%s_Y"%dname)
  cbf.add_row(("%s_X         translation detector %s_Y     1  0  0 . . . detector_arm"%(dname,dname)).split()) ; axis_names.append("%s_X"%dname)

  axis_settings.append(["AXIS_SOURCE" ,"FRAME1","0","0"])
  axis_settings.append(["AXIS_GRAVITY","FRAME1","0","0"])
  axis_settings.append([dname+"_X"    ,"FRAME1","0","0"])
  axis_settings.append([dname+"_Y"    ,"FRAME1","0","0"])
  axis_settings.append([dname+"_Z"    ,"FRAME1","0",str(-distance)])

  for key in sorted(metro):
    basis = metro[key]
    assert len(key) > 0 and len(key) <= 4

    cbf.add_frame_shift(basis, axis_settings)
    axis_names.append(basis.axis_name)

    if len(key) == 4:

      dim_pixel = basis.pixel_size
      dim_readout = basis.dimension

      # Add the two vectors for each asic that describe the fast and slow directions pixels should be laid out in real space
      offset_fast = -dim_pixel[0]*((dim_readout[1] - 1) / 2)
      offset_slow = +dim_pixel[1]*((dim_readout[0] - 1) / 2)

      aname = "D%dQ%dS%dA%d"%key

      cbf.add_row(["AXIS_"+ aname + "_S", "translation","detector",basis.axis_name    ,"0", "-1","0","%f"%offset_fast,"%f"%offset_slow,"0.0", "detector_asic"])
      cbf.add_row(["AXIS_"+ aname + "_F", "translation","detector","AXIS_"+aname +"_S","1","0","0","0","0","0.0", "detector_asic"])
      axis_names.append("AXIS_"+ aname + "_F"); axis_names.append("AXIS_"+ aname + "_S")
      axis_settings.append(["AXIS_"+ aname + "_F","FRAME1","0","0"])
      axis_settings.append(["AXIS_"+ aname + "_S","FRAME1","0","0"])

  """Data items in the DIFFRN_SCAN_AXIS category describe the settings of
     axes for particular scans.  Unspecified axes are assumed to be at
     their zero points."""
  # leave all the settings zero. Levels with settings are set below.
  cbf.add_category("diffrn_scan_axis",["axis_id","scan_id","angle_start","angle_range","angle_increment",
                                       "displacement_start","displacement_range","displacement_increment"])
  for name in axis_names:
    cbf.add_row([name,"SCAN1","0","0","0","0","0","0"])

  """Data items in the DIFFRN_SCAN_FRAME_AXIS category describe the
     settings of axes for particular frames.  Unspecified axes are
     assumed to be at their zero points."""
  cbf.add_category("diffrn_scan_frame_axis",["axis_id","frame_id","angle","displacement"])
  for row in axis_settings:
    cbf.add_row(row)

  """Data items in the ARRAY_STRUCTURE_LIST category record the size
     and organization of each array dimension.
     The relationship to physical axes may be given."""
  cbf.add_category("array_structure_list",["array_id","index","dimension","precedence","direction","axis_set_id"])
  for tilename,tilekey in zip(tilestrs,tilekeys):
    b = metro[tilekey]
    cbf.add_row(["ARRAY_"+tilename,"1","%d"%b.dimension[0],"1","increasing","AXIS_"+tilename+"_F"])
    cbf.add_row(["ARRAY_"+tilename,"2","%d"%b.dimension[1],"2","increasing","AXIS_"+tilename+"_S"])

  """Data items in the ARRAY_STRUCTURE_LIST_AXIS category describe
     the physical settings of sets of axes for the centres of pixels that
     correspond to data points described in the
     ARRAY_STRUCTURE_LIST category."""
  cbf.add_category("array_structure_list_axis",["axis_set_id","axis_id","displacement","displacement_increment"])
  for tilename,tilekey in zip(tilestrs,tilekeys):
    cbf.add_row(["AXIS_"+tilename+"_F","AXIS_"+tilename+"_F","0.0","%f"%(metro[tilekey].pixel_size[0])])
    cbf.add_row(["AXIS_"+tilename+"_S","AXIS_"+tilename+"_S","0.0","%f"%(metro[tilekey].pixel_size[1])])

  # rest of these involve the binary data
  if not header_only:

    """ Data items in the ARRAY_INTENSITIES category record the
     information required to recover the intensity data from
     the set of data values stored in the ARRAY_DATA category."""
    # More detail here: http://www.iucr.org/__data/iucr/cifdic_html/2/cif_img.dic/Carray_intensities.html
    cbf.add_category("array_intensities",["array_id","binary_id","linearity","gain","gain_esd","overload","undefined_value"])
    for i, (tilename, tilekey) in enumerate(zip(tilestrs,tilekeys)):
      cbf.add_row(["ARRAY_"+ tilename,str(i+1),"linear","1.0","0.1",str(metro[tilekey].saturation),"0.0"])

    """ Data items in the ARRAY_STRUCTURE category record the organization and
       encoding of array data in the ARRAY_DATA category."""
    cbf.add_category("array_structure",["id","encoding_type","compression_type","byte_order"])

    for tilename in tilestrs:
      cbf.add_row(["ARRAY_"+tilename,"signed 64-bit real IEEE","canonical","little_endian"])

    """ Data items in the ARRAY_STRUCTURE category record the organization and
       encoding of array data in the ARRAY_DATA category."""
    cbf.add_category("array_data",["array_id","binary_id","data"])

    if verbose:
      print "Compressing tiles...",

    for i, tilekey in enumerate(tilekeys):
      detector, quadrant, sensor, asic = tilekey
      focus = tiles[tilekey].focus()

      cbf.add_row(["ARRAY_" + tilestrs[i],str(i+1)])

      binary_id = i+1
      data = tiles[tilekey].as_double().copy_to_byte_str()
      elsize = 8
      elements = len(tiles[tilekey])
      byteorder = "little_endian"
      dimfast = focus[1]
      dimmid = focus[0]
      dimslow = 1
      padding = 0

      cbf.set_realarray_wdims_fs(\
        pycbf.CBF_CANONICAL,
        binary_id,
        data,
        elsize,
        elements,
        byteorder,
        dimfast,
        dimmid,
        dimslow,
        padding)

  cbf.write_widefile(destpath,pycbf.CBF,\
      pycbf.MIME_HEADERS|pycbf.MSG_DIGEST|pycbf.PAD_4K,0)

  if verbose:
    print "%s written"%destpath
