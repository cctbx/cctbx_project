from __future__ import absolute_import, division, print_function
from six.moves import range
#-*- Mode: Python; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
#
# $Id: cspad_cbf_tbx.py
#

import pycbf, os
from scitbx import matrix
from scitbx.array_family import flex
import six
from six.moves import zip
from six.moves import map

# need to define these here since it not defined in SLAC's metrology definitions
asic_dimension = (194,185)
asic_gap = 3
pixel_size = 0.10992
from xfel.cxi.cspad_ana.cspad_tbx import cspad_saturated_value, cspad_min_trusted_value

def get_psana_corrected_data(psana_det, evt, use_default=False, dark=True, common_mode=None, apply_gain_mask=True,
                             gain_mask_value=None, per_pixel_gain=False, gain_mask=None):
  """
  Given a psana Detector object, apply corrections as appropriate and return the data from the event
  @param psana_det psana Detector object
  @param evt psana event
  @param use_default If true, apply the default calibration only, using the psana algorithms. Otherise, use the corrections
  specified by the rest of the flags and values passed in.
  @param dark Whether to apply the detector dark, bool or numpy array
  @param common_mode Which common mode algorithm to apply. None: apply no algorithm. Default: use the algorithm specified
  in the calib folder. Otherwise should be a list as specified by the psana documentation for common mode customization
  @param apply_gain_mask Whether to apply the common mode gain mask correction
  @param gain_mask_value Multiplier to apply to the pixels, according to the gain mask
  @param per_pixel_gain If available, use the per pixel gain deployed to the calibration folder
  @param, gain_mask, gain mask showing which pixels to apply gain mask value
  @return Numpy array corrected as specified.
  """
  # order is pedestals, then common mode, then gain mask, then per pixel gain
  import numpy as np
  run = evt.run()

  if use_default:
    return psana_det.calib(evt)  # applies psana's complex run-dependent calibrations
  data = psana_det.raw_data(evt)
  if data is None:
    return

  data = data.astype(np.float64)
  if isinstance(dark, bool):
    if dark:
      data -= psana_det.pedestals(run)
  elif isinstance( dark, np.ndarray ):
    data -= dark

  if common_mode is not None and common_mode != "default":
    if common_mode == 'cspad_default':
      common_mode = (1,25,25,100,1)  # default parameters for CSPAD images
      psana_det.common_mode_apply(run, data, common_mode)
    elif common_mode == 'unbonded':
      common_mode = (5,0,0,0,0)  # unbonded pixels used for correction
      psana_det.common_mode_apply(run, data, common_mode)
    else:  # this is how it was before.. Though I think common_mode would need to be a tuple..
      psana_det.common_mode_apply(run, data, common_mode)
  if apply_gain_mask:
    if gain_mask is None:  # TODO: consider try/except here
      gain_mask = psana_det.gain_mask(run) == 1
    if gain_mask_value is None:
      try:
        gain_mask_value = psana_det._gain_mask_factor
      except AttributeError:
        print("No gain set for psana detector, using gain value of 1, consider disabling gain in your phil file")
        gain_mask_value = 1
    data[gain_mask] = data[gain_mask]*gain_mask_value
  if per_pixel_gain: # TODO: test this
    data *= psana_det.gain(run)
  return data

from dxtbx.format.FormatCBFMultiTile import cbf_wrapper as dxtbx_cbf_wrapper
class cbf_wrapper(dxtbx_cbf_wrapper):
  """ Wrapper class that provides convenience functions for working with cbflib"""

  def add_frame_shift(self, basis, axis_settings):
    """Add an axis representing a frame shift (a rotation axis with an offset)"""
    angle, axis = angle_and_axis(basis)

    if angle == 0:
      axis = (0,0,1)

    if basis.include_translation:
      translation = basis.translation
    else:
      translation = (0,0,0)

    self.add_row([basis.axis_name,"rotation","detector",basis.depends_on,
                  str(axis[0]),str(axis[1]),str(axis[2]),
                  str(translation[0]),
                  str(translation[1]),
                  str(translation[2]),
                  basis.equipment_component])

    axis_settings.append([basis.axis_name, "FRAME1", str(angle), "0"])

def angle_and_axis(basis):
  """Normalize a quaternion and return the angle and axis
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
  def __init__(self, orientation = None, translation = None, panelgroup = None, homogenous_transformation = None):
    """
    Provide only orientation + translation or a panelgroup or a homogenous_transformation.

    @param orientation rotation in the form of a quarternion
    @param translation vector translation in relation to the parent frame
    @param panelgroup dxtbx panelgroup object whose local d matrix will represent the
    basis shift
    @param homogenous_transformation 4x4 matrix.sqr object representing a translation
    and a rotation. Must not also contain a scale as this won't be decomposed properly.
    """
    self.include_translation = True

    if panelgroup is not None:
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

    elif orientation is not None or translation is not None:
      assert orientation is not None and translation is not None
      self.orientation = orientation
      self.translation = translation

    else:
      # Decompose the homegenous transformation assuming no scale factors were used
      h = homogenous_transformation
      self.orientation = matrix.sqr((h[0],h[1],h[2],
                                     h[4],h[5],h[6],
                                     h[8],h[9],h[10])).r3_rotation_matrix_as_unit_quaternion()
      self.translation = matrix.col((h[3],
                                     h[7],
                                     h[11]))
      assert h[12] == h[13] == h[14] == 0 and h[15] == 1

  def as_homogenous_transformation(self):
    """ Returns this basis change as a 4x4 transformation matrix in homogenous coordinates"""
    r3 = self.orientation.normalize().unit_quaternion_as_r3_rotation_matrix()
    return matrix.sqr((r3[0],r3[1],r3[2],self.translation[0],
                       r3[3],r3[4],r3[5],self.translation[1],
                       r3[6],r3[7],r3[8],self.translation[2],
                       0,0,0,1))

  def __mul__(self, other):
    """ Use homogenous matrices to multiply bases together """
    if hasattr(other, 'as_homogenous_transformation'):
      return basis(homogenous_transformation = self.as_homogenous_transformation() * other.as_homogenous_transformation())
    elif hasattr(other, 'n'):
      if other.n == (3,1):
        b = matrix.col((other[0], other[1], other[2], 1))
      elif other.n == (4,1):
        b = other
      else:
        raise TypeError(b, "Incompatible matrices")
      p = self.as_homogenous_transformation() * b
      if other.n == (3,1):
        return matrix.col(p[0:3])
      else:
        return p
    else:
      raise TypeError(b)


def basis_from_geo(geo, use_z = True):
  """ Given a psana GeometryObject, construct a basis object """
  rotx = matrix.col((1,0,0)).axis_and_angle_as_r3_rotation_matrix(
    geo.rot_x + geo.tilt_x, deg=True)
  roty = matrix.col((0,1,0)).axis_and_angle_as_r3_rotation_matrix(
    geo.rot_y + geo.tilt_y, deg=True)
  rotz = matrix.col((0,0,1)).axis_and_angle_as_r3_rotation_matrix(
    geo.rot_z + geo.tilt_z, deg=True)

  rot = (rotx*roty*rotz).r3_rotation_matrix_as_unit_quaternion()

  if use_z:
    trans = matrix.col((geo.x0/1000, geo.y0/1000, geo.z0/1000))
  else:
    trans = matrix.col((geo.x0/1000, geo.y0/1000, 0))

  return basis(orientation = rot, translation = trans)

def read_slac_metrology(path = None, geometry = None, plot=False, include_asic_offset=False):
  if path is None and geometry is None:
    raise Sorry("Need to provide a geometry object or a path to a geometry file")

  if path is not None and geometry is not None:
    raise Sorry("Cannot provide a geometry object and a geometry file. Ambiguous")

  if geometry is None:
    try:
      from PSCalib.GeometryAccess import GeometryAccess
      geometry = GeometryAccess(path)
    except Exception as e:
      raise Sorry("Can't parse this metrology file")

  metro = {}
  pixel_size = geometry.get_pixel_scale_size()/1000
  null_ori = matrix.col((0,0,1)).axis_and_angle_as_unit_quaternion(0, deg=True)

  # collapse any transformations above those of the quadrants into one X/Y offset,
  # but don't keep Z transformations, as those come from the XTC stream
  root = geometry.get_top_geo()
  root_basis = basis_from_geo(root, use_z=False)
  while len(root.get_list_of_children()) != 4 and len(root.get_list_of_children()) != 32:
    assert len(root.get_list_of_children()) == 1
    root = root.get_list_of_children()[0]
    root_basis *= basis_from_geo(root, use_z=False)

  metro[(0,)] = root_basis

  def add_sensor(quad_id, sensor_id, sensor):
    metro[(0,quad_id,sensor_id)] = basis_from_geo(sensor)

    x, y, z = sensor.get_pixel_coords()
    x/=1000; y/=1000; z/=1000
    assert x.shape == y.shape == z.shape
    sensor_px_slow = x.shape[0]
    sensor_px_fast = x.shape[1]
    assert sensor_px_fast % 2 == 0

    a0ul = sul = matrix.col((x[0,0],y[0,0],z[0,0]))
    a1ur = sur = matrix.col((x[0,sensor_px_fast-1],y[0,sensor_px_fast-1],z[0,sensor_px_fast-1]))
    a1lr = slr = matrix.col((x[sensor_px_slow-1,sensor_px_fast-1],y[sensor_px_slow-1,sensor_px_fast-1],z[sensor_px_slow-1,sensor_px_fast-1]))
    a0ll = sll = matrix.col((x[sensor_px_slow-1,0],y[sensor_px_slow-1,0],z[sensor_px_slow-1,0]))

    a0ur = matrix.col((x[0,sensor_px_fast//2-1],y[0,sensor_px_fast//2-1],z[0,sensor_px_fast//2-1]))
    a0lr = matrix.col((x[sensor_px_slow-1,sensor_px_fast//2-1],y[sensor_px_slow-1,sensor_px_fast//2-1],z[sensor_px_slow-1,sensor_px_fast//2-1]))

    a1ul = matrix.col((x[0,sensor_px_fast//2],y[0,sensor_px_fast//2],z[0,sensor_px_fast//2]))
    a1ll = matrix.col((x[sensor_px_slow-1,sensor_px_fast//2],y[sensor_px_slow-1,sensor_px_fast//2],z[sensor_px_slow-1,sensor_px_fast//2]))

    sensor_center = center([sul,sur,slr,sll])
    asic0_center = center([a0ul,a0ur,a0lr,a0ll])
    asic1_center = center([a1ul,a1ur,a1lr,a1ll])

    asic_trans0 = (asic0_center-sensor_center).length()
    asic_trans1 = (asic1_center-sensor_center).length()

    if include_asic_offset:
      rotated_ori = matrix.col((1,0,0)).axis_and_angle_as_unit_quaternion(180.0, deg=True)
      offset_fast = -pixel_size*((sensor_px_fast) / 4) # 4 because sensor_px_fast is for sensor
      offset_slow = +pixel_size*((sensor_px_slow) / 2) # Sensor is divided into 2 only in fast direction
      metro[(0,quad_id,sensor_id,0)] = basis(orientation=rotated_ori,translation=matrix.col((-asic_trans0,0,0)))
      metro[(0,quad_id,sensor_id,1)] = basis(orientation=rotated_ori,translation=matrix.col((+asic_trans1,0,0)))
      metro[(0,quad_id,sensor_id,0)].translation += matrix.col((offset_fast, offset_slow, 0))
      metro[(0,quad_id,sensor_id,1)].translation += matrix.col((offset_fast, offset_slow, 0))
    else:
      metro[(0,quad_id,sensor_id,0)] = basis(orientation=null_ori,translation=matrix.col((-asic_trans0,0,0)))
      metro[(0,quad_id,sensor_id,1)] = basis(orientation=null_ori,translation=matrix.col((+asic_trans1,0,0)))

  if len(root.get_list_of_children()) == 4:
    for quad_id, quad in enumerate(root.get_list_of_children()):
      metro[(0,quad_id)] = basis_from_geo(quad)
      for sensor_id, sensor in enumerate(quad.get_list_of_children()):
        add_sensor(quad_id, sensor_id, sensor)
  elif len(root.get_list_of_children()) == 32:
    for quad_id in range(4):
      metro[(0,quad_id)] = basis(orientation = null_ori, translation = matrix.col((0,0,0)))
      sensors = root.get_list_of_children()
      for sensor_id in range(8):
        add_sensor(quad_id, sensor_id, sensors[quad_id*4+sensor_id])
  else:
    assert False

  return metro

def get_calib_file_path(env, address, run):
  """ Findes the path to the SLAC metrology file stored in a psana environment
      object's calibration store
      @param env psana environment object
      @param address address string for a detector
      @param run psana run object or run number
  """
  from psana import Detector
  try:
    # try to get it from the detector interface
    psana_det = Detector(address, run.env())
    return psana_det.pyda.geoaccess(run.run()).path
  except Exception as e:
    pass

  # try to get it from the calib store directly
  from psana import ndarray_uint8_1, Source
  cls = env.calibStore()
  src = Source('DetInfo(%s)'%address)
  path_nda = cls.get(ndarray_uint8_1, src, 'geometry-calib')
  if path_nda is None:
    return None
  return ''.join(map(chr, path_nda))

def env_dxtbx_from_slac_metrology(run, address):
  """ Loads a dxtbx cspad cbf header only object from the metrology path stored
      in a psana run object's calibration store
      @param env psana run object
      @param address address string for a detector
  """
  from xfel.command_line.xtc_process import PSANA2_VERSION
  if PSANA2_VERSION:
    det = run.ds.Detector(address)
    geometry = det.geometry(run)
  else:
    from psana import Detector
    try:
      # try to load the geometry from the detector interface
      psana_det = Detector(address, run.env())
      geometry = psana_det.pyda.geoaccess(run.run())
    except Exception as e:
      geometry = None

  if geometry is None:
    metro_path = get_calib_file_path(run.env(), address, run)
  elif geometry.valid:
    metro_path = None
  else:
    from libtbx.utils import Sorry
    import socket, os
    raise Sorry("Could not read geometry, hostname: %s"%socket.gethostname())

  if metro_path is None and geometry is None:
    return None

  metro = read_slac_metrology(metro_path, geometry)

  cbf = get_cspad_cbf_handle(None, metro, 'cbf', None, "test", None, 100, verbose = True, header_only = True)

  from dxtbx.format.FormatCBFCspad import FormatCBFCspadInMemory
  return FormatCBFCspadInMemory(cbf)

def format_object_from_data(base_dxtbx, data, distance, wavelength, timestamp, address, round_to_int=True):
  """
  Given a preloaded dxtbx format object and raw data, assemble the tiles
  and set the distance.
  @param base_dxtbx A header only dxtbx format object
  @param data 32x185x388 CSPAD byte array from XTC stream
  @param distance Detector distance (mm)
  @param wavelength Shot wavelength (angstroms)
  @param timestamp Human readable timestamp
  @param address Detector address, put in CBF header
  """
  from dxtbx.format.FormatCBFCspad import FormatCBFCspadInMemory
  import copy
  import numpy as np
  cbf = copy_cbf_header(base_dxtbx._cbf_handle)
  cspad_img = FormatCBFCspadInMemory(cbf)
  cbf.set_datablockname(address + "_" + timestamp)

  if round_to_int:
    data = flex.double(data.astype(np.float64)).iround()
  else:
    data = flex.double(data.astype(np.float64))
  data.reshape(flex.grid((4,8,185,388)))

  n_asics = data.focus()[0] * data.focus()[1]
  add_frame_specific_cbf_tables(cbf, wavelength,timestamp,
    [(cspad_min_trusted_value,cspad_saturated_value)]*n_asics)

  # Set the distance, I.E., the length translated along the Z axis
  cbf.find_category("diffrn_scan_frame_axis")
  cbf.find_column("axis_id")
  cbf.find_row("AXIS_D0_Z") # XXX discover the Z axis somehow, don't use D0 here
  cbf.find_column("displacement")
  cbf.set_value(str(-distance))

  # Explicitly reset the detector object now that the distance is set correctly
  cspad_img._detector_instance = cspad_img._detector()

  # Explicitly set up the beam object now that the tables are all loaded correctly
  cspad_img._beam_instance = cspad_img._beam()

  # Get the data and add it to the cbf handle. Split it out by quads.
  tiles = {}
  for i in range(4):
    tiles[(0,i)] = data[i:i+1,:,:,:]
    tiles[(0,i)].reshape(flex.grid((8,185,388)))
  add_tiles_to_cbf(cbf,tiles)

  return cspad_img


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

  #if plot:
    #print "Showing un-adjusted, parsed metrology from flatfile"
    #import matplotlib.pyplot as plt
    #from matplotlib.patches import Polygon
    #fig = plt.figure()
    #ax = fig.add_subplot(111, aspect='equal')
    #cents = []

    #for q_id, quadrant in quadrants.iteritems():
      #for s_id, sensor in quadrant.iteritems():
        #ax.add_patch(Polygon([p[0:2] for p in sensor], closed=True, color='green', fill=False, hatch='/'))
        #cents.append(center(sensor)[0:2])

    #for i, c in enumerate(cents):
      #ax.annotate(i*2,c)
    #ax.set_xlim((0, 200000))
    #ax.set_ylim((0, 200000))
    #plt.show()

  quadrants_trans = {}
  if detector == "CxiDs1":
    # rotate sensors 6 and 7 180 degrees
    for q_id, quadrant in six.iteritems(quadrants):
      six = quadrant[6]
      svn = quadrant[7]
      assert len(six) == 4 and len(svn) == 4
      quadrant[6] = [six[2],six[3],six[0],six[1]]
      quadrant[7] = [svn[2],svn[3],svn[0],svn[1]]
      quadrants[q_id] = quadrant

    # apply transformations: bring to order (slow, fast) <=> (column,
     # row).  This takes care of quadrant rotations
    for (q, sensors) in six.iteritems(quadrants):
      quadrants_trans[q] = {}

      q_apa = q

      if q == 0:
        # Q0:
        #   x -> -slow
        #   y -> -fast
        for (s, vertices) in six.iteritems(sensors):
          quadrants_trans[q][s] = [matrix.col((-v[1]/1000, +v[0]/1000, v[2]/1000))
                                   for v in quadrants[q_apa][s]]
      elif q == 1:
        # Q1:
        #   x -> +fast
        #   y -> -slow
        for (s, vertices) in six.iteritems(sensors):
          quadrants_trans[q][s] = [matrix.col((+v[0]/1000, +v[1]/1000, v[2]/1000))
                                   for v in quadrants[q_apa][s]]
      elif q == 2:
        # Q2:
        #   x -> +slow
        #   y -> +fast
        for (s, vertices) in six.iteritems(sensors):
          quadrants_trans[q][s] = [matrix.col((+v[1]/1000, -v[0]/1000, v[2]/1000))
                                   for v in quadrants[q_apa][s]]
      elif q == 3:
        # Q3:
        #   x -> -fast
        #   y -> +slow
        for (s, vertices) in six.iteritems(sensors):
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
    slen = len(quadrants[list(quadrants.keys())[0]])
    for (q, sensors) in six.iteritems(quadrants):
      assert len(sensors) == slen
      for s, sensor in six.iteritems(sensors):
        o += center(sensor)
        N += 1
    o /= N

    rot_mat = matrix.col((0,0,1)).axis_and_angle_as_r3_rotation_matrix(180,deg=True)
    sensors_to_rotate = [4,5,10,11,12,13,14,15,16,17,18,19,22,23,24,25]

    for (q, quadrant) in six.iteritems(quadrants):
      quadrants_trans[q] = {}
      for (s, vertices) in six.iteritems(quadrant):
        # move to origin, rotate 180 degrees around origin, and scale
        vertices = [rot_mat*(v-o)/1000 for v in vertices]

        #rotate a subset of the sensors into position (determined emperically)
        if (q*8)+s in sensors_to_rotate:
          vertices = [vertices[2],vertices[3],vertices[0],vertices[1]]

        quadrants_trans[q][s] = vertices


  if plot:
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')

    a = []; b = []; c = []; d = []; cents = {}

    for q_id, q in six.iteritems(quadrants_trans):
      q_c = matrix.col((0,0,0))
      for s_id, s in six.iteritems(q):
        q_c += center(s)
      q_c /= len(q)
      cents["Q%d"%q_id] = q_c
      for s_id, s in six.iteritems(q):
        sensor = ((s[0][0], s[0][1]),
                  (s[1][0], s[1][1]),
                  (s[2][0], s[2][1]),
                  (s[3][0], s[3][1]))
        cents["S%d"%((q_id*8)+s_id)] = center([matrix.col(v) for v in sensor])
        ax.add_patch(Polygon(sensor, closed=True, color='green', fill=False, hatch='/'))

        a.append(s[0]); b.append(s[1]); c.append(s[2]); d.append(s[3])

    ax.set_xlim((-100, 100))
    ax.set_ylim((-100, 100))
    plt.scatter([v[0] for v in a], [v[1] for v in a], c = 'black')
    for i, v in six.iteritems(cents):
        ax.annotate(i, (v[0],v[1]))
    plt.scatter([v[0] for v in b], [v[1] for v in b], c = 'yellow')
    #plt.scatter([v[0] for v in c], [v[1] for v in c], c = 'yellow')
    plt.scatter([v[0] for v in d], [v[1] for v in d], c = 'yellow')
    plt.scatter([v[0] for v in cents.values()], [v[1] for v in cents.values()], c = 'red')
    plt.show()

  null_ori = matrix.col((0,0,1)).axis_and_angle_as_unit_quaternion(0, deg=True)
  metro = { (0,): basis(null_ori, matrix.col((0,0,0))) }

  for q_id, q in six.iteritems(quadrants_trans):
    # calculate the center of the quadrant
    q_c = matrix.col((0,0,0))
    for s_id, s in six.iteritems(q):
      q_c += center(s)
    q_c /= len(q)

    metro[(0,q_id)] = basis(null_ori,q_c)

    for s_id, s in six.iteritems(q):
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

  if plot:
    print("Validating transofmation matrices set up correctly")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')

    keys = [key for key in metro if len(key) == 4]

    alen = pixel_size*asic_dimension[0]/2
    awid = pixel_size*asic_dimension[1]/2

    for key in sorted(keys):
      asic = [matrix.col((-alen,+awid,0,1)),
              matrix.col((+alen,+awid,0,1)),
              matrix.col((+alen,-awid,0,1)),
              matrix.col((-alen,-awid,0,1))]

      transform = metro[key[0:1]].as_homogenous_transformation() * \
                  metro[key[0:2]].as_homogenous_transformation() * \
                  metro[key[0:3]].as_homogenous_transformation() * \
                  metro[key[0:4]].as_homogenous_transformation()

      asic = [(transform * point)[0:2] for point in asic]

      ax.add_patch(Polygon(asic, closed=True, color='green', fill=False, hatch='/'))

    ax.set_xlim((-100, 100))
    ax.set_ylim((-100, 100))
    plt.show()

  return metro

def cbf_file_to_basis_dict(path):
  """ Maps a cbf file to a dictionary of tuples and basis objects, in the same form as the above from
  read_optical_metrology_from_flat_file
  @param path cbf file path """
  from dxtbx.format.Registry import Registry
  reader = Registry.find(path)
  instance = reader(path)
  return map_detector_to_basis_dict(instance.get_detector())

def map_detector_to_basis_dict(detector):
  """ Maps a dxtbx detector object representing a CSPAD CBF to a dictionary of tuples and basis objects,
  in the same form as the above from read_optical_metrology_from_flat_file
  @param detector cbf file path """

  root = detector.hierarchy()

  d = 0 # only allow one detector for now
  metro = {(d,):basis(panelgroup=root)}

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

        offset = matrix.col((-pix_size[0]*(img_size[0])/2,
                             +pix_size[1]*(img_size[1])/2,0))

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

def add_frame_specific_cbf_tables(cbf, wavelength, timestamp, trusted_ranges, diffrn_id = "DS1", is_xfel = True):
  """ Adds tables to cbf handle that won't already exsist if the cbf file is just a header
  @ param wavelength Wavelength in angstroms
  @ param timestamp String formatted timestamp for the image
  @ param trusted_ranges Array of trusted range tuples (min, max), one for each element """

  """Data items in the DIFFRN_RADIATION category describe
   the radiation used for measuring diffraction intensities,
   its collimation and monochromatization before the sample.

   Post-sample treatment of the beam is described by data
   items in the DIFFRN_DETECTOR category."""
  cbf.add_category("diffrn_radiation", ["diffrn_id","wavelength_id","probe"])
  cbf.add_row([diffrn_id,"WAVELENGTH1","x-ray"])

  """ Data items in the DIFFRN_RADIATION_WAVELENGTH category describe
   the wavelength of the radiation used in measuring the diffraction
   intensities. Items may be looped to identify and assign weights
   to distinct wavelength components from a polychromatic beam."""
  cbf.add_category("diffrn_radiation_wavelength", ["id","wavelength","wt"])
  cbf.add_row(["WAVELENGTH1",str(wavelength),"1.0"])

  """Data items in the DIFFRN_MEASUREMENT category record details
   about the device used to orient and/or position the crystal
   during data measurement and the manner in which the
   diffraction data were measured."""
  cbf.add_category("diffrn_measurement",["diffrn_id","id","number_of_axes","method","details"])
  cbf.add_row([diffrn_id,
    "INJECTION" if is_xfel else "unknown","0",
    "electrospray" if is_xfel else "unknown"
    "crystals injected by electrospray" if is_xfel else "unknown"])

  """ Data items in the DIFFRN_SCAN category describe the parameters of one
     or more scans, relating axis positions to frames."""
  cbf.add_category("diffrn_scan",["id","frame_id_start","frame_id_end","frames"])
  cbf.add_row(["SCAN1","FRAME1","FRAME1","1"])

  """Data items in the DIFFRN_SCAN_FRAME category describe
   the relationships of particular frames to scans."""
  cbf.add_category("diffrn_scan_frame",["frame_id","frame_number","integration_time","scan_id","date"])
  cbf.add_row(["FRAME1","1","0.0","SCAN1",timestamp])

  """ Data items in the ARRAY_INTENSITIES category record the
   information required to recover the intensity data from
   the set of data values stored in the ARRAY_DATA category."""
  # More detail here: http://www.iucr.org/__data/iucr/cifdic_html/2/cif_img.dic/Carray_intensities.html
  array_names = []
  cbf.find_category("diffrn_data_frame")
  while True:
    try:
      cbf.find_column("array_id")
      array_names.append(cbf.get_value())
      cbf.next_row()
    except Exception as e:
      assert "CBF_NOTFOUND" in e.message
      break

  cbf.add_category("array_intensities",["array_id","binary_id","linearity","gain","gain_esd","overload","undefined_value"])
  for i, array_name in enumerate(array_names):
    cbf.add_row([array_name,str(i+1),"linear","1.0","0.1",str(trusted_ranges[i][1]),str(trusted_ranges[i][0])])

def add_tiles_to_cbf(cbf, tiles, verbose = False):
  """
  Given a cbf handle, add the raw data and the necessary tables to support it
  """
  array_names = []
  cbf.find_category("diffrn_data_frame")
  while True:
    try:
      cbf.find_column("array_id")
      array_names.append(cbf.get_value())
      cbf.next_row()
    except Exception as e:
      assert "CBF_NOTFOUND" in e.message
      break

  tileisint = flex.bool()
  for tilekey in sorted(tiles.keys()):
    assert len(tiles[tilekey].focus()) == 3
    if isinstance(tiles[tilekey],flex.int):
      tileisint.append(True)
    elif isinstance(tiles[tilekey],flex.double):
      tileisint.append(False)
    else:
      raise TypeError("Ints or doubles are required")

  """ Data items in the ARRAY_STRUCTURE category record the organization and
  encoding of array data in the ARRAY_DATA category."""
  cbf.add_category("array_structure",["id","encoding_type","compression_type","byte_order"])
  for i, array_name in enumerate(array_names):
    if tileisint[i]:
      cbf.add_row([array_name,"signed 32-bit integer","packed","little_endian"])
    else:
      cbf.add_row([array_name,"signed 64-bit real IEEE","packed","little_endian"])

  """ Data items in the ARRAY_DATA category are the containers for the array data
  items described in the category ARRAY_STRUCTURE. """
  cbf.add_category("array_data",["array_id","binary_id","data"])

  if verbose:
    print("Compressing tiles...", end=' ')

  for i, (tilekey, array_name) in enumerate(zip(sorted(tiles.keys()), array_names)):
    focus = tiles[tilekey].focus()

    cbf.add_row([array_name,str(i+1)])

    binary_id = i+1
    data = tiles[tilekey].copy_to_byte_str()
    elements = len(tiles[tilekey])
    byteorder = "little_endian"
    dimfast = focus[2]
    dimmid = focus[1]
    dimslow = focus[0]
    padding = 0

    if tileisint[i]:
      elsize = 4
      elsigned = 1

      cbf.set_integerarray_wdims_fs(\
        pycbf.CBF_PACKED,
        binary_id,
        data,
        elsize,
        elsigned,
        elements,
        byteorder,
        dimfast,
        dimmid,
        dimslow,
        padding)
    else:
      elsize = 8

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

def copy_cbf_header(src_cbf, skip_sections = False):
  """ Given a cbf handle, copy the header tables only to a new cbf handle instance
  @param src_cbf cbf_handle instance that has the header information
  @param skip_sections If True, don't copy array array_structure_list_section,
  which may not always be present
  @return cbf_wrapper instance with the header information from the source """
  dst_cbf = cbf_wrapper()
  dst_cbf.new_datablock("dummy")

  categories = ["diffrn",
                "diffrn_source",
                "diffrn_detector",
                "diffrn_detector_axis",
                "diffrn_detector_element",
                "diffrn_data_frame",
                "array_structure_list"]
  if not skip_sections:
    categories.append("array_structure_list_section")
  categories.extend([
                "array_structure_list_axis",
                "axis",
                "diffrn_scan_axis",
                "diffrn_scan_frame_axis"])

  for cat in categories:
    src_cbf.find_category(cat)
    columns = []
    for i in range(src_cbf.count_columns()):
      src_cbf.select_column(i)
      columns.append(src_cbf.column_name())

    dst_cbf.add_category(cat, columns)

    for i in range(src_cbf.count_rows()):
      src_cbf.select_row(i)
      row = []
      for j in range(src_cbf.count_columns()):
        src_cbf.select_column(j)
        row.append(src_cbf.get_value())
      dst_cbf.add_row(row)

  return dst_cbf

def write_cspad_cbf(tiles, metro, metro_style, timestamp, cbf_root, wavelength, distance, verbose = True, header_only = False):
  cbf = get_cspad_cbf_handle(tiles, metro, metro_style, timestamp, cbf_root, wavelength, distance, verbose, header_only)

  cbf.write_widefile(cbf_root,pycbf.CBF,\
      pycbf.MIME_HEADERS|pycbf.MSG_DIGEST|pycbf.PAD_4K,0)

  if verbose:
    print("%s written"%cbf_root)

def get_cspad_cbf_handle(tiles, metro, metro_style, timestamp, cbf_root, wavelength, distance, verbose = True, header_only = False):
  assert metro_style in ['calibdir','flatfile','cbf']
  if metro_style == 'calibdir':
    metro = metro_phil_to_basis_dict(metro)

  # set up the metrology dictionary to include axis names, pixel sizes, and so forth
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
      basis.include_translation = False # don't include the translation in the rotation axis offset below, instead it will be
                                        # included in the axis_settings table below
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
      basis.dimension = asic_dimension
      basis.trusted_range = (cspad_min_trusted_value, cspad_saturated_value)
    else:
      assert False # shouldn't be reached as it would indicate more than four levels of hierarchy for this detector
    basis.axis_name = detector_axes_names[-1]


  # the data block is the root cbf node
  cbf=cbf_wrapper()
  cbf.new_datablock(os.path.splitext(os.path.basename(cbf_root))[0])

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

  # create a series of strings representing each asic.  Q here is for quadrant instead of panel.
  tilestrs = []
  tilequads = []
  tilekeys = sorted([key for key in metro if len(key) == 4])
  for tilekey in tilekeys:
    tilestrs.append("D%dQ%dS%dA%d"%tilekey)
    tilequads.append("D%dQ%d"%(tilekey[0],tilekey[1]))

  # create a series of strings representing each quad.  Q here is for quadrant instead of panel.
  quadstrs = []
  quadkeys = sorted([key for key in metro if len(key) == 2])
  for quadkey in quadkeys:
    quadstrs.append("D%dQ%d"%quadkey)

  """Data items in the DIFFRN_DETECTOR_ELEMENT category record
   the details about spatial layout and other characteristics
   of each element of a detector which may have multiple elements."""
  cbf.add_category("diffrn_detector_element",["id","detector_id"])

  for quadname in quadstrs:
    cbf.add_row(["ELE_" + quadname, "CSPAD_FRONT"])

  """Data items in the DIFFRN_DATA_FRAME category record
   the details about each frame of data."""
  cbf.add_category("diffrn_data_frame",["id","detector_element_id","array_id","binary_id"])

  for i, quadname in enumerate(quadstrs):
    cbf.add_row(["FRAME1","ELE_"+quadname,"ARRAY_"+quadname,"%d"%(i+1)])

  if not header_only:
    add_frame_specific_cbf_tables(cbf, wavelength, timestamp, [metro[k].trusted_range for k in tilekeys])

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

  root_basis = metro[(0,)]

  axis_settings.append(["AXIS_SOURCE" ,"FRAME1","0","0"])
  axis_settings.append(["AXIS_GRAVITY","FRAME1","0","0"])
  axis_settings.append([dname+"_X"    ,"FRAME1","0",str(root_basis.translation[0])])
  axis_settings.append([dname+"_Y"    ,"FRAME1","0",str(root_basis.translation[1])])
  axis_settings.append([dname+"_Z"    ,"FRAME1","0",str(root_basis.translation[2])])

  for key in sorted(metro):
    basis = metro[key]
    assert len(key) > 0 and len(key) <= 4

    cbf.add_frame_shift(basis, axis_settings)
    axis_names.append(basis.axis_name)

    if len(key) == 4:

      dim_pixel = basis.pixel_size
      dim_readout = basis.dimension

      # Add the two vectors for each asic that describe the fast and slow directions pixels should be laid out in real space
      offset_fast = -dim_pixel[0]*((dim_readout[0]) / 2)
      offset_slow = +dim_pixel[1]*((dim_readout[1]) / 2)

      aname = "D%dQ%dS%dA%d"%key

      cbf.add_row(["AXIS_"+ aname + "_S", "translation","detector",basis.axis_name    ,"0", "-1","0",str(offset_fast),str(offset_slow),"0.0", "detector_asic"])
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

  # find the asic sizes
  for tilekey in tilekeys:
    b = metro[tilekey]
    if not "x_dim" in locals():
      x_dim = b.dimension[0]
      y_dim = b.dimension[1]
    else:
      assert x_dim == b.dimension[0] and y_dim == b.dimension[1]

  sensor_quad_keys = [(key[0],key[1]) for key in metro if len(key) == 3]
  for quadkey in quadkeys:
    if not "z_dim" in locals():
      z_dim = sensor_quad_keys.count(quadkey)
    else:
      assert z_dim == sensor_quad_keys.count(quadkey)

  cbf.add_category("array_structure_list",["array_id","array_section_id","index","dimension","precedence","direction","axis_set_id"])
  for quadname in quadstrs:
    cbf.add_row(["ARRAY_"+quadname,".","1","%d"%(2*x_dim),"1","increasing","."])
    cbf.add_row(["ARRAY_"+quadname,".","2","%d"%y_dim,"2","increasing","."])
    cbf.add_row(["ARRAY_"+quadname,".","3","%d"%z_dim,"3","increasing","."])
  for tilename,tilekey,tilequad in zip(tilestrs,tilekeys,tilequads):
    b = metro[tilekey]
    cbf.add_row(["ARRAY_"+tilequad,"ARRAY_"+tilename,"1","%d"%b.dimension[0],"1","increasing","AXIS_"+tilename+"_F"])
    cbf.add_row(["ARRAY_"+tilequad,"ARRAY_"+tilename,"2","%d"%b.dimension[1],"2","increasing","AXIS_"+tilename+"_S"])

  """Data items in the ARRAY_STRUCTURE_LIST_SECTION category identify
     the dimension-by-dimension start, end and stride of each section of an
     array that is to be referenced."""
  cbf.add_category("array_structure_list_section",["array_id","id","index","start","end"])
  for tilename,tilekey,tilequad in zip(tilestrs,tilekeys,tilequads):
    b = metro[tilekey]
    if tilekey[3] == 0:
      x_dim = (1, b.dimension[0])
    elif tilekey[3] == 1:
      x_dim = (b.dimension[0]+1, 2*b.dimension[0])
    cbf.add_row(["ARRAY_"+tilequad,"ARRAY_"+tilename,"1","%d"%x_dim[0],"%d"%x_dim[1]])
    cbf.add_row(["ARRAY_"+tilequad,"ARRAY_"+tilename,"2","1","%d"%b.dimension[1]])
    cbf.add_row(["ARRAY_"+tilequad,"ARRAY_"+tilename,"3","%d"%(tilekey[2]+1),"%d"%(tilekey[2]+1)])

  """Data items in the ARRAY_STRUCTURE_LIST_AXIS category describe
     the physical settings of sets of axes for the centres of pixels that
     correspond to data points described in the
     ARRAY_STRUCTURE_LIST category."""
  cbf.add_category("array_structure_list_axis",["axis_set_id","axis_id","displacement","displacement_increment"])
  for tilename,tilekey in zip(tilestrs,tilekeys):
    cbf.add_row(["AXIS_"+tilename+"_F","AXIS_"+tilename+"_F","0.0",str(metro[tilekey].pixel_size[0])])
    cbf.add_row(["AXIS_"+tilename+"_S","AXIS_"+tilename+"_S","0.0",str(metro[tilekey].pixel_size[1])])

  if not header_only:
    add_tiles_to_cbf(cbf, tiles)

  return cbf
