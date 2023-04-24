from __future__ import division

from libtbx.utils import Sorry
from scitbx import matrix
from scitbx.array_family import flex
import serialtbx.util.time
import os

# The CAMP and CSpad counters are both 14 bits wide (Strüder et al
# 2010; Philipp et al., 2007), which means the physical limit is 2**14 - 1.
# However, in practice, when the pixels are in the low gain mode, after
# correcting by a gain value of around 6.87, the pixels tend to saturate
# around 90000. See xpp experiment xppe0314, run 184 as evidence.
cspad_saturated_value = 90000

# The dark average for the CSPAD detector is around 1100-1500. A pixel
# histogram of a minimum projection of an uncorrected (raw) light run shows
# a mostly flat tail up to ~800 ADU with a few bumps in the tail which
# represent true underloads. Assume a dark average of 1200 ADU. After dark
# subtraction, 800 - 1200 gives a minimum trusted value of -400. Reject
# pixels less than this.
cspad_min_trusted_value = -400

# The pixel size in mm.  The pixel size is fixed and square, with side
# length of 110 µm (Philipp et al., 2007).  XXX Should really clarify
# this with Sol and Chris.
#
# XXX Andor: 13.5 µm square, CAMP: 75 µm, square (Strüder et al.,
# 2010)
# Apr 14 2023: commenting this out and using the more accurate version below
#pixel_size = 110e-3


# need to define these here since it not defined in SLAC's metrology definitions
asic_dimension = (194,185)
asic_gap = 3
pixel_size = 0.10992

PSANA2_VERSION = 0
try:
  PSANA2_VERSION = os.environ.get('PSANA2_VERSION', 0)
except AttributeError:
  pass

from serialtbx.detector import basis, center
from serialtbx.detector.xtc import basis_from_geo

def read_slac_metrology(path = None, geometry = None, plot=False, include_asic_offset=False):
  if path is None and geometry is None:
    raise Sorry("Need to provide a geometry object or a path to a geometry file")

  if path is not None and geometry is not None:
    raise Sorry("Cannot provide a geometry object and a geometry file. Ambiguous")

  if geometry is None:
    try:
      from PSCalib.GeometryAccess import GeometryAccess
      geometry = GeometryAccess(path)
    except Exception:
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

def get_psana_corrected_data(psana_det, evt, use_default=False, dark=True, common_mode=None, apply_gain_mask=True,
                             gain_mask_value=None, per_pixel_gain=False, gain_mask=None, additional_gain_factor=None):
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
  @param gain_mask gain mask showing which pixels to apply gain mask value
  @param additional_gain_factor Additional gain factor. Pixels counts are divided by this number after all other
  corrections.
  @return Numpy array corrected as specified.
  """
  # order is pedestals, then common mode, then gain mask, then per pixel gain
  import numpy as np

  if PSANA2_VERSION:
      # in psana2, data are stored as raw, fex, etc so the selection
      # has to be given here when the detector interface is used.
      # for now, assumes cctbx uses "raw".
      psana_det = psana_det.raw

  if use_default:
    return psana_det.calib(evt)  # applies psana's complex run-dependent calibrations
  data = psana_det.raw_data(evt)
  if data is None:
    return

  data = data.astype(np.float64)
  if isinstance(dark, bool):
    if dark:
      if PSANA2_VERSION:
        data -= psana_det.pedestals()
      else:
        data -= psana_det.pedestals(evt)
  elif isinstance( dark, np.ndarray ):
    data -= dark

  if common_mode is not None and common_mode != "default":
    if common_mode == 'cspad_default':
      common_mode = (1,25,25,100,1)  # default parameters for CSPAD images
      psana_det.common_mode_apply(data, common_mode)
    elif common_mode == 'unbonded':
      common_mode = (5,0,0,0,0)  # unbonded pixels used for correction
      psana_det.common_mode_apply(data, common_mode)
    else:  # this is how it was before.. Though I think common_mode would need to be a tuple..
      psana_det.common_mode_apply(data, common_mode)
  if apply_gain_mask:
    if gain_mask is None:  # TODO: consider try/except here
      gain_mask = psana_det.gain_mask(evt) == 1
    if gain_mask_value is None:
      try:
        gain_mask_value = psana_det._gain_mask_factor
      except AttributeError:
        print("No gain set for psana detector, using gain value of 1, consider disabling gain in your phil file")
        gain_mask_value = 1
    data[gain_mask] = data[gain_mask]*gain_mask_value
  if per_pixel_gain: # TODO: test this
    data *= psana_det.gain()
  if additional_gain_factor is not None:
    data /= additional_gain_factor
  return data





def dpack(active_areas=None,
          address=None,
          beam_center_x=None,
          beam_center_y=None,
          ccd_image_saturation=None,
          data=None,
          distance=None,
          pixel_size=pixel_size,
          saturated_value=None,
          timestamp=None,
          wavelength=None,
          xtal_target=None,
          min_trusted_value=None):
  """XXX Check completeness.  Should fill in sensible defaults."""

  # Must have data.
  if data is None:
    return None

  # Create a time stamp of the current time if none was supplied.
  if timestamp is None:
    timestamp = serialtbx.util.time.timestamp()

  # For unknown historical reasons, the dictionary must contain both
  # CCD_IMAGE_SATURATION and SATURATED_VALUE items.
  if ccd_image_saturation is None:
    if saturated_value is None:
      ccd_image_saturation = cspad_saturated_value
    else:
      ccd_image_saturation = saturated_value
  if saturated_value is None:
    saturated_value = ccd_image_saturation

  # Use a minimum value if provided for the pixel range
  if min_trusted_value is None:
    min_trusted_value = cspad_min_trusted_value

  # By default, the beam center is the center of the image.  The slow
  # (vertical) and fast (horizontal) axes correspond to x and y,
  # respectively.
  if beam_center_x is None:
    beam_center_x = pixel_size * data.focus()[1] / 2
  if beam_center_y is None:
    beam_center_y = pixel_size * data.focus()[0] / 2

  # By default, the entire detector image is an active area.  There is
  # no sensible default for distance nor wavelength.  XXX But setting
  # wavelength to zero may be disastrous?
  if active_areas is None:
    # XXX Verify order with non-square detector
    active_areas = flex.int((0, 0, data.focus()[0], data.focus()[1]))
  if distance is None:
    distance = 0
  if wavelength is None:
    wavelength = 0

  # The size must match the image dimensions.  The length along the
  # slow (vertical) axis is SIZE1, the length along the fast
  # (horizontal) axis is SIZE2.
  return {'ACTIVE_AREAS': active_areas,
          'BEAM_CENTER_X': beam_center_x,
          'BEAM_CENTER_Y': beam_center_y,
          'CCD_IMAGE_SATURATION': ccd_image_saturation,
          'DATA': data,
          'DETECTOR_ADDRESS': address,
          'DISTANCE': distance,
          'PIXEL_SIZE': pixel_size,
          'SATURATED_VALUE': saturated_value,
          'MIN_TRUSTED_VALUE': min_trusted_value,
          'SIZE1': data.focus()[0],
          'SIZE2': data.focus()[1],
          'TIMESTAMP': timestamp,
          'SEQUENCE_NUMBER': 0, # XXX Deprecated
          'WAVELENGTH': wavelength,
          'xtal_target': xtal_target}
