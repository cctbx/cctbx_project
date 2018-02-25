from __future__ import absolute_import, division
#!/usr/bin/env python
# detector.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# A model for the detector for the "updated experimental model" project
# documented in internal ticket #1555. This is not designed to be used outside
# of the XSweep classes. N.B. this should probably be generalized for non
# flat detectors, or composite detectors constructed from a number of flat
# elements.

import pycbf
from scitbx import matrix
from dxtbx_model_ext import Panel, Detector
from dxtbx_model_ext import SimplePxMmStrategy, ParallaxCorrectedPxMmStrategy
from dxtbx.model.detector_helpers import detector_helper_sensors
from dxtbx.model.detector_helpers import find_undefined_value
import libtbx.phil

import os
if 'DXTBX_OVERLOAD_SCALE' in os.environ:
  dxtbx_overload_scale = float(os.environ['DXTBX_OVERLOAD_SCALE'])
else:
  dxtbx_overload_scale = 1

detector_phil_scope = libtbx.phil.parse('''
  detector
    .expert_level = 1
    .short_caption = "Detector overrides"
  {
    panel
      .multiple = True
    {
      id = 0
        .type = int
        .help = "The panel number"
        .short_caption = "Panel ID"

      name = None
        .type = str
        .help = "Override the panel name"
        .short_caption = "Panel name"

      type = None
        .type = str
        .help = "Override the panel type"
        .short_caption = "Panel type"

      gain = None
        .type = float(value_min=0)
        .help = "The gain of the detector panel"
        .short_caption = "Gain value"

      pixel_size = None
        .type = floats(size=2)
        .help = "Override the panel pixel size"
        .short_caption = "Panel pixel size"

      image_size = None
        .type = ints(size=2)
        .help = "Override the panel image size"
        .short_caption= "Panel image size"

      trusted_range = None
        .type = floats(size=2)
        .help = "Override the panel trusted range"
        .short_caption = "Panel trusted range"

      thickness = None
        .type = float
        .help = "Override the panel thickness"
        .short_caption = "Panel thickness"

      material = None
        .type = str
        .help = "Override the panel material"
        .short_caption = "Panel material"

      fast_axis = None
        .type = floats(size=3)
        .help = "Override the panel fast axis. Requires slow_axis and origin."
        .short_caption = "Panel fast axis direction"

      slow_axis = None
        .type = floats(size=3)
        .help = "Override the panel slow axis. Requires fast_axis and origin."
        .short_caption = "Panel slow axis direction"

      origin = None
        .type = floats(size=3)
        .help = "Override the panel origin. Requires fast_axis and slow_axis."
        .short_caption = "Panel origin vector"

      parallax_correction = None
        .type = bool
        .help = "Enable parallax correction. By default in overwrite mode, the"
                "value of None does nothing."
        .short_caption = "Enable parallax correction"
    }

    hierarchy
      {

      name = None
        .type = str
        .help = "Override the group name"
        .short_caption = "group name"

      fast_axis = None
        .type = floats(size=3)
        .help = "Override the panel fast axis. Requires slow_axis and origin."
        .short_caption = "Panel fast axis direction"

      slow_axis = None
        .type = floats(size=3)
        .help = "Override the panel slow axis. Requires fast_axis and origin."
        .short_caption = "Panel slow axis direction"

      origin = None
        .type = floats(size=3)
        .help = "Override the panel origin. Requires fast_axis and slow_axis."
        .short_caption = "Panel origin vector"

      group
        .multiple = True
      {
        id = None
          .type = ints
          .help = "The group identifier specifying the place in the hierarchy"
          .short_caption = "Group ID"

        name = None
          .type = str
          .help = "Override the group name"
          .short_caption = "Group name"

        fast_axis = None
        .type = floats(size=3)
          .help = "Override the group fast axis. Requires slow_axis and origin."
          .short_caption = "Group fast axis direction"

        slow_axis = None
          .type = floats(size=3)
          .help = "Override the group slow axis. Requires fast_axis and origin."
          .short_caption = "Group slow axis direction"

        origin = None
          .type = floats(size=3)
          .help = "Override the group origin. Requires fast_axis and slow_axis."
          .short_caption = "Group origin vector"

        panel = None
          .multiple = True
          .type = int
          .help = "The panel id"
          .short_caption = "Panel ID"
      }
    }

    mosflm_beam_centre = None
      .type = floats(size=2)
      .help = "Override the beam centre from the image headers, following "
              "the mosflm convention."
      .short_caption = "Beam centre coordinates (mm, mm) using the Mosflm convention"

    distance = None
      .type = float
      .help = "The detector distance (used when mosflm_beam_centre is set)"
      .short_caption = "Detector distance"

    slow_fast_beam_centre = None
      .type = ints(size_min=2, size_max=3)
      .help = "Override the beam centre from the image headers, following "
              "the slow/fast pixel convention used by dials.image_viewer."
              "The first two values are the slow and fast pixel coordinate."
              "If the third is supplied it specifies a panel number."
      .short_caption = "Beam centre coordinates (px slow, px fast, [panel id])"
  }
''')


class DetectorFactory:
  '''A factory class for detector objects, which will encapsulate standard
  detector designs to make it a little easier to get started with these. In
  cases where a CBF image is provided a full description can be used, in
  other cases assumptions will be made about the experiment configuration.
  In all cases information is provided in the CBF coordinate frame.'''

  def __init__(self):
    pass

  @staticmethod
  def generate_from_phil(params, beam=None):
    '''
    Create a new detector model from phil parameters

    '''
    from cctbx.eltbx import attenuation_coefficient
    detector = Detector()

    # Create a list of panels
    panel_list = {}
    for panel_params in params.detector.panel:
      panel = Panel()
      if panel_params.name is not None:
        panel.set_name(panel_params.name)
      if panel_params.type is not None:
        panel.set_type(panel_params.type)
      if panel_params.gain is not None:
        panel.set_gain(panel_params.gain)
      if panel_params.pixel_size is not None:
        panel.set_pixel_size(panel_params.pixel_size)
      else:
        raise RuntimeError('No pixel size set')
      if panel_params.image_size is not None:
        panel.set_image_size(panel_params.image_size)
      else:
        raise RuntimeError('No image size set')
      if panel_params.trusted_range is not None:
        panel.set_trusted_range(panel_params.trusted_range)
      else:
        raise RuntimeError('No trusted range set')
      if panel_params.thickness is not None:
        panel.set_thickness(panel_params.thickness)
      if panel_params.material is not None:
        panel.set_material(panel_params.material)
      if panel_params.parallax_correction is True:
        if panel_params.material is None:
          raise RuntimeError("No material for parallax correction")
        if panel_params.thickness is None:
          raise RuntimeError("No thickness for parallax correction")
        if beam is None:
          raise RuntimeError("No beam for parallax correction")
        table = attenuation_coefficient.get_table(panel_params.material)
        mu = table.mu_at_angstrom(beam.get_wavelength()) / 10.0
        t0 = panel_params.thickness
        panel.set_px_mm_strategy(ParallaxCorrectedPxMmStrategy(mu, t0))
      if panel_params.fast_axis is None:
        panel_params.fast_axis = (1, 0, 0)
      if panel_params.slow_axis is None:
        panel_params.slow_axis = (0, 1, 0)
      if panel_params.origin is None:
        panel_params.origin = (0, 0, 0)
      panel.set_local_frame(
        panel_params.fast_axis,
        panel_params.slow_axis,
        panel_params.origin)
      if panel_params.id in panel_list:
        raise RuntimeError('Multiple panels with id=%d' % panel_params.id)
      panel_list[panel_params.id] = panel

    # Create the hierarchy
    panel_counter = 0
    root = detector.hierarchy()
    if params.detector.hierarchy.name is not None:
      root.set_name(params.detector.hierarchy.name)
    if params.detector.hierarchy.fast_axis is None:
      params.detector.hierarchy.fast_axis = (1, 0, 0)
    if params.detector.hierarchy.slow_axis is None:
      params.detector.hierarchy.slow_axis = (0, 1, 0)
    if params.detector.hierarchy.origin is None:
      params.detector.hierarchy.origin = (0, 0, 0)
    root.set_frame(
      params.detector.hierarchy.fast_axis,
      params.detector.hierarchy.slow_axis,
      params.detector.hierarchy.origin)
    def get_parent(node, index):
      if len(index) == 0:
        return node
      return get_parent(node[index[0]], index[1:])
    for group_params in params.detector.hierarchy.group:
      parent = get_parent(root, group_params.id[:-1])
      assert len(parent) == group_params.id[-1]
      group = parent.add_group()
      if group_params.name is not None:
        group.set_name(group_params.name)
      if group_params.fast_axis is None:
        group_params.fast_axis = (1, 0, 0)
      if group_params.slow_axis is None:
        group_params.slow_axis = (0, 1, 0)
      if group_params.origin is None:
        group_params.origin = (0, 0, 0)
      group.set_local_frame(
        group_params.fast_axis,
        group_params.slow_axis,
        group_params.origin)
      for panel_id in group_params.panel:
        assert panel_id == panel_counter
        group.add_panel(panel_list[panel_id])
        panel_counter += 1
    if panel_counter == 0:
      for panel_id in range(max(panel_list.keys())+1):
        detector.add_panel(panel_list[panel_id])
    elif panel_counter != len(panel_list):
      raise RuntimeError("Inconsistent number of panels in hierarchy")

    # Return detector
    return detector

  @staticmethod
  def overwrite_from_phil(params, detector, beam=None):
    '''
    Overwrite from phil parameters

    '''
    from cctbx.eltbx import attenuation_coefficient
    # Override any panel parameters
    for panel_params in params.detector.panel:
      panel = detector[panel_params.id]
      if panel_params.name is not None:
        panel.set_name(panel_params.name)
      if panel_params.type is not None:
        panel.set_type(panel_params.type)
      if panel_params.gain is not None:
        panel.set_gain(panel_params.gain)
      if panel_params.pixel_size is not None:
        panel.set_pixel_size(panel_params.pixel_size)
      if panel_params.image_size is not None:
        panel.set_image_size(panel_params.image_size)
      if panel_params.trusted_range is not None:
        panel.set_trusted_range(panel_params.trusted_range)
      if panel_params.thickness is not None:
        panel.set_thickness(panel_params.thickness)
      if panel_params.material is not None:
        panel.set_material(panel_params.material)
      if panel_params.parallax_correction is None:
        if isinstance(panel.get_px_mm_strategy(), ParallaxCorrectedPxMmStrategy):
          panel_params.parallax_correction = True
        else:
          panel_params.parallax_correction = False
      if panel_params.parallax_correction is True:
        if beam is None:
          raise RuntimeError("No beam for parallax correction")
        table = attenuation_coefficient.get_table(panel.get_material())
        mu = table.mu_at_angstrom(beam.get_wavelength()) / 10.0
        t0 = panel.get_thickness()
        panel.set_px_mm_strategy(ParallaxCorrectedPxMmStrategy(mu, t0))
      else:
        panel.set_px_mm_strategy(SimplePxMmStrategy())
      axes = [panel_params.fast_axis,
              panel_params.slow_axis,
              panel_params.origin]
      if axes.count(None) != 3:
        if panel_params.fast_axis is None:
          panel_params.fast_axis = panel.get_local_fast_axis()
        if panel_params.slow_axis is None:
          panel_params.slow_axis = panel.get_local_slow_axis()
        if panel_params.origin is None:
          panel_params.origin = panel.get_local_origin()
        panel.set_local_frame(
          panel_params.fast_axis,
          panel_params.slow_axis,
          panel_params.origin)

    # Create the hierarchy
    if params.detector.hierarchy is not None:
      root = detector.hierarchy()
      if params.detector.hierarchy.name is not None:
        root.set_name(params.detector.hierarchy.name)
      if (params.detector.hierarchy.fast_axis is not None or
          params.detector.hierarchy.slow_axis is not None or
          params.detector.hierarchy.origin is not None):
        if params.detector.hierarchy.fast_axis is None:
          params.detector.hierarchy.fast_axis = root.get_fast_axis()
        if params.detector.hierarchy.slow_axis is None:
          params.detector.hierarchy.slow_axis = root.get_slow_axis()
        if params.detector.hierarchy.origin is None:
          params.detector.hierarchy.origin = root.get_origin()
        root.set_frame(
          params.detector.hierarchy.fast_axis,
          params.detector.hierarchy.slow_axis,
          params.detector.hierarchy.origin)
      def get_group(node, index):
        if len(index) == 0:
          return node
        return get_group(node[index[0]], index[1:])
      for group_params in params.detector.hierarchy.group:
        group = get_group(root, group_params.id)
        if group_params.name is not None:
          group.set_name(group_params.name)
        if (group_params.fast_axis is not None or
            group_params.slow_axis is not None or
            group_params.origin is not None):
          if group_params.fast_axis is None:
            group_params.fast_axis = group.get_local_fast_axis()
          if group_params.slow_axis is None:
            group_params.slow_axis = group.get_local_slow_axis()
          if group_params.origin is None:
            group_params.origin = group.get_local_origin()
          group.set_local_frame(
            group_params.fast_axis,
            group_params.slow_axis,
            group_params.origin)
        if len(group_params.panel) != 0:
          raise RuntimeError("Can't reassign panels in groups")

    # Return the detector
    return detector

  @staticmethod
  def from_phil(params, reference=None, beam=None):
    '''
    Convert phil parameters into detector model

    '''
    from dxtbx.model.detector_helpers import set_detector_distance
    from dxtbx.model.detector_helpers import set_mosflm_beam_centre
    from dxtbx.model.detector_helpers import set_slow_fast_beam_centre_mm

    # Check the input. If no reference detector is provided then
    # Create the detector model from scratch from the parameters
    if reference is None:
      detector = DetectorFactory.generate_from_phil(params, beam)
    else:
      detector = DetectorFactory.overwrite_from_phil(params, reference, beam)

    # If the distance is set
    if params.detector.distance is not None:
      set_detector_distance(detector, params.detector.distance)

    # If the mosflm beam centre is set then update
    if params.detector.mosflm_beam_centre is not None:
      assert beam is not None
      set_mosflm_beam_centre(
        detector,
        beam,
        params.detector.mosflm_beam_centre)

    # If the slow fast beam centre is set then update
    if params.detector.slow_fast_beam_centre is not None:
      panel_id = 0
      if len(params.detector.slow_fast_beam_centre) > 2:
        panel_id = params.detector.slow_fast_beam_centre[2]
      if panel_id >= len(detector):
        raise Sorry('Detector does not have panel index {0}'.format(panel_id))
      px_size_f, px_size_s = detector[0].get_pixel_size()
      slow_fast_beam_centre_mm = (
          params.detector.slow_fast_beam_centre[0] * px_size_s,
          params.detector.slow_fast_beam_centre[1] * px_size_f)
      assert beam is not None
      set_slow_fast_beam_centre_mm(
        detector,
        beam,
        slow_fast_beam_centre_mm,
        panel_id=panel_id)

    # Return the model
    return detector

  @staticmethod
  def from_dict(d, t=None):
    ''' Convert the dictionary to a detector model

    Params:
        d The dictionary of parameters
        t The template dictionary to use

    Returns:
        The detector model

    '''
    from dxtbx.model import Detector

    # If None, return None
    if d == None:
      if t == None: return None
      else: return from_dict(t, None)
    elif t != None:
      if isinstance(d, list):
        d = { 'panels' : d }
      d2 = dict(t.items() + d.items())
    else:
      if isinstance(d, list):
        d = { 'panels' : d }

    # Create the model from the dictionary
    return Detector.from_dict(d)

  @staticmethod
  def make_detector(stype, fast_axis, slow_axis, origin,
                    pixel_size, image_size, trusted_range = (0.0, 0.0),
                    px_mm=None, name="Panel", thickness=0.0, material='',
                    mu=0.0, gain=None, identifier=""):
    """Ensure all types are correct before creating c++ detector class."""

    if px_mm is None:
      px_mm = SimplePxMmStrategy()
    try:
      d = Detector()
      p = d.add_panel()
      p.set_type(str(stype))
      p.set_name(str(name))
      p.set_local_frame(
          tuple(map(float, fast_axis)),
          tuple(map(float, slow_axis)),
          tuple(map(float, origin)))
      p.set_pixel_size(tuple(map(float, pixel_size)))
      p.set_image_size(tuple(map(int, image_size)))
      p.set_trusted_range(tuple(map(float, trusted_range)))
      p.set_thickness(thickness)
      p.set_material(material)
      p.set_px_mm_strategy(px_mm)
      p.set_identifier(identifier)
      if gain is not None:
        p.set_gain(gain)
    except Exception as e:
      print e
      raise e
    return d

  @staticmethod
  def simple(sensor, distance, beam_centre, fast_direction, slow_direction,
             pixel_size, image_size, trusted_range = (0.0, 0.0), mask = [],
             px_mm=None, mu=0.0, gain=None, identifier=""):
    '''Construct a simple detector at a given distance from the sample
    along the direct beam presumed to be aligned with -z, offset by the
    beam centre - the directions of which are given by the fast and slow
    directions, which are themselves given as +x, +y, -x, -y. The pixel
    size is given in mm in the fast and slow directions and the image size
    is given in pixels in the same order. Everything else is the same as
    for the main reference frame.'''

    assert(fast_direction in ['-x', '+y', '+x', '-y'])
    assert(slow_direction in ['-x', '+y', '+x', '-y'])

    assert(fast_direction[1] != slow_direction[1])

    direction_map = {
        '+x':(1.0, 0.0, 0.0),
        '-x':(-1.0, 0.0, 0.0),
        '+y':(0.0, 1.0, 0.0),
        '-y':(0.0, -1.0, 0.0)
        }

    fast = matrix.col(direction_map[fast_direction])
    slow = matrix.col(direction_map[slow_direction])

    origin = matrix.col((0, 0, -1)) * distance - \
             fast * beam_centre[0] - slow * beam_centre[1]

    detector = DetectorFactory.make_detector(
        DetectorFactory.sensor(sensor),
        fast,
        slow,
        origin,
        pixel_size,
        image_size,
        trusted_range,
        px_mm,
        mu=mu,
        gain=gain,
        identifier=identifier)
    detector[0].mask = mask
    return detector

  @staticmethod
  def two_theta(sensor, distance, beam_centre, fast_direction,
                slow_direction, two_theta_direction, two_theta_angle,
                pixel_size, image_size, trusted_range = (0.0, 0.0),
                mask = [], px_mm = None, mu=0.0, gain=None, identifier=""):
    '''Construct a simple detector at a given distance from the sample
    along the direct beam presumed to be aligned with -z, offset by the
    beam centre - the directions of which are given by the fast and slow
    directions, which are themselves given as +x, +y, -x, -y. The pixel
    size is given in mm in the fast and slow directions and the image size
    is given in pixels in the same order. Everything else is the same as
    for the main reference frame. Also given are the direction of the
    two-theta axis and the angle in degrees by which the detector is
    moved.'''

    assert(fast_direction in ['-x', '+y', '+x', '-y'])
    assert(slow_direction in ['-x', '+y', '+x', '-y'])
    assert(two_theta_direction in ['-x', '+y', '+x', '-y'])

    assert(fast_direction[1] != slow_direction[1])

    direction_map = {
        '+x':(1.0, 0.0, 0.0),
        '-x':(-1.0, 0.0, 0.0),
        '+y':(0.0, 1.0, 0.0),
        '-y':(0.0, -1.0, 0.0)
        }

    fast = matrix.col(direction_map[fast_direction])
    slow = matrix.col(direction_map[slow_direction])

    origin = matrix.col((0, 0, -1)) * distance - \
             fast * beam_centre[0] - slow * beam_centre[1]

    two_theta = matrix.col(direction_map[two_theta_direction])

    R = two_theta.axis_and_angle_as_r3_rotation_matrix(two_theta_angle,
                                                       deg = True)

    detector = DetectorFactory.make_detector(
        DetectorFactory.sensor(sensor),
        (R * fast), (R * slow), (R * origin), pixel_size,
        image_size, trusted_range, px_mm, mu=mu, gain=gain, identifier=identifier)

    detector.mask = mask
    return detector

  @staticmethod
  def complex(sensor, origin, fast, slow, pixel, size,
              trusted_range = (0.0, 0.0), px_mm = None, gain = None,
              identifier=""):
    '''A complex detector model, where you know exactly where everything
    is. This is useful for implementation of the Rigaku Saturn header
    format, as that is exactly what is in there. Origin, fast and slow are
    vectors in the CBF reference frame, pixel is the dimensions as a tuple
    as is size.'''

    assert(len(origin) == 3)
    assert(len(fast) == 3)
    assert(len(slow) == 3)
    assert(len(pixel) == 2)
    assert(len(size) == 2)

    return DetectorFactory.make_detector(
            DetectorFactory.sensor(sensor),
            fast, slow, origin, pixel, size, trusted_range, px_mm, gain=gain,
            identifier=identifier)

  @staticmethod
  def imgCIF(cif_file, sensor):
    '''Initialize a detector model from an imgCIF file.'''

    cbf_handle = pycbf.cbf_handle_struct()
    cbf_handle.read_file(cif_file, pycbf.MSG_DIGEST)

    return DetectorFactory.imgCIF_H(cbf_handle, sensor)

  @staticmethod
  def imgCIF_H(cbf_handle, sensor):
    '''Initialize a detector model from an imgCIF file handle, where it
    is assumed that the file has already been read.'''

    cbf_detector = cbf_handle.construct_detector(0)

    pixel = (cbf_detector.get_inferred_pixel_size(1),
             cbf_detector.get_inferred_pixel_size(2))

    origin = tuple(cbf_detector.get_pixel_coordinates(0, 0))
    fast = cbf_detector.get_detector_axis_fast()
    slow = cbf_detector.get_detector_axis_slow()
    size = tuple(reversed(cbf_handle.get_image_size(0)))

    try:
      underload = find_undefined_value(cbf_handle)
      overload = cbf_handle.get_overload(0)
      trusted_range = (underload, overload * dxtbx_overload_scale)
    except Exception:
      trusted_range = (0.0, 0.0)

    cbf_detector.__swig_destroy__(cbf_detector)
    del(cbf_detector)

    return DetectorFactory.make_detector(
                  DetectorFactory.sensor(sensor),
                  fast, slow, origin, pixel, size, trusted_range)

  @staticmethod
  def sensor(name):
    '''Return the correct sensor token for a given name, for example:

    ccd, CCD
    image_plate, IMAGE_PLATE
    pad, PAD

    to the appropriate static token which will be used as a handle
    everywhere else in this. Also allow existing token to be passed in.'''

    if detector_helper_sensors.check_sensor(name):
      return name

    if name.upper() == 'PAD':
      return detector_helper_sensors.SENSOR_PAD
    elif name.upper() == 'CCD':
      return detector_helper_sensors.SENSOR_CCD
    elif name.upper() == 'IMAGE_PLATE':
      return detector_helper_sensors.SENSOR_IMAGE_PLATE
    elif name.upper() == 'UNKNOWN':
      return detector_helper_sensors.SENSOR_UNKNOWN

    raise RuntimeError('name %s not known' % name)
