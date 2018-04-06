from __future__ import absolute_import, division
#!/usr/bin/env python
#
# dxtbx.serialize.xds.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Authors: James Parkhurst, Richard Gildea, Graeme Winter
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

import glob
import sys
from scitbx import matrix
from rstbx.cftbx.coordinate_frame_helpers import align_reference_frame
from dxtbx.model.detector_helpers_types import detector_helpers_types
from dxtbx.sweep_filenames import template_regex

def to_imageset(input_filename, extra_filename=None):
  '''Get an image set from the xds input filename plus an extra filename

  Params:
      input_filename The XDS.INP file
      extra_filename A (G)XPARM.XDS, INTGRATE.HKL or XDS_ASCII.HKL file

  Returns:
      The imageset

  '''
  from iotbx.xds import xds_inp
  from dxtbx.imageset import ImageSetFactory
  import dxtbx

  # Read the input filename
  handle = xds_inp.reader()
  handle.read_file(input_filename)

  # Get the template
  template = handle.name_template_of_data_frames[0].replace('?', '#')
  image_range = handle.data_range
  detector_name = handle.detector

  if extra_filename is not None:
    # we can get all the extra dxtbx models from extra_filename
    check_format = False
  else:
    # we need the image files present to get the dxtbx models
    check_format = True

  # If an extra filename has been specified, try to load models
  if extra_filename:
    models = dxtbx.load(extra_filename)
    detector = models.get_detector()
    if detector_name.strip() in ('PILATUS', 'EIGER') or handle.silicon is not None:
      from dxtbx.model import ParallaxCorrectedPxMmStrategy
      from cctbx.eltbx import attenuation_coefficient
      if handle.silicon is None:
        table = attenuation_coefficient.get_table("Si")
        wavelength = models.get_beam().get_wavelength()
        mu = table.mu_at_angstrom(wavelength) / 10.0
      else:
        mu = handle.silicon
      t0 = handle.sensor_thickness
      for panel in detector:
        panel.set_px_mm_strategy(ParallaxCorrectedPxMmStrategy(mu, t0))
        panel.set_trusted_range(
          (handle.minimum_valid_pixel_value, handle.overload))
    beam = models.get_beam()
    detector = models.get_detector()
    goniometer = models.get_goniometer()
    scan = models.get_scan()
    scan.set_image_range(image_range)
  else:
    beam = None
    detector = None
    goniometer = None
    scan = None

  # Create the imageset
  imageset = ImageSetFactory.from_template(
    template,
    image_range=image_range,
    check_format=check_format,
    beam = beam,
    detector = detector,
    goniometer = goniometer,
    scan = scan)[0]

  # Return the imageset
  return imageset

def to_crystal(filename):
  ''' Get the crystal model from the xparm file

  Params:
      filename The xparm/or integrate filename

  Return:
      The crystal model

  '''
  from rstbx.cftbx.coordinate_frame_converter import \
      coordinate_frame_converter
  from cctbx.sgtbx import space_group, space_group_symbols

  # Get the real space coordinate frame
  cfc = coordinate_frame_converter(filename)
  real_space_a = cfc.get('real_space_a')
  real_space_b = cfc.get('real_space_b')
  real_space_c = cfc.get('real_space_c')
  sg = cfc.get('space_group_number')
  space_group = space_group(space_group_symbols(sg).hall())
  mosaicity = cfc.get('mosaicity')

  # Return the crystal model
  if (mosaicity is None):
    from dxtbx.model import Crystal
    crystal = Crystal(
        real_space_a=real_space_a,
        real_space_b=real_space_b,
        real_space_c=real_space_c,
        space_group=space_group)
  else:
    from dxtbx.model import MosaicCrystalKabsch2010
    crystal = MosaicCrystalKabsch2010(
        real_space_a=real_space_a,
        real_space_b=real_space_b,
        real_space_c=real_space_c,
        space_group=space_group)
    crystal.set_mosaicity(mosaicity)
  return crystal

def xds_detector_name(dxtbx_name):
  '''Translate from a xia2 name from the detector library to an XDS detector
  name.'''
  # http://xds.mpimf-heidelberg.mpg.de/html_doc/xds_parameters.html#DETECTOR=

  if 'pilatus' in dxtbx_name:
    return 'PILATUS'
  if 'eiger' in dxtbx_name:
    return 'PILATUS'
  if 'rayonix' in dxtbx_name:
    return 'CCDCHESS'
  if 'adsc' in dxtbx_name:
    return 'ADSC'
  if 'holton' in dxtbx_name:
    return 'ADSC'
  if 'saturn' in dxtbx_name:
    return 'SATURN'
  if 'raxis' in dxtbx_name:
    return 'RAXIS'
  if 'mar-345' in dxtbx_name:
    return 'MAR345'
  if 'mar' in dxtbx_name:
    return 'MAR'

  raise RuntimeError('detector %s unknown' % dxtbx_name)


class to_xds(object):
  '''A class to export contents of a Sweep as XDS.INP or XPARM.XDS.'''

  def __init__(self, sweep):
    self._sweep = sweep

    # detector dimensions in pixels
    self.detector_size = map(
      int,
      (max(panel.get_raw_image_offset()[0]+panel.get_image_size()[0]
           for panel in self.get_detector()),
       max(panel.get_raw_image_offset()[1]+panel.get_image_size()[1]
           for panel in self.get_detector())))
    self.fast, self.slow = self.detector_size

    if len(self.get_detector()) > 1:
      fast = self.get_detector()[0].get_parent_fast_axis()
      slow = self.get_detector()[0].get_parent_slow_axis()
      Rd = align_reference_frame(fast, (1,0,0), slow, (0,1,0))
      origin = Rd * matrix.col(self.get_detector()[0].get_parent_origin())
    else:
      fast = self.get_detector()[0].get_fast_axis()
      slow = self.get_detector()[0].get_slow_axis()
      Rd = align_reference_frame(fast, (1,0,0), slow, (0,1,0))
      origin = Rd * matrix.col(self.get_detector()[0].get_origin())

    self.detector_x_axis = (Rd * matrix.col(fast)).elems
    self.detector_y_axis = (Rd * matrix.col(slow)).elems

    F = Rd * matrix.col(fast)
    S = Rd * matrix.col(slow)
    N = F.cross(S)
    self.detector_normal = N.elems

    self.pixel_size = self.get_detector()[0].get_pixel_size() # assume all panels same pixel size

    centre = -(origin - origin.dot(N) * N)
    x = centre.dot(F)
    y = centre.dot(S)

    f, s = self.pixel_size
    self.detector_distance = origin.dot(N)
    # Need to add 0.5 because XDS seems to do centroids in fortran coords
    self.detector_origin = (x/f + 0.5, y/f + 0.5)

    self.imagecif_to_xds_transformation_matrix = Rd

    self.panel_limits = []
    self.panel_x_axis = []
    self.panel_y_axis = []
    self.panel_origin = []
    self.panel_distance = []
    self.panel_normal = []

    for panel_id, panel in enumerate(self.get_detector()):

      f = Rd * matrix.col(panel.get_fast_axis())
      s = Rd * matrix.col(panel.get_slow_axis())
      n = f.cross(s)

      xmin, ymin = panel.get_raw_image_offset()
      xmax = xmin + panel.get_image_size()[0]
      ymax = ymin + panel.get_image_size()[1]
      self.panel_limits.append((xmin+1, xmax, ymin+1, ymax))

      o = Rd * matrix.col(panel.get_origin())
      op = o.dot(n) * n
      d0 = matrix.col((-x, -y, self.detector_distance))
      orgsx = (op - o + d0).dot(f) / self.pixel_size[0] + xmin
      orgsy = (op - o + d0).dot(s) / self.pixel_size[1] + ymin
      panel_distance = op.dot(n) - d0.dot(n)

      # axes in local (i.e. detector) frame
      fl = matrix.col(panel.get_local_fast_axis())
      sl = matrix.col(panel.get_local_slow_axis())
      nl = fl.cross(sl)

      self.panel_x_axis.append(fl.elems)
      self.panel_y_axis.append(sl.elems)
      self.panel_normal.append(nl.elems)
      self.panel_origin.append((orgsx, orgsy))
      self.panel_distance.append(panel_distance)

    # Beam stuff
    self.wavelength = self.get_beam().get_wavelength()
    self.beam_vector = Rd * matrix.col(self.get_beam().get_direction())
    # just to make sure it is the correct length
    self.beam_vector = self.beam_vector.normalize() #/ self.wavelength
    self.beam_vector = (- self.beam_vector).elems

    # Scan and goniometer stuff
    self.starting_frame = self.get_scan().get_image_range()[0]
    self.starting_angle = self.get_scan().get_oscillation()[0]
    self.oscillation_range = self.get_scan().get_oscillation()[1]
    self.rotation_axis = (
        Rd * matrix.col(self.get_goniometer().get_rotation_axis())).elems
    return

  def get_detector(self):
    return self._sweep.get_detector()

  def get_goniometer(self):
    return self._sweep.get_goniometer()

  def get_beam(self):
    return self._sweep.get_beam()

  def get_scan(self):
    return self._sweep.get_scan()

  def get_template(self):
    try:
      return self._sweep.get_template()
    except AttributeError:
      return 'FIXME####.h5'


  def XDS_INP(self, out=None,
              space_group_number=None,
              real_space_a=None, real_space_b=None, real_space_c=None,
              job_card="XYCORR INIT COLSPOT IDXREF DEFPIX INTEGRATE CORRECT"):
    if out is None:
      out = sys.stdout

    assert [real_space_a, real_space_b, real_space_c].count(None) in (0,3)

    sensor = self.get_detector()[0].get_type()
    fast, slow = self.detector_size
    f, s = self.pixel_size
    df = int(1000 * f)
    ds = int(1000 * s)

    # FIXME probably need to rotate by pi about the X axis

    detector = xds_detector_name(
        detector_helpers_types.get(sensor, fast, slow, df, ds))
    trusted = self.get_detector()[0].get_trusted_range()

    print >> out, 'DETECTOR=%s MINIMUM_VALID_PIXEL_VALUE=%d OVERLOAD=%d' % \
          (detector, trusted[0] + 1, trusted[1])

    if detector == 'PILATUS':
      print >> out, 'SENSOR_THICKNESS= %.3f' % \
        self.get_detector()[0].get_thickness()
      if self.get_detector()[0].get_material():
        from cctbx.eltbx import attenuation_coefficient
        material = self.get_detector()[0].get_material()
        table = attenuation_coefficient.get_table(material)
        mu = table.mu_at_angstrom(self.wavelength) / 10.0
        print >> out, '!SENSOR_MATERIAL / THICKNESS %s %.3f' % \
          (material, self.get_detector()[0].get_thickness())
        print >> out, '!SILICON= %f' % mu

    print >> out, 'DIRECTION_OF_DETECTOR_X-AXIS= %.5f %.5f %.5f' % \
          self.detector_x_axis

    print >> out, 'DIRECTION_OF_DETECTOR_Y-AXIS= %.5f %.5f %.5f' % \
          self.detector_y_axis

    print >> out, 'NX=%d NY=%d QX=%.4f QY=%.4f' % (fast, slow, f, s)

    print >> out, 'DETECTOR_DISTANCE= %.6f' % self.detector_distance
    print >> out, 'ORGX= %.2f ORGY= %.2f' % self.detector_origin
    print >> out, 'ROTATION_AXIS= %.5f %.5f %.5f' % \
          self.rotation_axis
    print >> out, 'STARTING_ANGLE= %.3f' % \
          self.starting_angle
    print >> out, 'OSCILLATION_RANGE= %.3f' % \
          self.oscillation_range
    print >> out, 'X-RAY_WAVELENGTH= %.5f' % \
          self.wavelength
    print >> out, 'INCIDENT_BEAM_DIRECTION= %.3f %.3f %.3f' % \
          tuple([b / self.wavelength for b in self.beam_vector])

    # FIXME LATER
    if hasattr(self.get_beam(), "get_polarization_fraction"):
      print >> out, 'FRACTION_OF_POLARIZATION= %.3f' % \
          self.get_beam().get_polarization_fraction()
      print >> out, 'POLARIZATION_PLANE_NORMAL= %.3f %.3f %.3f' % \
          self.get_beam().get_polarization_normal()
    template = self.get_template()
    if template.endswith('master.h5'):
      master_file = template
      g = glob.glob(template.split('master.h5')[0]+'data_*[0-9].h5')
      assert len(g), 'No associated data files found for %s' % master_file
      template = template_regex(g[0])[0]
      template = master_file.split('master.h5')[0] + template.split('data_')[-1]
    print >> out, 'NAME_TEMPLATE_OF_DATA_FRAMES= %s' % \
        template.replace('#', '?')
    print >> out, 'TRUSTED_REGION= 0.0 1.41'
    for f0, s0, f1, s1 in self.get_detector()[0].get_mask():
      print >> out, 'UNTRUSTED_RECTANGLE= %d %d %d %d' % \
            (f0, f1 + 1, s0, s1 + 1)

    start_end = self.get_scan().get_image_range()

    if start_end[0] == 0:
      start_end = (1, start_end[1])

    print >> out, 'DATA_RANGE= %d %d' % start_end
    print >> out, 'JOB=%s' %job_card
    if space_group_number is not None:
      print >> out, 'SPACE_GROUP_NUMBER= %i' %space_group_number
    if [real_space_a, real_space_b, real_space_c].count(None) == 0:
      R = self.imagecif_to_xds_transformation_matrix
      unit_cell_a_axis = R * matrix.col(real_space_a)
      unit_cell_b_axis = R * matrix.col(real_space_b)
      unit_cell_c_axis = R * matrix.col(real_space_c)
      print >> out, "UNIT_CELL_A-AXIS= %.6f %.6f %.6f" %unit_cell_a_axis.elems
      print >> out, "UNIT_CELL_B-AXIS= %.6f %.6f %.6f" %unit_cell_b_axis.elems
      print >> out, "UNIT_CELL_C-AXIS= %.6f %.6f %.6f" %unit_cell_c_axis.elems

    if len(self.panel_x_axis) > 1:
      for panel_id, panel_x_axis in enumerate(self.panel_x_axis):

        print >> out
        print >> out, "!"
        print >> out, "! SEGMENT %d" %(panel_id+1)
        print >> out, "!"
        print >> out, 'SEGMENT= %d %d %d %d' % self.panel_limits[panel_id]
        print >> out, 'DIRECTION_OF_SEGMENT_X-AXIS= %.5f %.5f %.5f' % \
              panel_x_axis

        print >> out, 'DIRECTION_OF_SEGMENT_Y-AXIS= %.5f %.5f %.5f' % \
              self.panel_y_axis[panel_id]

        print >> out, 'SEGMENT_DISTANCE= %.3f' % self.panel_distance[panel_id]

        print >> out, 'SEGMENT_ORGX= %.2f SEGMENT_ORGY= %.2f' % self.panel_origin[panel_id]
        print >> out


  def xparm_xds(self, real_space_a, real_space_b, real_space_c,
                space_group, out=None):
    from cctbx import uctbx
    R = self.imagecif_to_xds_transformation_matrix
    unit_cell_a_axis = R * matrix.col(real_space_a)
    unit_cell_b_axis = R * matrix.col(real_space_b)
    unit_cell_c_axis = R * matrix.col(real_space_c)
    A_inv = matrix.sqr(unit_cell_a_axis.elems +
                       unit_cell_b_axis.elems +
                       unit_cell_c_axis.elems)
    metrical_matrix = (A_inv * A_inv.transpose()).as_sym_mat3()
    unit_cell = uctbx.unit_cell(metrical_matrix=metrical_matrix)
    from iotbx.xds import xparm
    writer = xparm.writer(
        self.starting_frame,
        self.starting_angle,
        self.oscillation_range,
        self.rotation_axis,
        self.wavelength,
        self.beam_vector,
        space_group,
        unit_cell.parameters(),
        unit_cell_a_axis.elems,
        unit_cell_b_axis.elems,
        unit_cell_c_axis.elems,
        None, # num_segments
        self.detector_size,
        self.pixel_size,
        self.detector_origin,
        self.detector_distance,
        self.detector_x_axis,
        self.detector_y_axis,
        self.detector_normal,
        segments=None,
        orientation=None)
    writer.show(out=out)
