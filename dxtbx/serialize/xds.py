from __future__ import absolute_import, division, print_function

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
  if template.endswith("h5"):
    template = template.replace("######", "master")
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
              job_card="XYCORR INIT COLSPOT IDXREF DEFPIX INTEGRATE CORRECT",
              as_str=False):
    if out is None:
      out = sys.stdout

    # horrible hack to allow returning result as string; would be nice if the
    # structure of this was to make the XDS.INP in memory then print it...
    # see also show() method on things.
    str_result = []
    if as_str:
      def print(str, file=None):
        str_result.append(str)
    else:
      try:
        import __builtin__
      except ImportError:
        import builtins as __builtin__
      print = __builtin__.print

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

    print('DETECTOR=%s MINIMUM_VALID_PIXEL_VALUE=%d OVERLOAD=%d' % \
          (detector, trusted[0] + 1, trusted[1]), file=out)

    if detector == 'PILATUS':
      print('SENSOR_THICKNESS= %.3f' % \
        self.get_detector()[0].get_thickness(), file=out)
      if self.get_detector()[0].get_material():
        from cctbx.eltbx import attenuation_coefficient
        material = self.get_detector()[0].get_material()
        table = attenuation_coefficient.get_table(material)
        mu = table.mu_at_angstrom(self.wavelength) / 10.0
        print('!SENSOR_MATERIAL / THICKNESS %s %.3f' % \
          (material, self.get_detector()[0].get_thickness()), file=out)
        print('!SILICON= %f' % mu, file=out)

    print('DIRECTION_OF_DETECTOR_X-AXIS= %.5f %.5f %.5f' % \
          self.detector_x_axis, file=out)

    print('DIRECTION_OF_DETECTOR_Y-AXIS= %.5f %.5f %.5f' % \
          self.detector_y_axis, file=out)

    print('NX=%d NY=%d QX=%.4f QY=%.4f' % (fast, slow, f, s), file=out)

    print('DETECTOR_DISTANCE= %.6f' % self.detector_distance, file=out)
    print('ORGX= %.2f ORGY= %.2f' % self.detector_origin, file=out)
    print('ROTATION_AXIS= %.5f %.5f %.5f' % \
          self.rotation_axis, file=out)
    print('STARTING_ANGLE= %.3f' % \
          self.starting_angle, file=out)
    print('OSCILLATION_RANGE= %.3f' % \
          self.oscillation_range, file=out)
    print('X-RAY_WAVELENGTH= %.5f' % \
          self.wavelength, file=out)
    print('INCIDENT_BEAM_DIRECTION= %.3f %.3f %.3f' % \
          tuple([b / self.wavelength for b in self.beam_vector]), file=out)

    # FIXME LATER
    if hasattr(self.get_beam(), "get_polarization_fraction"):
      print('FRACTION_OF_POLARIZATION= %.3f' % \
          self.get_beam().get_polarization_fraction(), file=out)
      print('POLARIZATION_PLANE_NORMAL= %.3f %.3f %.3f' % \
          self.get_beam().get_polarization_normal(), file=out)
    template = self.get_template()
    if template.endswith('master.h5'):
      master_file = template
      g = glob.glob(template.split('master.h5')[0]+'data_*[0-9].h5')
      assert g, 'No associated data files found for %s' % master_file
      template = template_regex(g[0])[0]
      template = master_file.split('master.h5')[0] + template.split('data_')[-1]
    print('NAME_TEMPLATE_OF_DATA_FRAMES= %s' % \
        template.replace('#', '?'), file=out)
    print('TRUSTED_REGION= 0.0 1.41', file=out)
    for f0, s0, f1, s1 in self.get_detector()[0].get_mask():
      print('UNTRUSTED_RECTANGLE= %d %d %d %d' % \
            (f0, f1 + 1, s0, s1 + 1), file=out)

    start_end = self.get_scan().get_image_range()

    if start_end[0] == 0:
      start_end = (1, start_end[1])

    print('DATA_RANGE= %d %d' % start_end, file=out)
    print('JOB=%s' %job_card, file=out)
    if space_group_number is not None:
      print('SPACE_GROUP_NUMBER= %i' %space_group_number, file=out)
    if [real_space_a, real_space_b, real_space_c].count(None) == 0:
      R = self.imagecif_to_xds_transformation_matrix
      unit_cell_a_axis = R * matrix.col(real_space_a)
      unit_cell_b_axis = R * matrix.col(real_space_b)
      unit_cell_c_axis = R * matrix.col(real_space_c)
      print("UNIT_CELL_A-AXIS= %.6f %.6f %.6f" %unit_cell_a_axis.elems, file=out)
      print("UNIT_CELL_B-AXIS= %.6f %.6f %.6f" %unit_cell_b_axis.elems, file=out)
      print("UNIT_CELL_C-AXIS= %.6f %.6f %.6f" %unit_cell_c_axis.elems, file=out)

    if len(self.panel_x_axis) > 1:
      for panel_id, panel_x_axis in enumerate(self.panel_x_axis):

        print(file=out)
        print("!", file=out)
        print("! SEGMENT %d" %(panel_id+1), file=out)
        print("!", file=out)
        print('SEGMENT= %d %d %d %d' % self.panel_limits[panel_id], file=out)
        print('DIRECTION_OF_SEGMENT_X-AXIS= %.5f %.5f %.5f' % \
              panel_x_axis, file=out)

        print('DIRECTION_OF_SEGMENT_Y-AXIS= %.5f %.5f %.5f' % \
              self.panel_y_axis[panel_id], file=out)

        print('SEGMENT_DISTANCE= %.3f' % self.panel_distance[panel_id], file=out)

        print('SEGMENT_ORGX= %.2f SEGMENT_ORGY= %.2f' % self.panel_origin[panel_id], file=out)
        print(file=out)

    if as_str:
      return '\n'.join(str_result)

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
