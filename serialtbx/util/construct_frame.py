from __future__ import absolute_import, division, print_function
from six.moves import range

from scitbx.array_family import flex
from cctbx import crystal, miller
from cctbx.crystal_orientation import crystal_orientation
import cctbx

class ConstructFrame(object):
  def get_template_pickle(self):
    return {
                       'beam_s0': 0,
                       'beam_polarization_normal': 0,
                       'current_cb_op_to_primitive': 0,
                       'current_orientation':0,
                       'distance':0,
                       'effective_tiling':0,
                       'identified_isoform':None,
                       'mapped_predictions':[[]],
                       'max_signal':0,
                       'ML_domain_size_ang':[0],
                       'ML_half_mosaicity_deg':[0],
                       'mosaicity':0,
                       'model_partialities':[None],
                       'observations':[0],
                       'pointgroup':0,
                       'residual':0,
                       'sa_parameters':['None'],
                       's1_vec':[0],
                       'wavelength':0,
                       'xbeam':0,
                       'ybeam':0}

  def __init__(self, reflections, experiment):
    # assemble template and unpack files
    self.frame = self.get_template_pickle()
    self.pixel_size = experiment.detector[0].get_pixel_size()[0]

    if 'intensity.prf.value' in reflections:
      self.method = 'prf' # integration by profile fitting
    elif 'intensity.sum.value' in reflections:
      self.method = 'sum' # integration by simple summation
    self.reflections = reflections.select(reflections['intensity.' + self.method + '.variance'] > 0) # keep only spots with sigmas above zero

    self.xtal = experiment.crystal
    self.beam_obj = experiment.beam
    self.det = experiment.detector
    self.gonio = experiment.goniometer
    self.scan = experiment.scan
    self.img_sweep = experiment.imageset

  # experiment-dependent components ---------------------------------------------------------------------------

  # get wavelength
  def populate_wavelength(self):
    assert self.beam_obj.get_wavelength() is not None, "no wavelength"
    self.frame['wavelength'] = self.beam_obj.get_wavelength()
    self.frame['beam_s0'] = self.beam_obj.get_s0()
    self.frame['beam_polarization_normal'] = self.beam_obj.get_polarization_normal()

  # get detector distance in mm
  def populate_distance(self):
    assert self.det[0].get_distance() is not None, "no detector distance"
    self.frame['distance'] = self.det[0].get_distance()

  # get xbeam and ybeam in mm
  def populate_beam_dir(self):
    assert self.beam_obj.get_s0() is not None, "no beam direction"
    self.frame['xbeam'], self.frame['ybeam'] = self.det[0].get_beam_centre(self.beam_obj.get_s0())

  # get isoform
  def populate_isoform(self):
    if hasattr(self.xtal, "identified_isoform"):
      self.frame['identified_isoform'] = self.xtal.identified_isoform

  # get max signal
  def populate_max_signal(self):
    pass

  # get effective tiling
  def populate_effective_tiling(self):
    pass

  # indicate simulated annealing parameters, if present
  def populate_sa_params(self):
    pass

  # crystal-dependent components ------------------------------------------------------------------------------

  # generate a crystal orientation object from the A* matrix
  def populate_orientation(self):
    assert self.xtal.get_A() is not None, "no crystal orientation matrix"
    self.frame['current_orientation'] = [crystal_orientation(self.xtal.get_A(), True)]

  # generate change-of-basis operation for current to primitive cell
  def populate_op_to_primitive(self):
    assert self.xtal.get_space_group() is not None, "no space group"
    self.frame['current_cb_op_to_primitive'] = [self.xtal.get_space_group().z2p_op()]

  # fetch the point group associated with the crystal
  def populate_point_group(self):
    assert self.xtal.get_space_group() is not None, "no space group"
    self.frame['pointgroup'] = str(self.xtal.get_space_group().build_derived_point_group().info())

  # get mosaicity
  def populate_mosaicity(self):
    try:
      assert self.xtal.get_mosaicity() is not None, "no mosaicity"
      self.frame['mosaicity'] = self.xtal.get_mosaicity()
    except Exception:
      pass

  # get any available ML values
  def populate_ML_values(self):
    try:
      self.frame['ML_half_mosaicity_deg'] = [self.xtal.get_half_mosaicity_deg()]
    except AttributeError:
      pass
    try:
      self.frame['ML_domain_size_ang'] = [self.xtal.get_domain_size_ang()]
    except AttributeError:
      pass

  # observations-dependent components -------------------------------------------------------------------------

  # generate a miller array containing the Miller indices, intensities and variances for one frame
  def populate_observations(self):
    intensities = self.reflections['intensity.' + self.method + '.value']
    variances = self.reflections['intensity.' + self.method + '.variance']
    space_group = crystal.symmetry(self.xtal.get_unit_cell(), str(self.xtal.get_space_group().info()))
    miller_set = miller.set(space_group, self.reflections['miller_index'])
    self.frame['observations'][0] = cctbx.miller.array(miller_set, intensities, flex.sqrt(variances)).set_observation_type_xray_intensity()
    self.frame['s1_vec'][0] = self.reflections['s1']

  # collect predicted spot positions
  def populate_pixel_positions(self):
    assert 'xyzcal.px' in self.reflections, "no calculated spot positions"
    self.frame['mapped_predictions'][0] = flex.vec2_double()
    for i in range(len(self.reflections['xyzcal.px'])):
      self.frame['mapped_predictions'][0].append(tuple(self.reflections['xyzcal.px'][i][1::-1])) # 1::-1 reverses the order taking only the first two elements first.

  # generate a list of dictionaries containing a series of corrections for each predicted reflection
  def populate_corrections(self):
    assert 'xyzobs.px.value' in self.reflections and 'xyzcal.px' in self.reflections, "no calculated or observed spot positions"
    assert self.frame['xbeam'] != 0 and self.frame['ybeam'] != 0, "invalid beam center"
    self.frame['correction_vectors'] = [[]]
    for idx in range(len(self.reflections['xyzobs.px.value'])):
      if self.reflections['xyzcal.px'][idx][0:2] != self.reflections['xyzobs.px.value'][idx][0:2]:
        theoret_center = 1765/2, 1765/2
        refined_center = self.frame['xbeam']/self.pixel_size, self.frame['ybeam']/self.pixel_size # px to mm conversion
        hkl = self.reflections['miller_index'][idx]
        obsspot = tuple(self.reflections['xyzobs.px.value'][idx][0:2])
        predspot = tuple(self.reflections['xyzcal.px'][idx][0:2])
        self.frame['correction_vectors'][0].append({'refinedcenter':refined_center, 'hkl':hkl, 'setting_id':0, 'azimuthal':0, 'radial':0,
          'obsspot':obsspot, 'obscenter':theoret_center, 'predspot':predspot})

  # get partialities
  def populate_partialities(self):
    pass

  # produce residuals
  def populate_residuals(self):
    pass

  # collect kapton corrections if present
  def populate_kapton_corrections(self):
    if 'kapton_absorption_correction' in self.reflections:
      self.frame['fuller_kapton_absorption_correction'] = [self.reflections['kapton_absorption_correction']]
      if 'kapton_absorption_correction_sigmas' in self.reflections:
        self.frame['fuller_kapton_absorption_correction_sigmas'] = [self.reflections['kapton_absorption_correction_sigmas']]

  # combine all of the above
  def make_frame(self):
    self.populate_wavelength()
    self.populate_distance()
    self.populate_beam_dir()
    self.populate_isoform()
    self.populate_max_signal()
    self.populate_effective_tiling()
    self.populate_sa_params()
    self.populate_orientation()
    self.populate_op_to_primitive()
    self.populate_point_group()
    self.populate_mosaicity()
    self.populate_ML_values()
    self.populate_observations()
    self.populate_pixel_positions()
    # self.populate_corrections() # works, but unnecessary
    self.populate_partialities()
    self.populate_residuals()
    self.populate_kapton_corrections()
    return self.frame
