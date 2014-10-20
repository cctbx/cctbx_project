from __future__ import division
# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# LIBTBX_SET_DISPATCHER_NAME cxi.frame_extractor
# LIBTBX_SET_DISPATCHER_NAME xfel.frame_extractor
#
# $Id: frame_extractor.py idyoung $

from dials.array_family import flex
from dials.util.options import Importer, flatten_reflections, flatten_experiments, OptionParser
from cctbx import crystal, miller
from cctbx.crystal_orientation import crystal_orientation
import iotbx.phil
import cctbx

class construct_frame(object):
  def __init__(self, pickle_name, json_name, pixel_size):
    # assemble template and unpack files
    template_pickle = {'correction_vectors':[[]],
                       'current_cb_op_to_primitive': 0,
                       'current_orientation':0,
                       'distance':0,
                       'effective_tiling':0,
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
                       'wavelength':0,
                       'xbeam':0,
                       'ybeam':0}

    self.frame = template_pickle
    self.pixel_size = pixel_size

    importer = Importer([pickle_name, json_name], read_experiments=True, read_reflections=True)
    if importer.unhandled:
      print "unable to process:", importer.unhandled

    # load the integration.pickle file (reflection table) into memory
    self.reflections = flatten_reflections(importer.reflections)[0]

    # load the experiments.json file (json) into memory, piecewise
    experiments = flatten_experiments(importer.experiments)[0]

    self.xtal = experiments.crystal
    self.beam_obj = experiments.beam
    self.det = experiments.detector
    self.gonio = experiments.goniometer
    self.scan = experiments.scan
    self.img_sweep = experiments.imageset

  # experiment-dependent components ---------------------------------------------------------------------------

  # get wavelength
  def populate_wavelength(self):
    assert self.beam_obj.get_wavelength() is not None, "no wavelength"
    self.frame['wavelength'] = self.beam_obj.get_wavelength()

  # get detector distance in mm
  def populate_distance(self):
    assert self.det[0].get_distance() is not None, "no detector distance"
    self.frame['distance'] = self.det[0].get_distance()

  # get xbeam and ybeam in mm
  def populate_beam_dir(self):
    assert self.beam_obj.get_s0() is not None, "no beam direction"
    self.frame['xbeam'], self.frame['ybeam'] = self.det[0].get_beam_centre(self.beam_obj.get_s0())

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
    self.frame['current_orientation'] = [crystal_orientation(self.xtal.get_A().elems, True)]

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
    assert self.xtal.get_mosaicity() is not None, "no mosaicity"
    self.frame['mosaicity'] = self.xtal.get_mosaicity()

  # get any available ML values
  def populate_ML_values(self):
    try:
      self.frame['ML_half_mosaicity_deg'] = [self.xtal._ML_half_mosaicity_deg]
    except AttributeError:
      pass
    try:
      self.frame['ML_domain_size_ang'] = [self.xtal._ML_domain_size_ang]
    except AttributeError:
      pass

  # observations-dependent components -------------------------------------------------------------------------

  # generate a miller array containing the Miller indices, intensities and variances for one frame
  def populate_observations(self):
    if self.reflections.has_key('intensity.prf.value'):
      intensities = self.reflections['intensity.prf.value']
      variances = self.reflections['intensity.prf.variance']
    elif self.reflections.has_key('intensity.sum.value'):
      intensities = self.reflections['intensity.sum.value']
      variances = self.reflections['intensity.sum.variance']
    space_group = crystal.symmetry(self.xtal.get_unit_cell(), str(self.xtal.get_space_group().info()))
    miller_set = miller.set(space_group, self.reflections['miller_index'])
    self.frame['observations'][0] = cctbx.miller.array(miller_set, intensities, variances).set_observation_type_xray_intensity()

  # collect predicted spot positions
  def populate_pixel_positions(self):
    assert self.reflections.has_key('xyzcal.px'), "no calculated spot positions"
    self.frame['mapped_predictions'][0] = flex.vec2_double()
    for i in range(len(self.reflections['xyzcal.px'])):
      self.frame['mapped_predictions'][0].append(tuple(self.reflections['xyzcal.px'][i][0:2]))

  # generate a list of dictionaries containing a series of corrections for each predicted reflection
  def populate_corrections(self):
    assert self.reflections.has_key('xyzobs.px.value') and self.reflections.has_key('xyzcal.px'), "no calculated or observed spot positions"
    assert self.frame['xbeam'] is not 0 and self.frame['ybeam'] is not 0, "invalid beam center"
    for idx in xrange(len(self.reflections['xyzobs.px.value'])):
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

  # combine all of the above
  def make_frame(self):
    self.populate_wavelength()
    self.populate_distance()
    self.populate_beam_dir()
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
    self.populate_corrections()
    self.populate_partialities()
    self.populate_residuals()
    return self.frame

if __name__ == "__main__":
  master_phil_scope = iotbx.phil.parse("""
    pickle_name = None
      .type = path
      .help = path to a reflection table (integrated.pickle) file
    json_name = None
      .type = path
      .help = path to an experiments.json file
    pixel_size = 0.11
      .type = float
      .help = detector-specific parameter for pixel size in mm
      """)
  parser = OptionParser(phil=master_phil_scope)
  params, options = parser.parse_args(show_diff_phil=True)
  frame = construct_frame(params.pickle_name, params.json_name).make_frame()
