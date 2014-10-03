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
  def __init__(self, pickle_name, json_name):
    # assemble template and unpack files
    template_pickle = {'mapped_predictions':[[]],
                       'distance':0,
                       'ybeam':0,
                       'current_orientation':0,
                       'ML_half_mosaicity_deg':[0],
                       'current_cb_op_to_primitive': 0,
                       'effective_tiling':0,
                       'residual':0,
                       'sa_parameters':['None'],
                       'model_partialities':[None],
                       'ML_domain_size_ang':[0],
                       'mosaicity':0,
                       'observations':[0],
                       'wavelength':0,
                       'xbeam':0,
                       'pointgroup':0,
                       'max_signal':0,
                       'correction_vectors':[[0]]}

    self.frame = template_pickle

    importer = Importer([pickle_name, json_name], read_experiments=True, read_reflections=True)

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
    self.frame['observations'][0] = cctbx.miller.array(miller_set, intensities, variances)

  # convert a flex.vec3_double object with key 'xyzobs.px.value' to a flex.vec2_double object with key 'mapped_predictions'
  def populate_pixel_positions(self):
    assert self.reflections.has_key('xyzobs.px.value')
    self.frame['mapped_predictions'][0] = flex.vec2_double()
    for i in range(len(self.reflections['xyzobs.px.value'])):
      self.frame['mapped_predictions'][0].append((self.reflections['xyzobs.px.value'][i][0], self.reflections['xyzobs.px.value'][i][1]))


  # fetch the point group associated with the crystal
  def populate_point_group(self):
    assert self.xtal.get_space_group() is not None
    self.frame['pointgroup'] = str(self.xtal.get_space_group().build_derived_point_group().info())

  # get wavelength
  def populate_wavelength(self):
    assert self.beam_obj.get_wavelength() is not None
    self.frame['wavelength'] = self.beam_obj.get_wavelength()

  # get mosaicity
  def populate_mosaicity(self):
    assert self.xtal.get_mosaicity() is not None
    self.frame['mosaicity'] = self.xtal.get_mosaicity()

  # get distance in mm
  def populate_distance(self):
    assert self.det[0].get_distance() is not None
    self.frame['distance'] = self.det[0].get_distance()

  # get xbeam and ybeam in mm
  def populate_beam_dir(self):
    assert self.beam_obj.get_s0() is not None
    self.frame['xbeam'], self.frame['ybeam'] = self.det[0].get_beam_centre(self.beam_obj.get_s0())

  # generate a crystal orientation object from the A matrix
  def populate_orientation(self):
    assert self.xtal.get_A() is not None
    self.frame['current_orientation'] = [crystal_orientation(self.xtal.get_A().elems, True)] # if buggy, check here first

  # generate change-of-basis operation for current to primitive cell
  def populate_op_to_primitive(self):
    assert self.xtal.get_space_group() is not None
    self.frame['current_cb_op_to_primitive'] = [self.xtal.get_space_group().z2p_op()]

  # combine all of the above
  def make_frame(self):
    self.populate_observations()
    self.populate_pixel_positions()
    self.populate_point_group()
    self.populate_wavelength()
    self.populate_mosaicity()
    self.populate_distance()
    self.populate_beam_dir()
    self.populate_orientation()
    self.populate_op_to_primitive()
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
