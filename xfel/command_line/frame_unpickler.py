from __future__ import absolute_import, division, print_function
from six.moves import range
# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# LIBTBX_SET_DISPATCHER_NAME frame.unpickler
#
# $Id: frame_unpickler.py idyoung $

from dials.array_family import flex
from scitbx.array_family import flex as sciflex
from libtbx import easy_pickle
from libtbx.utils import Sorry
from dials.util.options import ArgumentParser
from dxtbx.model import BeamFactory, DetectorFactory
from dxtbx.format import FormatMultiImage
from dxtbx.model import Experiment, ExperimentList
from dxtbx import imageset
import dxtbx
import iotbx.phil
import os

def int_pickle_to_filename(int_name, prefix, suffix):
  parts = int_name.split("-")
  parts = parts[0:3] + parts[3].split(":")
  parts = parts[0:4] + parts[4].split(".")
  parts = parts[0:5] + parts[5].split("_")
  outname = prefix + parts[1] + parts[2] \
    + parts[3][0:2] + parts[3][3:5] + parts[4][0:2] + parts[4][3:5] + parts[5] + suffix
  return outname

def find_matching_img(pickle, img_location=None):
  # find the image associated with the pickle file to be processed, or use the pickle file itself as a placeholder
  if pickle is None:
    raise Sorry("Cannot find matching image for NoneType in place of pickle")
  if img_location is None:
    if os.path.exists(pickle.split(".pickle")[0] + ".cbf"):
      return pickle.split(".pickle")[0] + ".cbf"
    elif os.path.exists(os.path.join(os.path.dirname(pickle), "..", "out")):
      loc = os.path.join(os.path.dirname(pickle), "..", "out")
      name = os.path.basename(pickle).split(".pickle")[0]
      imgname = int_pickle_to_filename(name, "idx-", ".pickle")
      if os.path.exists(os.path.join(loc, imgname)):
        return os.path.join(loc, imgname)
      else:
        return None
    else:
      return None
  else:
    loc = img_location
    name = os.path.basename(pickle).split(".pickle")[0]
    imgname = int_pickle_to_filename(name, "idx-", ".pickle")
    if os.path.exists(os.path.join(loc, imgname)):
      return os.path.join(loc, imgname)
    else:
      return None

class construct_reflection_table_and_experiment_list(object):
  def __init__(self, pickle, img_location, pixel_size, proceed_without_image=False):
    # unpickle pickle file and keep track of image location
    if img_location is None:
      if proceed_without_image:
        self.img_location = []
      else:
        raise Sorry("No image found at specified location. Override by setting proceed_without_image to False"
        + "to produce experiment lists that may only be read when check_format is False.")
    else:
      self.img_location = [img_location]

    # pickle can be already unpickled, if so, loading it will fail with an AttributeError. A load
    # error will fail with EOFError
    try:
      self.data = easy_pickle.load(pickle)
    except EOFError:
      self.data = None
      self.pickle = None
    except AttributeError:
      self.data = pickle
      self.pickle = None
    else:
      self.pickle = pickle
    if self.data is not None:
      self.length = len(self.data['observations'][0].data())
      self.pixel_size = pixel_size

  # extract things from pickle file
  def unpack_pickle(self):
    """Extract all relevant information from an integration pickle file."""

    # crystal-dependent
    self.ori = self.data['current_orientation'][0]
    self.ucell = self.data['current_orientation'][0].unit_cell()

    # experiment-dependent
    self.wavelength = self.data['wavelength']
    self.det_dist = self.data['distance']

    # observation-dependent
    self.observations = self.data['observations'][0]
    self.predictions = self.data['mapped_predictions'][0]

    if 'fuller_kapton_absorption_correction' in self.data:
      self.fuller_correction = self.data['fuller_kapton_absorption_correction'][0]
      if 'fuller_kapton_absorption_correction_sigmas' in self.data:
        self.fuller_correction_sigmas = self.data['fuller_kapton_absorption_correction_sigmas'][0]

  # construct the experiments and experiment_list objects
  def expt_beam_maker(self):
    """Construct the beam object for the experiments file."""
    self.beam = BeamFactory.simple(self.wavelength)

  def expt_crystal_maker(self):
    """Construct the crystal object for the experiments file."""
    a, b, c = self.ucell.parameters()[0:3]
    direct_matrix = self.ori.direct_matrix()
    real_a = direct_matrix[0:3]
    real_b = direct_matrix[3:6]
    real_c = direct_matrix[6:9]
    lattice = self.ucell.lattice_symmetry_group()
    found_it = False
    if 'ML_half_mosaicity_deg' in self.data:
      assert 'ML_domain_size_ang' in self.data
      if d['ML_half_mosaicity_deg'][0] is None or d['ML_domain_size_ang'][0] is None:
        assert d['ML_half_mosaicity_deg'][0] is None and d['ML_domain_size_ang'][0] is None
      else:
        found_it = True
        if 'mosaicity' in self.data and self.data['mosaicity'] > 0:
          print("Warning, two kinds of mosaicity found. Using Sauter2014 model")
        from dxtbx.model import MosaicCrystalSauter2014
        self.crystal = MosaicCrystalSauter2014(real_a, real_b, real_c, space_group=lattice)
        self.crystal.set_half_mosaicity_deg(self.data['ML_half_mosaicity_deg'][0])
        self.crystal.set_domain_size_ang(self.data['ML_domain_size_ang'][0])
    if not found_it:
      if 'mosaicity' in self.data:
        from dxtbx.model import MosaicCrystalKabsch2010
        self.crystal = MosaicCrystalKabsch2010(real_a, real_b, real_c, space_group=lattice)
        self.crystal.set_mosaicity(self.data['mosaicity'])
      else:
        from dxtbx.model import Crystal
        self.crystal = Crystal(real_a, real_b, real_c, space_group=lattice)
    if 'identified_isoform' in self.data and self.data['identified_isoform'] is not None:
      self.crystal.identified_isoform = self.data['identified_isoform']

  def expt_detector_maker(self):
    """Construct the detector object for the experiments file. This function generates a monolithic flattening of the
    CSPAD detector if not supplied with an image file."""
    self.distance = self.data['distance']
    self.xbeam, self.ybeam = self.data['xbeam'], self.data['ybeam']
    if len(self.img_location) > 0 and not dxtbx.load(self.img_location[0])._image_file.endswith("_00000.pickle"):
      self.detector = dxtbx.load(self.img_location[0])._detector()
    else:
      self.detector = DetectorFactory.simple('SENSOR_UNKNOWN',self.distance,(self.xbeam, self.ybeam),'+x','-y',
      (self.pixel_size, self.pixel_size),(1765,1765))

  def expt_gonio_maker(self):
    """XFEL data consisting of stills is expected to have been generated by an experiment without a goniometer -- use placeholder None."""
    self.goniometer = None

  def expt_imageset_maker(self):
    """Construct the imageset object for the experiments file."""
    if len(self.img_location) == 0:
      self.imageset = None
      return
    self.filename = self.img_location
    self.format = FormatMultiImage.FormatMultiImage()
    self.reader = imageset.MultiFileReader(self.format, self.filename)
    self.imageset = imageset.ImageSet(self.reader)

  def expt_scan_maker(self):
    """XFEL data consisting of stills is expected not to contain scans -- use placeholder None."""
    self.scan = None

  def assemble_experiments(self):
    self.unpack_pickle()
    self.expt_beam_maker()
    self.expt_crystal_maker()
    self.expt_detector_maker()
    self.expt_gonio_maker()
    self.expt_imageset_maker()
    self.expt_scan_maker()
    self.experiment = Experiment(beam = self.beam,
                                                 crystal=self.crystal,
                                                 detector=self.detector,
                                                 goniometer=self.goniometer,
                                                 imageset=self.imageset,
                                                 scan=self.scan)
    self.experiment_list = ExperimentList([self.experiment])

  def experiments_to_json(self, path_name=None):
    if path_name is None:
      loc = os.path.dirname(self.pickle)
    else:
      loc = path_name
    name = os.path.basename(self.pickle).split(".pickle")[0]
    expt_name = int_pickle_to_filename(name, "idx-", ".expt")
    experiments = os.path.join(loc, expt_name)
    self.experiment_list.as_file(experiments)

  # construct the reflection table
  def refl_table_maker(self):
    self.reflections = flex.reflection_table()

  def refl_bkgd_maker(self):
    self.reflections['background.mean'] = sciflex.double(self.length)
    self.reflections['background.mse'] = sciflex.double(self.length)

  def refl_bbox_maker(self):
    self.reflections['bbox'] = flex.int6(self.length)

  def refl_correlation_maker(self):
    self.reflections['correlation.ideal.profile'] = sciflex.double(self.length)

  def refl_entering_maker(self):
    self.reflections['entering'] = flex.bool(self.length)

  def refl_flags_maker(self):
    self.reflections['flags'] = flex.size_t(self.length, 1)

  def refl_id_maker(self):
    self.reflections['id'] = sciflex.size_t(self.length)

  def refl_intensities_maker(self):
    self.reflections['intensity.sum.value'] = self.observations.data()
    self.reflections['intensity.sum.variance'] = self.observations.sigmas()**2

  def refl_lp_maker(self):
    self.reflections['lp'] = sciflex.double(self.length)

  def refl_millers_maker(self):
    self.reflections['miller_index'] = self.observations._indices

  def refl_nbackgroundforeground_maker(self):
    self.reflections['n_background'] = sciflex.size_t(self.length)
    self.reflections['n_foreground'] = sciflex.size_t(self.length)

  def refl_panel_maker(self):
    self.reflections['panel'] = sciflex.size_t(self.length, 0)

  def refl_profile_maker(self):
    self.reflections['profile.correlation'] = sciflex.double(self.length)

  def refl_s1_maker(self):
    from scitbx.matrix import col
    self.reflections['s1'] = sciflex.vec3_double(self.length)
    for idx in range(self.length):
      coords = col(self.detector[0].get_pixel_lab_coord(self.reflections['xyzobs.px.value'][idx][0:2])).normalize()
      self.reflections['s1'][idx] = tuple(coords)

  def refl_xyzcal_maker(self):
    self.reflections['xyzcal.px'] = sciflex.vec3_double(self.predictions.parts()[1], self.predictions.parts()[0], sciflex.double(self.length, 0.0))
    self.reflections['xyzcal.mm'] = self.pixel_size * self.reflections['xyzcal.px']

  def refl_xyzobs_maker(self):
    self.reflections['xyzobs.px.value'] = sciflex.vec3_double(self.predictions.parts()[1], self.predictions.parts()[0], sciflex.double(self.length, 0.5))
    self.reflections['xyzobs.px.variance'] = sciflex.vec3_double(self.length, (0.0,0.0,0.0))

  def refl_zeta_maker(self):
    self.reflections['zeta'] = sciflex.double(self.length)

  def refl_kapton_correction_maker(self):
    if hasattr(self, 'fuller_correction'):
      self.reflections['kapton_absorption_correction'] = self.fuller_correction
      if hasattr(self, 'fuller_correction_sigmas'):
        self.reflections['kapton_absorption_correction_sigmas'] = self.fuller_correction_sigmas

  def assemble_reflections(self):
    self.refl_table_maker()
    self.refl_bkgd_maker()
    self.refl_bbox_maker()
    self.refl_correlation_maker()
    self.refl_entering_maker()
    self.refl_flags_maker()
    self.refl_id_maker()
    self.refl_intensities_maker()
    self.refl_millers_maker()
    self.refl_nbackgroundforeground_maker()
    self.refl_panel_maker()
    self.refl_xyzcal_maker()
    self.refl_xyzobs_maker()
    self.refl_zeta_maker()
    self.refl_kapton_correction_maker()
    self.refl_s1_maker() # depends on successful completion of refl_xyz_obs_maker

  def reflections_to_pickle(self, path_name=None):
    if path_name is None:
      loc = os.path.dirname(self.pickle)
    else:
      loc = path_name
    name = os.path.basename(self.pickle).split(".refl")[0]
    refl_name = int_pickle_to_filename(name, "idx-", "_integrated.refl")
    reflections = os.path.join(loc, refl_name)
    self.reflections.as_pickle(reflections)

if __name__ == "__main__":
  master_phil_scope = iotbx.phil.parse("""
    pickle_location = None
      .type = path
      .help = supply complete path to integration pickle to be unpacked.
    img_location = None
      .type = path
      .help = supply complete path to image as either cbf or image pickle. If None, will attempt to locate cbfs
      .help = in the same directory or the sibling out directory
    json_location = None
      .type = path
      .help = supply complete path to new json directory. If None, will place all generated files in the
      .help = same directory as the integration pickles.
    refl_location = None
      .type = path
      .help = supply complete path to new reflection table directory. If None, will place all generated
      .help = files in the same directory as the integration pickles.
    pixel_size = 0.11
      .type = float
      .help = detector-specific parameter for pixel size in mm
    proceed_without_image = False
      .type = bool
      .help = override raising exceptions if no image can be found. Write a json with a dummy image (pickle name)
    """)
  parser = ArgumentParser(phil=master_phil_scope)
  params, options = parser.parse_args(show_diff_phil=False)
  for pickle_file in os.listdir(params.pickle_location):
    pickle_path = os.path.join(params.pickle_location, pickle_file)
    result = construct_reflection_table_and_experiment_list(pickle_path, find_matching_img(pickle_path, params.img_location), params.pixel_size, proceed_without_image=params.proceed_without_image)
    if result.data is not None:
      result.assemble_experiments()
      result.assemble_reflections()
      result.experiments_to_json(params.json_location)
      result.reflections_to_pickle(params.refl_location)
    else:
      print("Skipping unreadable pickle file at", pickle_path)
  print('Generated experiments.expt and integrated.refl files.')
