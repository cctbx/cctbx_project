"""Computes F_TAAM - F_IAM difference map based on a model."""
from __future__ import absolute_import, division, print_function
try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  from libtbx.program_template import ProgramTemplate
import os
from libtbx.utils import Sorry
from iotbx import crystal_symmetry_from_any
from scitbx.array_family import flex
from cctbx import miller
# Check if pydiscamb is installed
try:
  import pydiscamb
except ImportError:
  raise Sorry('PyDiscamb not installed.')
from pydiscamb import DiscambWrapper

master_phil_str = '''
include scope libtbx.phil.interface.tracking_params
resolution = None
  .type = float
scattering_table = *wk1995 it1992 electron
  .type = choice
  .short_caption = Scattering table
  .help = Scattering table for structure factors calculations
apply_scaling = True
  .type = bool
  .short_caption = Apply a scale factor for Fcalc(TAAM)
  .help = Apply a scale factor for Fcalc(TAAM)
'''

# =============================================================================

class Program(ProgramTemplate):

  description = '''

Computes F_TAAM - F_IAM map based on a model.

Inputs:
  - Model file
  and
  - Resolution limit
  or
  - Reflection file

Usage examples:
  1. phenix.TAAM_minus_IAM model.cif data.mtz
  2. phenix.TAAM_minus_IAM model.pdb data.mtz
  3. phenix.TAAM_minus_IAM model.pdb resolution=1.0

'''

  datatypes = ['model', 'phil', 'miller_array']

  master_phil_str = master_phil_str
  use_scattering_table_for_default_type = 'scattering_table'

  # ---------------------------------------------------------------------------

  def validate(self):
    print('Validating inputs:\n', file=self.logger)
    self.data_manager.has_models(
      raise_sorry = True,
      expected_n  = 1,
      exact_count = True)
    has_miller_arrays = self.data_manager.has_miller_arrays()
    if not has_miller_arrays and self.params.resolution is None:
      raise Sorry('Supply either resolution keyword or a reflection file')
    if has_miller_arrays and self.params.resolution is not None:
      raise Sorry('Resolution or a reflection file are needed.')
    if (len(self.data_manager.get_miller_array_names()) > 1):
      raise Sorry('Dont input more than one (1) reflection file.')
    if has_miller_arrays:
      self.fmodel = None
      try:
        fmodel_params = self.data_manager.get_fmodel_params()
        fmodel_params.xray_data.r_free_flags.required = False
        fmodel_params.xray_data.r_free_flags.ignore_r_free_flags = True
        self.data_manager.set_fmodel_params(fmodel_params)
        self.fmodel = self.data_manager.get_fmodel(
          scattering_table = self.params.scattering_table)
      #except Sorry as s:
      #  if 'previously used R-free flags are available run this command again' in str(s):
      except Exception as e:
        print(e)
      if self.fmodel is None:
        raise Sorry('Failed to create fmodel. Please submit a bug report.')

  # ---------------------------------------------------------------------------

  def check_crystal_symmetry(self):
    crystal_symmetries = []
    files = self.data_manager.get_model_names() + \
      self.data_manager.get_miller_array_names()
    for f in files:
      cs = crystal_symmetry_from_any.extract_from(f)
      if (cs is not None):
        crystal_symmetries.append(cs)
    if (len(crystal_symmetries) == 0): raise Sorry("No crystal symmetry found.")
    elif (len(crystal_symmetries) > 1):
      if (not crystal_symmetries[0].is_similar_symmetry(crystal_symmetries[1])):
        raise Sorry("Crystal symmetry mismatch between model and reflection file.")
      crystal_symmetry = crystal_symmetries[0]

  # ---------------------------------------------------------------------------

  def scale(self, x, y):
    x = flex.abs(x)
    y = flex.abs(y)
    return flex.sum(x*y)/flex.sum(y*y)

  # ---------------------------------------------------------------------------

  def r(self, x, y):
    x = flex.abs(x)
    y = flex.abs(y)
    scale = flex.sum(x*y)/flex.sum(y*y)
    num = flex.sum(flex.abs(x-scale*y))
    den = flex.sum(flex.abs(x+scale*y))
    return num/den*2*100.

  # ---------------------------------------------------------------------------

  def cctbx_direct(self, xrs, complete_set):
    return complete_set.structure_factors_from_scatterers(
      xray_structure = xrs,
      algorithm      = "direct",
      cos_sin_table  =  False).f_calc()

  # ---------------------------------------------------------------------------

  def discamb_taam(self, xrs, complete_set):
    wrapper = pydiscamb.DiscambWrapper(
      xrs,
      method=pydiscamb.FCalcMethod.TAAM,
      assignment_info="atom_type_assignment.log") # XXX Can't be None
    wrapper.show_atom_type_assignment(log=self.logger)
    wrapper.set_indices(complete_set.indices())
    return wrapper.f_calc()

  # ---------------------------------------------------------------------------

  def run(self):

    print('Using model file:', self.data_manager.get_default_model_name(),
      file=self.logger)
    print('\nUsing scattering table:', self.params.scattering_table)

    model = self.data_manager.get_model()
    xrs = model.get_xray_structure()
    xrs.scattering_type_registry(table = self.params.scattering_table)

    if self.params.resolution is not None:
      miller_array = miller.build_set(
        crystal_symmetry = xrs.crystal_symmetry(),
        anomalous_flag   = False,
        d_min            = self.params.resolution)
    else:
      miller_array = None
      print('Using reflection file:',
        self.data_manager.get_default_miller_array_name(), file=self.logger)
      self.check_crystal_symmetry()
      user_selected_labels = self.data_manager.get_miller_array_user_selected_labels()
      if not user_selected_labels: user_selected_labels = None
      miller_arrays = self.data_manager.get_miller_arrays(
        labels=user_selected_labels)
      data_sizes = flex.int([ma.data().size() for ma in miller_arrays])
      if(data_sizes.all_eq(data_sizes[0])): miller_array = miller_arrays[0]
      else:
        raise Sorry('Reflection file contains arrays of different lengths. \
                    Please select one using "labels.name=" keyword.')
      miller_array = miller_arrays[0]
      assert(miller_array is not None)
      miller_array.show_comprehensive_summary(f = self.logger, prefix="  ")
      miller_array = miller_array.map_to_asu().customized_copy(
        data = flex.double(miller_array.data().size(), 1))

    fc1 = self.cctbx_direct(xrs=xrs, complete_set=miller_array)
    fc2 = self.discamb_taam(xrs=xrs, complete_set=miller_array)
    if self.params.apply_scaling:
      sc = self.scale(x=fc1.data(), y=fc2)
      print("\nscale:", sc, file=self.logger)
      fc2 = fc1.customized_copy(data = fc2*sc)
      print("\nr:",self.r(fc1.data(), fc2.data()), file=self.logger)
    else:
      fc2 = fc1.customized_copy(data = fc2)
    diff = fc1.customized_copy(data = fc2.data()-fc1.data())

    mtz_dataset = diff.as_mtz_dataset(column_root_label = "TAAM-IAM")
    mtz_object = mtz_dataset.mtz_object()

    model_basename = os.path.basename(
      self.data_manager.get_default_model_name().split(".")[0])
    fn = model_basename + "_TAAM_minus_IAM.mtz"
    mtz_object.write(file_name = fn)
    print('\nFile %s was written.' % fn, file=self.logger)
