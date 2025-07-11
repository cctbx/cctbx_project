"""Calculate model structure factors"""
from __future__ import absolute_import, division, print_function
try:
  from phenix.program_template import ProgramTemplate
except ImportError:
  from libtbx.program_template import ProgramTemplate
import os
import mmtbx.utils
import iotbx.phil
import iotbx.pdb
import random
from libtbx.utils import Sorry
from cctbx.array_family import flex
from libtbx import group_args

fmodel_from_xray_structure_params_str = """\
fmodel
  .short_caption = F(model) options
  .expert_level = 1
  .style = auto_align box
{
  k_sol = 0.0
    .type = float
    .help = Bulk solvent k_sol values
    .short_caption=Bulk solvent K_sol value
  b_sol = 0.0
    .type = float
    .help = Bulk solvent b_sol values
    .short_caption=Bulk solvent B_sol value
  b_cart = 0 0 0 0 0 0
    .type = floats(6)
    .help = Anisotropic scale matrix
    .input_size = 200
    .short_caption = Anisotropic scale matrix
  scale = 1.0
    .type = float
    .help = Overall scale factor
}
structure_factors_accuracy
  .short_caption = Structure factors accuracy
  .style = auto_align box
{
  include scope mmtbx.f_model.sf_and_grads_accuracy_master_params
}
mask
  .short_caption = Bulk solvent mask
  .style = auto_align box
{
  include scope mmtbx.masks.mask_master_params
}
"""

master_phil_str = '''
include scope libtbx.phil.interface.tracking_params

high_resolution = None
  .type = float
  .expert_level=1
  .style = noauto bold
low_resolution = None
  .type = float
  .expert_level=1
  .style = noauto
r_free_flags_fraction = None
  .type = float
  .expert_level=1
  .style = noauto
add_sigmas = False
  .type = bool
  .expert_level=1
  .help = Adds calculated Sigma(F) column to output file.
  .style = noauto
add_random_error_to_amplitudes_percent = None
  .type = float
  .short_caption = Add random error (percent)
  .style = noauto
scattering_table = wk1995  it1992  *n_gaussian  neutron electron
  .type = choice
  .help = Choices of scattering table for structure factors calculations.  \
    n_gaussian is the standard set of X-ray scattering factors.
  .expert_level=1
  .style = noauto
custom_scattering_factors = None
  .type = path
  .help = Use custom scattering factors and replaces default values entirely
%s
random_seed=None
  .type = int
  .help = Random seed
  .expert_level=2
twin_law = None
  .type = str
  .help = Optional twin law if we want to generate a twinned dataset
  .input_size = 120
  .style = noauto
twin_fraction = None
  .type = float
  .help = Twin fraction, ignored if twin_law is not specified
  .style = noauto
wavelength = None
  .type = float
  .input_size = 80
  .help = Wavelength, sets all atoms to anomalous
  .style = noauto
generate_fake_p1_symmetry = False
  .type = bool
  .short_caption = Generate fake symmetry if necessary
  .help = Allows use of PDB files without CRYST1 records as input.  The \
    crystal symmetry will be assumed to be a P1 box.
output
  .short_caption = Reflection output
  .expert_level=0
  .style = noauto
{
  format = *mtz cns
    .type = choice
    .short_caption = File format
    .input_size = 100
  label = FMODEL
    .type = str
    .short_caption = Data label
    .input_size = 100
  type = real *complex
    .type = choice
    .short_caption = Output data type
    .help = Numeric type of output data.  'real' is amplitudes only, \
      'complex' is complete structure factors as complex numbers.
    .expert_level=1
    .style = bold
  obs_type = *amplitudes intensities
    .type = choice
    .help = Experimental observation type to output.  Certain restrictions \
      apply if intensities are selected.
    .expert_level = 2
  file_name = None
    .type = path
    .short_caption = Output file
    .style = bold noauto new_file
  include scope libtbx.phil.interface.tracking_params
}
anomalous_scatterers
  .short_caption = Anomalous sites
  .style = menu_item noauto
{
  group
    .optional = True
    .multiple = True
    .short_caption = Anomalous scatterer group
    .style = auto_align
  {
    selection = None
      .type = atom_selection
      .short_caption = Atom selection
      .input_size = 400
    f_prime = 0
      .type = float
      .short_caption = f'
    f_double_prime = 0
      .type = float
      .short_caption = f''
  }
}
gui
  .help = "GUI-specific parameter required for output directory"
{
  output_dir = None
  .type = path
  .style = output_dir

  data_column_label = None
  .type = str
  .style = noauto renderer:draw_any_label_widget
  .input_size = 300
}
'''%fmodel_from_xray_structure_params_str

fmodel_from_xray_structure_params = iotbx.phil.parse(
  fmodel_from_xray_structure_params_str, process_includes=True)

master_phil = iotbx.phil.parse(master_phil_str, process_includes = True)


def set_fp_fdp_for_anomalous_scatterers(pdb_hierarchy, xray_structure,
  anomalous_scatterer_groups):
  scatterers = xray_structure.scatterers()
  for group in anomalous_scatterer_groups:
    iselection = pdb_hierarchy.atom_selection_cache().selection(
      string = group.selection).iselection()
    if(iselection.size() == 0):
      raise Sorry(
        "Empty selection: selection string '%s' does not select any atom."%
        group.selection)
    for i_seq in iselection:
      scatterers[i_seq].fp = group.f_prime
      scatterers[i_seq].fdp = group.f_double_prime

# =============================================================================

class Program(ProgramTemplate):

  description = '''
phenix.fmodel: a tool to compute structure factors, Fmodel:

  Fmodel = scale * exp(AnisoScale) * (Fcalc + k_sol * exp(-b_sol*s^2/4) * Fmask)

  where:

  - Fmodel - total model structure factor (complex value)
  - AnisoScale = -ht*A(-1)*b_cart*A(-1)th/4
  - h - column vector with Miller indices
  - A - orthogonalization matrix
  - b_cart - anisotropic scale matrix
  - t and (-1) denotes transposition and inversion operations
  - scale - overall scale factor
  - Fcalc - structure factors calculated from atomic model
  - k_sol and b_sol - Flat Bulk solvent model parameters
  - Fmask - structure factors calculated from bulk solvent mask

Usage examples:

  1) phenix.fmodel model.pdb high_resolution=1.5

     will result in a file containing complete set of Fmodel = Fcalc computed
     from atomic model up to 1.5A resolution.

  2) phenix.fmodel model.pdb scale=2 k_sol=0.35 b_sol=50 b_cart="1 2 3 0 4 7" high_res=1.5 low_res=10

     will result in a file containing complete set of Fmodel computed using the
     above formula in resolution range 1.5-20.0A.

  3) phenix.fmodel model.pdb high_resolution=1.5 algorithm=direct

     is similar to "1)" but the Fcalc are computed using direct summation algorithm.

  4) phenix.fmodel model.pdb high_res=1.5 format=cns label=FOBS type=real r_free=0.1

     will result in CNS formatted file containing complete set of amplitudes of
     Fmodel = Fcalc computed up to 1.5A resolution, labelled as FOBS, and free-R
     flags with 10% of test reflections. This is a typical command to simulate Fobs.

  5) phenix.fmodel model.pdb high_res=1.5 scattering_table=neutron

     will result in a file containing complete set of Fmodel = Fcalc computed
     from atomic model up to 1.5A resolution using neutron scattering table.

  6) phenix.fmodel model.pdb parameters.txt

     will result in a structure factor file, where Fmodel were computed using
     parameters defined in parameters.txt file. The parameters.txt file can
     contain all or any subset of parameters listed below. Note, that each {
     must have a matching one }.

  7) phenix.fmodel model.pdb reflection_data.mtz

     will result in a file containing a set of Fmodel = Fcalc that will match
     the set of Miller indices of the data in reflection_data.mtz file.

  8) phenix.fmodel model.pdb reflection_data.mtz data_column_label="FOBS,SIGMA"

     similar to "7)", where the specific data array is selected.

  9) phenix.fmodel model.pdb reflection_data.mtz twin_law="l,-k,h" twin_fraction=0.3

     generates twin data set (real type) with given twin law and fraction.

See below for complete list of available parameters.
'''

  datatypes = ['model', 'phil', 'miller_array']

  master_phil_str = master_phil_str

#  show_data_manager_scope_by_default = True

  # ---------------------------------------------------------------------------

  def validate(self):
    print('Validating inputs:\n', file=self.logger)
    self.data_manager.has_models(
      raise_sorry = True,
      expected_n  = 1,
      exact_count = True)
    if len(self.data_manager.get_miller_array_names()) > 1:
      raise Sorry('Please supply at most one reflection file.')
    if self.data_manager.has_miller_arrays():
      if([self.params.high_resolution, self.params.low_resolution].count(None) != 2):
        raise Sorry("high_resolution and low_resolution must be undefined "+
                    "if reflection data file is given.")
      if len(self.data_manager.get_miller_array_user_selected_labels()) > 1:
        raise Sorry('Supply not more than one label name.')
    else:
      if (self.params.high_resolution is None):
        raise Sorry("Input data file or high_resolution has to be provided.")

    if (self.params.output.type == "complex") and (self.params.add_sigmas):
      raise Sorry("Sigma values only supported when the output type is 'real'.")
    if (self.params.low_resolution is not None
        and self.params.high_resolution is not None):
      if self.params.low_resolution < self.params.high_resolution :
        raise Sorry("Low-resolution cutoff must be larger than the high-"+
          "resolution cutoff.")
    if (self.params.output.obs_type == "intensities"):
      if (self.params.output.type == "complex"):
        raise Sorry("Output type must be 'real' when intensities specified "+
          "for obs_type.")
      if (not self.params.output.label.upper().startswith("I")):
        raise Sorry("Output label must start with 'I' (any case) when "+
          "intensities specified for obs_type (was: %s)." % self.params.output.label)
      if (self.params.output.format != "mtz"):
        raise Sorry("Output format must be 'mtz' when intensities specified.")
    if (self.params.wavelength is not None):
      if (self.params.scattering_table == "neutron"):
        raise Sorry("Wavelength parameter not supported when the neutron "+
          "scattering table is used.")

# For GUI
#def validate_params(params, callback=None):
#  if len(params.pdb_file) == 0 :
#    raise Sorry("You must provide at least one model file to use for "+
#      "F(model) calculations.")
#  if (params.high_resolution is None):
#    if (params.reference_file is None):
#      raise Sorry("Please specify a high-resolution cutoff.")
#  elif (params.reference_file is not None):
#    if (params.data_column_label is None):
#      raise Sorry("Please select a column label to use in the reference "+
#        "data file.")
#    elif ([params.high_resolution, params.low_resolution].count(None) != 2):
#      raise Sorry("High resolution and low resolution must be undefined "+
#                  "if reflection data file is given.")
#  if (params.output.file_name is None):
#    raise Sorry("Please specify an output file.")
#  validate_params_command_line(params)


  # ---------------------------------------------------------------------------

  def run(self):

    print('Parameters to compute Fmodel:', file=self.logger)
    # TODO print processed non defaults from phil

    print('Using model file:', self.data_manager.get_default_model_name(),
      file=self.logger)

    mo = self.data_manager.get_model()
    cs = mo.crystal_symmetry()

    miller_array = None
    if self.data_manager.has_miller_arrays():
      print('Using reflection file:',
        self.data_manager.get_default_miller_array_name(),
        file=self.logger)
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
      cs_from_ma = miller_array.crystal_symmetry()
      if (self.params.generate_fake_p1_symmetry and cs_from_ma is not None):
        raise Sorry("The reflection data already define crystal symmetry; "+
          "you may not use this in combination with the option "+
          "generate_fake_p1_symmetry=True.")
      if cs is None:
        cs = cs_from_ma

    if not self.params.generate_fake_p1_symmetry:
      msg = '''Symmetry information in model file is incomplete or missing.
If you want the program to generate P1 symmetry automatically, set
generate_fake_p1_symmetry=True.'''
      if not cs: raise Sorry(msg)
      elif cs.is_incomplete(): raise Sorry(msg)
    else:
      print('\nGenerating fake P1 symmetry', file=self.logger)

    pdb_hierarchy = mo.get_hierarchy()
    # need to preserve the order in the hierarchy in case we have to perform an
    # atom selection later
    # if cs is None, this will create a fake box, not sure how to get the same
    # box via model obj directly
    xray_structure = pdb_hierarchy.extract_xray_structure(crystal_symmetry = cs)
    if (cs is None): cs = xray_structure.crystal_symmetry()
    print('\nCrystal symmetry used: ', file=self.logger)
    cs.show_summary(f=self.logger)

    if (miller_array is not None):
      if (miller_array.crystal_symmetry() is None):
        miller_array = miller_array.customized_copy(crystal_symmetry=cs)
    xray_structure.show_summary(f = self.logger, prefix='  ')
    if(len(self.params.anomalous_scatterers.group) != 0):
      pdb_atoms = pdb_hierarchy.atoms()
      pdb_atoms.reset_i_seq()
      set_fp_fdp_for_anomalous_scatterers(
        pdb_hierarchy              = pdb_hierarchy,
        xray_structure             = xray_structure,
        anomalous_scatterer_groups = self.params.anomalous_scatterers.group)
    elif (self.params.wavelength is not None):
      print("Setting inelastic form factors for wavelength = %g" % \
        self.params.wavelength, file=self.logger)
      xray_structure.set_inelastic_form_factors(
        photon=self.params.wavelength,
        table="sasaki")

    if(self.params.random_seed is not None):
      random.seed(self.params.random_seed)
      flex.set_random_seed(self.params.random_seed)

    print("-"*79, file=self.logger)
    print("Computing model structure factors, Fmodel:", file=self.logger)
    if(self.params.output.format == "cns"): extension = ".hkl"
    elif(self.params.output.format == "mtz"): extension = ".mtz"
    ofn = self.params.output.file_name
    if(ofn is None):
      ofn = os.path.basename(self.data_manager.get_default_model_name())
      ofn = ofn + extension

    use_custom_scattering_dictionary = False
    if self.params.custom_scattering_factors:
      use_custom_scattering_dictionary = True
    mmtbx.utils.fmodel_from_xray_structure(
      xray_structure = xray_structure,
      f_obs          = miller_array,
      add_sigmas     = self.params.add_sigmas,
      params         = self.params,
      twin_law       = self.params.twin_law,
      twin_fraction  = self.params.twin_fraction,
      use_custom_scattering_dictionary = use_custom_scattering_dictionary,
      out            = self.logger).write_to_file(file_name = ofn,
        obs_type=self.params.output.obs_type)
    print("Output file name:", ofn, file=self.logger)
    print("All done.", file=self.logger)
    print("-"*79, file=self.logger)
    self.output_file = ofn

  def get_results(self):
    return group_args(
     output_file=self.output_file,
     model=self.data_manager.get_default_model_name())

