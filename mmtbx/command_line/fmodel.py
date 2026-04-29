"""Calculate structure factors from model"""
# LIBTBX_SET_DISPATCHER_NAME phenix.fmodel
from __future__ import absolute_import, division, print_function

from iotbx.cli_parser import run_program
from mmtbx.programs import fmodel

if __name__ == '__main__':
  run_program(program_class=fmodel.Program)

#from __future__ import absolute_import, division, print_function
## LIBTBX_SET_DISPATCHER_NAME phenix.fmodel
#
#import sys, os
#import mmtbx.utils
#import iotbx.phil
#import iotbx.pdb
#from scitbx.array_family import flex
#from libtbx import runtime_utils
#from libtbx.utils import Sorry
#import random
#
#legend = """
#phenix.fmodel: a tool to compute structure factors, Fmodel:
#
#  Fmodel = scale * exp(AnisoScale) * (Fcalc + k_sol * exp(-b_sol*s^2/4) * Fmask)
#
#  where:
#
#  - Fmodel - total model structure factor (complex value)
#  - AnisoScale = -ht*A(-1)*b_cart*A(-1)th/4
#  - h - column vector with Miller indices
#  - A - orthogonalization matrix
#  - b_cart - anisotropic scale matrix
#  - t and (-1) denotes transposition and inversion operations
#  - scale - overall scale factor
#  - Fcalc - structure factors calculated from atomic model
#  - k_sol and b_sol - Flat Bulk solvent model parameters
#  - Fmask - structure factors calculated from bulk solvent mask
#
#Usage examples:
#
#  1) phenix.fmodel model.pdb high_resolution=1.5
#
#     will result in a file containing complete set of Fmodel = Fcalc computed
#     from atomic model up to 1.5A resolution.
#
#  2) phenix.fmodel model.pdb scale=2 k_sol=0.35 b_sol=50 b_cart="1 2 3 0 4 7" high_res=1.5 low_res=10
#
#     will result in a file containing complete set of Fmodel computed using the
#     above formula in resolution range 1.5-20.0A.
#
#  3) phenix.fmodel model.pdb high_resolution=1.5 algorithm=direct
#
#     is similar to "1)" but the Fcalc are computed using direct summation algorithm.
#
#  4) phenix.fmodel model.pdb high_res=1.5 format=cns label=FOBS type=real r_free=0.1
#
#     will result in CNS formatted file containing complete set of amplitudes of
#     Fmodel = Fcalc computed up to 1.5A resolution, labelled as FOBS, and free-R
#     flags with 10% of test reflections. This is a typical command to simulate Fobs.
#
#  5) phenix.fmodel model.pdb high_res=1.5 scattering_table=neutron
#
#     will result in a file containing complete set of Fmodel = Fcalc computed
#     from atomic model up to 1.5A resolution using neutron scattering table.
#
#  6) phenix.fmodel model.pdb parameters.txt
#
#     will result in a structure factor file, where Fmodel were computed using
#     parameters defined in parameters.txt file. The parameters.txt file can
#     contain all or any subset of parameters listed below. Note, that each {
#     must have a matching one }.
#
#  7) phenix.fmodel model.pdb reflection_data.mtz
#
#     will result in a file containing a set of Fmodel = Fcalc that will match
#     the set of Miller indices of the data in reflection_data.mtz file.
#
#  8) phenix.fmodel model.pdb reflection_data.mtz data_column_label="FOBS,SIGMA"
#
#     similar to "7)", where the specific data array is selected.
#
#  9) phenix.fmodel model.pdb reflection_data.mtz twin_law="l,-k,h" twin_fraction=0.3
#
#     generates twin data set (real type) with given twin law and fraction.
#
#See below for complete list of available parameters.
#"""
#
#fmodel_from_xray_structure_params_str = """\
#fmodel
#  .short_caption = F(model) options
#  .expert_level = 1
#  .style = auto_align box
#{
#  k_sol = 0.0
#    .type = float
#    .help = Bulk solvent k_sol values
#    .short_caption=Bulk solvent K_sol value
#  b_sol = 0.0
#    .type = float
#    .help = Bulk solvent b_sol values
#    .short_caption=Bulk solvent B_sol value
#  b_cart = 0 0 0 0 0 0
#    .type = floats(6)
#    .help = Anisotropic scale matrix
#    .input_size = 200
#    .short_caption = Anisotropic scale matrix
#  scale = 1.0
#    .type = float
#    .help = Overall scale factor
#}
#structure_factors_accuracy
#  .short_caption = Structure factors accuracy
#  .style = auto_align box
#{
#  include scope mmtbx.f_model.sf_and_grads_accuracy_master_params
#}
#mask
#  .short_caption = Bulk solvent mask
#  .style = auto_align box
#{
#  include scope mmtbx.masks.mask_master_params
#}
#"""
#fmodel_from_xray_structure_params = iotbx.phil.parse(
#  fmodel_from_xray_structure_params_str, process_includes=True)
#
#fmodel_from_xray_structure_master_params_str = """\
#high_resolution = None
#  .type = float
#  .expert_level=1
#  .style = noauto bold
#low_resolution = None
#  .type = float
#  .expert_level=1
#  .style = noauto
#r_free_flags_fraction = None
#  .type = float
#  .expert_level=1
#  .style = noauto
#add_sigmas = False
#  .type = bool
#  .expert_level=1
#  .help = Adds calculated Sigma(F) column to output file.
#  .style = noauto
#add_random_error_to_amplitudes_percent = None
#  .type = float
#  .short_caption = Add random error (percent)
#  .style = noauto
#scattering_table = wk1995  it1992  *n_gaussian  neutron electron
#  .type = choice
#  .help = Choices of scattering table for structure factors calculations.  \
#    n_gaussian is the standard set of X-ray scattering factors.
#  .expert_level=1
#  .style = noauto
#custom_scattering_factors = None
#  .type = path
#  .help = Use custom scattering factors and replaces default values entirely
#pdb_file = None
#  .type = path
#  .multiple = True
#  .optional = True
#  .short_caption = Model file
#  .style = bold noauto file_type:pdb input_file OnChange:update_output_file_name
#reference_file = None
#  .type = path
#  .short_caption = Reference set
#  .help = Reflections file containing Miller indices (h,k,l) to use in output \
#    file.
#  .style = noauto input_file file_type:mtz OnChange:update_reference_column_labels
#data_column_label = None
#  .type = str
#  .short_caption = Reference file label
#  .style = noauto renderer:draw_any_label_widget
#%s
#random_seed=None
#  .type = int
#  .help = Random seed
#  .expert_level=2
#twin_law = None
#  .type = str
#  .help = Optional twin law if we want to generate a twinned dataset
#  .input_size = 120
#  .style = noauto
#twin_fraction = None
#  .type = float
#  .help = Twin fraction, ignored if twin_law is not specified
#  .style = noauto
#wavelength = None
#  .type = float
#  .input_size = 80
#  .help = Wavelength, sets all atoms to anomalous
#  .style = noauto
#generate_fake_p1_symmetry = False
#  .type = bool
#  .short_caption = Generate fake symmetry if necessary
#  .help = Allows use of PDB files without CRYST1 records as input.  The \
#    crystal symmetry will be assumed to be a P1 box.
#output
#  .short_caption = Reflection output
#  .expert_level=0
#  .style = noauto
#{
#  format = *mtz cns
#    .type = choice
#    .short_caption = File format
#  label = FMODEL
#    .type = str
#    .short_caption = Data label
#    .input_size = 100
#  type = real *complex
#    .type = choice
#    .short_caption = Output data type
#    .help = Numeric type of output data.  'real' is amplitudes only, \
#      'complex' is complete structure factors as complex numbers.
#    .expert_level=1
#    .style = bold
#  obs_type = *amplitudes intensities
#    .type = choice
#    .help = Experimental observation type to output.  Certain restrictions \
#      apply if intensities are selected.
#    .expert_level = 2
#  file_name = None
#    .type = path
#    .short_caption = Output file
#    .style = bold noauto new_file
#  include scope libtbx.phil.interface.tracking_params
#}
#anomalous_scatterers
#  .short_caption = Anomalous sites
#  .style = menu_item noauto
#{
#  group
#    .optional = True
#    .multiple = True
#    .short_caption = Anomalous scatterer group
#    .style = auto_align
#  {
#    selection = None
#      .type = atom_selection
#      .short_caption = Atom selection
#      .input_size = 400
#    f_prime = 0
#      .type = float
#      .short_caption = f'
#    f_double_prime = 0
#      .type = float
#      .short_caption = f''
#  }
#}
#"""%fmodel_from_xray_structure_params_str
#
#fmodel_from_xray_structure_master_params = iotbx.phil.parse(
#  fmodel_from_xray_structure_master_params_str, process_includes=True)
#master_phil = fmodel_from_xray_structure_master_params # XXX for phenix docs
#
#def set_fp_fdp_for_anomalous_scatterers(pdb_hierarchy, xray_structure,
#  anomalous_scatterer_groups):
#  scatterers = xray_structure.scatterers()
#  for group in anomalous_scatterer_groups:
#    iselection = pdb_hierarchy.atom_selection_cache().selection(
#      string = group.selection).iselection()
#    if(iselection.size() == 0):
#      raise Sorry(
#        "Empty selection: selection string '%s' does not select any atom."%
#        group.selection)
#    for i_seq in iselection:
#      scatterers[i_seq].fp = group.f_prime
#      scatterers[i_seq].fdp = group.f_double_prime
#
#def run(args, log = sys.stdout):
#  print(legend, file=log)
#  # XXX: pre-processing for GUI; duplicates some of mmtbx.utils
#  sources = []
#  for arg in args :
#    if os.path.isfile(arg):
#      try :
#        file_phil = iotbx.phil.parse(file_name=arg)
#      except KeyboardInterrupt :
#        raise
#      except RuntimeError :
#        pass
#      else :
#        if len(file_phil.objects) != 0 :
#          sources.append(file_phil)
#  if len(sources) > 0 :
#    cmdline_phil = fmodel_from_xray_structure_master_params.fetch(
#      sources=sources)
#    params = cmdline_phil.extract()
#    if len(params.pdb_file) > 0 :
#      args.extend(params.pdb_file)
#    if params.reference_file is not None :
#      args.append(params.reference_file)
#  # end of preprocessing
#  processed_args = mmtbx.utils.process_command_line_args(args = args, log = log,
#    master_params = fmodel_from_xray_structure_master_params)
#  pdb_combined = iotbx.pdb.combine_unique_pdb_files(
#    file_names = processed_args.pdb_file_names)
#  pdb_combined.report_non_unique(out = log)
#  print("-"*79, file=log)
#  print("\nParameters to compute Fmodel::\n", file=log)
#  processed_args.params.show(out = log, prefix=" ")
#  params = processed_args.params.extract()
#  if(params.random_seed is not None):
#    random.seed(params.random_seed)
#    flex.set_random_seed(params.random_seed)
#  pdb_file_names = processed_args.pdb_file_names
#  if len(pdb_file_names) == 0 :
#    pdb_file_names = params.pdb_file # for GUI
#  pdb_combined = iotbx.pdb.combine_unique_pdb_files(file_names=pdb_file_names)
#  pdb_combined.report_non_unique(out = log)
#  if(len(pdb_combined.unique_file_names) == 0):
#    raise Sorry("Model file is not provided.")
#  print("-"*79, file=log)
#  print("\nInput model file(s):", " ".join(processed_args.pdb_file_names), file=log)
#  pdb_inp = iotbx.pdb.input(source_info = None,
#    lines = flex.std_string(pdb_combined.raw_records))
#  # select miller array to use as a set of miller indices for f_model
#  miller_array = None
#  if(len(processed_args.reflection_files) > 1):
#    raise Sorry("Multiple reflection files found at input.")
#  # FIXME this does not pick up the reference_file parameter!  in fact, it
#  # appears to be ignored completely when run on the command line.
#  if(len(processed_args.reflection_files) == 1):
#    print("-"*79, file=log)
#    print("Input reflection data:", \
#      " ".join(processed_args.reflection_file_names), file=log)
#    if([params.high_resolution, params.low_resolution].count(None) != 2):
#      raise Sorry("high_resolution and low_resolution must be undefined "+
#                  "if reflection data file is given.")
#    miller_arrays = processed_args.reflection_files[0].as_miller_arrays()
#    data_sizes = flex.int([ma.data().size() for ma in miller_arrays])
#    if(data_sizes.all_eq(data_sizes[0])): miller_array = miller_arrays[0]
#    else:
#      all_labels = []
#      for ma in miller_arrays:
#        if(params.data_column_label is not None and
#           ma.info().label_string() == params.data_column_label):
#          miller_array = ma
#          break
#        all_labels.append(",".join(ma.info().labels))
#    if(miller_array is None):
#      raise Sorry("Multiple data available in input reflection file:\n%s\n%s"%(
#        "\n".join(all_labels),"Please select one using 'data_column_label=' keyword."))
#    else:
#      miller_array.show_comprehensive_summary(f = log, prefix="  ")
#  if(miller_array is not None):
#    miller_array = miller_array.map_to_asu().customized_copy(
#      data = flex.double(miller_array.data().size(), 1))
#  #
#  cryst1 = pdb_inp.crystal_symmetry_from_cryst1()
#  if(cryst1 is None and miller_array is not None):
#    cryst1 = miller_array.crystal_symmetry()
#    if (cryst1 is not None) and (params.generate_fake_p1_symmetry):
#      raise Sorry("The input reference data already define crystal symmetry; "+
#        "you may not use this in combination with the option "+
#        "generate_fake_p1_symmetry=True.")
#  if (not params.generate_fake_p1_symmetry):
#    if(cryst1 is None):
#      raise Sorry(
#        "CRYST1 record in input PDB file is incomplete or missing.  "+
#        "If you want the program to generate P1 symmetry automatically, set "+
#        "generate_fake_p1_symmetry=True.")
#    else:
#      if([cryst1.unit_cell(), cryst1.space_group_info()].count(None) != 0):
#        raise Sorry(
#          "CRYST1 record in input PDB file is incomplete or missing. "+
#          "If you want the program to generate P1 symmetry automatically, "+
#          "set generate_fake_p1_symmetry=True.")
#  pdb_hierarchy = pdb_inp.construct_hierarchy()
#  # need to preserve the order in the hierarchy in case we have to perform an
#  # atom selection later
#  xray_structure = pdb_hierarchy.extract_xray_structure(
#    crystal_symmetry = cryst1)
#  if (cryst1 is None):
#    cryst1 = xray_structure.crystal_symmetry()
#  if (miller_array is not None):
#    if (miller_array.crystal_symmetry() is None):
#      miller_array = miller_array.customized_copy(crystal_symmetry=cryst1)
#  xray_structure.show_summary(f = log, prefix="  ")
#  if(len(params.anomalous_scatterers.group) != 0):
#    pdb_atoms = pdb_hierarchy.atoms()
#    pdb_atoms.reset_i_seq()
#    set_fp_fdp_for_anomalous_scatterers(
#      pdb_hierarchy              = pdb_hierarchy,
#      xray_structure             = xray_structure,
#      anomalous_scatterer_groups = params.anomalous_scatterers.group)
#  elif (params.wavelength is not None):
#    if (params.scattering_table == "neutron"):
#      raise Sorry("Wavelength parameter not supported when the neutron "+
#        "scattering table is used.")
#    print("Setting inelastic form factors for wavelength = %g" % \
#      params.wavelength, file=log)
#    xray_structure.set_inelastic_form_factors(
#      photon=params.wavelength,
#      table="sasaki")
#  #
#  validate_params_command_line(params)
#  #
#  print("-"*79, file=log)
#  print("Computing model structure factors, Fmodel:", file=log)
#  if(params.output.format == "cns"): extension = ".hkl"
#  elif(params.output.format == "mtz"): extension = ".mtz"
#  ofn = params.output.file_name
#  if(ofn is None):
#    ofn = os.path.basename(processed_args.pdb_file_names[0])
#    if(len(processed_args.pdb_file_names)==1): ofn = ofn + extension
#    else: ofn = ofn + "_et_al" + extension
#  if([miller_array, params.high_resolution].count(None)==2):
#    raise Sorry("Input data file or high_resolution has to be provided.")
#  use_custom_scattering_dictionary = False
#  if params.custom_scattering_factors:
#    use_custom_scattering_dictionary = True
#  mmtbx.utils.fmodel_from_xray_structure(
#    xray_structure = xray_structure,
#    f_obs          = miller_array,
#    add_sigmas     = params.add_sigmas,
#    params         = params,
#    twin_law       = params.twin_law,
#    twin_fraction  = params.twin_fraction,
#    use_custom_scattering_dictionary = use_custom_scattering_dictionary,
#    out            = log).write_to_file(file_name = ofn,
#      obs_type=params.output.obs_type)
#  print("Output file name:", ofn, file=log)
#  print("All done.", file=log)
#  print("-"*79, file=log)
#  return ofn
#
#class launcher(runtime_utils.target_with_save_result):
#  def run(self):
#    return run(args=list(self.args), log=sys.stdout)
#
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
#
#def validate_params_command_line(params):
#  if (params.output.type == "complex") and (params.add_sigmas):
#    raise Sorry("Sigma values only supported when the output type is 'real'.")
#  if (    params.low_resolution is not None
#      and params.high_resolution is not None):
#    if params.low_resolution < params.high_resolution :
#      raise Sorry("Low-resolution cutoff must be larger than the high-"+
#        "resolution cutoff.")
#  if (params.output.obs_type == "intensities"):
#    if (params.output.type == "complex"):
#      raise Sorry("Output type must be 'real' when intensities specified "+
#        "for obs_type.")
#    if (not params.output.label.upper().startswith("I")):
#      raise Sorry("Output label must start with 'I' (any case) when "+
#        "intensities specified for obs_type (was: %s)." % params.output.label)
#    if (params.output.format != "mtz"):
#      raise Sorry("Output format must be 'mtz' when intensities specified.")
#  return True
#
#def finish_job(result):
#  output_files = []
#  if (result is not None) and (os.path.isfile(result)):
#    output_files.append((os.path.abspath(result), "MTZ file"))
#  return (output_files, [])
#
#if (__name__ == "__main__"):
#  run(args=sys.argv[1:])
#

