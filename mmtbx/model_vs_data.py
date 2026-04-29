"""Compare model with x-ray data"""
from __future__ import absolute_import, division, print_function
import sys, random
from cctbx.array_family import flex
from iotbx import pdb
from libtbx.utils import Sorry
from iotbx import reflection_file_utils
from libtbx.str_utils import format_value
import iotbx
from mmtbx import utils
from iotbx import pdb
from six.moves import cStringIO as StringIO
import mmtbx.model
import iotbx.phil
import mmtbx.f_model
from iotbx import extract_xtal_data

if(1):
  random.seed(0)
  flex.set_random_seed(0)

def show_header(l, log):
  print(l, "-"*(79-len(l)), file=log)

def reflection_file_server(crystal_symmetry, reflection_files):
  return reflection_file_utils.reflection_file_server(
    crystal_symmetry=crystal_symmetry,
    force_symmetry=True,
    reflection_files=reflection_files,
    err=StringIO())

msg="""\
Inputs:
  - File with reflection data (Fobs or Iobs), and R-free flags (optionally);
  - label(s) selecting which reflection data arrays should be used (in case
    there are multiple choices in input file, there is no need to provide labels
    otherwise);
  - PDB file with input model;
  - some other optional parameters.

Usage examples:
  1. phenix.model_vs_data model.pdb data.hkl
  2. phenix.model_vs_data model.pdb data.hkl f_obs_label="F" r_free_flags_label="FREE"
  3. phenix.model_vs_data model.pdb data.hkl scattering_table=neutron
  4. phenix.model_vs_data model.pdb data.hkl twin_law='h,-k,l+h'

"""

master_params_str="""\
f_obs_label = None
  .type = str
r_free_flags_label = None
  .type = str
scattering_table = wk1995  it1992  *n_gaussian  neutron
  .type = choice
high_resolution = None
  .type = float
twin_law = None
  .type = str
n_bins = 20
  .type = int
"""

def defaults(log, silent):
  if(not silent): print("Default params:\n", file=log)
  parsed = iotbx.phil.parse(master_params_str)
  if(not silent): parsed.show(prefix="  ", out=log)
  if(not silent): print(file=log)
  return parsed

def run(args,
        out = None,
        log = sys.stdout):
  if(len(args)==0) or (args == ["--help"]):
    print(msg, file=log)
    defaults(log=log, silent=False)
    return
  parsed = defaults(log=log, silent=True)
  #
  processed_args = utils.process_command_line_args(args = args,
    log = log, master_params = parsed)
  params = processed_args.params.extract()
  #
  reflection_files = processed_args.reflection_files
  if(len(reflection_files) == 0):
    raise Sorry("No reflection file found.")
  crystal_symmetry = processed_args.crystal_symmetry
  if(crystal_symmetry is None):
    raise Sorry("No crystal symmetry found.")
  if(len(processed_args.pdb_file_names) == 0):
    raise Sorry("No PDB file found.")
  pdb_file_names = processed_args.pdb_file_names
  #
  rfs = reflection_file_server(
    crystal_symmetry = crystal_symmetry,
    reflection_files = reflection_files)
  parameters = extract_xtal_data.data_and_flags_master_params().extract()
  if(params.f_obs_label is not None):
    parameters.labels = params.f_obs_label
  if(params.r_free_flags_label is not None):
    parameters.r_free_flags.label = params.r_free_flags_label
  if (params.high_resolution is not None):
    parameters.high_resolution = params.high_resolution
  determine_data_and_flags_result = extract_xtal_data.run(
    reflection_file_server = rfs,
    parameters             = parameters,
    keep_going             = True)
  f_obs = determine_data_and_flags_result.f_obs
  # Data
  show_header(l="Data:", log=log)
  f_obs.show_comprehensive_summary(prefix="  ", f = log)
  # R-free-flags
  show_header(l="R-free-flags:", log=log)
  r_free_flags = determine_data_and_flags_result.r_free_flags
  test_flag_value = determine_data_and_flags_result.test_flag_value
  if(r_free_flags is None):
    r_free_flags=f_obs.array(data=flex.bool(f_obs.data().size(), False))
    test_flag_value=None
    print("  not available", file=log)
  else:
    print("  flag value:", test_flag_value, file=log)
  # Model
  pdb_combined = iotbx.pdb.combine_unique_pdb_files(
    file_names=processed_args.pdb_file_names)
  pdb_combined.report_non_unique(out=log)
  if (len(pdb_combined.unique_file_names) == 0):
    raise Sorry("No coordinate file given.")
  raw_records = pdb_combined.raw_records
  try:
    pdb_inp = iotbx.pdb.input(source_info = None,
                              lines       = flex.std_string(raw_records))
  except ValueError as e :
    raise Sorry("Model format (PDB or mmCIF) error:\n%s" % str(e))
  model = mmtbx.model.manager(
    model_input      = pdb_inp,
    crystal_symmetry = crystal_symmetry,
    log              = StringIO())
  #
  scattering_table = params.scattering_table
  exptl_method = pdb_inp.get_experiment_type()
  if exptl_method.is_neutron():
    scattering_table = "neutron"
  model.setup_scattering_dictionaries(
    scattering_table = scattering_table,
    d_min            = f_obs.d_min())
  #
  # Model vs data
  #
  show_header(l="Model vs Data:", log=log)
  fmodel = mmtbx.f_model.manager(
    xray_structure = model.get_xray_structure(),
    f_obs          = f_obs,
    r_free_flags   = r_free_flags,
    twin_law       = params.twin_law)
  fmodel.update_all_scales(update_f_part1=True)
  fmodel.show(log=log, show_header=False, show_approx=False)
  print("  r_work: %6.4f"%fmodel.r_work(), file=log)
  if(test_flag_value is not None):
    print("  r_free: %6.4f"%fmodel.r_free(), file=log)
  else:
    print("  r_free: None", file=log)
  print(file=log)
  n_outl = f_obs.data().size() - fmodel.f_obs().data().size()
  print("  Number of F-obs outliers:", n_outl, file=log)
  #
  # Extract information from PDB file header and output (if any)
  #
  pub_r_work       = None
  pub_r_free       = None
  pub_high         = None
  pub_low          = None
  pub_sigma        = None
  pub_program_name = None
  pub_solv_cont    = None
  pub_matthews     = None
  published_results = pdb_inp.get_r_rfree_sigma(file_name=pdb_file_names[0])
  if(published_results is not None):
    pub_r_work = published_results.r_work
    pub_r_free = published_results.r_free
    pub_high   = published_results.high
    pub_low    = published_results.low
    pub_sigma  = published_results.sigma
  pub_program_name = pdb_inp.get_program_name()
  pub_solv_cont    = pdb_inp.get_solvent_content()
  pub_matthews     = pdb_inp.get_matthews_coeff()
  #
  show_header(l="Information extracted from PDB file header:", log=log)
  print("  program_name    : %-s"%format_value("%s",pub_program_name), file=log)
  print("  year            : %-s"%format_value("%s",pdb_inp.extract_header_year()), file=log)
  print("  r_work          : %-s"%format_value("%s",pub_r_work), file=log)
  print("  r_free          : %-s"%format_value("%s",pub_r_free), file=log)
  print("  high_resolution : %-s"%format_value("%s",pub_high), file=log)
  print("  low_resolution  : %-s"%format_value("%s",pub_low), file=log)
  print("  sigma_cutoff    : %-s"%format_value("%s",pub_sigma), file=log)
  print("  matthews_coeff  : %-s"%format_value("%s",pub_matthews), file=log)
  print("  solvent_cont    : %-s"%format_value("%s",pub_solv_cont), file=log)
  if(exptl_method is not None):
    print("  exptl_method    : %-s"%format_value("%s", exptl_method), file=log)
  #
  # Recompute R-factors using published cutoffs
  fmodel_cut = fmodel
  tmp_sel = flex.bool(fmodel.f_obs().data().size(), True)
  if(pub_sigma is not None and fmodel.f_obs().sigmas() is not None):
    tmp_sel &= fmodel.f_obs().data() > fmodel.f_obs().sigmas()*pub_sigma
  if(pub_high is not None and abs(pub_high-fmodel.f_obs().d_min()) > 0.03):
    tmp_sel &= fmodel.f_obs().d_spacings().data() > pub_high
  if(pub_low is not None and abs(pub_low-fmodel.f_obs().d_max_min()[0]) > 0.03):
    tmp_sel &= fmodel.f_obs().d_spacings().data() < pub_low
  if(tmp_sel.count(True) != tmp_sel.size() and tmp_sel.count(True) > 0):
    show_header(l="After applying resolution and sigma cutoffs:", log=log)
    fmodel = mmtbx.f_model.manager(
      xray_structure = model.get_xray_structure(),
      f_obs          = fmodel.f_obs().select(tmp_sel),
      r_free_flags   = fmodel.r_free_flags().select(tmp_sel),
      twin_law       = params.twin_law)
    fmodel.update_all_scales(update_f_part1=True)
    fmodel.show(log=log, show_header=False, show_approx=False)
    print("  r_work: %6.4f"%fmodel.r_work(), file=log)
    if(test_flag_value is not None):
      print("  r_free: %6.4f"%fmodel.r_free(), file=log)
    else:
      print("  r_free: None", file=log)
    print(file=log)
    n_outl = f_obs.data().size() - fmodel.f_obs().data().size()
    print("  Number of F-obs outliers:", n_outl, file=log)
  return fmodel

