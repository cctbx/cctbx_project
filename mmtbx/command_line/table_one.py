from __future__ import absolute_import, division, print_function
from mmtbx.validation import dummy_validation
from mmtbx.validation import molprobity
import mmtbx.command_line.molprobity
from iotbx import file_reader
import iotbx.table_one
from libtbx import object_oriented_patterns as oop
from libtbx.str_utils import make_header
from libtbx.utils import Sorry, null_out
from libtbx import runtime_utils
import libtbx.phil.command_line
from libtbx import easy_pickle
from libtbx import easy_mp
from libtbx import Auto
import libtbx.phil
import os.path
import time
import sys
from six.moves import zip
from six.moves import range

structure_params_str = """
  structure
    .multiple = True
    .optional = True
    .short_caption = Structure files
    .style = auto_align box
  {
    name = None
      .type = str
      .input_size = 400
      .short_caption = Structure name
      .help = Structure name (will become column label)
      .style = bold
    pdb_file = None
      .type = path
      .short_caption = PDB file
      .style = bold file_type:pdb
      .help = PDB file
    mtz_file = None
      .type = path
      .short_caption = MTZ file
      .style = bold file_type:mtz OnChange:extract_table_one_labels
      .help = MTZ file
    data_labels = None
      .type = str
      .input_size = 160
      .style = renderer:draw_table_one_label_widget
    r_free_flags_label = None
      .type = str
      .short_caption = R-free flags label
      .input_size = 160
      .style = renderer:draw_table_one_label_widget
      .help = R-free flags label
    wavelength = None
      .type = float
    cif_file = None
      .type = path
      .multiple = True
      .optional = True
      .short_caption = CIF file
      .help = Restraints file
    cif_directory = None
      .type = path
      .short_caption = CIF directory
      .style = directory
      .help = Directory containing restraints - all CIF files found will be \
        used.
    data_type = *xray neutron
      .type = choice
    unmerged_data = None
      .type = path
      .help = Secondary reflections file with unmerged intensities, for \
        calculation of merging statistics.
      .style = bold file_type:hkl OnChange:extract_unmerged_intensities \
        help_page:unmerged_data.htm
    unmerged_labels = None
      .type = str
      .input_size = 160
      .style = renderer:draw_unmerged_intensities_widget
    use_internal_variance = False
      .type = bool
      .help = Estimate intensity variance for unmerged data
    count_anomalous_pairs_separately = False
      .type = bool
      .short_caption = Count anomalous pairs separately
      .help = If true, the program will treat F+ and F- (if present) as \
        independent reflections when calculating data statistics.  (Not \
        recommended.)
    enable_twinning = False
      .type = bool
      .help = Enable twinning
    twin_law = Auto
      .type = str
      .short_caption = Twin law
      .help = Twin law if available, automatic detection if not
  }
"""

master_phil_str = """
table_one {
  %s
  processing {
    re_compute_r_factors = True
      .type = bool
      .short_caption = Always re-compute R-factors
      .style = bold
    n_bins = 10
      .type = int
      .short_caption = Number of resolution bins
    ligand_selection = None
      .type = str
      .short_caption = Ligand atom selection
      .help = If specified, will determine which atoms or residues are \
        counted as ligands (instead of automatic behavior).
      .input_size = 400
  }
  multiprocessing {
    include scope libtbx.easy_mp.parallel_phil_str_no_threading
  }
  output {
    directory = None
      .type = path
      .help = This is only used by the PHENIX GUI.
      .short_caption = Output directory
      .style = bold output_dir
    include scope libtbx.phil.interface.tracking_params
    show_missing_fields = True
      .type = bool
    format = txt csv *rtf
      .type = choice(multi=True)
      .caption = Text CSV RTF
      .short_caption = Output formats
      .style = bold
    base_name = Table1
      .type = str
      .short_caption = Base file name
      .style = bold
    verbose = True
      .type = str
    text_field_separation = 2
      .type = int
  }
}""" % structure_params_str

class _(oop.injector, molprobity.molprobity):
  """
  Injector dummy class to add as_table1_column() method to the main molprobity
  object.  This extracts statistics from the various validation objects and
  unmerged data, checks for consistency, and returns an iotbx.table_one.column
  object.
  """
  def as_table1_column(self,
      label,
      wavelength,
      log,
      re_compute_r_factors=Auto):
    """
    Extract information for display in the traditional 'Table 1' of
    crystallographic statistics in structure articles.
    """
    outer_shell = None
    data_stats = self.data_stats
    if (data_stats is None):
      data_stats = dummy_validation()
    merging_stats = dummy_validation()
    merging_outer = dummy_validation()
    n_refl_uniq = data_stats.n_refl
    n_refl_refine = data_stats.n_refl_refine
    n_free = data_stats.n_free
    completeness = data_stats.completeness
    completeness_outer = data_stats.completeness_outer
    d_max_min = self.d_max_min()
    d_max, d_min = d_max_min
    if (self.merging is not None):
      merging_stats = self.merging.overall
      merging_outer = self.merging.bins[-1]
      n_refl_uniq = merging_stats.n_uniq
      epsilon = 0.001
      if ((merging_stats.d_min > d_min + 2*epsilon) or
          (merging_stats.d_max < d_max - 2*epsilon)):
        raise Sorry(("Resolution limits for unmerged data in the structure "+
          "'%s' do not cover the "+
          "full range present in the merged data: %g - %g (merged) versus "+
          "%g - %g (unmerged)") % (label, d_max, d_min,
            merging_stats.d_max, merging_stats.d_min))
    r_work = self.r_work()
    r_free = self.r_free()
    n_tls_groups = None
    if (self.header_info is not None):
      if (self.header_info.n_tls_groups > 0):
        n_tls_groups = self.header_info.n_tls_groups
      use_header_values = (not re_compute_r_factors or
          (not self.header_info.is_phenix_refinement() and
           (re_compute_r_factors is Auto)))
      r_work, r_free, warned = rfactor_sanity_check(
        r_work_pdb=self.header_info.r_work,
        r_free_pdb=self.header_info.r_free,
        r_work_fmodel=r_work,
        r_free_fmodel=r_free,
        out=log,
        structure_name=label,
        re_compute_r_factors=not use_header_values)
      if (use_header_values):
        n_refl_refine = data_stats.n_refl
    adp_result = self.adp_stats.result()
    adp_mean = [None for i in range(4)]
    for i,prop in enumerate(['overall', 'protein', 'other', 'water']):
      if getattr(adp_result, prop) is not None:
        adp_mean[i] = getattr(adp_result, prop).mean
    return iotbx.table_one.column(
      label=label,
      space_group=self.space_group_info(),
      unit_cell=self.unit_cell().parameters(),
      # data properties
      wavelength=wavelength,
      d_max_min=d_max_min,
      n_refl_all=merging_stats.n_obs,
      n_refl=n_refl_uniq,
      multiplicity=merging_stats.mean_redundancy,
      completeness=completeness * 100.0,
      i_over_sigma=merging_stats.i_over_sigma_mean,
      wilson_b=data_stats.wilson_b,
      r_sym=merging_stats.r_merge,
      r_meas=merging_stats.r_meas,
      r_pim=merging_stats.r_pim,
      cc_one_half=merging_stats.cc_one_half,
      cc_star=merging_stats.cc_star,
      # refinement
      n_refl_refine=n_refl_refine,
      n_free=n_free,
      r_work=r_work,
      r_free=r_free,
      cc_work=merging_stats.cc_work,
      cc_free=merging_stats.cc_free,
      # model properties
      n_atoms=self.model_stats_new.result().n_atoms - self.model_stats_new.result().n_hd,
      n_macro_atoms=self.model_stats_new.result().n_protein_atoms + self.model_stats_new.result().n_nucleotide_atoms,
      n_ligand_atoms=self.model_stats_new.result().n_other_atoms,
      n_waters=self.model_stats_new.result().n_water_atoms,
      n_residues=self.model_stats_new.result().n_protein,
      bond_rmsd=self.rms_bonds(),
      angle_rmsd=self.rms_angles(),
      rama_favored=self.rama_favored(),
      rama_allowed=self.rama_allowed(),
      rama_outliers=self.rama_outliers(),
      rota_outliers=self.rota_outliers(),
      clashscore=self.clashscore(),
      adp_mean=adp_mean[0],
      adp_mean_mm=adp_mean[1],
      adp_mean_lig=adp_mean[2],
      adp_mean_wat=adp_mean[3],
      n_tls_groups=n_tls_groups,
      anomalous_flag=data_stats.anomalous_flag,
      ).add_outer_shell(
        # XXX we need a consistency check here as well
        d_max_min=(data_stats.d_max_outer, data_stats.d_min_outer),
        n_refl=data_stats.n_refl_outer,
        n_refl_all=merging_outer.n_obs,
        n_refl_refine=data_stats.n_refl_refine_outer,
        n_free=data_stats.n_free_outer,
        cc_one_half=merging_outer.cc_one_half,
        cc_star=merging_outer.cc_star,
        r_sym=merging_outer.r_merge,
        r_meas=merging_outer.r_meas,
        r_pim=merging_outer.r_pim,
        i_over_sigma=merging_outer.i_over_sigma_mean,
        multiplicity=merging_outer.mean_redundancy,
        completeness=completeness_outer * 100,
        cc_work=merging_outer.cc_work,
        cc_free=merging_outer.cc_free,
        r_work=data_stats.r_work_outer,
        r_free=data_stats.r_free_outer)

# XXX This is possibly problematic.  We should encourage the deposition of as
# much data as possible, regardless of what was actually used in the final
# refinement.
def resolution_sanity_check(
    d_max_pdb,
    d_min_pdb,
    d_max_f_obs,
    d_min_f_obs,
    structure_name,
    out):
  warned = False
  if (d_max_pdb is not None):
    (d_max, d_min) = (d_max_pdb, d_min_pdb)
    print("Using resolution limits reported in PDB file for structure %s" % \
      structure_name, file=out)
    if ((d_max != d_max_pdb) or (d_min != d_min_pdb)):
      warned = True
      print("""\
*** WARNING: Resolution limits in the PDB file for structure '%s'
             are inconsistent with resolution limits in the MTZ file:
             %.2f - %.2f (in PDB file)
             %.2f - %.2f (in MTZ file)
    This is not a fatal error, since the data processing may output more data
    than are used in refinement, but you should check to make sure that you
    have not accidentally specified the wrong model file.
""" % (structure_name, d_max_pdb, d_min_pdb, d_max_f_obs, d_min_f_obs), file=out)
  else :
    (d_max, d_min) = (d_max_f_obs, d_min_f_obs)
  return (d_max, d_min, warned)

def rfactor_sanity_check(
    r_work_pdb,
    r_free_pdb,
    r_work_fmodel,
    r_free_fmodel,
    out,
    structure_name,
    re_compute_r_factors,
    tolerance_ignore=0.001,
    tolerance_warn=0.004):
  warned = False
  if (r_work_pdb is not None) and (r_free_pdb is not None):
    if (not re_compute_r_factors):
      print("Using R-factors in PDB header", file=out)
      (r_work, r_free) = (r_work_pdb, r_free_pdb)
    else :
      print("Using re-computed R-factors", file=out)
      (r_work, r_free) = (r_work_fmodel, r_free_fmodel)
    delta_r_work = abs(r_work_pdb - r_work_fmodel)
    delta_r_free = abs(r_free_pdb - r_free_fmodel)
    if ((delta_r_work > tolerance_warn) or
        (delta_r_free > tolerance_warn)):
      warned = True
      print("""\
*** WARNING: R-factors reported in the PDB file do not match those calculated
             by phenix.model_vs_data:
             r_work=%.4f r_free=%.4f (in PDB file)
             r_work=%.4f r_free=%.4f (from phenix.model_vs_data)
    This is not a fatal error, but if we can't recalculate the reported
    R-factors, others probably won't be able to either.  This may be due to
    inconsistent use of a twin operator; phenix.model_vs_data may apply one
    automatically, but only if it significantly improves R-work.  Alternately,
    it can indicate that the model was further modified after the last round
    of refinement, that the wrong MTZ file was used, or that a different
    resolution limit was used for refinement.
    If applicable, use F-obs-filtered data from MTZ file from last refinement run.
""" % (r_work_pdb, r_free_pdb, r_work_fmodel, r_free_fmodel), file=out)
    elif ((delta_r_work > tolerance_ignore) or
          (delta_r_free > tolerance_ignore)):
      print("""\
    NOTE: small discrepancy in R-factors reported in PDB file versus those
          calculated by phenix.model_vs_data:
             r_work=%.4f r_free=%.4f (in PDB file)
             r_work=%.4f r_free=%.4f (from phenix.model_vs_data)
    This is probably not a big deal, but please double-check that you used the
    correct data file and have not changed the resolution limit for refinement.
    If applicable, use F-obs-filtered data from MTZ file from last refinement run.
""" % (r_work_pdb, r_free_pdb, r_work_fmodel, r_free_fmodel), file=out)
  else :
    print("Using re-computed R-factors for structure %s" % \
      structure_name, file=out)
    (r_work, r_free) = (r_work_fmodel, r_free_fmodel)
  return (r_work, r_free, warned)

def run_single_structure(params,
    n_bins,
    ligand_selection=None,
    log=None):
  if (log is None) : log = null_out()
  twin_law = Auto
  skip_twin_detection = not params.enable_twinning
  if (params.twin_law is not None):
    twin_law = params.twin_law
    skip_twin_detection = True
  molprobity_args = [
    "pdb.file_name=\"%s\"" % params.pdb_file,
    "xray_data.file_name=\"%s\"" % params.mtz_file,
    "xray_data.labels=\"%s\"" % params.data_labels,
    "xray_data.r_free_flags.file_name=\"%s\"" % params.mtz_file,
    "xray_data.r_free_flags.label=\"%s\"" % params.r_free_flags_label,
    "n_bins=%d" % n_bins,
    "count_anomalous_pairs_separately=%s" % \
    params.count_anomalous_pairs_separately,
    "coot=False",
    "maps=False",
    "probe_dots=False",
    "use_pdb_header_resolution_cutoffs=True",
    "output.quiet=True",
    "nonbonded_distance_threshold=None",
    "input.twin_law=%s" % twin_law,
    "input.skip_twin_detection=%s" % skip_twin_detection
  ] + [ "monomers.file_name=\"%s\"" % cif for cif in params.cif_file ]
  if (ligand_selection is not None):
    molprobity_args.append("ligand_selection=\"%s\"" % ligand_selection)
  if (params.unmerged_data is not None):
    molprobity_args.extend([
      "unmerged_data.file_name=\"%s\"" % params.unmerged_data,
      "unmerged_data.labels=\"%s\"" % params.unmerged_labels,
      "unmerged_data.use_internal_variance=\"%s\"" %
      params.use_internal_variance
    ])
  if (params.data_type == "neutron"):
    molprobity_args.extend(["scattering_table=neutron", "keep_hydrogens=True"])
  if (params.cif_directory is not None):
    files = os.listdir(params.cif_directory)
    for file_name in files :
      full_path = os.path.join(params.cif_directory, file_name)
      if (os.path.isdir(full_path)) : continue
      f = file_reader.any_file(full_path)
      if (f.file_type == "cif"):
        molprobity_args.append("monomers.file_name=\"%s\"" % f.file_name)
  return mmtbx.command_line.molprobity.run(args=molprobity_args, out=log)

class table_one(iotbx.table_one.table):
  __slots__ = iotbx.table_one.table.__slots__ + [
    "output_dir", "params", "output_files", ] #"n_warnings"]

  def __init__(self, params, out=sys.stdout):
    iotbx.table_one.table.__init__(self,
      text_field_separation=params.output.text_field_separation)
    self.output_dir = os.getcwd()
    self.params = params
    self.output_files = []
    make_header("Running data analysis and validation", out=out)
    results = easy_mp.parallel_map(
      iterable=list(range(len(self.params.structure))),
      func=self.run_single_structure,
      processes=params.multiprocessing.nproc,
      method=params.multiprocessing.technology,
      preserve_exception_message=True)
    for structure, result in zip(params.structure, results):
      print("", file=out)
      print("Collecting stats for structure %s" % structure.name, file=out)
      column = result.validation.as_table1_column(
        label=structure.name,
        wavelength=structure.wavelength,
        re_compute_r_factors=params.processing.re_compute_r_factors,
        log=out)
      self.add_column(column)

  def run_single_structure(self, i_struct):
    return run_single_structure(
      params=self.params.structure[i_struct],
      n_bins=self.params.processing.n_bins,
      ligand_selection=self.params.processing.ligand_selection)

  def save_multiple(self, file_base, formats):
    for format in formats :
      file_name = "%s.%s" % (file_base, format)
      method = getattr(self, "save_%s" % format)
      method(file_name)
      self.output_files.append((file_name, "Table as '%s'" % format))

  def finish_job(self, job=None):
    return (self.output_files, [])

def extract_labels(params, out, parameter_scope="structure"):
  """
  Guess MTZ file column labels for experimental data and R-free flags.  Only
  invoked when this program is run from the command line, but the Phenix GUI
  does something similar.
  """
  for i, structure in enumerate(params.structure):
    if (structure.mtz_file is None):
      raise Sorry("Missing MTZ file for structure #%d." % (i+1))
    if ([structure.data_labels, structure.r_free_flags_label].count(None)>0):
      mtz_file = file_reader.any_file(structure.mtz_file, force_type="hkl")
      mtz_file.assert_file_type("hkl")
      server = mtz_file.file_server
      file_name = mtz_file.file_name
      if (structure.data_labels is None):
        print("Attempting to guess labels for %s..." % file_name, file=out)
        data = server.get_xray_data(
          file_name=file_name,
          labels=None,
          ignore_all_zeros=True,
          parameter_scope=parameter_scope,
          parameter_name="data_labels")
        structure.data_labels = data.info().label_string()
      if (structure.r_free_flags_label is None):
        print("Attempting to guess R-free label for %s..." % file_name, file=out)
        rfree = server.get_r_free_flags(
          file_name=file_name,
          label=None,
          test_flag_value=None,
          disable_suitability_test=False,
          parameter_scope=parameter_scope+".r_free_flags")
        structure.r_free_flags_label = rfree[0].info().label_string()

def run(args,
    out=sys.stdout,
    auto_extract_labels=True,
    use_current_directory_if_not_specified=False,
    warn=True):
  master_params = libtbx.phil.parse(master_phil_str,
    process_includes=True)
  if (len(args) == 0):
    print("""\
************************************************************************
  phenix.table_one - statistics harvesting for publication
************************************************************************

  note: this is somewhat difficult to configure on the command line at
        present; you may find it more convenient to use the PHENIX GUI.

""", file=out)
    print("# Parameter template for phenix.table_one:", file=out)
    master_params.show(out=out)
    print("# (the 'structure' scope may be copied as many times as ", file=out)
    print("#  necessary to handle multiple datasets.)", file=out)
    print("# Alternate usage:", file=out)
    print("#   phenix.table_one model.pdb data.mtz [logfile]*", file=out)
    return None
  if (warn):
    print("""
  note: this is somewhat difficult to configure on the command line at
        present; you may find it more convenient to use the PHENIX GUI.
    """, file=out)
    time.sleep(2)
  master_parmas = libtbx.phil.parse(master_phil_str)
  interpreter = libtbx.phil.command_line.argument_interpreter(
    master_phil=master_params,
    home_scope="table_one")
  file_phil = []
  cmdline_phil = []
  pdb_file = None
  mtz_file = None
  unmerged_data = None
  log_files = []
  for arg in args :
    if os.path.isfile(arg):
      f = file_reader.any_file(arg)
      if (f.file_type == "phil"):
        file_phil.append(f.file_object)
      elif (f.file_type == "pdb"):
        pdb_file = f.file_name
      elif (f.file_type == "hkl"):
        mtz_file = f.file_name
      elif (f.file_type == "txt"):
        log_files.append(f.file_name)
    else :
      if arg.startswith("unmerged_data="):
        unmerged_data = os.path.abspath("=".join(arg.split("=")[1:]))
        continue
      if arg.startswith("--"):
        arg = arg[2:] + "=True"
      try :
        arg_phil = interpreter.process(arg=arg)
      except RuntimeError :
        print("Ignoring unknown argument %s" % arg, file=out)
      else :
        cmdline_phil.append(arg_phil)
  working_phil = master_params.fetch(sources=file_phil+cmdline_phil)
  params = working_phil.extract()
  if (pdb_file is not None):
    if (len(params.table_one.structure) > 0):
      raise Sorry("You already have a structure defined in the parameter "+
        "file; to add structures, you should edit the parameters instead of "+
        "specifying additional PDB and data files on the command line.")
    if (mtz_file is None):
      raise Sorry("You have supplied a PDB file, but no corresponding MTZ "+
                  "file.")
    log_file_str = "\n".join([ "log_file=%s" % f for f in log_files ])
    structure_params = libtbx.phil.parse(structure_params_str)
    new_structure = structure_params.extract().structure[0]
    new_structure.pdb_file = pdb_file
    new_structure.mtz_file = mtz_file
    new_structure.unmerged_data = unmerged_data
    params.table_one.structure.append(new_structure)
  if auto_extract_labels :
    extract_labels(params.table_one, out=out)
  if use_current_directory_if_not_specified :
    if (params.table_one.output.directory is None):
      params.table_one.output.directory = os.getcwd()
  validate_params(params)
  if (params.table_one.multiprocessing.nproc is None):
    params.table_one.multiprocessing.nproc = 1
  final_phil = master_params.format(python_object=params)
  if params.table_one.output.verbose :
    print("", file=out)
    print("#Final effective parameters:", file=out)
    final_phil.show(out=out)
    print("#---end", file=out)
    print("", file=out)
  with open("table_one.eff", "w") as f:
    final_phil.show(out=f)
  table1 = table_one(params.table_one, out=out)
  easy_pickle.dump("%s.pkl" % params.table_one.output.base_name, table1)
  table1.save_multiple(
    file_base=params.table_one.output.base_name,
    formats=params.table_one.output.format)
  return table1

def validate_params(params):
  if (len(params.table_one.structure) == 0):
    raise Sorry("No structures defined.")
  if (not os.path.isdir(params.table_one.output.directory)):
    raise Sorry("Please specify a valid output directory.")
  for i, struct in enumerate(params.table_one.structure):
    if (None in [struct.name, struct.pdb_file, struct.mtz_file]):
      raise Sorry(("Structure #%d is missing either a PDB file or an MTZ "+
        "file or a structure name.") % (i+1))
    for cfile in struct.cif_file:
      if (file_reader.any_file(cfile).file_type == "pdb"):
        raise Sorry("A restraint cif is file expected in the CIF file field.")
#    elif (None in [struct.data_labels, struct.r_free_flags_label]):
#      raise Sorry(("Need both data labels and R-free flag label for MTZ file "+
#        "%s (structure #%d).") % (struct.mtz_file, i+1))
  if (params.table_one.output.text_field_separation < 1):
    raise Sorry("Text field separation must be at least one character.")
  if (len(params.table_one.output.format) == 0):
    raise Sorry("No formats selected for output!")
  return True

class launcher(runtime_utils.target_with_save_result):
  def run(self):
    if not os.path.exists(self.output_dir):
      os.makedirs(self.output_dir)
    os.chdir(self.output_dir)
    return run(args=list(self.args), out=sys.stdout, warn=False)

if (__name__ == "__main__"):
  run(args=sys.argv[1:],
    auto_extract_labels=True,
    use_current_directory_if_not_specified=True)
