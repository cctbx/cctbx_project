from cctbx.array_family import flex
from libtbx.utils import Sorry, null_out
import iotbx.phil
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
import iotbx.symmetry
import iotbx.file_reader
from cStringIO import StringIO
import sys, os
import mmtbx.f_model
from iotbx.option_parser import iotbx_option_parser
from cctbx import miller
from mmtbx import utils
from iotbx.pdb import combine_unique_pdb_files
import iotbx.pdb
from libtbx import runtime_utils
import libtbx.callbacks # import dependency

fo_minus_fo_master_params_str = """\
f_obs_1_file_name = None
  .type = path
  .help = File with Fobs data
  .short_caption = Reflections file 1
  .style = bold file_type:hkl process_hkl child:fobs:f_obs_1_label force_data
f_obs_1_label = None
  .type = str
  .short_caption = Labels 1
  .input_size = 160
  .style = bold renderer:draw_fobs_label_widget
f_obs_2_file_name = None
  .type = path
  .help = File with Fobs data
  .short_caption = Reflections file 2
  .style = bold file_type:hkl process_hkl child:fobs:f_obs_2_label force_data
f_obs_2_label = None
  .type = str
  .short_caption = Labels 2
  .input_size = 160
  .style = bold renderer:draw_fobs_label_widget
high_resolution = None
  .type = float
  .help = High resolution data cutoff
  .style = bold resolution
low_resolution = None
  .type = float
  .help = Low resolution data cutoff
  .style = bold resolution
sigma_cutoff = None
  .type = float
  .help = Fobs sigma cutoff
  .short_caption = F(obs) sigma cutoff
phase_source = None
  .type = path
  .help = PDB file with a model or reflection file with the phases
  .short_caption = PDB file for phasing
  .style = bold OnChange:validate_phase_source file_type:pdb
scattering_table = *xray neutron
  .type = choice(multi=False)
  .help = Choices of scattering table for structure factors calculations
output_file = None
  .type = path
  .style = bold new_file file_type:mtz
ignore_non_isomorphous_unit_cells = False
  .type = bool
  .short_caption = Ignore non-isomorphous unit cells
include scope libtbx.phil.interface.tracking_params
find_peaks_holes = False
  .type = bool
  .short_caption = Find peaks and holes in map
map_cutoff = 3.0
  .type = float
  .short_caption = Sigma-level cutoff for peak search
peak_search
  .short_caption = Peak search parameters
  .style = auto_align box menu_item
{
include scope mmtbx.find_peaks.master_params
}
"""
def fo_minus_fo_master_params():
  return iotbx.phil.parse(fo_minus_fo_master_params_str, process_includes=True)

def compute_fo_minus_fo_map(data_arrays, xray_structure, log, silent,
    output_file=None, peak_search=False, map_cutoff=None,
    peak_search_params=None):
  fmodels = []
  for i_seq, d in enumerate(data_arrays):
    if(not silent):
      print >> log, "Data set: %d"%i_seq
    if(d.anomalous_flag()):
      d = d.average_bijvoet_mates()
    r_free_flags = d.array(data = flex.bool(d.data().size(), False))
    fmodel = mmtbx.f_model.manager(
      xray_structure = xray_structure,
      r_free_flags   = r_free_flags,
      target_name    = "ls_wunit_k1",
      f_obs          = d)
    fmodel.remove_outliers()
    fmodel.update_solvent_and_scale()
    if(not silent):
      fmodel.info().show_rfactors_targets_scales_overall()
      print >> log
    fmodels.append(fmodel)
  # prepare Fobs for map calculation (apply scaling):
  f_obss = []
  for fmodel in fmodels:
    obs = fmodel.f_obs()
    f_obs_scale   = 1.0 / fmodel.k_anisotropic() / fmodel.k_isotropic()
    obs = miller.array(miller_set = fmodel.f_model(),
                       data       = obs.data()*f_obs_scale)
    f_obss.append(obs)
  # given two Fobs sets, make them one-to-one matching, get phases and map coefficients
  # Note: f_calc below is just f_calc from atoms (no bulk solvent etc applied)
  fobs_1, f_model = f_obss[0].common_sets(other = fmodels[1].f_model())
  fobs_1, fobs_2 = fobs_1.common_sets(other = f_obss[1])
  fobs_1, f_model = fobs_1.common_sets(other = f_model)
  assert fobs_2.indices().all_eq(fobs_1.indices())
  assert f_model.indices().all_eq(fobs_1.indices())
  # scale again
  scale_k1 = 1
  den = flex.sum(flex.abs(fobs_2.data())*flex.abs(fobs_2.data()))
  if(den != 0):
    scale_k1 = flex.sum(flex.abs(fobs_1.data())*flex.abs(fobs_2.data())) / den
  #
  fobs_2 = fobs_2.array(data = fobs_2.data()*scale_k1)
  if 0: fobs_1 = fobs_2.multiscale(other = fobs_1, reflections_per_bin=250)
  if(not silent):
    print >> log, "Fobs1_vs_Fobs2 statistics:"
    print >> log, "Bin# Resolution range  Compl.  No.of refl. R-factor"
    fobs_1.setup_binner(reflections_per_bin = 500)
    fobs_2.use_binning_of(fobs_1)
    for i_bin in fobs_1.binner().range_used():
      sel = fobs_1.binner().selection(i_bin)
      f1  = fobs_1.select(sel)
      f2  = fobs_2.select(sel)
      d_max, d_min = fobs_1.d_max_min()
      compl = fobs_1.completeness(d_max = d_max)
      n_ref = sel.count(True)
      r = flex.sum(flex.abs(f1.data()-f2.data())) / \
        flex.sum(flex.abs(f1.data()+f2.data())/2)
      d_range = fobs_1.binner().bin_legend(
                     i_bin = i_bin, show_bin_number = False, show_counts = False)
      fmt = "%3d: %-17s   %4.2f %6d         %6.4f"
      print >> log, fmt % (i_bin, d_range, compl, n_ref, r)
  # map coefficients
  diff = miller.array(
    miller_set = f_model,
    data       = fobs_1.data()-fobs_2.data())
  def phase_transfer(miller_array, phase_source):
    tmp = miller.array(miller_set = miller_array,
      data = flex.double(miller_array.indices().size(), 1)
      ).phase_transfer(phase_source = phase_source)
    return miller.array(miller_set = miller_array,
      data = miller_array.data() * tmp.data() )
  map_coeff = phase_transfer(
    miller_array = diff,
    phase_source = f_model)
  # output MTZ file with map coefficients
  class map_coeffs_mtz_label_manager:
    def __init__(self, amplitudes, phases):
      self._amplitudes = amplitudes
      self._phases = phases
    def amplitudes(self):
      return self._amplitudes
    def phases(self, root_label, anomalous_sign=None):
      assert anomalous_sign is None or not anomalous_sign
      return self._phases
  mtz_history_buffer = flex.std_string()
  lbl_mgr = map_coeffs_mtz_label_manager(amplitudes = "FoFo", phases = "PHFc")
  if(map_coeff.anomalous_flag()):
    map_coeff = map_coeff.average_bijvoet_mates()
  mtz_dataset = map_coeff.as_mtz_dataset(
    column_root_label=lbl_mgr.amplitudes(),
    label_decorator=lbl_mgr)
  mtz_history_buffer.append("> column label %s = phenix %s" % (
      lbl_mgr.amplitudes(), "FoFoPHFc"))
  if output_file is not None :
    file_name = output_file
  else :
    file_name = "FoFoPHFc.mtz"
  mtz_history_buffer.append("file name %s"%file_name)
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.add_history(mtz_history_buffer)
  mtz_object.write(file_name=file_name)
  file_names = [ file_name ]
  if (peak_search) :
    from mmtbx.command_line import find_peaks_holes
    from mmtbx import find_peaks
    peak_search_log = log
    if (silent) : peak_search_log = null_out()
    peaks = find_peaks.manager(
      fmodel=fmodel,
      map_type=None,
      map_coeffs=map_coeff,
      map_cutoff=map_cutoff,
      params=peak_search_params,
      log=peak_search_log).peaks_mapped()
    peaks.sites = fmodel.xray_structure.unit_cell().orthogonalize(peaks.sites)
    holes = find_peaks.manager(
      fmodel=fmodel,
      map_type=None,
      map_coeffs=map_coeff,
      map_cutoff=-map_cutoff,
      params=peak_search_params,
      log=peak_search_log).peaks_mapped()
    holes.sites = fmodel.xray_structure.unit_cell().orthogonalize(holes.sites)
    result = find_peaks_holes.peaks_holes_container(
      peaks=peaks,
      holes=holes,
      map_cutoff=map_cutoff)
    pdb_out = os.path.splitext(file_name)[0] + "_peaks.pdb"
    result.save_pdb_file(
      file_name=pdb_out,
      include_anom=False,
      include_water=False,
      log=peak_search_log)
    file_names.append(pdb_out)
  return file_names

def run(args, command_name = "phenix.fobs_minus_fobs_map"):
  if(len(args) == 0): args = ["--help"]
  examples = """Examples:

phenix.fobs_minus_fobs_map f_obs_1_file=data1.mtz f_obs_2_file=data2.sca \
f_obs_1_label=FOBS1 f_obs_2_label=FOBS2 model.pdb

phenix.fobs_minus_fobs_map f_obs_1_file=data.mtz f_obs_2_file=data.mtz \
f_obs_1_label=FOBS1 f_obs_2_label=FOBS2 phase_source=model.pdb \
high_res=2.0 sigma_cutoff=2 scattering_table=neutron"""

  command_line = (iotbx_option_parser(
    usage="%s [options]" % command_name,
    description=examples)
    .option("--silent",
      action="store_true",
      help="Suppress output to the screen.")
    .enable_symmetry_comprehensive()
    ).process(args=args)
  #
  log = sys.stdout
  if(not command_line.options.silent):
    utils.print_header("phenix.fobs_minus_fobs_map", out = log)
    print >> log, "Command line arguments: "
    print >> log, args
    print >> log
  #
  processed_args = utils.process_command_line_args(
    args=command_line.args,
    cmd_cs=command_line.symmetry,
    master_params=fo_minus_fo_master_params(),
    log=log,
    suppress_symmetry_related_errors=True)
  working_phil = processed_args.params
  if(not command_line.options.silent):
    print >> log, "*** Parameters:"
    working_phil.show(out = log)
    print >> log
  params = working_phil.extract()
  consensus_symmetry = None
  if (params.ignore_non_isomorphous_unit_cells) :
    if (None in [params.f_obs_1_file_name, params.f_obs_2_file_name,
        params.phase_source]):
      raise Sorry("The file parameters (f_obs_1_file_name, f_obs_2_file_name, "+
        "phase_source) must be specified explicitly when "+
        "ignore_non_isomorphous_unit_cells=True.")
    symm_manager = iotbx.symmetry.manager()
    pdb_in = iotbx.file_reader.any_file(params.phase_source, force_type="pdb")
    symm_manager.process_pdb_file(pdb_in)
    hkl_in_1 = iotbx.file_reader.any_file(params.f_obs_1_file_name,
      force_type="hkl")
    sg_err_1, uc_err_1 = symm_manager.process_reflections_file(hkl_in_1)
    hkl_in_2 = iotbx.file_reader.any_file(params.f_obs_2_file_name,
      force_type="hkl")
    sg_err_2, uc_err_2 = symm_manager.process_reflections_file(hkl_in_2)
    out = StringIO()
    symm_manager.show(out=out)
    if (sg_err_1) or (sg_err_2) :
      raise Sorry(("Incompatible space groups in input files:\n%s\nAll files "+
        "must have the same point group (and ideally the same space group). "+
        "Please note that any symmetry information in the PDB file will be "+
        "used first.") % out.getvalue())
    elif (uc_err_1) or (uc_err_2) :
      libtbx.call_back(message="warn",
        data=("Crystal symmetry mismatch:\n%s\nCalculations will continue "+
          "using the symmetry in the PDB file (or if not available, the "+
          "first reflection file), but the maps should be treated with "+
          "extreme suspicion.") % out.getvalue())
    crystal_symmetry = symm_manager.as_symmetry_object()
  else :
    processed_args = utils.process_command_line_args(
      args=command_line.args,
      cmd_cs=command_line.symmetry,
      master_params=fo_minus_fo_master_params(),
      log=StringIO())
    crystal_symmetry = processed_args.crystal_symmetry
  #
  pdb_file_names = processed_args.pdb_file_names
  if(len(processed_args.pdb_file_names) == 0):
    if(params.phase_source is not None):
      pdb_file_names = [params.phase_source]
    else:
      raise Sorry("No PDB file found.")
  # Extaract Fobs1, Fobs2
  f_obss = []
  if(len(processed_args.reflection_files)==2):
    for reflection_file in processed_args.reflection_files:
      reflection_file_server = reflection_file_utils.reflection_file_server(
        crystal_symmetry = crystal_symmetry,
        force_symmetry   = True,
        reflection_files = [reflection_file],
        err              = null_out())
      determine_data_and_flags_result = utils.determine_data_and_flags(
        reflection_file_server  = reflection_file_server,
        keep_going              = True,
        log                     = null_out())
      f_obss.append(determine_data_and_flags_result.f_obs)
  else:
    if([params.f_obs_1_file_name,params.f_obs_2_file_name].count(None)==2):
      raise Sorry("No reflection data file found.")
    for file_name, label in zip([params.f_obs_1_file_name,params.f_obs_2_file_name],
                                [params.f_obs_1_label,params.f_obs_2_label]):
      reflection_file = reflection_file_reader.any_reflection_file(
        file_name = file_name, ensure_read_access = False)
      reflection_file_server = reflection_file_utils.reflection_file_server(
        crystal_symmetry = crystal_symmetry,
        force_symmetry   = True,
        reflection_files = [reflection_file],
        err              = null_out())
      parameters = utils.data_and_flags_master_params().extract()
      if(label is not None):
        parameters.labels = [label]
      determine_data_and_flags_result = utils.determine_data_and_flags(
          reflection_file_server  = reflection_file_server,
          parameters              = parameters,
          keep_going              = True,
          log                     = null_out())
      f_obss.append(determine_data_and_flags_result.f_obs)
  if(len(f_obss)!=2):
    raise Sorry(" ".join(errors))
  if(not command_line.options.silent):
    for ifobs, fobs in enumerate(f_obss):
      print >> log, "*** Summary for data set %d:"%ifobs
      fobs.show_comprehensive_summary(f = log)
      print >> log
  pdb_combined = combine_unique_pdb_files(file_names = pdb_file_names)
  pdb_combined.report_non_unique(out = log)
  if(len(pdb_combined.unique_file_names) == 0):
    raise Sorry("No coordinate file given.")
  #
  raw_recs = flex.std_string()
  for rec in pdb_combined.raw_records:
    if(rec.upper().count("CRYST1")==0):
      raw_recs.append(rec)
  raw_recs.append(iotbx.pdb.format_cryst1_record(
    crystal_symmetry = crystal_symmetry))
  #
  xray_structure = iotbx.pdb.input(source_info = None, lines =
    raw_recs).xray_structure_simple()
  if(not command_line.options.silent):
    print >> log, "*** Model summary:"
    xray_structure.show_summary(f = log)
    print >> log
  info0 = f_obss[0].info()
  info1 = f_obss[1].info()
  f_obss[0] = f_obss[0].resolution_filter(d_min = params.high_resolution,
    d_max = params.low_resolution).set_info(info0)
  f_obss[1] = f_obss[1].resolution_filter(d_min = params.high_resolution,
    d_max = params.low_resolution).set_info(info1)
  if(params.sigma_cutoff is not None):
    for i in [0,1]:
      if(f_obss[i].sigmas() is not None):
        sel = f_obss[i].data() > f_obss[i].sigmas()*params.sigma_cutoff
        f_obss[i] = f_obss[i].select(sel).set_info(info0)
  for k, f_obs in enumerate(f_obss) :
    if (f_obs.indices().size() == 0) :
      raise Sorry("No data left in array %d (labels=%s) after filtering!" % (k+1,
        f_obs.info().label_string()))
  output_files = compute_fo_minus_fo_map(
    data_arrays = f_obss,
    xray_structure = xray_structure,
    log = log,
    silent = command_line.options.silent,
    output_file = params.output_file,
    peak_search=params.find_peaks_holes,
    map_cutoff=params.map_cutoff,
    peak_search_params=params.peak_search)
  return output_files

class launcher (runtime_utils.target_with_save_result) :
  def run (self) :
    return run(args=list(self.args))

def validate_params (params, callback=None) :
  if (None in [params.f_obs_1_file_name, params.f_obs_2_file_name]) :
    raise Sorry("You must supply two files containing F(obs).")
  if (None in [params.f_obs_1_label, params.f_obs_2_label]) :
    raise Sorry("You must define the labels for both reflection files.")
  if (params.phase_source is None) :
    raise Sorry("You must specify a PDB file for phasing.")
  if (params.output_file is None) or (params.output_file == "") :
    raise Sorry("You must specify an output file.")
  output_dir, output_file = os.path.split(params.output_file)
  if os.path.isdir(params.output_file) :
    raise Sorry("Output file is a directory!")
  elif not output_file.endswith(".mtz") :
    raise Sorry("Output file must be an MTZ file.")
  elif not os.path.isdir(output_dir) :
    raise Sorry("Output directory does not exist.")

def finish_job (result) :
  output_files = []
  if (result is not None) :
    assert (isinstance(result, list)) and (len(result) > 0)
    output_files.append((result[0], "Map coefficients"))
    if (len(result) > 1) :
      output_files.append((result[1], "Map peaks and holes"))
  return (output_files, [])
