
from __future__ import absolute_import, division, print_function
from mmtbx import utils
import mmtbx.f_model
from iotbx.option_parser import iotbx_option_parser
from iotbx.pdb import combine_unique_pdb_files
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
import iotbx.file_reader
import iotbx.symmetry
import iotbx.phil
from cctbx.array_family import flex
import iotbx.pdb
from cctbx import miller
from libtbx.str_utils import format_value
from libtbx.utils import Sorry, null_out
from libtbx import adopt_init_args
from libtbx import runtime_utils
import libtbx.callbacks # import dependency
from six.moves import cStringIO as StringIO
import os
import sys
import mmtbx.model
from six.moves import zip
from iotbx import extract_xtal_data

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
scattering_table = wk1995 it1992 *n_gaussian electron neutron
  .type = choice(multi=False)
  .help = Choices of scattering table for structure factors calculations
output_file = None
  .type = path
  .style = hidden
output_dir = None
  .type = path
  .short_caption = Output folder
  .style = output_dir
  .expert_level = 4
file_name_prefix = None
  .type = str
  .input_size = 400
  .help = Base file name for output (GUI parameter only)
  .expert_level = 4
job_id = None
  .type = int
  .style = hidden
  .help = GUI parameter only
  .expert_level = 4
ignore_non_isomorphous_unit_cells = False
  .type = bool
  .short_caption = Ignore non-isomorphous unit cells
include scope libtbx.phil.interface.tracking_params
advanced
  .expert_level = 3
  .short_caption = Advanced options
  .style = auto_align menu_item
{
  multiscale = False
    .type = bool
    .expert_level = 2
    .short_caption = Use multiscaling method
  omit_selection = None
    .type = atom_selection
    .short_caption = Omit atom selection
    .help = The selected atoms will be left out of the phasing model.
    .expert_level = 1
  anomalous = False
    .type = bool
    .short_caption = Use anomalous differences instead of amplitudes
    .help = Experimental feature: generate a difference map using the \
      anomalous differences of each dataset instead of the amplitudes, \
      similar to the log-likelihood gradient difference map in Phaser \
      SAD phasing (or simple anomaloues residual map), only using two \
      experimental datasets instead of F-obs and F-calc.
    .expert_level = 3
}
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
structure_factors_accuracy
  .short_caption = Structure factors accuracy
  .style = auto_align box
{
  include scope mmtbx.f_model.sf_and_grads_accuracy_master_params
}
"""
def fo_minus_fo_master_params():
  return iotbx.phil.parse(fo_minus_fo_master_params_str, process_includes=True)

class compute_fo_minus_fo_map(object):
  def __init__(self,
      data_arrays,
      xray_structure,
      sf_accuracy_params=None,
      log=None,
      silent=False,
      output_file=None,
      peak_search=False,
      map_cutoff=None,
      peak_search_params=None,
      r_free_arrays=None,
      write_map=True,
      multiscale=False,
      anomalous=False):
    if (log is None) : log = sys.stdout
    adopt_init_args(self, locals())
    fmodels = []
    for i_seq, d in enumerate(data_arrays):
      if(not silent):
        print("Data set: %d"%i_seq, file=log)
      if(d.anomalous_flag()) and (not anomalous):
        d = d.average_bijvoet_mates()
      elif (anomalous):
        assert d.anomalous_flag()
      if (r_free_arrays is not None) and (i_seq < len(r_free_arrays)):
        r_free_flags = r_free_arrays[i_seq]
      else :
        r_free_flags = d.array(data = flex.bool(d.data().size(), False))
      fmodel = mmtbx.f_model.manager(
        xray_structure = xray_structure.deep_copy_scatterers(),
        r_free_flags   = r_free_flags,
        target_name    = "ls_wunit_k1",
        sf_and_grads_accuracy_params = sf_accuracy_params,
        f_obs          = d)
      fmodel.update_all_scales(log=None)

      #fmodel.export( out=open("taam_%d.mtz"%i_seq,"w") )

      if(not silent):
        fmodel.info().show_rfactors_targets_scales_overall(out=log)
        print(file=log)
        fmodel.show(show_header=False, show_approx=False, log=log)

      fmodels.append(fmodel)
    self.fmodel = fmodels[0]
    # prepare Fobs for map calculation (apply scaling):
    f_obss = []
    for fmodel in fmodels:
      f_obss.append(fmodel.f_obs_scaled(include_fom=True))
    # given two Fobs sets, make them one-to-one matching, get phases and map coefficients
    # Note: f_calc below is just f_calc from atoms (no bulk solvent etc applied)
    fobs_1, f_model = f_obss[0].common_sets(other = fmodels[1].f_model())
    fobs_1, fobs_2 = fobs_1.common_sets(other = f_obss[1])
    fobs_1, f_model = fobs_1.common_sets(other = f_model)
    self.f_model = f_model
    assert fobs_2.indices().all_eq(fobs_1.indices())
    assert f_model.indices().all_eq(fobs_1.indices())

    # scale again
    scale_k1 = 1
    den = flex.sum(flex.abs(fobs_2.data())*flex.abs(fobs_2.data()))
    if(den != 0):
      scale_k1 = flex.sum(flex.abs(fobs_1.data())*flex.abs(fobs_2.data())) / den
    #
    fobs_2 = fobs_2.array(data = fobs_2.data()*scale_k1)

    if multiscale:
      #fobs_1 = fobs_2.multiscale(other = fobs_1, reflections_per_bin=50)
      fobs_1 = fobs_2.multiscale(
        other = fobs_1, reflections_per_bin=250, use_exp_scale=True)
    if(not silent):
      print("", file=log)
      print("Fobs1_vs_Fobs2 statistics:", file=log)
      print("Bin# Resolution range  Compl.  No.of refl. CC   R-factor", file=log)
      fobs_1.setup_binner(reflections_per_bin = min(50, fobs_1.data().size()))
      fobs_2.use_binning_of(fobs_1)
      for i_bin in fobs_1.binner().range_used():
        sel = fobs_1.binner().selection(i_bin)
        f1  = fobs_1.select(sel)
        f2  = fobs_2.select(sel)
        d_max, d_min = fobs_1.d_max_min()
        compl = fobs_1.completeness(d_max = d_max)
        n_ref = sel.count(True)
        num = flex.sum(flex.abs(f1.data()-f2.data()))
        den = flex.sum(flex.abs(f1.data()+f2.data())/2)
        r = None
        if(den!=0):
          r = num/den
        cc = flex.linear_correlation(x=f1.data(), y=f2.data()).coefficient()
        d_range = fobs_1.binner().bin_legend(
                       i_bin = i_bin, show_bin_number = False, show_counts = False)
        fmt = "%3d: %-17s   %4.2f %6d         %5.3f  %6s"
        print(fmt % (i_bin, d_range, compl, n_ref, cc,
          format_value("%6.4f", r)), file=log)
    # overall statistics
    self.cc = flex.linear_correlation(
      x=fobs_1.data(),
      y=fobs_2.data()).coefficient()
    num = flex.sum(flex.abs(fobs_1.data()-fobs_2.data()))
    den = flex.sum(flex.abs(fobs_2.data()+fobs_2.data())/2)
    self.r_factor = None
    if (den != 0):
      self.r_factor = num / den
    # map coefficients
    def phase_transfer(miller_array, phase_source):
      tmp = miller.array(miller_set = miller_array,
        data = flex.double(miller_array.indices().size(), 1)
        ).phase_transfer(phase_source = phase_source)
      return miller.array(miller_set = miller_array,
        data = miller_array.data() * tmp.data() )
    if (not anomalous):
      diff = miller.array(
        miller_set = f_model,
        data       = fobs_1.data()-fobs_2.data())
      self.map_coeff = phase_transfer(
        miller_array = diff,
        phase_source = f_model)
    else :
      dano_1 = fobs_1.anomalous_differences()
      dano_2 = fobs_2.anomalous_differences()
      assert dano_1.indices().all_eq(dano_2.indices())
      diff = miller.array(
        miller_set = dano_1,
        data = dano_1.data() - dano_2.data())
      f_model_phases = f_model.average_bijvoet_mates().common_set(diff)
      map_coeffs = phase_transfer(
        miller_array = diff,
        phase_source = f_model_phases)
      self.map_coeff = map_coeffs.customized_copy(data=map_coeffs.data()/(2j))
    if(self.map_coeff.anomalous_flag()):
      self.map_coeff = map_coeff.average_bijvoet_mates()
    self.file_names = []
    if (write_map):
      self.file_names = self.write_map_file()

  def write_map_file(self):
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
    mtz_dataset = self.map_coeff.as_mtz_dataset(
      column_root_label=lbl_mgr.amplitudes(),
      label_decorator=lbl_mgr)
    mtz_history_buffer.append("> column label %s = phenix %s" % (
        lbl_mgr.amplitudes(), "FoFoPHFc"))
    if self.output_file is not None :
      file_name = self.output_file
    else :
      file_name = "FoFoPHFc.mtz"
    mtz_history_buffer.append("file name %s"%file_name)
    mtz_object = mtz_dataset.mtz_object()
    mtz_object.add_history(mtz_history_buffer)
    mtz_object.write(file_name=file_name)
    self.file_names = [ file_name ]
    if (self.peak_search):
      from mmtbx.command_line import find_peaks_holes
      from mmtbx import find_peaks
      peak_search_log = self.log
      if (self.silent) : peak_search_log = null_out()
      fmodel = self.fmodel
      peaks = find_peaks.manager(
        map_cutoff     = self.map_cutoff,
        xray_structure = fmodel.xray_structure,
        params         = self.peak_search_params,
        log            = peak_search_log,
        map_coeffs     = self.map_coeff).peaks_mapped()
      peaks.sites = fmodel.xray_structure.unit_cell().orthogonalize(peaks.sites)
      holes = find_peaks.manager(
        map_cutoff     = -self.map_cutoff,
        xray_structure = fmodel.xray_structure,
        params         = self.peak_search_params,
        log            = peak_search_log,
        map_coeffs     = self.map_coeff).peaks_mapped()
      holes.sites = fmodel.xray_structure.unit_cell().orthogonalize(holes.sites)
      result = find_peaks_holes.peaks_holes_container(
        peaks=peaks,
        holes=holes,
        map_cutoff=self.map_cutoff)
      pdb_out = os.path.splitext(file_name)[0] + "_peaks.pdb"
      result.save_pdb_file(
        file_name=pdb_out,
        include_anom=False,
        include_water=False,
        log=peak_search_log)
      self.file_names.append(pdb_out)
    return self.file_names

def run(args, command_name = "phenix.fobs_minus_fobs_map", log=None):
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
  if (log is None):
    log = sys.stdout
  if(not command_line.options.silent):
    utils.print_header("phenix.fobs_minus_fobs_map", out = log)
    print("Command line arguments: ", file=log)
    print(args, file=log)
    print(file=log)
  #
  processed_args = utils.process_command_line_args(
    args=command_line.args,
    cmd_cs=command_line.symmetry,
    master_params=fo_minus_fo_master_params(),
    absolute_angle_tolerance=5,
    absolute_length_tolerance=1,
    log=log,
    suppress_symmetry_related_errors=True)
  working_phil = processed_args.params
  if(not command_line.options.silent):
    print("*** Parameters:", file=log)
    working_phil.show(out = log)
    print(file=log)
  params = working_phil.extract()
  consensus_symmetry = None
  if (params.ignore_non_isomorphous_unit_cells):
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
    if (sg_err_1) or (sg_err_2):
      raise Sorry(("Incompatible space groups in input files:\n%s\nAll files "+
        "must have the same point group (and ideally the same space group). "+
        "Please note that any symmetry information in the PDB file will be "+
        "used first.") % out.getvalue())
    elif (uc_err_1) or (uc_err_2):
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
      suppress_symmetry_related_errors = False,
      absolute_angle_tolerance=5,
      absolute_length_tolerance=1,
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
      # XXX UGLY !!!
      try:
        parameters = extract_xtal_data.data_and_flags_master_params().extract()
        if(params.f_obs_1_label is not None):
          parameters.labels = [params.f_obs_1_label]
        determine_data_and_flags_result = extract_xtal_data.run(
          reflection_file_server = reflection_file_server,
          keep_going             = True,
          parameters             = parameters)
      except: # intentional
        parameters = extract_xtal_data.data_and_flags_master_params().extract()
        if(params.f_obs_2_label is not None):
          parameters.labels = [params.f_obs_2_label]
        determine_data_and_flags_result = extract_xtal_data.run(
          reflection_file_server = reflection_file_server,
          keep_going             = True,
          parameters             = parameters)
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
      parameters = extract_xtal_data.data_and_flags_master_params().extract()
      if(label is not None):
        parameters.labels = [label]
      determine_data_and_flags_result = extract_xtal_data.run(
          reflection_file_server = reflection_file_server,
          parameters             = parameters,
          keep_going             = True)
      f_obss.append(determine_data_and_flags_result.f_obs)
  if(len(f_obss)!=2):
    raise Sorry(" ".join(errors))
  if(not command_line.options.silent):
    for ifobs, fobs in enumerate(f_obss):
      print("*** Summary for data set %d:"%ifobs, file=log)
      fobs.show_comprehensive_summary(f = log)
      print(file=log)
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
  pdb_in = iotbx.pdb.input(source_info = None, lines = raw_recs)
  model = mmtbx.model.manager(model_input = pdb_in)
  d_min = min(f_obss[0].d_min(), f_obss[1].d_min())
  model.setup_scattering_dictionaries(
    scattering_table = params.scattering_table,
    d_min            = d_min)
  xray_structure = model.get_xray_structure()
  hierarchy = model.get_hierarchy()
  #
  omit_sel = flex.bool(hierarchy.atoms_size(), False)
  if (params.advanced.omit_selection is not None):
    print("Will omit selection from phasing model:", file=log)
    print("  " + params.advanced.omit_selection, file=log)
    omit_sel = hierarchy.atom_selection_cache().selection(
      params.advanced.omit_selection)
    print("%d atoms selected for removal" % omit_sel.count(True), file=log)
  del hierarchy
  xray_structure = xray_structure.select(~omit_sel)
  if(not command_line.options.silent):
    print("*** Model summary:", file=log)
    xray_structure.show_summary(f = log)
    print(file=log)
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
  for k, f_obs in enumerate(f_obss):
    if (f_obs.indices().size() == 0):
      raise Sorry("No data left in array %d (labels=%s) after filtering!" % (k+1,
        f_obs.info().label_string()))
  output_file_name = params.output_file
  if (output_file_name is None) and (params.file_name_prefix is not None):
    output_file_name = "%s_%s.mtz" % (params.file_name_prefix, params.job_id)
  output_files = compute_fo_minus_fo_map(
    data_arrays = f_obss,
    xray_structure = xray_structure,
    sf_accuracy_params = params.structure_factors_accuracy,
    log = log,
    silent = command_line.options.silent,
    output_file = output_file_name,
    peak_search=params.find_peaks_holes,
    map_cutoff=params.map_cutoff,
    peak_search_params=params.peak_search,
    multiscale=params.advanced.multiscale,
    anomalous=params.advanced.anomalous).file_names
  return output_files

class launcher(runtime_utils.target_with_save_result):
  def run(self):
    os.makedirs(self.output_dir)
    os.chdir(self.output_dir)
    return run(args=list(self.args))

def validate_params(params, callback=None):
  if (None in [params.f_obs_1_file_name, params.f_obs_2_file_name]):
    raise Sorry("You must supply two files containing F(obs).")
  if (None in [params.f_obs_1_label, params.f_obs_2_label]):
    raise Sorry("You must define the labels for both reflection files.")
  if (params.phase_source is None):
    raise Sorry("You must specify a PDB file for phasing.")
  if not os.path.isdir(params.output_dir):
    raise Sorry("Output directory does not exist.")

def finish_job(result):
  output_files = []
  if (result is not None):
    assert (isinstance(result, list)) and (len(result) > 0)
    output_files.append((result[0], "Map coefficients"))
    if (len(result) > 1):
      output_files.append((result[1], "Map peaks and holes"))
  return (output_files, [])
