# LIBTBX_SET_DISPATCHER_NAME phenix.cc_star

from __future__ import division
from libtbx.str_utils import make_sub_header, format_value
from libtbx.utils import Sorry, Usage
from libtbx import runtime_utils
from libtbx import Auto
import libtbx.phil
import sys

def merging_and_model_statistics (
    f_obs,
    f_model,
    r_free_flags,
    unmerged_i_obs,
    n_bins=20,
    sigma_filtering=Auto) :
  from iotbx import merging_statistics
  free_sel = r_free_flags
  # very important: must use original intensities for i_obs, not squared f_obs,
  # because French-Wilson treatment is one-way
  assert (unmerged_i_obs.sigmas() is not None)
  info = unmerged_i_obs.info()
  assert (info is not None)
  unmerged_i_obs = unmerged_i_obs.customized_copy(crystal_symmetry=f_obs)
  unmerged_i_obs = unmerged_i_obs.select(
    unmerged_i_obs.sigmas() >= 0).set_info(info)
  filter = merging_statistics.filter_intensities_by_sigma(
    array=unmerged_i_obs,
    sigma_filtering=sigma_filtering)
  i_obs = filter.array_merged
  unmerged_i_obs = filter.array
  if (i_obs.anomalous_flag()) :
    i_obs = i_obs.average_bijvoet_mates()
  if (f_model.anomalous_flag()) :
    f_model = f_model.average_bijvoet_mates()
  if (free_sel.anomalous_flag()) :
    free_sel = free_sel.average_bijvoet_mates()
  if (free_sel.data().count(True) == 0) :
    raise Sorry("R-free array does not select any reflections.  To calculate "+
      "CC* and related statistics, a valid set of R-free flags must be used.")
  work_sel = free_sel.customized_copy(data=~free_sel.data())
  i_obs, f_model = i_obs.common_sets(other=f_model)
  i_obs, f_obs = i_obs.common_sets(other=f_obs)
  i_obs, work_sel = i_obs.common_sets(other=work_sel)
  i_obs, free_sel = i_obs.common_sets(other=free_sel)
  i_calc = abs(f_model).f_as_f_sq()
  d_max, d_min = i_calc.d_max_min()
  model_arrays = merging_statistics.model_based_arrays(
    f_obs=f_obs,
    i_obs=i_obs,
    i_calc=i_calc,
    work_sel=work_sel,
    free_sel=free_sel)
  return merging_statistics.dataset_statistics(
    i_obs=unmerged_i_obs,
    crystal_symmetry=i_calc,
    d_min=d_min,
    d_max=d_max,
    n_bins=n_bins,
    model_arrays=model_arrays,
    sigma_filtering=None) # no need, since it was done here

master_phil = libtbx.phil.parse("""
data = None
  .type = path
  .help = Data file (usually MTZ) containing R-free flags and either the \
    pre-calculated F(model) array or experimental amplitudes or intensities.
  .style = file_type:hkl input_file bold process_hkl child:ampl:f_obs_labels \
           child:rfree:r_free_flags.label child:fmodel:f_model_labels \
           force_data
f_obs_labels = None
  .type = str
  .help = Column labels for experimental data array
  .short_caption = F(obs) labels
  .input_size = 150
  .style = bold renderer:draw_fobs_label_widget
f_model_labels = None
  .type = str
  .short_caption = F(model) labels
  .style = renderer:draw_fmodel_label_widget
  .input_size = 150
r_free_flags.label = None
  .type = str
  .help = Column label for R-free flags
  .short_caption = Free R label
  .style = bold renderer:draw_rfree_label_widget
  .input_size = 150
r_free_flags.test_flag_value = None
  .type = int
  .help = Test flag value.  Not normally required.
model = None
  .type = path
  .help = Model file, required if F(model) is not pre-calculated.
  .style = file_type:pdb input_file
unmerged_data = None
  .type = path
  .help = File containing scaled, unmerged intensities
  .style = bold file_type:hkl OnChange:extract_unmerged_intensities input_file
unmerged_labels = None
  .type = str
  .help = Labels for unmerged intensity array
  .style = bold renderer:draw_unmerged_intensities_widget
  .input_size = 150
n_bins = 20
  .type = int(value_min=5, value_max=50)
  .help = Number of resolution bins
  .input_size = 64
  .style = spinner
include scope iotbx.merging_statistics.sigma_filtering_phil_str
include scope libtbx.phil.interface.tracking_params
loggraph = False
  .type = bool
""", process_includes=True)
master_params = master_phil # for phenix GUI

def show_symmetry_error (file1, file2, symm1, symm2) :
  import cStringIO
  symm_out1 = cStringIO.StringIO()
  symm_out2 = cStringIO.StringIO()
  symm1.show_summary(f=symm_out1, prefix="  ")
  symm2.show_summary(f=symm_out2, prefix="  ")
  raise Sorry("Incompatible symmetry definitions:\n%s:\n%s\n%s\n%s" %
    (file1, symm_out1.getvalue(), file2, symm_out2.getvalue()))

def load_and_validate_unmerged_data (f_obs, file_name, data_labels,
    log=sys.stdout) :
  from iotbx import merging_statistics
  unmerged_i_obs = merging_statistics.select_data(
    file_name=file_name,
    data_labels=data_labels,
    log=log)
  if ((unmerged_i_obs.space_group() is not None) and
      (unmerged_i_obs.unit_cell() is not None)) :
    if (not unmerged_i_obs.is_similar_symmetry(f_obs)) :
      show_symmetry_error("Data file", "Unmerged data", unmerged_i_obs, f_obs)
  return unmerged_i_obs

def run (args=None, params=None, out=sys.stdout) :
  assert [args, params].count(None) == 1
  if args is not None:
    if (len(args) == 0) or ("--help" in args) :
      raise Usage("""
    phenix.cc_star model.pdb data.mtz unmerged_data=data.hkl [n_bins=X] [options]
    phenix.cc_star model_refine_001.mtz unmerged_data=data.hkl [...]

  Implementation of the method for assessing data and model quality described in:
  Karplus PA & Diederichs K (2012) Science 336:1030-3.

  Full parameters:
  %s
  """ % master_phil.as_str(prefix=" ", attributes_level=1))
    import iotbx.phil
    cmdline = iotbx.phil.process_command_line_with_files(
      args=args,
      master_phil=master_phil,
      pdb_file_def="model",
      reflection_file_def="data")
    params = cmdline.work.extract()
  from iotbx import merging_statistics
  from iotbx import file_reader
  if (params.data is None) :
    raise Sorry("Please specify a data file (usually MTZ format).")
  if (params.unmerged_data is None) :
    raise Sorry("Please specify unmerged_data file")
  hkl_in = file_reader.any_file(params.data, force_type="hkl")
  hkl_in.check_file_type("hkl")
  f_model = f_obs = r_free_flags = None
  f_models = []
  data_arrays = []
  f_model_labels = []
  if (params.f_model_labels is None) :
    for array in hkl_in.file_server.miller_arrays :
      labels = array.info().label_string()
      if (array.is_complex_array()) :
        if (labels.startswith("F-model") or labels.startswith("FMODEL")) :
          f_models.append(array)
          f_model_labels.append(labels)
    if (len(f_models) > 1) :
      raise Sorry(("Multiple F(model) arrays found:\n%s\nPlease specify the "+
        "'labels' parameter.") % "\n".join(f_model_labels))
    elif (len(f_models) == 1) :
      f_model = f_models[0]
      if (f_model.anomalous_flag()) :
        info = f_model.info()
        f_model = f_model.average_bijvoet_mates().set_info(info)
      print >> out, "F(model):"
      f_model.show_summary(f=out, prefix="  ")
    else :
      data_array = hkl_in.file_server.get_xray_data(
        file_name=params.data,
        labels=params.f_obs_labels,
        ignore_all_zeros=True,
        parameter_scope="")
      if (data_array.is_xray_intensity_array()) :
        from cctbx import french_wilson
        f_obs = french_wilson.french_wilson_scale(
          miller_array=data_array,
          out=out)
      else :
        f_obs = data_array
      if (f_obs.anomalous_flag()) :
        info = f_obs.info()
        f_obs = f_obs.average_bijvoet_mates().set_info(info)
      print >> out, "F(obs):"
      f_obs.show_summary(f=out, prefix="  ")
      print >> out, ""
  else :
    for array in hkl_in.file_server.miller_arrays :
      array_labels = array.info().label_string()
      if (array_labels == params.f_model_labels) :
        if (array.is_complex_array()) :
          f_model = array
          break
        else :
          raise Sorry("The data in %s are not of the required type." %
            array_labels)
  if (f_model is not None) :
    assert (f_obs is None)
    for array in hkl_in.file_server.miller_arrays :
      labels = array.info().label_string()
      if (labels == params.f_obs_labels) :
        f_obs = array
        break
    else :
      try :
        f_obs = hkl_in.file_server.get_amplitudes(
          file_name=params.f_obs_labels,
          labels=None,
          convert_to_amplitudes_if_necessary=False,
          parameter_name="f_obs_labels",
          parameter_scope="",
          strict=True)
      except Sorry :
        raise Sorry("You must supply a file containing both F-obs and F-model "+
          "if you want to use a pre-calculated F-model array.")
  assert (f_obs.is_xray_amplitude_array())
  if (f_obs.anomalous_flag()) :
    info = f_obs.info()
    f_obs = f_obs.average_bijvoet_mates().set_info(info)
  print >> out, "F(obs):"
  f_obs.show_summary(f=out, prefix="  ")
  print >> out, ""
  r_free_flags, test_flag_value = hkl_in.file_server.get_r_free_flags(
    file_name=params.data,
    label=params.r_free_flags.label,
    test_flag_value=params.r_free_flags.test_flag_value,
    disable_suitability_test=False,
    parameter_scope="")
  info = r_free_flags.info()
  r_free_flags = r_free_flags.customized_copy(
    data=r_free_flags.data()==test_flag_value).set_info(info)
  if (r_free_flags.anomalous_flag()) :
    r_free_flags = r_free_flags.average_bijvoet_mates().set_info(info)
  print >> out, "R-free flags:"
  r_free_flags.show_summary(f=out, prefix="  ")
  print >> out, ""
  unmerged_i_obs = load_and_validate_unmerged_data(
    f_obs=f_obs,
    file_name=params.unmerged_data,
    data_labels=params.unmerged_labels,
    log=out)
  print >> out, "Unmerged intensities:"
  unmerged_i_obs.show_summary(f=out, prefix="  ")
  print >> out, ""
  if (f_model is None) :
    assert (f_obs is not None)
    if (params.model is None) :
      raise Sorry("A PDB file is required if F(model) is not pre-calculated.")
    make_sub_header("Calculating F(model)", out=out)
    pdb_in = file_reader.any_file(params.model, force_type="pdb")
    pdb_in.check_file_type("pdb")
    pdb_symm = pdb_in.file_object.crystal_symmetry()
    if (pdb_symm is None) :
      pdb_symm = f_obs
    else :
      if (not pdb_symm.is_similar_symmetry(f_obs)) :
        show_symmetry_error(
          file1="PDB file",
          file2="data file",
          symm1=pdb_symm,
          symm2=f_obs)
    xray_structure = pdb_in.file_object.xray_structure_simple(
      crystal_symmetry=pdb_symm)
    from mmtbx.utils import fmodel_simple
    # XXX this gets done anyway later, but they need to be consistent before
    # creating the fmodel manager
    if (f_obs.anomalous_flag()) :
      f_obs = f_obs.average_bijvoet_mates()
    f_obs, r_free_flags = f_obs.common_sets(other=r_free_flags)
    fmodel = fmodel_simple(
      update_f_part1_for="refinement",
      f_obs=f_obs,
      r_free_flags=r_free_flags,
      xray_structures=[xray_structure],
      scattering_table="n_gaussian")
    fmodel.show(log=out)
    f_model = fmodel.f_model()
    r_free_flags = f_model.customized_copy(data=fmodel.arrays.free_sel)
  else :
    if (f_model.anomalous_flag()) :
      f_model = f_model.average_bijvoet_mates()
    f_model, r_free_flags = f_model.common_sets(other=r_free_flags)
  stats = merging_and_model_statistics(
    f_model=f_model,
    f_obs=f_obs,
    r_free_flags=r_free_flags,
    unmerged_i_obs=unmerged_i_obs,
    n_bins=params.n_bins,
    sigma_filtering=params.sigma_filtering)
  stats.show_cc_star(out=out)
  if (params.loggraph) :
    stats.show_loggraph(out=out)
  print >> out, ""
  print >> out, "Reference:"
  print >> out, "  Karplus PA & Diederichs K (2012) Science 336:1030-3."
  print >> out, ""
  return stats

def validate_params (params) :
  if (params.data is None) or (params.f_obs_labels is None) :
    raise Sorry("No experimental data supplied!")
  if (params.f_model_labels is None) and (params.model is None) :
    raise Sorry("You must supply either a pre-calculated F(model) array, "+
      "or the current refined model.")
  return True

class launcher (runtime_utils.target_with_save_result) :
  def run (self) :
    return run(args=list(self.args), out=sys.stdout)

def finish_job (result) :
  stats = []
  if (result is not None) :
    stats = [
      ("High resolution", format_value("%.3g", result.overall.d_min)),
      ("Redundancy", format_value("%.1f", result.overall.mean_redundancy)),
      ("<I/sigma>", format_value("%.2g", result.overall.i_over_sigma_mean)),
      ("<I/sigma> (high-res)", format_value("%.2g",
        result.bins[-1].i_over_sigma_mean)),
      ("Completeness", format_value("%.1f%%", result.overall.completeness*100)),
      ("Completeness (high-res)", format_value("%.1f%%",
        result.bins[-1].completeness*100)),
      ("CC*", format_value("%.3f", result.overall.cc_star)),
      ("CC* (high-res)", format_value("%.3f", result.bins[-1].cc_star)),
      ("CC(work)", format_value("%.3f", result.overall.cc_work)),
      ("CC(work) (high-res)", format_value("%.3f", result.bins[-1].cc_work)),
      ("CC(free)", format_value("%.3f", result.overall.cc_free)),
      ("CC(free) (high-res)", format_value("%.3f", result.bins[-1].cc_free)),
    ]
  return ([], stats)

if (__name__ == "__main__") :
  run(sys.argv[1:])
