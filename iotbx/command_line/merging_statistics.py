"""Analyze merging statistics for a reflection file"""
# LIBTBX_SET_DISPATCHER_NAME phenix.merging_statistics
# LIBTBX_SET_DISPATCHER_NAME iotbx.merging_statistics

from __future__ import absolute_import, division, print_function
import iotbx.merging_statistics
import iotbx.phil
from libtbx.str_utils import format_value
from libtbx.utils import Sorry
from libtbx import runtime_utils
import sys

citations_str = iotbx.merging_statistics.citations_str

master_phil = """
file_name = None
  .type = path
  .short_caption = Unmerged data
  .style = file_type:hkl OnChange:extract_unmerged_intensities bold
labels = None
  .type = str
  .input_size = 200
  .style = renderer:draw_unmerged_intensities_widget
space_group = None
  .type = space_group
  .input_size = 120
unit_cell = None
  .type = unit_cell
symmetry_file = None
  .type = path
  .style = file_type:pdb,hkl OnChange:extract_symmetry
%s
debug = False
  .type = bool
loggraph = False
  .type = bool
estimate_cutoffs = False
  .type = bool
include scope libtbx.phil.interface.tracking_params
json {
  file_name = None
    .type = path
  indent = None
    .type = int(value_min=0)
}
mmcif {
  file_name = None
    .type = path
  data_name = data
    .type = str
}
""" % iotbx.merging_statistics.merging_params_str

# Hack for handling SHELX files
class cmdline_processor(iotbx.phil.process_command_line_with_files):
  def process_other(self, arg):
    if ("=" in arg):
      fields = arg.split("=")
      if (len(fields) == 2) and (fields[1] in ["amplitudes", "intensities",
          "hklf3", "hklf4"]):
        from iotbx import reflection_file_reader
        hkl_in = reflection_file_reader.any_reflection_file(arg)
        if (hkl_in.file_type() is not None):
          return iotbx.phil.parse("%s=%s" % (self.reflection_file_def,
            arg))
    return False

def run(args, out=None, master_params=None,
    assume_shelx_observation_type_is="intensities"):
  if (out is None) : out = sys.stdout
  import iotbx.phil
  if (master_params is None):
    master_params = iotbx.phil.parse(master_phil, process_includes=True)
  cmdline = cmdline_processor(
    args=args,
    master_phil=master_params,
    reflection_file_def="file_name",
    pdb_file_def="symmetry_file",
    space_group_def="space_group",
    unit_cell_def="unit_cell",
    usage_string="""\
phenix.merging_statistics [data_file] [options...]

Calculate merging statistics for non-unique data, including R-merge, R-meas,
R-pim, and redundancy.  Any format supported by Phenix is allowed, including
MTZ, unmerged Scalepack, or XDS/XSCALE (and possibly others).  Data should
already be on a common scale, but with individual observations unmerged.
%s
""" % citations_str)
  params = cmdline.work.extract()
  i_obs = iotbx.merging_statistics.select_data(
    file_name=params.file_name,
    data_labels=params.labels,
    log=out,
    assume_shelx_observation_type_is=assume_shelx_observation_type_is,
    anomalous=params.anomalous,
  )
  params.labels = i_obs.info().label_string()
  validate_params(params)
  symm = sg = uc = None
  if (params.symmetry_file is not None):
    from iotbx import crystal_symmetry_from_any
    symm = crystal_symmetry_from_any.extract_from(
      file_name=params.symmetry_file)
    if (symm is None):
      raise Sorry("No symmetry records found in %s." % params.symmetry_file)
  else :
    sg = i_obs.space_group()
    if (params.space_group is not None):
      sg = params.space_group.group()
    elif (sg is None):
      raise Sorry("Missing space group information.")
    uc = i_obs.unit_cell()
    if (params.unit_cell is not None):
      uc = params.unit_cell
    elif (uc is None):
      raise Sorry("Missing unit cell information.")
    from cctbx import crystal
    symm = crystal.symmetry(
      space_group=sg,
      unit_cell=uc)
  if (i_obs.sigmas() is None):
    raise Sorry("Sigma(I) values required for this application.")
  result = iotbx.merging_statistics.dataset_statistics(
    i_obs=i_obs,
    crystal_symmetry=symm,
    d_min=params.high_resolution,
    d_max=params.low_resolution,
    n_bins=params.n_bins,
    reflections_per_bin=params.reflections_per_bin,
    binning_method=params.binning_method,
    anomalous=params.anomalous,
    debug=params.debug,
    file_name=params.file_name,
    sigma_filtering=params.sigma_filtering,
    use_internal_variance=params.use_internal_variance,
    eliminate_sys_absent=params.eliminate_sys_absent,
    extend_d_max_min=params.extend_d_max_min,
    cc_one_half_significance_level=params.cc_one_half_significance_level,
    cc_one_half_method=params.cc_one_half_method,
    log=out)
  result.show(out=out)
  if (getattr(params, "loggraph", False)):
    result.show_loggraph(out=out)
  if (params.estimate_cutoffs):
    result.show_estimated_cutoffs(out=out)
  if params.json.file_name is not None:
    result.as_json(file_name=params.json.file_name, indent=params.json.indent)
  if params.mmcif.file_name is not None:
    import iotbx.cif.model
    cif = iotbx.cif.model.cif()
    cif[params.mmcif.data_name] = result.as_cif_block()
    with open(params.mmcif.file_name, 'w') as f:
      print(cif, file=f)
  print("", file=out)
  print("References:", file=out)
  print(citations_str, file=out)
  print("", file=out)
  return result

#-----------------------------------------------------------------------
# Phenix GUI stuff
def validate_params(params):
  if (params.file_name is None):
    raise Sorry("No data file specified!")
  elif (params.labels is None):
    raise Sorry("No data labels selected!")
  if (not None in [params.high_resolution, params.low_resolution]):
    if (params.low_resolution < params.high_resolution):
      raise Sorry("Resolution limits flipped - high resolution must be a "+
        "smaller number than low resolution.")
  elif (params.extend_d_max_min):
    raise Sorry("High and low resolution limits must be explicitly given "+
      "when calculating statistics relative to user-defined resolution "+
      "range (extend_d_max_min=True).")
  return True

class launcher(runtime_utils.target_with_save_result):
  def run(self):
    return run(args=list(self.args),
               out=sys.stdout,
               assume_shelx_observation_type_is="intensities")

def finish_job(result):
  stats = []
  if (result is not None):
    stats = [
      ("High resolution", format_value("%.3g", result.overall.d_min)),
      ("Redundancy", format_value("%.1f", result.overall.mean_redundancy)),
      ("R-meas", format_value("%.3g", result.overall.r_meas)),
      ("R-meas (high-res)", format_value("%.3g", result.bins[-1].r_meas)),
      ("<I/sigma>", format_value("%.2g", result.overall.i_over_sigma_mean)),
      ("<I/sigma> (high-res)", format_value("%.2g",
        result.bins[-1].i_over_sigma_mean)),
      ("Completeness", format_value("%.1f%%", result.overall.completeness*100)),
      ("Completeness (high-res)", format_value("%.1f%%",
        result.bins[-1].completeness*100)),
    ]
  return ([], stats)

if (__name__ == "__main__"):
  run(sys.argv[1:])
