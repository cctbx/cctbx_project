# LIBTBX_SET_DISPATCHER_NAME phenix.merging_statistics

from __future__ import division
import iotbx.merging_statistics
from libtbx.str_utils import format_value
from libtbx.utils import Sorry, Usage
from libtbx import runtime_utils
import sys

citations_str = """\
  Diederichs K & Karplus PA (1997) Nature Structural Biology 4:269-275
    (with erratum in: Nat Struct Biol 1997 Jul;4(7):592)
  Weiss MS (2001) J Appl Cryst 34:130-135.
  Karplus PA & Diederichs K (2012) Science 336:1030-3."""

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
include scope libtbx.phil.interface.tracking_params
""" % iotbx.merging_statistics.merging_params_str

def run (args, out=None) :
  if (out is None) : out = sys.stdout
  import iotbx.phil
  master_params = iotbx.phil.parse(master_phil, process_includes=True)
  if (len(args) == 0) :
    raise Usage("""\
phenix.merging_statistics [data_file] [options...]

Calculate merging statistics for non-unique data, including R-merge, R-meas,
R-pim, and redundancy.  Any format supported by Phenix is allowed, including
MTZ, unmerged Scalepack, or XDS/XSCALE (and possibly others).  Data should
already be on a common scale, but with individual observations unmerged.
%s

Full parameters:
%s
""" % (citations_str, master_params.as_str(prefix="  ")))
  import iotbx.phil
  cmdline = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil=master_params,
    reflection_file_def="file_name",
    pdb_file_def="symmetry_file")
  params = cmdline.work.extract()
  i_obs = iotbx.merging_statistics.select_data(
    file_name=params.file_name,
    data_labels=params.labels,
    log=out)
  symm = None
  if (params.symmetry_file is not None) :
    from iotbx import crystal_symmetry_from_any
    symm = crystal_symmetry_from_any.extract_from(
      file_name=params.symmetry_file)
    if (symm is None) :
      raise Sorry("No symmetry records found in %s." % params.symmetry_file)
  if (symm is None) :
    sg = i_obs.space_group()
    if (sg is None) :
      if (params.space_group is not None) :
        sg = params.space_group.group()
      else :
        raise Sorry("Missing space group information.")
    uc = i_obs.unit_cell()
    if (uc is None) :
      if (params.unit_cell is not None) :
        uc = params.unit_cell
      else :
        raise Sorry("Missing unit cell information.")
    from cctbx import crystal
    symm = crystal.symmetry(
      space_group=sg,
      unit_cell=uc)
  if (i_obs.sigmas() is None) :
    raise Sorry("Sigma(I) values required for this application.")
  result = iotbx.merging_statistics.dataset_statistics(
    i_obs=i_obs,
    crystal_symmetry=symm,
    d_min=params.high_resolution,
    d_max=params.low_resolution,
    n_bins=params.n_bins,
    anomalous=params.anomalous,
    debug=params.debug,
    file_name=params.file_name,
    log=out)
  result.show(out=out)
  if (params.loggraph) :
    result.show_loggraph(out=out)
  print >> out, ""
  print >> out, "References:"
  print >> out, citations_str
  print >> out, ""
  return result

#-----------------------------------------------------------------------
# Phenix GUI stuff
def validate_params (params) :
  if (params.file_name is None) :
    raise Sorry("No data file specified!")
  elif (params.labels is None) :
    raise Sorry("No data labels selected!")
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

if (__name__ == "__main__") :
  run(sys.argv[1:])
