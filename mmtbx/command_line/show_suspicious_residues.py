
from __future__ import division
from __future__ import print_function
import libtbx.load_env
from libtbx import Auto
import sys
import os

def master_phil():
  from mmtbx.command_line import generate_master_phil_with_inputs
  return generate_master_phil_with_inputs(
    enable_automatic_twin_detection=True,
    phil_string="""
hetatms_only = True
  .type = bool
skip_xtal_solution_mols = False
  .type = bool
  .help = If True, common crystallization components such as sulfate or \
    glycerol will be ignored.
skip_single_atoms = True
  .type = bool
skip_alt_confs = True
  .type = bool
min_acceptable_cc = 0.8
  .type = float
min_acceptable_2fofc = 1.0
  .type = float
max_frac_atoms_below_min = 0.5
  .type = float
write_coot_script = True
  .type = bool
write_maps = Auto
  .type = bool
""")

common_xtal_mols = [
  "SO4", "GOL", "EDO", "PEG", "ACT", "BME",
]

def run(args, out=sys.stdout):
  usage_string = """\
mmtbx.show_suspicious_residues [model] [data] [options...]

Flag bad residues based on fit to electron density.  By default, only HETATM
records will be inspected.
"""
  from mmtbx import real_space_correlation
  from mmtbx.command_line import load_model_and_data
  cmdline = load_model_and_data(
    args=args,
    master_phil=master_phil(),
    out=out,
    process_pdb_file=False,
    usage_string=usage_string)
  params = cmdline.params
  ignore_list = []
  if (params.skip_xtal_solution_mols):
    ignore_list = common_xtal_mols
  outliers = real_space_correlation.find_suspicious_residues(
    fmodel=cmdline.fmodel,
    pdb_hierarchy=cmdline.pdb_hierarchy,
    hetatms_only=params.hetatms_only,
    skip_single_atoms=params.skip_single_atoms,
    skip_alt_confs=params.skip_alt_confs,
    ignore_resnames=ignore_list,
    min_acceptable_cc=params.min_acceptable_cc,
    min_acceptable_2fofc=params.min_acceptable_2fofc,
    max_frac_atoms_below_min=params.max_frac_atoms_below_min,
    log=out)
  map_file = None
  if ((params.write_maps == True) or
      ((len(outliers) > 0) and (params.write_maps in [Auto, True]))):
    import mmtbx.maps.utils
    import iotbx.map_tools
    f_map, diff_map = mmtbx.maps.utils.get_maps_from_fmodel(cmdline.fmodel)
    anom_map = None
    if (cmdline.fmodel.f_obs().anomalous_flag()):
      anom_map = mmtbx.maps.utils.get_anomalous_map(cmdline.fmodel)
    base_name = os.path.basename(
      os.path.splitext(params.input.xray_data.file_name)[0])
    map_file = base_name + "_maps.mtz"
    iotbx.map_tools.write_map_coeffs(
      fwt_coeffs=f_map,
      delfwt_coeffs=diff_map,
      file_name=map_file,
      anom_coeffs=anom_map)
    print("Wrote maps to %s" % map_file)
  if (len(outliers) > 0) and (params.write_coot_script):
    zoom_list_base = libtbx.env.find_in_repositories(
      relative_path="cctbx_project/cootbx/simple_zoom_list.py",
      test=os.path.isfile)
    script = open("coot_bad_residues.py", "w")
    script.write(open(zoom_list_base).read())
    script.write("\n")
    for file_name in params.input.pdb.file_name :
      script.write("""read_pdb("%s")\n""" % file_name)
    if (map_file is not None):
      script.write("auto_read_make_and_draw_maps(\"%s\")\n" % map_file)
    script.write("""
draw_simple_zoom_list(
  title="Residues in suspicious density",
  items=%s)
""" % str(outliers))
    script.close()
    print("Coot script is coot_bad_residues.py")

if (__name__ == "__main__"):
  run(sys.argv[1:])
