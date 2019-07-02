
from __future__ import absolute_import, division, print_function
from libtbx import Auto
import os.path
import sys

master_phil_str = """
map_type = *2mFo-DFc mFo-DFc anom anom_residual llg
  .type = choice
exclude_free_r_reflections = True
  .type = bool
fill_missing_f_obs = False
  .type = bool
b_sharp = None
  .type = float
output_file = Auto
  .type = path
"""

map_type_labels = {
  "2mFo-DFc" : "2FOFCWT",
  "mFo-DFc" : "FOFCWT",
  "anom" : "ANOM",
  "anom_residual" : "ANOM_DIFF",
  "llg" : "LLG",
}

def master_phil():
  from mmtbx.command_line import generate_master_phil_with_inputs
  return generate_master_phil_with_inputs(
    phil_string=master_phil_str,
    enable_automatic_twin_detection=True)

def run(args, out=sys.stdout):
  master_params = master_phil()
  usage_string = """\
mmtbx.compute_map_coefficients model.pdb data.mtz map_type=MAP_TYPE [options]

Utility to compute a single set of map coefficients with minimal input.
"""
  import mmtbx.command_line
  from mmtbx import map_tools
  import iotbx.mtz
  cmdline = mmtbx.command_line.load_model_and_data(
    args=args,
    master_phil=master_params,
    process_pdb_file=False,
    prefer_anomalous=True,
    usage_string=usage_string,
    set_wavelength_from_model_header=True,
    set_inelastic_form_factors="sasaki",
    out=out)
  fmodel = cmdline.fmodel
  xray_structure = fmodel.xray_structure
  params = cmdline.params
  map_coeffs = fmodel.map_coefficients(
    map_type=params.map_type,
    exclude_free_r_reflections=params.exclude_free_r_reflections,
    fill_missing=params.fill_missing_f_obs)
  if (params.b_sharp is not None):
    if (params.b_sharp is Auto):
      params.b_sharp = - map_coeffs.d_min() * 10
    map_coeffs = map_tools.sharp_map(
      sites_frac=None,
      map_coeffs=map_coeffs,
      b_sharp=b_sharp)
  dec = iotbx.mtz.label_decorator(phases_prefix="PH")
  mtz_dataset = map_coeffs.as_mtz_dataset(
    column_root_label=map_type_labels[params.map_type],
    label_decorator=dec)
  if (params.output_file is Auto):
    pdb_file = os.path.basename(params.input.pdb.file_name[0])
    params.output_file = os.path.splitext(pdb_file)[0] + "_%s.mtz" % \
      params.map_type
  mtz_dataset.mtz_object().write(params.output_file)
  print("Wrote %s" % params.output_file, file=out)

if (__name__ == "__main__"):
  run(sys.argv[1:])
