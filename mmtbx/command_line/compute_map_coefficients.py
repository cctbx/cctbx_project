
from __future__ import division
from libtbx.utils import Usage
from libtbx import Auto
import libtbx.phil
import os.path
import sys

master_phil_str = """
include scope mmtbx.utils.cmdline_input_phil_str
map_type = *2mFo-DFc mFo-DFc anom anom_residual llg
  .type = choice
exclude_free_r_reflections = True
  .type = bool
fill_missing_f_obs = False
  .type = bool
wavelength = None
  .type = float
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

def master_phil () :
  return libtbx.phil.parse(master_phil_str, process_includes=True)

def run (args, out=sys.stdout) :
  master_params = master_phil()
  if (len(args) == 0) :
    raise Usage("""\
mmtbx.compute_map_coefficients model.pdb data.mtz map_type=MAP_TYPE [options]

Utility to compute a single set of map coefficients with minimal input.

Full parameters:
%s
""" % master_params.as_str(prefix="  "))
  import mmtbx.utils
  from mmtbx import map_tools
  import iotbx.mtz
  cmdline = mmtbx.utils.cmdline_load_pdb_and_data(
    update_f_part1_for="map",
    args=args,
    master_phil=master_params,
    process_pdb_file=False,
    prefer_anomalous=True,
    out=out)
  fmodel = cmdline.fmodel
  xray_structure = fmodel.xray_structure
  params = cmdline.params
  if (params.wavelength is not None) :
    xray_structure.set_inelastic_form_factors(
      photon=params.wavelength,
      table="sasaki")
    fmodel.update_xray_structure(xray_structure, update_f_calc=True)
  map_coeffs = fmodel.map_coefficients(
    map_type=params.map_type,
    exclude_free_r_reflections=params.exclude_free_r_reflections,
    fill_missing=params.fill_missing_f_obs)
  if (params.b_sharp is not None) :
    if (params.b_sharp is Auto) :
      params.b_sharp = - map_coeffs.d_min() * 10
    map_coeffs = map_tools.sharp_map(
      sites_frac=None,
      map_coeffs=map_coeffs,
      b_sharp=b_sharp)
  dec = iotbx.mtz.label_decorator(phases_prefix="PH")
  mtz_dataset = map_coeffs.as_mtz_dataset(
    column_root_label=map_type_labels[params.map_type],
    label_decorator=dec)
  if (params.output_file is Auto) :
    pdb_file = os.path.basename(params.input.pdb.file_name[0])
    params.output_file = os.path.splitext(pdb_file)[0] + "_%s.mtz" % \
      params.map_type
  mtz_dataset.mtz_object().write(params.output_file)
  print >> out, "Wrote %s" % params.output_file

if (__name__ == "__main__") :
  run(sys.argv[1:])
