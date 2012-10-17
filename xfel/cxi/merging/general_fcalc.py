from __future__ import division
import mmtbx.command_line.fmodel
import mmtbx.utils
from iotbx import file_reader
import libtbx.phil.command_line
import sys
import math

def run (params) :
  pdb_in = file_reader.any_file(params.model, force_type="pdb")
  pdb_in.assert_file_type("pdb")
  xray_structure = pdb_in.file_object.xray_structure_simple()
  xray_structure.show_summary()
  phil2 = mmtbx.command_line.fmodel.fmodel_from_xray_structure_master_params
  params2 = phil2.extract()
  # adjust the cutoff of the generated intensities to assure that
  # statistics will be reported to the desired high-resolution limit
  # even if the observed unit cell differs slightly from the reference.
  ISO_ALLOWANCE = 0.1 # isomorphous recip cell volume changes no more than 10%
  params2.high_resolution = params.d_min / math.pow( (1.+ISO_ALLOWANCE),(1./3.) )
  params2.output.type = "real"
  if (params.include_bulk_solvent) :
    params2.fmodel.k_sol = params.k_sol
    params2.fmodel.b_sol = params.b_sol
  f_model = mmtbx.utils.fmodel_from_xray_structure(
    xray_structure = xray_structure,
    f_obs          = None,
    add_sigmas     = True,
    params         = params2).f_model
  if not params.merge_anomalous:
    f_model = f_model.generate_bijvoet_mates()
  i_model = f_model.as_intensity_array()

  return i_model

if (__name__ == "__main__") :
  run(sys.argv[1:])
