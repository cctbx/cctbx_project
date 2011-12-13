import mmtbx.command_line.fmodel
import mmtbx.utils
from iotbx import file_reader
import libtbx.phil.command_line
import sys

def run (params) :
  pdb_in = file_reader.any_file(params.model, force_type="pdb")
  pdb_in.assert_file_type("pdb")
  xray_structure = pdb_in.file_object.xray_structure_simple()
  xray_structure.show_summary()
  phil2 = mmtbx.command_line.fmodel.fmodel_from_xray_structure_master_params
  params2 = phil2.extract()
  params2.high_resolution = params.d_min
  params2.output.type = "real"
  if (params.include_bulk_solvent) :
    params2.fmodel.k_sol = params.k_sol
    params2.fmodel.b_sol = params.b_sol
  f_model = mmtbx.utils.fmodel_from_xray_structure(
    xray_structure = xray_structure,
    f_obs          = None,
    add_sigmas     = False,
    params         = params2).f_model.generate_bijvoet_mates()
  i_model = f_model.as_intensity_array()

  return i_model

if (__name__ == "__main__") :
  run(sys.argv[1:])
