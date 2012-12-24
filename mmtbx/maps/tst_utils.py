
from __future__ import division
from libtbx.utils import null_out
import os

def exercise () :
  from mmtbx.regression import make_fake_anomalous_data
  import mmtbx.maps.utils
  import mmtbx.utils
  from iotbx.file_reader import any_file
  mtz_file, pdb_file = make_fake_anomalous_data.generate_calcium_inputs(
    "tst_map_utils")
  # create_map_from_pdb_and_mtz
  mfn1 = "tst_map_utils_1.mtz"
  mmtbx.maps.utils.create_map_from_pdb_and_mtz(
    pdb_file=pdb_file,
    mtz_file=mtz_file,
    output_file=mfn1,
    fill=False,
    out=null_out())
  assert (os.path.isfile(mfn1))
  maps_in = any_file(mfn1)
  assert (len(maps_in.file_server.miller_arrays) == 3)
  # generate_water_omit_map
  pdb_in = any_file(pdb_file)
  hierarchy = pdb_in.file_object.construct_hierarchy()
  xrs = pdb_in.file_object.xray_structure_simple()
  mtz_in = any_file(mtz_file)
  f_obs = mtz_in.file_server.miller_arrays[0]
  flags = mtz_in.file_server.miller_arrays[1]
  flags = flags.customized_copy(data=flags.data()==1)
  fmodel = mmtbx.utils.fmodel_simple(
    f_obs=f_obs,
    r_free_flags=flags,
    scattering_table="n_gaussian",
    xray_structures=[xrs],
    bulk_solvent_correction=False,
    skip_twin_detection=True)
  omit = mmtbx.maps.utils.generate_water_omit_map(
    fmodel=fmodel,
    pdb_hierarchy=hierarchy,
    log=null_out())
  assert (omit.n_waters == 3)
  mfn2 = "tst_map_utils_2.mtz"
  pdbfn2 = "tst_map_utils_2.pdb"
  omit.write_map_coeffs(mfn2)
  omit.write_pdb_file(pdbfn2)
  for fn in [mfn1, mfn2, pdbfn2] :
    os.remove(fn)

if (__name__ == "__main__") :
  exercise()
