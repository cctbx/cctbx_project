
from __future__ import division
from builtins import object
from libtbx.utils import null_out
import os

def exercise_1 () :
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
  hierarchy = pdb_in.file_object.hierarchy
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

def exercise_2():
  import iotbx.pdb
  import mmtbx.maps.utils
  pdb_str1="""
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1
HETATM    1  C    C      1       4.271   0.000   0.000  1.00  5.00           C
HETATM    2  ?    ?      2       5.729   0.000   0.000  1.00  5.00           ?
HETATM    1  X    D      1       4.271   0.000   0.000  1.00  5.00           X
HETATM    1  Z    E      1       4.271   0.000   0.000  1.00  5.00           Z
END
"""
  pdb_str2="""
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1
HETATM    1  C    C      1       4.271   0.000   0.000  1.00  5.00           C
END
"""
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str1)
  pdb_inp.write_pdb_file(file_name = "tst_exercise_2_map_utils.pdb")
  fc = iotbx.pdb.input(source_info=None,
    lines=pdb_str2).xray_structure_simple().structure_factors(d_min=2).f_calc()
  class dummy(object):
    def amplitudes(self): return "2FOFCWT"
    def phases(self,root_label=None): return "PH2FOFCWT"
  mtz_dataset = fc.as_mtz_dataset(column_root_label=dummy().amplitudes(),
    label_decorator=dummy())
  mtz_dataset.add_miller_array(miller_array=abs(fc), column_root_label="FOBS_X")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = "tst_exercise_2_map_utils.mtz")
  mfn1 = "tst_exercise_2_map_utils_output.mtz"
  mmtbx.maps.utils.create_map_from_pdb_and_mtz(
    pdb_file="tst_exercise_2_map_utils.pdb",
    mtz_file="tst_exercise_2_map_utils.mtz",
    output_file=mfn1,
    out=null_out())

if (__name__ == "__main__") :
  exercise_1()
  exercise_2()
