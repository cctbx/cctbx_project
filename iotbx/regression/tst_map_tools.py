
from __future__ import absolute_import, division, print_function
from iotbx import file_reader
from iotbx import map_tools
import iotbx.pdb.hierarchy
from iotbx import mtz
from scitbx.array_family import flex
import os

def exercise_map_tools():
  prefix = "tmp_iotbx_map_tools"
  pdb_file = prefix + ".pdb"
  mtz_file = prefix + ".mtz"
  pdb_in = iotbx.pdb.input(source_info=None, lines="""\
ATOM      1  N   GLY P  -1     -22.866  -2.627  15.217  1.00  0.00           N
ATOM      2  CA  GLY P  -1     -22.714  -3.068  16.621  1.00  0.00           C
ATOM      3  C   GLY P  -1     -21.276  -3.457  16.936  1.00  0.00           C
ATOM      4  O   GLY P  -1     -20.538  -3.887  16.047  1.00  0.00           O
ATOM      5  H1  GLY P  -1     -22.583  -3.364  14.590  1.00  0.00           H
ATOM      6  H2  GLY P  -1     -22.293  -1.817  15.040  1.00  0.00           H
ATOM      7  H3  GLY P  -1     -23.828  -2.392  15.027  1.00  0.00           H
ATOM      8  HA2 GLY P  -1     -23.016  -2.261  17.288  1.00  0.00           H
ATOM      9  HA3 GLY P  -1     -23.352  -3.933  16.803  1.00  0.00           H
""")
  xrs = pdb_in.xray_structure_simple()
  hierarchy = pdb_in.construct_hierarchy()
  with open(pdb_file, "w") as f:
    f.write(hierarchy.as_pdb_string(xrs))
  fc = xrs.structure_factors(d_min=1.5).f_calc()
  dec = mtz.label_decorator(phases_prefix="PH")
  # part 1: phenix.refine style
  mtz_data = fc.as_mtz_dataset(
    column_root_label="2FOFCWT",
    label_decorator=dec)
  mtz_data.add_miller_array(fc,
    column_root_label="2FOFCWT_no_fill",
    label_decorator=dec)
  mtz_data.add_miller_array(fc,
    column_root_label="FOFCWT",
    label_decorator=dec)
  mtz_data.add_miller_array(fc,
    column_root_label="ANOM",
    label_decorator=dec)
  mtz_data.mtz_object().write(mtz_file)
  converted = map_tools.auto_convert_map_coefficients(
    mtz_file=mtz_file,
    pdb_file=pdb_file)
  assert (not None in [converted.f_map,converted.diff_map,converted.anom_map])
  assert (converted.f_map == "tmp_iotbx_map_tools_2mFo-DFc.ccp4")
  assert (converted.f_map_type == "2mFo-DFc")
  server = map_tools.server(mtz_file)
  files = server.convert_phenix_maps(file_base=prefix)
  files = [ os.path.basename(file_name) for file_name in files ]
  assert (files == ['tmp_iotbx_map_tools_anomalous.ccp4',
    'tmp_iotbx_map_tools_mFo-DFc.ccp4', 'tmp_iotbx_map_tools_2mFo-DFc.ccp4',
    'tmp_iotbx_map_tools_2mFo-DFc_no_fill.ccp4'])
  for fn in files :
    assert os.path.isfile(fn)
    os.remove(fn)
  file_name = server.convert_any_map(
    f_label="2FOFCWT,PH2FOFCWT",
    phi_label=None,
    fom_label=None,
    use_standard_map_names=True)
  assert (file_name == "tmp_iotbx_map_tools_2mFo-DFc.ccp4")
  assert os.path.isfile(file_name)
  # part 2: Phaser/Refmac style
  mtz_data = fc.as_mtz_dataset(
    column_root_label="FWT",
    label_decorator=mtz.ccp4_label_decorator())
  mtz_data.add_miller_array(fc,
    column_root_label="DELFWT",
    label_decorator=mtz.ccp4_label_decorator())
  mtz_data.mtz_object().write(mtz_file)
  converted = map_tools.auto_convert_map_coefficients(
    mtz_file=mtz_file,
    pdb_file=pdb_file)
  assert (not None in [converted.f_map, converted.diff_map])
  assert (converted.f_map == 'tmp_iotbx_map_tools_FWT.ccp4')
  assert (converted.f_map_type == "2mFo-DFc")
  server = map_tools.server(mtz_file)
  map_files = server.convert_ccp4_map(pdb_file=pdb_file)
  map_files = [ os.path.basename(file_name) for file_name in map_files ]
  assert (map_files == ['tmp_iotbx_map_tools_FWT.ccp4',
    'tmp_iotbx_map_tools_DELFWT.ccp4'])
  # part 3: Resolve style
  mtz_data = fc.as_mtz_dataset(
    column_root_label="FWT",
    label_decorator=mtz.ccp4_label_decorator())
  ampl = abs(fc)
  phases = fc.phases()
  fom = fc.customized_copy(
    data=flex.double(fc.data().size(), 1),
    observation_type=None)
  mtz_data.add_miller_array(ampl, column_root_label="FP")
  mtz_data.add_miller_array(phases, column_root_label="PHIM")
  mtz_data.add_miller_array(fom, column_root_label="FOMM")
  mtz_data.mtz_object().write(mtz_file)
  server = map_tools.server(mtz_file)
  map_file = server.convert_resolve_map(pdb_file=None,
    force=False)
  assert (map_file == 'tmp_iotbx_map_tools.ccp4')
  assert os.path.exists(map_file)
  # write_map_coefficients_generic
  map_tools.write_map_coefficients_generic(
    map_coeffs=[fc, fc.generate_bijvoet_mates(), fc, fc],
    map_types=["2mFo-DFc", "mFo-DFc", "anom", "other"],
    file_name="tmp_iotbx_map_tools2.mtz")
  mtz_in = file_reader.any_file("tmp_iotbx_map_tools2.mtz")
  labels = [a.info().label_string() for a in mtz_in.file_server.miller_arrays]
  assert (labels == ['2FOFCWT,PH2FOFCWT', 'FOFCWT,PHFOFCWT', 'ANOM,PHANOM',
                     'other,PHother'])

if (__name__ == "__main__"):
  exercise_map_tools()
  print("OK")
