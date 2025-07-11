from __future__ import absolute_import, division, print_function
from six.moves import cStringIO as StringIO
import os.path
import iotbx.pdb.hierarchy
from mmtbx.regression import model_1yjp
from iotbx.data_manager import DataManager
from libtbx.utils import null_out
import mmtbx.model
from mmtbx.building.extend_sidechains import master_params
import iotbx.pdb
from mmtbx.validation import rotalyze

from six.moves import zip

def exercise_model_only():
  from mmtbx.building import extend_sidechains
  import mmtbx.monomer_library
  pdb_in = iotbx.pdb.input(source_info=None, lines="""
ATOM     65  N   LYS A   7       6.033   4.704   1.582  1.00 17.49           N
ATOM     66  CA  LYS A   7       5.159   5.427   2.499  1.00 18.23           C
ATOM     67  C   LYS A   7       4.673   4.437   3.507  1.00 14.78           C
ATOM     68  O   LYS A   7       4.777   3.208   3.297  1.00 15.83           O
ATOM     69  CB  LYS A   7       3.959   6.057   1.760  1.00 23.56           C
ATOM     70  CG  LYS A   7       4.345   7.215   0.830  1.00 33.58           C
ATOM     71  CD  LYS A   7       3.213   7.570  -0.123  1.00 41.39           C
ATOM     72  CE  LYS A   7       2.976   6.471  -1.165  1.00 48.81           C
""")
  h = pdb_in.construct_hierarchy()
  extend_sidechains.extend_protein_model(
    pdb_hierarchy=h,
    mon_lib_srv=mmtbx.monomer_library.server.server())
  assert (h.as_pdb_string() == """\
ATOM      1  N   LYS A   7       6.033   4.704   1.582  1.00 17.49           N
ATOM      2  CA  LYS A   7       5.159   5.427   2.499  1.00 18.23           C
ATOM      3  C   LYS A   7       4.673   4.437   3.507  1.00 14.78           C
ATOM      4  O   LYS A   7       4.777   3.208   3.297  1.00 15.83           O
ATOM      5  CB  LYS A   7       3.980   6.048   1.757  1.00 23.56           C
ATOM      6  CG  LYS A   7       4.366   7.205   0.850  1.00 33.58           C
ATOM      7  CD  LYS A   7       3.241   7.559  -0.109  1.00 41.39           C
ATOM      8  CE  LYS A   7       3.011   6.457  -1.129  1.00 48.81           C
ATOM      9  NZ  LYS A   7       1.911   6.790  -2.075  0.10 26.71           N
TER
""")

def exercise_cmdline_cif():
  #
  params = iotbx.phil.parse(input_string = master_params).extract()
  #
  pdb_in = iotbx.pdb.input(source_info=None, lines=model_1yjp)
  m1 = mmtbx.model.manager(model_input = pdb_in, log = null_out())
  h_in = m1.get_hierarchy()
  #
  mtz_file = "tst_extend_sidechains.mtz"
  xrs = m1.get_xray_structure()
  f_calc = abs(xrs.structure_factors(d_min=1.5).f_calc())
  flags = f_calc.generate_r_free_flags(fraction=0.1)
  mtz = f_calc.as_mtz_dataset(column_root_label="F")
  mtz.add_miller_array(flags, column_root_label="FreeR_flag")
  mtz.mtz_object().write(mtz_file)
  sel_str = "not (resname TYR and not (name c or name o or name n or name oxt or name ca or name cb))"
  sel = m1.selection(sel_str)
  h_trimmed = h_in.select(sel)

  # Convert to mmcif_only:
  h_trimmed.rename_chain_id('A','AXLONG')

  pdb_file = "tst_extend_sidechains.cif"
  h_trimmed.write_pdb_or_mmcif_file(pdb_file)

  pdb_out = "tst_extend_sidechains_out.cif"
  prefix=os.path.splitext(os.path.basename(pdb_out))[0]

  dm = DataManager()
  m3 = dm.get_model(pdb_file)
  ma = dm.get_miller_arrays(filename = mtz_file)
  fmo3 = dm.get_fmodel(scattering_table="n_gaussian")
  out1 = StringIO()
  mmtbx.building.extend_sidechains.extend_and_refine(
    pdb_hierarchy=m3.get_hierarchy(),
    xray_structure=m3.get_xray_structure(),
    fmodel=fmo3,
    params=params,
    prefix=prefix,
    out=out1,
    output_model=pdb_out)

  assert ("1 sidechains extended." in out1.getvalue()), out1.getvalue()

  pdb_new = iotbx.pdb.input(file_name=pdb_out)
  h_new = pdb_new.construct_hierarchy()
  r1 = rotalyze.rotalyze(pdb_hierarchy=h_in, outliers_only=False)
  r2 = rotalyze.rotalyze(pdb_hierarchy=h_new, outliers_only=False)
  for o1, o2 in zip(r1.results, r2.results):
    assert o1.rotamer_name == o2.rotamer_name
  # Cleanup
  if os.path.isfile(pdb_out): os.remove(pdb_out)
  if os.path.isfile(pdb_file): os.remove(pdb_file)
  if os.path.isfile(prefix+'_maps.mtz'): os.remove(prefix+'_maps.mtz')
  #
  # Part 2: with sequence corrections
  #
  out2 = StringIO()
  #seq_file = "tst_extend_sidechains.fa"
  #with open(seq_file, "w") as f:
  #  f.write(">1yjp_new\nGNDQQNY")

  sequences = [iotbx.bioinformatics.sequence("GNDQQNY")]
  h_in_mod = h_trimmed.deep_copy()
  n_changed = mmtbx.building.extend_sidechains.correct_sequence(
    pdb_hierarchy=h_in_mod,
    sequences=sequences,
    out=out2)
  #print(h_in.as_sequence())
  #print(h_in_mod.as_sequence())

  pdb_file = "tst_extend_sidechains_2.cif"
  h_in_mod.write_pdb_or_mmcif_file(pdb_file)

  dm = DataManager()
  m4 = dm.get_model(pdb_file)
  ma = dm.get_miller_arrays(filename = mtz_file)
  fmo4 = dm.get_fmodel(scattering_table="n_gaussian")

  pdb_out = "tst_extend_sidechains_out2.cif"
  prefix=os.path.splitext(os.path.basename(pdb_out))[0]

  mmtbx.building.extend_sidechains.extend_and_refine(
    pdb_hierarchy=m4.get_hierarchy(),
    xray_structure=m4.get_xray_structure(),
    fmodel=fmo4,
    params=params,
    prefix=prefix,
    out=out2,
    output_model=pdb_out)

  assert ("2 sidechains extended." in out2.getvalue()), out2.getvalue()
  # Cleanup
  if os.path.isfile(pdb_out): os.remove(pdb_out)
  if os.path.isfile(pdb_file): os.remove(pdb_file)
  if os.path.isfile(prefix+'_maps.mtz'): os.remove(prefix+'_maps.mtz')
  if os.path.isfile(mtz_file): os.remove(mtz_file)

def exercise_cmdline():
  #
  params = iotbx.phil.parse(input_string = master_params).extract()
  #
  pdb_in = iotbx.pdb.input(source_info=None, lines=model_1yjp)
  m1 = mmtbx.model.manager(model_input = pdb_in, log = null_out())
  h_in = m1.get_hierarchy()
  #
  mtz_file = "tst_extend_sidechains.mtz"
  xrs = m1.get_xray_structure()
  f_calc = abs(xrs.structure_factors(d_min=1.5).f_calc())
  flags = f_calc.generate_r_free_flags(fraction=0.1)
  mtz = f_calc.as_mtz_dataset(column_root_label="F")
  mtz.add_miller_array(flags, column_root_label="FreeR_flag")
  mtz.mtz_object().write(mtz_file)
  sel_str = "not (resname TYR and not (name c or name o or name n or name oxt or name ca or name cb))"
  sel = m1.selection(sel_str)
  h_trimmed = h_in.select(sel)

  pdb_file = "tst_extend_sidechains.pdb"
  h_trimmed.write_pdb_file(pdb_file)

  pdb_out = "tst_extend_sidechains_out.pdb"
  prefix=os.path.splitext(os.path.basename(pdb_out))[0]

  dm = DataManager()
  m3 = dm.get_model(pdb_file)
  ma = dm.get_miller_arrays(filename = mtz_file)
  fmo3 = dm.get_fmodel(scattering_table="n_gaussian")
  out1 = StringIO()
  mmtbx.building.extend_sidechains.extend_and_refine(
    pdb_hierarchy=m3.get_hierarchy(),
    xray_structure=m3.get_xray_structure(),
    fmodel=fmo3,
    params=params,
    prefix=prefix,
    out=out1,
    output_model=pdb_out)

  assert ("1 sidechains extended." in out1.getvalue()), out1.getvalue()

  pdb_new = iotbx.pdb.input(file_name=pdb_out)
  h_new = pdb_new.construct_hierarchy()
  r1 = rotalyze.rotalyze(pdb_hierarchy=h_in, outliers_only=False)
  r2 = rotalyze.rotalyze(pdb_hierarchy=h_new, outliers_only=False)
  for o1, o2 in zip(r1.results, r2.results):
    assert o1.rotamer_name == o2.rotamer_name
  # Cleanup
  if os.path.isfile(pdb_out): os.remove(pdb_out)
  if os.path.isfile(pdb_file): os.remove(pdb_file)
  if os.path.isfile(prefix+'_maps.mtz'): os.remove(prefix+'_maps.mtz')
  #
  # Part 2: with sequence corrections
  #
  out2 = StringIO()
  #seq_file = "tst_extend_sidechains.fa"
  #with open(seq_file, "w") as f:
  #  f.write(">1yjp_new\nGNDQQNY")

  sequences = [iotbx.bioinformatics.sequence("GNDQQNY")]
  h_in_mod = h_trimmed.deep_copy()
  n_changed = mmtbx.building.extend_sidechains.correct_sequence(
    pdb_hierarchy=h_in_mod,
    sequences=sequences,
    out=out2)
  #print(h_in.as_sequence())
  #print(h_in_mod.as_sequence())

  pdb_file = "tst_extend_sidechains_2.pdb"
  h_in_mod.write_pdb_file(pdb_file)

  dm = DataManager()
  m4 = dm.get_model(pdb_file)
  ma = dm.get_miller_arrays(filename = mtz_file)
  fmo4 = dm.get_fmodel(scattering_table="n_gaussian")

  pdb_out = "tst_extend_sidechains_out2.pdb"
  prefix=os.path.splitext(os.path.basename(pdb_out))[0]

  mmtbx.building.extend_sidechains.extend_and_refine(
    pdb_hierarchy=m4.get_hierarchy(),
    xray_structure=m4.get_xray_structure(),
    fmodel=fmo4,
    params=params,
    prefix=prefix,
    out=out2,
    output_model=pdb_out)

  assert ("2 sidechains extended." in out2.getvalue()), out2.getvalue()
  # Cleanup
  if os.path.isfile(pdb_out): os.remove(pdb_out)
  if os.path.isfile(pdb_file): os.remove(pdb_file)
  if os.path.isfile(prefix+'_maps.mtz'): os.remove(prefix+'_maps.mtz')
  if os.path.isfile(mtz_file): os.remove(mtz_file)

def exercise_correct_sequence():
  from mmtbx.building import extend_sidechains
  from mmtbx.regression import model_1yjp
  import iotbx.bioinformatics
  hierarchy = iotbx.pdb.input(source_info=None, lines=model_1yjp).construct_hierarchy()
  sequences = [ iotbx.bioinformatics.sequence("GNDQQNY") ]
  out = StringIO()
  n_changed = extend_sidechains.correct_sequence(
    pdb_hierarchy=hierarchy.deep_copy(),
    sequences=sequences,
    out=out)
  assert (n_changed == 1)
  assert ("  chain 'A'    3  ASP --> ASP (1 atoms removed)" in out.getvalue())
  out = StringIO()
  n_changed = extend_sidechains.correct_sequence(
    pdb_hierarchy=hierarchy.deep_copy(),
    sequences=sequences,
    truncate_to_cbeta=True,
    out=out)
  assert (n_changed == 1)
  assert ("  chain 'A'    3  ASP --> ASP (3 atoms removed)" in out.getvalue())

if (__name__ == "__main__"):
  exercise_cmdline_cif()
  exercise_model_only()
  exercise_correct_sequence()
  exercise_cmdline()
  print("OK")
