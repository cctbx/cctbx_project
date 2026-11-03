"""reduce2 model_id is the MODEL number (PDB MODEL record, mmCIF
pdbx_PDB_model_num), checked before the run: a number that matches no model,
or one <= 0, stops with a Sorry listing the available numbers. A single model
without a number is MODEL 1. It used to be a 0-based position, falling back
to the first model when out of range. The selected model keeps the user's
restraint CIFs. Models: 1yjp Gly 1 - Asn 2, and ligand 04 that only its CIF
describes; each further model is shifted by 5 A in x.
"""
from __future__ import absolute_import, division, print_function
import os
import iotbx.pdb
from libtbx.utils import null_out, Sorry
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2
from mmtbx.regression.reduce2_neutron_utils import heavy_pdb, zq4_atoms, zq4_cif

pdb_str = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
ATOM      1  N   GLY A   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      3  C   GLY A   1      -8.015   3.140   4.419  1.00 16.16           C
ATOM      4  O   GLY A   1      -7.523   2.521   5.381  1.00 16.78           O
ATOM      5  N   ASN A   2      -7.656   2.923   3.155  1.00 15.02           N
ATOM      6  CA  ASN A   2      -6.522   2.038   2.831  1.00 14.10           C
ATOM      7  C   ASN A   2      -5.241   2.537   3.427  1.00 13.13           C
ATOM      8  O   ASN A   2      -4.978   3.742   3.426  1.00 11.91           O
ATOM      9  CB  ASN A   2      -6.346   1.881   1.341  1.00 15.38           C
ATOM     10  CG  ASN A   2      -7.584   1.342   0.692  1.00 14.08           C
ATOM     11  OD1 ASN A   2      -8.025   0.227   1.016  1.00 17.46           O
ATOM     12  ND2 ASN A   2      -8.204   2.155  -0.169  1.00 11.72           N
END
"""

x_first = -9.009   # N GLY 1 in the first model

def make_input(prefix, numbers, cif=False):
  """input file with one model per MODEL number (None: no MODEL record)"""
  inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  m1 = inp.construct_hierarchy().models()[0]
  root = iotbx.pdb.hierarchy.root()
  for k, number in enumerate(numbers):
    m = m1.detached_copy()
    m.id = "" if number is None else str(number)
    for a in m.atoms():
      x, y, z = a.xyz
      a.xyz = (x + 5.0 * k, y, z)
    root.append_model(m)
  cs = inp.crystal_symmetry()
  fn = prefix + (".cif" if cif else ".pdb")
  with open(fn, "w") as f:
    if cif: f.write(root.as_mmcif_string(crystal_symmetry=cs))
    else:   f.write(root.as_pdb_string(crystal_symmetry=cs))
  return fn

def run(fn_in, prefix, model_id):
  """-> (number of models written, shift of the first model's N GLY 1)"""
  fn_out = prefix + "_out.pdb"
  if os.path.exists(fn_out): os.remove(fn_out)
  run_program(program_class=reduce2.Program, logger=null_out(),
    args=[fn_in, "model_id=%d" % model_id, "output.filename=" + fn_out,
          "output.description_file_name=" + prefix + "_out.txt",
          "overwrite=True"])
  models = iotbx.pdb.input(fn_out).construct_hierarchy().models()
  n = [a for a in models[0].atoms() if a.name.strip() == "N" and
       a.parent().resname == "GLY"][0]
  return len(models), round(n.xyz[0] - x_first, 3)

def check_sorry(fn_in, prefix, model_id, available):
  try:
    run(fn_in, prefix, model_id)
  except Sorry as e:
    expected = "model_id=%d matches no model. Available MODEL numbers: %s" % (
      model_id, available)
    assert expected in str(e), str(e)
    assert not os.path.exists(prefix + "_out.pdb")
  else:
    raise AssertionError("model_id=%d: no Sorry" % model_id)

def test_two_models():
  fn = make_input("tst_reduce2_model_id_1", [1, 2])
  assert run(fn, "tst_reduce2_model_id_1_1", 1) == (1, 0.0)
  assert run(fn, "tst_reduce2_model_id_1_2", 2) == (1, 5.0)
  for model_id in (3, 0, -1):
    check_sorry(fn, "tst_reduce2_model_id_1_x", model_id, "1, 2")

def test_non_sequential_numbers():
  fn = make_input("tst_reduce2_model_id_2", [3, 7])
  assert run(fn, "tst_reduce2_model_id_2_7", 7) == (1, 5.0)
  check_sorry(fn, "tst_reduce2_model_id_2_x", 1, "3, 7")

def test_single_model_without_number():
  fn = make_input("tst_reduce2_model_id_3", [None])
  with open(fn) as f: assert "MODEL" not in f.read()
  assert run(fn, "tst_reduce2_model_id_3_1", 1) == (1, 0.0)
  check_sorry(fn, "tst_reduce2_model_id_3_x", 2, "1")

def test_mmcif_model_numbers():
  fn = make_input("tst_reduce2_model_id_4", [1, 2], cif=True)
  assert run(fn, "tst_reduce2_model_id_4_2", 2) == (1, 5.0)
  assert run(fn, "tst_reduce2_model_id_4_1", 1) == (1, 0.0)

def make_ligand_input(prefix, numbers):
  """ligand 04 and its CIF; one model per MODEL number (None: no MODEL)"""
  src = heavy_pdb("04", zq4_atoms)
  cryst1 = src.splitlines()[0]
  atoms = [l for l in src.splitlines() if l.startswith("HETATM")]
  lines = [cryst1]
  for k, number in enumerate(numbers):
    if number is not None: lines.append("MODEL %8d" % number)
    lines += [l[:30] + "%8.3f" % (float(l[30:38]) + 5.0 * k) + l[38:]
              for l in atoms]
    if number is not None: lines.append("ENDMDL")
  with open(prefix + ".pdb", "w") as f: f.write("\n".join(lines) + "\nEND\n")
  with open(prefix + "_04.cif", "w") as f: f.write(zq4_cif)
  return prefix + ".pdb", prefix + "_04.cif"

zq4_h = sorted(x[0] for x in zq4_atoms if x[1] == "H")
zq4_c1 = [x[3] for x in zq4_atoms if x[0] == "C1"][0] + 10.0  # heavy_pdb offset

def run_ligand(fn_in, fn_cif, prefix, model_id):
  """-> [(H names, shift of C1) per model written]"""
  fn_out = prefix + "_out.pdb"
  if os.path.exists(fn_out): os.remove(fn_out)
  run_program(program_class=reduce2.Program, logger=null_out(),
    args=[fn_in, fn_cif, "model_id=%d" % model_id, "output.filename=" + fn_out,
          "output.description_file_name=" + prefix + "_out.txt",
          "overwrite=True"])
  result = []
  for m in iotbx.pdb.input(fn_out).construct_hierarchy().models():
    c1 = [a for a in m.atoms() if a.name.strip() == "C1"][0]
    result.append((sorted(a.name.strip() for a in m.atoms()
                          if a.element.strip() == "H"),
                   round(c1.xyz[0] - zq4_c1, 3)))
  return result

def test_user_cif_single_model():
  # used to stop: Restraints were not found ... 04
  fn, cif = make_ligand_input("tst_reduce2_model_id_5", [None])
  assert run_ligand(fn, cif, "tst_reduce2_model_id_5_1", 1) == [
    (zq4_h, 0.0)]

def test_user_cif_two_models():
  fn, cif = make_ligand_input("tst_reduce2_model_id_6", [1, 2])
  assert run_ligand(fn, cif, "tst_reduce2_model_id_6_2", 2) == [
    (zq4_h, 5.0)]

if __name__ == "__main__":
  test_two_models()
  test_non_sequential_numbers()
  test_single_model_without_number()
  test_mmcif_model_numbers()
  test_user_cif_single_model()
  test_user_cif_two_models()
  print("OK")
