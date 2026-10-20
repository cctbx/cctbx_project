"""reduce2 with qmready=True calculates the total charge (as qr.finalise does
with calculate_charge=True) and appends it to the default output file name:
model.pdb -> modelH_charge-2.pdb.
"""
from __future__ import absolute_import, division, print_function
import os
from libtbx.utils import null_out
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2

# Asp-Glu: Asp -1, Glu -1, N-terminal NH3+ +1, C-terminal COO- -1 -> -2
pdb_str_neg = """\
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
ATOM      1  N   ASP A   1       5.000   5.000   5.000  1.00 10.00           N
ATOM      2  CA  ASP A   1       6.458   5.000   5.000  1.00 10.00           C
ATOM      3  C   ASP A   1       7.009   6.420   5.000  1.00 10.00           C
ATOM      4  O   ASP A   1       6.251   7.390   5.000  1.00 10.00           O
ATOM      5  CB  ASP A   1       6.994   4.250   3.780  1.00 10.00           C
ATOM      6  CG  ASP A   1       8.510   4.250   3.780  1.00 10.00           C
ATOM      7  OD1 ASP A   1       9.120   4.900   4.650  1.00 10.00           O
ATOM      8  OD2 ASP A   1       9.080   3.580   2.890  1.00 10.00           O
ATOM      9  N   GLU A   2       8.330   6.540   5.000  1.00 10.00           N
ATOM     10  CA  GLU A   2       8.960   7.850   5.000  1.00 10.00           C
ATOM     11  C   GLU A   2      10.470   7.700   5.000  1.00 10.00           C
ATOM     12  O   GLU A   2      11.000   6.600   5.000  1.00 10.00           O
ATOM     13  CB  GLU A   2       8.520   8.640   6.240  1.00 10.00           C
ATOM     14  CG  GLU A   2       8.950  10.100   6.250  1.00 10.00           C
ATOM     15  CD  GLU A   2       8.500  10.850   7.490  1.00 10.00           C
ATOM     16  OE1 GLU A   2       7.910  10.220   8.390  1.00 10.00           O
ATOM     17  OE2 GLU A   2       8.740  12.070   7.560  1.00 10.00           O
END
"""

# Ala-Ser-Lys: Lys +1, N-terminal NH3+ +1, C-terminal COO- -1 -> +1
pdb_str_pos = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
ATOM      1  N   ALA A   1      11.204   2.315  10.312  1.00 20.00           N
ATOM      2  CA  ALA A   1      12.555   2.835  10.528  1.00 20.00           C
ATOM      3  C   ALA A   1      13.588   1.726  10.381  1.00 20.00           C
ATOM      4  O   ALA A   1      13.325   0.567  10.721  1.00 20.00           O
ATOM      5  CB  ALA A   1      12.866   3.965   9.556  1.00 20.00           C
ATOM      6  N   SER A   2      14.767   2.090   9.884  1.00 20.00           N
ATOM      7  CA  SER A   2      15.853   1.130   9.701  1.00 20.00           C
ATOM      8  C   SER A   2      16.968   1.370  10.713  1.00 20.00           C
ATOM      9  O   SER A   2      17.281   2.516  11.041  1.00 20.00           O
ATOM     10  CB  SER A   2      16.411   1.224   8.278  1.00 20.00           C
ATOM     11  OG  SER A   2      15.456   0.826   7.311  1.00 20.00           O
ATOM     12  N   LYS A   3      17.560   0.282  11.199  1.00 20.00           N
ATOM     13  CA  LYS A   3      18.643   0.370  12.172  1.00 20.00           C
ATOM     14  C   LYS A   3      19.945  -0.148  11.570  1.00 20.00           C
ATOM     15  O   LYS A   3      19.964  -1.176  10.892  1.00 20.00           O
ATOM     16  CB  LYS A   3      18.290  -0.424  13.432  1.00 20.00           C
ATOM     17  CG  LYS A   3      19.308  -0.301  14.556  1.00 20.00           C
ATOM     18  CD  LYS A   3      18.913  -1.137  15.762  1.00 20.00           C
ATOM     19  CE  LYS A   3      19.919  -1.002  16.894  1.00 20.00           C
ATOM     20  NZ  LYS A   3      19.531  -1.812  18.059  1.00 20.00           N
ATOM     21  OXT LYS A   3      20.995   0.451  11.763  1.00 20.00           O
END
"""

# Peptide fragments from 3vt9 (143-144, 147, 150, 232-233): NH3+ and COO- at
# each fragment, Asp -1 -> -1
pdb_str_frag = """\
CRYST1   27.973   30.241   36.301  90.00  90.00  90.00 P 1
ATOM      1  N   TYR A 143      18.189  17.846   5.000  1.00 23.64           N
ATOM      2  CA  TYR A 143      18.357  18.668   6.193  1.00 24.34           C
ATOM      3  C   TYR A 143      18.140  20.148   5.894  1.00 25.18           C
ATOM      4  O   TYR A 143      18.982  20.801   5.278  1.00 26.15           O
ATOM      5  CB  TYR A 143      19.753  18.468   6.795  1.00 23.00           C
ATOM      6  CG  TYR A 143      19.906  19.085   8.170  1.00 26.16           C
ATOM      7  CD1 TYR A 143      19.546  18.377   9.315  1.00 26.15           C
ATOM      8  CD2 TYR A 143      20.374  20.392   8.326  1.00 29.44           C
ATOM      9  CE1 TYR A 143      19.644  18.953  10.580  1.00 23.54           C
ATOM     10  CE2 TYR A 143      20.475  20.979   9.587  1.00 19.77           C
ATOM     11  CZ  TYR A 143      20.108  20.254  10.711  1.00 29.80           C
ATOM     12  OH  TYR A 143      20.202  20.824  11.966  1.00 26.06           O
ATOM     13  N   ASP A 144      17.005  20.671   6.345  1.00 26.39           N
ATOM     14  CA  ASP A 144      16.656  22.072   6.145  1.00 26.41           C
ATOM     15  C   ASP A 144      16.935  22.856   7.424  1.00 27.66           C
ATOM     16  O   ASP A 144      16.168  22.794   8.387  1.00 27.41           O
ATOM     17  CB  ASP A 144      15.177  22.187   5.774  1.00 26.75           C
ATOM     18  CG  ASP A 144      14.663  23.604   5.876  1.00 32.49           C
ATOM     19  OD1 ASP A 144      15.387  24.533   5.455  1.00 28.31           O
ATOM     20  OD2 ASP A 144      13.534  23.787   6.382  1.00 41.39           O
ATOM     21  N   TYR A 147      14.226  25.192   9.139  1.00 35.91           N
ATOM     22  CA  TYR A 147      13.181  24.667  10.017  1.00 36.98           C
ATOM     23  C   TYR A 147      11.773  24.971   9.506  1.00 37.02           C
ATOM     24  O   TYR A 147      10.821  25.037  10.285  1.00 36.71           O
ATOM     25  CB  TYR A 147      13.343  25.241  11.430  1.00 40.78           C
ATOM     26  CG  TYR A 147      14.497  24.655  12.217  1.00 42.98           C
ATOM     27  CD1 TYR A 147      15.745  24.462  11.627  1.00 44.97           C
ATOM     28  CD2 TYR A 147      14.340  24.301  13.555  1.00 50.93           C
ATOM     29  CE1 TYR A 147      16.806  23.931  12.347  1.00 50.90           C
ATOM     30  CE2 TYR A 147      15.396  23.769  14.285  1.00 54.15           C
ATOM     31  CZ  TYR A 147      16.625  23.587  13.674  1.00 54.34           C
ATOM     32  OH  TYR A 147      17.672  23.064  14.392  1.00 56.08           O
ATOM     33  N   PHE A 150       8.896  22.695  10.613  1.00 32.01           N
ATOM     34  CA  PHE A 150       8.615  22.605  12.043  1.00 33.65           C
ATOM     35  C   PHE A 150       7.281  23.232  12.427  1.00 34.16           C
ATOM     36  O   PHE A 150       6.780  23.014  13.527  1.00 33.81           O
ATOM     37  CB  PHE A 150       9.737  23.269  12.847  1.00 36.30           C
ATOM     38  CG  PHE A 150      10.953  22.404  13.022  1.00 35.86           C
ATOM     39  CD1 PHE A 150      11.504  21.724  11.937  1.00 36.69           C
ATOM     40  CD2 PHE A 150      11.548  22.267  14.271  1.00 32.00           C
ATOM     41  CE1 PHE A 150      12.628  20.920  12.095  1.00 37.39           C
ATOM     42  CE2 PHE A 150      12.672  21.466  14.442  1.00 32.96           C
ATOM     43  CZ  PHE A 150      13.214  20.791  13.352  1.00 37.79           C
ATOM     59  N   TYR A 232       9.572  13.044  12.058  1.00 22.04           N
ATOM     60  CA  TYR A 232      10.003  13.787  10.879  1.00 19.32           C
ATOM     61  C   TYR A 232      11.468  13.494  10.564  1.00 19.71           C
ATOM     62  O   TYR A 232      11.830  13.220   9.418  1.00 18.34           O
ATOM     63  CB  TYR A 232       9.805  15.286  11.125  1.00 16.65           C
ATOM     64  CG  TYR A 232      10.406  16.197  10.076  1.00 17.99           C
ATOM     65  CD1 TYR A 232       9.680  16.574   8.949  1.00 22.93           C
ATOM     66  CD2 TYR A 232      11.691  16.710  10.231  1.00 23.97           C
ATOM     67  CE1 TYR A 232      10.220  17.444   8.000  1.00 24.15           C
ATOM     68  CE2 TYR A 232      12.241  17.578   9.293  1.00 21.34           C
ATOM     69  CZ  TYR A 232      11.504  17.942   8.183  1.00 27.57           C
ATOM     70  OH  TYR A 232      12.055  18.799   7.259  1.00 25.99           O
ATOM     71  N   SER A 233      12.301  13.542  11.598  1.00 19.07           N
ATOM     72  CA  SER A 233      13.732  13.341  11.442  1.00 19.75           C
ATOM     73  C   SER A 233      14.047  11.941  10.931  1.00 19.48           C
ATOM     74  O   SER A 233      14.990  11.756  10.163  1.00 17.69           O
ATOM     75  CB  SER A 233      14.449  13.574  12.776  1.00 20.71           C
ATOM     76  OG  SER A 233      14.283  14.908  13.217  1.00 13.41           O
END
"""

def run(model_text, prefix, expected_files, in_suffix=".pdb", extra_args=()):
  fn_in = prefix + in_suffix
  with open(fn_in, "w") as f:
    f.write(model_text)
  for fn in expected_files:
    if os.path.exists(fn): os.remove(fn)
  results = run_program(program_class=reduce2.Program, logger=null_out(),
    args=[fn_in, "overwrite=True"] + list(extra_args))
  for fn in expected_files:
    assert os.path.isfile(fn), fn
  return results

def test_qmready_negative_charge():
  prefix = "tst_reduce2_qmready_charge_1"
  r = run(pdb_str_neg, prefix,
    expected_files=[prefix+"H_charge-2.pdb", prefix+"H_charge-2.txt"],
    extra_args=["qmready=True"])
  assert r.total_charge == -2, r.total_charge

def test_qmready_positive_charge():
  # Same format as qr.finalise: no '+' sign
  prefix = "tst_reduce2_qmready_charge_2"
  r = run(pdb_str_pos, prefix,
    expected_files=[prefix+"H_charge1.pdb", prefix+"H_charge1.txt"],
    extra_args=["qmready=True"])
  assert r.total_charge == 1, r.total_charge

def test_qmready_cif():
  # mmCIF in -> mmCIF out, default name with charge
  import iotbx.pdb
  prefix = "tst_reduce2_qmready_charge_3"
  inp = iotbx.pdb.input(lines=pdb_str_neg.split("\n"), source_info=None)
  cif_str = inp.construct_hierarchy().as_mmcif_string(
    crystal_symmetry=inp.crystal_symmetry())
  r = run(cif_str, prefix,
    expected_files=[prefix+"H_charge-2.cif", prefix+"H_charge-2.txt"],
    in_suffix=".cif", extra_args=["qmready=True"])
  assert r.total_charge == -2, r.total_charge

def test_qmready_explicit_filename():
  prefix = "tst_reduce2_qmready_charge_4"
  fn_out = prefix + "_custom.pdb"
  r = run(pdb_str_neg, prefix, expected_files=[fn_out],
    extra_args=["qmready=True", "output.filename="+fn_out])
  assert r.total_charge == -2, r.total_charge
  assert not os.path.exists(prefix + "_custom_charge-2.pdb")

def test_no_qmready_no_charge():
  prefix = "tst_reduce2_qmready_charge_5"
  r = run(pdb_str_neg, prefix,
    expected_files=[prefix+"H.pdb", prefix+"H.txt"])
  assert r.total_charge is None, r.total_charge
  assert not os.path.exists(prefix+"H_charge-2.pdb")

def test_qmready_chain_breaks():
  # NH3 at each fragment start and OXT at each fragment end, not only at the
  # first and last residue of the chain
  import iotbx.pdb
  prefix = "tst_reduce2_qmready_charge_6"
  fn_out = prefix+"H_charge-1.pdb"
  r = run(pdb_str_frag, prefix,
    expected_files=[fn_out, prefix+"H_charge-1.txt"],
    extra_args=["qmready=True"])
  assert r.total_charge == -1, r.total_charge
  nh3, oxt = [], []
  for rg in iotbx.pdb.input(fn_out).construct_hierarchy().residue_groups():
    names = [a.name.strip() for a in rg.atoms()]
    if 'H1' in names and 'H2' in names and 'H3' in names:
      nh3.append(rg.resseq.strip())
    if 'OXT' in names:
      oxt.append(rg.resseq.strip())
  assert nh3 == ['143', '147', '150', '232'], nh3
  assert oxt == ['144', '147', '150', '233'], oxt

def test_qmready_no_cryst1():
  # OXT added to a model without symmetry: reduce2 crashed boxing it
  import iotbx.pdb
  prefix = "tst_reduce2_qmready_charge_7"
  fn_out = prefix+"H_charge-2.pdb"
  pdb_str = "".join(l + "\n" for l in pdb_str_neg.splitlines()
    if not l.startswith("CRYST1"))
  r = run(pdb_str, prefix, expected_files=[fn_out, prefix+"H_charge-2.txt"],
    extra_args=["qmready=True"])
  assert r.total_charge == -2, r.total_charge
  inp = iotbx.pdb.input(fn_out)
  assert inp.crystal_symmetry() is None
  names = [a.name.strip() for a in inp.construct_hierarchy().atoms()]
  assert 'OXT' in names

def ccd_ligand_pdb(code, offset=(-5.0, 20.0, 30.0)):
  """heavy atoms of a CCD entry at its ideal coordinates, chain C"""
  from mmtbx.chemical_components import get_cif_dictionary
  lines = []
  for i, a in enumerate(get_cif_dictionary(code)["_chem_comp_atom"]):
    if a.type_symbol == "H": continue
    xyz = [float(v) + o for v, o in zip((a.pdbx_model_Cartn_x_ideal,
      a.pdbx_model_Cartn_y_ideal, a.pdbx_model_Cartn_z_ideal), offset)]
    name = a.atom_id if len(a.atom_id) == 4 else " " + a.atom_id
    lines.append("HETATM%5d %-4s %3s C   1    %8.3f%8.3f%8.3f  1.00 20.00"
      "          %2s" % (900 + i, name, code, xyz[0], xyz[1], xyz[2],
      a.type_symbol))
  return "\n".join(lines) + "\n"

def test_qmready_ccd_only_ligand():
  # BGM and 6EM: in neither the monomer library nor GeoStd. The charge step
  # had only the user CIFs, not reduce2's CCD objects, and stopped: "Number
  # of atoms with unknown nonbonded energy type symbols". 1yjp is neutral
  # (NH3+, COO-); BGM: 0, all 13 CCD H incl. HOP2/HOP3; 6EM: +1 on N4
  from mmtbx.regression.reduce2_neutron_utils import pdb_1yjp
  for i, code, n_h, charge in ((8, "BGM", 13, 0), (9, "6EM", 17, 1)):
    prefix = "tst_reduce2_qmready_charge_%d" % i
    model_text = pdb_1yjp.replace("END\n", "") + ccd_ligand_pdb(code) + "END\n"
    r = run(model_text, prefix, expected_files=[
      prefix + "H_charge%d.pdb" % charge, prefix + "H_charge%d.txt" % charge],
      extra_args=["qmready=True"])
    assert r.total_charge == charge, (code, r.total_charge)
    h = [a for a in r.model.get_hierarchy().atoms()
         if a.parent().resname.strip() == code and a.element_is_hydrogen()]
    assert len(h) == n_h, (code, len(h))

if __name__ == "__main__":
  test_qmready_negative_charge()
  test_qmready_positive_charge()
  test_qmready_cif()
  test_qmready_explicit_filename()
  test_no_qmready_no_charge()
  test_qmready_chain_breaks()
  test_qmready_no_cryst1()
  test_qmready_ccd_only_ligand()
  print("OK")
