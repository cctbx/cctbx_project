
"""
Test for mmtbx.ncs.ligands module and associated command-line tool
mmtbx.apply_ncs_to_ligands (used in phenix.ligand_pipeline)
"""

from __future__ import absolute_import, division, print_function
from libtbx import easy_run
from libtbx.utils import null_out
import os.path as op
import os

def make_inputs(prefix="tst_ligand_ncs"):
  pdb_raw = """\
CRYST1   25.000   10.000   25.477  90.00 107.08  90.00 P 1
ATOM      1  N   GLY A   1       9.010   0.535  -6.096  1.00 16.23           N
ATOM      2  CA  GLY A   1       9.052   0.053  -4.707  1.00 16.20           C
ATOM      3  C   GLY A   1       8.015  -1.161  -4.448  1.00 15.91           C
ATOM      4  O   GLY A   1       7.563  -1.905  -5.385  1.00 16.11           O
ATOM      5  N   ASN A   2       7.642  -1.374  -3.185  1.00 15.02           N
ATOM      6  CA  ASN A   2       6.537  -2.383  -2.848  1.00 13.92           C
ATOM      7  C   ASN A   2       5.231  -1.825  -3.428  1.00 12.24           C
ATOM      8  O   ASN A   2       4.965  -0.471  -3.432  1.00 11.42           O
ATOM      9  CB  ASN A   2       6.391  -2.550  -1.340  1.00 14.42           C
ATOM     10  CG  ASN A   2       7.639  -3.142  -0.697  1.00 13.92           C
ATOM     11  ND2 ASN A   2       8.253  -2.228   0.168  1.00 12.78           N
ATOM     12  OD1 ASN A   2       8.062  -4.416  -0.978  1.00 14.39           O
ATOM     13  N   ASN A   3       4.414  -2.880  -3.904  1.00 12.20           N
ATOM     14  CA  ASN A   3       3.168  -2.511  -4.594  1.00 11.81           C
ATOM     15  C   ASN A   3       1.927  -3.147  -3.910  1.00 10.59           C
ATOM     16  O   ASN A   3       1.840  -4.501  -3.653  1.00 10.32           O
ATOM     17  CB  ASN A   3       3.235  -3.076  -6.046  1.00 12.51           C
ATOM     18  CG  ASN A   3       1.973  -2.713  -6.848  1.00 12.59           C
ATOM     19  ND2 ASN A   3       1.258  -3.883  -7.268  1.00 12.39           N
ATOM     20  OD1 ASN A   3       1.661  -1.397  -7.106  1.00 13.64           O
ATOM     21  N   GLN A   4       0.973  -2.158  -3.608  1.00 10.34           N
ATOM     22  CA  GLN A   4      -0.369  -2.634  -3.208  1.00 10.00           C
ATOM     23  C   GLN A   4      -1.409  -1.846  -4.085  1.00 10.21           C
ATOM     24  O   GLN A   4      -1.521  -0.467  -4.070  1.00  8.99           O
ATOM     25  CB  GLN A   4      -0.659  -2.286  -1.736  1.00 10.00           C
ATOM     26  CG  GLN A   4      -1.933  -3.051  -1.229  1.00 10.50           C
ATOM     27  CD  GLN A   4      -2.488  -2.371   0.060  1.00 11.36           C
ATOM     28  NE2 GLN A   4      -2.690  -3.328   1.055  1.00 10.43           N
ATOM     29  OE1 GLN A   4      -2.751  -1.015   0.151  1.00 12.29           O
ATOM     30  N   GLN A   5      -2.163  -2.714  -4.857  1.00 10.48           N
ATOM     31  CA  GLN A   5      -3.260  -2.063  -5.630  1.00 11.24           C
ATOM     32  C   GLN A   5      -4.607  -2.687  -5.178  1.00 11.40           C
ATOM     33  O   GLN A   5      -4.805  -4.059  -5.092  1.00 11.94           O
ATOM     34  CB  GLN A   5      -3.037  -2.282  -7.129  1.00 11.14           C
ATOM     35  CG  GLN A   5      -1.863  -1.378  -7.653  1.00 10.65           C
ATOM     36  CD  GLN A   5      -1.348  -1.972  -8.965  1.00 10.73           C
ATOM     37  NE2 GLN A   5      -1.523  -1.092 -10.027  1.00 11.31           N
ATOM     38  OE1 GLN A   5      -0.803  -3.209  -9.028  1.00 10.14           O
ATOM     39  N   ASN A   6      -5.522  -1.659  -4.872  1.00 11.56           N
ATOM     40  CA  ASN A   6      -6.830  -2.073  -4.359  1.00 12.07           C
ATOM     41  C   ASN A   6      -7.880  -1.587  -5.353  1.00 13.18           C
ATOM     42  O   ASN A   6      -8.282  -0.283  -5.362  1.00 13.64           O
ATOM     43  CB  ASN A   6      -7.073  -1.296  -3.017  1.00 12.12           C
ATOM     44  CG  ASN A   6      -5.983  -1.631  -1.998  1.00 12.31           C
ATOM     45  ND2 ASN A   6      -5.212  -0.477  -1.645  1.00 11.88           N
ATOM     46  OD1 ASN A   6      -5.849  -2.909  -1.517  1.00 13.43           O
ATOM     47  N   TYR A   7      -8.324  -2.630  -6.193  1.00 14.34           N
ATOM     48  CA  TYR A   7      -9.191  -2.233  -7.317  1.00 15.00           C
ATOM     49  C   TYR A   7     -10.635  -2.021  -6.893  1.00 15.64           C
ATOM     50  O   TYR A   7     -11.045  -2.569  -5.838  1.00 15.68           O
ATOM     51  CB  TYR A   7      -9.093  -3.430  -8.412  1.00 15.31           C
ATOM     52  CG  TYR A   7      -7.683  -3.606  -8.931  1.00 15.06           C
ATOM     53  CD1 TYR A   7      -6.809  -4.645  -8.368  1.00 15.24           C
ATOM     54  CD2 TYR A   7      -7.219  -2.707  -9.960  1.00 14.96           C
ATOM     55  CE1 TYR A   7      -5.510  -4.797  -8.830  1.00 14.94           C
ATOM     56  CE2 TYR A   7      -5.928  -2.849 -10.429  1.00 15.13           C
ATOM     57  CZ  TYR A   7      -5.075  -3.892  -9.861  1.00 14.97           C
ATOM     58  OH  TYR A   7      -3.785  -4.029 -10.335  1.00 14.93           O
ATOM     59  OXT TYR A   7     -11.413  -1.296  -7.601  1.00 15.89           O1-
TER      59      TYR A   7
HETATM   60  C   ACT A   8      -0.915  -6.207  -5.572  1.00 18.56           C
HETATM   61  O   ACT A   8       0.141  -6.081  -6.241  1.00 19.26           O
HETATM   62  CH3 ACT A   8      -1.548  -7.692  -5.404  1.00 18.08           C
HETATM   63  OXT ACT A   8      -1.457  -5.087  -5.051  1.00 18.36           O
HETATM   64  O   HOH A   9     -10.495  -2.647  -3.168  1.00 17.57           O
HETATM   65  O   HOH A  10       6.480   1.219  -7.070  1.00 21.27           O
HETATM   66  O   HOH A  11     -11.846   0.122  -9.956  1.00 27.52           O
HETATM   67  O   HOH A  12       1.576  -3.897 -11.035  1.00 44.76           O
TER      68      HOH A  12
ATOM     69  N   GLY B   1      -9.010   5.535   6.096  1.00 16.23           N
ATOM     70  CA  GLY B   1      -9.052   5.053   4.707  1.00 16.20           C
ATOM     71  C   GLY B   1      -8.015   3.839   4.448  1.00 15.91           C
ATOM     72  O   GLY B   1      -7.563   3.095   5.385  1.00 16.11           O
ATOM     73  N   ASN B   2      -7.642   3.626   3.185  1.00 15.02           N
ATOM     74  CA  ASN B   2      -6.537   2.617   2.848  1.00 13.92           C
ATOM     75  C   ASN B   2      -5.231   3.175   3.428  1.00 12.24           C
ATOM     76  O   ASN B   2      -4.965   4.529   3.432  1.00 11.42           O
ATOM     77  CB  ASN B   2      -6.391   2.450   1.340  1.00 14.42           C
ATOM     78  CG  ASN B   2      -7.639   1.858   0.697  1.00 13.92           C
ATOM     79  ND2 ASN B   2      -8.253   2.772  -0.168  1.00 12.78           N
ATOM     80  OD1 ASN B   2      -8.062   0.584   0.978  1.00 14.39           O
ATOM     81  N   ASN B   3      -4.414   2.120   3.904  1.00 12.20           N
ATOM     82  CA  ASN B   3      -3.168   2.489   4.594  1.00 11.81           C
ATOM     83  C   ASN B   3      -1.927   1.853   3.910  1.00 10.59           C
ATOM     84  O   ASN B   3      -1.840   0.499   3.653  1.00 10.32           O
ATOM     85  CB  ASN B   3      -3.235   1.924   6.046  1.00 12.51           C
ATOM     86  CG  ASN B   3      -1.973   2.287   6.848  1.00 12.59           C
ATOM     87  ND2 ASN B   3      -1.258   1.117   7.268  1.00 12.39           N
ATOM     88  OD1 ASN B   3      -1.661   3.603   7.106  1.00 13.64           O
ATOM     89  N   GLN B   4      -0.973   2.842   3.608  1.00 10.34           N
ATOM     90  CA  GLN B   4       0.369   2.366   3.208  1.00 10.00           C
ATOM     91  C   GLN B   4       1.409   3.154   4.085  1.00 10.21           C
ATOM     92  O   GLN B   4       1.521   4.533   4.070  1.00  8.99           O
ATOM     93  CB  GLN B   4       0.659   2.714   1.736  1.00 10.00           C
ATOM     94  CG  GLN B   4       1.933   1.949   1.229  1.00 10.50           C
ATOM     95  CD  GLN B   4       2.488   2.629  -0.060  1.00 11.36           C
ATOM     96  NE2 GLN B   4       2.690   1.672  -1.055  1.00 10.43           N
ATOM     97  OE1 GLN B   4       2.751   3.985  -0.151  1.00 12.29           O
ATOM     98  N   GLN B   5       2.163   2.286   4.857  1.00 10.48           N
ATOM     99  CA  GLN B   5       3.260   2.937   5.630  1.00 11.24           C
ATOM    100  C   GLN B   5       4.607   2.313   5.178  1.00 11.40           C
ATOM    101  O   GLN B   5       4.805   0.941   5.092  1.00 11.94           O
ATOM    102  CB  GLN B   5       3.037   2.718   7.129  1.00 11.14           C
ATOM    103  CG  GLN B   5       1.863   3.622   7.653  1.00 10.65           C
ATOM    104  CD  GLN B   5       1.348   3.028   8.965  1.00 10.73           C
ATOM    105  NE2 GLN B   5       1.523   3.908  10.027  1.00 11.31           N
ATOM    106  OE1 GLN B   5       0.803   1.791   9.028  1.00 10.14           O
ATOM    107  N   ASN B   6       5.522   3.341   4.872  1.00 11.56           N
ATOM    108  CA  ASN B   6       6.830   2.927   4.359  1.00 12.07           C
ATOM    109  C   ASN B   6       7.880   3.413   5.353  1.00 13.18           C
ATOM    110  O   ASN B   6       8.282   4.717   5.362  1.00 13.64           O
ATOM    111  CB  ASN B   6       7.073   3.704   3.017  1.00 12.12           C
ATOM    112  CG  ASN B   6       5.983   3.369   1.998  1.00 12.31           C
ATOM    113  ND2 ASN B   6       5.212   4.523   1.645  1.00 11.88           N
ATOM    114  OD1 ASN B   6       5.849   2.091   1.517  1.00 13.43           O
ATOM    115  N   TYR B   7       8.324   2.370   6.193  1.00 14.34           N
ATOM    116  CA  TYR B   7       9.191   2.767   7.317  1.00 15.00           C
ATOM    117  C   TYR B   7      10.635   2.979   6.893  1.00 15.64           C
ATOM    118  O   TYR B   7      11.045   2.431   5.838  1.00 15.68           O
ATOM    119  CB  TYR B   7       9.093   1.570   8.412  1.00 15.31           C
ATOM    120  CG  TYR B   7       7.683   1.394   8.931  1.00 15.06           C
ATOM    121  CD1 TYR B   7       6.809   0.355   8.368  1.00 15.24           C
ATOM    122  CD2 TYR B   7       7.219   2.293   9.960  1.00 14.96           C
ATOM    123  CE1 TYR B   7       5.510   0.203   8.830  1.00 14.94           C
ATOM    124  CE2 TYR B   7       5.928   2.151  10.429  1.00 15.13           C
ATOM    125  CZ  TYR B   7       5.075   1.108   9.861  1.00 14.97           C
ATOM    126  OH  TYR B   7       3.785   0.971  10.335  1.00 14.93           O
ATOM    127  OXT TYR B   7      11.413   3.704   7.601  1.00 15.89           O1-
TER     128      TYR B   7
HETATM  133  C   ACT B   8       1.032  -1.120   5.864  1.00 23.56           C
HETATM  134  O   ACT B   8      -0.050  -0.928   6.472  1.00 24.26           O
HETATM  135  CH3 ACT B   8       1.698  -2.599   5.925  1.00 23.08           C
HETATM  136  OXT ACT B   8       1.573  -0.070   5.213  1.00 23.36           O
HETATM  129  O   HOH B   9      10.495   2.353   3.168  1.00 17.57           O
HETATM  130  O   HOH B  10      -6.480   6.219   7.070  1.00 21.27           O
HETATM  131  O   HOH B  11      11.846   5.122   9.956  1.00 27.52           O
HETATM  132  O   HOH B  12      -1.576   1.103  11.035  1.00 44.76           O
END
"""
  pdb_in = prefix + "_in.pdb"
  with open(pdb_in, "w") as f:
    f.write(pdb_raw)
  mtz_in = prefix + "_in.mtz"
  params = """
    high_resolution = 1.75
    r_free_flags_fraction = 0.1
    add_sigmas = True
    random_seed = 12345
    output {
      label = F
      type = *real complex
      file_name = %s
    }
  """ % mtz_in
  with open("%s_fmodel.eff" % prefix, "w") as f:
    f.write(params)
  assert (easy_run.fully_buffered(
      "phenix.fmodel %s %s_fmodel.eff" % (pdb_in, prefix)
    ).raise_if_errors().return_code == 0)
  return pdb_in, mtz_in

def exercise():
  #
  # Test command-line program
  #
  pdb_in, mtz_in = make_inputs()
  import iotbx.pdb
  pdb_file = iotbx.pdb.input(pdb_in)
  hierarchy = pdb_file.construct_hierarchy()
  old_ligand = None
  for chain in hierarchy.only_model().chains():
    if (chain.id != "B") : continue
    for residue_group in chain.residue_groups():
      atom_group = residue_group.only_atom_group()
      if (atom_group.resname == "ACT"):
        old_ligand = atom_group.detached_copy()
        residue_group.remove_atom_group(atom_group)
        break
  assert old_ligand is not None
  with open("tst_ligand_ncs_start.pdb", "w") as f:
    f.write(hierarchy.as_pdb_string(
      crystal_symmetry=pdb_file.crystal_symmetry()))
  args = [
    "tst_ligand_ncs_start.pdb",
    mtz_in,
    "ligand_code=ACT",
  ]
  from mmtbx.command_line import apply_ncs_to_ligand
  if op.isfile("ncs_ligands.pdb"):
    os.remove("ncs_ligands.pdb")
  result = apply_ncs_to_ligand.run(args=args, out=null_out())
  assert result.n_ligands_new == 1
  assert op.isfile("ncs_ligands.pdb")
  pdb_out = iotbx.pdb.input("ncs_ligands.pdb")
  hierarchy_new = pdb_out.construct_hierarchy()
  new_ligand = None
  for chain in hierarchy_new.only_model().chains():
    if (chain.id != "B") : continue
    for residue_group in chain.residue_groups():
      atom_group = residue_group.only_atom_group()
      if (atom_group.resname == "ACT"):
        new_ligand = atom_group.detached_copy()
  assert new_ligand is not None
  rmsd = old_ligand.atoms().extract_xyz().rms_difference(
    new_ligand.atoms().extract_xyz())
  assert (rmsd < 0.5)
  #
  # Unit tests
  #

  # Excessive testing of NCS search machinery although with ligands. May be removed.
  import iotbx.ncs
  ncs_obj = iotbx.ncs.input(hierarchy=hierarchy, log=null_out())
  nrgl = ncs_obj.get_ncs_restraints_group_list()
  assert len(nrgl) == 1, len(nrgl)
  assert len(nrgl[0].master_iselection) == 59, len(nrgl[0].master_iselection)

if (__name__ == "__main__"):
  exercise()
  print("OK")
