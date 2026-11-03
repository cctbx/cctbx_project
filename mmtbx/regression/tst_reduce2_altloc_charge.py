"""The total charge reduce2 qmready=True reports for a model with alternate
conformations (mmtbx.ligands.electrons, which counts each conformer, altloc
blank + X, on its own and reports the charge they agree on). Counted on the
whole model the bonds of the shared atoms to both conformers made the DNA
trimer below, -2, come out as -9 and the protein fragment, -1, as -4. The
charge is measured on the finished model, all H placed and the dangling
atoms capped (qmready=True capping=True).
"""
from __future__ import absolute_import, division, print_function
import os
import mmtbx.model
from iotbx.cli_parser import run_program
from libtbx.utils import null_out
from mmtbx.ligands import electrons
from mmtbx.programs import reduce2

# 1BNA A 1-3: a 5'-OH end, two phosphates, a 3' end; charge -2
pdb_dna3 = """\
CRYST1   24.870   40.390   66.200  90.00  90.00  90.00 P 21 21 21
ATOM      1  O5'  DC A   1      18.935  34.195  25.617  1.00 64.35           O
ATOM      2  C5'  DC A   1      19.130  33.921  24.219  1.00 44.69           C
ATOM      3  C4'  DC A   1      19.961  32.668  24.100  1.00 31.28           C
ATOM      4  O4'  DC A   1      19.360  31.583  24.852  1.00 37.45           O
ATOM      5  C3'  DC A   1      20.172  32.122  22.694  1.00 46.72           C
ATOM      6  O3'  DC A   1      21.350  31.325  22.681  1.00 48.89           O
ATOM      7  C2'  DC A   1      18.948  31.223  22.647  1.00 30.88           C
ATOM      8  C1'  DC A   1      19.231  30.482  23.944  1.00 36.58           C
ATOM      9  N1   DC A   1      18.070  29.661  24.380  1.00 40.51           N
ATOM     10  C2   DC A   1      18.224  28.454  25.015  1.00 16.62           C
ATOM     11  O2   DC A   1      19.360  28.014  25.214  1.00 27.75           O
ATOM     12  N3   DC A   1      17.143  27.761  25.377  1.00 20.55           N
ATOM     13  C4   DC A   1      15.917  28.226  25.120  1.00 34.72           C
ATOM     14  N4   DC A   1      14.828  27.477  25.444  1.00 40.31           N
ATOM     15  C5   DC A   1      15.719  29.442  24.471  1.00 30.78           C
ATOM     16  C6   DC A   1      16.843  30.171  24.101  1.00 25.90           C
ATOM     17  P    DG A   2      22.409  31.286  21.483  1.00 58.85           P
ATOM     18  OP1  DG A   2      23.536  32.157  21.851  1.00 57.82           O
ATOM     19  OP2  DG A   2      21.822  31.459  20.139  1.00 78.33           O
ATOM     20  O5'  DG A   2      22.840  29.751  21.498  1.00 40.36           O
ATOM     21  C5'  DG A   2      23.543  29.175  22.594  1.00 47.19           C
ATOM     22  C4'  DG A   2      23.494  27.709  22.279  1.00 47.81           C
ATOM     23  O4'  DG A   2      22.193  27.252  22.674  1.00 38.76           O
ATOM     24  C3'  DG A   2      23.693  27.325  20.807  1.00 28.58           C
ATOM     25  O3'  DG A   2      24.723  26.320  20.653  1.00 40.44           O
ATOM     26  C2'  DG A   2      22.273  26.885  20.416  1.00 21.14           C
ATOM     27  C1'  DG A   2      21.721  26.304  21.716  1.00 33.95           C
ATOM     28  N9   DG A   2      20.237  26.470  21.780  1.00 34.00           N
ATOM     29  C8   DG A   2      19.526  27.584  21.429  1.00 36.47           C
ATOM     30  N7   DG A   2      18.207  27.455  21.636  1.00 32.37           N
ATOM     31  C5   DG A   2      18.083  26.212  22.142  1.00 15.06           C
ATOM     32  C6   DG A   2      16.904  25.525  22.545  1.00 11.88           C
ATOM     33  O6   DG A   2      15.739  25.916  22.518  1.00 21.30           O
ATOM     34  N1   DG A   2      17.197  24.279  23.037  1.00 15.44           N
ATOM     35  C2   DG A   2      18.434  23.717  23.155  1.00  9.63           C
ATOM     36  N2   DG A   2      18.508  22.456  23.668  1.00 16.69           N
ATOM     37  N3   DG A   2      19.537  24.360  22.770  1.00 30.98           N
ATOM     38  C4   DG A   2      19.290  25.594  22.274  1.00 18.56           C
ATOM     39  P    DC A   3      25.064  25.621  19.252  1.00 44.67           P
ATOM     40  OP1  DC A   3      26.506  25.316  19.220  1.00 53.89           O
ATOM     41  OP2  DC A   3      24.559  26.412  18.115  1.00 57.79           O
ATOM     42  O5'  DC A   3      24.260  24.246  19.327  1.00 35.42           O
ATOM     43  C5'  DC A   3      24.584  23.285  20.335  1.00 45.75           C
ATOM     44  C4'  DC A   3      23.523  22.233  20.245  1.00 43.02           C
ATOM     45  O4'  DC A   3      22.256  22.844  20.453  1.00 36.85           O
ATOM     46  C3'  DC A   3      23.424  21.557  18.903  1.00 40.14           C
ATOM     47  O3'  DC A   3      24.121  20.309  18.928  1.00 49.62           O
ATOM     48  C2'  DC A   3      21.930  21.406  18.661  1.00 53.79           C
ATOM     49  C1'  DC A   3      21.278  21.966  19.909  1.00 22.18           C
ATOM     50  N1   DC A   3      20.196  22.889  19.521  1.00 25.44           N
ATOM     51  C2   DC A   3      18.909  22.584  19.816  1.00 19.81           C
ATOM     52  O2   DC A   3      18.685  21.512  20.382  1.00 29.92           O
ATOM     53  N3   DC A   3      17.935  23.447  19.502  1.00 21.59           N
ATOM     54  C4   DC A   3      18.217  24.603  18.897  1.00 14.01           C
ATOM     55  N4   DC A   3      17.221  25.499  18.629  1.00 26.88           N
ATOM     56  C5   DC A   3      19.526  24.945  18.571  1.00 27.59           C
ATOM     57  C6   DC A   3      20.537  24.048  18.899  1.00 27.05           C
TER
END
"""

# LYS A 1 truncated at CD (capped by capping=True; no NZ, so 0), a chain
# break, ASP A 5 - GLY A 6: NH3+ +1 and COO- -1 on each fragment, ASP -1:
# charge -1
pdb_lys_asp_gly = """\
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1
ATOM      1  N   LYS A   1      27.221  24.554  25.195  1.00  5.98           N
ATOM      2  CA  LYS A   1      25.851  24.605  24.700  1.00  6.43           C
ATOM      3  C   LYS A   1      25.466  26.033  24.327  1.00  6.30           C
ATOM      4  O   LYS A   1      25.212  26.865  25.198  1.00  6.91           O
ATOM      5  CB  LYS A   1      24.884  24.064  25.755  1.00  7.86           C
ATOM      6  CG  LYS A   1      23.430  24.044  25.312  1.00 11.05           C
ATOM      7  CD  LYS A   1      22.527  23.473  26.393  1.00 16.83           C
ATOM      8  N   ASP A   5       9.580   7.432   8.699  1.00 37.74           N
ATOM      9  CA  ASP A   5      10.692   8.369   8.630  1.00 38.33           C
ATOM     10  C   ASP A   5      10.828   9.231   9.884  1.00 37.17           C
ATOM     11  O   ASP A   5      11.938   9.511  10.322  1.00 36.50           O
ATOM     12  CB  ASP A   5      10.625   9.197   7.354  1.00 39.43           C
ATOM     13  CG  ASP A   5      11.281   8.467   6.194  1.00 43.72           C
ATOM     14  OD1 ASP A   5      10.549   7.831   5.406  1.00 49.92           O
ATOM     15  OD2 ASP A   5      12.531   8.382   6.067  1.00 51.10           O
ATOM     16  N   GLY A   6       9.691   9.552  10.506  1.00 36.16           N
ATOM     17  CA  GLY A   6       9.668  10.262  11.762  1.00 35.61           C
ATOM     18  C   GLY A   6      10.374   9.397  12.850  1.00 34.61           C
ATOM     19  O   GLY A   6      11.228   9.895  13.566  1.00 33.43           O
END
"""

def split(text, resseq, names, shift=0.3):
  """text with the named atoms of residue resseq split into altlocs A and
  B (occupancy 0.5 each, B shifted along x)"""
  out = []
  for l in text.splitlines():
    if (l.startswith("ATOM") and l[22:26].strip() == str(resseq)
        and l[12:16].strip() in names):
      a = l[:16] + "A" + l[17:54] + "  0.50" + l[60:]
      x = float(l[30:38]) + shift
      b = l[:16] + "B" + l[17:30] + "%8.3f" % x + l[38:54] + "  0.50" + l[60:]
      out += [a, b]
    else:
      out.append(l)
  return "\n".join(out) + "\n"

def run(text, prefix):
  fn = prefix + ".pdb"
  with open(fn, "w") as f:
    f.write(text)
  out = prefix + "_out.pdb"
  if os.path.exists(out): os.remove(out)
  return run_program(program_class=reduce2.Program, logger=null_out(),
    args=[fn, "qmready=True", "capping=True", "overwrite=True",
          "output.filename=" + out])

def conformer_charge(model, altloc):
  """electrons on the altloc-free model of one conformer"""
  sel = model.selection("altloc ' ' or altloc %s" % altloc)
  conformer = model.select(sel)
  conformer.process(make_restraints=True)
  return electrons.run(model=conformer,
                       return_formal_charges=True).get_total_charge()

def test_dna_split_base():
  """The base of DC A 3 in altlocs A and B; capping adds nothing, so the
  finished model can be interpreted again and each conformer alone gives
  electrons the right charge, -2, as reduce2 reports."""
  base = ["N1", "C2", "O2", "N3", "C4", "N4", "C5", "C6"]
  r = run(split(pdb_dna3, 3, base), "tst_reduce2_altloc_charge_dna")
  model = mmtbx.model.manager(
    model_input      = None,
    pdb_hierarchy    = r.model.get_hierarchy(),
    crystal_symmetry = r.model.crystal_symmetry(),
    log              = null_out())
  for altloc in ("A", "B"):
    assert conformer_charge(model, altloc) == -2, altloc
  print("DNA trimer, base of DC A 3 split: reduce2 total charge",
        r.total_charge, "expected -2")
  assert r.total_charge == -2, r.total_charge

def test_protein_split_side_chain():
  """The side chain of ASP A 5 in altlocs A and B next to a truncated,
  capped LYS; -1."""
  r = run(split(pdb_lys_asp_gly, 5, ["CB", "CG", "OD1", "OD2"]),
          "tst_reduce2_altloc_charge_protein")
  print("LYS-ASP-GLY, side chain of ASP A 5 split: reduce2 total charge",
        r.total_charge, "expected -1")
  assert r.total_charge == -1, r.total_charge

if __name__ == "__main__":
  test_protein_split_side_chain()
  test_dna_split_base()
  print("OK")
