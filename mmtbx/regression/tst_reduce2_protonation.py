"""reduce2 protonation: the protonation state of selected residues (I1 step
1). Fixtures:

- frag: 136 heavy atoms cut from 1PQ7 (P 1): Ala 55 - Ala 62 (His 56, Cys
  57, Tyr 61), Asn 98 - Leu 100 (Asp 99), Ile 102 - Leu 105 (Lys 104), Ser
  142 - Gly 144 (Glu 143). Asp 99 and Glu 143 are named as IUPAC wants
  (chi2 13 deg, chi3 -33 deg); frag_swapped exchanges their O names, so
  the IUPAC flip would rename them.
- his_cd: 1XMK His A 354, NE2 2.3 A from the Cd A 397 of the mate
  -x-1/2,-y+2,z-1/2 (P 21 21 21).
- his_zn: the same His with a Zn 2.05 A from NE2 in the asymmetric unit;
  his_zn_altlocs puts the His in conformer A and the Zn in conformer B.
- frag_his_split: frag with His 56 split, A complete, B truncated after CG
  (or complete).
- frag_two_his: frag with a copy of His 56 (His 256) in empty space.
- His 56 before the NAG-Asn glycan of tst_reduce2_blank_parent_split (whose
  split Asn makes placement share H between covalent partners).

States are given in an .eff file: on the command line each word of a
repeated scope starts a new entry.
"""
import os
from io import StringIO
import iotbx.pdb
import libtbx.load_env
import mmtbx.model
from libtbx.utils import Sorry, null_out
from iotbx.cli_parser import run_program
from scitbx import matrix
from mmtbx.programs import reduce2
from mmtbx.hydrogens import reduce_hydrogen
from mmtbx.regression.tst_reduce2_blank_parent_split import pdb_glycan

frag_str = """\
CRYST1   32.869   37.022   39.780 102.89 104.59 102.37 P 1           1
ATOM    560  N   ALA A  55       7.590  -3.374   2.133  1.00  3.91           N
ATOM    561  CA  ALA A  55       8.954  -3.495   2.628  1.00  3.93           C
ATOM    562  C   ALA A  55       9.877  -2.455   2.017  1.00  3.99           C
ATOM    563  O   ALA A  55      10.763  -1.938   2.716  1.00  4.54           O
ATOM    564  CB  ALA A  55       9.495  -4.901   2.403  1.00  4.35           C
ATOM    570  N   HIS A  56       9.728  -2.142   0.726  1.00  4.02           N
ATOM    571  CA  HIS A  56      10.687  -1.249   0.118  1.00  4.42           C
ATOM    572  C   HIS A  56      10.576   0.171   0.697  1.00  4.38           C
ATOM    573  O   HIS A  56      11.534   0.956   0.573  1.00  5.25           O
ATOM    574  CB  HIS A  56      10.592  -1.249  -1.407  1.00  4.99           C
ATOM    575  CG  HIS A  56       9.480  -0.462  -1.999  1.00  5.02           C
ATOM    576  ND1 HIS A  56       8.294  -0.988  -2.415  1.00  5.14           N
ATOM    577  CD2 HIS A  56       9.435   0.875  -2.263  1.00  5.76           C
ATOM    578  CE1 HIS A  56       7.560  -0.010  -2.923  1.00  5.72           C
ATOM    579  NE2 HIS A  56       8.216   1.121  -2.838  1.00  5.92           N
ATOM    587  N   CYS A  57       9.445   0.516   1.315  1.00  4.44           N
ATOM    588  CA  CYS A  57       9.313   1.815   1.946  1.00  4.72           C
ATOM    589  C   CYS A  57      10.241   1.988   3.148  1.00  4.55           C
ATOM    590  O   CYS A  57      10.553   3.132   3.509  1.00  5.19           O
ATOM    591  CB  CYS A  57       7.875   2.020   2.439  1.00  5.19           C
ATOM    592  SG  CYS A  57       6.674   2.294   1.132  1.00  6.27           S
ATOM    597  N   VAL A  58      10.599   0.888   3.802  1.00  4.24           N
ATOM    598  CA  VAL A  58      11.264   0.950   5.096  1.00  4.21           C
ATOM    599  C   VAL A  58      12.556   0.168   5.196  1.00  4.23           C
ATOM    600  O   VAL A  58      13.370   0.456   6.092  1.00  5.18           O
ATOM    601  CB  VAL A  58      10.295   0.532   6.237  1.00  4.70           C
ATOM    602  CG1 VAL A  58       9.187   1.577   6.396  1.00  5.80           C
ATOM    603  CG2 VAL A  58       9.702  -0.860   6.006  1.00  5.12           C
ATOM    613  N   SER A  59      12.778  -0.854   4.356  1.00  4.26           N
ATOM    614  CA  SER A  59      13.878  -1.784   4.577  1.00  4.25           C
ATOM    615  C   SER A  59      15.218  -1.082   4.583  1.00  4.34           C
ATOM    616  O   SER A  59      15.584  -0.413   3.585  1.00  4.78           O
ATOM    617  CB  SER A  59      13.848  -2.830   3.468  1.00  4.63           C
ATOM    618  OG  SER A  59      14.845  -3.824   3.670  1.00  5.00           O
ATOM    624  N   GLY A  60      15.958  -1.240   5.666  1.00  4.87           N
ATOM    625  CA  GLY A  60      17.316  -0.715   5.786  1.00  5.12           C
ATOM    626  C   GLY A  60      17.442   0.783   6.005  1.00  4.78           C
ATOM    627  O   GLY A  60      18.577   1.250   6.182  1.00  5.77           O
ATOM    631  N   TYR A  61      16.349   1.550   5.974  1.00  4.94           N
ATOM    632  CA  TYR A  61      16.423   2.984   6.193  1.00  4.82           C
ATOM    633  C   TYR A  61      16.454   3.256   7.702  1.00  4.90           C
ATOM    634  O   TYR A  61      16.013   2.471   8.506  1.00  6.74           O
ATOM    635  CB  TYR A  61      15.220   3.667   5.574  1.00  4.94           C
ATOM    636  CG  TYR A  61      15.121   3.593   4.077  1.00  4.69           C
ATOM    637  CD1 TYR A  61      16.129   4.014   3.237  1.00  5.43           C
ATOM    638  CD2 TYR A  61      13.935   3.177   3.486  1.00  5.20           C
ATOM    639  CE1 TYR A  61      15.979   3.988   1.850  1.00  5.86           C
ATOM    640  CE2 TYR A  61      13.755   3.183   2.105  1.00  5.81           C
ATOM    641  CZ  TYR A  61      14.779   3.580   1.287  1.00  5.43           C
ATOM    642  OH  TYR A  61      14.591   3.579  -0.066  1.00  7.78           O
ATOM    652  N   ALA A  62      16.906   4.474   8.057  1.00  4.72           N
ATOM    653  CA  ALA A  62      16.882   4.852   9.472  1.00  5.12           C
ATOM    654  C   ALA A  62      15.428   4.953   9.936  1.00  5.55           C
ATOM    655  O   ALA A  62      14.640   5.767   9.424  1.00  5.94           O
ATOM    656  CB  ALA A  62      17.574   6.176   9.689  1.00  5.53           C
ATOM   1178  N   ASN A  98       7.463 -10.332  -2.229  1.00  4.08           N
ATOM   1179  CA  ASN A  98       7.087 -10.379  -0.817  1.00  4.07           C
ATOM   1180  C   ASN A  98       7.305  -9.020  -0.138  1.00  4.19           C
ATOM   1181  O   ASN A  98       7.889  -8.901   0.936  1.00  4.42           O
ATOM   1182  CB  ASN A  98       7.782 -11.518  -0.089  1.00  4.36           C
ATOM   1183  CG  ASN A  98       7.540 -12.804  -0.862  1.00  4.57           C
ATOM   1184  OD1 ASN A  98       6.421 -13.101  -1.249  1.00  4.98           O
ATOM   1185  ND2 ASN A  98       8.615 -13.563  -1.080  1.00  5.57           N
ATOM   1192  N   ASP A  99       6.781  -7.997  -0.822  1.00  4.39           N
ATOM   1193  CA  ASP A  99       7.059  -6.598  -0.493  1.00  4.18           C
ATOM   1194  C   ASP A  99       6.109  -6.069   0.574  1.00  4.05           C
ATOM   1195  O   ASP A  99       5.331  -5.131   0.339  1.00  4.55           O
ATOM   1196  CB  ASP A  99       6.985  -5.767  -1.770  1.00  4.35           C
ATOM   1197  CG  ASP A  99       7.501  -4.341  -1.626  1.00  4.17           C
ATOM   1198  OD1 ASP A  99       8.150  -4.039  -0.619  1.00  4.39           O
ATOM   1199  OD2 ASP A  99       7.252  -3.544  -2.578  1.00  5.58           O
ATOM   1204  N   LEU A 100       6.207  -6.657   1.782  1.00  4.36           N
ATOM   1205  CA  LEU A 100       5.439  -6.208   2.931  1.00  4.02           C
ATOM   1206  C   LEU A 100       6.363  -6.088   4.143  1.00  4.04           C
ATOM   1207  O   LEU A 100       7.309  -6.858   4.310  1.00  4.40           O
ATOM   1208  CB  LEU A 100       4.282  -7.134   3.342  1.00  5.00           C
ATOM   1209  CG  LEU A 100       3.097  -7.175   2.389  1.00  4.94           C
ATOM   1210  CD1 LEU A 100       2.024  -8.133   2.937  1.00  6.39           C
ATOM   1211  CD2 LEU A 100       2.474  -5.822   2.145  1.00  5.40           C
ATOM   1223  N   ALA A 101       5.997  -5.151   5.013  1.00  4.16           N
ATOM   1224  CA  ALA A 101       6.631  -4.986   6.301  1.00  4.24           C
ATOM   1225  C   ALA A 101       5.603  -4.632   7.342  1.00  4.06           C
ATOM   1226  O   ALA A 101       4.598  -3.960   7.077  1.00  4.70           O
ATOM   1227  CB  ALA A 101       7.691  -3.874   6.255  1.00  4.90           C
ATOM   1233  N   ILE A 102       5.887  -5.091   8.578  1.00  4.29           N
ATOM   1234  CA  ILE A 102       5.137  -4.764   9.778  1.00  4.41           C
ATOM   1235  C   ILE A 102       5.921  -3.725  10.564  1.00  4.51           C
ATOM   1236  O   ILE A 102       7.104  -3.938  10.901  1.00  5.47           O
ATOM   1237  CB  ILE A 102       4.902  -6.034  10.607  1.00  5.19           C
ATOM   1238  CG1 ILE A 102       3.862  -6.950   9.923  1.00  5.91           C
ATOM   1239  CG2 ILE A 102       4.517  -5.746  12.052  1.00  6.20           C
ATOM   1240  CD1 ILE A 102       2.424  -6.481  10.069  1.00  7.08           C
ATOM   1252  N   LEU A 103       5.265  -2.615  10.869  1.00  4.53           N
ATOM   1253  CA  LEU A 103       5.831  -1.546  11.686  1.00  4.65           C
ATOM   1254  C   LEU A 103       5.131  -1.566  13.014  1.00  4.74           C
ATOM   1255  O   LEU A 103       3.886  -1.626  13.072  1.00  5.84           O
ATOM   1256  CB  LEU A 103       5.616  -0.181  11.040  1.00  5.12           C
ATOM   1257  CG  LEU A 103       5.819  -0.085   9.526  1.00  5.52           C
ATOM   1258  CD1 LEU A 103       5.629   1.356   9.078  1.00  5.96           C
ATOM   1259  CD2 LEU A 103       7.158  -0.614   9.101  1.00  6.25           C
ATOM   1271  N   LYS A 104       5.896  -1.473  14.095  1.00  5.34           N
ATOM   1272  CA  LYS A 104       5.312  -1.279  15.428  1.00  5.61           C
ATOM   1273  C   LYS A 104       5.467   0.193  15.783  1.00  5.77           C
ATOM   1274  O   LYS A 104       6.506   0.817  15.526  1.00  6.99           O
ATOM   1275  CB  LYS A 104       5.966  -2.167  16.462  1.00  6.79           C
ATOM   1276  CG  LYS A 104       5.735  -3.629  16.196  1.00  7.88           C
ATOM   1277  CD  LYS A 104       4.269  -4.013  16.397  1.00  8.60           C
ATOM   1278  CE  LYS A 104       4.073  -5.431  15.965  1.00 10.43           C
ATOM   1279  NZ  LYS A 104       2.923  -6.097  16.468  1.00 24.27           N
ATOM   1293  N   LEU A 105       4.402   0.736  16.384  1.00  5.91           N
ATOM   1294  CA  LEU A 105       4.298   2.145  16.681  1.00  6.13           C
ATOM   1295  C   LEU A 105       4.719   2.430  18.119  1.00  6.31           C
ATOM   1296  O   LEU A 105       4.362   1.747  19.054  1.00  8.03           O
ATOM   1297  CB  LEU A 105       2.847   2.588  16.523  1.00  6.42           C
ATOM   1298  CG  LEU A 105       2.254   2.378  15.127  1.00  6.45           C
ATOM   1299  CD1 LEU A 105       0.847   2.970  15.077  1.00  7.95           C
ATOM   1300  CD2 LEU A 105       3.116   3.001  14.024  1.00  7.45           C
ATOM   1798  N   SER A 142       1.182  12.172 -11.782  0.73  5.79           N
ATOM   1800  CA  SER A 142       1.739  11.819 -13.076  0.74  5.55           C
ATOM   1802  C   SER A 142       3.059  11.110 -12.827  0.68  5.95           C
ATOM   1804  O   SER A 142       3.708  11.379 -11.781  0.70  6.00           O
ATOM   1806  CB  SER A 142       1.929  13.065 -13.950  0.69  6.71           C
ATOM   1808  OG  SER A 142       2.820  13.971 -13.380  0.67  9.34           O
ATOM   1820  N   GLU A 143       3.533  10.281 -13.765  0.72  5.88           N
ATOM   1822  CA  GLU A 143       4.788   9.637 -13.603  0.70  6.14           C
ATOM   1824  C   GLU A 143       5.921  10.674 -13.486  0.69  6.00           C
ATOM   1826  O   GLU A 143       5.983  11.547 -14.346  0.69  6.28           O
ATOM   1828  CB  GLU A 143       4.965   8.732 -14.821  0.70  6.24           C
ATOM   1830  CG  GLU A 143       6.246   7.964 -14.687  0.66  7.23           C
ATOM   1832  CD  GLU A 143       6.522   7.033 -15.935  0.71  6.99           C
ATOM   1834  OE1 GLU A 143       5.499   6.536 -16.508  0.72  7.01           O
ATOM   1836  OE2 GLU A 143       7.727   6.837 -16.226  0.68 10.48           O
ATOM   1850  N   GLY A 144       6.677  10.594 -12.480  1.00  9.01           N
ATOM   1851  CA  GLY A 144       7.778  11.566 -12.287  1.00  9.65           C
ATOM   1852  C   GLY A 144       7.312  12.943 -11.827  1.00  9.10           C
ATOM   1853  O   GLY A 144       8.140  13.824 -11.816  1.00 10.11           O
END
"""

his_cd_str = """\
CRYST1   35.557   43.526   45.471  90.00  90.00  90.00 P 21 21 21    4
ATOM    548  N   HIS A 354      -1.093  10.746  39.862  1.00  6.48           N
ATOM    549  CA  HIS A 354      -0.970   9.793  40.975  1.00  6.84           C
ATOM    550  C   HIS A 354       0.400   9.117  40.895  1.00  6.70           C
ATOM    551  O   HIS A 354       0.877   8.812  39.808  1.00  7.52           O
ATOM    552  CB  HIS A 354      -2.004   8.651  40.833  1.00  8.27           C
ATOM    553  CG  HIS A 354      -3.476   9.048  40.932  1.00  8.55           C
ATOM    554  ND1 HIS A 354      -4.160   9.785  40.003  1.00  8.74           N
ATOM    555  CD2 HIS A 354      -4.350   8.756  41.946  1.00  8.81           C
ATOM    556  CE1 HIS A 354      -5.447   9.936  40.440  1.00  8.37           C
ATOM    557  NE2 HIS A 354      -5.605   9.331  41.604  1.00  8.44           N
HETATM  666 CD    CD A 397      10.242  12.438  48.077  1.00  9.61          CD
END
"""

def his_zn_str():
  # his_cd without the Cd, a Zn 2.05 A from NE2 on the line from the ring
  # centre through NE2
  lines = [l for l in his_cd_str.splitlines() if not l.startswith("HETATM")
           and not l.startswith("END")]
  xyz = dict((l[12:16].strip(), matrix.col((float(l[30:38]), float(l[38:46]),
    float(l[46:54])))) for l in lines if l.startswith("ATOM"))
  ring = ["CG", "ND1", "CD2", "CE1", "NE2"]
  centre = sum([xyz[n] for n in ring], matrix.col((0, 0, 0))) / 5
  zn = xyz["NE2"] + (xyz["NE2"] - centre).normalize() * 2.05
  lines.append("HETATM  666 ZN    ZN A 397    %8.3f%8.3f%8.3f  1.00  9.61"
    "          ZN" % tuple(zn))
  return "\n".join(lines + ["END", ""])

def his_zn_altlocs_str():
  # his_zn with the His (all atoms) in conformer A, the Zn in conformer B
  out = []
  for l in his_zn_str().splitlines():
    if l.startswith(("ATOM", "HETATM")):
      l = l[:16] + ("A" if l.startswith("ATOM") else "B") + l[17:54] + \
        "  0.50" + l[60:]
    out.append(l)
  return "\n".join(out) + "\n"

def frag_his_split_str(complete=False):
  # His 56: N, CA, C, O blank; conformer A complete, conformer B CB and CG only
  # (complete: all side-chain atoms), 0.1 A away
  out = []
  for l in frag_str.splitlines():
    if l.startswith("ATOM") and l[17:26] == "HIS A  56" and \
        l[12:16].strip() not in ("N", "CA", "C", "O"):
      out.append(l[:16] + "A" + l[17:54] + "  0.50" + l[60:])
      if complete or l[12:16].strip() in ("CB", "CG"):
        out.append(l[:16] + "B" + l[17:30] + "%8.3f" % (float(l[30:38]) + 0.1)
          + l[38:54] + "  0.50" + l[60:])
      continue
    out.append(l)
  return "\n".join(out) + "\n"

def his_56_glycan_str():
  # His 56 of frag, then the glycan fixture's atoms (its cell)
  his = [l for l in frag_str.splitlines()
         if l.startswith("ATOM") and l[17:26] == "HIS A  56"]
  glycan = pdb_glycan.splitlines()
  return "\n".join(glycan[:1] + his + glycan[1:]) + "\n"

def frag_swapped_str():
  # Asp 99 OD1/OD2 and Glu 143 OE1/OE2 exchange names (and so chi2/chi3 is
  # beyond 90 deg)
  swap = {(" 99", " OD1"): " OD2", (" 99", " OD2"): " OD1",
          ("143", " OE1"): " OE2", ("143", " OE2"): " OE1"}
  out = []
  for l in frag_str.splitlines():
    key = (l[23:26], l[12:16])
    if l.startswith("ATOM") and key in swap:
      l = l[:12] + swap[key] + l[16:]
    out.append(l)
  return "\n".join(out) + "\n"

def frag_two_his_str():
  # frag with a copy of His 56 as His 256, moved by (10, -6, -10) A: 13.5 A
  # from every atom and lattice mate, so it gets a HisFlip and a His flipkin
  lines = frag_str.splitlines()
  copy = []
  for l in lines:
    if l.startswith("ATOM") and l[17:26] == "HIS A  56":
      x, y, z = float(l[30:38]) + 10, float(l[38:46]) - 6, float(l[46:54]) - 10
      copy.append(l[:22] + " 256" + l[26:30] + "%8.3f%8.3f%8.3f" % (x, y, z)
        + l[54:])
  return "\n".join(lines[:-1] + copy + lines[-1:]) + "\n"

def write(name, text):
  with open(name, "w") as f:
    f.write(text)

def eff(prefix, entries):
  # entries: [(selection, state)]
  name = prefix + "_states.eff"
  write(name, "".join('protonation {\n  selection = "%s"\n  state = %s\n}\n'
    % e for e in entries))
  return name

def run(prefix, pdb_text, entries=(), extra_args=()):
  """reduce2 on pdb_text with the protonation entries; returns (output
  model text, description, log)."""
  write(prefix + ".pdb", pdb_text)
  args = [prefix + ".pdb", "output.filename=" + prefix + "_out.pdb",
    "output.description_file_name=" + prefix + ".txt", "overwrite=True"]
  if entries: args.append(eff(prefix, entries))
  log = StringIO()
  run_program(program_class=reduce2.Program, logger=log,
    args=args + list(extra_args))
  with open(prefix + "_out.pdb") as f:
    model_text = f.read()
  with open(prefix + ".txt") as f:
    description = f.read()
  return model_text, description, log.getvalue()

def sorry(prefix, pdb_text, entries=(), extra_args=()):
  """The Sorry text of a reduce2 run that must stop."""
  try:
    run(prefix, pdb_text, entries, extra_args)
  except Sorry as e:
    return str(e)
  raise AssertionError("no Sorry")

def hierarchy(text):
  return iotbx.pdb.input(lines=text.split("\n"), source_info=None
    ).construct_hierarchy()

def atoms_of(text, resseq):
  """{name: atom} of a residue (one conformer, as in these fixtures)."""
  h = hierarchy(text)
  result = dict((a.name.strip(), a) for a in h.atoms()
    if a.parent().parent().resseq_as_int() == resseq)
  # the atoms keep their hierarchy (their parents) while in use
  atoms_of.keep.append(h)
  return result
atoms_of.keep = []

def h_names(text, resseq):
  return sorted(n for n, a in atoms_of(text, resseq).items()
    if a.element.strip() in ("H", "D"))

def place(pdb_text, states):
  """place_hydrogens with {resseq: state}; returns the placed model."""
  model = mmtbx.model.manager(model_input=iotbx.pdb.input(
    lines=pdb_text.split("\n"), source_info=None), log=null_out())
  keyed = {}
  for rg in model.get_hierarchy().residue_groups():
    if rg.resseq_as_int() in states:
      keyed[reduce_hydrogen.protonation_key(rg)] = states[rg.resseq_as_int()]
  p = reduce_hydrogen.place_hydrogens(model=model, protonation_states=keyed)
  p.run()
  return p

# ------------------------------------------------------------------------------

def test_states_placement():
  # Every state on one residue: the H present and left out, the riding type
  # of each new H, and for the acids the ACID-ASP/ACID-GLU modification
  # (HD2/HE2 in the carboxyl plane and the syn torsion).
  expected = {
    "his_hd1": {56: (["HD1"], ["HE2"])},
    "his_he2": {56: (["HE2"], ["HD1"])},
    "his_both": {56: (["HD1", "HE2"], [])},
    "protonated": {99: (["HD2"], []), 143: (["HE2"], [])},
    "deprotonated": {57: ([], ["HG"]), 61: ([], ["HH"])},
  }
  riding = {(56, "HD1"): "flat_2neigbs", (56, "HE2"): "flat_2neigbs",
            (99, "HD2"): "alg1b", (143, "HE2"): "alg1b"}
  acid = {99: ("CB", "CG", "OD1", "OD2", "HD2"),
          143: ("CG", "CD", "OE1", "OE2", "HE2")}
  for state, residues in expected.items():
    p = place(frag_str, dict((r, state) for r in residues))
    model = p.get_model()
    h = model.get_hierarchy()
    text = h.as_pdb_string()
    for resseq, (present, absent) in residues.items():
      names = h_names(text, resseq)
      for n in present: assert n in names, (state, resseq, n, names)
      for n in absent: assert n not in names, (state, resseq, n, names)
    assert p.protonation_conflicts == [], p.protonation_conflicts
    model.setup_riding_h_manager()
    hp = model.get_riding_h_manager().h_parameterization
    grm = model.get_restraints_manager().geometry
    for a in h.atoms():
      key = (a.parent().parent().resseq_as_int(), a.name.strip())
      if key in riding and a.parent().parent().resseq_as_int() in residues:
        assert hp[a.i_seq] is not None, key
        assert hp[a.i_seq].htype == riding[key], (key, hp[a.i_seq].htype)
    for resseq, names in acid.items():
      if resseq not in residues: continue
      i = dict((a.name.strip(), a.i_seq) for a in h.atoms()
        if a.parent().parent().resseq_as_int() == resseq)
      planes = [set(pp.i_seqs) for pp in grm.planarity_proxies]
      assert set(i[n] for n in names) in planes, (state, resseq)
      torsions = [tuple(dp.i_seqs) for dp in grm.dihedral_proxies]
      t = tuple(i[n] for n in (names[0], names[1], names[3], names[4]))
      assert t in torsions or t[::-1] in torsions, (state, resseq)
  # control: without states the dictionary H
  text = place(frag_str, {}).get_model().get_hierarchy().as_pdb_string()
  assert "HD1" in h_names(text, 56) and "HE2" in h_names(text, 56)
  assert "HD2" not in h_names(text, 99) and "HE2" not in h_names(text, 143)
  assert "HG" in h_names(text, 57) and "HH" in h_names(text, 61)

def test_report_and_optimizer():
  # One line per residue in the description; the His gets no Mover; the
  # acid H rotate and end in the carboxyl plane.
  out, desc, log = run("tst_reduce2_protonation_1", frag_str, [
    ("resseq 56", "his_he2"), ("resseq 99 or resseq 143", "protonated"),
    ("resseq 57 or resseq 61", "deprotonated")])
  for line in [
      "Protonation states set by the user (protonation):",
      "  HIS A 56: HE2 only (his_he2): HE2 added, HD1 left out",
      "  ASP A 99: protonated on OD2 (protonated): HD2 added",
      "  GLU A 143: protonated on OE2 (protonated): HE2 added",
      "  CYS A 57: deprotonated (deprotonated): HG left out",
      "  TYR A 61: deprotonated (deprotonated): HH left out"]:
    assert line + "\n" in desc, line
    assert line in log, line
  assert "Conflicts with the automatic rules" not in desc
  assert ("Did not add MoverHisFlip to chain A HIS 56 (protonation set by the "
    "user; no ion lock-down)") in desc
  assert "SingleHydrogenRotator at chain A ASP 99 HD2" in desc
  assert "SingleHydrogenRotator at chain A GLU 143 HE2" in desc
  assert h_names(out, 56).count("HE2") == 1 and "HD1" not in h_names(out, 56)
  # control: no entries, no lines
  out0, desc0, _ = run("tst_reduce2_protonation_2", frag_str)
  assert "Protonation states" not in desc0
  assert sorted(h_names(out0, 56)) != h_names(out, 56)

def test_his_with_flips():
  # add_flip_movers=True: the His keeps the set state and is not flipped.
  inp = atoms_of(frag_str, 56)
  for state, ring_h in (("his_hd1", ["HD1"]), ("his_he2", ["HE2"]),
                        ("his_both", ["HD1", "HE2"])):
    out, desc, _ = run("tst_reduce2_protonation_3", frag_str,
      [("resseq 56", state)], ["add_flip_movers=True"])
    a = atoms_of(out, 56)
    assert [n for n in h_names(out, 56) if n in ("HD1", "HE2")] == ring_h, (
      state, h_names(out, 56))
    for n in ("ND1", "CD2", "CE1", "NE2"):
      assert (matrix.col(a[n].xyz) - matrix.col(inp[n].xyz)).length() < 1.e-3
    assert "HisFlip at chain A HIS 56" not in desc
    assert "HisPlace at chain A HIS 56" not in desc
    assert "protonation set by the user" in desc

def kin_ring_h(file_name, resseq):
  """Ring H names of His A resseq in a flipkin file."""
  import re
  with open(file_name) as f:
    text = f.read()
  pattern = r"\{\s*(\S+)\s+his A\s+%d\s" % resseq
  return sorted(set(n.upper() for n in re.findall(pattern, text)
    if n.upper() in ("HD1", "HE2")))

def test_flipkins_keep_state():
  # add_flip_movers=True with flipkins: His 56 has a state, His 256 not (it
  # makes the His flipkin). Both flipkin files carry the state, and in every
  # flipkin pass His 56 gets no Mover and no ion lock. his_both: the
  # comparison pass (no Optimizer: the main result) lists no Mover there.
  for state, ring_h in (("his_he2", ["HE2"]), ("his_both", ["HD1", "HE2"])):
    prefix = "tst_reduce2_protonation_25_" + state
    if not os.path.isdir(prefix): os.mkdir(prefix)
    extra = ["add_flip_movers=True", "output.flipkin_directory=" + prefix]
    if state == "his_both":
      extra.append("comparison_file=" + prefix + ".pdb")
    out, desc, log = run(prefix, frag_two_his_str(), [("resseq 56", state)],
      extra)
    assert [n for n in h_names(out, 56) if n in ("HD1", "HE2")] == ring_h
    assert "HisFlip at chain A HIS 256" in desc
    for kin in ("flipnq", "fliphis"):
      name = os.path.join(prefix, "%sFH-%s.kin" % (prefix, kin))
      assert os.path.isfile(name), name
      assert kin_ring_h(name, 56) == ring_h, (state, kin, kin_ring_h(name, 56))
    passes = log.split("Results of optimization:")[1:]
    assert len(passes) == 4, len(passes)
    for p in passes:
      assert ("Did not add MoverHisFlip to chain A HIS 56 (protonation set by "
        "the user; no ion lock-down)") in p, state
      assert "MoverHisFlip 1 to chain A HIS 56\n" not in p, state
    if state == "his_both":
      with open(prefix + "_out.pdb_comparison.csv") as f:
        rows = f.read().splitlines()
      assert len(rows) > 1 and not [r for r in rows if "HIS 56" in r], rows

def test_his_next_to_zn():
  # The state is kept with the Zn at NE2 and the conflict is listed; without
  # a state the ion lock takes HE2 away (control).
  out, desc, _ = run("tst_reduce2_protonation_4", his_zn_str(),
    [("resseq 354", "his_he2")])
  assert h_names(out, 354).count("HE2") == 1, h_names(out, 354)
  assert "HD1" not in h_names(out, 354)
  assert ("Conflicts with the automatic rules (the state is kept):\n"
    "  HIS A 354: HE2 only, as set; NE2 2.05 A from ZN A 397\n") in desc, desc
  assert "ionic bond" not in desc
  out0, desc0, _ = run("tst_reduce2_protonation_5", his_zn_str())
  assert "HE2" not in h_names(out0, 354), h_names(out0, 354)
  assert "ionic bond to ZN" in desc0

def test_his_and_metal_in_different_conformers():
  # His in conformer A, the Zn in B: not in one conformer, no conflict. The
  # control with both blank is test_his_next_to_zn.
  out, desc, _ = run("tst_reduce2_protonation_21", his_zn_altlocs_str(),
    [("resseq 354", "his_he2")])
  assert "Conflicts with the automatic rules" not in desc, desc
  assert "HIS A 354: HE2 only (his_he2, conformers A): HE2 added, HD1 left out" \
    in desc, desc

def test_presence_per_conformer():
  # A complete, B truncated after CG: HE2 only in A, and the report says so
  out, desc, _ = run("tst_reduce2_protonation_22", frag_his_split_str(),
    [("resseq 56", "his_he2")])
  assert ("  HIS A 56: HE2 only (his_he2, conformers A B): HE2 added in A; HE2 "
    "not placed in B, HD1 left out\n") in desc, desc
  h = hierarchy(out)
  he2 = [a.parent().altloc for a in h.atoms()
         if a.parent().parent().resseq_as_int() == 56 and a.name.strip() == "HE2"]
  assert he2 == ["A"], he2

def test_his_with_glycan_link_split():
  # placement shares H between covalent partners for the split NAG-Asn (its
  # link_rule) while a His has a state: both work (the state's table had
  # taken the rule's name)
  out, desc, _ = run("tst_reduce2_protonation_23", his_56_glycan_str(),
    [("resseq 56", "his_he2")], ["approach=place"])
  assert [n for n in h_names(out, 56) if n in ("HD1", "HE2")] == ["HE2"]
  assert "HIS A 56: HE2 only (his_he2): HE2 added, HD1 left out" in desc
  # the glycan's H as without the state (the split H1 per Asn conformer)
  out0, _, _ = run("tst_reduce2_protonation_24", his_56_glycan_str(), (),
    ["approach=place"])
  def nag_h(text):
    h = hierarchy(text)
    return sorted((a.name.strip(), a.parent().altloc) for a in h.atoms()
      if a.parent().resname == "NAG" and a.element.strip() == "H")
  assert nag_h(out) == nag_h(out0) and len(nag_h(out)) > 10, nag_h(out)

def test_his_next_to_symmetry_mate_metal():
  # 1XMK His A 354 set to HE2 only, the Cd in a mate (symmetry copies are on
  # by default): kept, the conflict names the operator; no ion lock.
  out, desc, _ = run("tst_reduce2_protonation_6", his_cd_str,
    [("chain A and resseq 354", "his_he2")])
  assert [n for n in h_names(out, 354) if n in ("HD1", "HE2")] == ["HE2"]
  lines = [l for l in desc.splitlines() if l.startswith("  HIS A 354: HE2 only, "
    "as set; NE2 ")]
  assert len(lines) == 1, desc
  assert "from CD A 397 (sym " in lines[0], lines[0]
  assert "ionic bond" not in desc
  # control: defaults, the lock leaves HD1
  out0, desc0, _ = run("tst_reduce2_protonation_7", his_cd_str)
  assert [n for n in h_names(out0, 354) if n in ("HD1", "HE2")] == ["HD1"]

def test_acid_chi_beyond_90():
  # protonated Asp/Glu named against IUPAC: the H is on the input's
  # OD2/OE2 and the O keep their names (A38). Control: without a state the
  # names are exchanged.
  swapped = frag_swapped_str()
  inp99, inp143 = atoms_of(swapped, 99), atoms_of(swapped, 143)
  out, desc, _ = run("tst_reduce2_protonation_8", swapped,
    [("resseq 99 or resseq 143", "protonated")])
  for resseq, inp, o, h in ((99, inp99, "OD2", "HD2"),
                            (143, inp143, "OE2", "HE2")):
    a = atoms_of(out, resseq)
    assert (matrix.col(a[o].xyz) - matrix.col(inp[o].xyz)).length() < 1.e-3, (
      resseq, o)
    assert a[h].distance(a[o]) < 1.1, (resseq, a[h].distance(a[o]))
  assert "ASP A 99: OD1/OD2" not in desc and "GLU A 143: OE1/OE2" not in desc
  out0, desc0, _ = run("tst_reduce2_protonation_9", swapped)
  assert "ASP A 99: OD1/OD2" in desc0 and "GLU A 143: OE1/OE2" in desc0
  # an unrelated restraint file (B23: its processing does not flip): the H on
  # the input's OD2/OE2 coordinates as well, for place and add
  cif = libtbx.env.find_in_repositories(
    relative_path="chem_data/geostd/l/data_LYS.cif", test=os.path.isfile)
  if cif is None:
    print("Skipping the restraint-file part: data_LYS.cif not found")
    return
  for approach in ("place", "add"):
    out, desc, _ = run("tst_reduce2_protonation_9c", swapped,
      [("resseq 99 or resseq 143", "protonated")], ["approach=" + approach, cif])
    for resseq, inp, o, h in ((99, inp99, "OD2", "HD2"),
                              (143, inp143, "OE2", "HE2")):
      a = atoms_of(out, resseq)
      assert (matrix.col(a[o].xyz) - matrix.col(inp[o].xyz)).length() < 1.e-3, (
        approach, resseq, o)
      assert a[h].distance(a[o]) < 1.1, (approach, resseq)

def test_sorry():
  # Each entry error stops with a Sorry naming the entry.
  e = sorry("tst_reduce2_protonation_10", frag_str, [("resseq 999", "his_hd1")])
  assert e == ("protonation entry 1 (selection=resseq 999, state=his_hd1): the "
    "selection matches no atoms."), e
  e = sorry("tst_reduce2_protonation_11", frag_str,
    [("resseq 104", "protonated")])
  assert e == ("protonation entry 1 (selection=resseq 104, state=protonated): "
    "state protonated does not fit LYS A 104 (it is for ASP, GLU)."), e
  e = sorry("tst_reduce2_protonation_12", frag_str,
    [("resname HIS", "his_hd1"), ("resseq 56 and name CA", "his_he2")])
  assert e == ("protonation entry 2 (selection=resseq 56 and name CA, "
    "state=his_he2): HIS A 56 is also in protonation entry 1."), e
  e = sorry("tst_reduce2_protonation_13", frag_str, [("resseq 56", "his_hd1")],
    ["approach=optimize"])
  assert e == ("protonation entries need approach=add or place "
    "(approach=optimize): protonation entry 1 (selection=resseq 56, "
    "state=his_hd1). To keep the input H and set the state, use approach=add "
    "keep_existing_H=True."), e
  # supplied H that contradict the state, kept with keep_existing_H=True
  lines = frag_str.splitlines()
  last = max(i for i, l in enumerate(lines) if l[17:26] == "HIS A  56")
  with_he2 = "\n".join(lines[:last + 1] + ["ATOM    999  HE2 HIS A  56    "
    "%8.3f%8.3f%8.3f  1.00  5.00           H" % (8.9, 1.7, -3.4)] + lines[last + 1:]) + "\n"
  e = sorry("tst_reduce2_protonation_14", with_he2, [("resseq 56", "his_hd1")],
    ["keep_existing_H=True"])
  assert e == ("protonation entry 1 (selection=resseq 56, state=his_hd1): HIS A "
    "56 has supplied HE2, which contradicts state his_hd1 (the supplied H are "
    "kept: keep_existing_H=True). Remove it or change the state."), e
  # control: without keep_existing_H the supplied H are dropped
  out, _, _ = run("tst_reduce2_protonation_15", with_he2,
    [("resseq 56", "his_hd1")])
  assert "HE2" not in h_names(out, 56)

ALL_STATES = [("resseq 56", "his_he2"), ("resseq 99 or resseq 143", "protonated"),
  ("resseq 57 or resseq 61", "deprotonated")]

def without(text, pairs):
  """text without the atoms {(resseq, name)}."""
  return "\n".join(l for l in text.splitlines()
    if not (l.startswith("ATOM") and (int(l[22:26]), l[12:16].strip()) in pairs)
    ) + "\n"

def h_counts(text, resseq):
  """{(name, altloc): count} of the H/D of a residue."""
  result = {}
  h = hierarchy(text)
  for a in h.atoms():
    if a.parent().parent().resseq_as_int() == resseq and \
        a.element.strip() in ("H", "D"):
      k = (a.name.strip(), a.parent().altloc)
      result[k] = result.get(k, 0) + 1
  return result

def test_supplied_h():
  # keep_existing_H=True with supplied H: those that agree with the state are
  # kept (not added twice), missing state H are added. Supplied D: the line
  # names them. (Contradicting H: test_sorry.)
  full, _, _ = run("tst_reduce2_protonation_26", frag_str, ALL_STATES)
  out, desc, _ = run("tst_reduce2_protonation_27", full, ALL_STATES,
    ["keep_existing_H=True"])
  for resseq in (56, 57, 61, 99, 143):
    assert h_counts(out, resseq) == h_counts(full, resseq), resseq
  he2_in, he2_out = atoms_of(full, 56)["HE2"], atoms_of(out, 56)["HE2"]
  assert he2_in.distance(he2_out) < 1.e-3
  for line in ["  HIS A 56: HE2 only (his_he2): HE2 kept (supplied), HD1 left out",
               "  ASP A 99: protonated on OD2 (protonated): HD2 kept (supplied)",
               "  GLU A 143: protonated on OE2 (protonated): HE2 kept (supplied)"]:
    assert line + "\n" in desc, line
  # HE2 of His 56 and HD2 of Asp 99 missing: added, once
  out, desc, _ = run("tst_reduce2_protonation_28",
    without(full, {(56, "HE2"), (99, "HD2")}), ALL_STATES,
    ["keep_existing_H=True"])
  assert h_counts(out, 56) == h_counts(full, 56)
  assert h_counts(out, 99) == h_counts(full, 99)
  for line in ["  HIS A 56: HE2 only (his_he2): HE2 added, HD1 left out",
               "  ASP A 99: protonated on OD2 (protonated): HD2 added",
               "  GLU A 143: protonated on OE2 (protonated): HE2 kept (supplied)"]:
    assert line + "\n" in desc, line
  # all D (neutron_option=all_d): kept and named; a contradicting D stops
  full_d, _, _ = run("tst_reduce2_protonation_29", frag_str, ALL_STATES,
    ["neutron_option=all_d"])
  out, desc, _ = run("tst_reduce2_protonation_30", full_d, ALL_STATES,
    ["keep_existing_H=True", "neutron_option=all_d"])
  assert h_counts(out, 56) == h_counts(full_d, 56), h_counts(out, 56)
  assert "  HIS A 56: HE2 only (his_he2): DE2 kept (supplied), HD1 left out\n" \
    in desc, desc
  # H/D pairs (hd_and_h, written as altlocs A/B of the H): merged into one
  # site, kept, one phrase; A/B are layout letters, no conformers (B16)
  full_hd, _, _ = run("tst_reduce2_protonation_36", frag_str, ALL_STATES,
    ["neutron_option=hd_and_h"])
  _, desc, _ = run("tst_reduce2_protonation_37", full_hd, ALL_STATES,
    ["keep_existing_H=True", "neutron_option=hd_and_h"])
  assert ("  HIS A 56: HE2 only (his_he2): HE2/DE2 kept (supplied), HD1 left "
    "out\n") in desc, desc
  e = sorry("tst_reduce2_protonation_31", full_d, [("resseq 56", "his_hd1")],
    ["keep_existing_H=True", "neutron_option=all_d"])
  assert e == ("protonation entry 1 (selection=resseq 56, state=his_hd1): HIS A "
    "56 has supplied DE2, which contradicts state his_hd1 (the supplied H are "
    "kept: keep_existing_H=True). Remove it or change the state."), e

def test_supplied_half_pair():
  # B16 (I1 step B-2 review): hd_and_h, His 56 with DD1 alone (its HD1
  # removed) and a full HE2/DE2 pair: the pair is one site, named HE2/DE2;
  # the lone DD1 stays in its H/D-only group B, the only letter left
  full, _, _ = run("tst_reduce2_protonation_38", frag_str,
    [("resseq 56", "his_both")], ["neutron_option=hd_and_h"])
  half = without(full, {(56, "HD1")})
  assert h_counts(half, 56).get(("DD1", "B")) == 1
  _, desc, _ = run("tst_reduce2_protonation_39", half,
    [("resseq 56", "his_both")], ["keep_existing_H=True",
    "neutron_option=hd_and_h"])
  assert ("  HIS A 56: HD1 and HE2 (his_both, conformers B): DD1 kept "
    "(supplied), HE2/DE2 kept (supplied)\n") in desc, desc

def test_supplied_h_only_conformers():
  # B16 A-2 review, control: no D (no normalization), blank heavy atoms,
  # supplied HE2 in A and HB3 in B: the H-only altlocs are the conformers
  plain, _, _ = run("tst_reduce2_protonation_40", frag_str,
    [("resseq 56", "his_he2")])
  def alt(l):
    if l.startswith("ATOM") and int(l[22:26]) == 56:
      if l[12:16].strip() == "HE2": return l[:16] + "A" + l[17:]
      if l[12:16].strip() == "HB3": return l[:16] + "B" + l[17:]
    return l
  supplied = "\n".join(alt(l) for l in plain.splitlines()) + "\n"
  _, desc, _ = run("tst_reduce2_protonation_41", supplied,
    [("resseq 56", "his_he2")], ["keep_existing_H=True"])
  assert ("  HIS A 56: HE2 only (his_he2, conformers A B): HE2 kept (supplied) "
    "in A; HE2 added in B, HD1 left out\n") in desc, desc
  # B16 C-2/C-3 reviews: input D is normalized, but His 56's H-only altlocs
  # are no H/D layout and stay conformers. The D in another residue (Ala 55
  # DA); in His 56 itself, its amide H as a coincident H/D pair in C/D (0.5
  # each), also with HE2 A / HB3 B at 0.6/0.4
  def da(l):
    if l.startswith("ATOM") and int(l[22:26]) == 55 and \
       l[12:16].strip() == "HA":
      return l[:12] + " DA " + l[16:76] + " D" + l[78:]
    return l
  def amide_pair(l):
    if l.startswith("ATOM") and int(l[22:26]) == 56 and \
       l[12:16].strip() == "H":
      return (l[:16] + "C" + l[17:54] + "  0.50" + l[60:] + "\n" +
              l[:12] + " D  D" + l[17:54] + "  0.50" + l[60:76] + " D" +
              l[78:])
    return l
  def occ_64(l):
    if l.startswith("ATOM") and int(l[22:26]) == 56:
      if l[12:17] == " HE2A": return l[:54] + "  0.60" + l[60:]
      if l[12:17] == " HB3B": return l[:54] + "  0.40" + l[60:]
    return l
  with_d = "\n".join(da(l) for l in supplied.splitlines()) + "\n"
  assert h_counts(with_d, 55).get(("DA", "")) == 1
  with_pair = "\n".join(amide_pair(l) for l in supplied.splitlines()) + "\n"
  assert h_counts(with_pair, 56).get(("D", "D")) == 1
  with_pair_64 = "\n".join(occ_64(l) for l in with_pair.splitlines()) + "\n"
  assert with_pair_64.count("  0.60") == 1 and with_pair_64.count("  0.40") == 1
  for prefix, text in (("tst_reduce2_protonation_44", with_d),
                       ("tst_reduce2_protonation_45", with_pair),
                       ("tst_reduce2_protonation_46", with_pair_64)):
    model, desc, _ = run(prefix, text, [("resseq 56", "his_he2")],
      ["keep_existing_H=True", "neutron_option=hd_and_h"])
    assert ("  HIS A 56: HE2 only (his_he2, conformers A B): HE2 kept "
      "(supplied) in A; HE2 added in B, HD1 left out\n") in desc, (prefix,
      desc)
    # both conformers' HE2 sit on one site (blank ring): one HE2/DE2 pair
    sites = {"HE2": [], "DE2": []}
    out = hierarchy(model)   # the atoms keep their parents while in use
    for a in out.atoms():
      if a.parent().parent().resseq_as_int() == 56 and \
         a.name.strip() in sites:
        sites[a.name.strip()].append(a.xyz)
    assert len(sites["HE2"]) == 1 and sites["HE2"] == sites["DE2"], (prefix,
      sites)

def test_supplied_pair_relabelled():
  # B16 A-2 review (P3), repair (d): His 56 CG split A/B (0.6/0.4), the
  # supplied H/D pairs in A/B at 0.5/0.5 then sit in conformer groups; the
  # repair moves them to two free letters and the normalization merges them:
  # one supplied pair, in both conformers
  full_hd, _, _ = run("tst_reduce2_protonation_42", frag_str,
    [("resseq 56", "his_he2")], ["neutron_option=hd_and_h"])
  lines = []
  for l in full_hd.splitlines():
    if l.startswith("ATOM") and int(l[22:26]) == 56 and \
       l[12:16].strip() == "CG":
      lines.append(l[:16] + "A" + l[17:54] + "  0.60" + l[60:])
      lines.append(l[:16] + "B" + l[17:54] + "  0.40" + l[60:])
    else:
      lines.append(l)
  split = "\n".join(lines) + "\n"
  assert h_counts(split, 56).get(("DE2", "B")) == 1
  _, desc, log = run("tst_reduce2_protonation_43", split,
    [("resseq 56", "his_he2")], ["keep_existing_H=True",
    "neutron_option=hd_and_h"])
  assert "H/D pair in conformer groups, moved to H/D-only groups" in desc
  assert ("  HIS A 56: HE2 only (his_he2, conformers A B): HE2/DE2 kept "
    "(supplied), HD1 left out\n") in desc, desc

def test_acid_out_of_plane_line():
  # A46, report only: the acid H on frag_swapped's Asp 99 ends 67 deg from
  # the carboxyl plane, and its line gives the angle (no cause: none is
  # checked). frag with preference_magnitude=0: about 89 deg. Glu 143, and
  # both in frag with the default preference (control), end in the plane.
  _, desc, _ = run("tst_reduce2_protonation_32", frag_swapped_str(),
    [("resseq 99 or resseq 143", "protonated")])
  assert ("  ASP A 99: protonated on OD2 (protonated): HD2 added, HD2 67 deg "
    "out of the carboxyl plane\n") in desc, desc
  assert "  GLU A 143: protonated on OE2 (protonated): HE2 added\n" in desc
  _, desc0, _ = run("tst_reduce2_protonation_33", frag_str,
    [("resseq 99 or resseq 143", "protonated")])
  assert "carboxyl plane" not in desc0
  assert "  ASP A 99: protonated on OD2 (protonated): HD2 added\n" in desc0
  _, desc1, _ = run("tst_reduce2_protonation_34", frag_str,
    [("resseq 99", "protonated")], ["preference_magnitude=0"])
  assert ("  ASP A 99: protonated on OD2 (protonated): HD2 added, HD2 89 deg "
    "out of the carboxyl plane\n") in desc1, desc1

def test_supplied_h_per_conformer():
  # His 56 split A/B, HE2 supplied in A only, his_he2, keep_existing_H=True:
  # one phrase per conformer. B complete: HE2 added there. B truncated after
  # CG: HE2 not placed there.
  for complete, b_phrase in ((True, "HE2 added in B"),
                             (False, "HE2 not placed in B")):
    prefix = "tst_reduce2_protonation_35_%s" % complete
    full, _, _ = run(prefix + "_full", frag_his_split_str(complete),
      [("resseq 56", "his_he2")])
    supplied = "\n".join(l for l in full.splitlines() if not (
      l.startswith("ATOM") and l[17:26] == "HIS A  56" and
      l[12:16].strip() == "HE2" and l[16] == "B")) + "\n"
    h = hierarchy(supplied)
    assert [a.parent().altloc for a in h.atoms()
      if a.name.strip() == "HE2" and a.parent().resname == "HIS"] == ["A"]
    out, desc, _ = run(prefix, supplied, [("resseq 56", "his_he2")],
      ["keep_existing_H=True"])
    assert ("  HIS A 56: HE2 only (his_he2, conformers A B): HE2 kept (supplied) "
      "in A; %s, HD1 left out\n" % b_phrase) in desc, (complete, desc)

def test_neutron_options():
  # The new acid H is exchangeable: H/D pair with hd_and_h, D with all_d.
  out, _, _ = run("tst_reduce2_protonation_16", frag_str,
    [("resseq 99 or resseq 143", "protonated")], ["neutron_option=hd_and_h"])
  assert "HD2" in h_names_all(out, 99) and "DD2" in h_names_all(out, 99), \
    h_names_all(out, 99)
  assert "HE2" in h_names_all(out, 143) and "DE2" in h_names_all(out, 143)
  out, _, _ = run("tst_reduce2_protonation_17", frag_str,
    [("resseq 99 or resseq 143", "protonated")], ["neutron_option=all_d"])
  assert "DD2" in h_names_all(out, 99) and "HD2" not in h_names_all(out, 99)
  assert "DE2" in h_names_all(out, 143) and "HE2" not in h_names_all(out, 143)

def h_names_all(text, resseq):
  """H/D names of a residue in all its conformers."""
  h = hierarchy(text)
  return sorted(set(a.name.strip() for a in h.atoms()
    if a.parent().parent().resseq_as_int() == resseq
    and a.element.strip() in ("H", "D")))

def test_qmready_charge():
  # The qmready total charge follows the states: +1 per protonated acid,
  # -1 for a His with one ring H instead of both, -1 per deprotonated
  # Cys or Tyr.
  def charge(prefix, entries):
    _, _, log = run(prefix, frag_str, entries, ["qmready=True"])
    lines = [l for l in log.splitlines() if l.startswith("Total charge:")]
    assert len(lines) == 1, lines
    return int(lines[0].split()[-1])
  q0 = charge("tst_reduce2_protonation_18", [])
  q1 = charge("tst_reduce2_protonation_19",
    [("resseq 99 or resseq 143", "protonated")])
  assert q1 == q0 + 2, (q0, q1)
  q2 = charge("tst_reduce2_protonation_20", [("resseq 56", "his_hd1"),
    ("resseq 57 or resseq 61", "deprotonated")])
  assert q2 == q0 - 3, (q0, q2)

if __name__ == "__main__":
  test_states_placement()
  test_report_and_optimizer()
  test_his_with_flips()
  test_flipkins_keep_state()
  test_his_next_to_zn()
  test_his_and_metal_in_different_conformers()
  test_presence_per_conformer()
  test_his_with_glycan_link_split()
  test_his_next_to_symmetry_mate_metal()
  test_acid_chi_beyond_90()
  test_sorry()
  test_supplied_h()
  test_supplied_half_pair()
  test_supplied_h_only_conformers()
  test_supplied_pair_relabelled()
  test_acid_out_of_plane_line()
  test_supplied_h_per_conformer()
  test_neutron_options()
  test_qmready_charge()
  print("OK")
