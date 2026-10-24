from __future__ import absolute_import, division, print_function
import mmtbx.model
import iotbx.pdb
from scitbx import matrix
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out

# Gln-Asn-Tyr; the TYR OXT is the deposited one, kept for reference
tripeptide = """\
ATOM      1  N   GLN A   1      17.191   7.640  11.852  1.00 10.38           N
ATOM      2  CA  GLN A   1      18.307   8.180  12.621  1.00 11.39           C
ATOM      3  C   GLN A   1      19.631   7.587  12.153  1.00 11.52           C
ATOM      4  O   GLN A   1      19.805   6.365  12.035  1.00 12.05           O
ATOM      5  CB  GLN A   1      18.093   8.002  14.128  1.00 11.96           C
ATOM      6  CG  GLN A   1      16.866   8.769  14.628  1.00 10.81           C
ATOM      7  CD  GLN A   1      16.381   8.233  15.935  1.00 13.10           C
ATOM      8  OE1 GLN A   1      15.811   7.144  15.983  1.00 10.65           O
ATOM      9  NE2 GLN A   1      16.586   9.006  17.020  1.00 12.30           N
ATOM     10  N   ASN A   2      20.551   8.483  11.837  1.00 11.99           N
ATOM     11  CA  ASN A   2      21.868   8.129  11.299  1.00 12.30           C
ATOM     12  C   ASN A   2      22.891   8.580  12.305  1.00 13.40           C
ATOM     13  O   ASN A   2      23.256   9.762  12.355  1.00 13.92           O
ATOM     14  CB  ASN A   2      22.102   8.835   9.974  1.00 12.13           C
ATOM     15  CG  ASN A   2      20.998   8.554   8.984  1.00 12.77           C
ATOM     16  OD1 ASN A   2      20.835   7.423   8.532  1.00 14.27           O
ATOM     17  ND2 ASN A   2      20.232   9.566   8.660  1.00 10.07           N
ATOM     18  N   TYR A   3      23.329   7.636  13.128  1.00 14.70           N
ATOM     19  CA  TYR A   3      24.196   7.963  14.280  1.00 15.18           C
ATOM     20  C   TYR A   3      25.640   8.150  13.866  1.00 15.91           C
ATOM     21  O   TYR A   3      26.078   7.630  12.836  1.00 15.76           O
ATOM     22  CB  TYR A   3      24.098   6.884  15.350  1.00 15.35           C
ATOM     23  CG  TYR A   3      22.702   6.748  15.883  1.00 14.45           C
ATOM     24  CD1 TYR A   3      22.247   7.575  16.901  1.00 14.80           C
ATOM     25  CD2 TYR A   3      21.808   5.840  15.308  1.00 15.68           C
ATOM     26  CE1 TYR A   3      20.941   7.468  17.397  1.00 14.33           C
ATOM     27  CE2 TYR A   3      20.517   5.725  15.777  1.00 13.46           C
ATOM     28  CZ  TYR A   3      20.084   6.548  16.812  1.00 15.09           C
ATOM     29  OH  TYR A   3      18.803   6.408  17.272  1.00 14.39           O
ATOM     30  OXT TYR A   3      26.395   8.818  14.593  1.00 17.49           O
"""

def edit(lines, chain_id=None, shift=(0,0,0), keep=None):
  result = []
  for l in lines.splitlines():
    if keep is not None and not keep(l): continue
    if chain_id is not None: l = l[:21] + chain_id + l[22:]
    x, y, z = [float(l[30+8*i:38+8*i]) + shift[i] for i in range(3)]
    l = l[:30] + "%8.3f%8.3f%8.3f" % (x, y, z) + l[54:]
    result.append(l)
  return "\n".join(result) + "\n"

def is_oxt(l): return l[12:16] == " OXT"

def get_model(pdb_str):
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  return mmtbx.model.manager(model_input=pdb_inp, log=null_out())

def oxt_atoms(model):
  return [a for a in model.get_hierarchy().atoms() if a.name.strip()=="OXT"]

def check_geometry(model, skip_chain_ids=()):
  for oxt in oxt_atoms(model):
    if oxt.parent().parent().parent().id in skip_chain_ids: continue
    ag = oxt.parent()
    rg = ag.parent()
    names = {}
    for a in rg.atoms():
      if a.parent().altloc in ["", ag.altloc]: names[a.name.strip()] = a
    c, ca, o = [matrix.col(names[n].xyz) for n in ["C", "CA", "O"]]
    x = matrix.col(oxt.xyz)
    assert approx_equal((x-c).length(), 1.231, 1.e-3)
    assert 115 < (x-c).angle(o-c, deg=True) < 130
    assert 110 < (x-c).angle(ca-c, deg=True) < 130
    # in the CA-C-O plane
    assert approx_equal((x-c).dot((ca-c).cross(o-c).normalize()), 0, 1.e-3)

def exercise_chain_ends_and_breaks():
  """
  Chain A: complete tripeptide, followed by a water (same chain).
  Chain B: tripeptide whose residue 2 is deleted -> break after residue 1.
  Chain C: tripeptide that already has OXT.
  """
  pdb_str = (
    edit(tripeptide, keep=lambda l: not is_oxt(l)) +
    "HETATM   31  O   HOH A   4      30.000  10.000  10.000  1.00 20.00           O\n" +
    edit(tripeptide, chain_id="B", shift=(30,0,0),
      keep=lambda l: not is_oxt(l) and l[17:20]!="ASN") +
    edit(tripeptide, chain_id="C", shift=(0,30,0)))
  model = get_model(pdb_str)
  n_atoms = model.size()
  assert model.add_terminal_oxt() == 3
  assert model.size() == n_atoms + 3
  ids = sorted([(a.parent().parent().parent().id,
                 a.parent().parent().resseq.strip()) for a in oxt_atoms(model)])
  assert ids == [("A","3"), ("B","1"), ("B","3"), ("C","3")], ids
  check_geometry(model, skip_chain_ids=["C"]) # C has deposited OXT
  # compare with the deposited OXT
  deposited = matrix.col([float(l[30+8*i:38+8*i])
    for l in tripeptide.splitlines() if is_oxt(l) for i in range(3)])
  oxt_a = [a for a in oxt_atoms(model)
           if a.parent().parent().parent().id=="A"][0]
  assert (matrix.col(oxt_a.xyz) - deposited).length() < 0.2
  # new atoms inherit O properties; model is consistent
  assert approx_equal(oxt_a.b, 15.76)
  assert model.get_xray_structure().scatterers().size() == model.size()
  # second call adds nothing
  assert model.add_terminal_oxt() == 0
  # model with OXT can be processed (OXT is recognized by restraints)
  model.process(make_restraints=True)

def exercise_alt_confs():
  """
  C-terminal TYR with O in two alternative conformations; O of B is O of A
  rotated by 20 degrees around CA-C axis (C-O and CA-C-O are preserved).
  """
  lines = edit(tripeptide, keep=lambda l: not is_oxt(l)).splitlines()
  xyz = {}
  for l in lines:
    if l[17:20]=="TYR":
      xyz[l[12:16].strip()] = [float(l[30+8*i:38+8*i]) for i in range(3)]
  o_b = matrix.rotate_point_around_axis(
    axis_point_1=xyz["CA"], axis_point_2=xyz["C"], point=xyz["O"], angle=20,
    deg=True)
  result = []
  for l in lines:
    if l[17:20]=="TYR" and l[12:16]==" O  ":
      result.append(l[:16] + "A" + l[17:54] + "  0.50" + l[60:])
      result.append(l[:16] + "B" + l[17:30] + "%8.3f%8.3f%8.3f" % tuple(o_b) +
        "  0.50" + l[60:])
    else:
      result.append(l)
  model = get_model("\n".join(result)+"\n")
  assert model.add_terminal_oxt() == 2
  altlocs = sorted([a.parent().altloc for a in oxt_atoms(model)])
  assert altlocs == ["A", "B"], altlocs
  for a in oxt_atoms(model):
    assert approx_equal(a.occ, 0.5)
  check_geometry(model)

def exercise_missing_o():
  """ No OXT if O (or C, CA) is absent """
  model = get_model(edit(tripeptide,
    keep=lambda l: not is_oxt(l) and not (l[17:20]=="TYR" and l[12:16]==" O  ")))
  assert model.add_terminal_oxt() == 0

def exercise_bonded_c():
  """
  No OXT if C is bonded to another residue: an NH2 cap (next residue, not an
  amino acid; OXT was put 0.12 A from its N), or an atom of another chain
  (ring closure of a cyclic peptide, link). N sits at the deposited OXT site.
  """
  oxt = [l for l in tripeptide.splitlines() if is_oxt(l)][0]
  cap = "HETATM   31  N   NH2 A   4" + oxt[26:76] + " N\n"
  link = "HETATM   31  N   NH2 B   1" + oxt[26:76] + " N\n"
  for extra in [cap, link]:
    model = get_model(edit(tripeptide, keep=lambda l: not is_oxt(l)) + extra)
    assert model.add_terminal_oxt() == 0, extra

def run():
  exercise_chain_ends_and_breaks()
  exercise_alt_confs()
  exercise_missing_o()
  exercise_bonded_c()
  print("OK")

if (__name__ == "__main__"):
  run()
