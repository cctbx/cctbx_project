from __future__ import absolute_import, division, print_function
import mmtbx.model
import iotbx.pdb
from libtbx.utils import null_out
from libtbx import group_args
import pickle

# Gln-Asn-Tyr, no H
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
"""

# Cys 2 and Cys 14 form disulfide; His 9 is doubly protonated.
# Fragment of 2EFZ, H atoms added with reduce2.
cys_his_cys_h = """\
CRYST1   27.938   22.382   30.741  90.00  90.00  90.00 P 1
ATOM      1  N   CYS A   2      17.575  12.598  22.054  1.00  0.00           N
ATOM      2  CA  CYS A   2      16.798  11.919  21.025  1.00  0.00           C
ATOM      3  C   CYS A   2      15.621  12.795  20.605  1.00  0.00           C
ATOM      4  O   CYS A   2      15.258  13.724  21.326  1.00  0.00           O
ATOM      5  CB  CYS A   2      16.283  10.579  21.552  1.00  0.00           C
ATOM      6  SG  CYS A   2      16.060   9.307  20.297  1.00  0.00           S
ATOM      7  H   CYS A   2      17.728  13.428  21.887  1.00  0.00           H
ATOM      8  HA  CYS A   2      17.355  11.745  20.250  1.00  0.00           H
ATOM      9  HB3 CYS A   2      15.422  10.725  21.975  1.00  0.00           H
ATOM     10  HB2 CYS A   2      16.918  10.239  22.202  1.00  0.00           H
ATOM     12  N   HIS A   9      11.030  12.441  12.955  1.00  0.00           N
ATOM     13  CA  HIS A   9      10.744  12.272  11.547  1.00  0.00           C
ATOM     14  C   HIS A   9      11.998  11.616  10.970  1.00  0.00           C
ATOM     15  O   HIS A   9      12.993  11.480  11.684  1.00  0.00           O
ATOM     16  CB  HIS A   9      10.493  13.645  10.914  1.00  0.00           C
ATOM     17  CG  HIS A   9       9.210  14.294  11.370  1.00  0.00           C
ATOM     18  ND1 HIS A   9       8.969  14.837  12.612  1.00  0.00           N
ATOM     19  CD2 HIS A   9       8.059  14.421  10.637  1.00  0.00           C
ATOM     20  CE1 HIS A   9       7.699  15.273  12.625  1.00  0.00           C
ATOM     21  NE2 HIS A   9       7.105  15.050  11.441  1.00  0.00           N
ATOM     22  H   HIS A   9      10.412  12.147  13.475  1.00  0.00           H
ATOM     23  HA  HIS A   9       9.955  11.739  11.363  1.00  0.00           H
ATOM     24  HB2 HIS A   9      11.225  14.236  11.150  1.00  0.00           H
ATOM     25  HB3 HIS A   9      10.447  13.541   9.951  1.00  0.00           H
ATOM     26  HD1 HIS A   9       9.530  14.888  13.262  1.00  0.00           H
ATOM     27  HD2 HIS A   9       7.935  14.138   9.760  1.00  0.00           H
ATOM     28  HE1 HIS A   9       7.286  15.675  13.355  1.00  0.00           H
ATOM     29  HE2 HIS A   9       6.299  15.254  11.222  1.00  0.00           H
ATOM     30  N   CYS A  14      18.239   9.363  15.780  1.00  0.00           N
ATOM     31  CA  CYS A  14      18.026   8.685  17.055  1.00  0.00           C
ATOM     32  C   CYS A  14      16.614   8.108  17.168  1.00  0.00           C
ATOM     33  O   CYS A  14      16.451   6.911  17.375  1.00  0.00           O
ATOM     34  CB  CYS A  14      18.346   9.614  18.227  1.00  0.00           C
ATOM     35  SG  CYS A  14      18.031   8.821  19.820  1.00  0.00           S
ATOM     36  H   CYS A  14      18.875   9.028  15.307  1.00  0.00           H
ATOM     37  HA  CYS A  14      18.636   7.933  17.111  1.00  0.00           H
ATOM     38  HB3 CYS A  14      17.790  10.407  18.167  1.00  0.00           H
ATOM     39  HB2 CYS A  14      19.283   9.862  18.191  1.00  0.00           H
HETATM   40  O   HOH A 101       5.000   5.000   5.000  1.00  0.00           O
"""

def get_model(pdb_str):
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  return mmtbx.model.manager(model_input=pdb_inp, log=null_out())

def remove(pdb_str, resname, names):
  return "\n".join([l for l in pdb_str.splitlines()
    if not (l[17:20]==resname and l[12:16].strip() in names)]) + "\n"

def exercise_heavy_atoms():
  m = get_model(tripeptide)
  assert m.is_atom_complete()
  assert m.missing_atoms_by_residue() == {}
  # no H in the model: complete in non-H sense only
  assert not m.is_atom_complete(include_hydrogens=True)
  # missing side chain atom
  m = get_model(remove(tripeptide, "GLN", ["CD"]))
  assert not m.is_atom_complete()
  assert m.missing_atoms_by_residue() == {" GLN A   1 ": ["CD"]}
  # truncated side chain
  m = get_model(remove(tripeptide, "TYR", ["CG","CD1","CD2","CE1","CE2","CZ"]))
  assert m.missing_atoms_by_residue() == {
    " TYR A   3 ": ["CD1","CD2","CE1","CE2","CG","CZ"]}

def exercise_alt_confs():
  # TYR OH in A and B: complete (and no crash on alt. conformations)
  lines = []
  for l in tripeptide.splitlines():
    if l[17:20]=="TYR" and l[12:16]==" OH ":
      lines.append(l[:16]+"A"+l[17:54]+"  0.50"+l[60:])
      lines.append(l[:16]+"B"+l[17:54]+"  0.50"+l[60:])
    else:
      lines.append(l)
  m = get_model("\n".join(lines)+"\n")
  assert m.is_atom_complete()
  # CZ in A and B, OH only in A: conformer B is incomplete
  lines = []
  for l in tripeptide.splitlines():
    if l[17:20]=="TYR" and l[12:16] in [" CZ ", " OH "]:
      lines.append(l[:16]+"A"+l[17:54]+"  0.50"+l[60:])
      if l[12:16]==" CZ ":
        lines.append(l[:16]+"B"+l[17:54]+"  0.50"+l[60:])
    else:
      lines.append(l)
  m = get_model("\n".join(lines)+"\n")
  assert m.missing_atoms_by_residue() == {"BTYR A   3 ": ["OH"]}, \
    m.missing_atoms_by_residue()

def exercise_hydrogens():
  m = get_model(cys_his_cys_h)
  # SG-SG bonded CYS do not need HG; waters do not need H
  assert m.is_atom_complete(include_hydrogens=True)
  assert m.is_atom_complete()
  # HIS: either of HD1 or HE2 is enough
  for name in ["HD1", "HE2"]:
    m = get_model(remove(cys_his_cys_h, "HIS", [name]))
    assert m.is_atom_complete(include_hydrogens=True)
  m = get_model(remove(cys_his_cys_h, "HIS", ["HD1", "HE2"]))
  assert m.missing_atoms_by_residue(include_hydrogens=True) == {
    " HIS A   9 ": ["HD1", "HE2"]}
  # missing H, but complete in non-H sense
  m = get_model(remove(cys_his_cys_h, "HIS", ["HA"]))
  assert m.missing_atoms_by_residue(include_hydrogens=True) == {
    " HIS A   9 ": ["HA"]}
  assert m.is_atom_complete()
  # CYS without disulfide needs HG
  m = get_model(cys_his_cys_h)
  for a in m.get_hierarchy().atoms():
    if a.parent().parent().resseq.strip()=="14":
      a.set_xyz([a.xyz[0]+10, a.xyz[1], a.xyz[2]])
  assert m.missing_atoms_by_residue(include_hydrogens=True) == {
    " CYS A   2 ": ["HG"], " CYS A  14 ": ["HG"]}

def exercise_current_hierarchy_used():
  m = get_model(tripeptide)
  m.process(make_restraints=True)
  n_atoms = m.size()
  assert m.is_atom_complete()
  # model is not changed
  assert m.size() == n_atoms
  assert m.restraints_manager_available()
  # change hierarchy in place: result reflects the change
  ag = [ag for ag in m.get_hierarchy().atom_groups() if ag.resname=="GLN"][0]
  ag.remove_atom(ag.get_atom("CD"))
  assert m.missing_atoms_by_residue() == {" GLN A   1 ": ["CD"]}

def exercise_get_missing_atoms_not_stale():
  def heavy(m):
    return [k for k, v in m.get_missing_atoms().items()
            if v.get("missing", {}).get("heavy")]
  m = get_model(tripeptide)
  m.process()
  assert heavy(m) == []
  ag = [ag for ag in m.get_hierarchy().atom_groups() if ag.resname=="GLN"][0]
  ag.remove_atom(ag.get_atom("CD"))
  m.reset_after_changing_hierarchy()
  assert heavy(m) == [" GLN A   1 "] # lazy re-process after reset
  # re-process with atom put back: stale result is not kept
  m = get_model(remove(tripeptide, "GLN", ["CD"]))
  m.process()
  assert heavy(m) == [" GLN A   1 "]
  ag = [ag for ag in m.get_hierarchy().atom_groups() if ag.resname=="GLN"][0]
  cd = iotbx.pdb.hierarchy.atom()
  cd.set_name(" CD ")
  cd.set_element(" C")
  cd.set_xyz((16.381, 8.233, 15.935))
  ag.append_atom(cd)
  m.process()
  assert heavy(m) == []

MARKER = {" XXX X   1 ": ["MARKER"]}

def plant_marker(m):
  """ Replace cached records: getting MARKER back means cache was used """
  assert m._missing_atoms_records is not None
  m._missing_atoms_records = [group_args(
    key=" XXX X   1 ", legacy_key=" XXX X   1 ", resname="XXX",
    classification="other", incomplete=None,
    missing={"heavy": {"MARKER": None}}, sg_bonded=False)]

def exercise_cache():
  # process() fills cache; same result as computed from scratch
  for pdb_str in [tripeptide, remove(tripeptide, "GLN", ["CD"]),
                  cys_his_cys_h, remove(cys_his_cys_h, "HIS", ["HD1","HE2"])]:
    m1 = get_model(pdb_str)
    m1.process()
    m2 = get_model(pdb_str)
    assert m2._missing_atoms_records is None
    for h in [False, True]:
      assert m1.missing_atoms_by_residue(include_hydrogens=h) == \
             m2.missing_atoms_by_residue(include_hydrogens=h)
  # cache is used after process() and after first (fallback) call
  for do_process in [True, False]:
    m = get_model(tripeptide)
    if do_process: m.process()
    else:          m.missing_atoms_by_residue()
    plant_marker(m)
    assert m.missing_atoms_by_residue() == MARKER
    assert m.missing_atoms_by_residue(include_hydrogens=True) == MARKER
    assert not m.is_atom_complete()
  # returned dictionary is a copy
  m = get_model(remove(tripeptide, "GLN", ["CD"]))
  m.missing_atoms_by_residue()[" GLN A   1 "].append("XX")
  assert m.missing_atoms_by_residue() == {" GLN A   1 ": ["CD"]}
  # coordinate change does not invalidate cache
  m = get_model(tripeptide)
  m.process()
  plant_marker(m)
  m.set_sites_cart(m.get_sites_cart() + (1, 0, 0))
  assert m.missing_atoms_by_residue() == MARKER
  # reset_after_changing_hierarchy() clears cache
  m.reset_after_changing_hierarchy()
  assert m._missing_atoms_records is None
  assert m.is_atom_complete()

def exercise_fingerprint():
  """ In-place hierarchy changes invalidate cache """
  def gln(m):
    return [ag for ag in m.get_hierarchy().atom_groups() if ag.resname=="GLN"][0]
  def remove_atom(m):
    ag = gln(m)
    ag.remove_atom(ag.get_atom("CD"))
  def rename_atom(m):
    gln(m).get_atom("CD").set_name(" CX ")
  def rename_residue(m): # mutation
    gln(m).resname = "GLU"
  def change_altloc(m):
    gln(m).altloc = "A"
  def renumber(m):
    gln(m).parent().resseq = "  10"
  def add_atom(m):
    a = iotbx.pdb.hierarchy.atom()
    a.set_name(" OXT")
    a.set_element(" O")
    [ag for ag in m.get_hierarchy().atom_groups()
       if ag.resname=="TYR"][0].append_atom(a)
  for change in [remove_atom, rename_atom, rename_residue, change_altloc,
                 renumber, add_atom]:
    m = get_model(tripeptide)
    m.process()
    plant_marker(m)
    change(m)
    result = m.missing_atoms_by_residue()
    assert result != MARKER, change.__name__
  # stale result is not returned
  m = get_model(tripeptide)
  m.process()
  remove_atom(m)
  assert m.missing_atoms_by_residue() == {" GLN A   1 ": ["CD"]}

def exercise_ncs():
  """
  Identical chains: pdb_interpretation interprets only NCS master copy
  (use_ncs_to_build_restraints), still all copies must be checked.
  """
  def two_copies(pdb_str):
    result = []
    for l in pdb_str.splitlines():
      if not l.startswith("ATOM"): continue
      result.append(l)
    for l in list(result):
      x, y, z = [float(l[30+8*i:38+8*i]) for i in range(3)]
      result.append(l[:21]+"B"+l[22:30]+"%8.3f%8.3f%8.3f"%(x, y+30, z)+l[54:])
    return "\n".join(result)+"\n"
  m = get_model(two_copies(remove(tripeptide, "GLN", ["CD"])))
  m.process(make_restraints=True)
  expected = {" GLN A   1 ": ["CD"], " GLN B   1 ": ["CD"]}
  assert m.missing_atoms_by_residue() == expected, m.missing_atoms_by_residue()
  m = get_model(two_copies(tripeptide))
  m.process(make_restraints=True)
  assert m.is_atom_complete()

def exercise_get_missing_atoms():
  """ get_missing_atoms(): raw monomer library view, same records as above """
  # format: incomplete residue, raw H (CYS HG, water H)
  m = get_model(remove(tripeptide, "TYR",
    ["CB","CG","CD1","CD2","CE1","CE2","CZ","OH"]))
  r = m.get_missing_atoms()
  assert r[" TYR A   3 "] == {"resname": "TYR", "class": "peptide",
    "incomplete": "backbone_only"}
  assert sorted(r[" GLN A   1 "]["missing"].keys()) == ["hydrogens"]
  r = get_model(cys_his_cys_h).get_missing_atoms()
  assert list(r[" CYS A   2 "]["missing"]["hydrogens"].keys()) == ["HG"]
  assert sorted(r[" HOH A 101 "]["missing"]["hydrogens"].keys()) == ["H1","H2"]
  # no side effects if not processed: no restraints made
  m = get_model(tripeptide)
  m.get_missing_atoms()
  assert not m.restraints_manager_available()
  # in-place change without reset: result is current
  m = get_model(tripeptide)
  m.process(make_restraints=True)
  ag = [ag for ag in m.get_hierarchy().atom_groups() if ag.resname=="GLN"][0]
  ag.remove_atom(ag.get_atom("CD"))
  assert list(m.get_missing_atoms()[" GLN A   1 "]["missing"]["heavy"].keys())\
    == ["CD"]
  # NCS copies are reported
  lines = [l for l in remove(tripeptide, "GLN", ["CD"]).splitlines()
           if l.startswith("ATOM")]
  for l in list(lines):
    y = float(l[38:46]) + 30
    lines.append(l[:21] + "B" + l[22:38] + "%8.3f" % y + l[46:])
  m = get_model("\n".join(lines)+"\n")
  m.process(make_restraints=True)
  heavy = [k for k, v in m.get_missing_atoms().items()
           if v.get("missing", {}).get("heavy")]
  assert sorted(heavy) == [" GLN A   1 ", " GLN B   1 "], heavy

def exercise_pickle():
  m = get_model(remove(tripeptide, "GLN", ["CD"]))
  m.process()
  m2 = pickle.loads(pickle.dumps(m))
  assert m2.missing_atoms_by_residue() == {" GLN A   1 ": ["CD"]}

def run():
  exercise_heavy_atoms()
  exercise_alt_confs()
  exercise_hydrogens()
  exercise_current_hierarchy_used()
  exercise_get_missing_atoms_not_stale()
  exercise_cache()
  exercise_fingerprint()
  exercise_ncs()
  exercise_get_missing_atoms()
  exercise_pickle()
  print("OK")

if (__name__ == "__main__"):
  run()
