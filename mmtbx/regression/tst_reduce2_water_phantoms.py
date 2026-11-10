"""reduce2 Optimizer: no Phantom Hydrogens on a water Oxygen that carries H or D
in the pass. Its own H/D are in the model; phantoms on top of them gave the
water up to four donors. Since #1183 the Optimizer logged that it skipped such
a water but still made the phantoms.

- HOH 1 (H1, H2) and HOH 7 (D1, D2): no phantoms; HOH 2 and HOH 8, bare and
  2.8 A from them: phantoms.
- HOH 10 (H1 only) and HOH 12 (D1 only): no phantoms, so the second donor is
  not represented; HOH 11 and HOH 13, bare and 2.8 A from them: phantoms.
- HOH 4, bare, bonded to MG 3 (link_all=True): phantoms. Only an H or D
  neighbour stops them.
- HOH 5, O blank, H1 and H2 in altloc A only: no phantoms in pass A, phantoms
  in pass B (CL 9 in A and B gives the B pass).
- P 1, a = 20 A: HOH 1 (H1, H2) at x = 0.5 and HOH 2 at x = 17.7, 2.8 A apart
  across the cell edge. The symmetry copy of HOH 1 gets no phantoms, the copy
  of HOH 2 does.
"""
from __future__ import absolute_import, division, print_function
import re
import iotbx.pdb
import mmtbx.model
from libtbx import group_args
from libtbx.utils import null_out
from mmtbx.hydrogens import reduce_hydrogen, default_probe_phil
from mmtbx.reduce import Optimizers

pdb_str = """\
CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1
HETATM    1  O   HOH A   1      10.000  10.000  10.000  1.00 10.00           O
HETATM    2  H1  HOH A   1      10.960  10.000  10.000  1.00 10.00           H
HETATM    3  H2  HOH A   1       9.760  10.929  10.000  1.00 10.00           H
HETATM    4  O   HOH A   2      10.000   7.200  10.000  1.00 10.00           O
HETATM    5 MG    MG A   3      20.000  10.000  10.000  1.00 10.00          Mg
HETATM    6  O   HOH A   4      22.100  10.000  10.000  1.00 10.00           O
HETATM    7  O   HOH A   5      30.000  10.000  10.000  1.00 10.00           O
HETATM    8  H1 AHOH A   5      30.960  10.000  10.000  0.50 10.00           H
HETATM    9  H2 AHOH A   5      29.760  10.929  10.000  0.50 10.00           H
HETATM   10  O   HOH A   6      30.000   7.200  10.000  1.00 10.00           O
HETATM   11  O   HOH A   7      40.000  10.000  10.000  1.00 10.00           O
HETATM   12  D1  HOH A   7      40.960  10.000  10.000  1.00 10.00           D
HETATM   13  D2  HOH A   7      39.760  10.929  10.000  1.00 10.00           D
HETATM   14  O   HOH A   8      40.000   7.200  10.000  1.00 10.00           O
HETATM   15 CL  ACL  A   9      50.000  50.000  50.000  0.50 10.00          Cl
HETATM   16 CL  BCL  A   9      50.500  50.000  50.000  0.50 10.00          Cl
HETATM   17  O   HOH A  10      10.000  30.000  10.000  1.00 10.00           O
HETATM   18  H1  HOH A  10      10.960  30.000  10.000  1.00 10.00           H
HETATM   19  O   HOH A  11      10.000  27.200  10.000  1.00 10.00           O
HETATM   20  O   HOH A  12      20.000  30.000  10.000  1.00 10.00           O
HETATM   21  D1  HOH A  12      20.960  30.000  10.000  1.00 10.00           D
HETATM   22  O   HOH A  13      20.000  27.200  10.000  1.00 10.00           O
END
"""

pdb_edge = """\
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
HETATM    1  O   HOH A   1       0.500  10.000  10.000  1.00 10.00           O
HETATM    2  H1  HOH A   1       1.460  10.000  10.000  1.00 10.00           H
HETATM    3  H2  HOH A   1       0.260  10.929  10.000  1.00 10.00           H
HETATM    4  O   HOH A   2      17.700  10.000  10.000  1.00 10.00           O
END
"""

def process(pdb_text, link_all=False):
  m = mmtbx.model.manager(
    model_input=iotbx.pdb.input(source_info=None, lines=pdb_text), log=null_out())
  p = reduce_hydrogen.get_reduce_pdb_interpretation_params(False)
  p.pdb_interpretation.automatic_linking.link_all = link_all
  m.process(make_restraints=True, pdb_interpretation_params=p)
  return m

def symmetry_copies(model):
  """Copies of the residue groups within 3.6 A, as reduce2 builds them."""
  cs = reduce_hydrogen.contact_symmetry(model.crystal_symmetry())
  h = model.get_hierarchy()
  originals = h.atoms()
  copies, hierarchies = [], []
  for entry in h.sym_equiv_residue_groups(cs, 3.6):
    hierarchies.append(entry.hierarchy)
    for k, c in zip(entry.i_seqs, entry.hierarchy.atoms()):
      copies.append((c, originals[k], entry.rt_mx))
  return group_args(crystal_symmetry=cs, models={0: copies},
    hierarchies=hierarchies)

def phantoms(model, copies=None):
  """{(pass altloc, water): number of phantoms}, from the verbosity 3 log;
  water is 'chain A HOH 1', or 'copy of chain A HOH 1' for symmetry copies"""
  opt = Optimizers.Optimizer(default_probe_phil(), False, model,
    verbosity=3, fillAtomDump=False, symmetryCopies=copies)
  out, alt = {}, None
  for l in opt.getInfo().split('\n'):
    m = re.search(r"Running Reduce optimization on model index \d+, "
      r"alternate '(.*)'", l)
    if m:
      alt = m.group(1)
      continue
    m = re.match(r"\s*Added (\d+) phantom Hydrogens on (the )?(copy of )?"
      r"(chain \S+ \S*HOH \d+)", l)
    if m:
      key = (alt, (m.group(3) or '') + m.group(4))
      out[key] = out.get(key, 0) + int(m.group(1))
  return out

def test_water_with_h_or_d():
  '''HOH 1 (H) and HOH 7 (D) get no phantoms in any pass'''
  p = phantoms(process(pdb_str))
  for alt in ('A', 'B'):
    for w in ('chain A HOH 1', 'chain A HOH 7'):
      assert p.get((alt, w), 0) == 0, (alt, w, p)

def test_water_with_one_h_or_d():
  '''HOH 10 (one H) and HOH 12 (one D) get no phantoms: the model's H as given;
  HOH 11 and HOH 13, bare and 2.8 A from them, do'''
  p = phantoms(process(pdb_str))
  for alt in ('A', 'B'):
    for w in ('chain A HOH 10', 'chain A HOH 12'):
      assert p.get((alt, w), 0) == 0, (alt, w, p)
    for w in ('chain A HOH 11', 'chain A HOH 13'):
      assert p.get((alt, w), 0) > 0, (alt, w, p)

def test_bare_water():
  '''HOH 2 and HOH 8, bare, still get phantoms'''
  p = phantoms(process(pdb_str))
  for alt in ('A', 'B'):
    for w in ('chain A HOH 2', 'chain A HOH 8'):
      assert p.get((alt, w), 0) > 0, (alt, w, p)

def test_bonded_non_h():
  '''HOH 4 bonded to MG 3 and carrying no H still gets phantoms'''
  m = process(pdb_str, link_all=True)
  atoms = m.get_atoms()
  bp = m.get_restraints_manager().geometry.get_all_bond_proxies(
    sites_cart=m.get_sites_cart())[0]
  bonded = set()
  for b in bp:
    bonded.add(tuple(sorted(atoms[i].id_str() for i in b.i_seqs)))
  assert ('pdb=" O   HOH A   4 "', 'pdb="MG    MG A   3 "') in bonded, bonded
  p = phantoms(m)
  for alt in ('A', 'B'):
    assert p.get((alt, 'chain A HOH 4'), 0) > 0, (alt, p)

def test_altloc_pass():
  '''HOH 5: H only in altloc A; no phantoms in pass A, phantoms in pass B'''
  p = phantoms(process(pdb_str))
  assert p.get(('A', 'chain A HOH 5'), 0) == 0, p
  assert p.get(('B', 'chain A HOH 5'), 0) > 0, p

def test_symmetry_copies():
  '''the copy of HOH 1 (H1, H2) gets no phantoms, the copy of HOH 2 does;
  counted on the Optimizer's environment atoms (copies and their phantoms)'''
  m = process(pdb_edge)
  opt = Optimizers.Optimizer(default_probe_phil(), False, m, verbosity=0,
    fillAtomDump=False, symmetryCopies=symmetry_copies(m))
  p = {}
  for a in opt._environmentAtoms:
    if a.name == ' H?':
      resseq = a.parent().parent().resseq_as_int()
      p[resseq] = p.get(resseq, 0) + 1
  assert p.get(1, 0) == 0, p
  assert p.get(2, 0) > 0, p

if __name__ == "__main__":
  test_water_with_h_or_d()
  test_water_with_one_h_or_d()
  test_bare_water()
  test_bonded_non_h()
  test_altloc_pass()
  test_symmetry_copies()
  print("OK")
