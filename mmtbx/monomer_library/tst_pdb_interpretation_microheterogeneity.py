"""EXPECTED FAILURE (tst_list_expected_failures in mmtbx/run_tests.py) until
pdb_interpretation handles an atom shared by the residues of a
microheterogeneity site.

1ejg residue 22 is Pro in conformer A and Ser in B and C, with one N without
altloc, labelled PRO. That N belongs to every conformer, but
pdb_interpretation gives it to the residue it is labelled with only, so the
Ser conformers have no N:
- two residue groups (the site as read from a PDB file): N-CA of Ser is a
  bond of automatic linking (1.430 A, sigma 0.01; the model has 1.39 and
  1.60 A) instead of the Ser dictionary's, and the angles C(21)-N-CA,
  N-CA-C and N-CA-CB are missing;
- one residue group (mmCIF): no N-CA bond at all.
Refinement restrains these residues accordingly. The amide H of the Ser
conformers gets no N-H bond either: mmtbx.reduce2 adds what placement needs
itself (place_hydrogens.add_shared_parent_h_restraints(),
mmtbx/regression/tst_reduce2_microheterogeneity.py), and the model it writes
has the same gap when it is interpreted again.

When this test passes: move it to the regular tests and remove
add_shared_parent_h_restraints() from mmtbx/hydrogens/reduce_hydrogen.py.
"""
from __future__ import absolute_import, division, print_function
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
from mmtbx.regression.tst_reduce2_pass_environment import m_1ejg_22

def restraints(merge):
  inp = iotbx.pdb.input(source_info=None, lines=m_1ejg_22.splitlines())
  h = inp.construct_hierarchy()
  if merge:
    for chain in h.chains():
      primary = None
      for rg in chain.residue_groups():
        if primary is not None and primary.resid() == rg.resid():
          chain.merge_residue_groups(primary=primary, secondary=rg)
        else:
          primary = rg
    h.atoms().reset_serial()
  model = mmtbx.model.manager(model_input=None, pdb_hierarchy=h,
    crystal_symmetry=inp.crystal_symmetry(), log=null_out())
  model.process(make_restraints=True)
  grm = model.get_restraints_manager().geometry
  atoms = model.get_atoms()
  label = lambda i: (atoms[i].parent().parent().resseq.strip(),
    atoms[i].parent().altloc, atoms[i].name.strip())
  bonds = {}
  simple, asu = grm.get_all_bond_proxies(sites_cart=model.get_sites_cart())
  for p in simple:
    bonds[frozenset(label(i) for i in p.i_seqs)] = p.distance_ideal
  angles = {}
  for p in grm.angle_proxies:
    i, j, k = [label(x) for x in p.i_seqs]
    angles[(frozenset([i, k]), j)] = p.angle_ideal
  return bonds, angles

def exercise(merge):
  bonds, angles = restraints(merge)
  n, c_prev = ("22", "", "N"), ("21", "", "C")
  for altloc in ("B", "C"):
    ca, c, cb = [("22", altloc, x) for x in ("CA", "C", "CB")]
    # the Ser dictionary's N-CA (1.458 A)
    d = bonds.get(frozenset([n, ca]))
    assert d is not None, (merge, altloc, "no N-CA bond")
    assert abs(d - 1.458) < 0.01, (merge, altloc, d)
    for ends, centre in (((c_prev, ca), n), ((n, c), ca), ((n, cb), ca)):
      assert (frozenset(ends), centre) in angles, (merge, altloc, ends, centre)

def run():
  exercise(merge=False)
  exercise(merge=True)
  print("OK")

if (__name__ == "__main__"):
  run()
