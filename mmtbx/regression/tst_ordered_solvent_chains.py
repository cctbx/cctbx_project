from __future__ import absolute_import, division, print_function
import time
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
from scitbx import matrix
import mmtbx.refinement.refinement_flags
import mmtbx.solvent.ordered_solvent as ordered_solvent
from mmtbx.regression import tst_ordered_solvent_altlocs as t

"""
mmtbx.solvent.ordered_solvent.assign_water_chains: water goes next to the
nearest chain (symmetry copies included) and takes its ID.
"""

CS = "CRYST1   25.000   25.000   25.000  90.00  90.00  90.00 P 1 21 1"

def gly(chain, shift, resseq=1, i0=0):
  s = matrix.col(shift)
  return [
    t.atom_line(i0+1, " N  ", matrix.col((5.0, 5.0, 5.0))+s, resname="GLY",
      chain=chain, resseq=resseq),
    t.atom_line(i0+2, " CA ", matrix.col((6.4, 5.0, 5.0))+s, resname="GLY",
      chain=chain, resseq=resseq),
    t.atom_line(i0+3, " C  ", matrix.col((7.0, 6.4, 5.0))+s, resname="GLY",
      chain=chain, resseq=resseq),
    t.atom_line(i0+4, " O  ", matrix.col((6.3, 7.4, 5.0))+s, resname="GLY",
      chain=chain, resseq=resseq)]

def hoh(i, xyz, chain="S", resseq=None, altloc=" ", occ=1.0, b=20.,
        name=" O  ", element="O"):
  if(resseq is None): resseq = i
  return t.atom_line(100+i, name, xyz, altloc=altloc, occ=occ, resname="HOH",
    chain=chain, resseq=resseq, element=element, b=b, het=True)

def symmetry_copy(xyz):
  """ 2(1) along b of CS """
  return matrix.col(((-xyz[0])%25, (xyz[1]+12.5)%25, (-xyz[2])%25))

A_N, A_O = matrix.col((5., 5., 5.)), matrix.col((6.3, 7.4, 5.))
B_N, B_O = matrix.col((13., 13., 13.)), matrix.col((14.3, 15.4, 13.))
W1 = A_O + matrix.col((0, 2.8, 0))
W2 = B_O + matrix.col((0, 2.8, 0))
W3 = B_N + matrix.col((-2.8, 0, 0))   # where it has to be after the move
W4 = A_N + matrix.col((0, 0, -2.8))
W5 = A_N + matrix.col((-2.8, 0, 0))
W6A = matrix.col((17.2, 13., 13.))
W6B = matrix.col((17.9, 13.5, 13.))

def get_lines():
  aniso = "ANISOU%5d  O   HOH S%4d     1500   2500   3500    100    200" \
          "    300       O"
  return gly("A", (0,0,0)) + gly("B", (8,8,8), i0=10) + [
    hoh(1, W1, occ=0.7, b=22.),
    hoh(2, W2), aniso%(102, 2),
    hoh(3, symmetry_copy(W3)), aniso%(103, 3),
    hoh(4, W4, chain="A", resseq=50),
    hoh(5, W5),
    hoh(5, W5+matrix.col((-0.59, 0.76, 0)), name=" H1 ", element="H"),
    hoh(5, W5+matrix.col((-0.59,-0.76, 0)), name=" H2 ", element="H"),
    hoh(6, W6A, altloc="A", occ=0.6),
    hoh(6, W6B, altloc="B", occ=0.4)]

def water_info(model):
  """ (chain id, resseq, altloc, name) -> atom, for water """
  result = {}
  for a in model.get_hierarchy().atoms():
    ag = a.parent()
    if(ag.resname != "HOH"): continue
    rg = ag.parent()
    result[(rg.parent().id, rg.resseq_as_int(), ag.altloc, a.name.strip())] = a
  return result

def f_calc(model):
  return abs(model.get_xray_structure().structure_factors(
    d_min=2.0, algorithm="direct").f_calc()).data()

def exercise_main():
  model = t.get_model(get_lines(), cs=CS)
  model.setup_scattering_dictionaries(scattering_table="n_gaussian")
  model.process(make_restraints=True)
  n = model.size()
  assert n == 8+9
  atoms = model.get_hierarchy().atoms()
  def i_of(resseq, name="O", altloc=""):
    for a in atoms:
      if(a.parent().resname == "HOH" and a.name.strip() == name and
         a.parent().parent().resseq_as_int() == resseq and
         a.parent().altloc == altloc): return a.i_seq
  # refinement flags, to see that they follow the atoms
  aniso = model.get_xray_structure().use_u_aniso()
  assert aniso.count(True) == 2
  occ_groups = [[flex.size_t([i_of(1)])],
                [flex.size_t([i_of(6, altloc="A")]),
                 flex.size_t([i_of(6, altloc="B")])]]
  model.set_refinement_flags(mmtbx.refinement.refinement_flags.manager(
    individual_sites     = True,
    individual_adp       = True,
    occupancies          = True,
    sites_individual     = flex.bool(n, True),
    adp_individual_iso   = ~aniso,
    adp_individual_aniso = aniso,
    s_occupancies        = occ_groups))
  g = model.get_restraints_manager().geometry
  n_bonds = g.pair_proxies(sites_cart=model.get_sites_cart()
    ).bond_proxies.simple.size()
  n_angles = g.angle_proxies.size()
  ws = ordered_solvent.water_state(model=model, keep_input_water=True)
  assert ws.fixed_selection(model).count(True) == 9
  f_before = f_calc(model)
  #
  model = ordered_solvent.assign_water_chains(model=model, water_state=ws)
  #
  assert model.size() == n
  assert approx_equal(model.get_xray_structure().sites_cart(),
    model.get_hierarchy().atoms().extract_xyz())
  assert approx_equal(f_calc(model), f_before, 1.e-6)
  w = water_info(model)
  # chain IDs and residue numbers: staying water keeps its number, others
  # get the next ones
  assert sorted(w.keys()) == sorted([
    ("A", 50, "", "O"), ("A", 51, "", "O"),
    ("A", 52, "", "O"), ("A", 52, "", "H1"), ("A", 52, "", "H2"),
    ("B", 2, "", "O"), ("B", 3, "", "O"),
    ("B", 4, "A", "O"), ("B", 4, "B", "O")]), sorted(w.keys())
  # symmetry copy next to its chain
  a = w[("B", 3, "", "O")]
  assert approx_equal(a.xyz, W3, 1.e-3)
  assert approx_equal(w[("B", 2, "", "O")].xyz, W2, 1.e-3)
  # occupancy and ADP stay
  a = w[("A", 51, "", "O")]
  assert approx_equal(a.occ, 0.7) and approx_equal(a.b, 22.)
  sc = model.get_xray_structure().scatterers()
  assert approx_equal(sc[a.i_seq].occupancy, 0.7)
  for key in [("B", 2, "", "O"), ("B", 3, "", "O")]:
    assert sc[w[key].i_seq].flags.use_u_aniso()
  assert model.get_xray_structure().use_u_aniso().count(True) == 2
  assert approx_equal(w[("B", 2, "", "O")].uij,
    (0.15, 0.25, 0.35, 0.01, 0.02, 0.03), 1.e-4)
  # H are with their O
  for name in ["H1", "H2"]:
    d = w[("A", 52, "", name)].distance(w[("A", 52, "", "O")])
    assert approx_equal(d, 0.962, 0.01), d
  # one water chain per chain ID, after everything else
  chains = [(c.id, c.residue_groups()[0].atom_groups()[0].resname)
            for c in model.get_hierarchy().chains()]
  assert chains == [("A","GLY"), ("B","GLY"), ("A","HOH"), ("B","HOH")], \
    chains
  # restraints
  g = model.get_restraints_manager().geometry
  assert g.pair_proxies(sites_cart=model.get_sites_cart()
    ).bond_proxies.simple.size() == n_bonds
  assert g.angle_proxies.size() == n_angles
  # refinement flags
  rf = model.refinement_flags
  aniso = model.get_xray_structure().use_u_aniso()
  assert rf.adp_individual_aniso.all_eq(aniso)
  assert rf.adp_individual_iso.all_eq(~aniso)
  assert rf.sites_individual.count(True) == n
  got = [[list(s) for s in group] for group in rf.s_occupancies]
  assert got == [[[w[("A", 51, "", "O")].i_seq]],
                 [[w[("B", 4, "A", "O")].i_seq],
                  [w[("B", 4, "B", "O")].i_seq]]], got
  # water to keep is still known
  assert ws.fixed_selection(model).count(True) == 9
  assert ("A", "  51 ") in ws.input_water_keys, ws.input_water_keys
  # nothing to do the second time: the same object comes back
  model2 = ordered_solvent.assign_water_chains(model=model, water_state=ws)
  assert model2 is model
  t.no_relabel_needed(model)

def exercise_no_polymer():
  """ Without polymer water goes to the nearest other chain, or stays """
  lines = [
    t.atom_line(1, " C1 ", (5., 5., 5.), resname="LIG", chain="L",
      element="C", het=True),
    t.atom_line(2, " O1 ", (6.4, 5., 5.), resname="LIG", chain="L",
      element="O", het=True),
    hoh(1, (9.2, 5., 5.))]
  model = t.get_model(lines, cs=CS)
  model = ordered_solvent.assign_water_chains(model=model)
  assert list(water_info(model).keys()) == [("L", 2, "", "O")]
  # water only
  model = t.get_model([hoh(1, (9.2, 5., 5.)), hoh(2, (12., 5., 5.))], cs=CS)
  model2 = ordered_solvent.assign_water_chains(model=model)
  assert model2 is model
  assert sorted(water_info(model).keys()) == [
    ("S", 1, "", "O"), ("S", 2, "", "O")]
  # water farther than the largest cutoff from anything: stays
  cs = "CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1"
  model = t.get_model(gly("A", (0,0,0)) + [hoh(1, (50., 50., 50.))], cs=cs)
  model = ordered_solvent.assign_water_chains(model=model)
  assert list(water_info(model).keys()) == [("S", 1, "", "O")]
  # 10 A from the chain: found
  model = t.get_model(gly("A", (0,0,0)) + [hoh(1, (5., 5., 15.))], cs=cs)
  model = ordered_solvent.assign_water_chains(model=model)
  assert list(water_info(model).keys()) == [("A", 2, "", "O")]

def exercise_numbering():
  """ Residue numbers beyond 9999; chain ID of two characters """
  cs = "CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1"
  lines = gly("A", (0,0,0), resseq=9998) + [
    hoh(1, W1), hoh(2, W4), hoh(3, W5)]
  model = t.get_model(lines, cs=cs)
  model = ordered_solvent.assign_water_chains(model=model)
  resseqs = [rg.resseq for rg in model.get_hierarchy().residue_groups()
             if rg.atom_groups()[0].resname == "HOH"]
  assert resseqs == ["9999", "A000", "A001"], resseqs
  assert sorted(water_info(model).keys()) == [
    ("A", 9999, "", "O"), ("A", 10000, "", "O"), ("A", 10001, "", "O")]
  lines = [l[:20]+"AA"+l[22:] for l in gly("A", (0,0,0))] + [hoh(1, W1)]
  model = t.get_model(lines, cs=cs)
  model = ordered_solvent.assign_water_chains(model=model)
  assert list(water_info(model).keys()) == [("AA", 2, "", "O")]

def exercise_bonded_water():
  """ Water bonded to something (link to metal) is not touched """
  lines = gly("A", (0,0,0)) + [hoh(1, W1), hoh(2, symmetry_copy(W4))]
  model = t.get_model(lines, cs=CS)
  model.process(make_restraints=True)
  g = model.get_restraints_manager().geometry
  from cctbx import geometry_restraints
  i_w1 = 4
  g.add_new_bond_restraints_in_place(
    proxies = [geometry_restraints.bond_simple_proxy(
      i_seqs=(3, i_w1), distance_ideal=2.8, weight=1.)],
    sites_cart = model.get_sites_cart())
  model = ordered_solvent.assign_water_chains(model=model)
  assert sorted(water_info(model).keys()) == [
    ("A", 2, "", "O"), ("S", 1, "", "O")], sorted(water_info(model).keys())

def exercise_one_chain_id():
  """
  Model with one chain ID, water of the model is in the middle of it, with
  TLS group and occupancy group that include water: water becomes a chain
  object of its own after the chain, flags follow the atoms, and the second
  call has nothing to do.
  """
  cs = "CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1"
  lines = gly("A", (0,0,0)) + [
    hoh(1, W1, chain="A", resseq=10),
    hoh(2, W4, chain="A", resseq=11, occ=0.6),
    t.atom_line(20, " N  ", (6.4, 5.0, 2.0), resname="GLY", chain="A",
      resseq=12)]
  model = t.get_model(lines, cs=cs)
  n = model.size()
  model.set_refinement_flags(mmtbx.refinement.refinement_flags.manager(
    individual_sites = True,
    individual_adp   = True,
    tls              = True,
    occupancies      = True,
    sites_individual = flex.bool(n, True),
    adp_individual_iso = flex.bool([True]*5+[False]+[True]),
    adp_tls          = [flex.size_t([0,1,2,3,4,5]), flex.size_t([6])],
    s_occupancies    = [[flex.size_t([5])]]))
  model = ordered_solvent.assign_water_chains(model=model)
  atoms = model.get_hierarchy().atoms()
  names = [(a.parent().resname, a.parent().parent().resseq_as_int())
           for a in atoms]
  assert names == [("GLY",1)]*4 + [("GLY",12), ("HOH",10), ("HOH",11)], names
  chains = [(c.id, c.atoms_size()) for c in model.get_hierarchy().chains()]
  assert chains == [("A", 5), ("A", 2)], chains
  rf = model.refinement_flags
  assert list(rf.adp_individual_iso) == [True]*6+[False]
  assert [list(x) for x in rf.adp_tls] == [[0,1,2,3,5,6], [4]], \
    [list(x) for x in rf.adp_tls]
  assert [[list(y) for y in x] for x in rf.s_occupancies] == [[[6]]]
  assert approx_equal(atoms[6].occ, 0.6)
  model2 = ordered_solvent.assign_water_chains(model=model)
  assert model2 is model

def exercise_margin():
  """
  Water about as far from two chains stays in the chain it is in; it moves
  to the other one only if that one is clearly nearer.
  """
  cs = "CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1"
  shift = (1.3, 8.2, 0) # N of chain B is 5.8 A from O of chain A
  assert approx_equal((A_N+matrix.col(shift)-A_O).length(), 5.8, 1.e-3)
  for d_a, expected in [(2.85, "B"), (2.6, "A")]:
    w = A_O + matrix.col((0, d_a, 0))
    lines = gly("A", (0,0,0)) + gly("B", shift, i0=10) + [
      hoh(1, w, chain="B", resseq=7)]
    model = t.get_model(lines, cs=cs)
    model = ordered_solvent.assign_water_chains(model=model)
    keys = list(water_info(model).keys())
    assert len(keys) == 1 and keys[0][0] == expected, (d_a, keys)
    if(expected == "B"): assert keys[0][1] == 7 # and keeps its number

def run():
  exercise_one_chain_id()
  exercise_margin()
  exercise_main()
  exercise_no_polymer()
  exercise_numbering()
  exercise_bonded_water()

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %6.2f"%(time.time()-t0))
