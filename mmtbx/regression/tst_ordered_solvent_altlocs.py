from __future__ import absolute_import, division, print_function
import time
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out
from scitbx import matrix
import iotbx.pdb
import mmtbx.model
import mmtbx.f_model
import mmtbx.solvent.ordered_solvent as ordered_solvent
from iotbx.pdb import water_altlocs

"""
Water in alternative conformations in mmtbx.solvent.ordered_solvent.
"""

cryst1 = "CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1"

N   = matrix.col((10.0, 10.0, 10.0))
CA  = matrix.col((11.4, 10.0, 10.0))
C   = matrix.col((12.0, 11.4, 10.0))
O   = matrix.col((11.3, 12.4, 10.0))
CB  = matrix.col((12.0,  9.1, 11.1))
OGA = matrix.col((13.4,  9.1, 11.1))
OGB = matrix.col((11.5,  7.8, 11.2))

def atom_line(i, name, xyz, altloc=" ", occ=1.0, resname="SER", chain="A",
              resseq=1, element=None, b=10.0, het=False):
  if(element is None): element = name.strip()[0]
  rec = "HETATM" if het else "ATOM  "
  return "%s%5d %-4s%s%3s %s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s"%(
    rec, i, name, altloc, resname, chain, resseq, xyz[0], xyz[1], xyz[2],
    occ, b, element)

def ser_lines(altlocs=True, occ_a=0.7):
  result = [
    atom_line(1, " N  ", N),
    atom_line(2, " CA ", CA),
    atom_line(3, " C  ", C),
    atom_line(4, " O  ", O),
    atom_line(5, " CB ", CB)]
  if(altlocs):
    result.append(atom_line(6, " OG ", OGA, altloc="A", occ=occ_a))
    result.append(atom_line(7, " OG ", OGB, altloc="B", occ=1-occ_a))
  else:
    result.append(atom_line(6, " OG ", OGA))
  return result

def hoh(i, xyz, altloc=" ", occ=1.0, chain="A", b=10.0):
  return atom_line(100+i, " O  ", xyz, altloc=altloc, occ=occ, resname="HOH",
    chain=chain, resseq=100+i, b=b, het=True)

def get_model(lines, cs=cryst1):
  return mmtbx.model.manager(model_input = iotbx.pdb.input(
    source_info=None, lines=[cs]+list(lines)), log=null_out())

def assign(model, sites, heights=None, d_min=1.0, fixed=None):
  if(heights is None): heights = [10.-i for i in range(len(sites))]
  return ordered_solvent.assign_conformers(
    model                 = model,
    sites_cart            = flex.vec3_double([tuple(s) for s in sites]),
    heights               = flex.double(heights),
    params                = ordered_solvent.master_params().extract(),
    d_min                 = d_min,
    fixed_water_selection = fixed)

def no_relabel_needed(model):
  """
  Water in the model is such that iotbx.pdb.water_altlocs has nothing to
  label.
  """
  h = model.get_hierarchy().deep_copy()
  h.atoms().reset_i_seq()
  r = water_altlocs.run(h, model.crystal_symmetry())
  assert r.n_labelled == 0, r.labelled

# site next to OG of conformer B, 0.9 A from it, away from CB
P_B = OGB + (OGB-CB).normalize()*0.9
# the same but with OG of conformer A farther than dist_max
P_B2 = OGB + (OGB-OGA).normalize()*1.0
FREE = matrix.col((18., 18., 2.))

def exercise_no_letters():
  """ No altlocs in the model: no water in alternative conformations """
  m = get_model(ser_lines(altlocs=False))
  assert ordered_solvent.model_altlocs(m) == []
  r = assign(m, [FREE, FREE+matrix.col((1.2,0,0))])
  assert list(r.keep) == [True, False]
  assert r.altlocs == ["", ""]
  assert approx_equal(r.occupancies[0], 1)
  assert r.relabel == {}
  # far enough from each other: both, no altlocs
  r = assign(m, [FREE, FREE+matrix.col((2.5,0,0))])
  assert list(r.keep) == [True, True] and r.altlocs == ["", ""]
  # letters of water do not count
  m = get_model(ser_lines(altlocs=False)+[hoh(1, FREE, altloc="A", occ=0.5)])
  assert ordered_solvent.model_altlocs(m) == []

def exercise_driven():
  """ Site in the place of conformer B gets altloc and occupancy of A """
  m = get_model(ser_lines())
  assert ordered_solvent.model_altlocs(m) == ["A", "B"]
  assert approx_equal((P_B-OGB).length(), 0.9)
  assert 2.2 < (P_B-OGA).length() < 3.2
  r = assign(m, [P_B])
  assert list(r.keep) == [True]
  assert r.altlocs == ["A"], r.altlocs
  assert approx_equal(r.occupancies[0], 0.7)
  # and of B if it is in the place of A
  p_a = OGA + (OGA-CB).normalize()*0.9
  r = assign(m, [p_a])
  assert list(r.keep) == [True] and r.altlocs == ["B"], r.altlocs
  assert approx_equal(r.occupancies[0], 0.3)

def exercise_blank_clash():
  """ Site too close to atom without altloc is rejected """
  m = get_model(ser_lines())
  r = assign(m, [CA + matrix.col((0,0,-1.2))])
  assert list(r.keep) == [False]
  # too close to atom with altloc
  r = assign(m, [OGB + (OGB-CB).normalize()*0.3])
  assert list(r.keep) == [False]
  # excluded by A and by B: no altloc available
  r = assign(m, [(OGA+OGB)*0.5 + matrix.col((0,0,0.9))])
  assert list(r.keep) == [False]

def exercise_free_pair():
  """ Two sites away from atoms with altlocs """
  m = get_model(ser_lines())
  sites = [FREE, FREE+matrix.col((1.2,0,0))]
  r = assign(m, sites, d_min=1.0)
  assert list(r.keep) == [True, True]
  assert r.altlocs == ["A", "B"], r.altlocs
  assert approx_equal(r.occupancies, [0.5, 0.5])
  # stronger one is processed first whatever the order
  r = assign(m, sites, heights=[5, 10], d_min=1.0)
  assert r.altlocs == ["A", "B"] or r.altlocs == ["B", "A"]
  assert list(r.keep) == [True, True]
  # resolution does not support two sites 1.2 A apart: 0.7*2.0 = 1.4 A
  r = assign(m, sites, d_min=2.0)
  assert list(r.keep) == [True, False]
  assert r.altlocs == ["", ""]
  # anisotropic water: two sites closer than 1 A are one water
  sites_close = [FREE, FREE+matrix.col((0.9,0,0))]
  r = assign(m, sites_close, d_min=1.0)
  assert list(r.keep) == [True, True] and r.altlocs == ["A", "B"]
  r = ordered_solvent.assign_conformers(model=m,
    sites_cart=flex.vec3_double([tuple(s) for s in sites_close]),
    heights=flex.double([10, 9]),
    params=ordered_solvent.master_params().extract(), d_min=1.0,
    anisotropic=True)
  assert list(r.keep) == [True, False] and r.altlocs == ["", ""]
  # three sites, two altlocs: the weakest one is rejected
  r = assign(m, sites+[FREE+matrix.col((0.6,1.0,0))], d_min=1.0)
  assert list(r.keep) == [True, True, False], list(r.keep)

def exercise_order():
  """
  A site that has to be B (it excludes an atom of A) is next to a stronger
  site that excludes nothing else: the stronger one becomes A (used to be
  given B first, and the other site was rejected). Same if the stronger one
  is water of the model.
  """
  m = get_model(ser_lines())
  p_a = OGA + (OGA-CB).normalize()*0.9   # in place of OG A: has to be B
  q = p_a + (OGA-CB).normalize()*1.5     # 2.4 A from OG A: excludes nothing
  assert (q-OGA).length() > 2.2 and (q-OGB).length() > 3.2
  r = assign(m, [p_a, q], heights=[5, 10])
  assert list(r.keep) == [True, True], list(r.keep)
  assert r.altlocs == ["B", "A"], r.altlocs
  m = get_model(ser_lines() + [hoh(1, q)])
  r = assign(m, [p_a])
  assert list(r.keep) == [True] and r.altlocs == ["B"], r.altlocs
  assert r.relabel == {m.size()-1: "A"}, r.relabel

def exercise_three_letters():
  """ Three conformers: choice by H-bond partner, then by occupancy """
  assert (P_B2-OGA).length() > 3.2
  partner = P_B2 + matrix.col((0,-2.8,0))
  far     = matrix.col((3., 17., 17.))
  for xyz, expected, occ in [(partner, "C", 0.25), (far, "A", 0.7)]:
    lines = ser_lines() + [
      atom_line(8, " O  ", xyz, altloc="C", occ=0.25, resname="GLY",
        resseq=5)]
    m = get_model(lines)
    assert ordered_solvent.model_altlocs(m) == ["A", "B", "C"]
    r = assign(m, [P_B2])
    assert list(r.keep) == [True]
    assert r.altlocs == [expected], r.altlocs
    assert approx_equal(r.occupancies[0], occ)

def exercise_existing_water():
  """ Site next to existing water without altloc """
  lines = ser_lines() + [hoh(1, FREE)]
  m = get_model(lines)
  i_w = m.size()-1
  site = FREE + matrix.col((1.5,0,0))
  r = assign(m, [site])
  assert list(r.keep) == [True]
  assert r.relabel == {i_w: "A"}, r.relabel
  assert r.altlocs == ["B"]
  # water that must not change: site is rejected
  fixed = flex.bool(m.size(), False)
  fixed[i_w] = True
  r = assign(m, [site], fixed=fixed)
  assert list(r.keep) == [False] and r.relabel == {}
  # existing water has altloc: site gets another one
  m = get_model(ser_lines() + [hoh(1, FREE, altloc="B", occ=0.5)])
  r = assign(m, [site])
  assert list(r.keep) == [True] and r.altlocs == ["A"] and r.relabel == {}

def exercise_metal():
  """ Site 2.1 A from partially occupied metal: coordination """
  lines = ser_lines() + [atom_line(8, "MG  ", FREE, occ=0.6, resname=" MG",
    resseq=5, element="MG", het=True)]
  m = get_model(lines)
  r = assign(m, [FREE + matrix.col((2.1,0,0))])
  assert list(r.keep) == [True] and r.altlocs == [""]
  # 2.1 A from partially occupied atom that is not metal: rejected
  lines = ser_lines() + [atom_line(8, " O  ", FREE, occ=0.6, resname="GLY",
    resseq=5)]
  r = assign(get_model(lines), [FREE + matrix.col((2.1,0,0))])
  assert list(r.keep) == [False]

def exercise_exclusion_distance():
  p = ordered_solvent.master_params().extract()
  lines = ser_lines() + [
    hoh(1, FREE),
    atom_line(8, "MG  ", FREE+matrix.col((5,0,0)), occ=0.6, altloc="A",
      resname=" MG", resseq=5, element="MG", het=True),
    atom_line(9, " N  ", FREE+matrix.col((0,5,0)), occ=0.5, resname="GLY",
      resseq=6)]
  m = get_model(lines)
  d = {}
  for a in m.get_hierarchy().atoms():
    key = (a.parent().resname.strip(), a.name.strip(), a.parent().altloc)
    d[key] = ordered_solvent.exclusion_distance(a, p)
  assert water_altlocs.DISTANCE_CUTOFF == 2.2
  assert approx_equal(d[("HOH","O","")], 2.2)
  assert approx_equal(d[("SER","OG","A")], 2.2)
  assert approx_equal(d[("SER","N","")], 1.8)
  assert approx_equal(d[("MG","MG","A")], 1.8)
  assert approx_equal(d[("GLY","N","")], 2.2) # partial occupancy

def get_fmodel(model_true, model, d_min):
  xrs = model_true.get_xray_structure()
  f_obs = abs(xrs.structure_factors(d_min=d_min).f_calc())
  r_free_flags = f_obs.generate_r_free_flags()
  fmodel = mmtbx.f_model.manager(
    f_obs          = f_obs,
    r_free_flags   = r_free_flags,
    xray_structure = model.get_xray_structure())
  fmodel.update_all_scales()
  return fmodel

class manager_no_map_filter(ordered_solvent.manager):
  """
  A model of a few atoms in a big box is far from the data without one
  water, and sigma-scaled maps have no meaningful scale: map filter is not
  what these tests are about.
  """
  def _filter_map(self): pass

def run_manager(lines_true, lines, d_min, params=None, process=True,
                water_state=None, map_filter=True):
  model_true = get_model(lines_true)
  model = get_model(lines)
  for m in [model_true, model]:
    m.setup_scattering_dictionaries(scattering_table="n_gaussian")
  if(process): model.process(make_restraints=True)
  fmodel = get_fmodel(model_true, model, d_min)
  if(params is None):
    params = ordered_solvent.master_params().extract()
    params.add_hydrogens.d_min = None
  manager = ordered_solvent.manager
  if(not map_filter): manager = manager_no_map_filter
  o = manager(fmodel=fmodel, fmodels=None, model=model,
    is_neutron_scat_table=False, params=params, water_state=water_state,
    log=null_out())
  assert approx_equal(o.model.get_xray_structure().sites_cart(),
    o.model.get_hierarchy().atoms().extract_xyz())
  no_relabel_needed(o.model)
  return o

def waters(model):
  result = []
  for ag in model.get_hierarchy().atom_groups():
    if(ag.resname != "HOH"): continue
    for a in ag.atoms():
      if(a.element_is_hydrogen()): continue
      result.append(a)
  return result

def find_water(model, xyz, tolerance=0.3):
  result = [a for a in waters(model)
            if (matrix.col(a.xyz)-matrix.col(xyz)).length() < tolerance]
  assert len(result) == 1, len(result)
  return result[0]

def exercise_split_ser_end_to_end():
  """ Water that is there only with conformer A is built with altloc A """
  w_o = O + matrix.col((0,2.8,0))
  lines = ser_lines() + [hoh(1, w_o)]
  lines_true = lines + [hoh(2, P_B, altloc="A", occ=0.7)]
  for process in [True, False]:
    o = run_manager(lines_true, lines, d_min=1.0, process=process)
    assert len(waters(o.model)) == 2, len(waters(o.model))
    a = find_water(o.model, P_B)
    assert a.parent().altloc == "A"
    assert 0.4 < a.occ <= 1.0, a.occ
    assert find_water(o.model, w_o).parent().altloc == ""
    assert o.fmodel.r_work() < 0.1, o.fmodel.r_work()

def exercise_free_pair_end_to_end():
  """
  Two alternative positions of water away from conformers of the model are
  built as two sites if resolution allows it, and as one otherwise.
  """
  w1 = O + matrix.col((0,2.8,0))
  w2 = w1 + matrix.col((1.3,0,0))
  lines = ser_lines()
  lines_true = lines + [hoh(1, w1, altloc="A", occ=0.6),
                        hoh(2, w2, altloc="B", occ=0.4)]
  o = run_manager(lines_true, lines, d_min=0.9)
  ws = waters(o.model)
  assert len(ws) == 2, len(ws)
  a1, a2 = find_water(o.model, w1), find_water(o.model, w2)
  assert sorted([a1.parent().altloc, a2.parent().altloc]) == ["A", "B"]
  assert a1.occ > a2.occ, (a1.occ, a2.occ)
  o = run_manager(lines_true, lines, d_min=2.0)
  ws = waters(o.model)
  assert len(ws) == 1, len(ws)
  assert ws[0].parent().altloc == ""

def exercise_existing_water_preserved():
  """
  Water of the model keeps its ADP (also anisotropic), occupancy, H when new
  water is added, also when model has alternative conformations.
  """
  w1 = O + matrix.col((0,2.8,0))
  w2 = N + matrix.col((0,0,-2.8))
  w3 = w2 + matrix.col((0,2.8,0))
  w4 = w1 + matrix.col((0,0,2.8))
  aniso = "ANISOU  101  O   HOH A 101     1500   2500   3500    100    200" \
          "    300       O"
  lines = ser_lines() + [hoh(1, w1, b=20), aniso,
    hoh(2, w2, b=15),
    atom_line(103, " H1 ", w2+matrix.col((0.76,0,-0.59)), resname="HOH",
      resseq=102, element="H", het=True, b=15),
    atom_line(104, " H2 ", w2+matrix.col((-0.76,0,-0.59)), resname="HOH",
      resseq=102, element="H", het=True, b=15),
    hoh(3, w3, occ=0.6, b=33)]
  lines_true = lines + [hoh(4, w4, b=20)]
  o = run_manager(lines_true, lines, d_min=1.0, map_filter=False)
  m = o.model
  assert len(waters(m)) >= 4, len(waters(m))
  xrs = m.get_xray_structure()
  sc = xrs.scatterers()
  a = find_water(m, w1)
  assert sc[a.i_seq].flags.use_u_aniso()
  assert approx_equal(a.uij, (0.15, 0.25, 0.35, 0.01, 0.02, 0.03), 1.e-4)
  a = find_water(m, w2)
  assert [x.name.strip() for x in a.parent().atoms()] == ["O", "H1", "H2"]
  assert approx_equal(a.b, 15, 1.e-2)
  a = find_water(m, w3)
  assert approx_equal(a.occ, 0.6, 1.e-3) and approx_equal(a.b, 33, 1.e-2)
  assert approx_equal(sc[a.i_seq].occupancy, 0.6, 1.e-3)
  find_water(m, w4)
  # all water is in chain A now
  assert set([a.parent().parent().parent().id for a in waters(m)]) == \
    set(["A"])

def exercise_orphan_reverts():
  """
  Water that got altloc because of another water that was then removed has
  no altloc at the end.
  """
  w1 = O + matrix.col((0,2.8,0))
  w2 = w1 + matrix.col((1.3,0,0))
  lines = ser_lines()
  lines_true = lines + [hoh(1, w1, altloc="A", occ=0.6),
                        hoh(2, w2, altloc="B", occ=0.4)]
  params = ordered_solvent.master_params().extract()
  params.add_hydrogens.d_min = None
  params.occupancy_min = 0.58 # weaker site does not make it
  o = run_manager(lines_true, lines, d_min=0.9, params=params)
  ws = waters(o.model)
  assert len(ws) == 1, len(ws)
  assert ws[0].parent().altloc == "", ws[0].parent().altloc

def exercise_filter_only():
  """ mode=filter_only: nothing is added """
  w1 = O + matrix.col((0,2.8,0))
  lines = ser_lines()
  lines_true = lines + [hoh(1, w1)]
  params = ordered_solvent.master_params().extract()
  params.mode = "filter_only"
  o = run_manager(lines_true, lines, d_min=1.0, params=params)
  assert len(waters(o.model)) == 0

def run():
  exercise_no_letters()
  exercise_driven()
  exercise_blank_clash()
  exercise_free_pair()
  exercise_order()
  exercise_three_letters()
  exercise_existing_water()
  exercise_metal()
  exercise_exclusion_distance()
  exercise_split_ser_end_to_end()
  exercise_free_pair_end_to_end()
  exercise_existing_water_preserved()
  exercise_orphan_reverts()
  exercise_filter_only()

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %6.2f"%(time.time()-t0))
