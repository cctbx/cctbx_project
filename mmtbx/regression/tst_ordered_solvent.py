from __future__ import absolute_import, division, print_function
import time
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
import iotbx.pdb
import mmtbx.model
import mmtbx.solvent.ordered_solvent as ordered_solvent
import mmtbx.f_model
import math
from scitbx import matrix
from libtbx.utils import null_out

pdb_str = """\
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
ATOM      1  N   GLY A   1       5.000   5.000   5.000  1.00 10.00           N
ATOM      2  CA  GLY A   1       6.000   5.000   5.000  1.00 10.00           C
ATOM      3  C   GLY A   1       7.000   5.000   5.000  1.00 10.00           C
ATOM      4  O   GLY A   1       7.500   6.000   5.000  1.00 10.00           O
TER
END
"""

def exercise_added_solvent_atom_name():
  """
  Water added by the ordered-solvent machinery must have the oxygen
  right-justified in the PDB atom-name field, i.e. the element in column 14
  (" O  "), as the PDB spec requires for single-character elements. This calls
  add_solvent_to_model_inplace directly -- the exact call the ordered_solvent
  manager (phenix.refine / ligand pipeline) makes in _add_new_solvent, and
  which model.add_solvent also routes through. Regression test for the oxygen
  ending up in column 13 ("O   "), which coot and other parsers reject.
  """
  model = mmtbx.model.manager(
    model_input = iotbx.pdb.input(source_info=None, lines=pdb_str.splitlines()))
  model.setup_scattering_dictionaries(scattering_table="wk1995")
  sites_frac = flex.vec3_double([(0.40, 0.40, 0.40), (0.60, 0.60, 0.60)])
  params = ordered_solvent.master_params().extract()
  ordered_solvent.add_solvent_to_model_inplace(
    sites=sites_frac, model=model, params=params)
  # print(model.model_as_pdb())
  water_lines = [l for l in model.model_as_pdb().splitlines()
                 if l.startswith(("ATOM", "HETATM")) and l[17:20] == "HOH"]
  assert len(water_lines) == 2, water_lines
  for line in water_lines:
    # PDB columns 13-16 (0-based 12:16) are the atom-name field; a
    # single-character element must sit in column 14 with a leading blank.
    assert line[12:16] == " O  ", \
      "water oxygen must be in column 14, got name field %r in:\n%s" % (
        line[12:16], line)

pdb_str_with_water = pdb_str.replace("END\n", """\
HETATM    5  O   HOH S   1      10.000  10.000  10.000  1.00 30.00           O
HETATM    6  O   HOH S   2      12.000  10.000  10.000  1.00 30.00           O
HETATM    7  O   HOH S   3      14.000  10.000  10.000  1.00 30.00           O
HETATM    8  O   HOH S   4      16.000  10.000  10.000  1.00 30.00           O
END
""")

def exercise_added_solvent_unique_serials():
  """
  Atom serial numbers must stay unique after the filter-then-add sequence the
  ordered_solvent manager runs every macro-cycle. Removing waters with
  model.select() leaves the survivors with their old serials, so new waters
  must not be numbered from model.size(). Regression test for duplicate
  PDB serials / _atom_site.id on waters in phenix.refine output (UR-16031).
  """
  model = mmtbx.model.manager(
    model_input = iotbx.pdb.input(
      source_info=None, lines=pdb_str_with_water.splitlines()))
  model.setup_scattering_dictionaries(scattering_table="wk1995")
  # Emulate the "Filter" steps: drop two of the four waters.
  model = model.select(model.selection("not (chain S and resseq 2:3)"))
  # Emulate "Add new water".
  sites_frac = flex.vec3_double([(0.40, 0.40, 0.40), (0.60, 0.60, 0.60)])
  params = ordered_solvent.master_params().extract()
  ordered_solvent.add_solvent_to_model_inplace(
    sites=sites_frac, model=model, params=params)
  serials = [a.serial.strip() for a in model.get_hierarchy().atoms()]
  assert len(set(serials)) == len(serials), \
    "duplicate atom serials after add_solvent_to_model_inplace: %s" % serials
  pdb_serials = [l[6:11] for l in model.model_as_pdb().splitlines()
                 if l.startswith(("ATOM", "HETATM"))]
  assert len(set(pdb_serials)) == len(pdb_serials), \
    "duplicate serials in PDB output: %s" % pdb_serials

def exercise_added_solvent_merges_into_existing_chain():
  """
  New waters added to a model that already has a water chain with the same
  id must join that chain object, not be appended as a second chain object
  with the same id. A second chain object shows up as an extra label_asym_id
  for the same auth_asym_id in mmCIF output (one per solvent round), which
  confuses downstream tools (UR-16031, PDB-REDO).
  """
  model = mmtbx.model.manager(
    model_input = iotbx.pdb.input(
      source_info=None, lines=pdb_str_with_water.splitlines()))
  model.setup_scattering_dictionaries(scattering_table="wk1995")
  sites_frac = flex.vec3_double([(0.40, 0.40, 0.40), (0.60, 0.60, 0.60)])
  params = ordered_solvent.master_params().extract()
  ordered_solvent.add_solvent_to_model_inplace(
    sites=sites_frac, model=model, params=params)
  chain_ids = [c.id for c in model.get_hierarchy().only_model().chains()]
  assert chain_ids == ["A", "S"], chain_ids
  water_chain = model.get_hierarchy().only_model().chains()[-1]
  assert water_chain.residue_groups_size() == 6, \
    water_chain.residue_groups_size()
  # hierarchy and xray_structure must stay in the same atom order
  assert approx_equal(
    model.get_hierarchy().atoms().extract_xyz(),
    model.get_xray_structure().sites_cart())
  # no chain break inside the merged water chain
  assert "BREAK" not in model.model_as_pdb(), model.model_as_pdb()
  # one label_asym_id for all waters of chain S in mmCIF
  cif_block = model.get_hierarchy().as_cif_block()
  water_label_asym_ids = set(
    la for la, comp in zip(cif_block["_atom_site.label_asym_id"],
                           cif_block["_atom_site.label_comp_id"])
    if comp == "HOH")
  assert len(water_label_asym_ids) == 1, water_label_asym_ids

def exercise_added_solvent_water_chain_not_last():
  """
  When the existing water chain is not the last chain, new waters must still
  land at the end of the atom sequence (as a separate chain object), because
  add_solvent_to_model_inplace extends the xray_structure, the new-solvent
  selection and the refinement flags at the end.
  """
  pdb_str_water_then_ligand = pdb_str_with_water.replace("END\n", """\
HETATM    9 NA    NA L   1      18.000  18.000  18.000  1.00 20.00          NA
END
""")
  model = mmtbx.model.manager(
    model_input = iotbx.pdb.input(
      source_info=None, lines=pdb_str_water_then_ligand.splitlines()))
  model.setup_scattering_dictionaries(scattering_table="wk1995")
  sites_frac = flex.vec3_double([(0.40, 0.40, 0.40), (0.60, 0.60, 0.60)])
  params = ordered_solvent.master_params().extract()
  new_sel = ordered_solvent.add_solvent_to_model_inplace(
    sites=sites_frac, model=model, params=params)
  atoms = model.get_hierarchy().atoms()
  assert list(new_sel) == [False]*9 + [True]*2, list(new_sel)
  assert [a.parent().resname for a in atoms.select(new_sel)] == ["HOH"]*2
  assert approx_equal(
    atoms.extract_xyz(), model.get_xray_structure().sites_cart())

# Gln-Asn-Tyr with H atoms
peptide_h_str = """\
ATOM      1  N   GLN A   1      17.191   7.640  11.852  1.00 10.38           N
ATOM      2  CA  GLN A   1      18.307   8.180  12.621  1.00 11.39           C
ATOM      3  C   GLN A   1      19.631   7.587  12.153  1.00 11.52           C
ATOM      4  O   GLN A   1      19.805   6.365  12.035  1.00 12.05           O
ATOM      5  CB  GLN A   1      18.093   8.002  14.128  1.00 11.96           C
ATOM      6  CG  GLN A   1      16.866   8.769  14.628  1.00 10.81           C
ATOM      7  CD  GLN A   1      16.381   8.233  15.935  1.00 13.10           C
ATOM      8  OE1 GLN A   1      15.811   7.144  15.983  1.00 10.65           O
ATOM      9  NE2 GLN A   1      16.586   9.006  17.020  1.00 12.30           N
ATOM     10  H1  GLN A   1      17.464   7.449  11.027  1.00 10.38           H
ATOM     11  H2  GLN A   1      16.891   6.900  12.244  1.00 10.38           H
ATOM     12  H3  GLN A   1      16.537   8.242  11.816  1.00 10.38           H
ATOM     13  HA  GLN A   1      18.358   9.137  12.474  1.00 11.39           H
ATOM     14  HB2 GLN A   1      18.872   8.333  14.602  1.00 11.96           H
ATOM     15  HB3 GLN A   1      17.962   7.061  14.322  1.00 11.96           H
ATOM     16  HG2 GLN A   1      17.098   9.703  14.748  1.00 10.81           H
ATOM     17  HG3 GLN A   1      16.150   8.686  13.979  1.00 10.81           H
ATOM     18 HE21 GLN A   1      16.324   8.739  17.794  1.00 12.30           H
ATOM     19 HE22 GLN A   1      16.979   9.767  16.938  1.00 12.30           H
ATOM     20  N   ASN A   2      20.551   8.483  11.837  1.00 11.99           N
ATOM     21  CA  ASN A   2      21.868   8.129  11.299  1.00 12.30           C
ATOM     22  C   ASN A   2      22.891   8.580  12.305  1.00 13.40           C
ATOM     23  O   ASN A   2      23.256   9.762  12.355  1.00 13.92           O
ATOM     24  CB  ASN A   2      22.102   8.835   9.974  1.00 12.13           C
ATOM     25  CG  ASN A   2      20.998   8.554   8.984  1.00 12.77           C
ATOM     26  OD1 ASN A   2      20.835   7.423   8.532  1.00 14.27           O
ATOM     27  ND2 ASN A   2      20.232   9.566   8.660  1.00 10.07           N
ATOM     28  H   ASN A   2      20.440   9.331  11.925  1.00 11.99           H
ATOM     29  HA  ASN A   2      21.943   7.176  11.133  1.00 12.30           H
ATOM     30  HB2 ASN A   2      22.138   9.793  10.124  1.00 12.13           H
ATOM     31  HB3 ASN A   2      22.939   8.527   9.591  1.00 12.13           H
ATOM     32 HD21 ASN A   2      19.589   9.456   8.100  1.00 10.07           H
ATOM     33 HD22 ASN A   2      20.372  10.340   9.008  1.00 10.07           H
ATOM     34  N   TYR A   3      23.329   7.636  13.128  1.00 14.70           N
ATOM     35  CA  TYR A   3      24.196   7.963  14.280  1.00 15.18           C
ATOM     36  C   TYR A   3      25.640   8.150  13.866  1.00 15.91           C
ATOM     37  O   TYR A   3      26.078   7.630  12.836  1.00 15.76           O
ATOM     38  CB  TYR A   3      24.098   6.884  15.350  1.00 15.35           C
ATOM     39  CG  TYR A   3      22.702   6.748  15.883  1.00 14.45           C
ATOM     40  CD1 TYR A   3      22.247   7.575  16.901  1.00 14.80           C
ATOM     41  CD2 TYR A   3      21.808   5.840  15.308  1.00 15.68           C
ATOM     42  CE1 TYR A   3      20.941   7.468  17.397  1.00 14.33           C
ATOM     43  CE2 TYR A   3      20.517   5.725  15.777  1.00 13.46           C
ATOM     44  CZ  TYR A   3      20.084   6.548  16.812  1.00 15.09           C
ATOM     45  OH  TYR A   3      18.803   6.408  17.272  1.00 14.39           O
ATOM     46  OXT TYR A   3      26.395   8.818  14.593  1.00 17.49           O
ATOM     47  H   TYR A   3      23.145   6.799  13.051  1.00 14.70           H
ATOM     48  HA  TYR A   3      23.890   8.801  14.661  1.00 15.18           H
ATOM     49  HB2 TYR A   3      24.684   7.113  16.088  1.00 15.35           H
ATOM     50  HB3 TYR A   3      24.362   6.032  14.969  1.00 15.35           H
ATOM     51  HD1 TYR A   3      22.820   8.213  17.261  1.00 14.80           H
ATOM     52  HD2 TYR A   3      22.088   5.306  14.600  1.00 15.68           H
ATOM     53  HE1 TYR A   3      20.655   8.002  18.103  1.00 14.33           H
ATOM     54  HE2 TYR A   3      19.938   5.100  15.403  1.00 13.46           H
ATOM     55  HH  TYR A   3      18.433   7.160  17.325  1.00 14.39           H
"""

def get_xyz(lines, resname, name):
  for l in lines.splitlines():
    if l[17:20]==resname and l[12:16].strip()==name:
      return matrix.col([float(l[30+8*i:38+8*i]) for i in range(3)])
  raise RuntimeError(name)

def hoh_line(i, xyz, name=" O  ", element="O", chain="S", altloc=" ",
             occ=1.0):
  return "HETATM%5d %s%sHOH %s%4d    %8.3f%8.3f%8.3f%6.2f 20.00          %2s" % (
    900+i, name, altloc, chain, i, xyz[0], xyz[1], xyz[2], occ, element)

def water_with_h(i, o, away_from, chain="S"):
  """ Water O at o with two H pointing away from site away_from """
  u = (o - away_from).normalize()
  ref = matrix.col((0,0,1))
  if abs(u.dot(ref)) > 0.9: ref = matrix.col((1,0,0))
  v = u.cross(ref).normalize()
  a = math.radians(104.5/2)
  h1 = o + 0.96*(u*math.cos(a) + v*math.sin(a))
  h2 = o + 0.96*(u*math.cos(a) - v*math.sin(a))
  return [hoh_line(i, o, chain=chain),
          hoh_line(i, h1, name=" H1 ", element="H", chain=chain),
          hoh_line(i, h2, name=" H2 ", element="H", chain=chain)]

def get_model_h(water_lines, cryst1="CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1"):
  pdb = cryst1 + "\n" + peptide_h_str + "\n".join(water_lines) + "\n"
  return mmtbx.model.manager(
    model_input = iotbx.pdb.input(source_info=None, lines=pdb), log=null_out())

def kept_waters(model):
  return sorted([int(rg.resseq) for rg in model.get_hierarchy().residue_groups()
    if rg.atom_groups()[0].resname=="HOH"])

def exercise_filter_by_distance_h():
  """
  Distance filter with H atoms: water too close to H bonded to carbon is
  clash; polar H can H-bond with water; own H of water are ignored; water too
  close to non-H atom is removed also without altlocs.
  """
  P = peptide_h_str
  g = get_xyz
  def along(a, b, d): return b + (b-a).normalize()*d
  w1 = along(g(P,"TYR","CD1"), g(P,"TYR","HD1"), 2.0)
  lines = [
    # 1: 2.0 A from non-polar HD1 (2.9 A from CD1): removed
    hoh_line(1, w1),
    # 2: partner of 1, nothing close: kept
    hoh_line(2, along(g(P,"TYR","HD1"), w1, 2.8)),
    # 3: H-bond with polar HD21 (1.8 A), 2.8 A from ND2: kept
    hoh_line(3, along(g(P,"ASN","ND2"), g(P,"ASN","HD21"), 1.8)),
    # 4: 1.4 A from polar HE21: removed
    hoh_line(4, along(g(P,"GLN","NE2"), g(P,"GLN","HE21"), 1.4)),
    # 6: 1.5 A from OE1 (non-H): removed
    hoh_line(6, along(g(P,"GLN","CD"), g(P,"GLN","OE1"), 1.5))]
  # 5: water with own H, H-bond to OD1 (2.8 A): kept
  od1 = g(P,"ASN","OD1")
  lines.extend(water_with_h(5, along(g(P,"ASN","CG"), od1, 2.8),
    away_from=od1))
  m = get_model_h(lines)
  for altlocs in [False, True]:
    r = ordered_solvent.filter_by_distance(m.deep_copy(),
      fix_altlocs_and_filter_was_run=altlocs)
    assert kept_waters(r) == [2, 3, 5], kept_waters(r)
  # thresholds are parameters
  r = ordered_solvent.filter_by_distance(m.deep_copy(), False,
    dist_min_h_nonpolar=1.9, dist_min_h_polar=1.3, dist_min=1.4)
  assert kept_waters(r) == [1, 2, 3, 4, 5, 6], kept_waters(r)
  # fix_altlocs_and_filter: own H of water are not neighbors (used to get
  # altloc assigned or be removed)
  r = ordered_solvent.fix_altlocs_and_filter(m.deep_copy(), dist_min=1.8)
  assert 5 in kept_waters(r)

def exercise_no_duplicate_water_h_restraints():
  """ Adding water must not add restraints again for existing water with H """
  P = peptide_h_str
  oe1, cd = get_xyz(P,"GLN","OE1"), get_xyz(P,"GLN","CD")
  m = get_model_h(water_with_h(1, oe1 + (oe1-cd).normalize()*2.8,
    away_from=oe1))
  m.process(make_restraints=True)
  n_angles = m.get_restraints_manager().geometry.angle_proxies.size()
  params = ordered_solvent.master_params().extract()
  for k in range(2):
    ordered_solvent.add_solvent_to_model_inplace(model=m, params=params,
      sites=flex.vec3_double([(0.8+k*0.05, 0.8, 0.8)]))
  assert m.get_restraints_manager().geometry.angle_proxies.size() == n_angles

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

def get_manager_stub(model, fmodel, params):
  """ manager without running the pipeline, to call its steps directly """
  o = ordered_solvent.manager.__new__(ordered_solvent.manager)
  o.model, o.fmodel, o.params, o.log = model, fmodel, params, null_out()
  o._maps, o._model_state, o.n_water = None, None, None
  o.new_solvent_selection = None
  o.model_size_init = model.size()
  return o

def exercise_build_h():
  """
  H atoms are built for water in different chains (used to crash), water with
  alternative conformations is not touched (used to be duplicated).
  """
  cryst1 = "CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1"
  w1 = water_with_h(1, matrix.col((5,5,5)), matrix.col((5,5,0)), chain="A")
  w2 = water_with_h(2, matrix.col((5,14,5)), matrix.col((5,14,0)), chain="S")
  w3 = [hoh_line(3, matrix.col((14,4,4)), chain="S", altloc="A", occ=0.5),
        hoh_line(3, matrix.col((14,4.8,4)), chain="S", altloc="B", occ=0.5)]
  lines_true = w1 + w3 + w2
  lines = [l for l in lines_true if l[76:78].strip() != "H"]
  model_true = get_model_h(lines_true, cryst1=cryst1)
  model = get_model_h(lines, cryst1=cryst1)
  for m in [model_true, model]:
    m.setup_scattering_dictionaries(scattering_table="n_gaussian")
  fmodel = get_fmodel(model_true, model, d_min=0.9)
  params = ordered_solvent.master_params().extract()
  # at this resolution H peaks are closer to O and H-O-H angle is ~77
  params.add_hydrogens.min_angle = 60
  o = get_manager_stub(model, fmodel, params)
  o._get_maps()
  o._build_h()
  h = o.model.get_hierarchy()
  waters = {}
  for c in h.chains():
    for rg in c.residue_groups():
      if rg.atom_groups()[0].resname != "HOH": continue
      waters.setdefault(c.id, []).append(
        sorted([(ag.altloc, a.name.strip()) for ag in rg.atom_groups()
                for a in ag.atoms()]))
  assert waters["A"] == [[("", "H1"), ("", "H2"), ("", "O")]], waters
  assert sorted(waters["S"]) == sorted([
    [("", "H1"), ("", "H2"), ("", "O")],
    [("A", "O"), ("B", "O")]]), waters
  assert o.model.size() == model.size() + 4

def exercise_refine_new_water():
  """ New water is refined even if model size did not change """
  P = peptide_h_str
  m = get_model_h([])
  params = ordered_solvent.master_params().extract()
  o = get_manager_stub(m, None, params)
  calls = []
  o.refine_oat = lambda: calls.append(1)
  o.n_water = 1
  o.new_solvent_selection = flex.bool(m.size(), False)
  o._refine()
  assert calls == []
  o.new_solvent_selection[-1] = True
  o._refine()
  assert calls == [1]*params.n_cycles

def exercise_manager():
  """
  Run whole pipeline (find_peaks_params=None used to crash); maps are computed
  only when needed; no water clashes with H.
  """
  cryst1 = "CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1"
  P = peptide_h_str
  oe1, cd = get_xyz(P,"GLN","OE1"), get_xyz(P,"GLN","CD")
  oh, cz = get_xyz(P,"TYR","OH"), get_xyz(P,"TYR","CZ")
  lines_true = [hoh_line(1, oe1 + (oe1-cd).normalize()*2.8),
                hoh_line(2, oh + (oh-cz).normalize()*2.7)]
  model_true = get_model_h(lines_true, cryst1=cryst1)
  model = get_model_h(lines_true[:1], cryst1=cryst1)
  model.process(make_restraints=True)
  fmodel = get_fmodel(model_true, model, d_min=1.5)
  n_maps = []
  init = ordered_solvent.maps.__init__
  def counting_init(self, *args, **kwargs):
    n_maps.append(1)
    init(self, *args, **kwargs)
  ordered_solvent.maps.__init__ = counting_init
  try:
    o = ordered_solvent.manager(fmodel=fmodel, fmodels=None, model=model,
      is_neutron_scat_table=False, find_peaks_params=None, log=null_out())
  finally:
    ordered_solvent.maps.__init__ = init
  assert len(n_maps) <= 3, len(n_maps)
  assert o.n_water >= 1
  # no water closer than allowed to H
  m = o.model
  atoms = m.get_hierarchy().atoms()
  nb = ordered_solvent.get_water_neighbors(m, distance_cutoff=2.2)
  for i, neighbors in nb.items():
    for j, d in neighbors:
      if atoms[j].element_is_hydrogen():
        if ordered_solvent.is_polar_h(atoms[j]): assert d >= 1.6
        else:                                    assert d >= 2.2

def run():
  exercise_added_solvent_atom_name()
  exercise_added_solvent_unique_serials()
  exercise_added_solvent_merges_into_existing_chain()
  exercise_added_solvent_water_chain_not_last()
  exercise_filter_by_distance_h()
  exercise_no_duplicate_water_h_restraints()
  exercise_build_h()
  exercise_refine_new_water()
  exercise_manager()

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %6.2f" % (time.time() - t0))
