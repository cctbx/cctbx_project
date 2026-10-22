"""Tests of iotbx.pdb.hierarchy.root.sym_equiv_residue_groups."""
import iotbx.pdb.hierarchy
from cctbx import adptbx, crystal, sgtbx
from cctbx.array_family import flex
from cctbx.crystal import super_cell
from libtbx.test_utils import Exception_expected, approx_equal
from libtbx.utils import Sorry
from scitbx import matrix

# (space group, unit cell, fractional sites, selection): site 0 lies on a
# symmetry element.
special_position_fixtures = [
  ("C 1 2 1", (6, 7, 28, 90, 90, 90),                 # two-fold
   [(0.0, 0.2, 0.0), (0.13, 0.31, 0.07), (0.21, 0.44, 0.19)], [0, 1]),
  ("F m m 2", (9.673, 12.574, 16.444, 90, 90, 90),    # mirror, F centred
   [(0.0, 1 / 6., 0.0), (0.113, 0.207, 0.331), (0.717, 0.311, 0.613)], [1, 2]),
  ("R 3 :H", (12.1, 12.1, 16.7, 90, 90, 120),         # three-fold
   [(0.0, 0.0, 0.11), (0.213, 0.331, 0.07), (0.41, 0.62, 0.28)], [1, 2]),
  ("P 43 21 2", (12.1, 12.1, 16.7, 90, 90, 90),       # two-fold at (x, x, 0)
   [(0.19, 0.19, 0.0), (0.113, 0.407, 0.231), (0.617, 0.311, 0.513)], [1, 2]),
  ("P -1", (7, 8, 9, 90, 90, 90),                     # inversion centre
   [(0.0, 0.0, 0.0), (0.21, 0.33, 0.14), (0.44, 0.19, 0.37)], [1, 2]),
  ("P m -3 m", (9, 9, 9, 90, 90, 90),                 # site symmetry m-3m
   [(0.0, 0.0, 0.0), (0.213, 0.331, 0.07), (0.41, 0.28, 0.62)], [1, 2]),
]

def hierarchy_from_sites(sites_cart, residue_sizes, elements=None):
  """A one-chain hierarchy of the sites, residue_sizes[i] atoms in residue i,
  with i_seq left unset.  Every atom is carbon unless elements says otherwise.
  """
  if (elements is None):
    elements = ["C"] * sites_cart.size()
  root = iotbx.pdb.hierarchy.root()
  model = iotbx.pdb.hierarchy.model()
  chain = iotbx.pdb.hierarchy.chain(id="A")
  i_site = 0
  for i_rg, n_atoms in enumerate(residue_sizes):
    rg = iotbx.pdb.hierarchy.residue_group(resseq=f"{i_rg + 1:4d}")
    ag = iotbx.pdb.hierarchy.atom_group(resname="UNK")
    for i_atom in range(n_atoms):
      atom = iotbx.pdb.hierarchy.atom()
      atom.name = f"{elements[i_site]:>2}{i_atom + 1:<2d}"
      atom.element = f"{elements[i_site]:>2}"
      atom.xyz = sites_cart[i_site]
      ag.append_atom(atom)
      i_site += 1
    rg.append_atom_group(ag)
    chain.append_residue_group(rg)
  assert i_site == sites_cart.size()
  model.append_chain(chain)
  root.append_model(model)
  return root

def general_position_structure(space_group_symbol="P 21 21 21"):
  """Twelve sites, all in general positions, in a cell of 3000 A**3."""
  space_group_info = sgtbx.space_group_info(space_group_symbol)
  cs = crystal.symmetry(
    unit_cell=space_group_info.any_compatible_unit_cell(volume=3000.),
    space_group_info=space_group_info)
  sites_frac = flex.vec3_double()
  for i in range(12):
    sites_frac.append((
      (0.0417 + 0.0713 * i) % 0.83,
      (0.1103 + 0.0891 * i) % 0.79,
      (0.0629 + 0.1097 * i) % 0.87))
  return cs, cs.unit_cell().orthogonalize(sites_frac)

def brute_force_copies(sites_cart, crystal_symmetry, radius,
      primary_selection, residue_sizes, n_shifts=2):
  """Every copy of a residue with an atom within *radius* of a selected site,
  from all operators over a block of lattice translations, without
  asu_mappings.

  Returns a list of (residue index, operator as_xyz, copied positions), the
  identity left out.
  """
  unit_cell = crystal_symmetry.unit_cell()
  sites_frac = unit_cell.fractionalize(sites_cart=sites_cart)
  selected_fracs = [sites_frac[i_seq] for i_seq in primary_selection]
  identity_as_xyz = sgtbx.rt_mx().as_xyz()
  residues = []
  first = 0
  for n_atoms in residue_sizes:
    residues.append(range(first, first + n_atoms))
    first += n_atoms
  result = []
  for rotation_op in crystal_symmetry.space_group().all_ops():
    for i in range(-n_shifts, n_shifts + 1):
      for j in range(-n_shifts, n_shifts + 1):
        for k in range(-n_shifts, n_shifts + 1):
          op = rotation_op + (i, j, k)
          as_xyz = op.as_xyz()
          if (as_xyz == identity_as_xyz): continue
          for i_rg, residue in enumerate(residues):
            copy_fracs = [op * sites_frac[i_seq] for i_seq in residue]
            in_range = False
            for copy_frac in copy_fracs:
              for selected_frac in selected_fracs:
                if (unit_cell.distance(copy_frac, selected_frac) <= radius):
                  in_range = True
            if (in_range):
              result.append((i_rg, as_xyz,
                [unit_cell.orthogonalize(copy_frac)
                 for copy_frac in copy_fracs]))
  return result

def rotated_uij(unit_cell, rt_mx, uij):
  """uij under rt_mx, rotated as u_star in fractional coordinates."""
  r = matrix.sqr(rt_mx.r().as_double())
  u_star = matrix.sym(sym_mat3=adptbx.u_cart_as_u_star(unit_cell, uij))
  return adptbx.u_star_as_u_cart(
    unit_cell, (r * u_star * r.transpose()).as_sym_mat3())

def index_of_position(positions, position, tolerance):
  """Index of the first of *positions* within *tolerance* of *position* in
  every coordinate, or None."""
  for i, candidate in enumerate(positions):
    if (approx_equal(candidate, position, eps=tolerance, out=None)):
      return i
  return None

def exercise_match_brute_force(space_group_symbol):
  """A residue group with an atom whose equivalent is in range is copied
  whole under that operator, with its anisotropic U rotated, and nothing
  else is copied.

  The sites are in general positions, so no copied atom lands on the model.
  Every third atom is isotropic.
  """
  residue_sizes = [2, 3, 3, 4]
  residue_of = []
  for i_rg, n_atoms in enumerate(residue_sizes):
    residue_of.extend([i_rg] * n_atoms)
  primary_selection = [0, 3, 7]
  radius = 5.0
  cs, sites_cart = general_position_structure(space_group_symbol)
  hierarchy = hierarchy_from_sites(sites_cart, residue_sizes)
  uijs = flex.sym_mat3_double()
  for i_seq in range(sites_cart.size()):
    if (i_seq % 3 == 2):
      uijs.append((-1, -1, -1, -1, -1, -1))
    else:
      uijs.append((0.30, 0.20, 0.25, 0.05 + 0.01 * i_seq, -0.04, 0.03))
  hierarchy.atoms().set_uij(uijs)
  copies = hierarchy.sym_equiv_residue_groups(cs, radius,
    selection=flex.size_t(primary_selection))
  obtained = set()
  for copy in copies:
    i_rgs = set([residue_of[i_seq] for i_seq in copy.i_seqs])
    whole = [i_seq for i_seq in range(sites_cart.size())
             if residue_of[i_seq] in i_rgs]
    assert list(copy.i_seqs) == whole, (
      space_group_symbol, copy.rt_mx.as_xyz(), list(copy.i_seqs))
    assert approx_equal(
      copy.hierarchy.atoms().extract_xyz(),
      super_cell.sym_equiv_sites_cart(
        sites_cart, cs.unit_cell(), copy.rt_mx, copy.i_seqs),
      out=None), (space_group_symbol, copy.rt_mx.as_xyz())
    for i_seq, atom in zip(copy.i_seqs, copy.hierarchy.atoms()):
      info = (space_group_symbol, copy.rt_mx.as_xyz(), i_seq)
      if (i_seq % 3 == 2):
        assert not atom.uij_is_defined(), info
      else:
        assert approx_equal(atom.uij,
          rotated_uij(cs.unit_cell(), copy.rt_mx, uijs[i_seq]), out=None), info
    obtained.update([(i_rg, copy.rt_mx.as_xyz()) for i_rg in i_rgs])
  expected = set()
  for i_rg, as_xyz, _ in brute_force_copies(
      sites_cart, cs, radius, primary_selection, residue_sizes):
    expected.add((i_rg, as_xyz))
  assert len(expected) > 0, space_group_symbol
  assert obtained == expected, (
    space_group_symbol, sorted(expected - obtained),
    sorted(obtained - expected))

def exercise_leave_out_atoms_on_the_model():
  """A copied atom that lands on an atom of the model is left out.

  The two-fold of C 1 2 1 runs along (0, y, 0).  It swaps the two atoms of
  residue 1 and holds the first atom of residue 2.
  """
  cs = crystal.symmetry(unit_cell=(20, 24, 28, 90, 90, 90),
                        space_group_symbol="C 1 2 1")
  sites_cart = cs.unit_cell().orthogonalize(flex.vec3_double([
    (0.03, 0.2, 0.02), (-0.03, 0.2, -0.02),
    (0.0, 0.3, 0.0), (0.08, 0.31, 0.05),
    (0.1, 0.25, 0.06)]))
  copies = hierarchy_from_sites(sites_cart, [2, 2, 1]) \
    .sym_equiv_residue_groups(cs, 5.0)
  i_seqs_by_op = {}
  for copy in copies:
    i_seqs_by_op[copy.rt_mx.as_xyz()] = list(copy.i_seqs)
  assert i_seqs_by_op.get("-x,y,-z") == [3, 4], i_seqs_by_op
  n_copied = 0
  for copy in copies:
    for xyz in copy.hierarchy.atoms().extract_xyz():
      assert (flex.min(flex.sqrt((sites_cart - xyz).dot()))
              >= super_cell.min_distance_sym_equiv_coincident), (
        copy.rt_mx.as_xyz(), xyz)
      n_copied += 1
  assert n_copied > 0

def exercise_keep_atoms_of_another_element():
  """A copied atom next to an atom of another element is kept; next to one of
  its own element it is that atom.

  Under the two-fold of P 1 2 1, along (0, y, 0), each atom lands 0.22 A from
  the other.
  """
  cs = crystal.symmetry(unit_cell=(20, 15, 20, 90, 90, 90),
                        space_group_symbol="P 1 2 1")
  sites_cart = flex.vec3_double([(1.0, 5.0, 0.3), (-0.8, 5.0, -0.4)])
  for elements, expected in [(["C", "O"], [[0, 1]]), (["C", "C"], [])]:
    copies = hierarchy_from_sites(sites_cart, [1, 1], elements) \
      .sym_equiv_residue_groups(cs, 3.0, min_distance_sym_equiv=0.5)
    on_two_fold = [list(copy.i_seqs) for copy in copies
                   if copy.rt_mx.as_xyz() == "-x,y,-z"]
    assert on_two_fold == expected, (elements, on_two_fold)

def exercise_on_an_element_get_every_copy():
  """A residue group with an atom on a symmetry element is copied under every
  operator that places that atom in range, and the atom is copied once.

  Atom 0 lies on the two-fold of P 1 2 1 and is the only atom of residue 1
  in range of atom 2.  x,y+1,z and -x,y+1,-z both place it at (0, 20, 0),
  and place atom 1 apart.
  """
  cs = crystal.symmetry(unit_cell=(20, 15, 20, 90, 90, 90),
                        space_group_symbol="P 1 2 1")
  sites_cart = flex.vec3_double([(0, 5, 0), (1.4, 5, 0.3), (0.3, 18.5, 0.2)])
  copies = hierarchy_from_sites(sites_cart, [2, 1]) \
    .sym_equiv_residue_groups(cs, 1.6, selection=flex.size_t([2]))
  residue_1 = {}
  for copy in copies:
    i_seqs = [i_seq for i_seq in copy.i_seqs if i_seq < 2]
    if (i_seqs):
      residue_1[copy.rt_mx.as_xyz()] = i_seqs
  assert residue_1 == {"-x,y+1,-z": [0, 1], "x,y+1,z": [1]}, residue_1

def exercise_on_special_positions_match_brute_force():
  """Every copied atom a direct search finds is found once, and no others,
  for a residue group with an atom on a symmetry element.

  The direct search keeps each position once and none on a site of the
  model; all atoms are of one element.
  """
  residue_sizes = [2, 1]
  tolerance = super_cell.min_distance_sym_equiv_coincident
  for fixture in special_position_fixtures:
    space_group_symbol, unit_cell, fracs, selection = fixture
    cs = crystal.symmetry(unit_cell=unit_cell,
                          space_group_symbol=space_group_symbol)
    sites_cart = cs.unit_cell().orthogonalize(flex.vec3_double(fracs))
    hierarchy = hierarchy_from_sites(sites_cart, residue_sizes)
    for radius in [4.0, 6.0]:
      info = (space_group_symbol, radius)
      obtained = []
      for copy in hierarchy.sym_equiv_residue_groups(
          cs, radius, selection=flex.size_t(selection)):
        obtained.extend(copy.hierarchy.atoms().extract_xyz())
      expected = []
      for _, _, positions in brute_force_copies(
          sites_cart, cs, radius, selection, residue_sizes):
        for position in positions:
          if (flex.min(flex.sqrt((sites_cart - position).dot()))
              < tolerance):
            continue
          if (index_of_position(expected, position, tolerance) is None):
            expected.append(position)
      assert len(expected) > 0, info
      assert len(obtained) == len(expected), (
        info, len(expected), len(obtained))
      unmatched = list(obtained)
      for position in expected:
        i = index_of_position(unmatched, position, 1.e-4)
        assert i is not None, (info, position, unmatched)
        del unmatched[i]

def exercise_tolerance():
  """min_distance_sym_equiv decides when a copied atom is the atom it lands
  next to.

  Site 0 sits 0.2 A off the two-fold of C 1 2 1, 0.4 A from its copy under
  -x,y,-z.
  """
  cs = crystal.symmetry(unit_cell=(20, 24, 28, 90, 90, 90),
                        space_group_symbol="C 1 2 1")
  sites_cart = cs.unit_cell().orthogonalize(flex.vec3_double([
    (0.01, 0.2, 0.0), (0.13, 0.31, 0.07)]))
  hierarchy = hierarchy_from_sites(sites_cart, [1, 1])
  for tolerance, copied in [
      (super_cell.min_distance_sym_equiv_coincident, True), (0.5, False)]:
    copies = hierarchy.sym_equiv_residue_groups(
      cs, 5.0, selection=flex.size_t([0]),
      min_distance_sym_equiv=tolerance)
    on_two_fold = [list(copy.i_seqs) for copy in copies
                   if copy.rt_mx.as_xyz() == "-x,y,-z"]
    assert len(on_two_fold) == 1, (tolerance, on_two_fold)
    assert (0 in on_two_fold[0]) == copied, (tolerance, on_two_fold)

def exercise_contract():
  """A bool and an index selection give the same copies, an empty one gives
  none, and more than one model is rejected.  Operators come in as_xyz order
  on denominators (1, 12), each copy keeps the chain ID, and i_seq is not
  relied on.
  """
  cs, sites_cart = general_position_structure()
  hierarchy = hierarchy_from_sites(sites_cart, [3, 3, 3, 3])
  assert hierarchy.atoms().extract_i_seq().all_eq(0)
  primary_selection = flex.size_t([0, 4, 9])
  by_index = hierarchy.sym_equiv_residue_groups(
    cs, 5.0, selection=primary_selection)
  by_bool = hierarchy.sym_equiv_residue_groups(
    cs, 5.0, selection=flex.bool(sites_cart.size(), primary_selection))
  assert len(by_index) > 0
  assert ([(copy.rt_mx.as_xyz(), list(copy.i_seqs)) for copy in by_bool]
          == [(copy.rt_mx.as_xyz(), list(copy.i_seqs)) for copy in by_index])
  as_xyz = [copy.rt_mx.as_xyz() for copy in by_index]
  assert as_xyz == sorted(as_xyz), as_xyz
  for copy in by_index:
    assert copy.rt_mx.r().den() == 1, copy.rt_mx.as_xyz()
    assert copy.rt_mx.t().den() == 12, copy.rt_mx.as_xyz()
    assert copy.hierarchy.models_size() == 1
    assert [chain.id for chain in copy.hierarchy.chains()] == ["A"]
  assert hierarchy.sym_equiv_residue_groups(
    cs, 5.0, selection=flex.size_t()) == []
  two_models = hierarchy.deep_copy()
  two_models.append_model(hierarchy.models()[0].detached_copy())
  try:
    two_models.sym_equiv_residue_groups(cs, 5.0)
  except Sorry as e:
    assert str(e) == "Symmetry copies need one model, not 2.", str(e)
  else:
    raise Exception_expected

def run():
  for space_group_symbol in ["P 21 21 21", "I 41/a c d", "P 63", "R 3 :H"]:
    exercise_match_brute_force(space_group_symbol)
  exercise_leave_out_atoms_on_the_model()
  exercise_keep_atoms_of_another_element()
  exercise_on_an_element_get_every_copy()
  exercise_on_special_positions_match_brute_force()
  exercise_tolerance()
  exercise_contract()
  print("OK")

if (__name__ == "__main__"):
  run()
