"""Tests of iotbx.pdb.hierarchy.root.sym_equiv_residue_groups."""
import iotbx.pdb
from cctbx import adptbx, crystal, sgtbx
from cctbx.array_family import flex
from cctbx.crystal import super_cell
from cctbx.crystal.sym_equiv_test_utils import (
  assert_same_positions, brute_force_sym_equiv, general_position_structure,
  special_position_fixtures, unique_positions)
from cctbx.development import debug_utils
from libtbx.test_utils import Exception_expected, approx_equal, show_diff
from libtbx.utils import Sorry
from scitbx import matrix
import sys

def residue_of_atoms(residue_sizes):
  """The residue index of each atom, residue_sizes[i] atoms in residue i."""
  result = []
  for i_rg, n_atoms in enumerate(residue_sizes):
    result.extend([i_rg] * n_atoms)
  return result

def hierarchy_from_sites(sites_cart, residue_sizes, elements=None):
  """A one-chain hierarchy of the sites, residue_sizes[i] atoms in residue i.
  Every atom is carbon unless elements says otherwise.
  """
  import mmtbx.model
  if (elements is None):
    elements = ["C"] * sites_cart.size()
  resseq_list = [iotbx.pdb.resseq_encode(i_rg + 1)
                 for i_rg in residue_of_atoms(residue_sizes)]
  return mmtbx.model.manager.from_sites_cart(
    sites_cart     = sites_cart,
    atom_name_list = [f"{element:>2}  " for element in elements],
    scatterer_list = elements,
    resseq_list    = resseq_list).get_hierarchy()

def brute_force_copies(sites_cart, crystal_symmetry, radius,
      primary_selection, residue_sizes):
  """Every copy of a residue with an atom that brute_force_sym_equiv finds.

  Returns a list of (residue index, operator as_xyz, copied positions), the
  identity left out.
  """
  unit_cell = crystal_symmetry.unit_cell()
  sites_frac = unit_cell.fractionalize(sites_cart=sites_cart)
  residue_of = residue_of_atoms(residue_sizes)
  copies = set()
  for j_seq, ops in brute_force_sym_equiv(
      sites_cart, crystal_symmetry, radius, primary_selection).items():
    for as_xyz in ops:
      copies.add((residue_of[j_seq], as_xyz))
  result = []
  for i_rg, as_xyz in sorted(copies):
    op = sgtbx.rt_mx(as_xyz)
    result.append((i_rg, as_xyz,
      [unit_cell.orthogonalize(op * sites_frac[i_seq])
       for i_seq in range(len(residue_of)) if residue_of[i_seq] == i_rg]))
  return result

def rotated_uij(unit_cell, rt_mx, uij):
  """uij under rt_mx, rotated as u_star in fractional coordinates."""
  r = matrix.sqr(rt_mx.r().as_double())
  u_star = matrix.sym(sym_mat3=adptbx.u_cart_as_u_star(unit_cell, uij))
  return adptbx.u_star_as_u_cart(
    unit_cell, (r * u_star * r.transpose()).as_sym_mat3())

def exercise_simple_example():
  """A small model and its sym_equiv copies within 4 A, as PDB strings.

  The two-fold of P 1 2 1 places the copy of N 3.2 A from N, so the
  residue is copied whole, C included; no equivalent of the water comes
  within 4 A.
  """
  pdb_inp = iotbx.pdb.input(source_info=None, lines="""\
CRYST1   20.000   15.000   20.000  90.00  90.00  90.00 P 1 2 1
ATOM      1  N   ALA A   1       1.500   5.000   0.500  1.00 10.00           N
ATOM      2  CA  ALA A   1       2.500   6.000   1.500  1.00 10.00           C
ATOM      3  C   ALA A   1       4.000   6.000   1.000  1.00 10.00           C
HETATM    4  O   HOH B   2       8.000   6.000   7.000  1.00 10.00           O
""")
  # All atoms are selected by default.
  copies = pdb_inp.construct_hierarchy().sym_equiv_residue_groups(
    pdb_inp.crystal_symmetry(), 4.0)
  # One entry per operator that copies a residue; the identity is never one.
  assert len(copies) == 1
  assert copies[0].rt_mx.as_xyz() == "-x,y,-z"
  # The atoms copied, as indices into the model's atoms(): all of ALA.
  assert list(copies[0].i_seqs) == [0, 1, 2]
  # Their copies, moved by rt_mx, keep chain A and residue 1.
  assert not show_diff(copies[0].hierarchy.as_pdb_string(), """\
ATOM      1  N   ALA A   1      -1.500   5.000  -0.500  1.00 10.00           N
ATOM      2  CA  ALA A   1      -2.500   6.000  -1.500  1.00 10.00           C
ATOM      3  C   ALA A   1      -4.000   6.000  -1.000  1.00 10.00           C
TER
""")

def exercise_match_brute_force(space_group_info):
  """A residue group with an atom whose equivalent is in range is copied
  whole under that operator, with its anisotropic U rotated and its sigxyz
  and siguij cleared, and nothing else is copied.

  The sites are in general positions, so no copied atom lands on the model.
  Every third atom is isotropic.
  """
  residue_sizes = [2, 3, 3, 4]
  residue_of = residue_of_atoms(residue_sizes)
  primary_selection = [0, 3, 7]
  radius = 5.0
  space_group_symbol = str(space_group_info)
  cs, sites_cart = general_position_structure(space_group_info)
  hierarchy = hierarchy_from_sites(sites_cart, residue_sizes)
  uijs = flex.sym_mat3_double()
  for i_seq in range(sites_cart.size()):
    if (i_seq % 3 == 2):
      uijs.append((-1, -1, -1, -1, -1, -1))
    else:
      uijs.append((0.30, 0.20, 0.25, 0.05 + 0.01 * i_seq, -0.04, 0.03))
  hierarchy.atoms().set_uij(uijs)
  hierarchy.atoms().set_sigxyz(
    flex.vec3_double(sites_cart.size(), (0.01, 0.02, 0.03)))
  if (iotbx.pdb.hierarchy.atom.has_siguij()):
    for atom in hierarchy.atoms():
      atom.set_siguij(new_siguij=(.001, .002, .003, .0004, .0005, .0006))
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
      assert atom.sigxyz == (0, 0, 0), info
      assert not atom.siguij_is_defined(), info
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
          expected.append(position)
      expected = unique_positions(expected, tolerance)
      assert len(expected) > 0, info
      assert_same_positions({0: expected}, {0: obtained}, info)

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
  none, and a Python list and more than one model are rejected.  Operators
  come in as_xyz order on denominators (1, 12), each copy keeps the chain
  ID, and i_seq is not relied on.
  """
  cs, sites_cart = general_position_structure()
  # copy_atoms gives every atom i_seq 0.
  hierarchy = hierarchy_from_sites(sites_cart, [3, 3, 3, 3]).select(
    flex.bool(sites_cart.size(), True), copy_atoms=True)
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
  for selection in [[0, 4, 9], []]:
    try:
      hierarchy.sym_equiv_residue_groups(cs, 5.0, selection=selection)
    except TypeError:
      pass
    else:
      raise Exception_expected
  two_models = hierarchy.deep_copy()
  two_models.append_model(hierarchy.models()[0].detached_copy())
  try:
    two_models.sym_equiv_residue_groups(cs, 5.0)
  except Sorry as e:
    assert str(e) == "Symmetry copies need one model, not 2.", str(e)
  else:
    raise Exception_expected

def run_call_back(flags, space_group_info):
  """The exercises that run once per space group."""
  exercise_match_brute_force(space_group_info)

def run():
  exercise_simple_example()
  exercise_leave_out_atoms_on_the_model()
  exercise_keep_atoms_of_another_element()
  exercise_on_an_element_get_every_copy()
  exercise_on_special_positions_match_brute_force()
  exercise_tolerance()
  exercise_contract()
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)
  print("OK")

if (__name__ == "__main__"):
  run()
