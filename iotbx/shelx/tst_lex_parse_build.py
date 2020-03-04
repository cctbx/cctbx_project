from __future__ import absolute_import, division, print_function
from cctbx import crystal
from cctbx import uctbx
from cctbx import sgtbx
from cctbx import adptbx
from cctbx import xray
from iotbx import shelx
from iotbx.shelx import crystal_symmetry_from_ins
from iotbx.shelx import tokens
import iotbx.builders
from libtbx.test_utils import approx_equal, Exception_expected
from libtbx.math_utils import are_equivalent
from six.moves import cStringIO as StringIO
from six.moves import range
from six.moves import zip

def exercise_lexing():
  stream = shelx.command_stream(file=StringIO(ins_mundane_tiny))
  i = iter(stream)
  try:
    cmd, line = next(i)
    assert cmd == ('TITL', ('in Pbca',))
    cmd, line =  next(i)
    assert cmd == ('CELL', (0.71073, 7.35, 9.541, 12.842, 90, 90, 90))
    cmd, line = next(i)
    assert cmd == ('ZERR', (4, 0.002, 0.002, 0.003, 0, 0, 0))
    cmd, line = next(i)
    assert cmd == ('LATT', (1,))
    cmd, line = next(i)
    assert cmd == ('SYMM', ('0.5-X, -Y, 0.5+Z',))
    cmd, line = next(i)
    assert cmd == ('SYMM', ('-X, 0.5+Y, 0.5-Z',))
    cmd, line = next(i)
    assert cmd == ('SYMM', ('1/2+X, 0.5-Y, -Z',))
    cmd, line =  next(i)
    assert cmd == ('SFAC', ('C', 'H', 'O', 'N',))
    cmd, line = next(i)
    assert cmd == ('UNIT', (32, 40, 16, 8))
    cmd, line = next(i)
    assert cmd == ('TEMP', (-153,))
    cmd, line = next(i)
    assert cmd == ('L.S.', (4,))
    cmd, line = next(i)
    assert cmd == ('BOND', (tokens.element_token(element='H'),))
    cmd, line = next(i) # FMAP
    cmd, line = next(i) # PLAN
    cmd, line = next(i) # WGHT
    cmd, line = next(i) # EXTI
    cmd, line = next(i) # FVAR
    cmd, line = next(i)
    assert cmd == ('REM', ())
    cmd, line = next(i)
    assert cmd == ('+', '/path/to/filename.ins')
    cmd, line = next(i)
    assert cmd == ('REM', ('Protracted example of residues on command',))
    cmd, line = next(i)
    assert cmd == ('HFIX', (tokens.residue_number_tok, 1), (23,))
    cmd, line =  next(i)
    assert cmd == ('HFIX', (tokens.residue_class_tok, 'N'), (43,))
    cmd, line = next(i)
    assert cmd == ('EQIV', (1, '1-X, -Y, -Z'))
    cmd, line = next(i)
    assert cmd == ('CONF', (tokens.atomname_token(name='C4'),
                            tokens.atomname_token(name='N'),
                            tokens.atomname_token(name='H'),
                            tokens.atomname_token(name='O2', symmetry=1) ) )
    cmd, line = next(i)
    assert cmd == ('DFIX', (tokens.residue_number_tok, 1),
                   (1.5, tokens.atomname_token(name='C2'),
                    tokens.atomname_token(name='C3')))
    cmd, line = next(i)
    assert cmd == ('__ATOM__',
                   ('O2', 3, 0.362893, 0.160589, -0.035913, 11,
                          0.03926, 0.02517, 0.02140,
                          -0.00415, -0.00810, 0.01009))
    cmd, line = next(i)
    assert cmd == ('__ATOM__',
                   ('O3', 3, 0.696722, 0.119176, 0.260657, 11,
                          0.02838, 0.02133, 0.02918,
                          0.00011, -0.01030, -0.00048))
    cmd, line = next(i) # C1
    cmd, line = next(i) # C4
    cmd, line =  next(i)
    assert cmd == ('RESI', (1,))
    cmd, line = next(i) # C2
    cmd, line = next(i) # C3
    cmd, line =  next(i)
    assert cmd == ('RESI', ('N',))
    cmd, line = next(i) # N
    cmd, line = next(i) # HKLF
    try:
      cmd, line = next(i)
      raise AssertionError
    except StopIteration:
      pass
  except StopIteration:
    raise AssertionError

def exercise_lexing_bis():
  stream = shelx.command_stream(file=StringIO(ins_equal_sign_in_rem))
  i = iter(stream)
  try:
    cmd, line = next(i)
    assert cmd == ('REM',
                   ('Solution 1  R1  0.100,  Alpha = 0.0015  in P2(1)',))
    cmd, line = next(i)
    assert cmd == ('REM', ('C13 O10',))
    cmd, line = next(i)
    assert cmd == ('TITL', ('SUCROSE IN P2(1)',))
  except StopIteration:
    raise AssertionError

def exercise_crystal_symmetry_parsing():
  stream = shelx.command_stream(file=StringIO(ins_mundane_tiny))
  builder = iotbx.builders.crystal_symmetry_builder()
  stream = shelx.crystal_symmetry_parser(stream, builder)
  stream = shelx.wavelength_parser(stream.filtered_commands(), builder)
  stream.parse()
  assert builder.crystal_symmetry.is_similar_symmetry(
    crystal.symmetry(
      unit_cell=uctbx.unit_cell((7.350, 9.541, 12.842, 90, 90, 90)),
      space_group_symbol='Pbca'),
    relative_length_tolerance=1e-15,
    absolute_angle_tolerance=1e-15)
  cs = crystal_symmetry_from_ins.extract_from(
    file=StringIO(ins_mundane_tiny))
  assert cs.is_similar_symmetry(builder.crystal_symmetry)
  assert approx_equal(builder.wavelength_in_angstrom, 0.71073, eps=5e-6)

  stream = shelx.command_stream(file=StringIO(ins_P1))
  l = shelx.crystal_symmetry_parser(
    stream,
    builder=iotbx.builders.crystal_symmetry_builder())
  l.parse()
  assert l.builder.crystal_symmetry.is_similar_symmetry(
    crystal.symmetry(
      unit_cell=uctbx.unit_cell((1,2,3,99,100,101)),
      space_group_symbol='P1'),
    relative_length_tolerance=1e-15,
    absolute_angle_tolerance=1e-15)
  cs = crystal_symmetry_from_ins.extract_from(
    file=StringIO(ins_P1))
  assert cs.is_similar_symmetry(l.builder.crystal_symmetry)

def exercise_instruction_parsing():
  alternatives = (None,)
  try:
    _ = iotbx.builders.mixin_builder_class(
      "builder_to_test_instruction_parsing",
      iotbx.builders.reflection_data_source_builder,
      iotbx.builders.weighting_scheme_builder,
      iotbx.builders.twinning_builder)
    alternatives += (_(),)
  except AttributeError:
    pass
  for builder in alternatives:
    stream = shelx.command_stream(file=StringIO(ins_aspirin))
    l = shelx.instruction_parser(stream, builder)
    l.parse()
    ins = l.instructions
    assert ins['hklf']['s'] == 1
    assert ins['hklf']['matrix'].as_xyz() == 'z,x+y,x'
    assert ins['hklf']['n'] == 4
    assert ins['omit_hkl'] == [[2, 3, 4], [-1, -3, 2]]
    assert ins['omit']['s'] == -2
    assert ins['omit']['two_theta'] == 56
    assert ins['wght'] == {'a':0.0687,'b':0.4463,
                           'c':0, 'd':0, 'e':0, 'f':1/3}
    assert ins['merg'] == 2
    assert ins['twin']['matrix'].as_xyz() == '-x,y,-z'
    assert ins['twin']['n'] == 2
    assert ins['basf'] == (0.352,)
    assert ins['temp'] == -60
    if builder:
      assert builder.temperature_in_celsius == ins['temp']
      ws = builder.weighting_scheme
      assert isinstance(
        ws,
        iotbx.builders.least_squares.mainstream_shelx_weighting)
      assert ws.a == ins['wght']['a']
      assert ws.b == ins['wght']['b']
      assert builder.reflection_file_format == "hklf4"
      assert builder.data_change_of_basis_op.as_hkl()== "l,h+k,h"

  if len(alternatives) != 2: return
  builder = alternatives[-1]

  ins = StringIO(
    "HKLF 4 1  "
    "0.0000  0.0000  0.3330  1.0000  0.0000  0.0000  0.0000  1.0000 -0.3330")
  stream = shelx.command_stream(file=ins)
  stream = shelx.instruction_parser(stream, builder)
  stream.parse()
  assert builder.data_change_of_basis_op.as_xyz() == "y+3*z,x,y"

  ins = StringIO("HKLF 4 1 -1 2 0 -1 0 0 0 -1 1")
  stream = shelx.command_stream(file=ins)
  stream = shelx.instruction_parser(stream, builder)
  stream.parse()
  assert builder.data_change_of_basis_op.as_abc() == '-a+2*b,-a,-b+c'


def exercise_xray_structure_parsing():
  exercise_special_positions()
  exercise_aspirin()
  exercise_disordered()
  exercise_invalid()
  exercise_atom_with_peaks()
  exercise_q_peaks()

def exercise_atom_with_peaks():
  builder = iotbx.builders.crystal_structure_builder(set_grad_flags=True)
  stream = shelx.command_stream(
    file=StringIO(ins_with_atom_peak_heights))
  stream = shelx.crystal_symmetry_parser(stream, builder)
  stream = shelx.atom_parser(stream.filtered_commands(), builder=builder,
                             strictly_shelxl=False)
  stream.parse()
  sc = builder.structure.scatterers()

  assert sc[0].label == 'O001'
  assert approx_equal(sc[0].site, (0.60914, 0.62292, 0.82801), eps=1e-5)
  assert sc[0].occupancy == 1
  assert sc[0].flags.use_u_iso()
  assert approx_equal(sc[0].u_iso, 0.01991)

  assert sc[-1].label == 'C00N'
  assert approx_equal(sc[-1].site, (0.95213, 0.88631, 0.71108), eps=1e-5)
  assert sc[-1].occupancy == 1
  assert sc[-1].flags.use_u_iso()
  assert approx_equal(sc[-1].u_iso, 0.02971)

def exercise_q_peaks():
  builder = iotbx.builders.crystal_structure_builder(set_grad_flags=True)
  stream = shelx.command_stream(
    file=StringIO(ins_with_q_peaks))
  stream = shelx.crystal_symmetry_parser(stream, builder)
  stream = shelx.atom_parser(stream.filtered_commands(), builder=builder,
                             strictly_shelxl=False)
  stream.parse()
  uc = builder.structure.unit_cell()
  sc = builder.structure.scatterers()
  q = builder.electron_density_peaks

  assert len(sc) == 3

  assert sc[-1].label == 'O2'
  assert approx_equal(sc[-1].site, (0.645870, 0.883988, 0.489582), eps=1e-5)
  assert shelx_u_cif(uc, sc[-1].u_star) == (
    "0.01055   0.00848   0.00946   0.00007   0.00427   -0.00235")

  assert len(q) == 3

  assert approx_equal(q[0].site, (0.615600, 0.695100, 0.623900), eps=1e-6)
  assert approx_equal(q[0].height, 0.64, eps=1e-2)

  assert approx_equal(q[1].site, (0.526800, 0.678600, 0.431500), eps=1e-6)
  assert approx_equal(q[1].height, 0.63, eps=1e-2)

  assert approx_equal(q[2].site, (0.620800, 1.025100, 0.211700), eps=1e-6)
  assert approx_equal(q[2].height, 0.46, eps=1e-2)

def exercise_special_positions():
  structure = xray.structure.from_shelx(
    file=StringIO(ins_special_positions))
  occupancies = [ sc.occupancy for sc in structure.scatterers() ]
  multiplicities = [ sc.multiplicity() for sc in structure.scatterers() ]
  assert multiplicities == [ 2, 2, 2, 2, 6 ]
  assert approx_equal(occupancies, [ 1, 0.9999, 1, 1, 1 ], eps=1e-15)

def exercise_aspirin():
  for set_grad_flags in (False, True):
    structure = xray.structure.from_shelx(
      file=StringIO(ins_aspirin),
      set_grad_flags=set_grad_flags)
    isinstance(structure, xray.structure)
    assert structure.crystal_symmetry().is_similar_symmetry(
      crystal.symmetry(
        unit_cell=uctbx.unit_cell((11.492, 6.621, 11.453, 90, 95.615, 90)),
        space_group_symbol='P2(1)/c'),
      relative_length_tolerance=1e-15,
      absolute_angle_tolerance=1e-15)
    scatterers = structure.scatterers()
    unit_cell = structure.unit_cell()
    assert len(scatterers) == 21
    scatterer_labelled = dict([ (s.label, s) for s in scatterers ])
    o2 = scatterer_labelled['O2']
    assert o2.site == (0.879181, 0.140375, 0.051044)
    assert o2.occupancy == 1
    assert not o2.flags.use_u_iso() and o2.flags.use_u_aniso()
    assert shelx_u_cif(unit_cell, o2.u_star) == (
      "0.05893   0.05202   0.05656   0.01670   0.01559   0.01156")
      # copied from the .lst file
    h2 = scatterer_labelled['H2']
    assert h2.flags.use_u_iso() and not h2.flags.use_u_aniso()
    assert approx_equal(h2.u_iso, 0.08277, eps=1e-5)
    c3 = scatterer_labelled['C3']
    assert c3.flags.use_u_iso() and not c3.flags.use_u_aniso()
    h3 = scatterer_labelled['H3']
    assert approx_equal(h3.u_iso, 0.06952, eps=1e-5)
    h9a, h9b, h9c = [ scatterer_labelled[lbl]
                      for lbl in ('H9A', 'H9B', 'H9C') ]
    assert [ h.flags.use_u_iso() for h in (h9a, h9b, h9c) ]
    assert h9a.u_iso == h9b.u_iso == h9c.u_iso
    c9 = scatterer_labelled['C9']
    assert approx_equal(h9a.u_iso, 1.5*c9.u_iso_or_equiv(unit_cell))
    if not set_grad_flags:
      for scatt in scatterers:
        f = scatt.flags
        assert not f.grad_site()
        assert not f.grad_occupancy()
        assert not f.grad_u_iso()
        assert not f.grad_u_aniso()
        assert not f.grad_fp()
        assert not f.grad_fdp()
    else:
      for scatt in scatterers:
        f = scatt.flags
        assert f.grad_site() or scatt.label == 'C7'
        assert not f.grad_occupancy()
        assert are_equivalent(f.use_u_iso(),
                              f.grad_u_iso() and not f.grad_u_aniso())
        assert are_equivalent(f.use_u_aniso(),
                              not f.grad_u_iso() and f.grad_u_aniso())
        assert not f.grad_fp()
        assert not f.grad_fdp()

def exercise_disordered():
  for set_grad_flags in (False, True):
    builder = iotbx.builders.crystal_structure_builder(
      set_grad_flags=set_grad_flags)
    stream = shelx.command_stream(file=StringIO(ins_disordered))
    cs_parser = shelx.crystal_symmetry_parser(stream, builder)
    xs_parser = shelx.atom_parser(cs_parser.filtered_commands(), builder)
    xs_parser.parse()
    structure = builder.structure
    assert structure.crystal_symmetry().is_similar_symmetry(
      crystal.symmetry(
        unit_cell=uctbx.unit_cell((6.033, 6.830, 7.862,
                                   107.70, 103.17, 95.81)),
        space_group_symbol='P-1'),
      relative_length_tolerance=1e-15,
      absolute_angle_tolerance=1e-15)
    scatterers = structure.scatterers()
    unit_cell = structure.unit_cell()
    assert len(scatterers) == 15
    scatterer_labelled = dict([ (s.label, s) for s in scatterers ])
    br2, br3 = [ scatterer_labelled[x] for x in ('BR2', 'BR3') ]
    assert br2.occupancy == br3.occupancy
    assert approx_equal(br2.occupancy, 0.25054)
    if set_grad_flags:
      assert br2.flags.grad_occupancy()
      assert br3.flags.grad_occupancy()
    c12, h121, h122, h123 = [ scatterer_labelled[x]
                              for x in ('C12', 'H121', 'H122', 'H123') ]
    assert c12.occupancy == h121.occupancy == h122.occupancy == h123.occupancy
    assert approx_equal(c12.occupancy, 1-br2.occupancy)
    if set_grad_flags:
      for a in (c12, h121, h122, h123):
        assert a.flags.grad_occupancy()
    assert approx_equal(builder.conformer_indices,
                        [2, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2])
    assert builder.sym_excl_indices.count(0) == structure.scatterers().size()

def exercise_invalid():
  for set_grad_flags in (False, True):
    try:
      structure = xray.structure.from_shelx(
        file=StringIO(ins_invalid_scatt),
        set_grad_flags=set_grad_flags)
      raise Exception_expected
    except RuntimeError as e:
      assert str(e) == "ShelX: illegal argument '0.3.' at line 3"

    try:
      structure = xray.structure.from_shelx(
        file=StringIO(ins_invalid_scatt_1),
        set_grad_flags=set_grad_flags)
      raise Exception_expected
    except RuntimeError as e:
      assert str(e) == ("ShelX: wrong number of parameters "
                        "for scatterer at line 3")

    try:
      structure = xray.structure.from_shelx(
        file=StringIO(ins_missing_sfac),
        set_grad_flags=set_grad_flags)
      raise Exception_expected
    except RuntimeError as e:
      assert e.args[0].startswith('ShelX:')

  structure = xray.structure.from_shelx(
    file=StringIO(ins_disordered_with_part_sof))
  occ = 0.89064
  for sc in structure.scatterers():
    if sc.label in ('CL2', 'C28', 'H28A', 'H28B'):
      assert approx_equal(sc.occupancy, occ)
    elif sc.label in ('CL2B', 'C28B'):
      assert approx_equal(sc.occupancy, 1-occ)

def exercise_afix_parsing():
  import smtbx.refinement.constraints.geometrical.hydrogens as _
  builder = iotbx.builders.constrained_crystal_structure_builder()
  stream = shelx.command_stream(file=StringIO(ins_aspirin))
  l_cs = shelx.crystal_symmetry_parser(stream, builder)
  l_afix = shelx.afix_parser(l_cs.filtered_commands(), builder)
  l_xs = shelx.atom_parser(l_afix.filtered_commands(), builder)
  l_xs.parse()
  expected_geometrical_constraints = [
    _.staggered_terminal_tetrahedral_xh_site(
      constrained_site_indices=(1,),
      pivot=0),
    _.secondary_planar_xh_site(
      constrained_site_indices=(6,),
      pivot=5),
    _.secondary_planar_xh_site(
      constrained_site_indices=(10,),
      pivot=9),
    _.secondary_planar_xh_site(
      constrained_site_indices=(13,),
      pivot=12),
    _.secondary_planar_xh_site(
      constrained_site_indices=(16,),
      pivot=15),
    _.terminal_tetrahedral_xh3_sites(
      constrained_site_indices=(18,19,20),
      pivot=17,
      rotating=True)
    ]
  geometrical_constraints = [
    c for c in builder.constraints
    if c.__module__ == 'smtbx.refinement.constraints.geometrical.hydrogens' ]
  assert len(geometrical_constraints) == len(expected_geometrical_constraints)
  for result, expected in zip(geometrical_constraints,
                              expected_geometrical_constraints):
    assert (result == expected)

def exercise_u_iso_proportional_to_u_eq_parsing():
  import smtbx.refinement.constraints.adp as _
  builder = iotbx.builders.constrained_crystal_structure_builder()
  stream = shelx.command_stream(file=StringIO(ins_aspirin))
  l_cs = shelx.crystal_symmetry_parser(stream, builder)
  l_xs = shelx.atom_parser(l_cs.filtered_commands(), builder)
  l_xs.parse()
  expected_u_iso_constraints = [
    _.u_iso_proportional_to_pivot_u_eq(
      u_iso_scatterer_idx=1,
      u_eq_scatterer_idx=0,
      multiplier=1.5),
    _.u_iso_proportional_to_pivot_u_eq(
      u_iso_scatterer_idx=6,
      u_eq_scatterer_idx=5,
      multiplier=1.5),
    _.u_iso_proportional_to_pivot_u_eq(
      u_iso_scatterer_idx=10,
      u_eq_scatterer_idx=9,
      multiplier=1.5),
    _.u_iso_proportional_to_pivot_u_eq(
      u_iso_scatterer_idx=13,
      u_eq_scatterer_idx=12,
      multiplier=1.5),
    _.u_iso_proportional_to_pivot_u_eq(
      u_iso_scatterer_idx=16,
      u_eq_scatterer_idx=15,
      multiplier=1.5),
    _.u_iso_proportional_to_pivot_u_eq(
      u_iso_scatterer_idx=18,
      u_eq_scatterer_idx=17,
      multiplier=1.5),
    _.u_iso_proportional_to_pivot_u_eq(
      u_iso_scatterer_idx=19,
      u_eq_scatterer_idx=17,
      multiplier=1.5),
    _.u_iso_proportional_to_pivot_u_eq(
      u_iso_scatterer_idx=20,
      u_eq_scatterer_idx=17,
      multiplier=1.5),
  ]
  u_iso_constraints = [
    c for c in builder.constraints
    if c.__module__ == 'smtbx.refinement.constraints.adp' ]
  assert len(u_iso_constraints) == len(expected_u_iso_constraints)
  for result, expected in zip(u_iso_constraints, expected_u_iso_constraints):
    assert result == expected

def exercise_restraint_parsing():
  import smtbx.refinement.restraints
  def parse_restraints(ins_name):
    builder = iotbx.builders.restrained_crystal_structure_builder()
    stream = shelx.command_stream(file=StringIO(ins_name))
    l_cs = shelx.crystal_symmetry_parser(stream, builder)
    l_afix = shelx.afix_parser(l_cs.filtered_commands(), builder)
    l_xs = shelx.atom_parser(l_afix.filtered_commands(), builder)
    l_restraints = shelx.restraint_parser(l_xs.filtered_commands(), builder)
    l_restraints.parse()
    assert isinstance(builder.restraints_manager,
                      smtbx.refinement.restraints.manager)
    return builder.proxies()
  # exercise DFIX, DANG
  proxies = parse_restraints(ins_dfix_across_symm)
  shared_bond_proxy = proxies['bond']
  assert len(shared_bond_proxy) == 3
  assert approx_equal(shared_bond_proxy[0].distance_ideal, 1.75)
  assert approx_equal(shared_bond_proxy[1].distance_ideal, 1.75)
  assert approx_equal(shared_bond_proxy[2].distance_ideal, 1.75)
  assert shared_bond_proxy[0].i_seqs == (3,0)
  assert shared_bond_proxy[1].i_seqs == (3,2)
  assert shared_bond_proxy[2].i_seqs == (1,3)
  assert approx_equal(shared_bond_proxy[0].weight, 2500.0)
  assert approx_equal(shared_bond_proxy[1].weight, 1/(0.03*0.03))
  assert approx_equal(shared_bond_proxy[2].weight, 2500.0)
  assert shared_bond_proxy[0].rt_mx_ji == sgtbx.rt_mx('-x+1,y,-z+1/2')
  assert shared_bond_proxy[1].rt_mx_ji == sgtbx.rt_mx('-x+1,y,-z+1/2')
  assert shared_bond_proxy[2].rt_mx_ji is None
  proxies = parse_restraints(ins_dfix_multiple)
  shared_bond_proxy = proxies['bond']
  assert shared_bond_proxy.size() == 2
  # invalid DFIX instructions
  try:
    proxies = parse_restraints(ins_invalid_dfix)
  except RuntimeError as e:
    assert str(e) == "ShelX: Invalid DFIX instruction at line 3"
  try:
    proxies = parse_restraints(ins_invalid_dfix_2)
  except RuntimeError as e:
    assert str(e) == "ShelX: Invalid DFIX instruction at line 3"
  # exercise FLAT
  proxies = parse_restraints(ins_flat)
  shared_planarity_proxy = proxies['planarity']
  assert shared_planarity_proxy.size() == 2
  assert approx_equal(shared_planarity_proxy[0].i_seqs,
                      (0, 1, 2, 3, 4, 5, 6))
  assert approx_equal(shared_planarity_proxy[1].i_seqs,
                      (7, 8, 9, 10, 11, 12, 13))
  assert shared_planarity_proxy[0].sym_ops is None
  assert shared_planarity_proxy[1].sym_ops is None
  assert approx_equal(shared_planarity_proxy[0].weights,
                      (100, 100, 100, 100, 100, 100, 100))
  assert approx_equal(shared_planarity_proxy[1].weights,
                      (100, 100, 100, 100, 100, 100, 100))
  # invalid FLAT
  try:
    proxies = parse_restraints(ins_invalid_flat)
  except RuntimeError as e:
    assert str(e) == "ShelX: Invalid FLAT instruction at line 3"
  # SADI simple
  proxies = parse_restraints(ins_sadi)
  shared_bond_similarity_proxy\
      = proxies['bond_similarity']
  assert shared_bond_similarity_proxy.size() == 1
  assert approx_equal(
    shared_bond_similarity_proxy[0].i_seqs, ((0,2), (1,3)))
  assert shared_bond_similarity_proxy[0].sym_ops is None
  assert approx_equal(
    shared_bond_similarity_proxy[0].weights, (625.0, 625.0))
  # SADI with symmetry
  proxies = parse_restraints(ins_sadi_with_sym)
  shared_bond_similarity_proxy\
      = proxies['bond_similarity']
  assert shared_bond_similarity_proxy.size() == 1
  assert approx_equal(
    shared_bond_similarity_proxy[0].i_seqs, ((2,0), (1,3)))
  assert shared_bond_similarity_proxy[0].sym_ops\
    == (sgtbx.rt_mx('1-X,+Y,0.5-Z'),sgtbx.rt_mx('1-X,+Y,0.5-Z'))
  assert approx_equal(
    shared_bond_similarity_proxy[0].weights, (2500.0, 2500.0))
  # SIMU simple
  proxies = parse_restraints(ins_simu_simple)
  shared_apd_similarity_proxy \
      = proxies['adp_similarity']
  assert shared_apd_similarity_proxy.size() == 3
  assert shared_apd_similarity_proxy[0].i_seqs == (0,1)
  assert shared_apd_similarity_proxy[1].i_seqs == (1,2)
  assert shared_apd_similarity_proxy[2].i_seqs == (2,3)
  assert approx_equal(shared_apd_similarity_proxy[0].weight, 156.25)
  assert approx_equal(shared_apd_similarity_proxy[1].weight, 625.0)
  assert approx_equal(shared_apd_similarity_proxy[2].weight, 156.25)
  # SIMU with sigma given
  proxies = parse_restraints(ins_simu_s)
  shared_apd_similarity_proxy \
      = proxies['adp_similarity']
  assert shared_apd_similarity_proxy.size() == 3
  assert approx_equal(shared_apd_similarity_proxy[0].weight, 1/(0.06*0.06))
  assert approx_equal(shared_apd_similarity_proxy[1].weight, 1/(0.03*0.03))
  assert approx_equal(shared_apd_similarity_proxy[2].weight, 1/(0.06*0.06))
  # SIMU with sigma and sigma terminal
  proxies = parse_restraints(ins_simu_s_st)
  shared_apd_similarity_proxy \
      = proxies['adp_similarity']
  assert shared_apd_similarity_proxy.size() == 3
  assert approx_equal(shared_apd_similarity_proxy[0].weight, 1/(0.07*0.07))
  assert approx_equal(shared_apd_similarity_proxy[1].weight, 1/(0.03*0.03))
  assert approx_equal(shared_apd_similarity_proxy[2].weight, 1/(0.07*0.07))
  # SIMU with sigma and sigma terminal and dmax
  proxies = parse_restraints(ins_simu_s_st_dmax)
  shared_apd_similarity_proxy \
      = proxies['adp_similarity']
  assert shared_apd_similarity_proxy.size() == 3
  assert approx_equal(shared_apd_similarity_proxy[0].weight, 1/(0.07*0.07))
  assert approx_equal(shared_apd_similarity_proxy[1].weight, 1/(0.03*0.03))
  assert approx_equal(shared_apd_similarity_proxy[2].weight, 1/(0.07*0.07))
  # SIMU with just atoms
  proxies = parse_restraints(ins_simu_atoms)
  shared_apd_similarity_proxy \
      = proxies['adp_similarity']
  assert shared_apd_similarity_proxy.size() == 2
  assert shared_apd_similarity_proxy[0].i_seqs == (0,1)
  assert shared_apd_similarity_proxy[1].i_seqs == (1,2)
  assert approx_equal(shared_apd_similarity_proxy[0].weight, 156.25)
  assert approx_equal(shared_apd_similarity_proxy[1].weight, 625.0)
  # SIMU with sigma and atoms
  proxies = parse_restraints(ins_simu_s_atoms)
  shared_apd_similarity_proxy \
      = proxies['adp_similarity']
  assert shared_apd_similarity_proxy.size() == 2
  assert shared_apd_similarity_proxy[0].i_seqs == (0,1)
  assert shared_apd_similarity_proxy[1].i_seqs == (1,2)
  assert approx_equal(shared_apd_similarity_proxy[0].weight, 1/(0.06*0.06))
  assert approx_equal(shared_apd_similarity_proxy[1].weight, 1/(0.03*0.03))
  # SIMU with sigma, sigma terminal and atoms
  proxies = parse_restraints(ins_simu_s_st_atoms)
  shared_apd_similarity_proxy \
      = proxies['adp_similarity']
  assert shared_apd_similarity_proxy.size() == 2
  assert shared_apd_similarity_proxy[0].i_seqs == (0,1)
  assert shared_apd_similarity_proxy[1].i_seqs == (1,2)
  assert approx_equal(shared_apd_similarity_proxy[0].weight, 1/(0.07*0.07))
  assert approx_equal(shared_apd_similarity_proxy[1].weight, 1/(0.03*0.03))
  # SIMU with sigma, sigma terminal dmax and atoms
  proxies = parse_restraints(ins_simu_s_st_dmax_atoms)
  shared_apd_similarity_proxy \
      = proxies['adp_similarity']
  assert shared_apd_similarity_proxy.size() == 2
  assert shared_apd_similarity_proxy[0].i_seqs == (0,1)
  assert shared_apd_similarity_proxy[1].i_seqs == (1,2)
  assert approx_equal(shared_apd_similarity_proxy[0].weight, 1/(0.07*0.07))
  assert approx_equal(shared_apd_similarity_proxy[1].weight, 1/(0.03*0.03))
  # DELU simple
  proxies = parse_restraints(ins_delu_simple)
  shared_rigid_bond_proxy \
      = proxies['rigid_bond']
  assert shared_rigid_bond_proxy.size() == 5
  expected_i_seqs = ((0,1),(0,2),(1,2),(1,3),(2,3))
  for i in range(5):
    assert shared_rigid_bond_proxy[i].i_seqs == expected_i_seqs[i]
    assert approx_equal(shared_rigid_bond_proxy[i].weight, 10000)
  # DELU with s1
  proxies = parse_restraints(ins_delu_s1)
  shared_rigid_bond_proxy \
      = proxies['rigid_bond']
  assert shared_rigid_bond_proxy.size() == 5
  for i in range(5):
    assert shared_rigid_bond_proxy[i].i_seqs == expected_i_seqs[i]
    assert approx_equal(shared_rigid_bond_proxy[i].weight, 2500)
  # DELU with s1 and s2
  proxies = parse_restraints(ins_delu_s1_s2)
  shared_rigid_bond_proxy \
      = proxies['rigid_bond']
  assert shared_rigid_bond_proxy.size() == 5
  expected_weights = (10000,2500,10000,2500,10000)
  for i in range(5):
    assert shared_rigid_bond_proxy[i].i_seqs == expected_i_seqs[i]
    assert approx_equal(shared_rigid_bond_proxy[i].weight,
                        expected_weights[i])
  # DELU with atoms
  proxies = parse_restraints(ins_delu_atoms)
  shared_rigid_bond_proxy \
      = proxies['rigid_bond']
  assert shared_rigid_bond_proxy.size() == 3
  for i in range(3):
    assert shared_rigid_bond_proxy[i].i_seqs == expected_i_seqs[i]
    assert approx_equal(shared_rigid_bond_proxy[i].weight,10000)
  # DELU with s1 and atoms
  proxies = parse_restraints(ins_delu_s1_atoms)
  shared_rigid_bond_proxy \
      = proxies['rigid_bond']
  assert shared_rigid_bond_proxy.size() == 3
  for i in range(3):
    assert shared_rigid_bond_proxy[i].i_seqs == expected_i_seqs[i]
    assert approx_equal(shared_rigid_bond_proxy[i].weight,2500)
  # DELU with s1, s2 and atoms
  proxies = parse_restraints(ins_delu_s1_s2_atoms)
  shared_rigid_bond_proxy \
      = proxies['rigid_bond']
  assert shared_rigid_bond_proxy.size() == 3
  for i in range(3):
    assert shared_rigid_bond_proxy[i].i_seqs == expected_i_seqs[i]
    assert approx_equal(shared_rigid_bond_proxy[i].weight,
                        expected_weights[i])

  # RIGU simple
  proxies = parse_restraints(ins_rigu_simple)
  shared_rigu_proxy \
      = proxies['rigu']
  assert shared_rigu_proxy.size() == 5
  expected_i_seqs = ((0,1),(0,2),(1,2),(1,3),(2,3))
  for i in range(5):
    assert shared_rigu_proxy[i].i_seqs == expected_i_seqs[i]
    assert approx_equal(shared_rigu_proxy[i].weight, 62500)
  # RIGU with s1
  proxies = parse_restraints(ins_rigu_s1)
  shared_rigu_proxy \
      = proxies['rigu']
  assert shared_rigu_proxy.size() == 5
  for i in range(5):
    assert shared_rigu_proxy[i].i_seqs == expected_i_seqs[i]
    assert approx_equal(shared_rigu_proxy[i].weight, 2500)
  # RIGU with s1 and s2
  proxies = parse_restraints(ins_rigu_s1_s2)
  shared_rigu_proxy \
      = proxies['rigu']
  assert shared_rigu_proxy.size() == 5
  expected_weights = (10000,2500,10000,2500,10000)
  for i in range(5):
    assert shared_rigu_proxy[i].i_seqs == expected_i_seqs[i]
    assert approx_equal(shared_rigu_proxy[i].weight,
                        expected_weights[i])
  # RIGU with atoms
  proxies = parse_restraints(ins_rigu_atoms)
  shared_rigu_proxy \
      = proxies['rigu']
  assert shared_rigu_proxy.size() == 3
  for i in range(3):
    assert shared_rigu_proxy[i].i_seqs == expected_i_seqs[i]
    assert approx_equal(shared_rigu_proxy[i].weight,62500)
  # RIGU with s1 and atoms
  proxies = parse_restraints(ins_rigu_s1_atoms)
  shared_rigu_proxy \
      = proxies['rigu']
  assert shared_rigu_proxy.size() == 3
  for i in range(3):
    assert shared_rigu_proxy[i].i_seqs == expected_i_seqs[i]
    assert approx_equal(shared_rigu_proxy[i].weight,2500)
  # RIGU with s1, s2 and atoms
  proxies = parse_restraints(ins_rigu_s1_s2_atoms)
  shared_rigu_proxy \
      = proxies['rigu']
  assert shared_rigu_proxy.size() == 3
  for i in range(3):
    assert shared_rigu_proxy[i].i_seqs == expected_i_seqs[i]
    assert approx_equal(shared_rigu_proxy[i].weight,
                        expected_weights[i])

  # ISOR simple
  proxies = parse_restraints(ins_isor_simple)
  shared_isotropic_adp_proxy \
      = proxies['isotropic_adp']
  assert shared_isotropic_adp_proxy.size() == 4
  expected_i_seqs = (0,1,2,3)
  expected_weights = (25,100,100,25)
  for i in range(4):
    assert shared_isotropic_adp_proxy[i].i_seqs[0] == expected_i_seqs[i]
    assert approx_equal(shared_isotropic_adp_proxy[i].weight,
                        expected_weights[i])
  # ISOR with s
  proxies = parse_restraints(ins_isor_s)
  shared_isotropic_adp_proxy \
      = proxies['isotropic_adp']
  assert shared_isotropic_adp_proxy.size() == 4
  expected_weights = (6.25,25,25,6.25)
  for i in range(4):
    assert shared_isotropic_adp_proxy[i].i_seqs[0] == expected_i_seqs[i]
    assert approx_equal(shared_isotropic_adp_proxy[i].weight,
                        expected_weights[i])
  # ISOR with s and st
  proxies = parse_restraints(ins_isor_s_st)
  shared_isotropic_adp_proxy \
      = proxies['isotropic_adp']
  assert shared_isotropic_adp_proxy.size() == 4
  expected_weights = (1/(0.3*0.3), 25, 25, 1/(0.3*0.3))
  for i in range(4):
    assert shared_isotropic_adp_proxy[i].i_seqs[0] == expected_i_seqs[i]
    assert approx_equal(shared_isotropic_adp_proxy[i].weight,
                        expected_weights[i])
  # ISOR with atoms
  proxies = parse_restraints(ins_isor_atoms)
  shared_isotropic_adp_proxy \
      = proxies['isotropic_adp']
  assert shared_isotropic_adp_proxy.size() == 3
  expected_i_seqs = (0,1,2)
  expected_weights = (25,100,100)
  for i in range(3):
    assert shared_isotropic_adp_proxy[i].i_seqs[0] == expected_i_seqs[i]
    assert approx_equal(shared_isotropic_adp_proxy[i].weight,
                        expected_weights[i])
  # ISOR with s and atoms
  proxies = parse_restraints(ins_isor_s_atoms)
  shared_isotropic_adp_proxy \
      = proxies['isotropic_adp']
  assert shared_isotropic_adp_proxy.size() == 3
  expected_i_seqs = (0,1,2)
  expected_weights = (6.25,25,25)
  for i in range(3):
    assert shared_isotropic_adp_proxy[i].i_seqs[0] == expected_i_seqs[i]
    assert approx_equal(shared_isotropic_adp_proxy[i].weight,
                        expected_weights[i])
  # ISOR with s, st and atoms
  proxies = parse_restraints(ins_isor_s_st_atoms)
  shared_isotropic_adp_proxy \
      = proxies['isotropic_adp']
  assert shared_isotropic_adp_proxy.size() == 3
  expected_i_seqs = (0,1,2)
  expected_weights = (1/(0.3*0.3), 25, 25)
  for i in range(3):
    assert shared_isotropic_adp_proxy[i].i_seqs[0] == expected_i_seqs[i]
    assert approx_equal(shared_isotropic_adp_proxy[i].weight,
                        expected_weights[i])

def exercise_residues():
  proxies = parse_restraints(ins_residues)
  planarity_proxies = proxies['planarity']
  assert planarity_proxies.size() == 2
  p = planarity_proxies[1]
  assert approx_equal(p.i_seqs, [3, 1, 4, 2, 5])
  assert approx_equal(p.weights, [4]*len(p.weights))
  p = planarity_proxies[0]
  assert approx_equal(p.i_seqs, [15, 16, 17, 18, 12, 13, 14, 11, 10, 9, 8])
  assert approx_equal(p.weights, [100]*len(p.weights))
  assert p.sym_ops is None
  bond_proxies = proxies['bond']
  assert len(bond_proxies) == 10
  p = bond_proxies[0]
  assert p.distance_ideal == 2.24
  assert p.i_seqs == (9, 11)
  assert approx_equal(p.weight, 100)
  p = bond_proxies[2]
  assert p.distance_ideal == 2.425
  assert p.i_seqs == (1, 4)
  assert approx_equal(p.weight, 1479.2899408284025)
  p = bond_proxies[3]
  assert p.distance_ideal == 2.250
  assert p.i_seqs == (3, 4)
  assert approx_equal(p.weight, 3460.207612456747)
  p = bond_proxies[6]
  assert p.i_seqs == (2, 5)
  assert approx_equal(p.distance_ideal, 2.435)
  assert approx_equal(p.weight, 625.0)
  p = bond_proxies[9]
  assert p.i_seqs == (2, 4)
  assert approx_equal(p.distance_ideal, 1.329)
  assert approx_equal(p.weight, 2500.0)
  adp_sim_proxies = proxies['adp_similarity']
  assert len(adp_sim_proxies) == 19
  p = adp_sim_proxies[0]
  assert p.i_seqs == (0, 1)
  assert approx_equal(p.weight, 25)
  rigid_bond = proxies['rigid_bond']
  assert rigid_bond.size() == 43
  p = rigid_bond[1]
  assert p.i_seqs == (0, 2)
  assert approx_equal(p.weight, 10000)
  isot = proxies['isotropic_adp']
  assert isot.size() == 4
  p = isot[0]
  assert p.i_seqs == (19,)
  assert approx_equal(p.weight, 25)
  p = isot[-1]
  assert p.i_seqs == (22,)
  assert approx_equal(p.weight, 25)




def parse_restraints(ins_name):
  import smtbx.refinement.restraints
  builder = iotbx.builders.restrained_crystal_structure_builder()
  stream = shelx.command_stream(file=StringIO(ins_name))
  l_cs = shelx.crystal_symmetry_parser(stream, builder)
  l_afix = shelx.afix_parser(l_cs.filtered_commands(), builder)
  l_xs = shelx.atom_parser(l_afix.filtered_commands(), builder)
  l_restraints = shelx.restraint_parser(
    l_xs.filtered_commands(), builder)
  l_restraints.parse()
  builder = l_restraints.builder
  assert isinstance(builder.restraints_manager,
                    smtbx.refinement.restraints.manager)
  return builder.proxies()

def shelx_u_cif(unit_cell, u_star):
  u_cif = adptbx.u_star_as_u_cif(unit_cell, u_star)
  u_cif = u_cif[0:3] + (u_cif[-1], u_cif[-2], u_cif[-3])
  return (" "*3).join([ "%.5f" % x for x in  u_cif ])

def exercise_constrained_occupancies():
  from smtbx.refinement.constraints.occupancy import \
    occupancy_pair_affine_constraint
  builder = iotbx.builders.constrained_crystal_structure_builder()
  stream = shelx.command_stream(file=StringIO(ins_thpp))
  l_cs = shelx.crystal_symmetry_parser(stream, builder)
  l_xs = shelx.atom_parser(l_cs.filtered_commands(), builder)
  l_xs.parse()
  assert len(builder.constraints) == 2
  sc = builder.structure.scatterers()
  dsu = [(sc[c.scatterer_indices[0]].label, c) for c in builder.constraints]
  dsu.sort()
  constraints = [c for _,c in dsu]
  c0, c1 = constraints
  assert isinstance(c0, occupancy_pair_affine_constraint)
  assert c0.scatterer_indices == (len(sc) - 4, len(sc) - 2)
  assert c0.linear_form == ((1,1), 1)
  assert isinstance(c1, occupancy_pair_affine_constraint)
  assert c1.scatterer_indices == (3, len(sc) - 1)
  assert c1.linear_form == ((1,1), 1)

def exercise_wavelength():
  from iotbx.shelx.parsers import wavelength_parser
  builder = iotbx.builders.crystal_structure_builder()
  stream = shelx.command_stream(file=StringIO(ins_aspirin))
  stream = wavelength_parser(stream, builder)
  stream.parse()
  assert approx_equal(builder.wavelength_in_angstrom, 0.71073)

def run():
  exercise_lexing()
  exercise_lexing_bis()
  exercise_instruction_parsing()
  import libtbx.load_env
  if (not libtbx.env.has_module(name="smtbx")):
    print("Skipping some tests: smtbx module is not available.")
  else:
    exercise_restraint_parsing()
    exercise_constrained_occupancies()
    exercise_u_iso_proportional_to_u_eq_parsing()
    exercise_afix_parsing()
    exercise_residues()
  exercise_xray_structure_parsing()
  exercise_crystal_symmetry_parsing()
  exercise_wavelength()
  print('OK')

ins_mundane_tiny = (
"TITL in Pbca\n"
"CELL0.71073   7.350   9.541  12.842  90.000  90.000  90.000\n"
"ZERR    4.00   0.002   0.002   0.003   0.000   0.000   0.000\n"
"LATT  1\n"
"\n"
"SYMM 0.5-X, -Y, 0.5+Z\n"
"SYMM -X, 0.5+Y, 0.5-Z\n"
"SYMM 1/2+X, 0.5-y, -Z\n"
"SFAC C H O N\n"
"UNIT 32 40 16 8\n"
"TEMP -153\n"
" skipped\n"
"L.S. 4\n"
"BOND $H\n"
"FMAP 2 ! Fo - Fc\n"
"PLAN 5\n"
"WGHT    0.054500    0.220600\n"
"EXTI    0.034016\n"
"FVAR       0.89220\n"
"REM\n"
"+/path/to/filename.ins\n"
"REM  Protracted example of residues on command\n"
"HFIX_1 23\n"
"HFIX_N 43\n"
"EQIV $1 1-X, -Y, -Z\n"
"CONF C4 N H O2_$1\n"
"DFIX_1 1.5 C2 C3\n"
"O2    3    0.362893    0.160589   -0.035913    11.00000    0.03926    0.02517 =\n"
"         0.02140   -0.00415   -0.00810    0.01009\n"
"O3    3    0.696722    0.119176    0.260657    = some garbage\n"
" 11.00000    0.02838    0.02133 =\n"
"         0.02918    0.00011   -0.01030   -0.00048\n"
"C1    1    0.432360    0.196651    0.046472    11.00000    0.02226    0.01792 =\n"
"         0.01857    0.00058    0.00123    0.00249\n"
"C4    1    0.601261    0.176489    0.196615    11.00000    0.01844    0.01705 =\n"
"         0.02123   -0.00033   -0.00071   -0.00284\n"
"RESI 1\n"
"C2    1    0.411568    0.336421    0.100062    11.00000    0.02867    0.01555 =\n"
"         0.02069    0.00053   -0.00070    0.00268\n"
"C3    1    0.522448    0.322965    0.200394    11.00000    0.02158    0.01556 =\n"
"         0.02255   -0.00189   -0.00049   -0.00115\n"
"RESI N\n"
"N     4    0.544830    0.114061    0.104838    11.00000    0.02004    0.01586 =\n"
"         0.01943   -0.00094   -0.00096    0.00161\n"
" \n"
"HKLF 4 \n"
"END\n"
"     \n"
"WGHT      0.0536      0.2259 \n"
"REM Highest difference peak  0.330,  deepest hole -0.210,  1-sigma level  0.046\n"
"Q1    1   0.5622  0.2496  0.2006  11.00000  0.05    0.33\n"
"Q2    1   0.4189  0.2682  0.0749  11.00000  0.05    0.33\n"
"Q3    1   0.5706  0.1441  0.1462  11.00000  0.05    0.29\n"
"Q4    1   0.4966  0.1552  0.0716  11.00000  0.05    0.29\n"
"Q5    1   0.4666  0.3316  0.1506  11.00000  0.05    0.26\n")

ins_special_positions = """TITL cr4 in P6(3)
CELL 0.71073   5.1534   5.1534   8.6522  90.000  90.000 120.000
ZERR    2.00   0.0007   0.0007   0.0017   0.000   0.000   0.000
LATT -1
SYMM -Y, X-Y, Z
SYMM -X+Y, -X, Z
SYMM -X, -Y, 0.5+Z
SYMM Y, -X+Y, 0.5+Z
SYMM X-Y, X, 0.5+Z
SFAC LI O S K
UNIT 2 8 2 2

K1    4    0.000000    0.000000   -0.001950    10.33333    0.02427    0.02427 =
         0.02379    0.00000    0.00000    0.01214
S1    3    0.333333    0.666667    0.204215    10.33330    0.01423    0.01423 =
         0.01496    0.00000    0.00000    0.00712
LI1   1    0.333333    0.666667    0.815681    10.33333    0.02132    0.02132 =
         0.02256    0.00000    0.00000    0.01066
O1    2    0.333333    0.666667    0.035931    10.33333    0.06532    0.06532 =
         0.01669    0.00000    0.00000    0.03266
O2    2    0.343810    0.941658    0.258405    11.00000    0.02639    0.02079 =
         0.05284   -0.01180   -0.00053    0.01194
HKLF 4
"""

ins_P1 = (
"TITL P1\n"
"CELL 0.71073   1 2 3 99 100 101\n"
"ZERR    4.00   0.002   0.002   0.003   0.000   0.000   0.000\n"
"SFAC C H O N\n")

ins_aspirin = """TITL aspirin in P2(1)/c
CELL 0.71073 11.492 6.621 11.453 90.0 95.615 90.0
ZERR 4 0.002 0.001 0.002 0.0 0.0 0.0
LATT 1
SYMM -X,0.500+Y,0.500-Z
SFAC C H O
UNIT 36 28 16
TEMP -60.000
L.S. 20
WGHT    0.068700    0.446300
MERG 2
OMIT -2 56
OMIT 2 3 4
OMIT -1 -3 2
BASF 0.352
REM this twin command is nonsense...
TWIN -1 0 0 0 1 0 0 0 -1 2
FVAR       0.91641
O2    3    0.879181    0.140375    0.051044    11.00000    0.05893    0.05202 =
         0.05656    0.01670    0.01559    0.01156
AFIX 87
H2    2    0.925036    0.046729    0.065472    11.00000   -1.50000
AFIX   0
O3    3    0.714314    0.411976    0.087940    11.00000    0.04930    0.05061 =
         0.03685   -0.00419    0.00783   -0.00427
C2    1    0.846429    0.436128   -0.067477    11.00000    0.03922    0.03869 =
         0.03570   -0.00268    0.00014   -0.00150
O1    3    0.989649    0.187859   -0.096470    11.00000    0.06532    0.05656 =
         0.05359    0.00743    0.02230    0.01823
C3    1    0.883500    0.550173   -0.159841    11.00000    0.04635
AFIX  43
H3    2    0.945193    0.503717   -0.199359    11.00000   -1.50000
AFIX   0
C7    1    10.910013    0.242812   -0.037585    11.00000    0.04340    0.04037 =
         0.03632   -0.00344    0.00282   -0.00157
C1    1    0.754688    0.511841   -0.009478    11.00000    0.04076
C6    1    0.701446    0.694836   -0.041608    11.00000    0.05003    0.04601 =
         0.05788   -0.00711    0.00597    0.00688
AFIX  43
H6    2    0.640709    0.743710   -0.001563    11.00000   -1.50000
AFIX   0
O4    3    0.596499    0.219939   -0.034410    11.00000    0.07409    0.08015 =
         0.04201   -0.00406    0.00247   -0.02748
C5    1    0.739822    0.803914   -0.134050    11.00000    0.05726    0.04061 =
         0.06325    0.00538    0.00340    0.00642
AFIX  43
H5    2    0.704633    0.926467   -0.155929    11.00000   -1.50000
AFIX   0
C8    1    0.634563    0.262008    0.063872    11.00000    0.04357    0.05241 =
         0.04294   -0.00177    0.00752   -0.00018
C4    1    0.829521    0.731819   -0.193418    11.00000    0.05444    0.04790 =
         0.05321    0.00900    0.00503   -0.00209
AFIX  43
H4    2    0.854077    0.804532   -0.256010    11.00000   -1.50000
AFIX   0
C9    1    0.602602    0.163875    0.172952    11.00000    0.06769    0.07924 =
         0.04765    0.01171    0.00887   -0.00967
AFIX 137
H9A   2    0.571182    0.263181    0.222267    11.00000   -1.50000
H9C   2    0.670926    0.103865    0.213827    11.00000   -1.50000
H9B   2    0.545047    0.061174    0.153241    11.00000   -1.50000
HKLF 4 1 0 0 1 1 1 0 1 0 0

REM  aspirin in P2(1)/c
REM R1 =  0.0455 for   1038 Fo > 4sig(Fo)  and  0.0990 for all   1806 data
REM    110 parameters refined using      0 restraints

END

WGHT      0.0534      0.3514
"""

ins_disordered = """\
TITL 97srv101 in P-1
CELL 0.71073   6.033   6.830   7.862 107.70 103.17  95.81
ZERR    1.00   0.002   0.001   0.004   0.02   0.01   0.01
LATT  1
SFAC C H S BR
UNIT 9 9 4 1
TEMP -123
SIZE 0.25 0.15 0.04
OMIT -3 58
ACTA
L.S. 5
DFIX 1.5 C2 C12 C3 C13
BOND
FMAP 2
PLAN 10
WGHT    0.035100    0.787700
FVAR       0.70431   0.25054
PART 2
BR2   4    0.210409    0.205854    0.090628    21.00000    0.01590    0.03384 =
         0.03235    0.00702   -0.00500   -0.00273
PART   1
BR3   4    0.319774    0.727684    0.184240    21.00000    0.01979    0.03679 =
         0.03533    0.01864    0.00196    0.01154
PART   0
S1    3    0.695496    0.257529    0.325940    11.00000    0.01668    0.02391 =
         0.02776    0.00712    0.00046    0.00401
S2    3    0.795402    0.715230    0.403850    11.00000    0.01702    0.02605 =
         0.03256    0.01155    0.00143    0.00484
C1    1    0.893971    0.494238    0.443049    11.00000    0.01446    0.02402 =
         0.02423    0.01012    0.00393    0.00234
C2    1    0.470694    0.372115    0.231574    11.00000    0.01981    0.03196 =
         0.01877    0.00819    0.00393    0.00676
C3    1    0.516311    0.578190    0.267351    11.00000    0.01544    0.03197 =
         0.02369    0.01001    0.00228    0.00481
PART 1
C12   1    0.241589    0.221595    0.109681   -21.00000    0.03323
AFIX   7
H121  2    0.184843    0.141885    0.180774   -21.00000   -1.50000
H122  2    0.269551    0.125107   -0.000992   -21.00000   -1.50000
H123  2    0.125620    0.301975    0.072102   -21.00000   -1.50000
AFIX   0
PART   2
C13   1    0.349281    0.712664    0.197083   -21.00000    0.03570
AFIX   7
H131  2    0.336393    0.828845    0.303383   -21.00000   -1.50000
H132  2    0.194195    0.626606    0.130621   -21.00000   -1.50000
H133  2    0.407982    0.771230    0.112116   -21.00000   -1.50000
"""

ins_disordered_with_part_sof = """\
TITL 02srv053 SADABS in P2(1)/n
CELL 0.71073  12.823  13.422  18.550  90.000 103.91  90.000
ZERR    4.00   0.003   0.003   0.004   0.000   0.01   0.000
LATT  1
SYMM 0.5-X, 0.5+Y, 0.5-Z
SFAC C H S CL
REM ................ skipping stuff I don't test here ..............
FVAR       0.31540   0.05337
FVAR       0.03075   0.89064
S1    3    0.366783    0.363927    0.347492    11.00000    0.02050    0.01721 =
         0.02652    0.00339    0.00979   -0.00170
S2    3    0.151555    0.283625    0.288656    11.00000    0.02118    0.01442 =
         0.01699   -0.00030    0.00306    0.00146
S3    3    0.294777    0.579144    0.332456    11.00000    0.03635    0.01521 =
         0.02250    0.00137    0.01462   -0.00354
REM ............. skipping many atoms .........................
CL1   4   -0.126497    0.262053    0.433117    11.00000    0.02755    0.03548 =
         0.05670    0.00439    0.01834    0.00260
PART 1 41
CL2   4    0.096994    0.328954    0.480493    31.00000    0.03167    0.04296 =
         0.03718    0.00702    0.00390   -0.00524
C28   1    0.008586    0.236935    0.433750    31.00000    0.02306    0.02292 =
         0.04054    0.00599    0.01598    0.00352
AFIX   3
H28A  2    0.029496    0.171785    0.457990    21.00000    0.03492
H28B  2    0.015946    0.232245    0.381920    21.00000    0.03736
AFIX   0
PART 2 -41
CL2B  4    0.068342    0.365253    0.494288   -31.00000    0.07647
C28B  1    0.028119    0.260724    0.451358   -31.00000    0.08827
PART 0
HKLF 4
"""

ins_missing_sfac = """\
CELL 0 1 2 3 90 90 90
O 4 0.1 0.2 0.3 11 0.04
"""

ins_invalid_scatt = """\
CELL 0 1 2 3 90 90 90
SFAC C O
O 4 0.1 0.2 0.3. 11 0.04
"""

ins_invalid_scatt_1 = """\
CELL 0 1 2 3 90 90 90
SFAC C O
O 4 0.1 0.2 0.3 11
"""

ins_dfix_across_symm="""\
TITL s031 in C2/c #15
CELL 0.71073 28.3148 15.8506 18.226 90 118.092 90
ZERR 4 0.0006 0.0003 0.0004 0 0.001 0
LATT 7
SYMM -X,+Y,0.5-Z
SFAC C H Cl
UNIT 282 382 6
FVAR 0.1818
WGHT 0.1
EQIV $1 1-X,+Y,0.5-Z
DFIX 1.75 0.02 CL1_$1 C0AA
DFIX 1.75 0.03 C0AA CL3_$1
DFIX 1.75 CL2 C0AA
CL1   3     0.55003  0.45683  0.24917  10.25000  0.05233
CL2   3     0.54956  0.33497  0.36400  10.25000  0.05483
CL3   3     0.53790  0.33648  0.30419  10.25000  0.08474
C0AA  1     0.50786  0.40217  0.27444  10.25000  0.03965
HKLF 4
"""

ins_sadi = """\
CELL 0 1 2 3 90 90 90
SFAC C
SADI 0.04 C1 C3 C2 C4
C1    1     0.50000  0.24000  0.25000  11.00000  0.05233
C2    1     0.50000  0.26000  0.25000  11.00000  0.05233
C3    1     0.00000  0.25000  0.25000  11.00000  0.05233
C4    1     0.00000  0.25000  0.25000  11.00000  0.05233
"""

ins_sadi_with_sym = """\
CELL 0 1 2 3 90 90 90
SYMM -X,+Y,0.5-Z
EQIV $1 1-X,+Y,0.5-Z
SFAC C
SADI C1_$1 C3 C2 C4_$1
C1    1     0.50000  0.24000  0.25000  11.00000  0.05233
C2    1     0.50000  0.26000  0.25000  11.00000  0.05233
C3    1     0.00000  0.25000  0.25000  11.00000  0.05233
C4    1     0.00000  0.25000  0.25000  11.00000  0.05233
"""

ins_flat="""\
TITL fred in P2(1)/c
CELL 0.71073 15.149 11.4915 16.2458 90 99.462 90
ZERR 3 0.0054 0.0043 0.006 0 0.007 0
LATT 1
SYMM -X,0.5+Y,0.5-Z
SFAC C H
UNIT 134 142
FMAP 2.0
WGHT 0.1
FLAT 0.1 C31 C32 C33 C34 C35 C36 C37
FLAT C41 C42 C43 C44 C45 C46 C47
FVAR 0.16431
PART 1
C31   1     0.28967  0.31613 -0.02773  10.75000  0.05888  0.04514  0.10878 =
 -0.01962 -0.00471 -0.00143
C32   1     0.30768  0.41649  0.03308  10.75000  0.03521  0.04076  0.05711 =
 0.01504  0.01485  0.01552
C33   1     0.28897  0.41092  0.11289  10.75000  0.06042  0.09876  0.06719 =
 0.05095  0.03457  0.04942
C34   1     0.30734  0.50229  0.16767  10.75000  0.08796  0.14534  0.06011 =
 0.03723  0.02945  0.07980
C35   1     0.34782  0.60688  0.14285  10.75000  0.09111  0.10467  0.04177 =
 -0.01494  0.00086  0.05687
C36   1     0.36857  0.61562  0.06124  10.75000  0.08041  0.06114  0.05218 =
 0.00249  0.01982  0.01339
C37   1     0.34923  0.52118  0.01013  10.75000  0.04695  0.03550  0.04216 =
 -0.00583  0.00719 -0.00031
PART 2
C41   1     0.40898  0.68460  0.08429  10.25000  0.08117  0.05294  0.24136 =
 -0.04489  0.04281 -0.05242
C42   1     0.35781  0.57012  0.08520  10.25000  0.02162  0.06254  0.03143 =
 -0.00552  0.00755 -0.01371
C43   1     0.34017  0.49729  0.00895  10.25000  0.04695  0.03550  0.04216 =
 -0.00583  0.00719 -0.00031
C44   1     0.30214  0.37738  0.00610  10.25000  0.04686  0.18231  0.09967 =
 0.09419  0.03785  0.07323
C45   1     0.27044  0.34405  0.07795  10.25000  0.03487  0.04084  0.25569 =
 -0.01301 -0.02616  0.01005
C46   1     0.28348  0.40879  0.14914  10.25000  0.02746  0.08814  0.07379 =
 0.00346 -0.00267  0.02126
C47   1     0.32637  0.52368  0.16222  10.25000  0.03065  0.05189  0.02130 =
 -0.00058  0.00299 -0.00724
HKLF 4
"""

def template_ins(instruction):
  return """\
CELL 0 5 10 15 90 90 90
SFAC C
%s
C1    1     0.28967  0.31613 -0.02773  10.75000  0.05888  0.04514  0.10878 =
 -0.01962 -0.00471 -0.00143
C2    1     0.30768  0.41649  0.03308  10.75000  0.03521  0.04076  0.05711 =
 0.01504  0.01485  0.01552
C3    1     0.28897  0.41092  0.11289  10.75000  0.06042  0.09876  0.06719 =
 0.05095  0.03457  0.04942
C4    1     0.30734  0.50229  0.16767  10.75000  0.08796  0.14534  0.06011 =
 0.03723  0.02945  0.07980
""" %instruction

ins_simu_simple = template_ins("SIMU")
ins_simu_s = template_ins("SIMU 0.03")
ins_simu_s_st = template_ins("SIMU 0.03 0.07")
ins_simu_s_st_dmax = template_ins("SIMU 0.03 0.07 1.8")
ins_simu_atoms = template_ins("SIMU C1 C2 C3")
ins_simu_s_atoms = template_ins("SIMU 0.03 C1 C2 C3")
ins_simu_s_st_atoms = template_ins(
  "SIMU 0.03 0.07 C1 C2 C3")
ins_simu_s_st_dmax_atoms = template_ins(
  "SIMU 0.03 0.07 1.7 C1 C2 C3")

ins_delu_simple = template_ins("DELU")
ins_delu_s1 = template_ins("DELU 0.02")
ins_delu_s1_s2 = template_ins("DELU 0.01 0.02")
ins_delu_atoms = template_ins("DELU C1 C2 C3")
ins_delu_s1_atoms = template_ins("DELU 0.02 C1 C2 C3")
ins_delu_s1_s2_atoms = template_ins(
  "DELU 0.01 0.02 C1 C2 C3")

ins_rigu_simple = template_ins("RIGU")
ins_rigu_s1 = template_ins("RIGU 0.02")
ins_rigu_s1_s2 = template_ins("RIGU 0.01 0.02")
ins_rigu_atoms = template_ins("RIGU C1 C2 C3")
ins_rigu_s1_atoms = template_ins("RIGU 0.02 C1 C2 C3")
ins_rigu_s1_s2_atoms = template_ins(
  "RIGU 0.01 0.02 C1 C2 C3")

ins_isor_simple = template_ins("ISOR")
ins_isor_s = template_ins("ISOR 0.2")
ins_isor_s_st = template_ins("ISOR 0.2 0.3")
ins_isor_atoms = template_ins("ISOR C1 C2 C3")
ins_isor_s_atoms = template_ins("ISOR 0.2 C1 C2 C3")
ins_isor_s_st_atoms = template_ins(
  "ISOR 0.2 0.3 C1 C2 C3")

ins_dfix_multiple = template_ins("DFIX 1.7 C1 C2 C3 C4")
ins_invalid_dfix = template_ins("DFIX C1 C2")
ins_invalid_dfix_2 = template_ins("DFIX 1.7 C1 C2 C3")
ins_invalid_flat = template_ins("FLAT C1 C2 C3")

ins_equal_sign_in_rem = """\
REM Solution 1  R1  0.100,  Alpha = 0.0015  in P2(1)
REM C13 O10
TITL SUCROSE IN P2(1)
"""

ins_with_atom_peak_heights = """\
TITL SUCROSE IN P2(1)
CELL  0.71073   7.7830   8.7364  10.9002   90.000  102.984   90.000
ZERR 2 0.001 0.0012 0.0015 0 0.009 0
LATT -1
SYMM -X, 1/2+Y, -Z
SFAC C H O
UNIT 24 44 22
L.S. 10
BOND
LIST 6
FMAP 2
PLAN 20
ANIS
DELU
O001  3  0.60914  0.62292  0.82801 11.00000  0.01991   18.33
O002  3  0.63134  0.57137  0.62239 11.00000  0.02191   8.21
O003  3  0.68506  0.87525  0.78808 11.00000  0.02117   8.01
O004  3  0.25196  0.53353  0.76983 11.00000  0.02850   7.89
O005  3  0.37848  0.73377  0.96996 11.00000  0.02967   7.87
O006  3  0.79443  0.65096  1.07392 11.00000  0.02911   7.82
O007  3  1.08974  0.87253  1.02226 11.00000  0.02929   7.77
O008  3  0.96030  0.73231  0.67382 11.00000  0.03335   7.69
O009  3  0.29687  0.22403  0.69006 11.00000  0.03867   7.41
O00A  3  0.71468  0.42497  0.41792 11.00000  0.03487   7.28
C00B  1  0.64414  0.15501  0.65110 11.00000  0.02137   6.85
C00C  1  0.51386  0.61153  0.69938 11.00000  0.01950   6.36
C00D  1  0.78440  0.77821  0.99237 11.00000  0.02002   6.06
C00E  1  0.36586  0.49695  0.68863 11.00000  0.01963   6.03
C00F  1  0.94543  0.80269  0.93471 11.00000  0.02056   5.97
C00G  1  0.63136  0.77807  0.87607 11.00000  0.02041   5.94
C00H  1  0.55937  0.29776  0.62561 11.00000  0.02281   5.87
C00I  1  0.45580  0.83872  0.89621 11.00000  0.02767   5.77
C00J  1  0.87184  0.90817  0.82300 11.00000  0.02382   5.76
C00K  1  0.43548  0.33306  0.71436 11.00000  0.02549   5.76
C00L  1  0.81534  0.39973  0.54440 11.00000  0.02952   5.73
C00M  1  0.70459  0.41948  0.64084 11.00000  0.02597   5.72
C00N  1  0.95213  0.88631  0.71108 11.00000  0.02971   5.35
HKLF 4
END
"""

ins_with_q_peaks = """\
TITL 01mrh009 in P2(1)/c
CELL 0.71073   8.7925   8.2663   8.5561  90.000 111.141  90.000
ZERR    2.00   0.0018   0.0017   0.0017   0.000   0.030   0.000
LATT 1
SYMM -X, 0.5+Y, 0.5-Z
SFAC C  H  N  O  P  V
UNIT 8 24 4 20 4 4
L.S. 30
ACTA
BOND
FMAP 2
PLAN 5
HTAB
WGHT  0.032500  0.937600
FVAR  0.671230
MOLE 1
V1     6   0.657682  1.034967  0.322750 11.000000  0.005490  0.006760 =
           0.007550 -0.000030  0.002830 -0.000560
P1     5   0.557872  0.733904  0.520843 11.000000  0.006330  0.006410 =
           0.006990  0.000100  0.002880 -0.000130
O2     4   0.645870  0.883988  0.489582 11.000000  0.010550  0.008480 =
           0.009460  0.000070  0.004270 -0.002350
Q1     1   0.615600  0.695100  0.623900 11.000000  0.050000      0.64
Q2     1   0.526800  0.678600  0.431500 11.000000  0.050000      0.63
Q3     1   0.620800  1.025100  0.211700 11.000000  0.050000      0.46
HKLF 4
END
"""

ins_residues = """\
TITL 3odl
CELL  1.00000   59.140   59.140  186.740  90.00  90.00 120.00  ! lambda & Cell
ZERR        6    0.059    0.059    0.187   0.00   0.00   0.00  !Z & cell esds

LATT -1     !P, non centrosymmetric
SYMM -Y, X-Y, 2/3+Z
SYMM -X+Y, -X, 1/3+Z
SYMM Y, X, -Z
SYMM X-Y, -Y, 1/3-Z
SYMM -X, -X+Y, 2/3-Z
SFAC  C    H     N    O    S  P   NA   !scattering factor types and
UNIT  9264.63 13529.27 2472 6437.14 76.02 12   6   !unit cell contents

DELU $C_* $N_* $O_* $S_* $P_*     ! Rigid bond restraints -ignored for iso
SIMU 0.1 $C_* $N_* $O_* $S_* $P_* ! Similar U restraints - iso or anis.
ISOR 0.1 O_1001 > LAST            ! Approximate iso restraints for waters;

FLAT_FAD    AN1 AC2 AN3 AC4 AC5 AC6 AN6 AN7 AC8 AN9 AC1*
DFIX_FAD   1.362   AN1  AC6
DFIX_FAD   1.368   AN1  AC2
DANG_FAD 2.240 0.10 AN9  AN7
DANG_FAD 2.196 0.10 AN9  AC5
DANG_8 2.425 0.026 CA N_+
DANG_9 2.250 0.017 O_- N
DFIX_* 1.329 C_- N
DANG_* 2.425 CA_- N
DANG_* 2.250 O_- N
DANG_* 2.435 C_- CA
FLAT_8 0.5 O CA N_+ C CA_+

RESI    8   SER
N     3   -0.575274    0.141191    0.113592    11.00000    0.33808    0.75677 =
         1.11905   -0.27101    0.19945   -0.04973
CA    1   -0.551908    0.139429    0.115198    11.00000    0.48936    0.25522 =
         0.84266   -0.05734    0.13568   -0.12284
C     1   -0.525230    0.163359    0.114621    11.00000    0.36302    0.24661 =
         0.62351    0.04491    0.09919   -0.02604
O     4   -0.505790    0.161790    0.116768    11.00000    0.52134    0.27541 =
         0.56288    0.00938    0.07008    0.11741

RESI    9   LYS
N     3   -0.520798    0.186053    0.111792    11.00000    0.25817    0.26426 =
         0.24642   -0.03698    0.01196   -0.03745
CA    1   -0.494101    0.207878    0.111725    11.00000    0.16273    0.21215 =
         0.21593   -0.04270    0.01910    0.02736
C     1   -0.487071    0.220051    0.119169    11.00000    0.13133    0.19993 =
         0.20434   -0.03059    0.03663    0.03709
O     4   -0.503787    0.219454    0.123261    11.00000    0.14634    0.36185 =
         0.24067   -0.07422    0.04195    0.06255

RESI  415   FAD

AC1*  1    0.301769    0.460821    0.091581    11.00000    0.17934    0.39185 =
         0.16499   -0.02542   -0.00370    0.18623
AN9   3    0.329069    0.474690    0.094290    11.00000    0.16649    0.37106 =
         0.19437    0.02459    0.01637    0.16403
AC8   1    0.338965    0.493725    0.099532    11.00000    0.19389    0.42740 =
         0.20073   -0.01332   -0.03583    0.17001
AN7   3    0.363614    0.500512    0.100736    11.00000    0.18827    0.48315 =
         0.29375    0.00036   -0.03961    0.17524
AC5   1    0.369875    0.488232    0.095578    11.00000    0.20021    0.53248 =
         0.30473    0.04033   -0.01535    0.22652
AC6   1    0.394249    0.488600    0.094022    11.00000    0.20036    0.55609 =
         0.36639    0.15609    0.05588    0.23824
AN6   3    0.416744    0.505330    0.096997    11.00000    0.20718    0.56035 =
         0.55246    0.16445    0.00471    0.22417
AN1   3    0.392969    0.470604    0.089321    11.00000    0.25898    0.60811 =
         0.36195    0.17095    0.09827    0.30193
AC2   1    0.370280    0.453766    0.086063    11.00000    0.28075    0.52743 =
         0.37940    0.12160    0.10576    0.31394
AN3   3    0.347491    0.452243    0.087106    11.00000    0.24288    0.43252 =
         0.31698    0.06838    0.08486    0.24807
AC4   1    0.348480    0.469070    0.091969    11.00000    0.20789    0.48822 =
         0.24296    0.02751    0.01321    0.23800


RESI 1001   HOH
O     4   -0.446670    0.031582    0.038635    11.00000    0.16482    0.14614 =
         0.17802   -0.00150   -0.02394    0.01386

RESI 1002   HOH
O     4   -0.368462    0.037104    0.034398    11.00000    0.17923    0.12876 =
         0.18971   -0.03786   -0.01430    0.04562

RESI 1003   HOH
O     4   -0.422007    0.007732    0.032882    11.00000    0.19693    0.17138 =
         0.21503   -0.04115    0.01271    0.01223

RESI 1004   HOH
O     4   -0.350514    0.045091    0.020274    11.00000    0.17129    0.15651 =
         0.19856   -0.03983    0.01063    0.03671


"""

# This is one of the test structures shipped with Olex 2
# This version has the correct modeling of the disorder C7a/C7b
# plus a spurious split N3/C3 with occupancies adding up to 1.
ins_thpp = """\
TITL
CELL 0.71073 6.9196 14.5749 9.7248 90 90.637 90
ZERR 4 0.0001 0.0002 0.0001 0 0.001 0
LATT 1
SYMM 0.5-X,0.5+Y,0.5-Z

SFAC C H F N
UNIT 40 40 8 16
EXYZ N3 C3


L.S. 4
PLAN  20
more -3
fmap 2
REM <HKL>thpp.hkl</HKL>
WGHT 0.1
FVAR 0.35838 0.87977 0.5

F1    3     0.16726  0.42638 -0.23772  11.00000  0.03596  0.02944  0.01904 =
 -0.00761 -0.00298 -0.00009
F2    3     0.13095  0.31497 -0.01246  11.00000  0.04052  0.01691  0.03281 =
 -0.00250 -0.00502 -0.00279
N8    4     0.36220  0.66937  0.07644  11.00000  0.02801  0.01729  0.01820 =
 -0.00038  0.00220 -0.00146
N3    4     0.22193  0.43032  0.12983  31.00000  0.02131
C9    1     0.30101  0.58265  0.04215  11.00000  0.01747  0.01804  0.01766 =
 -0.00073  0.00094  0.00221
C4    1     0.27507  0.51672  0.15265  11.00000  0.02018  0.01935  0.01696 =
 0.00027  0.00016  0.00159
N5    4     0.30293  0.54252  0.28533  11.00000  0.03661  0.02191  0.01588 =
 0.00077 -0.00090 -0.00391
C2    1     0.18838  0.40270  0.00228  11.00000  0.02311  0.01671  0.02386 =
 -0.00235 -0.00203  0.00103
C10   1     0.26210  0.54987 -0.09061  11.00000  0.01942  0.02080  0.01702 =
 -0.00001  0.00153  0.00268
C1    1     0.20604  0.45735 -0.10937  11.00000  0.02242  0.02167  0.01778 =
 -0.00389 -0.00127  0.00281
C11   1     0.26071  0.60304 -0.21444  11.00000  0.03058  0.02565  0.01889 =
 -0.00137  0.00044  0.00016
C13   1     0.26803  0.47611  0.39401  11.00000  0.03989  0.02833  0.01867 =
 0.00388 -0.00072 -0.00350
C6    1     0.36454  0.63430  0.32066  11.00000  0.04896  0.02569  0.02000 =
 -0.00456 -0.00017 -0.00449
N12   4     0.25129  0.63888 -0.31966  11.00000  0.05489  0.03497  0.02311 =
 0.00384 -0.00181 -0.00450
PART 1
C7a   1     0.30537  0.70214  0.21212  21.00000  0.03415  0.02081  0.02014 =
 -0.00399  0.00271  0.00061
PART 0
C14   1     0.39584  0.73936 -0.02671  11.00000  0.03778  0.02073  0.02508 =
 0.00407  0.00124 -0.00460
PART 2
C7b   1     0.40368  0.69420  0.21920 -21.00000  0.02814  0.02120  0.02427 =
 -0.00252 -0.00560  0.00546
PART 0
C3    1     0.22193  0.43032  0.12983 -31.00000  0.02131
HKLF 4
END

"""


if __name__ == '__main__':
  run()
