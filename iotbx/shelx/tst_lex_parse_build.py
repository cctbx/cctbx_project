from __future__ import division
from cctbx import crystal
from cctbx import uctbx
from cctbx import sgtbx
from cctbx import adptbx
from cctbx import xray
from iotbx import shelx
from iotbx.shelx import crystal_symmetry_from_ins
import iotbx.builders
from libtbx.test_utils import approx_equal, Exception_expected
from libtbx.math_utils import are_equivalent
import cStringIO

def exercise_lexing():
  stream = shelx.command_stream(file=cStringIO.StringIO(ins_mundane_tiny))
  i = iter(stream)
  try:
    cmd, line = i.next()
    assert cmd, line == ('TITL', ('in Pbca',))
    cmd, line =  i.next()
    assert cmd, line == ('CELL', (0.71073, 7.35, 9.541, 12.842, 90, 90, 90))
    cmd, line = i.next()
    assert cmd, line == ('ZERR', (4, 0.002, 0.002, 0.003, 0, 0, 0))
    cmd, line = i.next()
    assert cmd, line == ('LATT', (1,))
    cmd, line = i.next()
    assert cmd, line == ('SYMM', ('0.5-X, -Y, 0.5+Z',))
    cmd, line = i.next()
    assert cmd, line == ('SYMM', ('-X, 0.5+Y, 0.5-Z',))
    cmd, line = i.next()
    assert cmd, line == ('SYMM', ('1/2+X, 0.5-Y, -Z',))
    cmd, line =  i.next()
    assert cmd, line == ('SFAC', ('C', 'H', 'O', 'N',))
    cmd, line = i.next() # UNIT
    cmd, line = i.next() # TEMP
    cmd, line = i.next()
    assert cmd, line == ('L.S.', (4,))
    cmd, line = i.next()
    assert cmd, line == ('BOND', ((stream.element_tok, 'H'),))
    cmd, line = i.next() # FMAP
    cmd, line = i.next() # PLAN
    cmd, line = i.next() # WGHT
    cmd, line = i.next() # EXTI
    cmd, line = i.next() # FVAR
    cmd, line = i.next()
    assert cmd, line == ('REM', ('Protracted example of residues on command',))
    cmd, line = i.next()
    assert cmd, line == ('HFIX', (stream.residue_number_tok, 1), (23,))
    cmd, line =  i.next()
    assert cmd, line == ('HFIX', (stream.residue_class_tok, 'N'), (43,))
    cmd, line = i.next()
    assert cmd, line == ('EQIV', (1, '1-X, -Y, -Z'))
    cmd, line = i.next()
    assert cmd, line == ('CONF', ( (stream.atom_tok, 'C4', None),
                             (stream.atom_tok, 'N', None),
                             (stream.atom_tok, 'H', None),
                             (stream.atom_tok, 'O2', 1) ) )
    cmd, line = i.next()
    assert cmd, line == ('__ATOM__',
                   ('O2', 3, 0.362893, 0.160589, -0.035913, 11,
                          0.03926, 0.02517, 0.02140,
                          -0.00415, -0.00810, 0.01009))
    cmd, line = i.next()
    assert cmd, line == ('__ATOM__',
                   ('O3', 3, 0.696722, 0.119176, 0.260657, 11,
                          0.02838, 0.02133, 0.02918,
                          0.00011, -0.01030, -0.00048))
    cmd, line = i.next() # C1
    cmd, line = i.next() # C4
    cmd, line =  i.next()
    assert cmd, line == ('RESI', (1,))
    cmd, line = i.next() # C2
    cmd, line = i.next() # C3
    cmd, line =  i.next()
    assert cmd, line == ('RESI', ('N',))
    cmd, line = i.next() # N
    cmd, line = i.next() # HKLF
    try:
      cmd, line = i.next()
      raise AssertionError
    except StopIteration:
      pass
  except StopIteration:
    raise AssertionError

def exercise_crystal_symmetry_parsing():
  stream = shelx.command_stream(file=cStringIO.StringIO(ins_mundane_tiny))
  l = shelx.crystal_symmetry_parser(
    stream,
    builder=iotbx.builders.crystal_symmetry_builder())
  l.parse()
  assert l.builder.crystal_symmetry.is_similar_symmetry(
    crystal.symmetry(
      unit_cell=uctbx.unit_cell((7.350, 9.541, 12.842, 90, 90, 90)),
      space_group_symbol='Pbca'),
    relative_length_tolerance=1e-15,
    absolute_angle_tolerance=1e-15)
  cs = crystal_symmetry_from_ins.extract_from(
    file=cStringIO.StringIO(ins_mundane_tiny))
  assert cs.is_similar_symmetry(l.builder.crystal_symmetry)

  stream = shelx.command_stream(file=cStringIO.StringIO(ins_P1))
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
    file=cStringIO.StringIO(ins_P1))
  assert cs.is_similar_symmetry(l.builder.crystal_symmetry)

def exercise_instruction_parsing():
  try:
    _ = iotbx.builders.weighting_scheme_builder
    alternatives = (None, _())
  except AttributeError:
    alternatives = (None, )
  for builder in alternatives:
    stream = shelx.command_stream(file=cStringIO.StringIO(ins_aspirin))
    l = shelx.instruction_parser(stream, builder)
    l.parse()
    ins = l.instructions
    assert ins['hklf']['s'] == 1
    assert ins['hklf']['matrix'].as_xyz() == 'x,y,z'
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
    if builder is not None:
      assert builder.temperature_in_celsius == ins['temp']
      ws = builder.weighting_scheme
      assert isinstance(
        ws,
        iotbx.builders.least_squares.mainstream_shelx_weighting)
      assert ws.a == ins['wght']['a']
      assert ws.b == ins['wght']['b']

def exercise_xray_structure_parsing():
  exercise_special_positions()
  exercise_aspirin()
  exercise_disordered()
  exercise_invalid()

def exercise_special_positions():
  structure = xray.structure.from_shelx(
    file=cStringIO.StringIO(ins_special_positions))
  occupancies = [ sc.occupancy for sc in structure.scatterers() ]
  multiplicities = [ sc.multiplicity() for sc in structure.scatterers() ]
  assert multiplicities == [ 2, 2, 2, 2, 6 ]
  assert approx_equal(occupancies, [ 1, 0.9999, 1, 1, 1 ], eps=1e-15)

def exercise_aspirin():
  for set_grad_flags in (False, True):
    structure = xray.structure.from_shelx(
      file=cStringIO.StringIO(ins_aspirin),
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
    stream = shelx.command_stream(file=cStringIO.StringIO(ins_disordered))
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
        file=cStringIO.StringIO(ins_invalid_scatt),
        set_grad_flags=set_grad_flags)
      raise Exception_expected
    except RuntimeError, e:
      assert str(e) == "ShelX: illegal argument '0.3.' at line 3"

    try:
      structure = xray.structure.from_shelx(
        file=cStringIO.StringIO(ins_invalid_scatt_1),
        set_grad_flags=set_grad_flags)
      raise Exception_expected
    except RuntimeError, e:
      assert str(e) == ("ShelX: wrong number of parameters "
                        "for scatterer at line 3")

    try:
      structure = xray.structure.from_shelx(
        file=cStringIO.StringIO(ins_missing_sfac),
        set_grad_flags=set_grad_flags)
      raise Exception_expected
    except RuntimeError, e:
      assert e.args[0].startswith('ShelX:')

  structure = xray.structure.from_shelx(
    file=cStringIO.StringIO(ins_disordered_with_part_sof))
  occ = 0.89064
  for sc in structure.scatterers():
    if sc.label in ('CL2', 'C28', 'H28A', 'H28B'):
      assert approx_equal(sc.occupancy, occ)
    elif sc.label in ('CL2B', 'C28B'):
      assert approx_equal(sc.occupancy, 1-occ)

def exercise_afix_parsing():
  import smtbx.refinement.constraints.geometrical.hydrogens as _
  try:
    builder = iotbx.builders.constrained_crystal_structure_builder()
  except AttributeError:
    print "Skipping AFIX parsing: : smtbx module not available"
    return
  stream = shelx.command_stream(file=cStringIO.StringIO(ins_aspirin))
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
  try:
    builder = iotbx.builders.constrained_crystal_structure_builder()
  except AttributeError:
    print "Skipping u_iso = k u_eq test: smtbx module not available"
    return
  stream = shelx.command_stream(file=cStringIO.StringIO(ins_aspirin))
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
  import libtbx.load_env
  if (not libtbx.env.has_module(name="smtbx")):
    print "Skipping exercise_restraint_parsing():" \
      " smtbx module is not available."
    return
  import smtbx.refinement.restraints
  def parse_restraints(ins_name):
    builder = iotbx.builders.restrained_crystal_structure_builder()
    stream = shelx.command_stream(file=cStringIO.StringIO(ins_name))
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
  except RuntimeError, e:
    assert str(e) == "ShelX: Invalid DFIX instruction at line 3"
  try:
    proxies = parse_restraints(ins_invalid_dfix_2)
  except RuntimeError, e:
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
  except RuntimeError, e:
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
  # ISOR simple
  proxies = parse_restraints(ins_isor_simple)
  shared_isotropic_adp_proxy \
      = proxies['isotropic_adp']
  assert shared_isotropic_adp_proxy.size() == 4
  expected_i_seqs = (0,1,2,3)
  expected_weights = (25,100,100,25)
  for i in range(4):
    assert shared_isotropic_adp_proxy[i].i_seq == expected_i_seqs[i]
    assert approx_equal(shared_isotropic_adp_proxy[i].weight,
                        expected_weights[i])
  # ISOR with s
  proxies = parse_restraints(ins_isor_s)
  shared_isotropic_adp_proxy \
      = proxies['isotropic_adp']
  assert shared_isotropic_adp_proxy.size() == 4
  expected_weights = (6.25,25,25,6.25)
  for i in range(4):
    assert shared_isotropic_adp_proxy[i].i_seq == expected_i_seqs[i]
    assert approx_equal(shared_isotropic_adp_proxy[i].weight,
                        expected_weights[i])
  # ISOR with s and st
  proxies = parse_restraints(ins_isor_s_st)
  shared_isotropic_adp_proxy \
      = proxies['isotropic_adp']
  assert shared_isotropic_adp_proxy.size() == 4
  expected_weights = (1/(0.3*0.3), 25, 25, 1/(0.3*0.3))
  for i in range(4):
    assert shared_isotropic_adp_proxy[i].i_seq == expected_i_seqs[i]
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
    assert shared_isotropic_adp_proxy[i].i_seq == expected_i_seqs[i]
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
    assert shared_isotropic_adp_proxy[i].i_seq == expected_i_seqs[i]
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
    assert shared_isotropic_adp_proxy[i].i_seq == expected_i_seqs[i]
    assert approx_equal(shared_isotropic_adp_proxy[i].weight,
                        expected_weights[i])

def shelx_u_cif(unit_cell, u_star):
  u_cif = adptbx.u_star_as_u_cif(unit_cell, u_star)
  u_cif = u_cif[0:3] + (u_cif[-1], u_cif[-2], u_cif[-3])
  return (" "*3).join([ "%.5f" % x for x in  u_cif ])

def run():
  exercise_instruction_parsing()
  exercise_restraint_parsing()
  exercise_u_iso_proportional_to_u_eq_parsing()
  exercise_afix_parsing()
  exercise_xray_structure_parsing()
  exercise_crystal_symmetry_parsing()
  exercise_lexing()
  print 'OK'

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
"REM  Protracted example of residues on command\n"
"HFIX_1 23\n"
"HFIX_N 43\n"
"EQIV $1 1-X, -Y, -Z\n"
"CONF C4 N H O2_$1\n"
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
HKLF 4 1 1 0 0 0 1 0 0 0 1

REM  aspirin in P2(1)/c
REM R1 =  0.0455 for   1038 Fo > 4sig(Fo)  and  0.0990 for all   1806 data
REM    110 parameters refined using      0 restraints

END

WGHT      0.0534      0.3514
"""

ins_disordered = """
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

ins_disordered_with_part_sof = """
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

ins_missing_sfac = """
CELL 0 1 2 3 90 90 90
O 4 0.1 0.2 0.3 11 0.04
"""

ins_invalid_scatt = """
CELL 0 1 2 3 90 90 90
SFAC C O
O 4 0.1 0.2 0.3. 11 0.04
"""

ins_invalid_scatt_1 = """
CELL 0 1 2 3 90 90 90
SFAC C O
O 4 0.1 0.2 0.3 11
"""

ins_dfix_across_symm="""
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

ins_sadi = """
CELL 0 1 2 3 90 90 90
SFAC C
SADI 0.04 C1 C3 C2 C4
C1    1     0.50000  0.24000  0.25000  11.00000  0.05233
C2    1     0.50000  0.26000  0.25000  11.00000  0.05233
C3    1     0.00000  0.25000  0.25000  11.00000  0.05233
C4    1     0.00000  0.25000  0.25000  11.00000  0.05233
"""

ins_sadi_with_sym = """
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

ins_flat="""
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
  return """
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

if __name__ == '__main__':
  run()
