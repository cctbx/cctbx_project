from cctbx import crystal
from cctbx import uctbx
from cctbx import sgtbx
from cctbx import adptbx
from cctbx import xray
from iotbx import shelx
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal, Exception_expected
from libtbx.math_utils import are_equivalent
import cStringIO

def exercise_lexing():
  stream = shelx.command_stream(file=cStringIO.StringIO(ins_mundane_tiny))
  i = iter(stream)
  try:
    cmd = i.next()
    assert cmd == ('TITL', ('in Pbca',))
    cmd =  i.next()
    assert cmd == ('CELL', (0.71073, 7.35, 9.541, 12.842, 90, 90, 90))
    cmd = i.next()
    assert cmd == ('ZERR', (4, 0.002, 0.002, 0.003, 0, 0, 0))
    cmd = i.next()
    assert cmd == ('LATT', (1,))
    cmd = i.next()
    assert cmd == ('SYMM', ('0.5-X, -Y, 0.5+Z',))
    cmd = i.next()
    assert cmd == ('SYMM', ('-X, 0.5+Y, 0.5-Z',))
    cmd = i.next()
    assert cmd == ('SYMM', ('1/2+X, 0.5-Y, -Z',))
    cmd =  i.next()
    assert cmd == ('SFAC', ('C', 'H', 'O', 'N',))
    cmd = i.next() # UNIT
    cmd = i.next() # TEMP
    cmd = i.next()
    assert cmd == ('L.S.', (4,))
    cmd = i.next()
    assert cmd == ('BOND', ((stream.element_tok, 'H'),))
    cmd = i.next() # FMAP
    cmd = i.next() # PLAN
    cmd = i.next() # WGHT
    cmd = i.next() # EXTI
    cmd = i.next() # FVAR
    cmd = i.next()
    assert cmd == ('REM ', ('Protracted example of residues on command',))
    cmd = i.next()
    assert cmd == ('HFIX', (stream.residue_number_tok, 1), (23,))
    cmd =  i.next()
    assert cmd == ('HFIX', (stream.residue_class_tok, 'N'), (43,))
    cmd = i.next()
    assert cmd == ('EQIV', (1, '1-X, -Y, -Z'))
    cmd = i.next()
    assert cmd == ('CONF', ( (stream.atom_tok, 'C4', None),
                             (stream.atom_tok, 'N', None),
                             (stream.atom_tok, 'H', None),
                             (stream.atom_tok, 'O2', 1) ) )
    cmd = i.next()
    assert cmd == ('O2', (3, 0.362893, 0.160589, -0.035913, 11,
                          0.03926, 0.02517, 0.02140,
                          -0.00415, -0.00810, 0.01009))
    cmd = i.next()
    assert cmd == ('O3', (3, 0.696722, 0.119176, 0.260657, 11,
                          0.02838, 0.02133, 0.02918,
                          0.00011, -0.01030, -0.00048))
    cmd = i.next() # C1
    cmd = i.next() # C4
    cmd =  i.next()
    assert cmd == ('RESI', (1,))
    cmd = i.next() # C2
    cmd = i.next() # C3
    cmd =  i.next()
    assert cmd == ('RESI', ('N',))
    cmd = i.next() # N
    cmd = i.next() # HKLF
    try:
      cmd = i.next()
      raise AssertionError
    except StopIteration:
      pass
  except StopIteration:
    raise AssertionError

def exercise_parsing():
  stream = shelx.command_stream(file=cStringIO.StringIO(ins_mundane_tiny))
  l = shelx.crystal_symmetry_parser(stream)
  l.parse()
  assert l.builder.crystal_symmetry.is_similar_symmetry(
    crystal.symmetry(
      unit_cell=uctbx.unit_cell((7.350, 9.541, 12.842, 90, 90, 90)),
      space_group_symbol='Pbca'),
    relative_length_tolerance=1e-15,
    absolute_angle_tolerance=1e-15)

  stream = shelx.command_stream(file=cStringIO.StringIO(ins_P1))
  l = shelx.crystal_symmetry_parser(stream)
  l.parse()
  assert l.builder.crystal_symmetry.is_similar_symmetry(
    crystal.symmetry(
      unit_cell=uctbx.unit_cell((1,2,3,99,100,101)),
      space_group_symbol='P1'),
    relative_length_tolerance=1e-15,
    absolute_angle_tolerance=1e-15)

  for set_grad_flags in (False, True):
    builder=shelx.crystal_structure_builder(set_grad_flags=set_grad_flags)
    stream = shelx.command_stream(file=cStringIO.StringIO(ins_aspirin))
    l_cs = shelx.crystal_symmetry_parser(stream, builder)
    l = shelx.atom_parser(l_cs.filtered_commands(), builder)
    l.parse()
    structure = l.builder.structure
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

  for set_grad_flags in (False, True):
    builder=shelx.crystal_structure_builder(set_grad_flags=set_grad_flags)
    stream = shelx.command_stream(file=cStringIO.StringIO(ins_disordered))
    l_cs = shelx.crystal_symmetry_parser(stream, builder)
    l = shelx.atom_parser(l_cs.filtered_commands(), builder)
    l.parse()
    structure = l.builder.structure
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

  builder=shelx.crystal_structure_builder()
  stream = shelx.command_stream(file=cStringIO.StringIO(ins_invalid_scatt))
  l_cs = shelx.crystal_symmetry_parser(stream, builder)
  l = shelx.atom_parser(l_cs.filtered_commands(), builder)
  try:
    l.parse()
    raise Exception_expected
  except shelx.illegal_argument_error, e:
    assert e.args[0] == 3 and e.args[-1] == '0.3.'

  builder=shelx.crystal_structure_builder()
  stream = shelx.command_stream(file=cStringIO.StringIO(ins_invalid_scatt_1))
  l_cs = shelx.crystal_symmetry_parser(stream, builder)
  l = shelx.atom_parser(l_cs.filtered_commands(), builder)
  try:
    l.parse()
    raise Exception_expected
  except shelx.illegal_scatterer_error, e:
    assert e.args[0] == 'O'

  builder=shelx.crystal_structure_builder()
  stream = shelx.command_stream(file=cStringIO.StringIO(ins_missing_sfac))
  l_cs = shelx.crystal_symmetry_parser(stream, builder)
  l = shelx.atom_parser(l_cs.filtered_commands(), builder)
  try:
    l.parse()
    raise Exception_expected
  except shelx.missing_sfac_error, e:
    pass

def shelx_u_cif(unit_cell, u_star):
  u_cif = adptbx.u_star_as_u_cif(unit_cell, u_star)
  u_cif = u_cif[0:3] + (u_cif[-1], u_cif[-2], u_cif[-3])
  return (" "*3).join([ "%.5f" % x for x in  u_cif ])

def run():
  exercise_lexing()
  exercise_parsing()
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
TEMP 20.000
L.S. 20
WGHT    0.068700    0.446300
FVAR       0.91641
O2    3    0.879181    0.140375    0.051044    11.00000    0.05893    0.05202 =
         0.05656    0.01670    0.01559    0.01156
AFIX 147
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
HKLF 4

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


if __name__ == '__main__':
  run()
