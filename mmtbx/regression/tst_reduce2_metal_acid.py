"""reduce2 placement at acids bound to metals: an O-H whose O has one heavy
neighbour, an oxoacid centre (carboxylic and carbonic acids, and the
oxoacids of P and S with a P=O or S=O), and is bound to a metal (within
0.55 A of the summed covalent radii) is taken as deprotonated, and reduce2
does not build its H. A supplied H is kept and listed; with a metal of lower
occupancy than the O the H is kept and listed as unresolved. A phenol is not
an acid.
Models: 6IHT, phosphoserine X 3 and Mn 303 (O3P 2.25 A, O2P 2.40 A);
bicarbonate D 401 of 3WU2 with a metal placed 2.0 to 2.6 A from O3, on the
C-O3 line, and a copy of it 10 A away without a metal; a Tyr with a Zn
placed 2.0 A from OH, on the CZ-OH line.
"""
from libtbx.utils import null_out
from scitbx import matrix
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2

pdb_6iht = """\
CRYST1   46.693   38.609   65.003  90.00 101.89  90.00 P 1 21 1
HETATM 1951  N   SEP X   3     -17.641   2.560  81.946  1.00 43.47           N
HETATM 1952  CA  SEP X   3     -18.182   1.212  82.103  1.00 50.69           C
HETATM 1953  C   SEP X   3     -19.574   1.130  82.798  1.00 76.16           C
HETATM 1954  O   SEP X   3     -20.191   0.056  82.870  1.00 81.18           O
HETATM 1955  CB  SEP X   3     -17.159   0.429  82.964  1.00 55.47           C
HETATM 1956  OG  SEP X   3     -16.390  -0.567  82.295  1.00 27.38           O
HETATM 1957  OXT SEP X   3     -20.144   2.095  83.333  1.00 75.42           O
HETATM 1958  O1P SEP X   3     -15.197   0.450  80.303  1.00 17.16           O
HETATM 1959  O2P SEP X   3     -14.274  -1.538  81.367  1.00 16.63           O
HETATM 1960  O3P SEP X   3     -14.216   0.696  82.681  1.00 20.00           O
HETATM 1961  P   SEP X   3     -14.985  -0.197  81.677  0.31 24.79           P
HETATM 1964 MN    MN A 303     -12.440  -0.678  82.647  1.00 15.70          MN
END
"""

pdb_bct_zn = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
HETATM    1  C   BCT A   1     -25.131 -58.646  19.754  1.00 20.00           C
HETATM    2  O1  BCT A   1     -25.291 -59.027  18.508  1.00 20.00           O
HETATM    3  O2  BCT A   1     -24.215 -59.207  20.449  1.00 20.00           O
HETATM    4  O3  BCT A   1     -25.936 -57.757  20.169  1.00 20.00           O
HETATM    5 ZN    ZN A   2     -27.205 -56.356  20.823  1.00 20.00          ZN
HETATM    6  C   BCT A   3     -25.131 -58.646  29.754  1.00 20.00           C
HETATM    7  O1  BCT A   3     -25.291 -59.027  28.508  1.00 20.00           O
HETATM    8  O2  BCT A   3     -24.215 -59.207  30.449  1.00 20.00           O
HETATM    9  O3  BCT A   3     -25.936 -57.757  30.169  1.00 20.00           O
END
"""

pdb_tyr_zn = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
ATOM    631  N   TYR A  61      16.349   1.550   5.974  1.00  4.94           N
ATOM    632  CA  TYR A  61      16.423   2.984   6.193  1.00  4.82           C
ATOM    633  C   TYR A  61      16.454   3.256   7.702  1.00  4.90           C
ATOM    634  O   TYR A  61      16.013   2.471   8.506  1.00  6.74           O
ATOM    635  CB  TYR A  61      15.220   3.667   5.574  1.00  4.94           C
ATOM    636  CG  TYR A  61      15.121   3.593   4.077  1.00  4.69           C
ATOM    637  CD1 TYR A  61      16.129   4.014   3.237  1.00  5.43           C
ATOM    638  CD2 TYR A  61      13.935   3.177   3.486  1.00  5.20           C
ATOM    639  CE1 TYR A  61      15.979   3.988   1.850  1.00  5.86           C
ATOM    640  CE2 TYR A  61      13.755   3.183   2.105  1.00  5.81           C
ATOM    641  CZ  TYR A  61      14.779   3.580   1.287  1.00  5.43           C
ATOM    642  OH  TYR A  61      14.591   3.579  -0.066  1.00  7.78           O
HETATM  643 ZN    ZN A 201      14.316   3.578  -2.047  1.00 10.00          ZN
END
"""

def bct_metal_str(metal, distance):
  '''pdb_bct_zn with the metal at distance from O3 of BCT 1, on the C-O3
  line'''
  c = matrix.col((-25.131, -58.646, 19.754))
  o3 = matrix.col((-25.936, -57.757, 20.169))
  x, y, z = o3 + (o3 - c).normalize() * distance
  zn = ('HETATM    5 ZN    ZN A   2     -27.205 -56.356  20.823  1.00 20.00'
    '          ZN')
  line = (f'HETATM    5 {metal:<2}    {metal:>2} A   2    {x:8.3f}{y:8.3f}'
    f'{z:8.3f}  1.00 20.00          {metal:>2}')
  assert zn in pdb_bct_zn
  return pdb_bct_zn.replace(zn, line)

# pdb_6iht with an H on O3P
pdb_6iht_hop3 = pdb_6iht.replace('HETATM 1964 MN ', (
  'HETATM 1962 HOP3 SEP X   3     -14.478   1.651  82.542  1.00 20.00'
  '           H\nHETATM 1964 MN '))
assert pdb_6iht_hop3 != pdb_6iht

# pdb_6iht with a D on O3P
pdb_6iht_dop3 = pdb_6iht_hop3.replace(
  'HOP3 SEP X   3     -14.478   1.651  82.542  1.00 20.00           H',
  'DOP3 SEP X   3     -14.478   1.651  82.542  1.00 20.00           D')
assert pdb_6iht_dop3 != pdb_6iht_hop3

NOT_PLACED = ('Acid H not placed: the O is bound to a metal, taken as '
  'deprotonated')
UNRESOLVED = 'Acid H placed at metal sites with unresolved protonation'
SUPPLIED = ('Supplied acid H at metal sites: kept as supplied; acids bound '
  'to metals are deprotonated')

def run(prefix, pdb_str, args=()):
  '''reduce2 on a model string; returns the (residue number, name) of the H
  of its model and the text of the description file'''
  fn = f'{prefix}_in.pdb'
  txt = f'{prefix}.txt'
  with open(fn, 'w') as f:
    f.write(pdb_str)
  result = run_program(program_class=reduce2.Program, logger=null_out(),
    args=[fn, f'output.filename={prefix}_out.pdb',
          f'output.description_file_name={txt}', '--overwrite'] + list(args))
  with open(txt) as f:
    description = f.read()
  names = set()
  for a in result.model.get_hierarchy().atoms():
    if a.element_is_hydrogen():
      names.add((a.parent().parent().resseq_as_int(), a.name.strip()))
  return names, description

def test_6iht_phosphoserine():
  '''Both P-OH on Mn 303 (O3P 2.25 A, O2P 2.40 A) lose their H and are
  listed; with Mn 303 at occupancy 0.50 both keep their H and are listed as
  unresolved'''
  names, description = run('tst_reduce2_metal_acid_6iht', pdb_6iht)
  assert (3, 'HOP3') not in names and (3, 'HOP2') not in names, sorted(names)
  lines = description.splitlines()
  assert NOT_PLACED in lines, description
  i = lines.index(NOT_PLACED)
  assert sorted(lines[i + 1:i + 3]) == [
    '  HOP2 SEP X   3  (O2P SEP X   3 - MN    MN A 303 2.40 A (geometric))',
    '  HOP3 SEP X   3  (O3P SEP X   3 - MN    MN A 303 2.25 A (geometric))'
    ], description
  mn = ('HETATM 1964 MN    MN A 303     -12.440  -0.678  82.647  1.00 15.70'
    '          MN')
  assert mn in pdb_6iht
  names, description = run('tst_reduce2_metal_acid_6iht_half',
    pdb_6iht.replace(mn, mn.replace(' 1.00 15.70', ' 0.50 15.70')))
  assert (3, 'HOP3') in names and (3, 'HOP2') in names, sorted(names)
  lines = description.splitlines()
  assert NOT_PLACED not in lines, description
  assert UNRESOLVED in lines, description
  i = lines.index(UNRESOLVED)
  assert sorted(lines[i + 1:i + 3]) == [
    f'  HOP{k} SEP X   3: O{k}P occupancy 1.00 exceeds coordinating MN '
    'occupancy 0.50. Protonation may differ between bound and unbound '
    f'populations. HOP{k} was added using the default acid state; review '
    'this site.' for k in (2, 3)], description

def test_bicarbonate():
  '''Bicarbonate on Zn or Fe through O3, up to 0.55 A past the summed
  covalent radii (Zn-O 2.43 A), loses HO3 and is listed with the acids;
  farther away, or without a metal, it keeps HO3'''
  for k, (metal, distance) in enumerate((('ZN', 2.0), ('FE', 2.0),
                                         ('ZN', 2.35), ('ZN', 2.6))):
    names, description = run(f'tst_reduce2_metal_acid_bct_{k}',
      bct_metal_str(metal, distance))
    assert (3, 'HO3') in names, (metal, distance, sorted(names))
    lines = description.splitlines()
    if distance > 2.43:
      assert (1, 'HO3') in names, (metal, distance, sorted(names))
      assert NOT_PLACED not in lines, (metal, distance, description)
    else:
      assert (1, 'HO3') not in names, (metal, distance, sorted(names))
      assert NOT_PLACED in lines, (metal, distance, description)
      assert lines[lines.index(NOT_PLACED) + 1] == ('  HO3 BCT A   1  (O3  '
        f'BCT A   1 - {metal}    {metal} A   2 {distance:.2f} A (geometric))'
        ), (metal, distance, description)

def test_phenol_not_an_acid():
  '''A Tyr OH on Zn keeps HH: the ring C carries no second O, so it is not
  an acid centre'''
  names, description = run('tst_reduce2_metal_acid_tyr_zn', pdb_tyr_zn)
  assert (61, 'HH') in names, sorted(names)
  assert NOT_PLACED not in description.splitlines(), description

def test_supplied_h_kept():
  '''A supplied HOP3 on the O on Mn stays and is listed, not as not placed,
  with add keep_existing_H=True and with optimize; a supplied DOP3 written
  with neutron_option=all_d is listed under its own name'''
  contact = ('  HOP3 SEP X   3  (O3P SEP X   3 - MN    MN A 303 2.25 A '
    '(geometric)); no metal-coordination proxy')
  for k, args in enumerate((['approach=add', 'keep_existing_H=True'],
                            ['approach=optimize'])):
    names, description = run(f'tst_reduce2_metal_acid_supplied_{k}',
      pdb_6iht_hop3, args)
    assert (3, 'HOP3') in names, (args, sorted(names))
    lines = description.splitlines()
    assert SUPPLIED in lines, (args, description)
    assert lines[lines.index(SUPPLIED) + 1] == contact, (args, description)
    assert ('  HOP3 SEP X   3  (O3P SEP X   3 - MN    MN A 303 2.25 A '
      '(geometric))') not in lines, (args, description)
    names, description = run(f'tst_reduce2_metal_acid_supplied_d_{k}',
      pdb_6iht_dop3, args + ['neutron_option=all_d'])
    assert (3, 'DOP3') in names and (3, 'HOP3') not in names, (args,
      sorted(names))
    lines = description.splitlines()
    assert lines[lines.index(SUPPLIED) + 1] == contact.replace('HOP3',
      'DOP3'), (args, description)

if __name__ == '__main__':
  test_6iht_phosphoserine()
  test_bicarbonate()
  test_phenol_not_an_acid()
  test_supplied_h_kept()
  print('OK')
