"""reduce2 placement at acids bound to metals: an O-H whose O has one heavy
neighbour, an oxoacid centre (carboxylic, phosphoric or phosphonic, sulfonic
or arsenic acid, bicarbonate), and is bound to a metal is taken as
deprotonated, and reduce2 does not build its H. A supplied H is kept and
listed.
Models: 6IHT, phosphoserine X 3 and Mn 303 (O3P 2.25 A, in the window; O2P
2.40 A, outside it); bicarbonate D 401 of 3WU2 with a Zn placed 2.0 A from
O3, on the C-O3 line, and a copy of it 10 A away without a metal.
"""
from libtbx.utils import null_out
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

# pdb_6iht with an H on O3P
pdb_6iht_hop3 = pdb_6iht.replace('HETATM 1964 MN ', (
  'HETATM 1962 HOP3 SEP X   3     -14.478   1.651  82.542  1.00 20.00'
  '           H\nHETATM 1964 MN '))
assert pdb_6iht_hop3 != pdb_6iht

NOT_PLACED = ('Acid H not placed: the O is bound to a metal, taken as '
  'deprotonated')
SUPPLIED = ('Supplied acid H at metal sites: kept as supplied; acids bound '
  'to metals are deprotonated')

def run(prefix, pdb_str, args=()):
  '''reduce2 on a model string; returns the names of the H of its model and
  the text of the description file'''
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
      names.add(a.name.strip())
  return names, description

def test_6iht_phosphoserine():
  '''The P-OH on Mn 303 (O3P) loses HOP3 and is listed; the P-OH outside
  the window (O2P) keeps HOP2'''
  names, description = run('tst_reduce2_metal_acid_6iht', pdb_6iht)
  assert 'HOP3' not in names and 'HOP2' in names, sorted(names)
  lines = description.splitlines()
  assert NOT_PLACED in lines, description
  assert lines[lines.index(NOT_PLACED) + 1] == ('  HOP3 SEP X   3  (O3P SEP X'
    '   3 - MN    MN A 303 2.25 A (geometric))'), description

def test_bicarbonate():
  '''Bicarbonate on Zn or Fe through O3 loses HO3 and is listed with the
  acids; the one without a metal keeps it'''
  for metal in ('ZN', 'FE'):
    pdb_str = pdb_bct_zn.replace('ZN    ZN', f'{metal}    {metal}').replace(
      '20.00          ZN', f'20.00          {metal}')
    names, description = run(f'tst_reduce2_metal_acid_bct_{metal}', pdb_str)
    assert 'HO3' in names, (metal, sorted(names))
    lines = description.splitlines()
    assert NOT_PLACED in lines, (metal, description)
    assert lines[lines.index(NOT_PLACED) + 1] == ('  HO3 BCT A   1  (O3  BCT'
      f' A   1 - {metal}    {metal} A   2 2.00 A (geometric))'), (metal,
      description)
    assert 'Phenol, enol or N-OH H' not in description, (metal,
      description)

def test_supplied_h_kept():
  '''A supplied HOP3 on the O on Mn stays and is listed, with add
  keep_existing_H=True and with optimize'''
  contact = ('  HOP3 SEP X   3  (O3P SEP X   3 - MN    MN A 303 2.25 A '
    '(geometric)); no metal-coordination proxy')
  for k, args in enumerate((['approach=add', 'keep_existing_H=True'],
                            ['approach=optimize'])):
    names, description = run(f'tst_reduce2_metal_acid_supplied_{k}',
      pdb_6iht_hop3, args)
    assert 'HOP3' in names, (args, sorted(names))
    lines = description.splitlines()
    assert SUPPLIED in lines, (args, description)
    assert lines[lines.index(SUPPLIED) + 1] == contact, (args, description)
    assert NOT_PLACED not in lines, (args, description)

if __name__ == '__main__':
  test_6iht_phosphoserine()
  test_bicarbonate()
  test_supplied_h_kept()
  print('OK')
