"""reduce2 placement at acids bound to metals: an O-H whose O has one heavy
neighbour, an oxoacid centre (carboxylic, phosphoric or phosphonic, sulfonic
or arsenic acid), and is bound to a metal is taken as deprotonated, and
reduce2 does not build its H. Bicarbonate keeps its H.
Models: 6IHT, phosphoserine X 3 and Mn 303 (O3P 2.25 A); bicarbonate D 401 of
3WU2 with a Zn placed 2.0 A from O3, on the C-O3 line.
"""
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
END
"""

NOT_PLACED = 'Acid H not placed'

def run(prefix, pdb_str):
  """Run reduce2 with its defaults on a model string in the current directory.

  Parameters
  ----------
  prefix : str
    Prefix of the input, output and description files written.
  pdb_str : str
    The input model in PDB format.

  Returns
  -------
  tuple
    The names of the H atoms of the model that reduce2 returns and the text
    of the description file.
  """
  from libtbx.utils import null_out
  from iotbx.cli_parser import run_program
  from mmtbx.programs import reduce2
  fn, txt = prefix + '_in.pdb', prefix + '.txt'
  with open(fn, 'w') as f:
    f.write(pdb_str)
  result = run_program(program_class=reduce2.Program, logger=null_out(),
    args=[fn, 'output.filename=' + prefix + '_out.pdb',
          'output.description_file_name=' + txt, '--overwrite'])
  with open(txt) as f:
    description = f.read()
  names = set()
  for a in result.model.get_hierarchy().atoms():
    if a.element_is_hydrogen():
      names.add(a.name.strip())
  return names, description

def test_6iht_phosphoserine():
  '''The P-OH on Mn 303 (O3P) loses HOP3; the free P-OH keeps HOP2'''
  names, description = run('tst_reduce2_metal_acid_6iht', pdb_6iht)
  assert 'HOP3' not in names and 'HOP2' in names, sorted(names)
  assert ('HOP3 SEP X   3  (O3P SEP X   3 - MN    MN A 303 2.25 A '
    '(geometric))') in description, description

def test_bicarbonate_kept():
  '''Bicarbonate on Zn through O3 keeps HO3'''
  names, description = run('tst_reduce2_metal_acid_bct', pdb_bct_zn)
  assert 'HO3' in names, sorted(names)
  assert NOT_PLACED not in description, description

if __name__ == '__main__':
  test_6iht_phosphoserine()
  test_bicarbonate_kept()
  print('OK')
