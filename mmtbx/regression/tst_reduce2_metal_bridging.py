"""reduce2 placement at O-H bridging two metals: an O-H whose O has one heavy
neighbour, a C or N, and is bound to two metals that coexist (blank or the
same altloc), neither alkali nor alkaline-earth, is taken as deprotonated, and
reduce2 does not build its H. Ca does not count, nor do alternate positions
of one metal.
Models: 8U7C, bis-tris D 207 with Tb D 202 (O4 2.21 A) and Er D 203 (O4
2.45 A, also O3 and O6); the same with Ca in place of Tb, and with Er in
altlocs A and B at the Er and Tb positions; 6LZB, hydroxamate EYF A 403 with
Co A 401 and A 402 (O1 1.91 and 1.99 A).
"""
pdb_8u7c = """\
CRYST1   37.514   40.892   50.150  92.62 106.14  98.87 P 1
HETATM 2307 TB    TB D 202       0.674  -0.668  -3.282  1.00 11.25          TB
HETATM 2308 ER   ER3 D 203      -0.352   0.798  -0.090  1.00 13.75          ER
HETATM 2330  N   BTB D 207      -0.304  -1.353   1.562  1.00  9.92           N
HETATM 2331  C1  BTB D 207       1.484  -3.198   1.769  1.00 11.03           C
HETATM 2332  C2  BTB D 207       1.055  -1.857   1.140  1.00  6.60           C
HETATM 2333  C3  BTB D 207       2.171  -0.798   1.403  1.00 11.00           C
HETATM 2334  C4  BTB D 207       1.063  -2.214  -0.328  1.00  8.09           C
HETATM 2335  C5  BTB D 207      -0.373  -0.766   2.931  1.00 12.64           C
HETATM 2336  C6  BTB D 207      -1.418   0.349   3.003  1.00 12.15           C
HETATM 2337  C7  BTB D 207      -1.370  -2.312   1.265  1.00 11.46           C
HETATM 2338  C8  BTB D 207      -2.353  -1.933   0.120  1.00 14.30           C
HETATM 2339  O1  BTB D 207       1.475  -3.088   3.181  1.00 16.04           O
HETATM 2340  O3  BTB D 207       1.824   0.527   1.045  1.00 11.81           O
HETATM 2341  O4  BTB D 207       0.634  -1.185  -1.135  1.00  8.52           O
HETATM 2342  O6  BTB D 207      -0.965   1.488   2.312  1.00 12.43           O
HETATM 2343  O8  BTB D 207      -2.446  -0.555  -0.164  1.00 15.12           O
END
"""

pdb_6lzb = """\
CRYST1   47.573   77.251   48.623  90.00  91.54  90.00 P 1 21 1
HETATM 2456 CO    CO A 401       0.630   0.402  12.184  1.00 16.09          CO
HETATM 2457 CO    CO A 402      -0.179  -0.701  15.125  1.00 12.80          CO
HETATM 2458  N   EYF A 403       0.789  -2.272  13.738  0.80 35.08           N
HETATM 2459  C   EYF A 403       0.718  -2.919  12.676  0.80 29.10           C
HETATM 2460  O   EYF A 403       0.804  -2.193  11.716  0.80 24.50           O
HETATM 2461  C1  EYF A 403       0.593  -4.418  12.484  0.80 30.77           C
HETATM 2462  C10 EYF A 403       4.195 -10.352  11.745  0.80 35.87           C
HETATM 2463  C11 EYF A 403       2.924  -9.788  11.789  0.80 31.05           C
HETATM 2464  C12 EYF A 403      -0.433  -8.758  12.698  0.80 28.23           C
HETATM 2465  C13 EYF A 403      -0.838  -7.633  13.210  0.80 26.04           C
HETATM 2466  C14 EYF A 403      -0.152  -6.663  12.779  0.80 25.90           C
HETATM 2467  C15 EYF A 403      -0.124  -5.348  13.265  0.80 27.35           C
HETATM 2468  C2  EYF A 403       1.262  -4.898  11.319  0.80 28.72           C
HETATM 2469  C3  EYF A 403       0.601  -7.112  11.628  0.80 29.07           C
HETATM 2470  C4  EYF A 403       0.770  -9.311  10.571  0.80 35.06           C
HETATM 2471  C5  EYF A 403       2.124  -9.882  10.663  0.80 34.57           C
HETATM 2472  C6  EYF A 403       2.556 -10.493   9.473  0.80 41.85           C
HETATM 2473  C7  EYF A 403       3.843 -11.035   9.413  0.80 48.39           C
HETATM 2474  C8  EYF A 403       3.648 -12.939   7.945  0.80 56.23           C
HETATM 2475  C9  EYF A 403       4.651 -10.971  10.560  0.80 41.23           C
HETATM 2476  N1  EYF A 403       1.307  -6.194  10.919  0.80 32.36           N
HETATM 2477  N2  EYF A 403       0.314  -8.417  11.592  0.80 33.14           N
HETATM 2478  O1  EYF A 403       0.989  -0.914  13.525  0.80 16.59           O
HETATM 2479  O2  EYF A 403       4.234 -11.665   8.252  0.80 51.44           O
END
"""

NOT_PLACED = 'Hydroxyl H not placed'

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

def test_8u7c_bis_tris():
  '''O4 of bis-tris, between Tb and Er, loses HO4; O3 and O6, on Er alone,
  keep HO3 and HO6'''
  names, description = run('tst_reduce2_metal_bridging_btb', pdb_8u7c)
  assert 'HO4' not in names and 'HO3' in names and 'HO6' in names, sorted(names)
  lines = [l for l in description.splitlines()
           if l.startswith('  HO4 BTB D 207  (')]
  assert len(lines) == 1, description
  assert 'TB    TB D 202 2.21 A' in lines[0], lines[0]
  assert 'ER   ER3 D 203 2.45 A' in lines[0], lines[0]

def test_8u7c_controls():
  '''O4 keeps HO4 with Ca in place of Tb, and with Er in altlocs A and B at
  the Er and Tb positions'''
  pdb_ca = pdb_8u7c.replace('TB    TB D 202', 'CA    CA D 202').replace(
    '          TB\n', '          CA\n')
  assert ' TB' not in pdb_ca
  names, description = run('tst_reduce2_metal_bridging_ca', pdb_ca)
  assert 'HO4' in names, sorted(names)
  assert 'Hydroxyl H' not in description, description
  metals = ('HETATM 2307 TB    TB D 202       0.674  -0.668  -3.282  1.00 11.25'
    '          TB\n'
    'HETATM 2308 ER   ER3 D 203      -0.352   0.798  -0.090  1.00 13.75'
    '          ER\n')
  alternates = ('HETATM 2307 ER  AER3 D 203      -0.352   0.798  -0.090  0.50 '
    '13.75          ER\n'
    'HETATM 2308 ER  BER3 D 203       0.674  -0.668  -3.282  0.50 11.25'
    '          ER\n')
  assert metals in pdb_8u7c
  names, description = run('tst_reduce2_metal_bridging_alt',
    pdb_8u7c.replace(metals, alternates))
  assert 'HO4' in names, sorted(names)
  assert 'Hydroxyl H' not in description, description

def test_6lzb_hydroxamate():
  '''The N-O oxygen of the hydroxamate, between two Co, loses H9; the N keeps
  H8'''
  names, description = run('tst_reduce2_metal_bridging_eyf', pdb_6lzb)
  assert 'H9' not in names and 'H8' in names, sorted(names)
  assert NOT_PLACED in description, description

if __name__ == '__main__':
  test_8u7c_bis_tris()
  test_8u7c_controls()
  test_6lzb_hydroxamate()
  print('OK')
