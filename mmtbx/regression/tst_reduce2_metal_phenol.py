"""reduce2 placement at phenols bound to metal ions of charge +3 or more: an
O-H whose O has one heavy neighbour, a C in a double, triple, aromatic or
delocalized bond (phenol, enol; acids: tst_reduce2_metal_acid.py) or an N
that is in such a bond or bonded to such an atom (hydroxamic acid, N-hydroxy
heterocycle), bound to a metal of PHENOLATE_METAL_ELEMENTS is taken as
deprotonated (Fe as Fe(III)), and reduce2 does not build its H. At other
metals the H stays (a phenol bound to Zn has a pKa of about 8.5). Supplied H
are kept with a note.
Models: 8RK1 (FutA with Fe(III), neutron), Tyr 13, 143, 199 and 200 on
Fe 401 (1.93-2.19 A); none of the four OH carries D in the neutron model.
4X1D, Tyr 426 and 517 on Yb 701 (2.37 and 2.40 A). 1Y93, acetohydroxamate on
Zn 264 (O 2.07 A), and the same with Fe in place of Zn.
"""
pdb_8rk1 = """\
CRYST1   39.542   78.316   47.893  90.00  97.36  90.00 P 1 21 1
ATOM    192  N   TYR A  13      23.356  -5.119  10.678  1.00 23.60           N
ATOM    193  CA  TYR A  13      22.877  -6.489  10.653  1.00 24.12           C
ATOM    194  C   TYR A  13      23.957  -7.423  10.098  1.00 27.53           C
ATOM    195  O   TYR A  13      23.844  -8.639  10.219  1.00 27.32           O
ATOM    196  CB  TYR A  13      21.596  -6.472   9.819  1.00 24.29           C
ATOM    197  CG  TYR A  13      20.430  -5.752  10.490  1.00 23.50           C
ATOM    198  CD1 TYR A  13      19.912  -6.221  11.699  1.00 23.00           C
ATOM    199  CD2 TYR A  13      19.855  -4.604   9.940  1.00 23.36           C
ATOM    200  CE1 TYR A  13      18.859  -5.581  12.336  1.00 21.91           C
ATOM    201  CE2 TYR A  13      18.802  -3.955  10.570  1.00 22.19           C
ATOM    202  CZ  TYR A  13      18.300  -4.454  11.762  1.00 21.50           C
ATOM    203  OH  TYR A  13      17.260  -3.814  12.362  1.00 21.87           O
ATOM   2580  N   TYR A 143      12.797   1.748  18.814  1.00 11.27           N
ATOM   2581  CA  TYR A 143      11.665   0.989  18.311  1.00 10.99           C
ATOM   2582  C   TYR A 143      10.446   1.285  19.180  1.00 10.74           C
ATOM   2583  O   TYR A 143       9.385   1.628  18.656  1.00 10.13           O
ATOM   2584  CB  TYR A 143      12.033  -0.499  18.356  1.00 11.05           C
ATOM   2585  CG  TYR A 143      13.160  -0.940  17.432  1.00 11.60           C
ATOM   2586  CD1 TYR A 143      13.594  -0.141  16.381  1.00 12.31           C
ATOM   2587  CD2 TYR A 143      13.808  -2.149  17.640  1.00 11.76           C
ATOM   2588  CE1 TYR A 143      14.622  -0.541  15.544  1.00 12.92           C
ATOM   2589  CE2 TYR A 143      14.845  -2.560  16.825  1.00 12.35           C
ATOM   2590  CZ  TYR A 143      15.251  -1.757  15.770  1.00 13.22           C
ATOM   2591  OH  TYR A 143      16.257  -2.180  14.951  1.00 13.34           O
ATOM   3536  N   TYR A 199      10.911  -8.879  11.791  1.00 12.96           N
ATOM   3537  CA  TYR A 199      12.039  -9.496  11.099  1.00 12.89           C
ATOM   3538  C   TYR A 199      12.918 -10.344  12.022  1.00 13.25           C
ATOM   3539  O   TYR A 199      13.518 -11.304  11.559  1.00 13.55           O
ATOM   3540  CB  TYR A 199      12.891  -8.477  10.337  1.00 12.84           C
ATOM   3541  CG  TYR A 199      13.479  -7.323  11.136  1.00 12.22           C
ATOM   3542  CD1 TYR A 199      14.685  -7.425  11.815  1.00 11.95           C
ATOM   3543  CD2 TYR A 199      12.811  -6.111  11.197  1.00 12.18           C
ATOM   3544  CE1 TYR A 199      15.192  -6.360  12.547  1.00 12.10           C
ATOM   3545  CE2 TYR A 199      13.296  -5.041  11.934  1.00 12.09           C
ATOM   3546  CZ  TYR A 199      14.492  -5.163  12.613  1.00 11.96           C
ATOM   3547  OH  TYR A 199      14.950  -4.091  13.318  1.00 12.22           O
ATOM   3557  N   TYR A 200      12.963 -10.037  13.324  1.00 13.46           N
ATOM   3558  CA  TYR A 200      13.689 -10.854  14.287  1.00 13.76           C
ATOM   3559  C   TYR A 200      13.044 -12.231  14.447  1.00 13.89           C
ATOM   3560  O   TYR A 200      13.746 -13.219  14.659  1.00 14.04           O
ATOM   3561  CB  TYR A 200      13.747 -10.185  15.663  1.00 14.03           C
ATOM   3562  CG  TYR A 200      14.561  -8.897  15.694  1.00 14.59           C
ATOM   3563  CD1 TYR A 200      15.905  -8.901  15.323  1.00 14.30           C
ATOM   3564  CD2 TYR A 200      13.999  -7.688  16.101  1.00 14.41           C
ATOM   3565  CE1 TYR A 200      16.661  -7.741  15.335  1.00 14.21           C
ATOM   3566  CE2 TYR A 200      14.757  -6.532  16.146  1.00 14.12           C
ATOM   3567  CZ  TYR A 200      16.083  -6.558  15.757  1.00 14.09           C
ATOM   3568  OH  TYR A 200      16.825  -5.413  15.775  1.00 14.13           O
HETATM 5556 FE    FE A 401      16.613  -4.028  14.445  1.00 11.87          FE
END
"""

pdb_4x1d = """\
CRYST1  156.502  136.863  107.843  90.00  90.41  90.00 C 1 2 1
ATOM   3166  N   TYR A 426     -55.552 -31.359 110.039  1.00 67.28           N
ATOM   3167  CA  TYR A 426     -56.245 -30.745 108.890  1.00 64.66           C
ATOM   3168  C   TYR A 426     -56.964 -31.817 108.022  1.00 62.14           C
ATOM   3169  O   TYR A 426     -56.929 -33.011 108.339  1.00 59.75           O
ATOM   3170  CB  TYR A 426     -55.307 -29.789 108.076  1.00 67.78           C
ATOM   3171  CG  TYR A 426     -54.161 -30.438 107.295  1.00 71.45           C
ATOM   3172  CD1 TYR A 426     -54.415 -31.137 106.136  1.00 73.02           C
ATOM   3173  CD2 TYR A 426     -52.825 -30.327 107.706  1.00 72.96           C
ATOM   3174  CE1 TYR A 426     -53.404 -31.743 105.420  1.00 74.39           C
ATOM   3175  CE2 TYR A 426     -51.797 -30.929 106.985  1.00 70.57           C
ATOM   3176  CZ  TYR A 426     -52.100 -31.646 105.834  1.00 73.64           C
ATOM   3177  OH  TYR A 426     -51.166 -32.304 105.047  1.00 62.65           O
ATOM   3862  N   TYR A 517     -50.433 -27.097  97.158  1.00 49.74           N
ATOM   3863  CA  TYR A 517     -51.879 -27.263  97.296  1.00 47.71           C
ATOM   3864  C   TYR A 517     -52.563 -26.381  96.267  1.00 43.89           C
ATOM   3865  O   TYR A 517     -53.418 -26.841  95.511  1.00 44.24           O
ATOM   3866  CB  TYR A 517     -52.408 -26.858  98.659  1.00 47.66           C
ATOM   3867  CG  TYR A 517     -51.804 -27.462  99.891  1.00 45.06           C
ATOM   3868  CD1 TYR A 517     -52.279 -28.677 100.418  1.00 44.74           C
ATOM   3869  CD2 TYR A 517     -50.830 -26.764 100.592  1.00 43.23           C
ATOM   3870  CE1 TYR A 517     -51.761 -29.194 101.598  1.00 43.94           C
ATOM   3871  CE2 TYR A 517     -50.300 -27.264 101.749  1.00 44.30           C
ATOM   3872  CZ  TYR A 517     -50.766 -28.469 102.257  1.00 48.05           C
ATOM   3873  OH  TYR A 517     -50.203 -28.910 103.436  1.00 56.79           O
HETATM10223 YB    YB A 701     -49.426 -31.118 103.955  1.00 49.88          YB
END
"""

pdb_1y93 = """\
CRYST1   50.913   59.552   53.497  90.00 115.14  90.00 C 1 2 1
HETATM 2416 ZN    ZN A 264       4.199  -1.587   4.821  1.00  6.40          ZN
HETATM 2421  N   HAE A 301       1.361  -2.059   5.185  1.00  8.42           N
HETATM 2422  O   HAE A 301       2.508  -1.997   5.938  1.00 11.90           O
HETATM 2423  C1  HAE A 301       0.171  -1.388   3.191  1.00 14.58           C
HETATM 2424  C2  HAE A 301       1.448  -1.511   3.991  1.00 12.17           C
HETATM 2425  O2  HAE A 301       2.529  -1.084   3.591  1.00 12.21           O
END
"""

NOT_PLACED = 'Phenol H not placed'
SUPPLIED = 'Supplied phenol H at metal sites'

def run(prefix, pdb_str, args=()):
  """Run reduce2 on a model string in the current directory.

  Parameters
  ----------
  prefix : str
    Prefix of the input, output and description files written.
  pdb_str : str
    The input model in PDB format.
  args : sequence of str
    Further reduce2 arguments.

  Returns
  -------
  tuple
    The hierarchy of the model that reduce2 returns and the text of the
    description file.
  """
  from libtbx.utils import null_out
  from iotbx.cli_parser import run_program
  from mmtbx.programs import reduce2
  fn, txt = prefix + '_in.pdb', prefix + '.txt'
  with open(fn, 'w') as f:
    f.write(pdb_str)
  result = run_program(program_class=reduce2.Program, logger=null_out(),
    args=[fn, 'output.filename=' + prefix + '_out.pdb',
          'output.description_file_name=' + txt, '--overwrite'] + list(args))
  with open(txt) as f:
    description = f.read()
  return result.model.get_hierarchy(), description

def section(text, title):
  """Return the indented lines under the report heading that starts with title.

  Parameters
  ----------
  text : str
    The description file.
  title : str
    Start of the heading.

  Returns
  -------
  list of str
    The lines of the section, stripped.
  """
  lines = text.splitlines()
  starts = [i for i, l in enumerate(lines) if l.startswith(title)]
  assert len(starts) == 1, (title, text)
  result = []
  for l in lines[starts[0] + 1:]:
    if not l.startswith('  ') or l.strip().startswith('Time to'):
      break
    result.append(l.strip())
  return result

def phenol_h(hierarchy):
  """Return the residue numbers of the Tyr that carry HH.

  Parameters
  ----------
  hierarchy : iotbx.pdb.hierarchy.root

  Returns
  -------
  list of int
  """
  return sorted(a.parent().parent().resseq_as_int() for a in hierarchy.atoms()
                if a.name.strip() == 'HH')

def test_8rk1_fe():
  '''Four Tyr on Fe(III): reduce2 places no HH and lists the four'''
  h, description = run('tst_reduce2_metal_phenol_fe', pdb_8rk1)
  assert phenol_h(h) == [], phenol_h(h)
  lines = section(description, NOT_PLACED)
  assert len(lines) == 4, lines
  assert lines[0].startswith('HH  TYR A  13'), lines

def test_8rk1_zn():
  '''The same site with Zn: reduce2 places the four HH'''
  pdb_str = pdb_8rk1.replace('FE    FE A 401', 'ZN    ZN A 401').replace(
    '11.87          FE', '11.87          ZN')
  h, description = run('tst_reduce2_metal_phenol_zn', pdb_str)
  assert phenol_h(h) == [13, 143, 199, 200], phenol_h(h)
  assert NOT_PLACED not in description, description

def test_8rk1_supplied():
  '''approach=optimize keeps supplied HH on Tyr bound to Fe, with a note'''
  pdb_str = pdb_8rk1.replace('FE    FE A 401', 'ZN    ZN A 401').replace(
    '11.87          FE', '11.87          ZN')
  h, description = run('tst_reduce2_metal_phenol_h', pdb_str)
  with_h = pdb_8rk1.splitlines()[0] + '\n' + h.as_pdb_string()
  with_h = with_h.replace('ZN    ZN A 401', 'FE    FE A 401')
  with_h = with_h.replace('          ZN', '          FE')
  h, description = run('tst_reduce2_metal_phenol_supplied', with_h,
    ['approach=optimize'])
  assert phenol_h(h) == [13, 143, 199, 200], phenol_h(h)
  assert len([l for l in section(description, SUPPLIED)
              if l.startswith('HH')]) == 4, description

def test_4x1d_yb():
  '''Two Tyr on Yb: reduce2 places no HH and lists both'''
  h, description = run('tst_reduce2_metal_phenol_yb', pdb_4x1d)
  assert phenol_h(h) == [], phenol_h(h)
  assert len(section(description, NOT_PLACED)) == 2, description

def test_1y93_hydroxamate():
  '''Acetohydroxamate on Fe loses HO and keeps HN; on Zn it keeps both'''
  pdb_fe = pdb_1y93.replace('ZN    ZN A 264', 'FE    FE A 264').replace(
    '          ZN\n', '          FE\n')
  assert ' ZN' not in pdb_fe
  for pdb_str, metal, expected in ((pdb_fe, 'fe', ['HN']),
                                   (pdb_1y93, 'zn', ['HN', 'HO'])):
    h, description = run('tst_reduce2_metal_phenol_hae_' + metal, pdb_str)
    names = sorted(a.name.strip() for a in h.atoms()
                   if a.name.strip() in ('HN', 'HO'))
    assert names == expected, (metal, names)

if __name__ == '__main__':
  test_8rk1_fe()
  test_8rk1_zn()
  test_8rk1_supplied()
  test_4x1d_yb()
  test_1y93_hydroxamate()
  print('OK')
