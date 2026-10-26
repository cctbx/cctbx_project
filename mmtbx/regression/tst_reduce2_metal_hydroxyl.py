"""reduce2 Optimizer: a rotatable H (OH, SH, SeH) whose parent is bound to a
metal is not turned towards it (metal-X-H below 90 degrees).
Model: 1MNZ, residues within 4 A of Mg 389 and Tris 392, whose O2 is on the
Mg (2.31 A). Without the rule, HO2 was placed 1.7 A from the Mg (Mg-O2-HO2
about 40 degrees).
"""
pdb_1mnz = """\
CRYST1   93.880   99.680  102.900  90.00  90.00  90.00 I 2 2 2
ATOM    120  N   TRP A  16      36.520  25.459  63.264  1.00  8.60           N
ATOM    121  CA  TRP A  16      37.638  25.635  62.338  1.00  8.12           C
ATOM    122  C   TRP A  16      38.738  26.600  62.804  1.00  8.25           C
ATOM    123  O   TRP A  16      39.761  26.755  62.105  1.00  9.70           O
ATOM    124  CB  TRP A  16      37.084  26.170  60.971  1.00  8.52           C
ATOM    125  CG  TRP A  16      36.544  27.597  61.087  1.00  8.54           C
ATOM    126  CD1 TRP A  16      37.277  28.742  60.854  1.00  9.24           C
ATOM    127  CD2 TRP A  16      35.220  28.039  61.456  1.00  8.62           C
ATOM    128  NE1 TRP A  16      36.508  29.882  61.042  1.00  9.58           N
ATOM    129  CE2 TRP A  16      35.245  29.482  61.405  1.00  8.72           C
ATOM    130  CE3 TRP A  16      33.981  27.430  61.774  1.00  8.63           C
ATOM    131  CZ2 TRP A  16      34.115  30.257  61.714  1.00  9.40           C
ATOM    132  CZ3 TRP A  16      32.888  28.197  62.058  1.00  9.20           C
ATOM    133  CH2 TRP A  16      32.920  29.596  62.014  1.00  9.44           C
ATOM    426  N   HIS A  54      30.941  25.089  60.159  1.00  9.05           N
ATOM    427  CA  HIS A  54      30.530  24.791  58.796  1.00  9.16           C
ATOM    428  C   HIS A  54      30.043  23.327  58.743  1.00  9.74           C
ATOM    429  O   HIS A  54      30.631  22.475  59.437  1.00 10.37           O
ATOM    430  CB  HIS A  54      31.691  24.907  57.775  1.00  9.64           C
ATOM    431  CG  HIS A  54      32.374  26.255  57.723  1.00  9.11           C
ATOM    432  ND1 HIS A  54      33.365  26.613  58.583  1.00  9.90           N
ATOM    433  CD2 HIS A  54      32.234  27.356  56.900  1.00  9.79           C
ATOM    434  CE1 HIS A  54      33.751  27.857  58.306  1.00  9.83           C
ATOM    435  NE2 HIS A  54      33.088  28.362  57.306  1.00 10.23           N
ATOM    744  N   PHE A  94      33.328  26.974  49.745  1.00  9.23           N
ATOM    745  CA  PHE A  94      34.600  27.490  49.253  1.00  9.19           C
ATOM    746  C   PHE A  94      35.763  26.504  49.171  1.00  9.36           C
ATOM    747  O   PHE A  94      36.681  26.756  48.350  1.00  9.54           O
ATOM    748  CB  PHE A  94      35.049  28.758  50.093  1.00  9.53           C
ATOM    749  CG  PHE A  94      35.035  28.506  51.606  1.00  9.27           C
ATOM    750  CD1 PHE A  94      36.116  27.953  52.238  1.00 10.23           C
ATOM    751  CD2 PHE A  94      33.892  28.842  52.342  1.00  9.56           C
ATOM    752  CE1 PHE A  94      36.079  27.734  53.619  1.00 11.09           C
ATOM    753  CE2 PHE A  94      33.859  28.593  53.703  1.00 10.69           C
ATOM    754  CZ  PHE A  94      34.935  28.035  54.335  1.00 11.01           C
ATOM   1090  N   TRP A 137      29.188  33.441  54.097  1.00  9.40           N
ATOM   1091  CA  TRP A 137      30.400  34.126  53.616  1.00  8.29           C
ATOM   1092  C   TRP A 137      30.666  33.642  52.211  1.00  8.42           C
ATOM   1093  O   TRP A 137      31.113  32.493  51.989  1.00  9.33           O
ATOM   1094  CB  TRP A 137      31.650  33.819  54.498  1.00  9.33           C
ATOM   1095  CG  TRP A 137      32.910  34.378  53.901  1.00  8.65           C
ATOM   1096  CD1 TRP A 137      33.215  35.679  53.579  1.00 10.09           C
ATOM   1097  CD2 TRP A 137      34.094  33.624  53.519  1.00  8.27           C
ATOM   1098  NE1 TRP A 137      34.479  35.791  53.039  1.00 10.29           N
ATOM   1099  CE2 TRP A 137      35.048  34.530  52.997  1.00  8.43           C
ATOM   1100  CE3 TRP A 137      34.452  32.270  53.595  1.00  9.92           C
ATOM   1101  CZ2 TRP A 137      36.318  34.103  52.538  1.00  8.43           C
ATOM   1102  CZ3 TRP A 137      35.690  31.839  53.147  1.00 10.06           C
ATOM   1103  CH2 TRP A 137      36.611  32.767  52.628  1.00  9.03           C
ATOM   1440  N   GLU A 181      28.262  37.889  55.661  1.00  8.56           N
ATOM   1441  CA  GLU A 181      29.599  38.261  56.138  1.00  8.67           C
ATOM   1442  C   GLU A 181      30.427  38.877  55.014  1.00  8.96           C
ATOM   1443  O   GLU A 181      30.889  38.110  54.138  1.00 10.81           O
ATOM   1444  CB  GLU A 181      30.326  37.050  56.714  1.00  9.47           C
ATOM   1445  CG  GLU A 181      31.717  37.397  57.258  1.00  9.70           C
ATOM   1446  CD  GLU A 181      32.337  36.296  58.084  1.00  9.82           C
ATOM   1447  OE1 GLU A 181      31.613  35.343  58.472  1.00 11.75           O
ATOM   1448  OE2 GLU A 181      33.586  36.381  58.336  1.00 10.73           O
ATOM   1731  N   GLU A 217      33.218  41.215  60.192  1.00  7.18           N
ATOM   1732  CA  GLU A 217      34.456  40.673  60.814  1.00  6.79           C
ATOM   1733  C   GLU A 217      35.479  41.824  60.901  1.00  6.90           C
ATOM   1734  O   GLU A 217      35.767  42.522  59.902  1.00  6.89           O
ATOM   1735  CB  GLU A 217      35.029  39.504  59.997  1.00  7.32           C
ATOM   1736  CG  GLU A 217      36.329  38.941  60.615  1.00  7.61           C
ATOM   1737  CD  GLU A 217      36.949  37.728  59.910  1.00  7.88           C
ATOM   1738  OE1 GLU A 217      36.220  36.767  59.606  1.00  9.28           O
ATOM   1739  OE2 GLU A 217      38.175  37.818  59.611  1.00  8.85           O
ATOM   1751  N   HIS A 220      38.227  41.947  58.203  1.00  7.03           N
ATOM   1752  CA  HIS A 220      37.775  42.108  56.809  1.00  6.92           C
ATOM   1753  C   HIS A 220      37.573  43.597  56.467  1.00  6.71           C
ATOM   1754  O   HIS A 220      38.047  44.056  55.424  1.00  7.52           O
ATOM   1755  CB  HIS A 220      36.504  41.292  56.555  1.00  6.93           C
ATOM   1756  CG  HIS A 220      36.780  39.797  56.611  1.00  7.14           C
ATOM   1757  ND1 HIS A 220      35.818  38.812  56.464  1.00  7.87           N
ATOM   1758  CD2 HIS A 220      37.947  39.102  56.832  1.00  7.84           C
ATOM   1759  CE1 HIS A 220      36.415  37.604  56.606  1.00  7.92           C
ATOM   1760  NE2 HIS A 220      37.714  37.711  56.810  1.00  8.47           N
ATOM   1948  N   ASP A 245      32.001  38.778  64.054  1.00  8.57           N
ATOM   1949  CA  ASP A 245      33.247  38.003  64.118  1.00  8.47           C
ATOM   1950  C   ASP A 245      34.364  38.927  64.643  1.00  7.54           C
ATOM   1951  O   ASP A 245      34.685  39.978  64.042  1.00  8.30           O
ATOM   1952  CB  ASP A 245      33.703  37.433  62.765  1.00  9.27           C
ATOM   1953  CG  ASP A 245      33.057  36.168  62.290  1.00  8.98           C
ATOM   1954  OD1 ASP A 245      32.175  35.603  62.915  1.00 10.07           O
ATOM   1955  OD2 ASP A 245      33.511  35.653  61.196  1.00 12.00           O
ATOM   2281  N   ASP A 287      36.488  34.483  65.233  1.00  8.36           N
ATOM   2282  CA  ASP A 287      36.923  33.831  63.989  1.00  8.10           C
ATOM   2283  C   ASP A 287      38.410  34.122  63.851  1.00  7.55           C
ATOM   2284  O   ASP A 287      38.799  35.169  63.306  1.00  8.85           O
ATOM   2285  CB  ASP A 287      36.124  34.323  62.777  1.00  8.12           C
ATOM   2286  CG  ASP A 287      36.583  33.565  61.504  1.00  8.51           C
ATOM   2287  OD1 ASP A 287      37.353  32.592  61.608  1.00  9.02           O
ATOM   2288  OD2 ASP A 287      36.146  33.998  60.394  1.00 10.30           O
TER
HETATM 3076 MG    MG A 389      34.719  35.310  59.621  1.00  8.72          MG
HETATM 3078  N   TRS A 392      33.373  33.550  59.189  1.00 10.61           N
HETATM 3079  C   TRS A 392      33.735  32.634  58.095  1.00 12.66           C
HETATM 3080  C1  TRS A 392      34.731  31.531  58.435  1.00 18.88           C
HETATM 3081  C2  TRS A 392      34.369  33.522  57.041  1.00 14.10           C
HETATM 3082  C3  TRS A 392      32.500  31.852  57.652  1.00 12.06           C
HETATM 3083  O1  TRS A 392      36.136  31.939  58.199  1.00 33.72           O
HETATM 3084  O2  TRS A 392      35.424  34.387  57.627  1.00 12.13           O
HETATM 3085  O3  TRS A 392      32.743  31.018  56.515  1.00 11.60           O
HETATM 3115  O   HOH A1024      38.402  34.599  58.475  1.00 11.37           O
HETATM 3141  O   HOH A1050      30.383  30.638  55.083  1.00 13.42           O
HETATM 3162  O   HOH A1071      31.894  32.988  61.504  1.00 16.22           O
HETATM 3234  O   HOH A1146      38.726  30.997  57.974  1.00 28.52           O
HETATM 3271  O   HOH A1186      38.206  34.527  55.348  1.00 22.20           O
END
"""

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
  iotbx.pdb.hierarchy.root
    The hierarchy of the model that reduce2 returns.
  """
  from libtbx.utils import null_out
  from iotbx.cli_parser import run_program
  from mmtbx.programs import reduce2
  fn = prefix + '_in.pdb'
  with open(fn, 'w') as f:
    f.write(pdb_str)
  result = run_program(program_class=reduce2.Program, logger=null_out(),
    args=[fn, 'output.filename=' + prefix + '_out.pdb',
          'output.description_file_name=' + prefix + '.txt', '--overwrite']
    + list(args))
  return result.model.get_hierarchy()

def metal_angle(hierarchy, resname, x_name, h_name, metal):
  """Return the metal-X-H angle of one hydrogen.

  Parameters
  ----------
  hierarchy : iotbx.pdb.hierarchy.root
    Model with one residue named resname and one metal atom.
  resname : str
    Name of the residue that holds X and the H.
  x_name, h_name : str
    Atom names of X and of the H on it.
  metal : str
    Residue and atom name of the metal.

  Returns
  -------
  float
    The angle in degrees.
  """
  from scitbx import matrix
  sites = {}
  for a in hierarchy.atoms():
    sites[(a.parent().resname.strip(), a.name.strip())] = matrix.col(a.xyz)
  m, x, h = (sites[(metal, metal)], sites[(resname, x_name)],
             sites[(resname, h_name)])
  return (m - x).angle(h - x, deg=True)

def test_1mnz_tris_o2():
  '''A hydroxyl on a metal keeps its H off the metal: reduce2 on 1MNZ, with
  and without the Tris amine placed as NH2, places HO2 with Mg-O2-HO2 of at
  least 90 degrees'''
  for args in ([], ['metal_bound_amines_as_nh2=True']):
    h = run('tst_reduce2_metal_hydroxyl_1mnz', pdb_1mnz, args)
    angle = metal_angle(h, 'TRS', 'O2', 'HO2', 'MG')
    assert angle >= 90, (args, angle)

def test_1mnz_mg_symmetry_mate():
  '''The metal can be in a symmetry mate: with Mg 389 moved to its equivalent
  under -x,-y,z, reduce2 places HO2 with at least 90 degrees between O2-HO2
  and O2 to the Mg image bound to O2'''
  moved = pdb_1mnz.replace('MG A 389      34.719  35.310  59.621',
                           'MG A 389     -34.719 -35.310  59.621')
  assert moved != pdb_1mnz
  h = run('tst_reduce2_metal_hydroxyl_1mnz_sym', moved)
  # put the Mg at its image bound to O2, where it was in 1MNZ
  for a in h.atoms():
    if a.name.strip() == 'MG':
      a.set_xyz((34.719, 35.310, 59.621))
  angle = metal_angle(h, 'TRS', 'O2', 'HO2', 'MG')
  assert angle >= 90, angle

if __name__ == '__main__':
  test_1mnz_tris_o2()
  test_1mnz_mg_symmetry_mate()
  print('OK')
