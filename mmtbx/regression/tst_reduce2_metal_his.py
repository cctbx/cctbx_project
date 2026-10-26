"""reduce2 Optimizer: a His ring N bound to a positive ion in a symmetry mate
carries no H, as for an ion in the model; the other ring N keeps its H.
Model: 6GTJ (neutron), His 306 and Cu 502, whose symmetry mate
(-y+1/2,-z,x-1/2) is 2.22 A from NE2. The neutron model has no D on NE2;
without the symmetry search reduce2 placed HE2 pointing at the Cu.
"""
pdb_6gtj = """\
CRYST1   97.980   97.980   97.980  90.00  90.00  90.00 P 21 3
ATOM   4587  N   HIS A 306      21.626  12.904 -13.586  1.00 10.46           N
ATOM   4588  CA  HIS A 306      21.568  11.718 -14.419  1.00  9.17           C
ATOM   4589  C   HIS A 306      22.952  11.160 -14.731  1.00 11.35           C
ATOM   4590  O   HIS A 306      23.112  10.463 -15.739  1.00  9.07           O
ATOM   4591  CB  HIS A 306      20.797  11.992 -15.705  1.00  8.50           C
ATOM   4592  CG  HIS A 306      20.121  10.770 -16.251  1.00  7.91           C
ATOM   4593  ND1 HIS A 306      20.798   9.658 -16.705  1.00  8.73           N
ATOM   4594  CD2 HIS A 306      18.804  10.470 -16.346  1.00 10.75           C
ATOM   4595  CE1 HIS A 306      19.926   8.741 -17.086  1.00  9.71           C
ATOM   4596  NE2 HIS A 306      18.708   9.209 -16.877  1.00  7.35           N
TER
HETATM 5111 CU    CU A 502      31.875  32.121  -7.996  1.00  9.68          CU
END
"""

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
  set of str
    Names of the H atoms of the model that reduce2 returns.
  """
  from libtbx.utils import null_out
  from iotbx.cli_parser import run_program
  from mmtbx.programs import reduce2
  fn = prefix + '_in.pdb'
  with open(fn, 'w') as f:
    f.write(pdb_str)
  result = run_program(program_class=reduce2.Program, logger=null_out(),
    args=[fn, 'output.filename=' + prefix + '_out.pdb',
          'output.description_file_name=' + prefix + '.txt', '--overwrite'])
  names = set()
  for a in result.model.get_hierarchy().atoms():
    if a.element_is_hydrogen():
      names.add(a.name.strip())
  return names

def test_6gtj_cu_symmetry_mate():
  '''His NE2 bound to the symmetry mate of Cu 502: reduce2 places HD1 and not
  HE2'''
  names = run('tst_reduce2_metal_his_6gtj', pdb_6gtj)
  assert 'HD1' in names and 'HE2' not in names, sorted(names)

def test_6gtj_no_cell():
  '''Without CRYST1 no symmetry mate is searched: the Cu is 30 A away and
  reduce2 places HD1 and HE2'''
  names = run('tst_reduce2_metal_his_6gtj_no_cell', pdb_6gtj.split('\n', 1)[1])
  assert 'HD1' in names and 'HE2' in names, sorted(names)

if __name__ == '__main__':
  test_6gtj_cu_symmetry_mate()
  test_6gtj_no_cell()
  print('OK')
