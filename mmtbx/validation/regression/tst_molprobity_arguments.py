from __future__ import absolute_import, division, print_function
from libtbx import easy_run
from libtbx.test_utils import approx_equal

pdb_str = '''
CRYST1   35.050   40.500   42.370  90.00  90.00  90.00 P 21 21 21
ATOM    868  N   THR A  58      12.684  27.992  19.917  1.00  4.88           N
ATOM    869  CA  THR A  58      12.449  29.168  20.750  1.00  5.30           C
ATOM    870  C   THR A  58      12.794  30.387  19.944  1.00  5.18           C
ATOM    871  O   THR A  58      13.829  30.407  19.238  1.00  7.53           O
ATOM    872  CB  THR A  58      13.340  29.099  21.992  1.00  6.13           C
ATOM    873  OG1 THR A  58      13.131  27.860  22.680  1.00  8.07           O
ATOM    874  CG2 THR A  58      13.019  30.214  22.977  1.00  8.02           C
ATOM    875  H   THR A  58      13.565  27.890  19.407  1.00  5.87           H
ATOM    876  HA  THR A  58      11.405  29.177  21.037  1.00  6.36           H
ATOM    877  HB  THR A  58      14.387  29.172  21.714  1.00  7.36           H
ATOM    878  HG1 THR A  58      13.227  27.203  21.983  1.00  9.69           H
ATOM    879 HG21 THR A  58      13.587  30.061  23.889  1.00  9.63           H
ATOM    880 HG22 THR A  58      13.295  31.165  22.548  1.00  9.63           H
ATOM    881 HG23 THR A  58      11.958  30.232  23.192  1.00  9.63           H
HETATM  931  O   HOH A  63      14.819  25.768  21.599  1.00  8.43           O
HETATM  932  H1  HOH A  63      15.358  25.054  22.012  1.00 10.12           H
HETATM  933  H2  HOH A  63      14.452  26.173  22.407  1.00 10.12           H
'''

def run_and_get_clashscore(cmd, out):
  rc = easy_run.fully_buffered(cmd)
  with open(out, 'r') as fo:
    for line in fo.readlines():
      if line.find('clashscore =')>-1:
        return float(line.split()[-1])
      if line.find('Clashscore            =')>-1:
        return float(line.split()[-1])
  assert 0

def run():
  with open('tst_mol_args.pdb', 'w') as fo:
    fo.write(pdb_str)
  for keep_hydrogens in [True, False]:
    for nuclear in [True, False]:
      cmd1 = " ".join([
        "phenix.clashscore tst_mol_args.pdb",
        "keep_hydrogens=%s"%str(keep_hydrogens),
        "nuclear=%s"%str(nuclear),
        ">& tst_molprobity_arguments.zlog1"])
      cmd2 = " ".join([
        "phenix.molprobity tst_mol_args.pdb",
        "keep_hydrogens=%s"%str(keep_hydrogens),
        "pdb_interpretation.use_neutron_distances=%s"%str(nuclear),
        ">& tst_molprobity_arguments.zlog2"])
      cs1=run_and_get_clashscore(cmd=cmd1, out="tst_molprobity_arguments.zlog1")
      cs2=run_and_get_clashscore(cmd=cmd2, out="tst_molprobity_arguments.zlog2")
      if(0):
        print()
        print(cmd1)
        print(cmd2)
        print(cs1, cs2)
        print()
      assert approx_equal(cs1, cs2, 1.e-3)

if __name__=='__main__':
  run()
