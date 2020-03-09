from __future__ import absolute_import, division, print_function
import sys

from iotbx import pdb
from six.moves import range

pdb_lines = '''
ATOM  18755  N   LYS N  56      48.225  62.706 -53.862  1.00 69.77           N
ATOM  18756  CA  LYS N  56      47.122  61.816 -53.629  1.00 70.02           C
ATOM  18757  C   LYS N  56      46.793  61.775 -52.151  1.00 70.49           C
ATOM  18758  O   LYS N  56      45.732  62.272 -51.762  1.00 70.82           O
ATOM  18759  CB  LYS N  56      47.445  60.405 -54.119  1.00 69.96           C
ATOM  18760  CG  LYS N  56      47.862  60.366 -55.558  1.00 69.09           C
ATOM  18761  CD  LYS N  56      46.738  60.829 -56.449  1.00 69.14           C
ATOM  18762  CE  LYS N  56      47.154  60.761 -57.900  1.00 68.65           C
ATOM  18763  NZ  LYS N  56      46.070  61.230 -58.799  1.00 67.64           N
ATOM  18764  N   LYS N  57      47.630  61.207 -51.295  1.00 70.59           N
ATOM  18765  CA  LYS N  57      47.259  61.174 -49.888  1.00 70.64           C
ATOM  18766  C   LYS N  57      48.216  61.977 -49.013  1.00 71.13           C
ATOM  18767  O   LYS N  57      47.781  62.990 -48.461  1.00 71.98           O
ATOM  18768  CB  LYS N  57      47.055  59.726 -49.452  1.00 70.28           C
ATOM  18769  CG  LYS N  57      45.619  59.278 -49.768  1.00 70.08           C
ATOM  18770  CD  LYS N  57      45.350  58.763 -51.161  1.00 68.27           C
ATOM  18771  CE  LYS N  57      43.872  58.439 -51.282  1.00 67.38           C
ATOM  18772  NZ  LYS N  57      43.719  57.886 -52.654  1.00 66.40           N
ATOM  18773  N   VAL N  58      49.522  61.811 -49.121  1.00 71.03           N
ATOM  18774  CA  VAL N  58      50.479  60.749 -49.246  1.00 70.78           C
ATOM  18775  C   VAL N  58      50.801  60.494 -47.743  1.00 70.09           C
ATOM  18776  O   VAL N  58      50.009  60.889 -46.880  1.00 70.44           O
ATOM  18777  CB  VAL N  58      51.665  61.148 -50.115  1.00 71.15           C
ATOM  18778  CG1 VAL N  58      51.233  61.264 -51.566  1.00 71.41           C
ATOM  18779  CG2 VAL N  58      52.289  62.461 -49.638  1.00 71.37           C
ATOM  18780  N   PRO N  59      51.836  59.722 -47.347  1.00 69.26           N
ATOM  18781  CA  PRO N  59      53.160  59.201 -46.949  1.00 68.80           C
ATOM  18782  C   PRO N  59      54.389  59.494 -47.837  1.00 68.33           C
ATOM  18783  O   PRO N  59      54.177  59.842 -48.977  1.00 68.60           O
ATOM  18784  CB  PRO N  59      52.909  57.680 -46.892  1.00 68.81           C
ATOM  18785  CG  PRO N  59      51.511  57.540 -46.598  1.00 69.03           C
ATOM  18786  CD  PRO N  59      50.868  58.605 -47.385  1.00 69.11           C
ATOM  18787  N   PHE N  60      55.619  59.124 -47.487  1.00 67.22           N
ATOM  18788  CA  PHE N  60      56.155  58.824 -46.158  1.00 65.75           C
ATOM  18789  C   PHE N  60      57.465  59.615 -46.163  1.00 65.03           C
ATOM  18790  O   PHE N  60      58.327  59.487 -45.292  1.00 64.85           O
ATOM  18791  CB  PHE N  60      56.325  57.317 -45.873  1.00 65.65           C
ATOM  18792  CG  PHE N  60      57.281  56.587 -46.787  1.00 65.73           C
ATOM  18793  CD1 PHE N  60      56.849  56.085 -48.007  1.00 66.68           C
ATOM  18794  CD2 PHE N  60      58.580  56.314 -46.380  1.00 66.05           C
ATOM  18795  CE1 PHE N  60      57.713  55.378 -48.834  1.00 67.64           C
ATOM  18796  CE2 PHE N  60      59.451  55.605 -47.199  1.00 67.13           C
ATOM  18797  CZ  PHE N  60      59.015  55.134 -48.427  1.00 67.63           C
ATOM  18798  N   LEU N  61      57.540  60.450 -47.206  1.00 64.08           N
ATOM  18799  CA  LEU N  61      58.603  61.375 -47.593  1.00 63.10           C
ATOM  18800  C   LEU N  61      57.880  62.245 -48.648  1.00 62.48           C
ATOM  18801  O   LEU N  61      56.647  62.255 -48.570  1.00 62.00           O
ATOM  18802  CB  LEU N  61      59.905  60.590 -47.871  1.00 62.99           C
ATOM  18803  CG  LEU N  61      60.257  59.711 -49.072  1.00 62.57           C
ATOM  18804  CD1 LEU N  61      61.732  59.308 -49.013  1.00 62.70           C
ATOM  18805  CD2 LEU N  61      59.353  58.526 -49.229  1.00 61.75           C
ATOM  18806  N   SER N  62      58.466  63.043 -49.565  1.00 61.91           N
ATOM  18807  CA  SER N  62      59.847  63.261 -50.010  1.00 61.11           C
ATOM  18808  C   SER N  62      60.085  64.721 -50.329  1.00 60.49           C
ATOM  18809  O   SER N  62      59.168  65.539 -50.236  1.00 60.11           O
ATOM  18810  CB  SER N  62      60.128  62.434 -51.262  1.00 60.97           C
ATOM  18811  OG  SER N  62      61.494  62.508 -51.606  1.00 61.00           O
'''

def run():
  filename = 'tst_multi.pdb'
  f=open(filename, 'w')
  f.write(pdb_lines)
  f.close()
  pdb_inp = pdb.input(filename)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  from mmtbx.conformation_dependent_library.testing_utils import \
      get_geometry_restraints_manager
  geometry_restraints_manager = get_geometry_restraints_manager(filename)
  pdb_hierarchy.reset_i_seq_if_necessary()
  refine = [False, # -179
            True,  # -44
            False, # 86
            True,  # -22
            False, # -179
            ]
  refine += [True]*5
  refine += [False]*14
  omegalyze = [False,
               False,
               False,
               True,
               False,
               ]
  omegalyze += [True]*3
  omegalyze += [False]*16
  from mmtbx.conformation_dependent_library import generate_protein_threes
  for i, threes in enumerate(generate_protein_threes(pdb_hierarchy,
                                        geometry_restraints_manager,
                                        cdl_class=True,
                                        #verbose=verbose,
                                        )):
    print(i, threes)
    print('  omega   %5.1f' % threes.get_omega_value())
    print("  cis?    %-5s %s" % (threes.cis_group(), threes.cis_group(limit=30)))
    print("  trans?  %-5s %s" % (threes.trans_group(), threes.trans_group(limit=30)))
    print("  rama    %s" % threes.get_ramalyze_key())
    print('  conf    %s' % threes.is_pure_main_conf())
    assert threes.cis_group()==refine[i], '%s!=%s' % (threes.cis_group(), refine[i])
    assert threes.cis_group(limit=30)==omegalyze[i]

  for j in range(0,181,10):
    i+=1
    print("  %3d %-5s %-8s %-5s" % (
      j,
      threes._define_omega_a_la_duke_using_limit(j)=='cis',
      threes._define_omega_a_la_duke_using_limit(j, limit=30),
      refine[i],
    ))
    assert (threes._define_omega_a_la_duke_using_limit(j)=='cis')==refine[i]



if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run()
