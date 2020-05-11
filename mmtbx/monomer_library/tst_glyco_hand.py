from __future__ import absolute_import, division, print_function
from mmtbx.monomer_library.glyco_chiral_values import alpha_beta

def main():
  print('codes',len(alpha_beta))
  tests = {"NAG" : "beta",
           'NDG' : 'alpha',
           'MAN' : "alpha",
           'BMA' : "beta",
           "FUC" : "alpha",
           "FUL" : "beta",
           'BDP' : 'beta',
           'GCU' : 'alpha',
           }
  for code in tests:
    print("  %-3s %-8s %s" % (code, tests[code], alpha_beta[code]))
    assert tests[code]==alpha_beta[code]

if __name__ == '__main__':
  main()
