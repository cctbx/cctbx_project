from __future__ import absolute_import, division, print_function
import iotbx.pdb

pdb_str1 = """
ATOM      1  N  AVAL A   1      -4.898   0.072  13.387  0.49  7.34           N
ATOM      2  CA AVAL A   1      -4.626   0.703  12.080  0.49  7.71           C
ATOM      3  C  AVAL A   1      -3.475   1.680  12.227  0.49  7.52           C
ATOM      4  O  AVAL A   1      -3.125   2.100  13.335  0.49  7.59           O
ATOM      5  CB AVAL A   1      -5.882   1.390  11.495  0.49  7.79           C
ATOM      6  CG1AVAL A   1      -6.979   0.367  11.258  0.49  7.20           C
ATOM      7  CG2AVAL A   1      -6.359   2.496  12.408  0.49  7.81           C
ATOM      8  H1 AVAL A   1      -4.794  -0.809  13.318  0.49  8.81           H
ATOM      9  H2 AVAL A   1      -4.331   0.391  13.995  0.49  8.81           H
ATOM     10  H3 AVAL A   1      -5.733   0.254  13.636  0.49  8.81           H
ATOM     11  HA AVAL A   1      -4.371   0.017  11.444  0.49  9.26           H
ATOM     12  HB AVAL A   1      -5.655   1.790  10.641  0.49  9.35           H
HETATM 2211  O   HOH S 216       9.105  -6.647  -4.343  0.25 56.37           O
HETATM 2219  O  BHOH S 224       6.977   3.045   9.044  0.31 55.60           O
"""

pdb_str2 = """
HETATM 2174  O   HOH S 179      -5.781   7.569   9.276  0.24 14.70           O
HETATM 2174  O  AHOH S 179      -5.781   7.569   9.276  0.24 14.70           O
HETATM 2174  O   HOH S 179      -4.781   7.569   9.276  0.24 14.70           H
HETATM 2174  O  CHOH S 179      -6.781   7.569   9.276  0.24 14.70           H
END
"""

pdb_str3 = """
HETATM 2174  O   HOH S 179      -5.781   7.569   9.276  0.24 14.70           O
HETATM 2174  O   HOH S 179      -5.781   7.569   9.276  0.24 14.70           O
HETATM 2174  O   HOH S 179      -4.781   7.569   9.276  0.24 14.70           H
HETATM 2174  O   HOH S 179      -6.781   7.569   9.276  0.24 14.70           H
END
"""

def get_h(s):
  return iotbx.pdb.input(source_info=None, lines=s).construct_hierarchy()

def run():
  #
  # This is current state of affairs. This was also the case in
  # Phenix 1.20 and 1.18.
  #
  altlocs = [ag.altloc for ag in get_h(pdb_str1).only_model().atom_groups()]
  assert altlocs == ['A', '', 'B'], altlocs
  #
  print("------------")
  altlocs = [ag.altloc for ag in get_h(pdb_str2).only_model().atom_groups()]
  assert altlocs == [' ', 'A', 'C'], altlocs

  print("------------")
  altlocs = [ag.altloc for ag in get_h(pdb_str3).only_model().atom_groups()]
  assert altlocs == [''], altlocs

if (__name__ == "__main__"):
  run()
  print("OK")
