from __future__ import absolute_import, division, print_function
import time
import mmtbx.model
import iotbx.pdb
from libtbx.utils import null_out
from six.moves import zip

#-----------------------------------------------------------------------------
# This test checks the parameterization of hydrogen atoms for planar X-H2 groups
# Aim is to test if ideal dihedral angles are correct for any configuration
#-----------------------------------------------------------------------------

def exercise(pdb_str):
  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  model = mmtbx.model.manager(
    model_input = pdb_inp,
    log         = null_out())
  model.process(make_restraints=True)
  pdb_hierarchy = model.get_hierarchy()
  sites_cart = model.get_sites_cart()
  atoms = pdb_hierarchy.atoms()

  model.setup_riding_h_manager()
  riding_h_manager = model.get_riding_h_manager()

  h_parameterization = riding_h_manager.h_parameterization

  diagnostics = riding_h_manager.diagnostics(
    sites_cart         = sites_cart,
    threshold          = 0.05)

  diagnostics = riding_h_manager.diagnostics(
    sites_cart         = sites_cart,
    threshold          = 0.05)
  h_distances   = diagnostics.h_distances

  number_h = model.get_hd_selection().count(True)
  number_h_para = len(h_parameterization) - h_parameterization.count(None)

# There are 90 H atoms in pdb_string, check if all of them are recognized
  assert (number_h_para == number_h), 'Not all H atoms are parameterized'

  for ih in h_distances:
    labels = atoms[ih].fetch_labels()
    assert (h_distances[ih] < 0.01), \
      'distance too large: %s  atom: %s (%s) residue: %s ' \
      % (h_parameterization[ih].htype, atoms[ih].name, ih, labels.resseq.strip())


pdb_str_00 = """
CRYST1   13.889   13.738   13.126  90.00  90.00  90.00 P 1
SCALE1      0.071999  0.000000  0.000000        0.00000
SCALE2      0.000000  0.072791  0.000000        0.00000
SCALE3      0.000000  0.000000  0.076185        0.00000
ATOM      1  N   ASN A   1       5.264   6.323   6.229  1.00 10.00           N
ATOM      2  CA  ASN A   1       6.108   7.341   6.843  1.00 10.00           C
ATOM      3  C   ASN A   1       5.944   8.680   6.132  1.00 10.00           C
ATOM      4  O   ASN A   1       5.516   8.733   4.979  1.00 10.00           O
ATOM      5  CB  ASN A   1       7.575   6.907   6.817  1.00 10.00           C
ATOM      6  CG  ASN A   1       7.857   5.750   7.754  1.00 10.00           C
ATOM      7  OD1 ASN A   1       7.361   5.713   8.880  1.00 10.00           O
ATOM      8  ND2 ASN A   1       8.659   4.797   7.294  1.00 10.00           N
ATOM      9  HA  ASN A   1       5.844   7.457   7.769  1.00 10.00           H
ATOM     10  HB2 ASN A   1       7.806   6.628   5.917  1.00 10.00           H
ATOM     11  HB3 ASN A   1       8.131   7.655   7.086  1.00 10.00           H
ATOM     12 HD21 ASN A   1       8.849   4.119   7.787  1.00 10.00           H
ATOM     13 HD22 ASN A   1       8.988   4.857   6.501  1.00 10.00           H
"""

pdb_str_01 = """
CRYST1   13.889   13.738   13.126  90.00  90.00  90.00 P 1
SCALE1      0.071999  0.000000  0.000000        0.00000
SCALE2      0.000000  0.072791  0.000000        0.00000
SCALE3      0.000000  0.000000  0.076185        0.00000
ATOM      1  N   ASN A   1       5.264   6.323   6.229  1.00 10.00           N
ATOM      2  CA  ASN A   1       6.108   7.341   6.843  1.00 10.00           C
ATOM      3  C   ASN A   1       5.944   8.680   6.132  1.00 10.00           C
ATOM      4  O   ASN A   1       5.516   8.733   4.979  1.00 10.00           O
ATOM      5  CB  ASN A   1       7.575   6.907   6.817  1.00 10.00           C
ATOM      6  CG  ASN A   1       7.857   5.750   7.754  1.00 10.00           C
ATOM      7  OD1 ASN A   1       7.361   5.713   8.880  1.00 10.00           O
ATOM      8  ND2 ASN A   1       8.659   4.797   7.294  1.00 10.00           N
ATOM      9  HA  ASN A   1       5.844   7.457   7.769  1.00 10.00           H
ATOM     10  HB2 ASN A   1       7.806   6.628   5.917  1.00 10.00           H
ATOM     11  HB3 ASN A   1       8.131   7.655   7.086  1.00 10.00           H
ATOM     12 HD21 ASN A   1       8.988   4.857   6.501  1.00 10.00           H
ATOM     13 HD22 ASN A   1       8.849   4.119   7.787  1.00 10.00           H
"""

pdb_str_02 = """
CRYST1   14.446   16.451   11.913  90.00  90.00  90.00 P 1
SCALE1      0.069223  0.000000  0.000000        0.00000
SCALE2      0.000000  0.060787  0.000000        0.00000
SCALE3      0.000000  0.000000  0.083942        0.00000
ATOM      1  N   ARG A  50       6.725  11.631   6.607  1.00 10.00           N
ATOM      2  CA  ARG A  50       7.490  10.491   6.117  1.00 10.00           C
ATOM      3  C   ARG A  50       8.914  10.908   5.763  1.00 10.00           C
ATOM      4  O   ARG A  50       9.701  11.267   6.639  1.00 10.00           O
ATOM      5  CB  ARG A  50       6.806   9.870   4.897  1.00 10.00           C
ATOM      6  CG  ARG A  50       5.402   9.344   5.166  1.00 10.00           C
ATOM      7  CD  ARG A  50       5.419   8.122   6.072  1.00 10.00           C
ATOM      8  NE  ARG A  50       6.111   6.992   5.457  1.00 10.00           N
ATOM      9  CZ  ARG A  50       6.529   5.916   6.119  1.00 10.00           C
ATOM     10  NH1 ARG A  50       6.327   5.812   7.427  1.00 10.00           N
ATOM     11  NH2 ARG A  50       7.150   4.940   5.472  1.00 10.00           N
ATOM     12  HA  ARG A  50       7.540   9.818   6.814  1.00 10.00           H
ATOM     13  HB2 ARG A  50       6.740  10.543   4.202  1.00 10.00           H
ATOM     14  HB3 ARG A  50       7.345   9.127   4.582  1.00 10.00           H
ATOM     15  HG2 ARG A  50       4.880  10.036   5.600  1.00 10.00           H
ATOM     16  HG3 ARG A  50       4.990   9.092   4.325  1.00 10.00           H
ATOM     17  HD2 ARG A  50       5.872   8.343   6.900  1.00 10.00           H
ATOM     18  HD3 ARG A  50       4.505   7.851   6.254  1.00 10.00           H
ATOM     19  HE  ARG A  50       6.259   7.024   4.611  1.00 10.00           H
ATOM     20 HH11 ARG A  50       5.925   6.440   7.855  1.00 10.00           H
ATOM     21 HH12 ARG A  50       6.600   5.114   7.847  1.00 10.00           H
ATOM     22 HH21 ARG A  50       7.283   5.003   4.624  1.00 10.00           H
ATOM     23 HH22 ARG A  50       7.420   4.245   5.899  1.00 10.00           H
"""

pdb_str_03 = """
CRYST1   14.446   16.451   11.913  90.00  90.00  90.00 P 1
SCALE1      0.069223  0.000000  0.000000        0.00000
SCALE2      0.000000  0.060787  0.000000        0.00000
SCALE3      0.000000  0.000000  0.083942        0.00000
ATOM      1  N   ARG A  50       6.725  11.631   6.607  1.00 10.00           N
ATOM      2  CA  ARG A  50       7.490  10.491   6.117  1.00 10.00           C
ATOM      3  C   ARG A  50       8.914  10.908   5.763  1.00 10.00           C
ATOM      4  O   ARG A  50       9.701  11.267   6.639  1.00 10.00           O
ATOM      5  CB  ARG A  50       6.806   9.870   4.897  1.00 10.00           C
ATOM      6  CG  ARG A  50       5.402   9.344   5.166  1.00 10.00           C
ATOM      7  CD  ARG A  50       5.419   8.122   6.072  1.00 10.00           C
ATOM      8  NE  ARG A  50       6.111   6.992   5.457  1.00 10.00           N
ATOM      9  CZ  ARG A  50       6.529   5.916   6.119  1.00 10.00           C
ATOM     10  NH1 ARG A  50       6.327   5.812   7.427  1.00 10.00           N
ATOM     11  NH2 ARG A  50       7.150   4.940   5.472  1.00 10.00           N
ATOM     12  HA  ARG A  50       7.540   9.818   6.814  1.00 10.00           H
ATOM     13  HB2 ARG A  50       6.740  10.543   4.202  1.00 10.00           H
ATOM     14  HB3 ARG A  50       7.345   9.127   4.582  1.00 10.00           H
ATOM     15  HG2 ARG A  50       4.880  10.036   5.600  1.00 10.00           H
ATOM     16  HG3 ARG A  50       4.990   9.092   4.325  1.00 10.00           H
ATOM     17  HD2 ARG A  50       5.872   8.343   6.900  1.00 10.00           H
ATOM     18  HD3 ARG A  50       4.505   7.851   6.254  1.00 10.00           H
ATOM     19  HE  ARG A  50       6.259   7.024   4.611  1.00 10.00           H
ATOM     21 HH11 ARG A  50       6.600   5.114   7.847  1.00 10.00           H
ATOM     20 HH12 ARG A  50       5.925   6.440   7.855  1.00 10.00           H
ATOM     23 HH21 ARG A  50       7.420   4.245   5.899  1.00 10.00           H
ATOM     22 HH22 ARG A  50       7.283   5.003   4.624  1.00 10.00           H
"""

pdb_str_04 = """
CRYST1   14.547   13.375   15.374  90.00  90.00  90.00 P 1
SCALE1      0.068743  0.000000  0.000000        0.00000
SCALE2      0.000000  0.074766  0.000000        0.00000
SCALE3      0.000000  0.000000  0.065045        0.00000
ATOM      1  N   GLN A  70       6.510   6.937   4.976  1.00 10.00           N
ATOM      2  CA  GLN A  70       6.727   6.312   6.275  1.00 10.00           C
ATOM      3  C   GLN A  70       5.406   6.141   7.018  1.00 10.00           C
ATOM      4  O   GLN A  70       4.992   5.022   7.317  1.00 10.00           O
ATOM      5  CB  GLN A  70       7.700   7.144   7.113  1.00 10.00           C
ATOM      6  CG  GLN A  70       8.079   6.506   8.440  1.00 10.00           C
ATOM      7  CD  GLN A  70       9.104   7.319   9.206  1.00 10.00           C
ATOM      8  OE1 GLN A  70       9.517   8.392   8.767  1.00 10.00           O
ATOM      9  NE2 GLN A  70       9.521   6.810  10.360  1.00 10.00           N
ATOM     10  HA  GLN A  70       7.116   5.433   6.144  1.00 10.00           H
ATOM     11  HB2 GLN A  70       8.516   7.275   6.605  1.00 10.00           H
ATOM     12  HB3 GLN A  70       7.292   8.003   7.304  1.00 10.00           H
ATOM     13  HG2 GLN A  70       7.287   6.427   8.993  1.00 10.00           H
ATOM     14  HG3 GLN A  70       8.456   5.628   8.273  1.00 10.00           H
ATOM     15 HE21 GLN A  70       9.210   6.057  10.635  1.00 10.00           H
ATOM     16 HE22 GLN A  70      10.103   7.233  10.832  1.00 10.00           H
"""

pdb_str_05 = """
CRYST1   14.547   13.375   15.374  90.00  90.00  90.00 P 1
SCALE1      0.068743  0.000000  0.000000        0.00000
SCALE2      0.000000  0.074766  0.000000        0.00000
SCALE3      0.000000  0.000000  0.065045        0.00000
ATOM      1  N   GLN A  70       6.510   6.937   4.976  1.00 10.00           N
ATOM      2  CA  GLN A  70       6.727   6.312   6.275  1.00 10.00           C
ATOM      3  C   GLN A  70       5.406   6.141   7.018  1.00 10.00           C
ATOM      4  O   GLN A  70       4.992   5.022   7.317  1.00 10.00           O
ATOM      5  CB  GLN A  70       7.700   7.144   7.113  1.00 10.00           C
ATOM      6  CG  GLN A  70       8.079   6.506   8.440  1.00 10.00           C
ATOM      7  CD  GLN A  70       9.104   7.319   9.206  1.00 10.00           C
ATOM      8  OE1 GLN A  70       9.517   8.392   8.767  1.00 10.00           O
ATOM      9  NE2 GLN A  70       9.521   6.810  10.360  1.00 10.00           N
ATOM     10  HA  GLN A  70       7.116   5.433   6.144  1.00 10.00           H
ATOM     11  HB2 GLN A  70       8.516   7.275   6.605  1.00 10.00           H
ATOM     12  HB3 GLN A  70       7.292   8.003   7.304  1.00 10.00           H
ATOM     13  HG2 GLN A  70       7.287   6.427   8.993  1.00 10.00           H
ATOM     14  HG3 GLN A  70       8.456   5.628   8.273  1.00 10.00           H
ATOM     16 HE21 GLN A  70      10.103   7.233  10.832  1.00 10.00           H
ATOM     15 HE22 GLN A  70       9.210   6.057  10.635  1.00 10.00           H
"""

pdb_list = [pdb_str_00, pdb_str_01, pdb_str_02, pdb_str_03,
  pdb_str_04, pdb_str_05]

pdb_list_name = ['pdb_str_00', 'pdb_str_01', 'pdb_str_02', 'pdb_str_03',
  'pdb_str_04', 'pdb_str_05']

def run():
  for pdb_str, str_name in zip(pdb_list,pdb_list_name):
    #print 'pdb_string:', str_name, 'use_ideal_bonds_angles =', True
    exercise(pdb_str=pdb_str)

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("OK. Time: %8.3f"%(time.time()-t0))
