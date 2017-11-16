from __future__ import division
import time
import mmtbx.model
import iotbx.pdb
from libtbx.utils import null_out
from mmtbx.monomer_library.pdb_interpretation import grand_master_phil_str
from libtbx.test_utils import approx_equal
from mmtbx.hydrogens.validate_H import validate_H
#from validate_H_cl_app import master_params_str

pdb_str = """
CRYST1   31.264   27.900   96.292  90.00  90.00  90.00 P 1
SCALE1      0.031986  0.000000  0.000000        0.00000
SCALE2      0.000000  0.035842  0.000000        0.00000
SCALE3      0.000000  0.000000  0.010385        0.00000
ATOM      1  N   ASN A   1      25.174 -24.761 -21.940  1.00 16.59           N
ATOM      2  CA  ASN A   1      25.971 -25.578 -22.850  1.00 17.25           C
ATOM      3  C   ASN A   1      26.037 -24.941 -24.225  1.00 16.66           C
ATOM      4  O   ASN A   1      27.077 -24.979 -24.874  1.00 18.42           O
ATOM      5  CB  ASN A   1      25.498 -27.040 -22.913  1.00 18.30           C
ATOM      6  CG  ASN A   1      24.187 -27.208 -23.672  1.00 20.42           C
ATOM      7  OD1 ASN A   1      23.167 -26.590 -23.335  1.00 20.87           O
ATOM      8  ND2 ASN A   1      24.200 -28.069 -24.685  1.00 18.15           N
ATOM      9  DA  ASN A   1      26.975 -25.581 -22.461  1.00 15.99           D
ATOM     10 DD21 ASN A   1      25.030 -28.545 -24.882  1.00 16.07           D
ATOM     11  D  AASN A   1      24.379 -25.105 -21.517  0.82 16.71           D
ATOM     12  DB2AASN A   1      25.358 -27.400 -21.898  0.65 18.82           D
ATOM     13  DB3AASN A   1      26.255 -27.641 -23.402  0.56 18.30           D
ATOM     14 DD22AASN A   1      23.370 -28.163 -25.192  0.72 16.98           D
ATOM     15  H  BASN A   1      24.379 -25.105 -21.517  0.18 16.71           H
ATOM     16  HB2BASN A   1      26.255 -27.641 -23.402  0.44 18.30           H
ATOM     17  HB3BASN A   1      25.358 -27.400 -21.898  0.35 18.82           H
ATOM     18 HD22BASN A   1      23.370 -28.163 -25.192  0.28 16.98           H
ATOM     19  N   ASN A   2      23.615 -22.281 -25.492  1.00 16.59           N
ATOM     20  CA  ASN A   2      24.412 -23.098 -26.402  1.00 17.25           C
ATOM     21  C   ASN A   2      24.478 -22.461 -27.777  1.00 16.66           C
ATOM     22  O   ASN A   2      25.518 -22.499 -28.426  1.00 18.42           O
ATOM     23  CB  ASN A   2      23.939 -24.560 -26.465  1.00 18.30           C
ATOM     24  CG  ASN A   2      22.628 -24.728 -27.224  1.00 20.42           C
ATOM     25  OD1 ASN A   2      21.608 -24.110 -26.887  1.00 20.87           O
ATOM     26  ND2 ASN A   2      22.641 -25.589 -28.237  1.00 18.15           N
ATOM     27  DA  ASN A   2      25.416 -23.101 -26.013  1.00 15.99           D
ATOM     28 HD21 ASN A   2      23.471 -26.065 -28.434  1.00 16.07           H
ATOM     29  D  AASN A   2      22.820 -22.625 -25.069  0.82 16.71           D
ATOM     30  DB2AASN A   2      23.799 -24.920 -25.450  0.65 18.82           D
ATOM     31  DB3AASN A   2      24.696 -25.161 -26.954  0.56 18.30           D
ATOM     32 DD22AASN A   2      21.811 -25.683 -28.744  0.72 16.98           D
ATOM     33  H  BASN A   2      22.820 -22.625 -25.069  0.18 16.71           H
ATOM     34  HB2BASN A   2      24.696 -25.161 -26.954  0.44 18.30           H
ATOM     35  HB3BASN A   2      23.799 -24.920 -25.450  0.35 18.82           H
ATOM     36 HD22BASN A   2      21.811 -25.683 -28.744  0.28 16.98           H
ATOM     37  N   LEU A   3      23.374 -21.865 -28.222  1.00 15.90           N
ATOM     38  CA  LEU A   3      23.342 -21.224 -29.530  1.00 15.47           C
ATOM     39  C   LEU A   3      23.637 -19.727 -29.491  1.00 18.18           C
ATOM     40  O   LEU A   3      24.174 -19.175 -30.462  1.00 20.24           O
ATOM     41  CB  LEU A   3      22.021 -21.521 -30.238  1.00 15.76           C
ATOM     42  CG  LEU A   3      21.790 -23.005 -30.533  1.00 19.30           C
ATOM     43  CD1 LEU A   3      20.460 -23.197 -31.290  1.00 20.74           C
ATOM     44  CD2 LEU A   3      22.957 -23.627 -31.336  1.00 17.02           C
ATOM     45  DA  LEU A   3      24.132 -21.653 -30.101  1.00 12.76           D
ATOM     46  DB3 LEU A   3      21.233 -21.196 -29.610  1.00 17.93           D
ATOM     47  D  ALEU A   3      22.567 -21.865 -27.663  0.91 16.63           D
ATOM     48  DB2ALEU A   3      21.964 -20.976 -31.155  0.67 16.51           D
ATOM     49  DG ALEU A   3      21.719 -23.530 -29.581  0.68 18.31           D
ATOM     50 DD11ALEU A   3      19.620 -22.925 -30.652  0.74 20.52           D
ATOM     51 DD12ALEU A   3      20.364 -24.239 -31.597  0.77 19.78           D
ATOM     52 DD13ALEU A   3      20.458 -22.568 -32.171  0.49 20.26           D
ATOM     53 DD21ALEU A   3      23.466 -22.855 -31.890  0.40 16.87           D
ATOM     54 DD22ALEU A   3      22.563 -24.357 -32.028  0.00 16.68           D
ATOM     55 DD23ALEU A   3      23.653 -24.109 -30.673  0.89 14.39           D
ATOM     56  H  BLEU A   3      22.567 -21.865 -27.663  0.09 16.63           H
ATOM     57  HB2BLEU A   3      21.964 -20.976 -31.155  0.33 14.51           H
ATOM     58  HG BLEU A   3      21.719 -23.530 -29.581  0.32 18.31           H
ATOM     59 HD11BLEU A   3      19.620 -22.925 -30.652  0.26 20.52           H
ATOM     60 HD12BLEU A   3      20.364 -24.239 -31.597  0.23 19.78           H
ATOM     61 HD13BLEU A   3      20.458 -22.568 -32.171  0.51 20.26           H
ATOM     62 HD21BLEU A   3      23.466 -22.855 -31.890  0.60 16.87           H
ATOM     63 HD22BLEU A   3      22.563 -24.357 -32.028  0.58 16.68           H
ATOM     64 HD23BLEU A   3      23.653 -24.109 -30.673  0.11 14.39           H
ATOM     65  N   SER A   4       9.243 -32.509 -36.471  1.00 38.02           N
ATOM     66  CA  SER A   4      10.063 -33.604 -36.981  1.00 39.55           C
ATOM     67  C   SER A   4      11.459 -33.153 -37.367  1.00 37.45           C
ATOM     68  O   SER A   4      11.633 -32.304 -38.239  1.00 38.82           O
ATOM     69  CB  SER A   4       9.384 -34.279 -38.180  1.00 41.66           C
ATOM     70  OG  SER A   4       8.375 -35.181 -37.753  1.00 44.91           O
ATOM     71  DA  SER A   4      10.165 -34.324 -36.187  1.00 38.72           D
ATOM     72  D  ASER A   4       8.730 -31.950 -37.103  0.62 38.48           D
ATOM     73  DB2ASER A   4      10.116 -34.818 -38.774  0.63 41.36           D
ATOM     74  DB3ASER A   4       8.918 -33.518 -38.793  0.36 41.73           D
ATOM     75  DG ASER A   4       8.699 -36.084 -37.842  0.36 43.61           D
ATOM     76  H  BSER A   4       8.730 -31.950 -37.103  0.28 38.48           H
ATOM     77  HB2BSER A   4       8.918 -33.518 -38.793  0.64 41.73           H
ATOM     78  HB3BSER A   4      10.116 -34.818 -38.774  0.37 41.36           H
ATOM     79  HG BSER A   4       8.699 -36.084 -37.842  0.64 43.61           H
ATOM     80  N   SER A   5      26.478 -19.398  45.776  1.00 19.78           N
ATOM     81  C   SER A   5      27.940 -19.562  43.818  1.00 21.99           C
ATOM     82  O   SER A   5      29.005 -19.770  43.306  1.00 21.20           O
ATOM     83  CA ASER A   5      27.849 -19.441  45.318  0.64 18.99           C
ATOM     84  CB ASER A   5      28.611 -20.590  45.987  0.64 21.66           C
ATOM     85  OG ASER A   5      28.819 -20.337  47.363  0.64 23.21           O
ATOM     86  DG ASER A   5      29.639 -20.155  47.499  0.64 22.10           D
ATOM     87  H  ASER A   5      26.064 -20.150  45.742  0.06 18.55           H
ATOM     88  HA ASER A   5      28.296 -18.610  45.561  0.64 20.89           H
ATOM     89  HB2ASER A   5      28.094 -21.408  45.903  0.00 22.17           H
ATOM     90  HB3ASER A   5      29.468 -20.686  45.544  0.64 22.38           H
ATOM     91  CA BSER A   5      27.851 -19.476  45.322  0.36 19.02           C
ATOM     92  CB BSER A   5      28.561 -20.680  45.932  0.36 21.98           C
ATOM     93  OG BSER A   5      28.005 -21.879  45.425  0.36 25.40           O
ATOM     94  D  BSER A   5      26.023 -20.128  45.688  0.94 18.55           D
ATOM     95  HA BSER A   5      28.327 -18.675  45.599  0.36 21.03           H
ATOM     96  HB2BSER A   5      29.498 -20.638  45.691  0.36 22.55           H
ATOM     97  HB3BSER A   5      28.462 -20.665  46.898  0.36 23.83           H
ATOM     98  HG BSER A   5      27.486 -22.196  46.000  0.36 22.10           H
TER
HETATM   99  O   HOH B   1      25.202 -21.329 -23.306  1.00 30.00           O
HETATM  100  D1  HOH B   1      25.963 -20.763 -23.255  1.00 30.00           D
HETATM  101  D2  HOH B   1      24.439 -20.764 -23.305  1.00 30.00           D
HETATM  102  O   HOH B   2      26.442 -21.454 -25.349  1.00 30.00           O
HETATM  103  D1  HOH B   2      27.249 -20.953 -25.339  1.00 30.00           D
HETATM  104  O   HOH B   3      24.970 -18.184 -24.460  1.00 30.00           O
HETATM  105  O   HOH B   4      21.573 -23.174 -24.023  1.00 30.00           O
HETATM  106  D2  HOH B   4      20.731 -22.736 -24.071  1.00 30.00           D
HETATM  107  D1 AHOH B   4      22.045 -22.784 -23.297  1.00 30.00           D
HETATM  108  D1 BHOH B   4      22.145 -21.784 -22.297  1.00 30.00           D
HETATM  109  D1  HOH B   5      21.772 -27.767 -28.880  1.00 30.00           D
TER
END
"""

def exercise():
  log = null_out()
  pdb_interpretation_phil = iotbx.phil.parse(
    input_string = grand_master_phil_str, process_includes = True)
  pi_params = pdb_interpretation_phil.extract()
  pi_params.pdb_interpretation.use_neutron_distances = True

  #validate_H_params = iotbx.phil.parse(master_params_str, process_includes=True)

  pdb_inp = iotbx.pdb.input(lines=pdb_str.split("\n"), source_info=None)
  model = mmtbx.model.manager(
      model_input = pdb_inp,
      process_input = True,
      pdb_interpretation_params = pi_params)

  c = validate_H(model)
  r = c.validate_inputs()
  if r != 0:
    return()
  c.run()
  results = c.get_results()

  oc = results.overall_counts_hd
  hd_atoms_with_occ_0 = oc.hd_atoms_with_occ_0
  assert (oc.count_h == 30)
  assert (oc.count_d == 36)
  assert (oc.count_h_protein == 30)
  assert (oc.count_d_protein == 29)
  assert (oc.count_h_water == 0)
  assert (oc.count_d_water == 7)
  assert (oc.count_water == 5)
  assert (oc.count_water_0h == 1)
  assert (oc.count_water_1h == 1)
  assert (oc.count_water_2h == 1)
  assert (oc.count_water_altconf == 1)
  assert (oc.count_water_no_oxygen == 1)

  answer = ['DD22ALEU A   3 ', ' HB2ASER A   5 ']
  for item, answer in zip(hd_atoms_with_occ_0, answer):
    assert(item[0][5:-1] == answer)

  renamed_answer = [('DB2', 'DB3'), ('DB3', 'DB2'),('DB2', 'DB3'),
    ('DB3', 'DB2'),('DB2', 'DB3'), ('DB3', 'DB2')]
  renamed = results.renamed
  for entry, answer in zip(renamed, renamed_answer):
    oldname = entry['oldname'].strip()
    newname = entry['atom'].name.strip()
    assert(oldname == answer[0])
    assert(newname == answer[1])

  hd_exchanged_sites = results.hd_exchanged_sites
  assert (len(results.hd_exchanged_sites.keys()) == 22)

  hd_sites_analysis  = results.hd_sites_analysis
  sites_different_xyz = hd_sites_analysis.sites_different_xyz
  xyz_answer = [0.0712]
  for item, answer in zip(sites_different_xyz, xyz_answer):
    assert approx_equal(item[2],answer, 1.e-2)

  sites_different_b   = hd_sites_analysis.sites_different_b
  b_answer = [-2]
  for item, answer in zip(sites_different_b, b_answer):
    assert approx_equal(item[2],answer, 1.e-1)

  sites_sum_occ_not_1 = hd_sites_analysis.sites_sum_occ_not_1
  sum_occ_answer=[0.58, 0.9]
  for item, answer in zip(sites_sum_occ_not_1, sum_occ_answer):
    assert approx_equal(item[2], answer, 1.e-2)

  sites_occ_sum_no_scattering = hd_sites_analysis.sites_occ_sum_no_scattering
  sum_scatt_answer=[0.4, 0.36]
  for item, answer in zip(sites_occ_sum_no_scattering, sum_scatt_answer):
    assert approx_equal(item[3], answer, 1.e-2)

  missing_HD_atoms   = results.missing_HD_atoms
  # XXX TODO

if (__name__ == "__main__"):
  t0 = time.time()
  exercise()
  print "OK. Time: %8.3f"%(time.time()-t0)

