
from __future__ import absolute_import, division, print_function
from mmtbx.secondary_structure import sec_str_master_phil_str, manager
from cctbx import geometry_restraints
import iotbx.pdb
import libtbx.load_env
import iotbx.pdb.secondary_structure as ioss
from libtbx.utils import null_out
from libtbx.test_utils import approx_equal
import os

def exercise_protein():
  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/1ywf.pdb",
    test=os.path.isfile)
  pdb_file_h = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/1ywf_h.pdb",
    test=os.path.isfile)
  if pdb_file is None :
    print("Skipping exercise(): input file not available.")
    return False
  if pdb_file_h is None :
    print("Skipping exercise(): input file not available.")
    return False
  run_ksdssp = True
  if (not libtbx.env.has_module(name="ksdssp")):
    print("Skipping KSDSSP tests: ksdssp module not available.")
    run_ksdssp = False
  log = null_out()
  import sys
  log = sys.stdout
  expected_distances = [2.9, 1.975]
  for k, file_name in enumerate([pdb_file, pdb_file_h]):
    pdb_in = iotbx.pdb.input(file_name)
    pdb_hierarchy = pdb_in.construct_hierarchy()
    sec_str_from_pdb_file = pdb_in.extract_secondary_structure()
    m = manager(pdb_hierarchy=pdb_hierarchy,
        sec_str_from_pdb_file=sec_str_from_pdb_file)
    m.params.secondary_structure.protein.remove_outliers = False
    proxies, hb_angles = m.create_protein_hbond_proxies(
        annotation=None, log=log)
    print(proxies.size(), hb_angles.size())
    assert proxies.size() == 109
    assert hb_angles.size() == 294+27, hb_angles.size() # helix+sheet angles
    assert isinstance(proxies, geometry_restraints.shared_bond_simple_proxy)
    (frac_alpha, frac_beta) = m.calculate_structure_content()
    assert approx_equal(frac_alpha, 0.643, eps=1e-3)
    assert approx_equal(frac_beta, 0.0539, eps=1e-3)

    # Make sure the hydrogen auto-detection override is working
    # Functionality is disabled
    # m.params.secondary_structure.protein.substitute_n_for_h = True
    # build_proxies = m.create_hbond_proxies(
    #     log=log,
    #     as_python_objects=True,
    #     as_regular_bond_proxies=False)
    # proxies = build_proxies.proxies
    # atom_ids = []
    # for i_seq in proxies[0].i_seqs :
    #   atom_ids.append(pdb_hierarchy.atoms()[i_seq].id_str())
    # print atom_ids
    # assert (atom_ids == ['pdb=" N   ARG A  41 "', 'pdb=" O   ASP A  37 "'])
    # assert (proxies[0].distance_ideal == 2.9)
    if (run_ksdssp):
      m = manager(pdb_hierarchy=pdb_hierarchy,
          sec_str_from_pdb_file=None)
      m.params.secondary_structure.protein.remove_outliers = False
      proxies, angle_p = m.create_protein_hbond_proxies(
          annotation=None, log=log)
      assert proxies.size() == 90
      assert angle_p.size() == 243+18, angle_p.size()

def exercise_sheet_ends():
  pdb_in = iotbx.pdb.input(source_info=None, lines="""\
SHEET    1   A 3 ALA A  38  ALA A  44  0
SHEET    2   A 3 ALA A  76  ALA A  83 -1  N  ALA A  83   O  ALA A  38
SHEET    3   A 3 ALA A  95  ALA A  97 -1  N  ALA A  97   O  ALA A  76
ATOM      1  N   ALA A  37      27.792  -1.121  13.875  1.00 20.00           N
ATOM      2  CA  ALA A  37      28.151   0.253  14.202  1.00 20.00           C
ATOM      3  C   ALA A  37      26.961   1.003  14.792  1.00 20.00           C
ATOM      4  O   ALA A  37      27.129   2.000  15.494  1.00 20.00           O
ATOM      5  CB  ALA A  37      28.674   0.973  12.969  1.00 20.00           C
ATOM      6  N   ALA A  38      25.758   0.515  14.502  1.00 20.00           N
ATOM      7  CA  ALA A  38      24.529   1.125  14.995  1.00 20.00           C
ATOM      8  C   ALA A  38      23.363   0.264  14.524  1.00 20.00           C
ATOM      9  O   ALA A  38      23.538  -0.687  13.756  1.00 20.00           O
ATOM     10  N   ALA A  39      22.168   0.611  14.996  1.00 20.00           N
ATOM     11  CA  ALA A  39      20.950  -0.103  14.648  1.00 20.00           C
ATOM     12  C   ALA A  39      19.925   0.874  14.089  1.00 20.00           C
ATOM     13  O   ALA A  39      19.845   2.028  14.517  1.00 20.00           O
ATOM     14  CB  ALA A  39      20.363  -0.841  15.859  1.00 20.00           C
ATOM     15  N   ALA A  40      19.142   0.396  13.127  1.00 20.00           N
ATOM     16  CA  ALA A  40      18.108   1.191  12.480  1.00 20.00           C
ATOM     17  C   ALA A  40      16.741   0.588  12.768  1.00 20.00           C
ATOM     18  O   ALA A  40      16.573  -0.635  12.745  1.00 20.00           O
ATOM     19  CB  ALA A  40      18.337   1.274  10.967  1.00 20.00           C
ATOM     20  N   ALA A  41      15.767   1.453  13.040  1.00 20.00           N
ATOM     21  CA  ALA A  41      14.401   1.043  13.338  1.00 20.00           C
ATOM     22  C   ALA A  41      13.470   1.644  12.296  1.00 20.00           C
ATOM     23  O   ALA A  41      13.452   2.865  12.106  1.00 20.00           O
ATOM     24  CB  ALA A  41      13.989   1.477  14.746  1.00 20.00           C
ATOM     25  N   ALA A  42      12.701   0.790  11.627  1.00 20.00           N
ATOM     26  CA  ALA A  42      11.760   1.215  10.599  1.00 20.00           C
ATOM     27  C   ALA A  42      10.357   1.223  11.194  1.00 20.00           C
ATOM     28  O   ALA A  42       9.857   0.180  11.630  1.00 20.00           O
ATOM     29  N   ALA A  43       9.727   2.394  11.209  1.00 20.00           N
ATOM     30  CA  ALA A  43       8.385   2.564  11.743  1.00 20.00           C
ATOM     31  C   ALA A  43       7.444   3.011  10.634  1.00 20.00           C
ATOM     32  O   ALA A  43       7.806   3.844   9.796  1.00 20.00           O
ATOM     33  CB  ALA A  43       8.367   3.583  12.887  1.00 20.00           C
ATOM     34  N   ALA A  44       6.235   2.453  10.633  1.00 20.00           N
ATOM     35  CA  ALA A  44       5.236   2.788   9.630  1.00 20.00           C
ATOM     36  C   ALA A  44       3.852   2.693  10.255  1.00 20.00           C
ATOM     37  O   ALA A  44       3.636   1.941  11.209  1.00 20.00           O
ATOM     38  CB  ALA A  44       5.332   1.870   8.406  1.00 20.00           C
ATOM     39  N   ALA A  45       2.919   3.463   9.705  1.00 20.00           N
ATOM     40  CA  ALA A  45       1.548   3.476  10.200  1.00 20.00           C
ATOM     41  C   ALA A  45       0.715   2.390   9.527  1.00 20.00           C
ATOM     42  O   ALA A  45       1.135   1.800   8.531  1.00 20.00           O
ATOM     43  CB  ALA A  45       0.917   4.842   9.979  1.00 20.00           C
ATOM     44  N   ALA A  74      -2.563   6.328  17.723  1.00 20.00           N
ATOM     45  CA  ALA A  74      -2.435   5.677  16.425  1.00 20.00           C
ATOM     46  C   ALA A  74      -1.470   4.499  16.498  1.00 20.00           C
ATOM     47  O   ALA A  74      -0.401   4.596  17.100  1.00 20.00           O
ATOM     48  CB  ALA A  74      -1.975   6.676  15.374  1.00 20.00           C
ATOM     49  N   ALA A  75      -1.856   3.385  15.882  1.00 20.00           N
ATOM     50  CA  ALA A  75      -1.023   2.191  15.881  1.00 20.00           C
ATOM     51  C   ALA A  75       0.069   2.303  14.825  1.00 20.00           C
ATOM     52  O   ALA A  75      -0.171   2.768  13.707  1.00 20.00           O
ATOM     53  CB  ALA A  75      -1.873   0.946  15.627  1.00 20.00           C
ATOM     54  N   ALA A  76       1.275   1.871  15.189  1.00 20.00           N
ATOM     55  CA  ALA A  76       2.411   1.917  14.283  1.00 20.00           C
ATOM     56  C   ALA A  76       3.331   0.742  14.576  1.00 20.00           C
ATOM     57  O   ALA A  76       3.469   0.314  15.726  1.00 20.00           O
ATOM     58  CB  ALA A  76       3.180   3.238  14.405  1.00 20.00           C
ATOM     59  N   ALA A  77       3.957   0.225  13.525  1.00 20.00           N
ATOM     60  CA  ALA A  77       4.871  -0.901  13.624  1.00 20.00           C
ATOM     61  C   ALA A  77       6.309  -0.431  13.446  1.00 20.00           C
ATOM     62  O   ALA A  77       6.578   0.615  12.849  1.00 20.00           O
ATOM     63  CB  ALA A  77       4.537  -1.975  12.583  1.00 20.00           C
ATOM     64  N   ALA A  78       7.237  -1.225  13.977  1.00 20.00           N
ATOM     65  CA  ALA A  78       8.654  -0.911  13.889  1.00 20.00           C
ATOM     66  C   ALA A  78       9.447  -2.207  13.806  1.00 20.00           C
ATOM     67  O   ALA A  78       8.994  -3.266  14.248  1.00 20.00           O
ATOM     68  CB  ALA A  78       9.122  -0.070  15.082  1.00 20.00           C
ATOM     69  N   ALA A  79      10.643  -2.108  13.230  1.00 20.00           N
ATOM     70  CA  ALA A  79      11.528  -3.254  13.071  1.00 20.00           C
ATOM     71  C   ALA A  79      12.965  -2.807  13.282  1.00 20.00           C
ATOM     72  O   ALA A  79      13.402  -1.815  12.692  1.00 20.00           O
ATOM     73  CB  ALA A  79      11.368  -3.894  11.687  1.00 20.00           C
ATOM     74  N   ALA A  80      13.693  -3.538  14.121  1.00 20.00           N
ATOM     75  CA  ALA A  80      15.083  -3.235  14.423  1.00 20.00           C
ATOM     76  C   ALA A  80      16.005  -4.091  13.566  1.00 20.00           C
ATOM     77  O   ALA A  80      15.755  -5.284  13.366  1.00 20.00           O
ATOM     78  CB  ALA A  80      15.383  -3.466  15.905  1.00 20.00           C
ATOM     79  N   ALA A  81      17.071  -3.475  13.060  1.00 20.00           N
ATOM     80  CA  ALA A  81      18.037  -4.171  12.222  1.00 20.00           C
ATOM     81  C   ALA A  81      19.414  -3.568  12.449  1.00 20.00           C
ATOM     82  O   ALA A  81      19.590  -2.353  12.320  1.00 20.00           O
ATOM     83  CB  ALA A  81      17.653  -4.084  10.741  1.00 20.00           C
ATOM     84  N   ALA A  82      20.381  -4.417  12.785  1.00 20.00           N
ATOM     85  CA  ALA A  82      21.741  -3.960  13.029  1.00 20.00           C
ATOM     86  C   ALA A  82      22.475  -3.742  11.713  1.00 20.00           C
ATOM     87  O   ALA A  82      22.392  -4.561  10.793  1.00 20.00           O
ATOM     88  CB  ALA A  82      22.499  -4.970  13.891  1.00 20.00           C
ATOM     89  N   ALA A  83      23.197  -2.628  11.628  1.00 20.00           N
ATOM     90  CA  ALA A  83      23.957  -2.270  10.443  1.00 20.00           C
ATOM     91  C   ALA A  83      25.415  -2.029  10.812  1.00 20.00           C
ATOM     92  O   ALA A  83      25.748  -1.734  11.963  1.00 20.00           O
ATOM     93  CB  ALA A  83      23.375  -1.024   9.761  1.00 20.00           C
ATOM     94  N   ALA A  84      26.284  -2.160   9.815  1.00 20.00           N
ATOM     95  CA  ALA A  84      27.714  -1.960  10.014  1.00 20.00           C
ATOM     96  C   ALA A  84      28.356  -1.340   8.778  1.00 20.00           C
ATOM     97  O   ALA A  84      27.800  -1.400   7.681  1.00 20.00           O
ATOM     98  CB  ALA A  84      28.391  -3.277  10.360  1.00 20.00           C
ATOM     99  N   ALA A  94       8.798  -7.831  13.977  1.00 20.00           N
ATOM    100  CA  ALA A  94       8.087  -6.560  13.921  1.00 20.00           C
ATOM    101  C   ALA A  94       7.265  -6.339  15.186  1.00 20.00           C
ATOM    102  O   ALA A  94       6.554  -7.235  15.641  1.00 20.00           O
ATOM    103  CB  ALA A  94       7.193  -6.507  12.691  1.00 20.00           C
ATOM    104  N   ALA A  95       7.368  -5.139  15.750  1.00 20.00           N
ATOM    105  CA  ALA A  95       6.649  -4.772  16.959  1.00 20.00           C
ATOM    106  C   ALA A  95       5.563  -3.754  16.637  1.00 20.00           C
ATOM    107  O   ALA A  95       5.679  -2.974  15.687  1.00 20.00           O
ATOM    108  CB  ALA A  95       7.601  -4.200  18.015  1.00 20.00           C
ATOM    109  N   ALA A  96       4.502  -3.771  17.439  1.00 20.00           N
ATOM    110  CA  ALA A  96       3.371  -2.869  17.277  1.00 20.00           C
ATOM    111  C   ALA A  96       3.194  -2.052  18.548  1.00 20.00           C
ATOM    112  O   ALA A  96       3.221  -2.603  19.654  1.00 20.00           O
ATOM    113  CB  ALA A  96       2.088  -3.642  16.958  1.00 20.00           C
ATOM    114  N   ALA A  97       3.014  -0.743  18.388  1.00 20.00           N
ATOM    115  CA  ALA A  97       2.831   0.162  19.512  1.00 20.00           C
ATOM    116  C   ALA A  97       1.807   1.223  19.144  1.00 20.00           C
ATOM    117  O   ALA A  97       1.829   1.760  18.033  1.00 20.00           O
ATOM    118  CB  ALA A  97       4.153   0.825  19.921  1.00 20.00           C
ATOM    119  N   ALA A  98       0.913   1.521  20.082  1.00 20.00           N
ATOM    120  CA  ALA A  98      -0.125   2.521  19.861  1.00 20.00           C
ATOM    121  C   ALA A  98       0.164   3.791  20.655  1.00 20.00           C
ATOM    122  O   ALA A  98       0.052   3.806  21.881  1.00 20.00           O
ATOM    123  CB  ALA A  98      -1.489   1.961  20.232  1.00 20.00           C
TER
END
""")
  pdb_hierarchy = pdb_in.construct_hierarchy()
  pdb_hierarchy.atoms().reset_i_seq()
  ss_from_file = pdb_in.extract_secondary_structure()
  if not libtbx.env.has_module(name="ksdssp"):
    print("Skipping KSDSSP tests: ksdssp module not available.")
    return
  m = manager(pdb_hierarchy=pdb_hierarchy,
    sec_str_from_pdb_file=None)
  # m.params.secondary_structure.use_ksdssp = False
  m.params.secondary_structure.protein.remove_outliers = False
  log = open("exercise_sheet_ends.log", "w")
  proxies, hb_angles = m.create_protein_hbond_proxies(annotation=ss_from_file, log=log)
  log.close()
  # probably the determination of SS by ksdssp is not stable.
  assert proxies.size() == 12, proxies.size()
  assert hb_angles.size() == 27, hb_angles.size()


def exercise_helix_bonding_pattern():
  alpha_h1_simple = iotbx.pdb.input(source_info=None, lines="""\
ATOM      1  N   ALA     1       1.643  -2.366  -1.408  1.00  0.00           N
ATOM      2  CA  ALA     1       1.280  -3.608  -2.069  1.00  0.00           C
ATOM      3  C   ALA     1      -0.114  -3.466  -2.684  1.00  0.00           C
ATOM      4  O   ALA     1      -0.327  -3.827  -3.840  1.00  0.00           O
ATOM      5  CB  ALA     1       1.361  -4.762  -1.068  1.00  0.00           C
ATOM      6  N   ALA     2      -1.028  -2.938  -1.882  1.00  0.00           N
ATOM      7  CA  ALA     2      -2.395  -2.743  -2.332  1.00  0.00           C
ATOM      8  C   ALA     2      -2.396  -1.855  -3.579  1.00  0.00           C
ATOM      9  O   ALA     2      -3.059  -2.167  -4.567  1.00  0.00           O
ATOM     10  CB  ALA     2      -3.228  -2.150  -1.194  1.00  0.00           C
ATOM     11  N   ALA     3      -1.646  -0.767  -3.491  1.00  0.00           N
ATOM     12  CA  ALA     3      -1.551   0.168  -4.599  1.00  0.00           C
ATOM     13  C   ALA     3      -1.044  -0.568  -5.841  1.00  0.00           C
ATOM     14  O   ALA     3      -1.601  -0.419  -6.927  1.00  0.00           O
ATOM     15  CB  ALA     3      -0.646   1.337  -4.205  1.00  0.00           C
ATOM     16  N   ALA     4       0.008  -1.348  -5.639  1.00  0.00           N
ATOM     17  CA  ALA     4       0.597  -2.109  -6.728  1.00  0.00           C
ATOM     18  C   ALA     4      -0.466  -3.023  -7.340  1.00  0.00           C
ATOM     19  O   ALA     4      -0.611  -3.085  -8.559  1.00  0.00           O
ATOM     20  CB  ALA     4       1.808  -2.887  -6.211  1.00  0.00           C
ATOM     21  N   ALA     5      -1.184  -3.711  -6.463  1.00  0.00           N
ATOM     22  CA  ALA     5      -2.231  -4.619  -6.901  1.00  0.00           C
ATOM     23  C   ALA     5      -3.253  -3.847  -7.737  1.00  0.00           C
ATOM     24  O   ALA     5      -3.647  -4.296  -8.813  1.00  0.00           O
ATOM     25  CB  ALA     5      -2.864  -5.294  -5.683  1.00  0.00           C
ATOM     26  N   ALA     6      -3.654  -2.699  -7.211  1.00  0.00           N
ATOM     27  CA  ALA     6      -4.623  -1.860  -7.896  1.00  0.00           C
ATOM     28  C   ALA     6      -4.090  -1.499  -9.284  1.00  0.00           C
ATOM     29  O   ALA     6      -4.809  -1.602 -10.276  1.00  0.00           O
ATOM     30  CB  ALA     6      -4.919  -0.623  -7.045  1.00  0.00           C
ATOM     31  N   ALA     7      -2.831  -1.084  -9.309  1.00  0.00           N
ATOM     32  CA  ALA     7      -2.192  -0.708 -10.559  1.00  0.00           C
ATOM     33  C   ALA     7      -2.243  -1.890 -11.529  1.00  0.00           C
ATOM     34  O   ALA     7      -2.600  -1.727 -12.695  1.00  0.00           O
ATOM     35  CB  ALA     7      -0.761  -0.243 -10.281  1.00  0.00           C
ATOM     36  N   ALA     8      -1.881  -3.055 -11.012  1.00  0.00           N
ATOM     37  CA  ALA     8      -1.882  -4.264 -11.817  1.00  0.00           C
ATOM     38  C   ALA     8      -3.285  -4.496 -12.382  1.00  0.00           C
ATOM     39  O   ALA     8      -3.442  -4.772 -13.570  1.00  0.00           O
ATOM     40  CB  ALA     8      -1.391  -5.441 -10.972  1.00  0.00           C
ATOM     41  N   ALA     9      -4.269  -4.376 -11.503  1.00  0.00           N
ATOM     42  CA  ALA     9      -5.653  -4.568 -11.898  1.00  0.00           C
ATOM     43  C   ALA     9      -6.000  -3.590 -13.022  1.00  0.00           C
ATOM     44  O   ALA     9      -6.590  -3.978 -14.029  1.00  0.00           O
ATOM     45  CB  ALA     9      -6.561  -4.400 -10.678  1.00  0.00           C
ATOM     46  N   ALA    10      -5.617  -2.338 -12.812  1.00  0.00           N
ATOM     47  CA  ALA    10      -5.879  -1.301 -13.795  1.00  0.00           C
ATOM     48  C   ALA    10      -5.242  -1.695 -15.130  1.00  0.00           C
ATOM     49  O   ALA    10      -5.880  -1.605 -16.177  1.00  0.00           O
ATOM     50  CB  ALA    10      -5.358   0.040 -13.274  1.00  0.00           C
TER
""").construct_hierarchy()
  alpha_h1_1ac = iotbx.pdb.input(source_info=None, lines="""\
ATOM      1  N  AALA     1       1.643  -2.366  -1.408  0.50  0.00           N
ATOM      2  CA AALA     1       1.280  -3.608  -2.069  0.50  0.00           C
ATOM      3  C  AALA     1      -0.114  -3.466  -2.684  0.50  0.00           C
ATOM      4  O  AALA     1      -0.327  -3.827  -3.840  0.50  0.00           O
ATOM      5  CB AALA     1       1.361  -4.762  -1.068  0.50  0.00           C
ATOM      6  N  BALA     1       1.743  -2.366  -1.408  0.50  0.00           N
ATOM      7  CA BALA     1       1.380  -3.608  -2.069  0.50  0.00           C
ATOM      8  C  BALA     1      -0.014  -3.466  -2.684  0.50  0.00           C
ATOM      9  O  BALA     1      -0.227  -3.827  -3.840  0.50  0.00           O
ATOM     10  CB BALA     1       1.461  -4.762  -1.068  0.50  0.00           C
ATOM     11  N   ALA     2      -1.028  -2.938  -1.882  1.00  0.00           N
ATOM     12  CA  ALA     2      -2.395  -2.743  -2.332  1.00  0.00           C
ATOM     13  C   ALA     2      -2.396  -1.855  -3.579  1.00  0.00           C
ATOM     14  O   ALA     2      -3.059  -2.167  -4.567  1.00  0.00           O
ATOM     15  CB  ALA     2      -3.228  -2.150  -1.194  1.00  0.00           C
ATOM     16  N   ALA     3      -1.646  -0.767  -3.491  1.00  0.00           N
ATOM     17  CA  ALA     3      -1.551   0.168  -4.599  1.00  0.00           C
ATOM     18  C   ALA     3      -1.044  -0.568  -5.841  1.00  0.00           C
ATOM     19  O   ALA     3      -1.601  -0.419  -6.927  1.00  0.00           O
ATOM     20  CB  ALA     3      -0.646   1.337  -4.205  1.00  0.00           C
ATOM     21  N   ALA     4       0.008  -1.348  -5.639  1.00  0.00           N
ATOM     22  CA  ALA     4       0.597  -2.109  -6.728  1.00  0.00           C
ATOM     23  C   ALA     4      -0.466  -3.023  -7.340  1.00  0.00           C
ATOM     24  O   ALA     4      -0.611  -3.085  -8.559  1.00  0.00           O
ATOM     25  CB  ALA     4       1.808  -2.887  -6.211  1.00  0.00           C
ATOM     26  N   ALA     5      -1.184  -3.711  -6.463  1.00  0.00           N
ATOM     27  CA  ALA     5      -2.231  -4.619  -6.901  1.00  0.00           C
ATOM     28  C   ALA     5      -3.253  -3.847  -7.737  1.00  0.00           C
ATOM     29  O   ALA     5      -3.647  -4.296  -8.813  1.00  0.00           O
ATOM     30  CB  ALA     5      -2.864  -5.294  -5.683  1.00  0.00           C
ATOM     31  N   ALA     6      -3.654  -2.699  -7.211  1.00  0.00           N
ATOM     32  CA  ALA     6      -4.623  -1.860  -7.896  1.00  0.00           C
ATOM     33  C   ALA     6      -4.090  -1.499  -9.284  1.00  0.00           C
ATOM     34  O   ALA     6      -4.809  -1.602 -10.276  1.00  0.00           O
ATOM     35  CB  ALA     6      -4.919  -0.623  -7.045  1.00  0.00           C
ATOM     36  N   ALA     7      -2.831  -1.084  -9.309  1.00  0.00           N
ATOM     37  CA  ALA     7      -2.192  -0.708 -10.559  1.00  0.00           C
ATOM     38  C   ALA     7      -2.243  -1.890 -11.529  1.00  0.00           C
ATOM     39  O   ALA     7      -2.600  -1.727 -12.695  1.00  0.00           O
ATOM     40  CB  ALA     7      -0.761  -0.243 -10.281  1.00  0.00           C
ATOM     41  N   ALA     8      -1.881  -3.055 -11.012  1.00  0.00           N
ATOM     42  CA  ALA     8      -1.882  -4.264 -11.817  1.00  0.00           C
ATOM     43  C   ALA     8      -3.285  -4.496 -12.382  1.00  0.00           C
ATOM     44  O   ALA     8      -3.442  -4.772 -13.570  1.00  0.00           O
ATOM     45  CB  ALA     8      -1.391  -5.441 -10.972  1.00  0.00           C
ATOM     46  N   ALA     9      -4.269  -4.376 -11.503  1.00  0.00           N
ATOM     47  CA  ALA     9      -5.653  -4.568 -11.898  1.00  0.00           C
ATOM     48  C   ALA     9      -6.000  -3.590 -13.022  1.00  0.00           C
ATOM     49  O   ALA     9      -6.590  -3.978 -14.029  1.00  0.00           O
ATOM     50  CB  ALA     9      -6.561  -4.400 -10.678  1.00  0.00           C
ATOM     51  N   ALA    10      -5.617  -2.338 -12.812  1.00  0.00           N
ATOM     52  CA  ALA    10      -5.879  -1.301 -13.795  1.00  0.00           C
ATOM     53  C   ALA    10      -5.242  -1.695 -15.130  1.00  0.00           C
ATOM     54  O   ALA    10      -5.880  -1.605 -16.177  1.00  0.00           O
ATOM     55  CB  ALA    10      -5.358   0.040 -13.274  1.00  0.00           C
TER      56      ALA    10
END
""").construct_hierarchy()

  alpha_h1_1_10ac = iotbx.pdb.input(source_info=None, lines="""\
ATOM      1  N  AALA     1       1.643  -2.366  -1.408  0.50  0.00           N
ATOM      2  CA AALA     1       1.280  -3.608  -2.069  0.50  0.00           C
ATOM      3  C  AALA     1      -0.114  -3.466  -2.684  0.50  0.00           C
ATOM      4  O  AALA     1      -0.327  -3.827  -3.840  0.50  0.00           O
ATOM      5  CB AALA     1       1.361  -4.762  -1.068  0.50  0.00           C
ATOM      6  N  BALA     1       1.743  -2.366  -1.408  0.50  0.00           N
ATOM      7  CA BALA     1       1.380  -3.608  -2.069  0.50  0.00           C
ATOM      8  C  BALA     1      -0.014  -3.466  -2.684  0.50  0.00           C
ATOM      9  O  BALA     1      -0.227  -3.827  -3.840  0.50  0.00           O
ATOM     10  CB BALA     1       1.461  -4.762  -1.068  0.50  0.00           C
ATOM     11  N   ALA     2      -1.028  -2.938  -1.882  1.00  0.00           N
ATOM     12  CA  ALA     2      -2.395  -2.743  -2.332  1.00  0.00           C
ATOM     13  C   ALA     2      -2.396  -1.855  -3.579  1.00  0.00           C
ATOM     14  O   ALA     2      -3.059  -2.167  -4.567  1.00  0.00           O
ATOM     15  CB  ALA     2      -3.228  -2.150  -1.194  1.00  0.00           C
ATOM     16  N   ALA     3      -1.646  -0.767  -3.491  1.00  0.00           N
ATOM     17  CA  ALA     3      -1.551   0.168  -4.599  1.00  0.00           C
ATOM     18  C   ALA     3      -1.044  -0.568  -5.841  1.00  0.00           C
ATOM     19  O   ALA     3      -1.601  -0.419  -6.927  1.00  0.00           O
ATOM     20  CB  ALA     3      -0.646   1.337  -4.205  1.00  0.00           C
ATOM     21  N   ALA     4       0.008  -1.348  -5.639  1.00  0.00           N
ATOM     22  CA  ALA     4       0.597  -2.109  -6.728  1.00  0.00           C
ATOM     23  C   ALA     4      -0.466  -3.023  -7.340  1.00  0.00           C
ATOM     24  O   ALA     4      -0.611  -3.085  -8.559  1.00  0.00           O
ATOM     25  CB  ALA     4       1.808  -2.887  -6.211  1.00  0.00           C
ATOM     26  N   ALA     5      -1.184  -3.711  -6.463  1.00  0.00           N
ATOM     27  CA  ALA     5      -2.231  -4.619  -6.901  1.00  0.00           C
ATOM     28  C   ALA     5      -3.253  -3.847  -7.737  1.00  0.00           C
ATOM     29  O   ALA     5      -3.647  -4.296  -8.813  1.00  0.00           O
ATOM     30  CB  ALA     5      -2.864  -5.294  -5.683  1.00  0.00           C
ATOM     31  N   ALA     6      -3.654  -2.699  -7.211  1.00  0.00           N
ATOM     32  CA  ALA     6      -4.623  -1.860  -7.896  1.00  0.00           C
ATOM     33  C   ALA     6      -4.090  -1.499  -9.284  1.00  0.00           C
ATOM     34  O   ALA     6      -4.809  -1.602 -10.276  1.00  0.00           O
ATOM     35  CB  ALA     6      -4.919  -0.623  -7.045  1.00  0.00           C
ATOM     36  N   ALA     7      -2.831  -1.084  -9.309  1.00  0.00           N
ATOM     37  CA  ALA     7      -2.192  -0.708 -10.559  1.00  0.00           C
ATOM     38  C   ALA     7      -2.243  -1.890 -11.529  1.00  0.00           C
ATOM     39  O   ALA     7      -2.600  -1.727 -12.695  1.00  0.00           O
ATOM     40  CB  ALA     7      -0.761  -0.243 -10.281  1.00  0.00           C
ATOM     41  N   ALA     8      -1.881  -3.055 -11.012  1.00  0.00           N
ATOM     42  CA  ALA     8      -1.882  -4.264 -11.817  1.00  0.00           C
ATOM     43  C   ALA     8      -3.285  -4.496 -12.382  1.00  0.00           C
ATOM     44  O   ALA     8      -3.442  -4.772 -13.570  1.00  0.00           O
ATOM     45  CB  ALA     8      -1.391  -5.441 -10.972  1.00  0.00           C
ATOM     46  N   ALA     9      -4.269  -4.376 -11.503  1.00  0.00           N
ATOM     47  CA  ALA     9      -5.653  -4.568 -11.898  1.00  0.00           C
ATOM     48  C   ALA     9      -6.000  -3.590 -13.022  1.00  0.00           C
ATOM     49  O   ALA     9      -6.590  -3.978 -14.029  1.00  0.00           O
ATOM     50  CB  ALA     9      -6.561  -4.400 -10.678  1.00  0.00           C
ATOM     51  N  AALA    10      -5.617  -2.338 -12.812  0.50  0.00           N
ATOM     52  CA AALA    10      -5.879  -1.301 -13.795  0.50  0.00           C
ATOM     53  C  AALA    10      -5.242  -1.695 -15.130  0.50  0.00           C
ATOM     54  O  AALA    10      -5.880  -1.605 -16.177  0.50  0.00           O
ATOM     55  CB AALA    10      -5.358   0.040 -13.274  0.50  0.00           C
ATOM     57  N  BALA    10      -5.517  -2.338 -12.812  0.50  0.00           N
ATOM     58  CA BALA    10      -5.779  -1.301 -13.795  0.50  0.00           C
ATOM     59  C  BALA    10      -5.142  -1.695 -15.130  0.50  0.00           C
ATOM     60  O  BALA    10      -5.780  -1.605 -16.177  0.50  0.00           O
ATOM     61  CB BALA    10      -5.258   0.040 -13.274  0.50  0.00           C
END
""").construct_hierarchy()

  alpha_h1_1_5_10ac = iotbx.pdb.input(source_info=None, lines="""\
ATOM      1  N  AALA     1       1.643  -2.366  -1.408  0.50  0.00           N
ATOM      2  CA AALA     1       1.280  -3.608  -2.069  0.50  0.00           C
ATOM      3  C  AALA     1      -0.114  -3.466  -2.684  0.50  0.00           C
ATOM      4  O  AALA     1      -0.327  -3.827  -3.840  0.50  0.00           O
ATOM      5  CB AALA     1       1.361  -4.762  -1.068  0.50  0.00           C
ATOM      6  N  BALA     1       1.743  -2.366  -1.408  0.50  0.00           N
ATOM      7  CA BALA     1       1.380  -3.608  -2.069  0.50  0.00           C
ATOM      8  C  BALA     1      -0.014  -3.466  -2.684  0.50  0.00           C
ATOM      9  O  BALA     1      -0.227  -3.827  -3.840  0.50  0.00           O
ATOM     10  CB BALA     1       1.461  -4.762  -1.068  0.50  0.00           C
ATOM     11  N   ALA     2      -1.028  -2.938  -1.882  1.00  0.00           N
ATOM     12  CA  ALA     2      -2.395  -2.743  -2.332  1.00  0.00           C
ATOM     13  C   ALA     2      -2.396  -1.855  -3.579  1.00  0.00           C
ATOM     14  O   ALA     2      -3.059  -2.167  -4.567  1.00  0.00           O
ATOM     15  CB  ALA     2      -3.228  -2.150  -1.194  1.00  0.00           C
ATOM     16  N   ALA     3      -1.646  -0.767  -3.491  1.00  0.00           N
ATOM     17  CA  ALA     3      -1.551   0.168  -4.599  1.00  0.00           C
ATOM     18  C   ALA     3      -1.044  -0.568  -5.841  1.00  0.00           C
ATOM     19  O   ALA     3      -1.601  -0.419  -6.927  1.00  0.00           O
ATOM     20  CB  ALA     3      -0.646   1.337  -4.205  1.00  0.00           C
ATOM     21  N   ALA     4       0.008  -1.348  -5.639  1.00  0.00           N
ATOM     22  CA  ALA     4       0.597  -2.109  -6.728  1.00  0.00           C
ATOM     23  C   ALA     4      -0.466  -3.023  -7.340  1.00  0.00           C
ATOM     24  O   ALA     4      -0.611  -3.085  -8.559  1.00  0.00           O
ATOM     25  CB  ALA     4       1.808  -2.887  -6.211  1.00  0.00           C
ATOM     26  N  AALA     5      -1.184  -3.711  -6.463  0.50  0.00           N
ATOM     27  CA AALA     5      -2.231  -4.619  -6.901  0.50  0.00           C
ATOM     28  C  AALA     5      -3.253  -3.847  -7.737  0.50  0.00           C
ATOM     29  O  AALA     5      -3.647  -4.296  -8.813  0.50  0.00           O
ATOM     30  CB AALA     5      -2.864  -5.294  -5.683  0.50  0.00           C
ATOM     31  N  BALA     5      -1.084  -3.711  -6.463  0.50  0.00           N
ATOM     32  CA BALA     5      -2.131  -4.619  -6.901  0.50  0.00           C
ATOM     33  C  BALA     5      -3.153  -3.847  -7.737  0.50  0.00           C
ATOM     34  O  BALA     5      -3.547  -4.296  -8.813  0.50  0.00           O
ATOM     35  CB BALA     5      -2.764  -5.294  -5.683  0.50  0.00           C
ATOM     36  N   ALA     6      -3.654  -2.699  -7.211  1.00  0.00           N
ATOM     37  CA  ALA     6      -4.623  -1.860  -7.896  1.00  0.00           C
ATOM     38  C   ALA     6      -4.090  -1.499  -9.284  1.00  0.00           C
ATOM     39  O   ALA     6      -4.809  -1.602 -10.276  1.00  0.00           O
ATOM     40  CB  ALA     6      -4.919  -0.623  -7.045  1.00  0.00           C
ATOM     41  N   ALA     7      -2.831  -1.084  -9.309  1.00  0.00           N
ATOM     42  CA  ALA     7      -2.192  -0.708 -10.559  1.00  0.00           C
ATOM     43  C   ALA     7      -2.243  -1.890 -11.529  1.00  0.00           C
ATOM     44  O   ALA     7      -2.600  -1.727 -12.695  1.00  0.00           O
ATOM     45  CB  ALA     7      -0.761  -0.243 -10.281  1.00  0.00           C
ATOM     46  N   ALA     8      -1.881  -3.055 -11.012  1.00  0.00           N
ATOM     47  CA  ALA     8      -1.882  -4.264 -11.817  1.00  0.00           C
ATOM     48  C   ALA     8      -3.285  -4.496 -12.382  1.00  0.00           C
ATOM     49  O   ALA     8      -3.442  -4.772 -13.570  1.00  0.00           O
ATOM     50  CB  ALA     8      -1.391  -5.441 -10.972  1.00  0.00           C
ATOM     51  N   ALA     9      -4.269  -4.376 -11.503  1.00  0.00           N
ATOM     52  CA  ALA     9      -5.653  -4.568 -11.898  1.00  0.00           C
ATOM     53  C   ALA     9      -6.000  -3.590 -13.022  1.00  0.00           C
ATOM     54  O   ALA     9      -6.590  -3.978 -14.029  1.00  0.00           O
ATOM     55  CB  ALA     9      -6.561  -4.400 -10.678  1.00  0.00           C
ATOM     56  N  AALA    10      -5.617  -2.338 -12.812  0.50  0.00           N
ATOM     57  CA AALA    10      -5.879  -1.301 -13.795  0.50  0.00           C
ATOM     58  C  AALA    10      -5.242  -1.695 -15.130  0.50  0.00           C
ATOM     59  O  AALA    10      -5.880  -1.605 -16.177  0.50  0.00           O
ATOM     60  CB AALA    10      -5.358   0.040 -13.274  0.50  0.00           C
ATOM     62  N  BALA    10      -5.517  -2.338 -12.812  0.50  0.00           N
ATOM     63  CA BALA    10      -5.779  -1.301 -13.795  0.50  0.00           C
ATOM     64  C  BALA    10      -5.142  -1.695 -15.130  0.50  0.00           C
ATOM     65  O  BALA    10      -5.780  -1.605 -16.177  0.50  0.00           O
ATOM     66  CB BALA    10      -5.258   0.040 -13.274  0.50  0.00           C
END
""").construct_hierarchy()

  alpha_h1_1_5_10ac_insertions = iotbx.pdb.input(source_info=None, lines="""\
ATOM      1  N  AALA     1       1.643  -2.366  -1.408  0.50  0.00           N
ATOM      2  CA AALA     1       1.280  -3.608  -2.069  0.50  0.00           C
ATOM      3  C  AALA     1      -0.114  -3.466  -2.684  0.50  0.00           C
ATOM      4  O  AALA     1      -0.327  -3.827  -3.840  0.50  0.00           O
ATOM      5  CB AALA     1       1.361  -4.762  -1.068  0.50  0.00           C
ATOM      6  N  BALA     1       1.743  -2.366  -1.408  0.50  0.00           N
ATOM      7  CA BALA     1       1.380  -3.608  -2.069  0.50  0.00           C
ATOM      8  C  BALA     1      -0.014  -3.466  -2.684  0.50  0.00           C
ATOM      9  O  BALA     1      -0.227  -3.827  -3.840  0.50  0.00           O
ATOM     10  CB BALA     1       1.461  -4.762  -1.068  0.50  0.00           C
ATOM     11  N   ALA     2      -1.028  -2.938  -1.882  1.00  0.00           N
ATOM     12  CA  ALA     2      -2.395  -2.743  -2.332  1.00  0.00           C
ATOM     13  C   ALA     2      -2.396  -1.855  -3.579  1.00  0.00           C
ATOM     14  O   ALA     2      -3.059  -2.167  -4.567  1.00  0.00           O
ATOM     15  CB  ALA     2      -3.228  -2.150  -1.194  1.00  0.00           C
ATOM     16  N   ALA     2A     -1.646  -0.767  -3.491  1.00  0.00           N
ATOM     17  CA  ALA     2A     -1.551   0.168  -4.599  1.00  0.00           C
ATOM     18  C   ALA     2A     -1.044  -0.568  -5.841  1.00  0.00           C
ATOM     19  O   ALA     2A     -1.601  -0.419  -6.927  1.00  0.00           O
ATOM     20  CB  ALA     2A     -0.646   1.337  -4.205  1.00  0.00           C
ATOM     21  N   ALA     2B      0.008  -1.348  -5.639  1.00  0.00           N
ATOM     22  CA  ALA     2B      0.597  -2.109  -6.728  1.00  0.00           C
ATOM     23  C   ALA     2B     -0.466  -3.023  -7.340  1.00  0.00           C
ATOM     24  O   ALA     2B     -0.611  -3.085  -8.559  1.00  0.00           O
ATOM     25  CB  ALA     2B      1.808  -2.887  -6.211  1.00  0.00           C
ATOM     26  N  AALA     2C     -1.184  -3.711  -6.463  0.50  0.00           N
ATOM     27  CA AALA     2C     -2.231  -4.619  -6.901  0.50  0.00           C
ATOM     28  C  AALA     2C     -3.253  -3.847  -7.737  0.50  0.00           C
ATOM     29  O  AALA     2C     -3.647  -4.296  -8.813  0.50  0.00           O
ATOM     30  CB AALA     2C     -2.864  -5.294  -5.683  0.50  0.00           C
ATOM     31  N  BALA     2C     -1.084  -3.711  -6.463  0.50  0.00           N
ATOM     32  CA BALA     2C     -2.131  -4.619  -6.901  0.50  0.00           C
ATOM     33  C  BALA     2C     -3.153  -3.847  -7.737  0.50  0.00           C
ATOM     34  O  BALA     2C     -3.547  -4.296  -8.813  0.50  0.00           O
ATOM     35  CB BALA     2C     -2.764  -5.294  -5.683  0.50  0.00           C
ATOM     36  N   ALA     2D     -3.654  -2.699  -7.211  1.00  0.00           N
ATOM     37  CA  ALA     2D     -4.623  -1.860  -7.896  1.00  0.00           C
ATOM     38  C   ALA     2D     -4.090  -1.499  -9.284  1.00  0.00           C
ATOM     39  O   ALA     2D     -4.809  -1.602 -10.276  1.00  0.00           O
ATOM     40  CB  ALA     2D     -4.919  -0.623  -7.045  1.00  0.00           C
ATOM     41  N   ALA     3      -2.831  -1.084  -9.309  1.00  0.00           N
ATOM     42  CA  ALA     3      -2.192  -0.708 -10.559  1.00  0.00           C
ATOM     43  C   ALA     3      -2.243  -1.890 -11.529  1.00  0.00           C
ATOM     44  O   ALA     3      -2.600  -1.727 -12.695  1.00  0.00           O
ATOM     45  CB  ALA     3      -0.761  -0.243 -10.281  1.00  0.00           C
ATOM     46  N   ALA     4      -1.881  -3.055 -11.012  1.00  0.00           N
ATOM     47  CA  ALA     4      -1.882  -4.264 -11.817  1.00  0.00           C
ATOM     48  C   ALA     4      -3.285  -4.496 -12.382  1.00  0.00           C
ATOM     49  O   ALA     4      -3.442  -4.772 -13.570  1.00  0.00           O
ATOM     50  CB  ALA     4      -1.391  -5.441 -10.972  1.00  0.00           C
ATOM     51  N   ALA     5      -4.269  -4.376 -11.503  1.00  0.00           N
ATOM     52  CA  ALA     5      -5.653  -4.568 -11.898  1.00  0.00           C
ATOM     53  C   ALA     5      -6.000  -3.590 -13.022  1.00  0.00           C
ATOM     54  O   ALA     5      -6.590  -3.978 -14.029  1.00  0.00           O
ATOM     55  CB  ALA     5      -6.561  -4.400 -10.678  1.00  0.00           C
ATOM     56  N  AALA     6      -5.617  -2.338 -12.812  0.50  0.00           N
ATOM     57  CA AALA     6      -5.879  -1.301 -13.795  0.50  0.00           C
ATOM     58  C  AALA     6      -5.242  -1.695 -15.130  0.50  0.00           C
ATOM     59  O  AALA     6      -5.880  -1.605 -16.177  0.50  0.00           O
ATOM     60  CB AALA     6      -5.358   0.040 -13.274  0.50  0.00           C
ATOM     62  N  BALA     6      -5.517  -2.338 -12.812  0.50  0.00           N
ATOM     63  CA BALA     6      -5.779  -1.301 -13.795  0.50  0.00           C
ATOM     64  C  BALA     6      -5.142  -1.695 -15.130  0.50  0.00           C
ATOM     65  O  BALA     6      -5.780  -1.605 -16.177  0.50  0.00           O
ATOM     66  CB BALA     6      -5.258   0.040 -13.274  0.50  0.00           C
END
""").construct_hierarchy()

  alpha_annot_1 = """\
HELIX    1   1 ALA      1  ALA     10  1                                  10
"""

  alpha_annot_2 = """\
HELIX    1   1 ALA      2  ALA      5  1                                  10
"""

  log = null_out()
  n_hbonds = []
  n_hangles = []
  for pdb_h, recs in [
                        (alpha_h1_simple, alpha_annot_1),
                        (alpha_h1_1ac, alpha_annot_1),
                        (alpha_h1_1_10ac, alpha_annot_1),
                        (alpha_h1_1_5_10ac, alpha_annot_1),
                        (alpha_h1_1_5_10ac_insertions, alpha_annot_1),
                        (alpha_h1_1_5_10ac_insertions, alpha_annot_2),
                        ]:
    ioss_annotation = ioss.annotation.from_records(records = recs.split('\n'))
    ann = ioss_annotation.as_restraint_groups(prefix_scope="secondary_structure")
    defpars = iotbx.phil.parse(sec_str_master_phil_str)
    custom_pars = defpars.fetch(iotbx.phil.parse(ann))
    custom_pars_ex = custom_pars.extract()
    ss_manager = manager(
                pdb_h,
                sec_str_from_pdb_file=None,
                params=custom_pars_ex.secondary_structure,
                verbose=-1)
    proxies_for_grm, angle_p = ss_manager.create_protein_hbond_proxies(
      annotation= None,
      log          = log)
    # pdb_h.write_pdb_file(file_name="alpha_h1_1_5_10ac.pdb")
    # print proxies_for_grm.as_pymol_dashes(pdb_inp.hierarchy)
    # print angle_p.as_pymol_dashes(pdb_hierarchy=pdb_inp.hierarchy)
    # print "hbonds, hangles:", proxies_for_grm.size(), angle_p.size()
    # STOP()
    n_hbonds.append(proxies_for_grm.size())
    n_hangles.append(angle_p.size())
  print(n_hbonds, n_hangles)
  assert n_hbonds ==  [ 6,  7,  8,  9,  9,  5]
  assert n_hangles == [18, 21, 24, 28, 28, 16]

def exercise_sheets_bonding_pattern():
  pdb_apar_h = iotbx.pdb.input(source_info=None, lines = """\
SCRYST1   46.460   46.460  193.210  90.00  90.00 120.00 P 31 2 1
SCALE1      0.021524  0.012427  0.000000        0.00000
SCALE2      0.000000  0.024854  0.000000        0.00000
SCALE3      0.000000  0.000000  0.005176        0.00000
ATOM      1  N   GLY A   1      28.066 -23.487   2.817  1.00 16.21           N
ATOM      2  CA  GLY A   1      27.219 -22.652   1.936  1.00 13.03           C
ATOM      3  C   GLY A   1      27.859 -21.336   1.546  1.00  9.56           C
ATOM      4  O   GLY A   1      28.868 -20.882   2.096  1.00 12.48           O
ATOM      5  N   TYR A   2      27.208 -20.701   0.590  1.00  7.29           N
ATOM      6  CA  TYR A   2      27.617 -19.424   0.052  1.00  7.96           C
ATOM      7  C   TYR A   2      26.483 -18.436   0.263  1.00  6.87           C
ATOM      8  O   TYR A   2      25.303 -18.771   0.249  1.00  6.97           O
ATOM      9  CB  TYR A   2      27.861 -19.541  -1.451  1.00  7.90           C
ATOM     10  CG  TYR A   2      28.902 -20.556  -1.857  1.00  9.09           C
ATOM     11  CD1 TYR A   2      30.255 -20.336  -1.592  1.00 11.43           C
ATOM     12  CD2 TYR A   2      28.566 -21.697  -2.545  1.00 10.59           C
ATOM     13  CE1 TYR A   2      31.227 -21.244  -1.987  1.00 13.93           C
ATOM     14  CE2 TYR A   2      29.518 -22.630  -2.915  1.00 11.76           C
ATOM     15  CZ  TYR A   2      30.847 -22.395  -2.659  1.00 13.67           C
ATOM     16  OH  TYR A   2      31.792 -23.309  -3.059  1.00 18.26           O
ATOM     17  N   SER A   3      26.854 -17.177   0.412  1.00  6.66           N
ATOM     18  CA  SER A   3      25.899 -16.083   0.383  0.53  6.83           C
ATOM     20  C   SER A   3      26.500 -14.946  -0.440  1.00  5.61           C
ATOM     21  O   SER A   3      27.729 -14.822  -0.565  1.00  8.14           O
ATOM     22  CB  SER A   3      25.569 -15.634   1.795  0.53  7.38           C
ATOM     24  OG  SER A   3      26.740 -15.136   2.390  0.53  9.79           O
ATOM     26  N   CYS A   4      25.627 -14.135  -0.995  1.00  5.11           N
ATOM     27  CA  CYS A   4      26.070 -13.062  -1.865  1.00  5.65           C
ATOM     28  C   CYS A   4      25.043 -11.934  -1.798  1.00  4.54           C
ATOM     29  O   CYS A   4      23.856 -12.180  -1.528  1.00  5.60           O
ATOM     30  CB  CYS A   4      26.253 -13.557  -3.295  1.00  7.00           C
ATOM     31  SG  CYS A   4      24.806 -14.269  -4.119  1.00  8.88           S
ATOM     32  N   ARG A   5      25.486 -10.691  -2.002  1.00  4.69           N
ATOM     33  CA  ARG A   5      24.558  -9.545  -2.064  1.00  4.87           C
ATOM     34  C   ARG A   5      25.196  -8.395  -2.796  1.00  4.57           C
ATOM     35  O   ARG A   5      26.416  -8.238  -2.793  1.00  5.50           O
ATOM     36  CB  ARG A   5      24.061  -9.108  -0.700  1.00  6.43           C
ATOM     37  CG  ARG A   5      25.121  -8.566   0.219  1.00  6.96           C
ATOM     38  CD  ARG A   5      24.461  -8.032   1.494  1.00  7.25           C
ATOM     39  NE  ARG A   5      25.452  -7.547   2.440  1.00  7.63           N
ATOM     40  CZ  ARG A   5      26.107  -8.341   3.280  1.00  9.10           C
ATOM     41  NH1 ARG A   5      25.867  -9.642   3.297  1.00  9.68           N
ATOM     42  NH2 ARG A   5      26.974  -7.836   4.146  1.00 10.30           N
ATOM     43  N   ALA A   6      24.325  -7.563  -3.358  1.00  4.39           N
ATOM     44  CA  ALA A   6      24.723  -6.362  -4.067  1.00  4.73           C
ATOM     45  C   ALA A   6      23.693  -5.275  -3.769  1.00  4.15           C
ATOM     46  O   ALA A   6      22.482  -5.458  -3.987  1.00  4.96           O
ATOM     47  CB  ALA A   6      24.831  -6.626  -5.558  1.00  5.96           C
ATOM     48  N   VAL A   7      24.165  -4.139  -3.284  1.00  4.97           N
ATOM     49  CA  VAL A   7      23.374  -2.917  -3.085  1.00  4.42           C
ATOM     50  C   VAL A   7      23.482  -2.046  -4.325  1.00  4.45           C
ATOM     51  O   VAL A   7      24.589  -1.663  -4.717  1.00  5.55           O
ATOM     52  CB  VAL A   7      23.830  -2.159  -1.806  1.00  5.09           C
ATOM     53  CG1 VAL A   7      23.111  -0.842  -1.686  1.00  5.65           C
ATOM     54  CG2 VAL A   7      23.612  -2.998  -0.570  1.00  6.88           C
ATOM    204  N   MET A  31      18.177  -3.966  -4.656  1.00  4.72           N
ATOM    205  CA  MET A  31      18.833  -4.887  -3.744  1.00  5.29           C
ATOM    206  C   MET A  31      18.765  -6.294  -4.322  1.00  4.60           C
ATOM    207  O   MET A  31      17.661  -6.738  -4.657  1.00  4.98           O
ATOM    208  CB  MET A  31      18.097  -4.868  -2.387  1.00  6.39           C
ATOM    209  CG  MET A  31      18.723  -5.755  -1.334  1.00  8.61           C
ATOM    210  SD  MET A  31      20.248  -5.074  -0.655  1.00 11.04           S
ATOM    211  CE  MET A  31      21.358  -6.427  -0.777  1.00  8.94           C
ATOM    212  N   ALA A  32      19.899  -6.986  -4.412  1.00  3.90           N
ATOM    213  CA  ALA A  32      19.934  -8.380  -4.864  1.00  3.55           C
ATOM    214  C   ALA A  32      20.720  -9.194  -3.858  1.00  3.90           C
ATOM    215  O   ALA A  32      21.762  -8.763  -3.393  1.00  5.01           O
ATOM    216  CB  ALA A  32      20.552  -8.500  -6.255  1.00  4.55           C
ATOM    217  N   SER A  33      20.230 -10.407  -3.567  1.00  4.02           N
ATOM    218  CA  SER A  33      20.980 -11.288  -2.645  1.00  3.79           C
ATOM    219  C   SER A  33      20.591 -12.727  -2.850  1.00  4.33           C
ATOM    220  O   SER A  33      19.532 -13.045  -3.412  1.00  4.84           O
ATOM    221  CB  SER A  33      20.830 -10.880  -1.167  1.00  4.55           C
ATOM    222  OG  SER A  33      19.498 -11.105  -0.710  1.00  5.25           O
ATOM    223  N   GLY A  34      21.415 -13.600  -2.283  1.00  4.96           N
ATOM    224  CA  GLY A  34      21.104 -14.997  -2.335  1.00  4.50           C
ATOM    225  C   GLY A  34      21.914 -15.837  -1.397  1.00  4.64           C
ATOM    226  O   GLY A  34      22.836 -15.343  -0.732  1.00  4.96           O
ATOM    227  N   THR A  35      21.521 -17.105  -1.329  1.00  4.35           N
ATOM    228  CA  THR A  35      22.277 -18.138  -0.654  1.00  4.02           C
ATOM    229  C   THR A  35      22.226 -19.392  -1.495  1.00  4.55           C
ATOM    230  O   THR A  35      21.221 -19.652  -2.170  1.00  4.35           O
ATOM    231  CB  THR A  35      21.715 -18.436   0.762  1.00  5.12           C
ATOM    232  OG1 THR A  35      20.356 -18.929   0.668  1.00  5.51           O
ATOM    233  CG2 THR A  35      21.733 -17.222   1.670  1.00  5.97           C
ATOM    234  N   SER A  36      23.294 -20.178  -1.426  1.00  4.63           N
ATOM    235  CA  SER A  36      23.402 -21.406  -2.221  1.00  4.58           C
ATOM    236  C   SER A  36      24.387 -22.368  -1.553  1.00  5.05           C
ATOM    237  O   SER A  36      24.929 -22.095  -0.497  1.00  6.22           O
ATOM    238  CB  SER A  36      23.881 -21.071  -3.639  1.00  5.69           C
ATOM    239  OG  SER A  36      25.213 -20.561  -3.633  1.00  7.12           O
""").construct_hierarchy()

  pdb_par_h = iotbx.pdb.input(source_info=None, lines = """\
CRYST1   46.460   46.460  193.210  90.00  90.00 120.00 P 31 2 1
SCALE1      0.021524  0.012427  0.000000        0.00000
SCALE2      0.000000  0.024854  0.000000        0.00000
SCALE3      0.000000  0.000000  0.005176        0.00000
ATOM     67  N   ALA A  15       5.011  -5.031  -8.967  1.00  5.73           N
ATOM     68  CA  ALA A  15       4.943  -6.455  -9.287  1.00  6.16           C
ATOM     69  C   ALA A  15       5.610  -7.252  -8.166  1.00  6.39           C
ATOM     70  O   ALA A  15       6.751  -7.007  -7.857  1.00  9.87           O
ATOM     71  CB  ALA A  15       5.684  -6.739 -10.604  1.00  7.46           C
ATOM     72  N   THR A  16       4.929  -8.263  -7.636  1.00  5.49           N
ATOM     73  C   THR A  16       5.316 -10.600  -7.084  1.00  5.07           C
ATOM     74  O   THR A  16       4.214 -11.002  -7.422  1.00  6.51           O
ATOM     75  CA  THR A  16       5.513  -9.172  -6.657  1.00  5.70           C
ATOM     76  CB  THR A  16       4.864  -9.001  -5.276  1.00  8.31           C
ATOM     77  N   GLY A  17       6.393 -11.375  -7.067  1.00  4.56           N
ATOM     78  CA  GLY A  17       6.325 -12.770  -7.439  1.00  4.26           C
ATOM     79  C   GLY A  17       7.219 -13.654  -6.621  1.00  4.41           C
ATOM     80  O   GLY A  17       8.263 -13.233  -6.114  1.00  5.01           O
ATOM     81  N   SER A  18       6.827 -14.921  -6.561  1.00  4.24           N
ATOM     82  CA  SER A  18       7.657 -15.945  -5.959  1.00  4.02           C
ATOM     83  C   SER A  18       7.539 -17.244  -6.724  1.00  3.64           C
ATOM     84  O   SER A  18       6.482 -17.526  -7.331  1.00  4.07           O
ATOM     85  CB  SER A  18       7.335 -16.157  -4.481  1.00  5.86           C
ATOM     86  N   ALA A  19       8.573 -18.049  -6.627  1.00  3.19           N
ATOM     87  CA  ALA A  19       8.578 -19.414  -7.195  1.00  3.31           C
ATOM     88  C   ALA A  19       9.370 -20.307  -6.238  1.00  3.20           C
ATOM     89  O   ALA A  19      10.484 -19.956  -5.817  1.00  4.21           O
ATOM     90  CB  ALA A  19       9.235 -19.415  -8.574  1.00  3.85           C
ATOM     91  N   THR A  20       8.825 -21.476  -5.940  1.00  3.77           N
ATOM     92  CA  THR A  20       9.478 -22.432  -5.047  1.00  3.70           C
ATOM     93  C   THR A  20       9.444 -23.827  -5.640  1.00  3.56           C
ATOM     94  O   THR A  20       8.383 -24.281  -6.108  1.00  4.14           O
ATOM     95  CB  THR A  20       8.787 -22.430  -3.673  1.00  4.76           C
ATOM     96  N   THR A  21      10.560 -24.542  -5.569  1.00  4.00           N
ATOM     97  CA  THR A  21      10.597 -25.962  -5.876  1.00  4.05           C
ATOM     98  C   THR A  21      10.984 -26.770  -4.636  1.00  4.53           C
ATOM     99  O   THR A  21      11.770 -26.361  -3.802  1.00  5.04           O
ATOM    100  CB  THR A  21      11.488 -26.293  -7.083  1.00  4.38           C
ATOM    189  N   GLN A  40       0.280  -6.099  -9.049  1.00  6.35           N
ATOM    190  CA  GLN A  40       0.087  -7.454  -9.580  1.00  6.35           C
ATOM    191  C   GLN A  40       0.964  -8.417  -8.788  1.00  6.09           C
ATOM    192  O   GLN A  40       2.080  -8.093  -8.393  1.00  6.87           O
ATOM    193  CB  GLN A  40       0.461  -7.523 -11.060  1.00  7.52           C
ATOM    194  N   THR A  41       0.419  -9.596  -8.544  1.00  6.66           N
ATOM    195  CA  THR A  41       1.108 -10.640  -7.800  1.00  6.93           C
ATOM    196  C   THR A  41       0.932 -12.005  -8.414  1.00  6.82           C
ATOM    197  O   THR A  41      -0.069 -12.258  -9.104  1.00  8.79           O
ATOM    198  CB  THR A  41       0.633 -10.636  -6.352  1.00 10.84           C
ATOM    199  N   ALA A  42       1.951 -12.847  -8.263  1.00  6.44           N
ATOM    200  CA  ALA A  42       1.923 -14.209  -8.797  1.00  6.59           C
ATOM    201  C   ALA A  42       2.829 -15.117  -7.992  1.00  5.51           C
ATOM    202  O   ALA A  42       3.835 -14.684  -7.420  1.00  5.94           O
ATOM    203  CB  ALA A  42       2.327 -14.218 -10.264  1.00  9.02           C
ATOM    204  N   LYS A  43       2.479 -16.398  -7.978  1.00  6.26           N
ATOM    205  CA  LYS A  43       3.247 -17.395  -7.256  1.00  6.48           C
ATOM    206  C   LYS A  43       3.186 -18.741  -7.955  1.00  5.78           C
ATOM    207  O   LYS A  43       2.206 -19.041  -8.623  1.00  9.40           O
ATOM    208  CB  LYS A  43       2.727 -17.535  -5.836  1.00  8.81           C
ATOM    209  N   SER A  44       4.243 -19.534  -7.818  1.00  4.43           N
ATOM    210  CA  SER A  44       4.241 -20.890  -8.325  1.00  4.28           C
ATOM    211  C   SER A  44       4.998 -21.811  -7.358  1.00  4.09           C
ATOM    212  O   SER A  44       5.865 -21.377  -6.584  1.00  4.53           O
ATOM    213  CB  SER A  44       4.831 -20.949  -9.731  1.00  5.33           C
ATOM    214  N   PHE A  45       4.660 -23.091  -7.444  1.00  4.39           N
ATOM    215  CA  PHE A  45       5.198 -24.135  -6.576  1.00  4.44           C
ATOM    216  C   PHE A  45       5.222 -25.415  -7.389  1.00  4.16           C
ATOM    217  O   PHE A  45       4.183 -25.754  -7.979  1.00  5.11           O
ATOM    218  CB  PHE A  45       4.254 -24.281  -5.370  1.00  5.22           C
ATOM    219  N   ALA A  46       6.347 -26.119  -7.403  1.00  3.62           N
ATOM    220  CA  ALA A  46       6.443 -27.338  -8.202  1.00  3.99           C
ATOM    221  C   ALA A  46       7.579 -28.205  -7.717  1.00  4.57           C
ATOM    222  O   ALA A  46       8.479 -27.750  -7.000  1.00  4.79           O
ATOM    223  CB  ALA A  46       6.607 -27.026  -9.678  1.00  4.52           C
TER
END
""").construct_hierarchy()

  pdb_par_ac_h = iotbx.pdb.input(source_info=None, lines = """\
CRYST1   46.460   46.460  193.210  90.00  90.00 120.00 P 31 2 1
SCALE1      0.021524  0.012427  0.000000        0.00000
SCALE2      0.000000  0.024854  0.000000        0.00000
SCALE3      0.000000  0.000000  0.005176        0.00000
ATOM      1  N  AALA A  15       5.011  -5.031  -8.967  0.50  5.73           N
ATOM      2  CA AALA A  15       4.943  -6.455  -9.287  0.50  6.16           C
ATOM      3  C  AALA A  15       5.610  -7.252  -8.166  0.50  6.39           C
ATOM      4  O  AALA A  15       6.751  -7.007  -7.857  0.50  9.87           O
ATOM      5  CB AALA A  15       5.684  -6.739 -10.604  0.50  7.46           C
ATOM      6  N  BALA A  15       5.111  -5.031  -8.967  0.50  5.73           N
ATOM      7  CA BALA A  15       5.043  -6.455  -9.287  0.50  6.16           C
ATOM      8  C  BALA A  15       5.710  -7.252  -8.166  0.50  6.39           C
ATOM      9  O  BALA A  15       6.851  -7.007  -7.857  0.50  9.87           O
ATOM     10  CB BALA A  15       5.784  -6.739 -10.604  0.50  7.46           C
ATOM     11  N  ATHR A  16       4.929  -8.263  -7.636  0.50  5.49           N
ATOM     12  C  ATHR A  16       5.316 -10.600  -7.084  0.50  5.07           C
ATOM     13  O  ATHR A  16       4.214 -11.002  -7.422  0.50  6.51           O
ATOM     14  CA ATHR A  16       5.513  -9.172  -6.657  0.50  5.70           C
ATOM     15  CB ATHR A  16       4.864  -9.001  -5.276  0.50  8.31           C
ATOM     16  N  BTHR A  16       5.029  -8.263  -7.636  0.50  5.49           N
ATOM     17  C  BTHR A  16       5.416 -10.600  -7.084  0.50  5.07           C
ATOM     18  O  BTHR A  16       4.314 -11.002  -7.422  0.50  6.51           O
ATOM     19  CA BTHR A  16       5.613  -9.172  -6.657  0.50  5.70           C
ATOM     20  CB BTHR A  16       4.964  -9.001  -5.276  0.50  8.31           C
ATOM     21  N   GLY A  17       6.393 -11.375  -7.067  1.00  4.56           N
ATOM     22  CA  GLY A  17       6.325 -12.770  -7.439  1.00  4.26           C
ATOM     23  C   GLY A  17       7.219 -13.654  -6.621  1.00  4.41           C
ATOM     24  O   GLY A  17       8.263 -13.233  -6.114  1.00  5.01           O
ATOM     25  N   SER A  18       6.827 -14.921  -6.561  1.00  4.24           N
ATOM     26  CA  SER A  18       7.657 -15.945  -5.959  1.00  4.02           C
ATOM     27  C   SER A  18       7.539 -17.244  -6.724  1.00  3.64           C
ATOM     28  O   SER A  18       6.482 -17.526  -7.331  1.00  4.07           O
ATOM     29  CB  SER A  18       7.335 -16.157  -4.481  1.00  5.86           C
ATOM     30  N   ALA A  19       8.573 -18.049  -6.627  1.00  3.19           N
ATOM     31  CA  ALA A  19       8.578 -19.414  -7.195  1.00  3.31           C
ATOM     32  C   ALA A  19       9.370 -20.307  -6.238  1.00  3.20           C
ATOM     33  O   ALA A  19      10.484 -19.956  -5.817  1.00  4.21           O
ATOM     34  CB  ALA A  19       9.235 -19.415  -8.574  1.00  3.85           C
ATOM     35  N   THR A  20       8.825 -21.476  -5.940  1.00  3.77           N
ATOM     36  CA  THR A  20       9.478 -22.432  -5.047  1.00  3.70           C
ATOM     37  C   THR A  20       9.444 -23.827  -5.640  1.00  3.56           C
ATOM     38  O   THR A  20       8.383 -24.281  -6.108  1.00  4.14           O
ATOM     39  CB  THR A  20       8.787 -22.430  -3.673  1.00  4.76           C
ATOM     40  N   THR A  21      10.560 -24.542  -5.569  1.00  4.00           N
ATOM     41  CA  THR A  21      10.597 -25.962  -5.876  1.00  4.05           C
ATOM     42  C   THR A  21      10.984 -26.770  -4.636  1.00  4.53           C
ATOM     43  O   THR A  21      11.770 -26.361  -3.802  1.00  5.04           O
ATOM     44  CB  THR A  21      11.488 -26.293  -7.083  1.00  4.38           C
ATOM     45  N   GLN A  40       0.280  -6.099  -9.049  1.00  6.35           N
ATOM     46  CA  GLN A  40       0.087  -7.454  -9.580  1.00  6.35           C
ATOM     47  C   GLN A  40       0.964  -8.417  -8.788  1.00  6.09           C
ATOM     48  O   GLN A  40       2.080  -8.093  -8.393  1.00  6.87           O
ATOM     49  CB  GLN A  40       0.461  -7.523 -11.060  1.00  7.52           C
ATOM     50  N   THR A  41       0.419  -9.596  -8.544  1.00  6.66           N
ATOM     51  CA  THR A  41       1.108 -10.640  -7.800  1.00  6.93           C
ATOM     52  C   THR A  41       0.932 -12.005  -8.414  1.00  6.82           C
ATOM     53  O   THR A  41      -0.069 -12.258  -9.104  1.00  8.79           O
ATOM     54  CB  THR A  41       0.633 -10.636  -6.352  1.00 10.84           C
ATOM     55  N   ALA A  42       1.951 -12.847  -8.263  1.00  6.44           N
ATOM     56  CA  ALA A  42       1.923 -14.209  -8.797  1.00  6.59           C
ATOM     57  C   ALA A  42       2.829 -15.117  -7.992  1.00  5.51           C
ATOM     58  O   ALA A  42       3.835 -14.684  -7.420  1.00  5.94           O
ATOM     59  CB  ALA A  42       2.327 -14.218 -10.264  1.00  9.02           C
ATOM     60  N   LYS A  43       2.479 -16.398  -7.978  1.00  6.26           N
ATOM     61  CA  LYS A  43       3.247 -17.395  -7.256  1.00  6.48           C
ATOM     62  C   LYS A  43       3.186 -18.741  -7.955  1.00  5.78           C
ATOM     63  O   LYS A  43       2.206 -19.041  -8.623  1.00  9.40           O
ATOM     64  CB  LYS A  43       2.727 -17.535  -5.836  1.00  8.81           C
ATOM     65  N   SER A  44       4.243 -19.534  -7.818  1.00  4.43           N
ATOM     66  CA  SER A  44       4.241 -20.890  -8.325  1.00  4.28           C
ATOM     67  C   SER A  44       4.998 -21.811  -7.358  1.00  4.09           C
ATOM     68  O   SER A  44       5.865 -21.377  -6.584  1.00  4.53           O
ATOM     69  CB  SER A  44       4.831 -20.949  -9.731  1.00  5.33           C
ATOM     70  N   PHE A  45       4.660 -23.091  -7.444  1.00  4.39           N
ATOM     71  CA  PHE A  45       5.198 -24.135  -6.576  1.00  4.44           C
ATOM     72  C   PHE A  45       5.222 -25.415  -7.389  1.00  4.16           C
ATOM     73  O   PHE A  45       4.183 -25.754  -7.979  1.00  5.11           O
ATOM     74  CB  PHE A  45       4.254 -24.281  -5.370  1.00  5.22           C
ATOM     75  N   ALA A  46       6.347 -26.119  -7.403  1.00  3.62           N
ATOM     76  CA  ALA A  46       6.443 -27.338  -8.202  1.00  3.99           C
ATOM     77  C   ALA A  46       7.579 -28.205  -7.717  1.00  4.57           C
ATOM     78  O   ALA A  46       8.479 -27.750  -7.000  1.00  4.79           O
ATOM     79  CB  ALA A  46       6.607 -27.026  -9.678  1.00  4.52           C
TER      80      ALA A  46
END
""").construct_hierarchy()


  pdb_apar2_h = iotbx.pdb.input(source_info=None, lines = """\
CRYST1  172.640  172.640  172.640  90.00  90.00  90.00 P 1
ATOM   2094  C   ASP N 271     109.854 129.638  88.152  1.00  0.00           C
ATOM   2095  CA  ASP N 271     108.849 130.718  87.818  1.00  0.00           C
ATOM   2096  CB  ASP N 271     109.517 131.848  87.048  1.00  0.00           C
ATOM   2097  CG  ASP N 271     108.534 132.901  86.590  1.00  0.00           C
ATOM   2098  N   ASP N 271     108.279 131.208  89.060  1.00  0.00           N
ATOM   2099  O   ASP N 271     110.669 129.804  89.047  1.00  0.00           O
ATOM   2100  OD1 ASP N 271     107.322 132.611  86.550  1.00  0.00           O
ATOM   2101  OD2 ASP N 271     108.973 134.031  86.298  1.00  0.00           O
ATOM   2102  C   LYS N 272     111.760 127.262  86.697  1.00  0.00           C
ATOM   2103  CA  LYS N 272     110.671 127.408  87.746  1.00  0.00           C
ATOM   2104  CB  LYS N 272     109.854 126.123  87.810  1.00  0.00           C
ATOM   2105  CD  LYS N 272     108.904 126.401  90.106  1.00  0.00           C
ATOM   2106  CE  LYS N 272     107.623 126.322  90.905  1.00  0.00           C
ATOM   2107  CG  LYS N 272     108.611 126.182  88.660  1.00  0.00           C
ATOM   2108  N   LYS N 272     109.816 128.535  87.424  1.00  0.00           N
ATOM   2109  NZ  LYS N 272     107.845 126.530  92.350  1.00  0.00           N
ATOM   2110  O   LYS N 272     111.665 126.402  85.830  1.00  0.00           O
ATOM   2111  C   ILE N 273     114.616 126.730  85.969  1.00  0.00           C
ATOM   2112  CA  ILE N 273     113.831 127.992  85.732  1.00  0.00           C
ATOM   2113  CB  ILE N 273     114.800 129.181  85.738  1.00  0.00           C
ATOM   2114  CD1 ILE N 273     116.646 130.274  87.071  1.00  0.00           C
ATOM   2115  CG1 ILE N 273     115.554 129.250  87.064  1.00  0.00           C
ATOM   2116  CG2 ILE N 273     114.064 130.470  85.480  1.00  0.00           C
ATOM   2117  N   ILE N 273     112.790 128.090  86.749  1.00  0.00           N
ATOM   2118  O   ILE N 273     114.607 126.208  87.067  1.00  0.00           O
ATOM   2119  C   SER N 274     117.460 125.414  84.290  1.00  0.00           C
ATOM   2120  CA  SER N 274     116.236 125.152  85.140  1.00  0.00           C
ATOM   2121  CB  SER N 274     115.642 123.795  84.795  1.00  0.00           C
ATOM   2122  N   SER N 274     115.266 126.216  84.935  1.00  0.00           N
ATOM   2123  O   SER N 274     117.336 125.793  83.133  1.00  0.00           O
ATOM   2124  OG  SER N 274     115.345 123.724  83.415  1.00  0.00           O
ATOM   2125  C   PHE N 275     120.234 124.460  83.057  1.00  0.00           C
ATOM   2126  CA  PHE N 275     119.856 125.505  84.086  1.00  0.00           C
ATOM   2127  CB  PHE N 275     121.036 125.723  85.016  1.00  0.00           C
ATOM   2128  CD1 PHE N 275     120.424 128.127  85.038  1.00  0.00           C
ATOM   2129  CD2 PHE N 275     121.668 127.271  86.864  1.00  0.00           C
ATOM   2130  CE1 PHE N 275     120.431 129.366  85.608  1.00  0.00           C
ATOM   2131  CE2 PHE N 275     121.671 128.508  87.443  1.00  0.00           C
ATOM   2132  CG  PHE N 275     121.048 127.068  85.652  1.00  0.00           C
ATOM   2133  CZ  PHE N 275     121.052 129.559  86.815  1.00  0.00           C
ATOM   2134  N   PHE N 275     118.643 125.218  84.853  1.00  0.00           N
ATOM   2135  O   PHE N 275     119.942 123.279  83.207  1.00  0.00           O
ATOM   2136  C   ALA N 276     122.598 123.230  81.591  1.00  0.00           C
ATOM   2137  CA  ALA N 276     121.463 124.047  80.996  1.00  0.00           C
ATOM   2138  CB  ALA N 276     121.959 124.853  79.811  1.00  0.00           C
ATOM   2139  N   ALA N 276     120.896 124.925  82.004  1.00  0.00           N
ATOM   2140  O   ALA N 276     123.300 123.693  82.488  1.00  0.00           O
ATOM   2141  C   GLY N 277     123.870 120.733  82.942  1.00  0.00           C
ATOM   2142  CA  GLY N 277     123.899 121.188  81.499  1.00  0.00           C
ATOM   2143  N   GLY N 277     122.761 122.003  81.113  1.00  0.00           N
ATOM   2144  O   GLY N 277     123.935 119.542  83.232  1.00  0.00           O
ATOM   2145  C   VAL N 278     122.701 120.500  85.709  1.00  0.00           C
ATOM   2146  CA  VAL N 278     123.826 121.425  85.269  1.00  0.00           C
ATOM   2147  CB  VAL N 278     123.748 122.729  86.052  1.00  0.00           C
ATOM   2148  CG1 VAL N 278     123.758 122.461  87.531  1.00  0.00           C
ATOM   2149  CG2 VAL N 278     124.889 123.634  85.653  1.00  0.00           C
ATOM   2150  N   VAL N 278     123.766 121.693  83.847  1.00  0.00           N
ATOM   2151  O   VAL N 278     121.576 120.625  85.240  1.00  0.00           O
ATOM   2152  C   LYS N 279     122.293 118.200  88.623  1.00  0.00           C
ATOM   2153  CA  LYS N 279     122.266 118.434  87.093  1.00  0.00           C
ATOM   2154  CB  LYS N 279     122.696 117.173  86.350  1.00  0.00           C
ATOM   2155  CD  LYS N 279     122.717 115.904  84.174  1.00  0.00           C
ATOM   2156  CE  LYS N 279     121.718 114.856  84.636  1.00  0.00           C
ATOM   2157  CG  LYS N 279     122.463 117.242  84.850  1.00  0.00           C
ATOM   2158  N   LYS N 279     123.053 119.597  86.639  1.00  0.00           N
ATOM   2159  NZ  LYS N 279     121.960 113.522  84.020  1.00  0.00           N
ATOM   2160  O   LYS N 279     123.105 118.821  89.310  1.00  0.00           O
ATOM   2161  C   PHE N 280     122.110 115.862  91.064  1.00  0.00           C
ATOM   2162  CA  PHE N 280     121.227 117.043  90.584  1.00  0.00           C
ATOM   2163  CB  PHE N 280     119.768 116.778  90.931  1.00  0.00           C
ATOM   2164  CD1 PHE N 280     118.146 118.596  91.497  1.00  0.00           C
ATOM   2165  CD2 PHE N 280     118.685 118.230  89.208  1.00  0.00           C
ATOM   2166  CE1 PHE N 280     117.295 119.619  91.132  1.00  0.00           C
ATOM   2167  CE2 PHE N 280     117.837 119.253  88.835  1.00  0.00           C
ATOM   2168  CG  PHE N 280     118.845 117.889  90.538  1.00  0.00           C
ATOM   2169  CZ  PHE N 280     117.141 119.950  89.799  1.00  0.00           C
ATOM   2170  N   PHE N 280     121.378 117.359  89.139  1.00  0.00           N
ATOM   2171  O   PHE N 280     122.105 114.774  90.467  1.00  0.00           O
ATOM   2602  C   ALA N 339     120.464 128.840  79.943  1.00  0.00           C
ATOM   2603  CA  ALA N 339     121.554 129.338  79.013  1.00  0.00           C
ATOM   2604  CB  ALA N 339     121.199 129.029  77.578  1.00  0.00           C
ATOM   2605  N   ALA N 339     122.841 128.757  79.352  1.00  0.00           N
ATOM   2606  O   ALA N 339     120.396 127.658  80.263  1.00  0.00           O
ATOM   2607  C   VAL N 340     117.492 128.599  80.322  1.00  0.00           C
ATOM   2608  CA  VAL N 340     118.447 129.411  81.169  1.00  0.00           C
ATOM   2609  CB  VAL N 340     117.734 130.666  81.702  1.00  0.00           C
ATOM   2610  CG1 VAL N 340     116.444 130.309  82.399  1.00  0.00           C
ATOM   2611  CG2 VAL N 340     118.643 131.413  82.642  1.00  0.00           C
ATOM   2612  N   VAL N 340     119.609 129.757  80.376  1.00  0.00           N
ATOM   2613  O   VAL N 340     117.373 128.821  79.119  1.00  0.00           O
ATOM   2614  C   ASN N 341     114.499 127.022  81.166  1.00  0.00           C
ATOM   2615  CA  ASN N 341     115.761 126.908  80.322  1.00  0.00           C
ATOM   2616  CB  ASN N 341     116.146 125.449  80.130  1.00  0.00           C
ATOM   2617  CG  ASN N 341     115.213 124.731  79.202  1.00  0.00           C
ATOM   2618  N   ASN N 341     116.838 127.636  80.954  1.00  0.00           N
ATOM   2619  ND2 ASN N 341     114.996 123.450  79.455  1.00  0.00           N
ATOM   2620  O   ASN N 341     114.354 126.315  82.161  1.00  0.00           O
ATOM   2621  OD1 ASN N 341     114.676 125.323  78.268  1.00  0.00           O
ATOM   2622  C   ILE N 342     111.425 126.940  81.377  1.00  0.00           C
ATOM   2623  CA  ILE N 342     112.363 128.130  81.513  1.00  0.00           C
ATOM   2624  CB  ILE N 342     111.639 129.408  81.044  1.00  0.00           C
ATOM   2625  CD1 ILE N 342     111.812 131.926  81.022  1.00  0.00           C
ATOM   2626  CG1 ILE N 342     112.389 130.643  81.531  1.00  0.00           C
ATOM   2627  CG2 ILE N 342     110.218 129.459  81.550  1.00  0.00           C
ATOM   2628  N   ILE N 342     113.615 127.943  80.794  1.00  0.00           N
ATOM   2629  O   ILE N 342     111.219 126.417  80.289  1.00  0.00           O
ATOM   2630  C   LEU N 343     108.457 125.498  82.528  1.00  0.00           C
ATOM   2631  CA  LEU N 343     109.977 125.339  82.513  1.00  0.00           C
ATOM   2632  CB  LEU N 343     110.373 124.459  83.691  1.00  0.00           C
ATOM   2633  CD1 LEU N 343     112.060 123.086  84.901  1.00  0.00           C
ATOM   2634  CD2 LEU N 343     111.765 122.813  82.425  1.00  0.00           C
ATOM   2635  CG  LEU N 343     111.736 123.787  83.591  1.00  0.00           C
ATOM   2636  N   LEU N 343     110.810 126.539  82.484  1.00  0.00           N
ATOM   2637  O   LEU N 343     107.761 124.516  82.764  1.00  0.00           O
END
""").construct_hierarchy()

  s_apar_records1 = """\
SHEET    1   A 2 TYR A   2  VAL A   7  0
SHEET    2   A 2 MET A  31  SER A  36 -1  O  ALA A  36   N  ALA A   2
"""
  s_apar_records2 = """\
SHEET    1   A 2 TYR A   2  VAL A   7  0
SHEET    2   A 2 MET A  31  SER A  36 -1  O  ALA A  32   N  ALA A   6
"""

  s_par_records1 = """\
SHEET    1   B 2 ALA A  15  THR A  21  0
SHEET    2   B 2 GLN A  40  ALA A  46  1  O  GLN A  40   N  THR A  16
"""
  s_par_records2 = """\
SHEET    1   B 2 ALA A  15  THR A  21  0
SHEET    2   B 2 GLN A  40  ALA A  46  1  O  GLN A  44   N  THR A  20
"""
  s_par_records3 = """\
SHEET    1   B 2 ALA A  15  THR A  21  0
SHEET    2   B 2 GLN A  40  ALA A  46  1  N  GLN A  42   O  THR A  16
"""
  s_par_records4 = """\
SHEET    1   B 2 ALA A  15  THR A  21  0
SHEET    2   B 2 GLN A  40  ALA A  46  1  N  GLN A  46   O  THR A  20
"""
  s_par_records5 = """\
SHEET    1   B 2 GLN A  40  ALA A  46  0
SHEET    2   B 2 ALA A  15  THR A  21  1  O  THR A  20   N  ALA A  46
"""
  s_par_records6 = """\
SHEET    1   B 2 GLN A  40  ALA A  46  0
SHEET    2   B 2 ALA A  15  THR A  21  1  N  THR A  20   O  ALA A  44
"""
  s_par_records7 = """\
SHEET    1   B 2 GLN A  40  ALA A  46  0
SHEET    2   B 2 ALA A  15  THR A  21  1  N  THR A  16   O  GLN A  40
"""
  s_par_records8 = """\
SHEET    1   B 2 GLN A  40  ALA A  46  0
SHEET    2   B 2 ALA A  15  THR A  21  1  O  THR A  16   N  ALA A  42
"""

  s_apar_records3 = """\
SHEET    1   D 2 LYS N 272  LYS N 279  0
SHEET    2   D 2 VAL N 340  ILE N 342 -1  N  ASN N 341   O  SER N 274
"""

  s_par_records9 = """\
SHEET    1   B 2 SER A  44  ALA A  46  0
SHEET    2   B 2 ALA A  15  THR A  21  1  N  THR A  20   O  SER A  44
"""

  s_par_records10 = """\
SHEET    1   B 2 GLN A  40  ALA A  46  0
SHEET    2   B 2 SER A  18  THR A  21  1  N  THR A  20   O  SER A  44
"""

  s_par_records11 = """\
SHEET    1   B 2 GLN A  40  ALA A  46  0
SHEET    2   B 2 ALA A  19  THR A  21  1  N  THR A  20   O  SER A  44
"""
  s_par_records12 = """\
SHEET    1   B 2 GLN A  40  ALA A  46  0
SHEET    2   B 2 THR A  20  THR A  21  1  O  THR A  20   N  ALA A  46
"""



  log = null_out()
  # defpars = sec_str_master_phil
  n_hbonds = []
  n_hangles = []
  for pdb_h, recs in [
                        (pdb_apar_h, s_apar_records1),
                        (pdb_apar_h, s_apar_records2),
                        (pdb_par_h,  s_par_records1),
                        (pdb_par_h,  s_par_records2),
                        (pdb_par_h,  s_par_records3),
                        (pdb_par_h,  s_par_records4),
                        (pdb_par_h,  s_par_records5),
                        (pdb_par_h,  s_par_records6),
                        (pdb_par_h,  s_par_records7),
                        (pdb_par_h,  s_par_records8),
                        (pdb_apar2_h, s_apar_records3),
                        (pdb_par_h, s_par_records9),
                        (pdb_par_h, s_par_records10),
                        (pdb_par_h, s_par_records11),
                        (pdb_par_h, s_par_records12),
                        (pdb_par_ac_h,  s_par_records1), # 8 hbonds
                        ]:
    ioss_annotation = ioss.annotation.from_records(records = recs.split('\n'))
    ann = ioss_annotation.as_restraint_groups(prefix_scope="secondary_structure")
    defpars = iotbx.phil.parse(sec_str_master_phil_str)
    custom_pars = defpars.fetch(iotbx.phil.parse(ann))
    custom_pars_ex = custom_pars.extract()
    ss_manager = manager(
                pdb_h,
                sec_str_from_pdb_file=None,
                params=custom_pars_ex.secondary_structure,
                verbose=-1)
    proxies_for_grm, h_angles = ss_manager.create_protein_hbond_proxies(
      annotation= None,
      log          = log)
    # print proxies_for_grm.size()
    n_hbonds.append(proxies_for_grm.size())
    n_hangles.append(h_angles.size())
  print(n_hbonds, n_hangles)
  assert n_hbonds == [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 2, 3, 4, 2, 2, 8]
  assert n_hangles ==[15, 15, 18, 18, 18, 18, 18, 18, 18, 18, 6, 6, 9, 6, 3, 24]

def exercise_segid():
  if (not libtbx.env.has_module(name="ksdssp")):
    print("KSDSSP not available, skipping exercise_segid()")
    return
  pdb_par_segid_h = iotbx.pdb.input(source_info=None, lines =  """\
CRYST1   46.460   46.460  193.210  90.00  90.00 120.00 P 31 2 1
SCALE1      0.021524  0.012427  0.000000        0.00000
SCALE2      0.000000  0.024854  0.000000        0.00000
SCALE3      0.000000  0.000000  0.005176        0.00000
ATOM     67  N   ALA A  15       5.011  -5.031  -8.967  1.00  5.73           N
ATOM     68  CA  ALA A  15       4.943  -6.455  -9.287  1.00  6.16           C
ATOM     69  C   ALA A  15       5.610  -7.252  -8.166  1.00  6.39           C
ATOM     70  O   ALA A  15       6.751  -7.007  -7.857  1.00  9.87           O
ATOM     71  CB  ALA A  15       5.684  -6.739 -10.604  1.00  7.46           C
ATOM     72  N   THR A  16       4.929  -8.263  -7.636  1.00  5.49           N
ATOM     73  C   THR A  16       5.316 -10.600  -7.084  1.00  5.07           C
ATOM     74  O   THR A  16       4.214 -11.002  -7.422  1.00  6.51           O
ATOM     75  CA  THR A  16       5.513  -9.172  -6.657  1.00  5.70           C
ATOM     76  CB  THR A  16       4.864  -9.001  -5.276  1.00  8.31           C
ATOM     77  N   GLY A  17       6.393 -11.375  -7.067  1.00  4.56           N
ATOM     78  CA  GLY A  17       6.325 -12.770  -7.439  1.00  4.26           C
ATOM     79  C   GLY A  17       7.219 -13.654  -6.621  1.00  4.41           C
ATOM     80  O   GLY A  17       8.263 -13.233  -6.114  1.00  5.01           O
ATOM     81  N   SER A  18       6.827 -14.921  -6.561  1.00  4.24           N
ATOM     82  CA  SER A  18       7.657 -15.945  -5.959  1.00  4.02           C
ATOM     83  C   SER A  18       7.539 -17.244  -6.724  1.00  3.64           C
ATOM     84  O   SER A  18       6.482 -17.526  -7.331  1.00  4.07           O
ATOM     85  CB  SER A  18       7.335 -16.157  -4.481  1.00  5.86           C
ATOM     86  N   ALA A  19       8.573 -18.049  -6.627  1.00  3.19           N
ATOM     87  CA  ALA A  19       8.578 -19.414  -7.195  1.00  3.31           C
ATOM     88  C   ALA A  19       9.370 -20.307  -6.238  1.00  3.20           C
ATOM     89  O   ALA A  19      10.484 -19.956  -5.817  1.00  4.21           O
ATOM     90  CB  ALA A  19       9.235 -19.415  -8.574  1.00  3.85       seg C
ATOM     91  N   THR A  20       8.825 -21.476  -5.940  1.00  3.77       seg N
ATOM     92  CA  THR A  20       9.478 -22.432  -5.047  1.00  3.70       seg C
ATOM     93  C   THR A  20       9.444 -23.827  -5.640  1.00  3.56       seg C
ATOM     94  O   THR A  20       8.383 -24.281  -6.108  1.00  4.14       seg O
ATOM     95  CB  THR A  20       8.787 -22.430  -3.673  1.00  4.76       seg C
ATOM     96  N   THR A  21      10.560 -24.542  -5.569  1.00  4.00       seg N
ATOM     97  CA  THR A  21      10.597 -25.962  -5.876  1.00  4.05       seg C
ATOM     98  C   THR A  21      10.984 -26.770  -4.636  1.00  4.53       seg C
ATOM     99  O   THR A  21      11.770 -26.361  -3.802  1.00  5.04       seg O
ATOM    100  CB  THR A  21      11.488 -26.293  -7.083  1.00  4.38       seg C
ATOM    189  N   GLN A  40       0.280  -6.099  -9.049  1.00  6.35       seg N
ATOM    190  CA  GLN A  40       0.087  -7.454  -9.580  1.00  6.35       seg C
ATOM    191  C   GLN A  40       0.964  -8.417  -8.788  1.00  6.09       seg C
ATOM    192  O   GLN A  40       2.080  -8.093  -8.393  1.00  6.87       seg O
ATOM    193  CB  GLN A  40       0.461  -7.523 -11.060  1.00  7.52       seg C
ATOM    194  N   THR A  41       0.419  -9.596  -8.544  1.00  6.66       seg N
ATOM    195  CA  THR A  41       1.108 -10.640  -7.800  1.00  6.93       seg C
ATOM    196  C   THR A  41       0.932 -12.005  -8.414  1.00  6.82       seg C
ATOM    197  O   THR A  41      -0.069 -12.258  -9.104  1.00  8.79       seg O
ATOM    198  CB  THR A  41       0.633 -10.636  -6.352  1.00 10.84       seg C
ATOM    199  N   ALA A  42       1.951 -12.847  -8.263  1.00  6.44       seg N
ATOM    200  CA  ALA A  42       1.923 -14.209  -8.797  1.00  6.59       seg C
ATOM    201  C   ALA A  42       2.829 -15.117  -7.992  1.00  5.51       seg C
ATOM    202  O   ALA A  42       3.835 -14.684  -7.420  1.00  5.94       seg O
ATOM    203  CB  ALA A  42       2.327 -14.218 -10.264  1.00  9.02       seg C
ATOM    204  N   LYS A  43       2.479 -16.398  -7.978  1.00  6.26       seg N
ATOM    205  CA  LYS A  43       3.247 -17.395  -7.256  1.00  6.48       seg C
ATOM    206  C   LYS A  43       3.186 -18.741  -7.955  1.00  5.78       seg C
ATOM    207  O   LYS A  43       2.206 -19.041  -8.623  1.00  9.40       seg O
ATOM    208  CB  LYS A  43       2.727 -17.535  -5.836  1.00  8.81       seg C
ATOM    209  N   SER A  44       4.243 -19.534  -7.818  1.00  4.43       seg N
ATOM    210  CA  SER A  44       4.241 -20.890  -8.325  1.00  4.28       seg C
ATOM    211  C   SER A  44       4.998 -21.811  -7.358  1.00  4.09       seg C
ATOM    212  O   SER A  44       5.865 -21.377  -6.584  1.00  4.53       seg O
ATOM    213  CB  SER A  44       4.831 -20.949  -9.731  1.00  5.33       seg C
ATOM    214  N   PHE A  45       4.660 -23.091  -7.444  1.00  4.39       seg N
ATOM    215  CA  PHE A  45       5.198 -24.135  -6.576  1.00  4.44       seg C
ATOM    216  C   PHE A  45       5.222 -25.415  -7.389  1.00  4.16       seg C
ATOM    217  O   PHE A  45       4.183 -25.754  -7.979  1.00  5.11       seg O
ATOM    218  CB  PHE A  45       4.254 -24.281  -5.370  1.00  5.22       seg C
ATOM    219  N   ALA A  46       6.347 -26.119  -7.403  1.00  3.62       seg N
ATOM    220  CA  ALA A  46       6.443 -27.338  -8.202  1.00  3.99       seg C
ATOM    221  C   ALA A  46       7.579 -28.205  -7.717  1.00  4.57       seg C
ATOM    222  O   ALA A  46       8.479 -27.750  -7.000  1.00  4.79       seg O
ATOM    223  CB  ALA A  46       6.607 -27.026  -9.678  1.00  4.52       seg C
TER
END
""").construct_hierarchy()
  log = null_out()
  defpars = iotbx.phil.parse(sec_str_master_phil_str)
  custom_pars = defpars.extract()
  custom_pars.secondary_structure.enabled=True
  # custom_pars = defpars.fetch(iotbx.phil.parse("secondary_structure_restraints=True"))

  ss_manager = manager(
              pdb_par_segid_h,
              sec_str_from_pdb_file=None,
              params=custom_pars.secondary_structure,
              verbose=-1)
  ss_manager.initialize()
  proxies_for_grm, h_angles = ss_manager.create_protein_hbond_proxies(
    annotation = None,
    log          = log)
  assert proxies_for_grm.size() == 4
  assert h_angles.size() == 12

def exercise_phil_generation():
  log = null_out()
  defpars = sec_str_master_phil.fetch()

if __name__ == "__main__" :
  exercise_protein()
  exercise_sheet_ends()

  exercise_helix_bonding_pattern()
  exercise_sheets_bonding_pattern()

  # exercise_segid() # appeared to be nonvalid test... now when we don't ignore
  # segid in SS, we cannot construct SHEET situated in two different segids...
  print("OK")
