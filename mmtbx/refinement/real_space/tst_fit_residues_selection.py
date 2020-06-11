from __future__ import absolute_import, division, print_function
import time
import mmtbx.refinement.real_space.fit_residues
import mmtbx.refinement.real_space
from mmtbx.rotamer.rotamer_eval import RotamerEval
rotamer_eval = RotamerEval()

pdb_poor = """\
CRYST1   36.039   33.133   33.907  90.00  90.00  90.00 P 1
ATOM    369  N   GLN A  23      12.323  18.703  12.553  1.00  0.00           N
ATOM    370  CA  GLN A  23      12.181  19.413  11.273  1.00  0.00           C
ATOM    371  C   GLN A  23      12.652  18.596  10.060  1.00  0.00           C
ATOM    372  O   GLN A  23      12.591  19.108   8.942  1.00  0.00           O
ATOM    373  CB  GLN A  23      12.912  20.768  11.283  1.00  0.00           C
ATOM    374  CG  GLN A  23      12.633  21.683  12.491  1.00  0.00           C
ATOM    375  CD  GLN A  23      11.920  22.989  12.156  1.00  0.00           C
ATOM    376  OE1 GLN A  23      10.721  23.090  12.304  1.00  0.00           O
ATOM    377  NE2 GLN A  23      12.465  23.786  11.261  1.00  0.00           N
ATOM    378  H   GLN A  23      13.024  19.067  13.190  1.00  0.00           H
ATOM    379  HA  GLN A  23      11.130  19.627  11.113  1.00  0.00           H
ATOM    380  HB2 GLN A  23      13.985  20.581  11.232  1.00  0.00           H
ATOM    381  HB3 GLN A  23      12.646  21.297  10.366  1.00  0.00           H
ATOM    382  HG2 GLN A  23      12.047  21.165  13.248  1.00  0.00           H
ATOM    383  HG3 GLN A  23      13.584  21.937  12.941  1.00  0.00           H
ATOM    384 HE21 GLN A  23      13.373  23.615  10.881  1.00  0.00           H
ATOM    385 HE22 GLN A  23      11.850  24.502  10.909  1.00  0.00           H
ATOM    386  N   LYS A  24      13.104  17.339  10.235  1.00  0.00           N
ATOM    387  CA  LYS A  24      13.177  16.380   9.120  1.00  0.00           C
ATOM    388  C   LYS A  24      12.362  15.085   9.207  1.00  0.00           C
ATOM    389  O   LYS A  24      12.509  14.250   8.322  1.00  0.00           O
ATOM    390  CB  LYS A  24      14.585  16.121   8.540  1.00  0.00           C
ATOM    391  CG  LYS A  24      15.222  17.335   7.850  1.00  0.00           C
ATOM    392  CD  LYS A  24      16.336  16.882   6.887  1.00  0.00           C
ATOM    393  CE  LYS A  24      16.513  17.834   5.699  1.00  0.00           C
ATOM    394  NZ  LYS A  24      16.975  19.185   6.105  1.00  0.00           N
ATOM    395  H   LYS A  24      13.162  16.988  11.180  1.00  0.00           H
ATOM    396  HA  LYS A  24      12.652  16.912   8.352  1.00  0.00           H
ATOM    397  HB2 LYS A  24      15.242  15.676   9.281  1.00  0.00           H
ATOM    398  HB3 LYS A  24      14.484  15.365   7.765  1.00  0.00           H
ATOM    399  HG2 LYS A  24      14.449  17.841   7.271  1.00  0.00           H
ATOM    400  HG3 LYS A  24      15.616  18.028   8.593  1.00  0.00           H
ATOM    401  HD2 LYS A  24      17.272  16.766   7.432  1.00  0.00           H
ATOM    402  HD3 LYS A  24      16.079  15.911   6.457  1.00  0.00           H
ATOM    403  HE2 LYS A  24      17.224  17.384   5.000  1.00  0.00           H
ATOM    404  HE3 LYS A  24      15.550  17.907   5.179  1.00  0.00           H
ATOM    405  HZ1 LYS A  24      17.842  19.131   6.619  1.00  0.00           H
ATOM    406  HZ2 LYS A  24      17.127  19.760   5.284  1.00  0.00           H
ATOM    407  HZ3 LYS A  24      16.271  19.638   6.677  1.00  0.00           H
ATOM    408  N   ILE A  25      11.585  14.851  10.259  1.00  0.00           N
ATOM    409  CA  ILE A  25      10.847  13.577  10.386  1.00  0.00           C
ATOM    410  C   ILE A  25       9.429  13.760  10.915  1.00  0.00           C
ATOM    411  O   ILE A  25       8.548  12.967  10.595  1.00  0.00           O
ATOM    412  CB  ILE A  25      11.688  12.555  11.204  1.00  0.00           C
ATOM    413  CG1 ILE A  25      12.819  11.934  10.342  1.00  0.00           C
ATOM    414  CG2 ILE A  25      10.821  11.423  11.770  1.00  0.00           C
ATOM    415  CD1 ILE A  25      14.186  12.604  10.544  1.00  0.00           C
ATOM    416  H   ILE A  25      11.460  15.580  10.945  1.00  0.00           H
ATOM    417  HA  ILE A  25      10.686  13.155   9.391  1.00  0.00           H
ATOM    418  HB  ILE A  25      12.120  13.047  12.071  1.00  0.00           H
ATOM    419 HG12 ILE A  25      12.951  10.880  10.585  1.00  0.00           H
ATOM    420 HG13 ILE A  25      12.546  11.962   9.288  1.00  0.00           H
ATOM    421 HG21 ILE A  25      10.166  11.080  10.978  1.00  0.00           H
ATOM    422 HG22 ILE A  25      11.417  10.590  12.132  1.00  0.00           H
ATOM    423 HG23 ILE A  25      10.183  11.789  12.572  1.00  0.00           H
ATOM    424 HD11 ILE A  25      14.120  13.685  10.454  1.00  0.00           H
ATOM    425 HD12 ILE A  25      14.564  12.366  11.533  1.00  0.00           H
ATOM    426 HD13 ILE A  25      14.886  12.233   9.796  1.00  0.00           H
ATOM    427  N   HIS A  26       9.285  14.608  11.927  1.00  0.00           N
ATOM    428  CA  HIS A  26       8.067  14.661  12.719  1.00  0.00           C
ATOM    429  C   HIS A  26       7.127  15.775  12.283  1.00  0.00           C
ATOM    430  O   HIS A  26       6.188  15.506  11.549  1.00  0.00           O
ATOM    431  CB  HIS A  26       8.452  14.761  14.176  1.00  0.00           C
ATOM    432  CG  HIS A  26       9.098  13.493  14.680  1.00  0.00           C
ATOM    433  ND1 HIS A  26       8.446  12.311  14.909  1.00  0.00           N
ATOM    434  CD2 HIS A  26      10.403  13.310  15.045  1.00  0.00           C
ATOM    435  CE1 HIS A  26       9.327  11.419  15.377  1.00  0.00           C
ATOM    436  NE2 HIS A  26      10.547  11.994  15.547  1.00  0.00           N
ATOM    437  H   HIS A  26      10.027  15.264  12.126  1.00  0.00           H
ATOM    438  HA  HIS A  26       7.503  13.734  12.601  1.00  0.00           H
ATOM    439  HB2 HIS A  26       9.108  15.619  14.347  1.00  0.00           H
ATOM    440  HB3 HIS A  26       7.511  14.924  14.690  1.00  0.00           H
ATOM    441  HD1 HIS A  26       7.449  12.149  14.789  1.00  0.00           H
ATOM    442  HD2 HIS A  26      11.175  14.064  14.986  1.00  0.00           H
ATOM    443  HE1 HIS A  26       9.055  10.402  15.640  1.00  0.00           H
ATOM    444  N   SER A  27       7.594  17.013  12.432  1.00  0.00           N
ATOM    445  CA  SER A  27       7.087  18.141  11.666  1.00  0.00           C
ATOM    446  C   SER A  27       8.277  18.930  11.142  1.00  0.00           C
ATOM    447  O   SER A  27       8.951  19.586  11.932  1.00  0.00           O
ATOM    448  CB  SER A  27       6.114  18.989  12.498  1.00  0.00           C
ATOM    449  OG  SER A  27       5.000  19.320  11.696  1.00  0.00           O
ATOM    450  H   SER A  27       8.440  17.164  12.960  1.00  0.00           H
ATOM    451  HA  SER A  27       6.539  17.768  10.798  1.00  0.00           H
ATOM    452  HB2 SER A  27       5.756  18.405  13.348  1.00  0.00           H
ATOM    453  HB3 SER A  27       6.592  19.896  12.873  1.00  0.00           H
ATOM    454  HG  SER A  27       5.170  20.114  11.183  1.00  0.00           H
"""

pdb_answer = """
CRYST1   36.039   33.133   33.907  90.00  90.00  90.00 P 1
SCALE1      0.027748  0.000000  0.000000        0.00000
SCALE2      0.000000  0.030181  0.000000        0.00000
SCALE3      0.000000  0.000000  0.029492        0.00000
ATOM      1  N   GLN A  23      12.323  18.703  12.553  1.00  0.00           N
ATOM      2  CA  GLN A  23      12.181  19.413  11.273  1.00  0.00           C
ATOM      3  C   GLN A  23      12.652  18.596  10.060  1.00  0.00           C
ATOM      4  O   GLN A  23      12.591  19.108   8.942  1.00  0.00           O
ATOM      5  CB  GLN A  23      12.912  20.768  11.283  1.00  0.00           C
ATOM      6  CG  GLN A  23      12.601  21.700  12.470  1.00  0.00           C
ATOM      7  CD  GLN A  23      11.817  22.958  12.110  1.00  0.00           C
ATOM      8  OE1 GLN A  23      10.614  22.913  11.968  1.00  0.00           O
ATOM      9  NE2 GLN A  23      12.443  23.921  11.465  1.00  0.00           N
ATOM     10  H   GLN A  23      13.024  19.067  13.190  1.00  0.00           H
ATOM     11  HA  GLN A  23      11.130  19.627  11.113  1.00  0.00           H
ATOM     12  HB2 GLN A  23      13.986  20.580  11.264  1.00  0.00           H
ATOM     13  HB3 GLN A  23      12.671  21.284  10.352  1.00  0.00           H
ATOM     14  HG2 GLN A  23      12.052  21.174  13.248  1.00  0.00           H
ATOM     15  HG3 GLN A  23      13.543  22.015  12.900  1.00  0.00           H
ATOM     16 HE21 GLN A  23      13.371  23.807  11.112  1.00  0.00           H
ATOM     17 HE22 GLN A  23      11.873  24.724  11.256  1.00  0.00           H
ATOM     18  N   LYS A  24      13.104  17.339  10.235  1.00  0.00           N
ATOM     19  CA  LYS A  24      13.177  16.380   9.120  1.00  0.00           C
ATOM     20  C   LYS A  24      12.362  15.085   9.207  1.00  0.00           C
ATOM     21  O   LYS A  24      12.509  14.250   8.322  1.00  0.00           O
ATOM     22  CB  LYS A  24      14.585  16.121   8.540  1.00  0.00           C
ATOM     23  CG  LYS A  24      15.232  17.342   7.872  1.00  0.00           C
ATOM     24  CD  LYS A  24      16.382  16.900   6.946  1.00  0.00           C
ATOM     25  CE  LYS A  24      16.502  17.777   5.695  1.00  0.00           C
ATOM     26  NZ  LYS A  24      16.826  19.191   6.009  1.00  0.00           N
ATOM     27  H   LYS A  24      13.162  16.988  11.180  1.00  0.00           H
ATOM     28  HA  LYS A  24      12.652  16.912   8.352  1.00  0.00           H
ATOM     29  HB2 LYS A  24      15.237  15.660   9.276  1.00  0.00           H
ATOM     30  HB3 LYS A  24      14.481  15.378   7.753  1.00  0.00           H
ATOM     31  HG2 LYS A  24      14.474  17.840   7.267  1.00  0.00           H
ATOM     32  HG3 LYS A  24      15.594  18.038   8.627  1.00  0.00           H
ATOM     33  HD2 LYS A  24      17.319  16.888   7.500  1.00  0.00           H
ATOM     34  HD3 LYS A  24      16.195  15.886   6.584  1.00  0.00           H
ATOM     35  HE2 LYS A  24      17.271  17.347   5.046  1.00  0.00           H
ATOM     36  HE3 LYS A  24      15.551  17.722   5.151  1.00  0.00           H
ATOM     37  HZ1 LYS A  24      16.876  19.746   5.167  1.00  0.00           H
ATOM     38  HZ2 LYS A  24      16.110  19.586   6.608  1.00  0.00           H
ATOM     39  HZ3 LYS A  24      17.710  19.246   6.503  1.00  0.00           H
ATOM     40  N   ILE A  25      11.585  14.851  10.259  1.00  0.00           N
ATOM     41  CA  ILE A  25      10.847  13.577  10.386  1.00  0.00           C
ATOM     42  C   ILE A  25       9.429  13.760  10.915  1.00  0.00           C
ATOM     43  O   ILE A  25       8.548  12.967  10.595  1.00  0.00           O
ATOM     44  CB  ILE A  25      11.688  12.555  11.204  1.00  0.00           C
ATOM     45  CG1 ILE A  25      12.838  11.957  10.351  1.00  0.00           C
ATOM     46  CG2 ILE A  25      10.826  11.407  11.744  1.00  0.00           C
ATOM     47  CD1 ILE A  25      14.168  12.711  10.496  1.00  0.00           C
ATOM     48  H   ILE A  25      11.460  15.580  10.945  1.00  0.00           H
ATOM     49  HA  ILE A  25      10.686  13.155   9.391  1.00  0.00           H
ATOM     50  HB  ILE A  25      12.102  13.042  12.083  1.00  0.00           H
ATOM     51 HG12 ILE A  25      13.033  10.925  10.641  1.00  0.00           H
ATOM     52 HG13 ILE A  25      12.547  11.919   9.302  1.00  0.00           H
ATOM     53 HG21 ILE A  25      11.483  10.708  12.248  1.00  0.00           H
ATOM     54 HG22 ILE A  25      10.067  11.754  12.439  1.00  0.00           H
ATOM     55 HG23 ILE A  25      10.352  10.863  10.930  1.00  0.00           H
ATOM     56 HD11 ILE A  25      14.948  12.270   9.880  1.00  0.00           H
ATOM     57 HD12 ILE A  25      14.035  13.744  10.191  1.00  0.00           H
ATOM     58 HD13 ILE A  25      14.489  12.685  11.537  1.00  0.00           H
ATOM     59  N   HIS A  26       9.285  14.608  11.927  1.00  0.00           N
ATOM     60  CA  HIS A  26       8.067  14.661  12.719  1.00  0.00           C
ATOM     61  C   HIS A  26       7.127  15.775  12.283  1.00  0.00           C
ATOM     62  O   HIS A  26       6.188  15.506  11.549  1.00  0.00           O
ATOM     63  CB  HIS A  26       8.452  14.761  14.176  1.00  0.00           C
ATOM     64  CG  HIS A  26       9.094  13.491  14.681  1.00  0.00           C
ATOM     65  ND1 HIS A  26       8.436  12.314  14.919  1.00  0.00           N
ATOM     66  CD2 HIS A  26      10.400  13.302  15.039  1.00  0.00           C
ATOM     67  CE1 HIS A  26       9.314  11.419  15.385  1.00  0.00           C
ATOM     68  NE2 HIS A  26      10.539  11.987  15.546  1.00  0.00           N
ATOM     69  H   HIS A  26      10.027  15.264  12.126  1.00  0.00           H
ATOM     70  HA  HIS A  26       7.503  13.734  12.601  1.00  0.00           H
ATOM     71  HB2 HIS A  26       9.111  15.617  14.346  1.00  0.00           H
ATOM     72  HB3 HIS A  26       7.511  14.928  14.690  1.00  0.00           H
ATOM     73  HD1 HIS A  26       7.437  12.157  14.805  1.00  0.00           H
ATOM     74  HD2 HIS A  26      11.176  14.051  14.973  1.00  0.00           H
ATOM     75  HE1 HIS A  26       9.038  10.404  15.654  1.00  0.00           H
ATOM     76  N   SER A  27       7.594  17.013  12.432  1.00  0.00           N
ATOM     77  CA  SER A  27       7.087  18.141  11.666  1.00  0.00           C
ATOM     78  C   SER A  27       8.277  18.930  11.142  1.00  0.00           C
ATOM     79  O   SER A  27       8.951  19.586  11.932  1.00  0.00           O
ATOM     80  CB  SER A  27       6.114  18.989  12.498  1.00  0.00           C
ATOM     81  OG  SER A  27       5.157  19.563  11.632  1.00  0.00           O
ATOM     82  H   SER A  27       8.440  17.164  12.960  1.00  0.00           H
ATOM     83  HA  SER A  27       6.539  17.768  10.798  1.00  0.00           H
ATOM     84  HB2 SER A  27       5.586  18.348  13.207  1.00  0.00           H
ATOM     85  HB3 SER A  27       6.637  19.770  13.054  1.00  0.00           H
ATOM     86  HG  SER A  27       4.531  20.102  12.122  1.00  0.00           H
TER
"""

pdb_for_map = pdb_answer


def count_outliers(h):
  n_outliers = 0
  for model in h.models():
    for chain in model.chains():
      for residue_group in chain.residue_groups():
        conformers = residue_group.conformers()
        for conformer in residue_group.conformers():
          residue = conformer.only_residue()
          fl = rotamer_eval.evaluate_residue(residue = residue)
          if fl.strip().upper()=="OUTLIER":
            n_outliers += 1
  return n_outliers

def exercise(i_pdb=0, d_min = 3, resolution_factor = 0.1):
  """
  Testing selections
  """
  t = mmtbx.refinement.real_space.setup_test(
    pdb_answer        = pdb_answer,
    pdb_poor          = pdb_poor,
    i_pdb             = i_pdb,
    d_min             = d_min,
    residues          = ["SER","HIS","ILE","LYS","GLN"],
    resolution_factor = resolution_factor,
    pdb_for_map       = pdb_for_map)
  #
  assert count_outliers(t.ph_poor) == 3
  assert count_outliers(t.ph_answer) == 0
  #
  result = mmtbx.refinement.real_space.fit_residues.run(
    pdb_hierarchy     = t.ph_poor,
    vdw_radii         = t.vdw,
    crystal_symmetry  = t.crystal_symmetry,
    map_data          = t.target_map,
    backbone_sample   = False,
    rotatable_hd      = t.rotatable_hd,
    rotamer_manager   = t.rotamer_manager,
    sin_cos_table     = t.sin_cos_table,
    mon_lib_srv       = t.mon_lib_srv)
  result.pdb_hierarchy.write_pdb_file(file_name = "refined_%s.pdb"%str(i_pdb),
    crystal_symmetry = t.crystal_symmetry)
  assert count_outliers(result.pdb_hierarchy) == 0

  # Now with selection
  t = mmtbx.refinement.real_space.setup_test(
    pdb_answer        = pdb_answer,
    pdb_poor          = pdb_poor,
    i_pdb             = i_pdb,
    d_min             = d_min,
    resolution_factor = resolution_factor,
    pdb_for_map       = pdb_for_map)

  bselection = t.ph_poor.atom_selection_cache().selection('resid 23:25')
  result = mmtbx.refinement.real_space.fit_residues.run(
    pdb_hierarchy     = t.ph_poor,
    vdw_radii         = t.vdw,
    crystal_symmetry  = t.crystal_symmetry,
    map_data          = t.target_map,
    bselection        = bselection,
    backbone_sample   = False,
    rotatable_hd      = t.rotatable_hd,
    rotamer_manager   = t.rotamer_manager,
    sin_cos_table     = t.sin_cos_table,
    mon_lib_srv       = t.mon_lib_srv)
  result.pdb_hierarchy.write_pdb_file(file_name = "refined_%s.pdb"%str(i_pdb),
    crystal_symmetry = t.crystal_symmetry)
  assert count_outliers(result.pdb_hierarchy) == 1


if(__name__ == "__main__"):
  t0 = time.time()
  exercise()
  print("Time: %6.4f"%(time.time()-t0))
