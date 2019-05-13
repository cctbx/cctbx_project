from __future__ import division
from __future__ import print_function
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

pdb_for_map = pdb_poor

pdb_answer = pdb_poor

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
    resolution_factor = resolution_factor,
    pdb_for_map       = pdb_for_map)
  #
  assert count_outliers(t.ph_answer) == 3
  #
  result = mmtbx.refinement.real_space.fit_residues.run(
    pdb_hierarchy     = t.ph_poor,
    vdw_radii         = t.vdw,
    crystal_symmetry  = t.crystal_symmetry,
    map_data          = t.target_map,
    do_all            = True,
    massage_map       = False,
    backbone_sample   = False,
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
    do_all            = True,
    massage_map       = False,
    backbone_sample   = False,
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
