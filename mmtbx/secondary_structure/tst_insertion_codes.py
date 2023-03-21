
from __future__ import absolute_import, division, print_function
from mmtbx.secondary_structure import sec_str_master_phil_str, manager
import iotbx.pdb
import iotbx.pdb.secondary_structure as ioss
from libtbx.utils import null_out


def exercise_helix_bonding_pattern_with_insertions():
  alpha_h1_ends = iotbx.pdb.input(source_info=None, lines="""\
ATOM      1  N   ALA     1       1.643  -2.366  -1.408  1.00  0.00           N
ATOM      2  CA  ALA     1       1.280  -3.608  -2.069  1.00  0.00           C
ATOM      3  C   ALA     1      -0.114  -3.466  -2.684  1.00  0.00           C
ATOM      4  O   ALA     1      -0.327  -3.827  -3.840  1.00  0.00           O
ATOM      5  CB  ALA     1       1.361  -4.762  -1.068  1.00  0.00           C
ATOM      6  N   ALA     1A     -1.028  -2.938  -1.882  1.00  0.00           N
ATOM      7  CA  ALA     1A     -2.395  -2.743  -2.332  1.00  0.00           C
ATOM      8  C   ALA     1A     -2.396  -1.855  -3.579  1.00  0.00           C
ATOM      9  O   ALA     1A     -3.059  -2.167  -4.567  1.00  0.00           O
ATOM     10  CB  ALA     1A     -3.228  -2.150  -1.194  1.00  0.00           C
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
ATOM     46  N   ALA     9A     -5.617  -2.338 -12.812  1.00  0.00           N
ATOM     47  CA  ALA     9A     -5.879  -1.301 -13.795  1.00  0.00           C
ATOM     48  C   ALA     9A     -5.242  -1.695 -15.130  1.00  0.00           C
ATOM     49  O   ALA     9A     -5.880  -1.605 -16.177  1.00  0.00           O
ATOM     50  CB  ALA     9A     -5.358   0.040 -13.274  1.00  0.00           C
TER
""").construct_hierarchy()

  alpha_h1_ac = iotbx.pdb.input(source_info=None, lines="""\
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
ATOM     11  N   ALA     1A     -1.028  -2.938  -1.882  1.00  0.00           N
ATOM     12  CA  ALA     1A     -2.395  -2.743  -2.332  1.00  0.00           C
ATOM     13  C   ALA     1A     -2.396  -1.855  -3.579  1.00  0.00           C
ATOM     14  O   ALA     1A     -3.059  -2.167  -4.567  1.00  0.00           O
ATOM     15  CB  ALA     1A     -3.228  -2.150  -1.194  1.00  0.00           C
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
ATOM     56  N  AALA     9A     -5.617  -2.338 -12.812  0.50  0.00           N
ATOM     57  CA AALA     9A     -5.879  -1.301 -13.795  0.50  0.00           C
ATOM     58  C  AALA     9A     -5.242  -1.695 -15.130  0.50  0.00           C
ATOM     59  O  AALA     9A     -5.880  -1.605 -16.177  0.50  0.00           O
ATOM     60  CB AALA     9A     -5.358   0.040 -13.274  0.50  0.00           C
ATOM     62  N  BALA     9A     -5.517  -2.338 -12.812  0.50  0.00           N
ATOM     63  CA BALA     9A     -5.779  -1.301 -13.795  0.50  0.00           C
ATOM     64  C  BALA     9A     -5.142  -1.695 -15.130  0.50  0.00           C
ATOM     65  O  BALA     9A     -5.780  -1.605 -16.177  0.50  0.00           O
ATOM     66  CB BALA     9A     -5.258   0.040 -13.274  0.50  0.00           C
END
""").construct_hierarchy()

  alpha_annot_1_1 = """\
HELIX    1   1 ALA      1  ALA      9A 1                                  10
"""
  alpha_annot_1_2 = """\
HELIX    1   1 ALA      1A ALA      9  1                                  10
"""
  alpha_annot_1_3 = """\
HELIX    1   1 ALA      1A ALA      9A 1                                  10
"""
  alpha_annot_1_4 = """\
HELIX    1   1 ALA      1  ALA      9  1                                  10
"""

  alpha_annot_2 = """\
HELIX    1   1 ALA      2  ALA      5  1                                  10
"""

  log = null_out()
  n_hbonds = []
  n_hangles = []
  for pdb_h, recs in [
                        (alpha_h1_ends, alpha_annot_1_1), #6
                        (alpha_h1_ends, alpha_annot_1_2), #4
                        (alpha_h1_ends, alpha_annot_1_3), #5
                        (alpha_h1_ends, alpha_annot_1_4), #5
                        (alpha_h1_ac, alpha_annot_1_1), #9
                        (alpha_h1_ac, alpha_annot_1_2), #5
                        (alpha_h1_ac, alpha_annot_1_3), #7
                        (alpha_h1_ac, alpha_annot_1_4), #7
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
    proxies_for_grm, hangles = ss_manager.create_protein_hbond_proxies(
      annotation= None,
      log          = log)
    n_hbonds.append(proxies_for_grm.size())
    n_hangles.append(hangles.size())
  print(n_hbonds, n_hangles)
  assert n_hbonds ==  [6,   4,  5,  5,  9,  5,  7,  7]
  assert n_hangles == [18, 12, 15, 15, 28, 16, 22, 22]

def exercise_sheets_bonding_pattern_with_insertions():
  pdb_apar_h = iotbx.pdb.input(source_info=None, lines = """\
SCRYST1   46.460   46.460  193.210  90.00  90.00 120.00 P 31 2 1
SCALE1      0.021524  0.012427  0.000000        0.00000
SCALE2      0.000000  0.024854  0.000000        0.00000
SCALE3      0.000000  0.000000  0.005176        0.00000
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
ATOM     17  N   SER A   2A     26.854 -17.177   0.412  1.00  6.66           N
ATOM     18  CA  SER A   2A     25.899 -16.083   0.383  0.53  6.83           C
ATOM     20  C   SER A   2A     26.500 -14.946  -0.440  1.00  5.61           C
ATOM     21  O   SER A   2A     27.729 -14.822  -0.565  1.00  8.14           O
ATOM     22  CB  SER A   2A     25.569 -15.634   1.795  0.53  7.38           C
ATOM     24  OG  SER A   2A     26.740 -15.136   2.390  0.53  9.79           O
ATOM     26  N   CYS A   2B     25.627 -14.135  -0.995  1.00  5.11           N
ATOM     27  CA  CYS A   2B     26.070 -13.062  -1.865  1.00  5.65           C
ATOM     28  C   CYS A   2B     25.043 -11.934  -1.798  1.00  4.54           C
ATOM     29  O   CYS A   2B     23.856 -12.180  -1.528  1.00  5.60           O
ATOM     30  CB  CYS A   2B     26.253 -13.557  -3.295  1.00  7.00           C
ATOM     31  SG  CYS A   2B     24.806 -14.269  -4.119  1.00  8.88           S
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
ATOM    234  N   SER A  35A     23.294 -20.178  -1.426  1.00  4.63           N
ATOM    235  CA  SER A  35A     23.402 -21.406  -2.221  1.00  4.58           C
ATOM    236  C   SER A  35A     24.387 -22.368  -1.553  1.00  5.05           C
ATOM    237  O   SER A  35A     24.929 -22.095  -0.497  1.00  6.22           O
ATOM    238  CB  SER A  35A     23.881 -21.071  -3.639  1.00  5.69           C
ATOM    239  OG  SER A  35A     25.213 -20.561  -3.633  1.00  7.12           O
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
ATOM     72  N   THR A  15A      4.929  -8.263  -7.636  1.00  5.49           N
ATOM     73  C   THR A  15A      5.316 -10.600  -7.084  1.00  5.07           C
ATOM     74  O   THR A  15A      4.214 -11.002  -7.422  1.00  6.51           O
ATOM     75  CA  THR A  15A      5.513  -9.172  -6.657  1.00  5.70           C
ATOM     76  CB  THR A  15A      4.864  -9.001  -5.276  1.00  8.31           C
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
ATOM     96  N   THR A  20A     10.560 -24.542  -5.569  1.00  4.00           N
ATOM     97  CA  THR A  20A     10.597 -25.962  -5.876  1.00  4.05           C
ATOM     98  C   THR A  20A     10.984 -26.770  -4.636  1.00  4.53           C
ATOM     99  O   THR A  20A     11.770 -26.361  -3.802  1.00  5.04           O
ATOM    100  CB  THR A  20A     11.488 -26.293  -7.083  1.00  4.38           C
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
ATOM    204  N   LYS A  42A      2.479 -16.398  -7.978  1.00  6.26           N
ATOM    205  CA  LYS A  42A      3.247 -17.395  -7.256  1.00  6.48           C
ATOM    206  C   LYS A  42A      3.186 -18.741  -7.955  1.00  5.78           C
ATOM    207  O   LYS A  42A      2.206 -19.041  -8.623  1.00  9.40           O
ATOM    208  CB  LYS A  42A      2.727 -17.535  -5.836  1.00  8.81           C
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
ATOM    219  N   ALA A  45A      6.347 -26.119  -7.403  1.00  3.62           N
ATOM    220  CA  ALA A  45A      6.443 -27.338  -8.202  1.00  3.99           C
ATOM    221  C   ALA A  45A      7.579 -28.205  -7.717  1.00  4.57           C
ATOM    222  O   ALA A  45A      8.479 -27.750  -7.000  1.00  4.79           O
ATOM    223  CB  ALA A  45A      6.607 -27.026  -9.678  1.00  4.52           C
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
ATOM     11  N  ATHR A  15A      4.929  -8.263  -7.636  0.50  5.49           N
ATOM     12  C  ATHR A  15A      5.316 -10.600  -7.084  0.50  5.07           C
ATOM     13  O  ATHR A  15A      4.214 -11.002  -7.422  0.50  6.51           O
ATOM     14  CA ATHR A  15A      5.513  -9.172  -6.657  0.50  5.70           C
ATOM     15  CB ATHR A  15A      4.864  -9.001  -5.276  0.50  8.31           C
ATOM     16  N  BTHR A  15A      5.029  -8.263  -7.636  0.50  5.49           N
ATOM     17  C  BTHR A  15A      5.416 -10.600  -7.084  0.50  5.07           C
ATOM     18  O  BTHR A  15A      4.314 -11.002  -7.422  0.50  6.51           O
ATOM     19  CA BTHR A  15A      5.613  -9.172  -6.657  0.50  5.70           C
ATOM     20  CB BTHR A  15A      4.964  -9.001  -5.276  0.50  8.31           C
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
ATOM     60  N   LYS A  42A      2.479 -16.398  -7.978  1.00  6.26           N
ATOM     61  CA  LYS A  42A      3.247 -17.395  -7.256  1.00  6.48           C
ATOM     62  C   LYS A  42A      3.186 -18.741  -7.955  1.00  5.78           C
ATOM     63  O   LYS A  42A      2.206 -19.041  -8.623  1.00  9.40           O
ATOM     64  CB  LYS A  42A      2.727 -17.535  -5.836  1.00  8.81           C
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

  s_apar_records1 = """\
SHEET    1   A 2 TYR A   2  VAL A   7  0
SHEET    2   A 2 MET A  31  SER A  35A-1  O  ALA A  35A  N  ALA A   2
"""
  s_apar_records2 = """\
SHEET    1   A 2 TYR A   2  VAL A   7  0
SHEET    2   A 2 MET A  31  SER A  35A-1  O  ALA A  32   N  ALA A   6
"""
  s_apar_records3 = """\
SHEET    1   A 2 SER A   2A VAL A   7  0
SHEET    2   A 2 MET A  31  SER A  35A-1  O  ALA A  32   N  ALA A   6
"""

  s_par_records1 = """\
SHEET    1   B 2 ALA A  15  THR A  20A 0
SHEET    2   B 2 GLN A  40  ALA A  45A 1  O  GLN A  40   N  THR A  15A
"""
  s_par_records2 = """\
SHEET    1   B 2 ALA A  15  THR A  20A 0
SHEET    2   B 2 GLN A  40  ALA A  45A 1  O  GLN A  44   N  THR A  20
"""
  s_par_records3 = """\
SHEET    1   B 2 ALA A  15A THR A  20A 0
SHEET    2   B 2 GLN A  40  ALA A  45  1  O  GLN A  44   N  THR A  20
"""

  log = null_out()
  # defpars = sec_str_master_phil
  n_hbonds = []
  n_hangles = []
  for pdb_h, recs in [
                        (pdb_apar_h, s_apar_records1), # 6
                        (pdb_apar_h, s_apar_records2), # 6
                        (pdb_apar_h, s_apar_records3), # 4
                        (pdb_par_h,  s_par_records1), # 6
                        (pdb_par_h,  s_par_records2), # 6
                        (pdb_par_h,  s_par_records3), # 5
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
    proxies_for_grm, hangles = ss_manager.create_protein_hbond_proxies(
      annotation= None,
      log          = log)
    # print proxies_for_grm.size()
    n_hbonds.append(proxies_for_grm.size())
    n_hangles.append(hangles.size())
  print(n_hbonds, n_hangles)
  assert n_hbonds  == [6, 6, 4, 6, 6, 5, 8]
  assert n_hangles == [15, 15, 12, 18, 18, 12, 24]


def exercise_phil_generation():
  log = null_out()
  defpars = sec_str_master_phil.fetch()

if __name__ == "__main__" :
  exercise_helix_bonding_pattern_with_insertions()
  exercise_sheets_bonding_pattern_with_insertions()
  print("OK")
