
from __future__ import division
from libtbx import easy_run
import os

def run (file_base="ca_frag") :
  mtz_file = file_base + ".mtz"
  if (os.path.isfile(mtz_file)) :
    os.remove(mtz_file)
  write_pdb_input(file_base=file_base)
  params = """
    high_resolution = 1.5
    r_free_flags_fraction = 0.1
    add_sigmas = True
    pdb_file = %s.pdb
    output {
      label = F
      type = *real complex
      file_name = %s.mtz
    }
    anomalous_scatterers {
      group {
        selection = element CA
        f_prime = 0.25
        f_double_prime = 0.5
      }
    }""" % (file_base, file_base)
  open("%s_fmodel.eff" % file_base, "w").write(params)
  assert (easy_run.fully_buffered("phenix.fmodel %s_fmodel.eff" % file_base,
    ).raise_if_errors().return_code == 0)
  return os.path.abspath(mtz_file), os.path.abspath("%s.pdb" % file_base)

def write_pdb_input (file_base="ca_frag") :
  import iotbx.pdb.hierarchy
  pdb_in = iotbx.pdb.hierarchy.input(source_info=None, pdb_string="""\
ATOM      1  N   ASP A  37      10.710  14.456   9.568  1.00 15.78           N
ATOM      2  CA  ASP A  37       9.318  14.587   9.999  1.00 18.38           C
ATOM      3  C   ASP A  37       8.402  13.523   9.395  1.00 15.46           C
ATOM      4  O   ASP A  37       7.295  13.309   9.892  1.00 16.65           O
ATOM      5  CB  ASP A  37       8.771  16.006   9.784  1.00 24.16           C
ATOM      6  CG  ASP A  37       8.900  16.478   8.357  1.00 28.30           C
ATOM      7  OD1 ASP A  37       9.145  15.643   7.457  1.00 30.98           O
ATOM      8  OD2 ASP A  37       8.747  17.700   8.132  1.00 36.64           O-1
ATOM      9  N   LEU A  38       8.862  12.840   8.346  1.00 16.14           N
ATOM     10  CA  LEU A  38       8.092  11.733   7.778  1.00 16.45           C
ATOM     11  C   LEU A  38       7.990  10.538   8.735  1.00 15.91           C
ATOM     12  O   LEU A  38       7.160   9.650   8.546  1.00 16.71           O
ATOM     13  CB  LEU A  38       8.675  11.305   6.425  1.00 17.74           C
ATOM     14  CG  LEU A  38       8.715  12.386   5.345  1.00 18.16           C
ATOM     15  CD1 LEU A  38       9.091  11.771   4.003  1.00 20.50           C
ATOM     16  CD2 LEU A  38       7.391  13.166   5.250  1.00 19.36           C
ATOM     17  N   ALA A  39       8.828  10.526   9.770  1.00 13.52           N
ATOM     18  CA  ALA A  39       8.825   9.454  10.756  1.00 15.05           C
ATOM     19  C   ALA A  39       7.758   9.662  11.829  1.00 16.22           C
ATOM     20  O   ALA A  39       7.489   8.754  12.609  1.00 15.36           O
ATOM     21  CB  ALA A  39      10.201   9.341  11.405  1.00 13.76           C
ATOM     22  N   ILE A  40       7.164  10.855  11.883  1.00 14.32           N
ATOM     23  CA  ILE A  40       6.164  11.156  12.915  1.00 15.53           C
ATOM     24  C   ILE A  40       4.885  11.819  12.386  1.00 18.72           C
ATOM     25  O   ILE A  40       4.178  12.497  13.136  1.00 21.36           O
ATOM     26  CB  ILE A  40       6.752  12.039  14.047  1.00 17.76           C
ATOM     27  CG1 ILE A  40       7.253  13.375  13.493  1.00 19.23           C
ATOM     28  CG2 ILE A  40       7.871  11.314  14.773  1.00 19.67           C
ATOM     29  CD1 ILE A  40       7.538  14.403  14.583  1.00 22.28           C
ATOM     30  N   ASP A  41       4.569  11.592  11.115  1.00 17.33           N
ATOM     31  CA  ASP A  41       3.461  12.305  10.474  1.00 20.42           C
ATOM     32  C   ASP A  41       2.132  11.539  10.466  1.00 23.63           C
ATOM     33  O   ASP A  41       1.117  12.047   9.976  1.00 24.40           O
ATOM     34  CB  ASP A  41       3.852  12.742   9.051  1.00 19.02           C
ATOM     35  CG  ASP A  41       4.135  11.566   8.114  1.00 20.05           C
ATOM     36  OD1 ASP A  41       4.047  10.386   8.529  1.00 17.79           O
ATOM     37  OD2 ASP A  41       4.465  11.830   6.938  1.00 19.57           O-1
ATOM     38  N   GLY A  42       2.132  10.328  11.012  1.00 19.52           N
ATOM     39  CA  GLY A  42       0.917   9.531  11.073  1.00 19.91           C
ATOM     40  C   GLY A  42       0.622   8.804   9.775  1.00 25.32           C
ATOM     41  O   GLY A  42      -0.390   8.104   9.659  1.00 25.73           O
ATOM     42  N   ASN A  43       1.505   8.967   8.794  1.00 20.57           N
ATOM     43  CA  ASN A  43       1.326   8.337   7.488  1.00 20.84           C
ATOM     44  C   ASN A  43       2.269   7.149   7.310  1.00 20.54           C
ATOM     45  O   ASN A  43       3.484   7.322   7.215  1.00 18.86           O
ATOM     46  CB  ASN A  43       1.563   9.367   6.382  1.00 20.63           C
ATOM     47  CG  ASN A  43       1.087   8.896   5.020  1.00 19.17           C
ATOM     48  ND2 ASN A  43       1.027   9.825   4.069  1.00 24.27           N
ATOM     49  OD1 ASN A  43       0.793   7.718   4.813  1.00 21.70           O
ATOM     50  N   PRO A  44       1.711   5.933   7.258  1.00 20.33           N
ATOM     51  CA  PRO A  44       2.542   4.737   7.081  1.00 20.08           C
ATOM     52  C   PRO A  44       3.221   4.664   5.716  1.00 20.25           C
ATOM     53  O   PRO A  44       4.166   3.891   5.562  1.00 26.14           O
ATOM     54  CB  PRO A  44       1.538   3.589   7.222  1.00 20.20           C
ATOM     55  CG  PRO A  44       0.211   4.201   6.901  1.00 25.57           C
ATOM     56  CD  PRO A  44       0.288   5.596   7.442  1.00 19.95           C
ATOM     57  N   ALA A  45       2.750   5.441   4.743  1.00 20.05           N
ATOM     58  CA  ALA A  45       3.316   5.406   3.393  1.00 23.75           C
ATOM     59  C   ALA A  45       4.525   6.323   3.231  1.00 25.35           C
ATOM     60  O   ALA A  45       5.211   6.281   2.205  1.00 26.63           O
ATOM     61  CB  ALA A  45       2.246   5.741   2.353  1.00 24.89           C
ATOM     62  N   THR A  46       4.783   7.153   4.239  1.00 19.37           N
ATOM     63  CA  THR A  46       5.952   8.027   4.226  1.00 19.58           C
ATOM     64  C   THR A  46       6.923   7.556   5.292  1.00 19.22           C
ATOM     65  O   THR A  46       6.516   6.968   6.298  1.00 16.37           O
ATOM     66  CB  THR A  46       5.588   9.477   4.547  1.00 17.37           C
ATOM     67  CG2 THR A  46       4.599  10.035   3.514  1.00 21.06           C
ATOM     68  OG1 THR A  46       4.999   9.537   5.849  1.00 18.86           O
ATOM     69  N   SER A  47       8.206   7.816   5.068  1.00 17.65           N
ATOM     70  CA  SER A  47       9.224   7.383   6.009  1.00 16.72           C
ATOM     71  C   SER A  47      10.503   8.190   5.885  1.00 17.50           C
ATOM     72  O   SER A  47      10.765   8.834   4.859  1.00 15.22           O
ATOM     73  CB  SER A  47       9.541   5.916   5.774  1.00 17.26           C
ATOM     74  OG  SER A  47      10.208   5.738   4.530  1.00 18.61           O
ATOM     75  N   TRP A  48      11.299   8.139   6.947  1.00 14.32           N
ATOM     76  CA  TRP A  48      12.634   8.723   6.955  1.00 13.15           C
ATOM     77  C   TRP A  48      13.652   7.627   6.683  1.00 13.64           C
ATOM     78  O   TRP A  48      13.756   6.646   7.433  1.00 12.78           O
ATOM     79  CB  TRP A  48      12.912   9.395   8.300  1.00 13.01           C
ATOM     80  CG  TRP A  48      14.327   9.861   8.482  1.00 13.43           C
ATOM     81  CD1 TRP A  48      14.924  10.924   7.876  1.00 14.36           C
ATOM     82  CD2 TRP A  48      15.313   9.291   9.357  1.00 13.44           C
ATOM     83  CE2 TRP A  48      16.485  10.060   9.221  1.00 12.26           C
ATOM     84  CE3 TRP A  48      15.315   8.202  10.239  1.00 12.49           C
ATOM     85  NE1 TRP A  48      16.228  11.048   8.307  1.00 13.42           N
ATOM     86  CZ2 TRP A  48      17.653   9.778   9.940  1.00 13.96           C
ATOM     87  CZ3 TRP A  48      16.476   7.922  10.948  1.00 11.99           C
ATOM     88  CH2 TRP A  48      17.628   8.705  10.790  1.00 13.64           C
TER      89      TRP A  48
ATOM     90  N   ILE A 153      12.422   4.326   8.671  1.00 13.05           N
ATOM     91  CA  ILE A 153      11.440   4.432   9.754  1.00 14.48           C
ATOM     92  C   ILE A 153      10.194   5.191   9.295  1.00 12.58           C
ATOM     93  O   ILE A 153      10.264   6.388   8.994  1.00 13.74           O
ATOM     94  CB  ILE A 153      12.019   5.189  10.975  1.00 15.38           C
ATOM     95  CG1 ILE A 153      13.320   4.550  11.466  1.00 15.63           C
ATOM     96  CG2 ILE A 153      10.986   5.280  12.090  1.00 15.20           C
ATOM     97  CD1 ILE A 153      13.950   5.283  12.647  1.00 15.77           C
ATOM     98  N   SER A 154       9.051   4.511   9.278  1.00 14.97           N
ATOM     99  CA  SER A 154       7.820   5.135   8.801  1.00 14.50           C
ATOM    100  C   SER A 154       7.001   5.838   9.874  1.00 15.23           C
ATOM    101  O   SER A 154       6.353   6.843   9.589  1.00 15.50           O
ATOM    102  CB  SER A 154       6.944   4.134   8.052  1.00 16.67           C
ATOM    103  OG  SER A 154       6.522   3.100   8.920  1.00 16.82           O
ATOM    104  N   GLU A 155       7.018   5.314  11.097  1.00 15.34           N
ATOM    105  CA  GLU A 155       6.218   5.885  12.176  1.00 15.01           C
ATOM    106  C   GLU A 155       6.825   5.575  13.533  1.00 12.95           C
ATOM    107  O   GLU A 155       7.234   4.443  13.789  1.00 15.35           O
ATOM    108  CB  GLU A 155       4.780   5.346  12.139  1.00 17.13           C
ATOM    109  CG  GLU A 155       3.862   6.037  11.124  1.00 22.41           C
ATOM    110  CD  GLU A 155       3.866   7.552  11.272  1.00 24.24           C
ATOM    111  OE1 GLU A 155       3.618   8.050  12.389  1.00 26.14           O
ATOM    112  OE2 GLU A 155       4.132   8.250  10.271  1.00 23.46           O-1
ATOM    113  N   ILE A 156       6.882   6.596  14.382  1.00 13.75           N
ATOM    114  CA  ILE A 156       7.220   6.440  15.788  1.00 12.29           C
ATOM    115  C   ILE A 156       6.066   7.007  16.602  1.00 11.05           C
ATOM    116  O   ILE A 156       5.584   8.113  16.320  1.00 13.14           O
ATOM    117  CB  ILE A 156       8.506   7.196  16.147  1.00 11.71           C
ATOM    118  CG1 ILE A 156       9.714   6.573  15.433  1.00 12.75           C
ATOM    119  CG2 ILE A 156       8.742   7.158  17.657  1.00 12.85           C
ATOM    120  CD1 ILE A 156      10.932   7.490  15.398  1.00 13.97           C
TER     121      ILE A 156
HETATM  122 CA   CA  S   1       5.334   8.357   8.032  1.00  5.89          CA+2
HETATM  123  O   HOH S   2       5.396  15.243  10.734  1.00 22.95           O
HETATM  124  O   HOH S   3       3.629   8.994  14.414  1.00 25.28           O
END""")
  xrs = pdb_in.input.xray_structure_simple(cryst1_substitution_buffer_layer=5)
  pdb_file = file_base + ".pdb"
  if (os.path.exists(pdb_file)) :
    os.remove(pdb_file)
  f = open(pdb_file, "w")
  f.write(pdb_in.hierarchy.as_pdb_string(crystal_symmetry=xrs))
  f.close()
