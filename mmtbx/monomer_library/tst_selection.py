from mmtbx import monomer_library
import mmtbx.monomer_library.pdb_interpretation
import mmtbx.monomer_library.server
import sys

test_pdb_file = """\
CRYST1  201.710   68.060   86.350  90.00  90.00  90.00 P 21 21 2     8
ATOM      1  N   SER A   1      79.022  -1.487  39.367  1.00 97.85           N
ATOM      2  CA  SER A   1      78.722  -1.348  37.919  1.00 98.11           C
ATOM      3  C   SER A   1      77.318  -0.760  37.705  1.00 97.34           C
ATOM      4  O   SER A   1      76.305  -1.440  37.908  1.00 98.14           O
ATOM      5  CB  SER A   1      78.826  -2.717  37.217  1.00 98.90           C
ATOM      6  OG  SER A   1      80.160  -3.231  37.252  1.00 99.99           O
ATOM      7  N   GLU A   2      77.245   0.508  37.315  1.00 95.57           N
ATOM      8  CA  GLU A   2      75.961   1.137  37.056  1.00 91.95           C
ATOM      9  C   GLU A   2      76.189   1.866  35.740  1.00 88.32           C
ATOM     10  O   GLU A   2      77.119   2.688  35.625  1.00 86.70           O
ATOM     11  CB  GLU A   2      75.615   2.131  38.169  1.00 92.62           C
ATOM     12  CG  GLU A   2      76.373   3.428  38.057  1.00 96.62           C
ATOM     13  CD  GLU A   2      76.044   4.392  39.120  1.00 99.90           C
ATOM     14  OE1 GLU A   2      74.977   4.253  39.777  1.00 99.99           O
ATOM     15  OE2 GLU A   2      76.841   5.347  39.313  1.00 99.99           O
ATOM     16  N   SER A   3      75.396   1.547  34.723  1.00 83.22           N
ATOM     17  CA  SER A   3      75.583   2.238  33.461  1.00 78.14           C
ATOM     18  C   SER A   3      75.232   3.680  33.734  1.00 73.03           C
ATOM     19  O   SER A   3      74.305   3.974  34.499  1.00 72.80           O
ATOM     20  CB  SER A   3      74.678   1.646  32.370  1.00 80.23           C
ATOM     21  OG  SER A   3      75.170   0.385  31.930  1.00 80.78           O
ATOM   2377  N   LEU A 313      82.669  33.333 -17.857  1.00 53.08           N
ATOM   2378  CA  LEU A 313      81.629  32.854 -16.955  1.00 53.14           C
ATOM   2379  C   LEU A 313      81.630  33.876 -15.842  1.00 54.15           C
ATOM   2380  O   LEU A 313      82.643  34.503 -15.597  1.00 55.52           O
ATOM   2381  CB  LEU A 313      82.006  31.491 -16.387  1.00 54.25           C
ATOM   2382  CG  LEU A 313      82.274  30.386 -17.411  1.00 56.18           C
ATOM   2383  CD1 LEU A 313      82.747  29.098 -16.751  1.00 52.25           C
ATOM   2384  CD2 LEU A 313      81.002  30.149 -18.180  1.00 59.92           C
ATOM   2385  N   VAL A 314      80.512  34.069 -15.166  1.00 56.43           N
ATOM   2386  CA  VAL A 314      80.488  35.043 -14.074  1.00 56.85           C
ATOM   2387  C   VAL A 314      80.903  34.401 -12.731  1.00 60.38           C
ATOM   2388  O   VAL A 314      81.023  33.150 -12.658  1.00 61.26           O
ATOM   2389  CB  VAL A 314      79.087  35.666 -13.946  1.00 52.88           C
ATOM   2390  CG1 VAL A 314      79.074  36.727 -12.868  1.00 47.39           C
ATOM   2391  CG2 VAL A 314      78.677  36.241 -15.293  1.00 52.66           C
ATOM   2392  OXT VAL A 314      81.110  35.156 -11.751  1.00 62.47           O
TER    2393      VAL A 314
ATOM   2394  N   SER B   1      85.404 -33.823  39.687  1.00 91.56           N
ATOM   2395  CA  SER B   1      84.670 -32.782  38.908  1.00 91.48           C
ATOM   2396  C   SER B   1      83.295 -33.255  38.414  1.00 92.36           C
ATOM   2397  O   SER B   1      83.109 -34.423  38.036  1.00 92.95           O
ATOM   2398  CB  SER B   1      85.529 -32.298  37.718  1.00 89.71           C
ATOM   2399  OG  SER B   1      85.912 -33.355  36.850  1.00 80.92           O
ATOM   2400  N   GLU B   2      82.331 -32.338  38.418  1.00 91.53           N
ATOM   2401  CA  GLU B   2      80.983 -32.666  37.980  1.00 89.92           C
ATOM   2402  C   GLU B   2      80.567 -31.880  36.725  1.00 86.18           C
ATOM   2403  O   GLU B   2      81.233 -30.933  36.292  1.00 84.19           O
ATOM   2404  CB  GLU B   2      79.993 -32.434  39.136  1.00 92.51           C
ATOM   2405  CG  GLU B   2      78.594 -33.015  38.906  1.00 97.08           C
ATOM   2406  CD  GLU B   2      77.710 -32.959  40.153  1.00 99.61           C
ATOM   2407  OE1 GLU B   2      77.422 -31.843  40.638  1.00 99.99           O
ATOM   2408  OE2 GLU B   2      77.305 -34.037  40.650  1.00 99.99           O
TER    4786      VAL B 314
ATOM   4788  P     C C   1      47.241  23.887  -1.304  1.00 99.99           P
ATOM   4789  O1P   C C   1      47.437  24.960  -0.241  1.00 99.99           O
ATOM   4790  O2P   C C   1      46.407  22.725  -0.791  1.00 99.99           O
ATOM   4791  O5'   C C   1      48.722  23.275  -1.633  1.00 99.99           O
ATOM   4792  C5'   C C   1      49.588  22.811  -0.567  1.00 99.99           C
ATOM   4793  C4'   C C   1      49.474  21.311  -0.424  1.00 99.98           C
ATOM   4794  O4'   C C   1      48.276  20.855  -1.121  1.00 99.99           O
ATOM   4795  C3'   C C   1      50.570  20.506  -1.099  1.00 99.99           C
ATOM   4796  O3'   C C   1      51.727  20.252  -0.373  1.00 99.99           O
ATOM   4797  C2'   C C   1      49.934  19.165  -1.338  1.00 99.99           C
ATOM   4798  O2'   C C   1      50.017  18.352  -0.186  1.00 99.99           O
ATOM   4799  C1'   C C   1      48.505  19.564  -1.689  1.00 99.99           C
ATOM   4800  N1    C C   1      48.330  19.637  -3.157  1.00 99.99           N
ATOM   4801  C2    C C   1      48.881  18.604  -3.984  1.00 99.99           C
ATOM   4802  O2    C C   1      49.493  17.647  -3.456  1.00 99.99           O
ATOM   4803  N3    C C   1      48.728  18.684  -5.330  1.00 99.99           N
ATOM   4804  C4    C C   1      48.065  19.714  -5.869  1.00 99.99           C
ATOM   4805  N4    C C   1      47.945  19.750  -7.196  1.00 99.99           N
ATOM   4806  C5    C C   1      47.497  20.758  -5.068  1.00 99.99           C
ATOM   4807  C6    C C   1      47.655  20.682  -3.732  1.00 99.99           C
ATOM   4808  P     G C   2      52.970  19.524  -1.111  1.00 99.99           P
ATOM   4809  O1P   G C   2      52.775  18.136  -1.648  1.00 99.99           O
ATOM   4810  O2P   G C   2      53.883  19.724  -0.024  1.00 99.99           O
ATOM   4811  O5*   G C   2      53.329  20.523  -2.299  1.00 97.90           O
ATOM   4812  C5*   G C   2      52.379  21.525  -2.667  1.00 92.61           C
ATOM   4813  C4*   G C   2      53.032  22.857  -2.853  1.00 88.12           C
ATOM   4814  O4*   G C   2      54.015  23.115  -1.816  1.00 87.30           O
ATOM   4815  C3*   G C   2      53.805  22.945  -4.141  1.00 86.97           C
ATOM   4816  O3*   G C   2      52.918  23.175  -5.222  1.00 86.86           O
ATOM   4817  C2*   G C   2      54.790  24.065  -3.839  1.00 87.48           C
ATOM   4818  O2*   G C   2      54.217  25.351  -3.939  1.00 86.43           O
ATOM   4819  C1*   G C   2      55.146  23.758  -2.382  1.00 86.25           C
ATOM   4820  N9    G C   2      56.247  22.803  -2.383  1.00 86.73           N
ATOM   4821  C8    G C   2      56.229  21.497  -1.947  1.00 87.03           C
ATOM   4822  N7    G C   2      57.342  20.861  -2.190  1.00 85.88           N
ATOM   4823  C5    G C   2      58.151  21.807  -2.801  1.00 85.19           C
ATOM   4824  C6    G C   2      59.462  21.698  -3.313  1.00 84.09           C
ATOM   4825  O6    G C   2      60.198  20.699  -3.353  1.00 82.37           O
ATOM   4826  N1    G C   2      59.905  22.912  -3.834  1.00 84.56           N
ATOM   4827  C2    G C   2      59.174  24.077  -3.864  1.00 84.51           C
ATOM   4828  N2    G C   2      59.774  25.156  -4.399  1.00 84.63           N
ATOM   4829  N3    G C   2      57.946  24.181  -3.406  1.00 85.45           N
ATOM   4830  C4    G C   2      57.499  23.018  -2.894  1.00 85.82           C
ATOM   4831  P     C C   3      53.092  22.329  -6.579  1.00 86.68           P
ATOM   4832  O1P   C C   3      51.825  22.390  -7.355  1.00 86.55           O
ATOM   4833  O2P   C C   3      53.667  21.001  -6.222  1.00 84.84           O
ATOM   4834  O5*   C C   3      54.167  23.202  -7.360  1.00 84.63           O
ATOM   4835  C5*   C C   3      53.965  24.621  -7.484  1.00 81.25           C
ATOM   4836  C4*   C C   3      55.164  25.277  -8.115  1.00 79.23           C
ATOM   4837  O4*   C C   3      56.262  25.305  -7.172  1.00 76.28           O
ATOM   4838  C3*   C C   3      55.737  24.591  -9.341  1.00 77.29           C
ATOM   4839  O3*   C C   3      55.051  24.887 -10.537  1.00 78.19           O
ATOM   4840  C2*   C C   3      57.158  25.115  -9.363  1.00 76.40           C
ATOM   4842  C1*   C C   3      57.483  25.148  -7.872  1.00 76.01           C
ATOM   4843  N1    C C   3      58.067  23.861  -7.480  1.00 73.97           N
ATOM   4844  C2    C C   3      59.445  23.727  -7.509  1.00 71.90           C
ATOM   4845  O2    C C   3      60.135  24.715  -7.840  1.00 66.31           O
ATOM   4846  N3    C C   3      59.996  22.533  -7.175  1.00 70.72           N
ATOM   4847  C4    C C   3      59.212  21.513  -6.810  1.00 69.96           C
ATOM   4848  N4    C C   3      59.793  20.365  -6.470  1.00 71.43           N
ATOM   4849  C5    C C   3      57.800  21.629  -6.772  1.00 69.13           C
ATOM   4850  C6    C C   3      57.272  22.810  -7.109  1.00 72.22           C
HETATM    1  N   PCA C   1      61.815  -9.826  47.209  1.00 22.41           N
HETATM    2  CA  PCA C   1      63.251  -9.543  47.447  1.00 21.95           C
HETATM    3  CB  PCA C   1      64.110 -10.479  46.604  1.00 22.40           C
HETATM    4  CG  PCA C   1      63.238 -10.886  45.438  1.00 22.83           C
HETATM    5  CD  PCA C   1      61.855 -10.549  45.934  1.00 22.46           C
HETATM    6  OE  PCA C   1      60.846 -10.859  45.286  1.00 24.62           O
HETATM    7  C   PCA C   1      63.573  -8.105  47.060  1.00 21.16           C
HETATM    8  O   PCA C   1      62.762  -7.428  46.430  1.00 21.46           O
HETATM 8099 MG    MG   584      69.634   3.963 -22.036  1.00 45.86          MG
HETATM 8100 MG    MG   585      72.985 -30.470 -21.424  1.00 39.36          MG
HETATM 8101  O   HOH   501      50.104 -12.783 -34.164  1.00 76.32           O
HETATM 8102  O   HOH   502      50.327  -5.724 -17.983  1.00 53.54           O
END
"""

def exercise_atom_selections():
  verbose = "--verbose" in sys.argv[1:]
  log = None
  if (verbose):
    log = sys.stdout
  open("tmp_selection.pdb", "w").write(test_pdb_file)
  try:
    mon_lib_srv = monomer_library.server.server()
  except monomer_library.server.MonomerLibraryServerError:
    print "Skipping exercise_atom_selections(): monomer library not available."
    return
  ener_lib = monomer_library.server.ener_lib()
  processed_pdb = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name="tmp_selection.pdb",
    log=log)
  isel = processed_pdb.all_chain_proxies.iselection
  assert list(isel("peptide")) == range(52) + range(114,122)
  assert list(isel("rna")) == range(52,95)
  assert list(isel("dna")) == range(95,114)
  assert list(isel("nucleotide")) == range(52,114)
  assert list(isel("water")) == [124,125]
  assert list(isel("peptide backbone")) \
      == list(isel("peptide and (name ca or name c or name o or name n)"))
  assert list(isel("nucleotide backbone")) \
      == [52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63,
          72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83,
          95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105]
  assert list(isel("backbone and (rna or dna)")) \
      == list(isel("nucleotide backbone"))
  assert list(isel("phosphate")) \
      == [52, 53, 54, 55, 60, 72, 73, 74, 75, 80, 95, 96, 97, 98, 103]
  assert list(isel("ribose")) \
      == [56, 57, 58, 59, 61, 62, 63,
          76, 77, 78, 79, 81, 82, 83,
          99, 100, 101, 102, 104, 105]
  assert list(isel("phosphate or ribose")) == list(isel("nuc backbone"))
  assert list(isel("backbone resname glu")) \
      == [6, 7, 8, 9, 43, 44, 45, 46]
  assert list(isel("within(0.5, backbone resname glu)")) \
      == [6, 7, 8, 9, 43, 44, 45, 46]
  assert list(isel("within(1.5, backbone resname glu)")) \
      == [2, 6, 7, 8, 9, 15, 39, 43, 44, 45, 46]
  assert list(isel("within(2.5, backbone resname glu)")) \
      == [1, 2, 3, 6, 7, 8, 9, 10, 15, 16, 38, 39, 40, 43, 44, 45, 46, 47]
  assert list(isel("resname ser sidechain")) == [4, 5, 19, 20, 41, 42]

def exercise():
  exercise_atom_selections()
  print "OK"

if (__name__ == "__main__"):
  exercise()
