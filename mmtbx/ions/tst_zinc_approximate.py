
from __future__ import division
from libtbx import easy_run
import os

def exercise () :
  from iotbx.file_reader import any_file
  if (os.path.isfile("zn_frag.mtz")) :
    os.remove("zn_frag.mtz")
  write_pdb_input()
  params = """
    high_resolution = 1.9
    r_free_flags_fraction = 0.1
    add_sigmas = True
    pdb_file = zn_frag.pdb
    output {
      label = F
      type = *real complex
      file_name = zn_frag.mtz
    }
    anomalous_scatterers {
      group {
        selection = element ZN
        f_prime = -1.3
        f_double_prime = 0.47
      }
    }"""
  open("zn_frag_fmodel.eff", "w").write(params)
  assert (easy_run.fully_buffered("phenix.fmodel zn_frag_fmodel.eff"
    ).raise_if_errors().return_code == 0)
  pdb_in = any_file("zn_frag.pdb")
  hierarchy = pdb_in.file_object.construct_hierarchy()
  for chain in hierarchy.models()[0].chains() :
    for residue_group in chain.residue_groups() :
      for atom_group in residue_group.atom_groups() :
        if (atom_group.resname == "ZN ") :
          atom_group.resname = "HOH"
          atom = atom_group.atoms()[0]
          atom.name = " O  "
          atom.element = " O"
          atom.charge = ""
          atom.segid = "ZN"
          break
  f = open("zn_frag_hoh.pdb", "w")
  f.write(hierarchy.as_pdb_string(pdb_in.file_object.crystal_symmetry()))
  f.close()
  args = ["zn_frag_hoh.pdb", "zn_frag.mtz", "wavelength=1.54", "nproc=1",
          "elements=CA,ZN"]
  result = easy_run.fully_buffered("mmtbx.water_screen %s" % " ".join(args)
    ).raise_if_errors()
  if (not "  Probable cation: ZN+2" in result.stdout_lines) :
    print "\n".join(result.stdout_lines)
    raise RuntimeError("Missing ZN")
  print "OK"

def write_pdb_input () :
  import iotbx.pdb.hierarchy
  # XXX extracted from 2whz (thermolysin, wavelength = 1.54A), but with an
  # incomplete coordination shell
  pdb_in = iotbx.pdb.hierarchy.input(source_info=None, pdb_string="""\
ATOM      1  OD1 ASN A 112      34.902  38.381  -2.709  1.00 12.10           O
ATOM      2  C   ALA A 113      31.533  40.581  -3.237  1.00 12.06           C
ATOM      3  O   ALA A 113      32.560  40.184  -3.842  1.00 12.75           O
ATOM      4  N   PHE A 114      30.551  41.251  -3.854  1.00 11.88           N
ATOM      5  CA  PHE A 114      30.707  41.621  -5.277  1.00 12.52           C
ATOM      6  C   PHE A 114      29.633  42.613  -5.705  1.00 12.63           C
ATOM      7  CB  PHE A 114      30.665  40.371  -6.202  1.00 12.72           C
ATOM      8  N   TRP A 115      29.933  43.311  -6.794  1.00 13.30           N
ATOM      9  CA  TRP A 115      28.968  44.084  -7.575  1.00 13.94           C
ATOM     10  O   TRP A 115      29.375  42.549  -9.378  1.00 15.06           O
ATOM     11  CB  TRP A 115      29.646  45.329  -8.137  1.00 14.57           C
ATOM     13  O   VAL A 139      33.773  47.006  -1.063  1.00 11.73           O
ATOM     15  N   HIS A 142      34.994  49.588  -2.455  1.00 13.21           N
ATOM     16  CA  HIS A 142      35.593  48.855  -3.582  1.00 11.38           C
ATOM     17  C   HIS A 142      34.455  48.358  -4.484  1.00 12.11           C
ATOM     18  O   HIS A 142      34.468  48.596  -5.693  1.00 11.57           O
ATOM     19  CB  HIS A 142      36.442  47.680  -3.082  1.00 11.93           C
ATOM     20  CG  HIS A 142      36.945  46.756  -4.150  1.00  9.54           C
ATOM     21  CD2 HIS A 142      36.303  45.839  -4.912  1.00  9.50           C
ATOM     22  ND1 HIS A 142      38.289  46.628  -4.440  1.00 11.23           N
ATOM     23  CE1 HIS A 142      38.446  45.706  -5.381  1.00 14.08           C
ATOM     24  NE2 HIS A 142      37.258  45.206  -5.676  1.00 14.47           N
ATOM     25  N   GLU A 143      33.448  47.732  -3.879  1.00 11.17           N
ATOM     26  CA  GLU A 143      32.363  47.129  -4.683  1.00 12.04           C
ATOM     27  C   GLU A 143      31.546  48.215  -5.395  1.00 12.79           C
ATOM     28  O   GLU A 143      31.237  48.078  -6.584  1.00 12.97           O
ATOM     29  CB  GLU A 143      31.461  46.221  -3.834  1.00 11.99           C
ATOM     30  CG  GLU A 143      32.199  45.096  -3.030  1.00 12.94           C
ATOM     31  CD  GLU A 143      32.974  44.075  -3.867  1.00 15.46           C
ATOM     32  OE1 GLU A 143      33.029  44.153  -5.112  1.00 16.48           O
ATOM     33  OE2 GLU A 143      33.563  43.175  -3.251  1.00 16.80           O1-
ATOM     35  C   THR A 145      33.305  50.876  -8.608  1.00 12.68           C
ATOM     36  N   HIS A 146      33.071  49.560  -8.540  1.00 12.36           N
ATOM     37  CA  HIS A 146      32.869  48.794  -9.756  1.00 13.33           C
ATOM     38  CB  HIS A 146      32.643  47.318  -9.470  1.00 12.80           C
ATOM     39  CG  HIS A 146      33.919  46.567  -9.214  1.00 14.00           C
ATOM     40  CD2 HIS A 146      34.290  45.747  -8.206  1.00 12.81           C
ATOM     41  ND1 HIS A 146      34.992  46.623 -10.084  1.00 12.00           N
ATOM     42  CE1 HIS A 146      35.970  45.863  -9.617  1.00 15.40           C
ATOM     43  NE2 HIS A 146      35.552  45.289  -8.506  1.00 11.49           N
ATOM     45  CG  TYR A 157      36.255  42.005 -14.146  1.00 19.06           C
ATOM     46  CD1 TYR A 157      36.607  43.119 -13.385  1.00 20.01           C
ATOM     47  CD2 TYR A 157      36.279  40.754 -13.542  1.00 19.69           C
ATOM     48  CE1 TYR A 157      36.970  42.993 -12.047  1.00 18.89           C
ATOM     49  CE2 TYR A 157      36.645  40.612 -12.193  1.00 20.96           C
ATOM     50  CZ  TYR A 157      36.971  41.732 -11.462  1.00 20.87           C
ATOM     51  OH  TYR A 157      37.304  41.604 -10.142  1.00 21.10           O
ATOM     53  O   GLY A 162      40.788  44.717 -13.313  1.00 13.40           O
ATOM     55  CA  ASN A 165      39.280  47.961 -13.045  1.00 11.01           C
ATOM     56  C   ASN A 165      39.422  47.647 -11.555  1.00 11.21           C
ATOM     57  O   ASN A 165      38.940  48.402 -10.709  1.00 10.60           O
ATOM     58  CB  ASN A 165      38.414  46.889 -13.730  1.00 12.11           C
ATOM     59  CG  ASN A 165      36.977  46.958 -13.330  1.00 13.51           C
ATOM     60  ND2 ASN A 165      36.082  46.882 -14.307  1.00 12.19           N
ATOM     61  OD1 ASN A 165      36.661  47.056 -12.155  1.00 13.28           O
ATOM     62  N   GLU A 166      40.119  46.539 -11.250  1.00 11.43           N
ATOM     63  CA  GLU A 166      40.413  46.155  -9.851  1.00 10.59           C
ATOM     64  C   GLU A 166      41.166  47.218  -9.055  1.00 10.24           C
ATOM     65  O   GLU A 166      40.744  47.581  -7.940  1.00  9.90           O
ATOM     66  CB  GLU A 166      41.140  44.803  -9.783  1.00 10.41           C
ATOM     67  CG  GLU A 166      40.169  43.654  -9.741  1.00 11.01           C
ATOM     68  CD  GLU A 166      39.267  43.743  -8.537  1.00 15.78           C
ATOM     69  OE1 GLU A 166      39.745  43.921  -7.380  1.00 14.70           O
ATOM     70  OE2 GLU A 166      38.063  43.637  -8.765  1.00 12.05           O1-
ATOM     71  N   ALA A 167      42.243  47.759  -9.628  1.00 11.39           N
ATOM     72  C   ILE A 168      39.558  51.420  -8.190  1.00 10.98           C
ATOM     73  N   SER A 169      39.076  50.192  -8.040  1.00 10.62           N
ATOM     74  CA  SER A 169      38.198  49.860  -6.904  1.00 11.23           C
ATOM     75  C   SER A 169      38.998  49.850  -5.589  1.00 11.38           C
ATOM     76  O   SER A 169      38.484  50.325  -4.547  1.00 11.01           O
ATOM     77  CB  SER A 169      37.462  48.529  -7.144  1.00 12.04           C
ATOM     78  OG  SER A 169      36.259  48.790  -7.899  1.00 11.91           O
ATOM     79  N   ASP A 170      40.233  49.353  -5.645  1.00 11.74           N
ATOM     80  CA  ASP A 170      41.142  49.411  -4.436  1.00 11.68           C
ATOM     81  CB  ASP A 170      42.424  48.581  -4.612  1.00 11.14           C
ATOM     82  CG  ASP A 170      42.166  47.082  -4.475  1.00 10.59           C
ATOM     83  OD1 ASP A 170      41.076  46.659  -3.990  1.00 11.83           O
ATOM     84  OD2 ASP A 170      43.010  46.278  -4.909  1.00 12.10           O1-
ATOM     86  NE  ARG A 203      42.877  43.716  -3.830  1.00 12.26           N
ATOM     87  CZ  ARG A 203      41.692  43.111  -3.748  1.00 11.93           C
ATOM     88  NH1 ARG A 203      41.631  41.782  -3.754  1.00 12.22           N1+
ATOM     89  NH2 ARG A 203      40.579  43.829  -3.675  1.00 11.99           N
ATOM     91  C   VAL A 230      43.021  41.431 -10.925  1.00 14.80           C
ATOM     92  O   VAL A 230      43.311  42.599 -10.716  1.00 14.05           O
ATOM     93  CB  VAL A 230      41.294  41.138 -12.775  1.00 14.48           C
ATOM     94  CG2 VAL A 230      40.404  40.176 -11.957  1.00 14.82           C
ATOM     95  N   HIS A 231      42.923  40.503  -9.959  1.00 15.61           N
ATOM     96  CA  HIS A 231      43.233  40.790  -8.542  1.00 15.70           C
ATOM     97  CB  HIS A 231      42.434  39.898  -7.566  1.00 15.19           C
ATOM     98  CG  HIS A 231      40.955  40.000  -7.734  1.00 16.56           C
ATOM     99  CD2 HIS A 231      40.056  40.844  -7.174  1.00 17.34           C
ATOM    100  ND1 HIS A 231      40.241  39.168  -8.573  1.00 16.47           N
ATOM    101  CE1 HIS A 231      38.960  39.505  -8.529  1.00 15.36           C
ATOM    102  NE2 HIS A 231      38.825  40.518  -7.687  1.00 15.29           N
HETATM  104  N   ILE A1001      35.271  41.187  -4.048  1.00 16.78           N
HETATM  105  CA  ILE A1001      36.473  41.816  -3.412  1.00 17.77           C
HETATM  106  C   ILE A1001      37.764  41.046  -3.756  1.00 18.36           C
HETATM  107  O   ILE A1001      38.832  41.640  -3.838  1.00 16.99           O
HETATM  108  CB  ILE A1001      36.323  42.040  -1.875  1.00 17.98           C
HETATM  109  CG1 ILE A1001      37.464  42.899  -1.320  1.00 16.76           C
HETATM  110  CG2 ILE A1001      36.360  40.753  -1.098  1.00 20.77           C
HETATM  111  CD1 ILE A1001      37.576  44.234  -1.913  1.00 17.79           C
HETATM  112  N   TYR A1002      37.653  39.739  -3.981  1.00 18.85           N
HETATM  113  CA  TYR A1002      38.834  38.984  -4.422  1.00 20.77           C
HETATM  114  C   TYR A1002      38.430  37.858  -5.360  1.00 21.31           C
HETATM  115  O   TYR A1002      39.302  37.219  -5.983  1.00 22.07           O
HETATM  116  CB  TYR A1002      39.632  38.473  -3.225  1.00 21.20           C
HETATM  117  OXT TYR A1002      37.226  37.594  -5.543  1.00 22.10           O
HETATM  118  O   HOH A2080      32.747  43.340 -10.935  1.00 44.97           O
HETATM  119  O   HOH A2097      33.150  44.796 -12.822  1.00 35.29           O
HETATM  120  O   HOH A2190      32.628  42.780  -7.708  1.00 17.31           O
HETATM  121  O   HOH A2206      36.319  49.625 -10.585  1.00 16.17           O
HETATM  122  O   HOH A2210      33.342  47.349 -12.959  1.00 24.18           O
HETATM  123  O   HOH A2284      42.614  44.224  -6.731  1.00 15.85           O
HETATM  125 ZN   ZN  A1317      36.762  44.115  -7.277  1.00 15.32          ZN+2
END
""")
# XXX this is the missing water
#HETATM  124  O   HOH A2351      36.007  42.038  -6.782  1.00 10.54           O
  xrs = pdb_in.input.xray_structure_simple(cryst1_substitution_buffer_layer=5)
  if (os.path.exists("zn_frag.pdb")) :
    os.remove("zn_frag.pdb")
  f = open("zn_frag.pdb", "w")
  f.write(pdb_in.hierarchy.as_pdb_string(crystal_symmetry=xrs))
  f.close()

if (__name__ == "__main__") :
  exercise()
