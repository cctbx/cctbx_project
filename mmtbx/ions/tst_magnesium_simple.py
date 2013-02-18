
from __future__ import division
from libtbx import easy_run
import os

def exercise () :
  if (os.path.isfile("mg_frag.mtz")) :
    os.remove("mg_frag.mtz")
  write_pdb_input()
  params = """
    high_resolution = 1.5
    r_free_flags_fraction = 0.1
    add_sigmas = True
    pdb_file = mg_frag.pdb
    output {
      label = F
      type = *real complex
      file_name = mg_frag.mtz
    }"""
  open("mg_frag_fmodel.eff", "w").write(params)
  assert (easy_run.fully_buffered("phenix.fmodel mg_frag_fmodel.eff"
    ).raise_if_errors().return_code == 0)
  assert os.path.isfile("mg_frag.mtz")
  from iotbx.file_reader import any_file
  pdb_in = any_file("mg_frag.pdb")
  hierarchy = pdb_in.file_object.construct_hierarchy()
  for chain in hierarchy.models()[0].chains() :
    for residue_group in chain.residue_groups() :
      for atom_group in residue_group.atom_groups() :
        if (atom_group.resname == "MG ") :
          atom_group.resname = "HOH"
          atom = atom_group.atoms()[0]
          atom.name = " O  "
          atom.element = " O"
          atom.charge = ""
          atom.segid = "MG"
          break
  f = open("mg_frag_hoh.pdb", "w")
  f.write(hierarchy.as_pdb_string(pdb_in.file_object.crystal_symmetry()))
  f.close()
  args = ["mg_frag_hoh.pdb", "mg_frag.mtz", "nproc=1", "use_phaser=False",
          "elements=MG"]
  result = easy_run.fully_buffered("mmtbx.water_screen %s" % " ".join(args)
    ).raise_if_errors()
  n_mg = 0
  for line in result.stdout_lines :
    print line
    if ("Probable cation: MG+2" in line) :
      n_mg += 1
  assert n_mg == 2
  print "OK"

def write_pdb_input () :
  import iotbx.pdb.hierarchy
  pdb_in = iotbx.pdb.hierarchy.input(source_info=None, pdb_string="""\
REMARK  fragment derived from 3ip0
ATOM      1  N   ASP A  95     -12.581   4.250  10.804  1.00  5.08           N
ATOM      2  CA  ASP A  95     -11.950   3.273   9.927  1.00  5.04           C
ATOM      3  C   ASP A  95     -11.505   2.089  10.784  1.00  4.97           C
ATOM      4  O   ASP A  95     -10.680   2.251  11.678  1.00  5.43           O
ATOM      5  CB  ASP A  95     -10.737   3.916   9.223  1.00  5.26           C
ATOM      6  CG  ASP A  95     -10.103   3.032   8.160  1.00  5.03           C
ATOM      7  OD1 ASP A  95     -10.879   2.437   7.372  1.00  5.23           O
ATOM      8  OD2 ASP A  95      -8.836   2.949   8.149  1.00  5.44           O1-
TER       9      ASP A  95
ATOM     10  N   ASP A  97      -9.724  -1.416  10.219  1.00  4.47           N
ATOM     11  CA  ASP A  97      -8.951  -2.179   9.250  1.00  4.64           C
ATOM     12  C   ASP A  97      -8.480  -3.469   9.910  1.00  4.67           C
ATOM     13  O   ASP A  97      -7.947  -3.449  11.042  1.00  5.63           O
ATOM     14  CB  ASP A  97      -7.746  -1.380   8.744  1.00  5.20           C
ATOM     15  CG  ASP A  97      -8.132  -0.223   7.828  1.00  5.06           C
ATOM     16  OD1 ASP A  97      -9.317  -0.118   7.445  1.00  5.31           O
ATOM     17  OD2 ASP A  97      -7.197   0.558   7.508  1.00  5.74           O1-
TER      18      ASP A  97
HETATM   19  PA  APC A 171     -11.494   1.693   3.101  1.00  5.93           P
HETATM   20  PB  APC A 171      -8.904   2.838   3.989  1.00  5.40           P
HETATM   21  PG  APC A 171      -6.189   1.883   3.824  1.00  5.44           P
HETATM   22  C5' APC A 171     -11.917  -0.062   1.209  1.00  6.23           C
HETATM   23  O5' APC A 171     -11.907   0.228   2.621  1.00  5.98           O
HETATM   24  C4' APC A 171     -12.005  -1.551   1.015  1.00  6.31           C
HETATM   25  O4' APC A 171     -10.811  -2.160   1.564  1.00  6.21           O
HETATM   26  C3' APC A 171     -13.142  -2.283   1.727  1.00  6.94           C
HETATM   27  O3' APC A 171     -14.374  -2.203   1.097  1.00  8.42           O
HETATM   28  C2' APC A 171     -12.584  -3.692   1.886  1.00  6.56           C
HETATM   29  O2' APC A 171     -13.250  -4.576   0.968  1.00  8.16           O
HETATM   30  C1' APC A 171     -11.120  -3.552   1.530  1.00  6.15           C
HETATM   31  N1  APC A 171      -9.299  -4.552   6.318  1.00  5.61           N
HETATM   32  O1A APC A 171     -12.241   2.714   2.317  1.00  6.59           O
HETATM   33  O1B APC A 171      -9.519   4.202   3.986  1.00  6.49           O
HETATM   34  O1G APC A 171      -5.931   1.987   5.301  1.00  5.83           O
HETATM   35  C2  APC A 171     -10.353  -3.827   5.903  1.00  5.64           C
HETATM   36  O2A APC A 171     -11.698   1.666   4.588  1.00  5.88           O
HETATM   37  O2B APC A 171      -8.855   2.079   5.313  1.00  5.21           O
HETATM   38  O2G APC A 171      -5.035   2.418   2.993  1.00  6.01           O
HETATM   39  C3A APC A 171      -9.734   1.834   2.758  1.00  5.70           C
HETATM   40  N3  APC A 171     -10.798  -3.610   4.683  1.00  5.70           N
HETATM   41  O3B APC A 171      -7.396   2.924   3.462  1.00  5.24           O
HETATM   42  O3G APC A 171      -6.675   0.503   3.400  1.00  6.98           O
HETATM   43  C4  APC A 171     -10.034  -4.246   3.763  1.00  5.51           C
HETATM   44  C5  APC A 171      -8.964  -5.056   4.033  1.00  5.40           C
HETATM   45  C6  APC A 171      -8.569  -5.207   5.373  1.00  5.37           C
HETATM   46  N6  APC A 171      -7.500  -5.912   5.734  1.00  5.71           N
HETATM   47  N7  APC A 171      -8.414  -5.587   2.862  1.00  5.73           N
HETATM   48  C8  APC A 171      -9.167  -5.072   1.920  1.00  6.02           C
HETATM   49  N9  APC A 171     -10.168  -4.254   2.383  1.00  5.71           N
HETATM   50  N1 AHHR A 181      -3.249   8.723   8.653  0.45  3.59           N
HETATM   51  C2 AHHR A 181      -3.170   8.473   9.936  0.45  3.27           C
HETATM   52  N2 AHHR A 181      -2.312   9.188  10.698  0.45  4.00           N
HETATM   53  N3 AHHR A 181      -3.822   7.446  10.541  0.45  3.91           N
HETATM   54  C4 AHHR A 181      -4.648   6.566   9.912  0.45  3.74           C
HETATM   55  O4 AHHR A 181      -5.221   5.635  10.511  0.45  4.44           O
HETATM   56  N5 AHHR A 181      -5.483   5.915   7.737  0.45  3.23           N
HETATM   57  C6 AHHR A 181      -5.582   6.218   6.448  0.45  3.92           C
HETATM   58  C6AAHHR A 181      -6.336   5.270   5.533  0.45  4.40           C
HETATM   59  O6AAHHR A 181      -7.145   4.378   6.216  0.45  4.54           O
HETATM   60  C7 AHHR A 181      -4.896   7.322   5.892  0.45  5.79           C
HETATM   61  N8 AHHR A 181      -4.134   8.150   6.611  0.45  5.53           N
HETATM   62  C9 AHHR A 181      -4.046   7.908   7.922  0.45  3.58           C
HETATM   63  C10AHHR A 181      -4.762   6.791   8.479  0.45  3.12           C
HETATM   64  N1 BHHS A 182      -3.604   8.719   8.748  0.55  6.49           N
HETATM   65  C2 BHHS A 182      -3.478   8.455  10.036  0.55  6.53           C
HETATM   66  N2 BHHS A 182      -2.673   9.243  10.809  0.55  6.01           N
HETATM   67  N3 BHHS A 182      -4.109   7.429  10.645  0.55  5.05           N
HETATM   68  C4 BHHS A 182      -4.896   6.564   9.963  0.55  5.93           C
HETATM   69  O4 BHHS A 182      -5.536   5.719  10.560  0.55  5.29           O
HETATM   70  N5 BHHS A 182      -5.756   5.984   7.764  0.55  7.30           N
HETATM   71  C6 BHHS A 182      -6.698   5.438   5.588  0.55  6.59           C
HETATM   72  C6ABHHS A 182      -5.887   6.282   6.483  0.55  5.89           C
HETATM   73  O6ABHHS A 182      -7.411   4.466   6.183  0.55  6.94           O
HETATM   74  O6BBHHS A 182      -6.769   5.796   4.398  0.55  7.23           O
HETATM   75  C7 BHHS A 182      -5.274   7.421   5.955  0.55  5.00           C
HETATM   76  N8 BHHS A 182      -4.530   8.228   6.695  0.55  4.94           N
HETATM   77  C9 BHHS A 182      -4.383   7.927   7.998  0.55  5.58           C
HETATM   78  C10BHHS A 182      -5.029   6.797   8.531  0.55  6.96           C
HETATM   79  O   HOH A 201     -12.056  -0.140   6.826  1.00  5.82           O
HETATM   80  O   HOH A 202     -10.078  -0.655   4.762  1.00  6.16           O
HETATM   81  O   HOH A 203      -5.851   3.057   8.110  1.00  6.43           O
HETATM   82  O   HOH A 205      -7.904   1.644  10.650  1.00  6.56           O
HETATM   83  O   HOH A 216     -13.376  -2.280   5.233  1.00  7.21           O
HETATM   84  O   HOH A 219      -5.353   2.683  10.917  1.00  8.33           O
HETATM   85  O   HOH A 243      -8.078   5.439   9.608  1.00  9.50           O
HETATM   86  O   HOH A 264     -14.344   0.622   4.874  1.00 14.43           O
HETATM   87  O  AHOH A 499      -7.735  -1.646   4.892  0.50  5.41           O
HETATM   88  O  BHOH A 499      -7.601  -1.751   4.259  0.50  8.56           O
HETATM   89  O  BHOH A 501     -13.472   3.920   6.947  0.56  7.54           O
HETATM   90  O  BHOH A 504      -9.634   5.594   6.224  0.25  5.98           O
HETATM   91  O  CHOH A 504      -9.119   6.099   7.366  0.25  7.90           O
HETATM   92 MG   MG  A 161     -10.440   0.924   6.059  1.00  5.29      MG  MG
HETATM   93 MG   MG  A 162      -7.326   2.422   6.749  1.00  5.40      MG  MG
""")
  xrs = pdb_in.input.xray_structure_simple(cryst1_substitution_buffer_layer=5)
  if (os.path.exists("mg_frag.pdb")) :
    os.remove("mg_frag.pdb")
  f = open("mg_frag.pdb", "w")
  f.write(pdb_in.hierarchy.as_pdb_string(crystal_symmetry=xrs))
  f.close()

if (__name__ == "__main__") :
  exercise()
