
from __future__ import division
from mmtbx.command_line import build_alt_confs
from mmtbx.command_line import fmodel
from mmtbx.validation import rotalyze
from iotbx import file_reader
from libtbx.utils import null_out
import libtbx.load_env

# derivative of 1yjp, with alternate conformation Asn3 in different rotamer,
# plus Asn2/Gln4 split.  the unit cell 'b' edge has been increased to 6A to
# compensate.
pdb_raw = """\
CRYST1   21.937    6.000   23.477  90.00 107.08  90.00 P 1 21 1      2
ATOM      1  N   GLY A   1      -9.047   4.634   6.066  1.00 16.37           N
ATOM      2  CA  GLY A   1      -9.040   4.191   4.677  1.00 16.17           C
ATOM      3  C   GLY A   1      -7.991   3.118   4.432  1.00 15.46           C
ATOM      4  O   GLY A   1      -7.507   2.493   5.371  1.00 16.71           O
ATOM      5  N  AASN A   2      -7.663   2.878   3.168  0.60 14.94           N
ATOM      6  CA AASN A   2      -6.533   2.012   2.845  0.60 14.28           C
ATOM      7  C  AASN A   2      -5.242   2.518   3.452  0.60 12.87           C
ATOM      8  O  AASN A   2      -5.000   3.723   3.505  0.60 12.55           O
ATOM      9  CB AASN A   2      -6.356   1.880   1.341  0.60 15.48           C
ATOM     10  CG AASN A   2      -7.500   1.158   0.693  0.60 13.90           C
ATOM     11  OD1AASN A   2      -8.052   1.622  -0.304  0.60 17.99           O
ATOM     12  ND2AASN A   2      -7.867   0.008   1.252  0.60 11.61           N
ATOM     13  N  BASN A   2      -7.626   2.899   3.175  0.40 15.11           N
ATOM     14  CA BASN A   2      -6.447   2.085   2.893  0.40 13.97           C
ATOM     15  C  BASN A   2      -5.209   2.731   3.488  0.40 13.04           C
ATOM     16  O  BASN A   2      -5.077   3.953   3.479  0.40 11.55           O
ATOM     17  CB BASN A   2      -6.250   1.906   1.394  0.40 15.07           C
ATOM     18  CG BASN A   2      -7.412   1.207   0.741  0.40 13.81           C
ATOM     19  OD1BASN A   2      -7.953   1.683  -0.257  0.40 17.17           O
ATOM     20  ND2BASN A   2      -7.805   0.067   1.298  0.40 11.52           N
ATOM     21  N  AASN A   3      -4.414   1.582   3.898  0.60 12.35           N
ATOM     22  CA AASN A   3      -3.182   1.908   4.599  0.60 11.68           C
ATOM     23  C  AASN A   3      -1.933   1.346   3.923  0.60 11.09           C
ATOM     24  O  AASN A   3      -1.858   0.142   3.653  0.60 10.18           O
ATOM     25  CB AASN A   3      -3.259   1.407   6.050  0.60 11.62           C
ATOM     26  CG AASN A   3      -2.011   1.746   6.852  0.60 13.00           C
ATOM     27  OD1AASN A   3      -1.704   2.921   7.070  0.60 15.44           O
ATOM     28  ND2AASN A   3      -1.287   0.720   7.295  0.60 12.87           N
ATOM     29  N  BASN A   3      -4.310   1.908   4.015  0.40 12.00           N
ATOM     30  CA BASN A   3      -3.078   2.406   4.614  0.40 11.53           C
ATOM     31  C  BASN A   3      -1.844   1.661   4.111  0.40 11.70           C
ATOM     32  O  BASN A   3      -1.810   0.428   4.100  0.40 10.95           O
ATOM     33  CB BASN A   3      -3.152   2.335   6.142  0.40 12.90           C
ATOM     34  CG BASN A   3      -4.074   3.380   6.732  0.40 13.32           C
ATOM     35  OD1BASN A   3      -4.367   4.398   6.099  0.40 15.04           O
ATOM     36  ND2BASN A   3      -4.528   3.143   7.958  0.40 12.97           N
ATOM     37  N  AGLN A   4      -0.970   2.230   3.649  0.60 10.54           N
ATOM     38  CA AGLN A   4       0.385   1.832   3.248  0.60 10.42           C
ATOM     39  C  AGLN A   4       1.438   2.482   4.154  0.60 11.71           C
ATOM     40  O  AGLN A   4       1.592   3.742   4.128  0.60  8.92           O
ATOM     41  CB AGLN A   4       0.671   2.161   1.771  0.60  9.88           C
ATOM     42  CG AGLN A   4       1.921   1.446   1.228  0.60 10.02           C
ATOM     43  CD AGLN A   4       2.481   2.048  -0.057  0.60 12.86           C
ATOM     44  OE1AGLN A   4       2.716   3.260  -0.143  0.60 14.16           O
ATOM     45  NE2AGLN A   4       2.719   1.195  -1.059  0.60  9.04           N
ATOM     46  N  BGLN A   4      -0.835   2.425   3.699  0.40 10.83           N
ATOM     47  CA BGLN A   4       0.428   1.863   3.233  0.40 10.25           C
ATOM     48  C  BGLN A   4       1.600   2.421   4.037  0.40 10.29           C
ATOM     49  O  BGLN A   4       2.033   3.551   3.826  0.40 10.50           O
ATOM     50  CB BGLN A   4       0.633   2.137   1.738  0.40  9.74           C
ATOM     51  CG BGLN A   4       1.808   1.376   1.129  0.40 10.05           C
ATOM     52  CD BGLN A   4       2.359   2.037  -0.120  0.40 12.79           C
ATOM     53  OE1BGLN A   4       2.503   3.262  -0.175  0.40 15.08           O
ATOM     54  NE2BGLN A   4       2.674   1.228  -1.135  0.40  8.91           N
ATOM     55  N  AGLN A   5       2.125   1.651   4.991  0.60 10.59           N
ATOM     56  CA AGLN A   5       3.274   2.224   5.659  0.60 11.43           C
ATOM     57  C  AGLN A   5       4.582   1.732   5.073  0.60 11.24           C
ATOM     58  O  AGLN A   5       4.749   0.545   4.827  0.60 11.99           O
ATOM     59  CB AGLN A   5       3.229   2.076   7.170  0.60 12.07           C
ATOM     60  CG AGLN A   5       2.235   3.000   7.859  0.60 10.78           C
ATOM     61  CD AGLN A   5       1.562   2.322   9.034  0.60 12.94           C
ATOM     62  OE1AGLN A   5       1.005   1.233   8.899  0.60 10.72           O
ATOM     63  NE2AGLN A   5       1.621   2.959  10.197  0.60 12.32           N
ATOM     55  N  BGLN A   5       2.125   1.651   4.991  0.40 10.59           N
ATOM     56  CA BGLN A   5       3.274   2.224   5.659  0.40 11.43           C
ATOM     57  C  BGLN A   5       4.582   1.732   5.073  0.40 11.24           C
ATOM     58  O  BGLN A   5       4.749   0.545   4.827  0.40 11.99           O
ATOM     59  CB BGLN A   5       3.229   2.076   7.170  0.40 12.07           C
ATOM     60  CG BGLN A   5       2.235   3.000   7.859  0.40 10.78           C
ATOM     61  CD BGLN A   5       1.562   2.322   9.034  0.40 12.94           C
ATOM     62  OE1BGLN A   5       1.005   1.233   8.899  0.40 10.72           O
ATOM     63  NE2BGLN A   5       1.621   2.959  10.197  0.40 12.32           N
ATOM     64  N   ASN A   6       5.508   2.663   4.851  1.00 11.72           N
ATOM     65  CA  ASN A   6       6.825   2.322   4.325  1.00 12.12           C
ATOM     66  C   ASN A   6       7.854   2.763   5.330  1.00 13.15           C
ATOM     67  O   ASN A   6       8.221   3.937   5.380  1.00 13.93           O
ATOM     68  CB  ASN A   6       7.061   3.013   2.994  1.00 11.96           C
ATOM     69  CG  ASN A   6       5.963   2.732   2.005  1.00 12.58           C
ATOM     70  OD1 ASN A   6       5.799   1.604   1.549  1.00 14.01           O
ATOM     71  ND2 ASN A   6       5.192   3.751   1.679  1.00  9.96           N
ATOM     72  N   TYR A   7       8.297   1.822   6.155  1.00 14.62           N
ATOM     73  CA  TYR A   7       9.162   2.146   7.291  1.00 15.04           C
ATOM     74  C   TYR A   7      10.611   2.329   6.888  1.00 15.56           C
ATOM     75  O   TYR A   7      11.046   1.810   5.854  1.00 15.52           O
ATOM     76  CB  TYR A   7       9.056   1.072   8.370  1.00 14.86           C
ATOM     77  CG  TYR A   7       7.657   0.941   8.898  1.00 14.38           C
ATOM     78  CD1 TYR A   7       6.767   0.030   8.334  1.00 15.46           C
ATOM     79  CD2 TYR A   7       7.206   1.753   9.930  1.00 14.37           C
ATOM     80  CE1 TYR A   7       5.476  -0.089   8.802  1.00 13.24           C
ATOM     81  CE2 TYR A   7       5.905   1.644  10.409  1.00 13.84           C
ATOM     82  CZ  TYR A   7       5.049   0.722   9.836  1.00 14.98           C
ATOM     83  OH  TYR A   7       3.767   0.593  10.299  1.00 14.28           O
ATOM     84  OXT TYR A   7      11.361   3.003   7.603  1.00 17.34           O
TER
HETATM   85  O   HOH S   8      -6.473   5.220   7.122  1.00 22.61           O
HETATM   86  O   HOH S   9      10.427   1.864   3.212  1.00 19.32           O
HETATM   87  O   HOH S  10     -11.288   1.762  -1.464  1.00 16.97           O
HETATM   88  O   HOH S  11      11.803   4.188   9.965  1.00 23.89           O
HETATM   89  O   HOH S  12      13.608   1.315   9.196  1.00 26.08           O
HETATM   90  O   HOH S  13      -2.736   3.452  10.015  1.00 38.68           O
HETATM   91  O   HOH S  14      -1.495   0.667  10.978  1.00 44.24           O
TER
END
"""

def exercise () :
  open("tst_build_alt_confs_in.pdb", "w").write(pdb_raw)
  args = [
    "tst_build_alt_confs_in.pdb",
    "high_resolution=1.2",
    "type=real",
    "label=F",
    "add_sigmas=True",
    "r_free_flags_fraction=0.1",
    "random_seed=12345",
    "output.file_name=tst_build_alt_confs.mtz",
  ]
  fmodel.run(args=args, log=null_out())
  pdb_in = file_reader.any_file("tst_build_alt_confs_in.pdb")
  hierarchy = pdb_in.file_object.construct_hierarchy()
  validate = rotalyze.rotalyze(pdb_hierarchy=hierarchy,
    outliers_only=False)
  rota_in = [ (r.id_str(), r.rotamer_name) for r in validate.results ]
  xrs = pdb_in.file_object.xray_structure_simple()
  for chain in hierarchy.only_model().chains() :
    for residue_group in chain.residue_groups() :
      atom_groups = residue_group.atom_groups()
      if (len(atom_groups) == 2) :
        residue_group.remove_atom_group(atom_groups[-1])
        for atom in residue_group.atoms() :
          atom.occ = 1.0
        atom_groups[0].altloc = ''
  assert hierarchy.atoms().extract_occ().all_eq(1.0)
  open("tst_build_alt_confs_start.pdb", "w").write(
    hierarchy.as_pdb_string(crystal_symmetry=xrs))
  args = [
    "tst_build_alt_confs_start.pdb",
    "tst_build_alt_confs.mtz",
    "selection='chain A and resseq 3'",
    "output.file_name=tst_build_alt_confs_out.pdb",
    "expected_occupancy=0.4",
    "window_size=2",
    "nproc=1",
    "rsr_fofc_map_target=False",
    "--verbose",
  ]
  build_alt_confs.run(args=args, out=null_out())
  pdb_out = file_reader.any_file("tst_build_alt_confs_out.pdb")
  hierarchy =  pdb_out.file_object.construct_hierarchy()
  validate = rotalyze.rotalyze(pdb_hierarchy=hierarchy,
    outliers_only=False)
  rota_out = [ (r.id_str(), r.rotamer_name) for r in validate.results ]
  #print rota_in
  #print rota_out
  assert (rota_out == rota_in)

if (__name__ == "__main__") :
  if (not libtbx.env.has_module("probe")) :
    print "probe not available, skipping this test"
  else :
    exercise()
    print "OK"
