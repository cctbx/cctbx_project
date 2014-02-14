
from __future__ import division
from mmtbx.monomer_library import pdb_interpretation
from mmtbx.utils import occupancy_selections
from mmtbx.refinement import occupancies
from mmtbx.command_line import fmodel
from iotbx import file_reader
from libtbx.test_utils import Exception_expected
from libtbx.utils import null_out, Sorry

pdb_raw = """\
CRYST1   21.937    4.866   23.477  90.00 107.08  90.00 P 1 21 1
SCALE1      0.045585  0.000000  0.014006        0.00000
SCALE2      0.000000  0.205508  0.000000        0.00000
SCALE3      0.000000  0.000000  0.044560        0.00000
ATOM      1  N   GLY A   1      -9.056   4.638   6.050  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.058   4.194   4.668  1.00 16.57           C
ATOM      3  C   GLY A   1      -7.993   3.144   4.430  1.00 16.16           C
ATOM      4  O   GLY A   1      -7.521   2.511   5.374  1.00 16.78           O
ATOM      5  N   ASN A   2      -7.616   2.953   3.169  1.00 15.02           N
ATOM      6  CA  ASN A   2      -6.526   2.044   2.840  1.00 14.10           C
ATOM      7  C   ASN A   2      -5.216   2.527   3.434  1.00 13.13           C
ATOM      8  O   ASN A   2      -4.943   3.727   3.466  1.00 11.91           O
ATOM      9  CB  ASN A   2      -6.382   1.888   1.330  1.00 15.38           C
ATOM     10  CG  ASN A   2      -7.632   1.344   0.685  1.00 14.08           C
ATOM     11  OD1 ASN A   2      -8.042   0.216   0.957  1.00 17.46           O
ATOM     12  ND2 ASN A   2      -8.247   2.142  -0.178  1.00 11.72           N
ATOM     13  N   ASN A   3      -4.405   1.583   3.898  1.00 12.26           N
ATOM     14  CA  ASN A   3      -3.172   1.915   4.595  1.00 11.74           C
ATOM     15  C   ASN A   3      -1.922   1.362   3.915  1.00 11.10           C
ATOM     16  O   ASN A   3      -1.816   0.158   3.672  1.00 10.42           O
ATOM     17  CB  ASN A   3      -3.243   1.409   6.039  1.00 12.15           C
ATOM     18  CG  ASN A   3      -2.000   1.749   6.841  1.00 12.82           C
ATOM     19  OD1 ASN A   3      -1.705   2.920   7.082  1.00 15.05           O
ATOM     20  ND2 ASN A   3      -1.272   0.724   7.270  1.00 13.48           N
ATOM     21  N   GLN A   4      -0.987   2.256   3.598  1.00 10.29           N
ATOM     22  CA  GLN A   4       0.361   1.860   3.201  1.00 10.53           C
ATOM     23  C   GLN A   4       1.398   2.605   4.031  1.00 10.24           C
ATOM     24  O   GLN A   4       1.454   3.834   4.025  1.00  8.86           O
ATOM     25  CB  GLN A   4       0.626   2.117   1.712  1.00  9.80           C
ATOM     26  CG  GLN A   4       1.924   1.459   1.221  1.00 10.25           C
ATOM     27  CD  GLN A   4       2.465   2.050  -0.073  1.00 12.43           C
ATOM     28  OE1 GLN A   4       2.674   3.260  -0.178  1.00 14.62           O
ATOM     29  NE2 GLN A   4       2.708   1.192  -1.059  1.00  9.05           N
ATOM     30  N  AGLN A   5       2.202   1.848   4.775  0.62 10.38           N
ATOM     31  CA AGLN A   5       3.288   2.419   5.569  0.62 11.39           C
ATOM     32  C  AGLN A   5       4.638   1.844   5.123  0.62 11.52           C
ATOM     33  O  AGLN A   5       4.824   0.625   5.095  0.62 12.05           O
ATOM     34  CB AGLN A   5       3.046   2.170   7.063  0.62 11.96           C
ATOM     35  CG AGLN A   5       1.854   2.946   7.622  0.62 10.81           C
ATOM     36  CD AGLN A   5       1.361   2.406   8.951  0.62 13.10           C
ATOM     37  OE1AGLN A   5       0.800   1.312   9.019  0.62 10.65           O
ATOM     38  NE2AGLN A   5       1.562   3.175  10.016  0.62 12.30           N
ATOM     39  N  BGLN A   5       2.239   1.858   4.725  0.38 10.38           N
ATOM     40  CA BGLN A   5       3.326   2.476   5.450  0.38 11.39           C
ATOM     41  C  BGLN A   5       4.639   1.850   5.057  0.38 11.52           C
ATOM     42  O  BGLN A   5       4.814   0.627   5.020  0.38 12.05           O
ATOM     43  CB BGLN A   5       3.110   2.331   6.919  0.38 11.96           C
ATOM     44  CG BGLN A   5       2.695   0.980   7.141  0.38 10.81           C
ATOM     45  CD BGLN A   5       2.882   0.618   8.479  0.38 13.10           C
ATOM     46  OE1BGLN A   5       2.538   1.369   9.406  0.38 10.65           O
ATOM     47  NE2BGLN A   5       3.380  -0.597   8.664  0.38 12.30           N
ATOM     48  N   ASN A   6       5.565   2.732   4.753  1.00 11.99           N
ATOM     49  CA  ASN A   6       6.868   2.339   4.280  1.00 12.30           C
ATOM     50  C   ASN A   6       7.881   2.785   5.302  1.00 13.40           C
ATOM     51  O   ASN A   6       8.262   3.954   5.351  1.00 13.92           O
ATOM     52  CB  ASN A   6       7.133   2.954   2.915  1.00 12.13           C
ATOM     53  CG  ASN A   6       5.988   2.721   1.955  1.00 12.77           C
ATOM     54  OD1 ASN A   6       5.795   1.608   1.466  1.00 14.27           O
ATOM     55  ND2 ASN A   6       5.211   3.764   1.690  1.00 10.07           N
ATOM     56  N  ATYR A   7       8.304   1.849   6.146  0.59 14.70           N
ATOM     57  CA ATYR A   7       9.167   2.166   7.280  0.59 15.18           C
ATOM     58  C  ATYR A   7      10.622   2.326   6.868  0.59 15.91           C
ATOM     59  O  ATYR A   7      11.054   1.799   5.844  0.59 15.76           O
ATOM     60  CB ATYR A   7       9.044   1.086   8.356  0.59 15.35           C
ATOM     61  CG ATYR A   7       7.640   0.946   8.887  0.59 14.45           C
ATOM     62  CD1ATYR A   7       6.759   0.027   8.335  0.59 15.68           C
ATOM     63  CD2ATYR A   7       7.187   1.750   9.924  0.59 14.80           C
ATOM     64  CE1ATYR A   7       5.469  -0.098   8.810  0.59 13.46           C
ATOM     65  CE2ATYR A   7       5.899   1.633  10.407  0.59 14.33           C
ATOM     66  CZ ATYR A   7       5.044   0.707   9.845  0.59 15.09           C
ATOM     67  OH ATYR A   7       3.759   0.583  10.319  0.59 14.39           O
ATOM     68  OXTATYR A   7      11.394   2.990   7.558  0.59 17.49           O
ATOM     70  N  BTYR A   7       8.323   1.843   6.116  0.41 14.70           N
ATOM     71  CA BTYR A   7       9.149   2.183   7.247  0.41 15.18           C
ATOM     72  C  BTYR A   7      10.629   2.316   6.861  0.41 15.91           C
ATOM     73  O  BTYR A   7      11.084   1.756   5.864  0.41 15.76           O
ATOM     74  CB BTYR A   7       8.954   1.147   8.348  0.41 15.35           C
ATOM     75  CG BTYR A   7       9.942   1.356   9.417  0.41 14.45           C
ATOM     76  CD1BTYR A   7       9.807   2.381  10.320  0.41 15.68           C
ATOM     77  CD2BTYR A   7      11.054   0.580   9.473  0.41 14.80           C
ATOM     78  CE1BTYR A   7      10.746   2.569  11.248  0.41 13.46           C
ATOM     79  CE2BTYR A   7      11.968   0.749  10.405  0.41 14.33           C
ATOM     80  CZ BTYR A   7      11.858   1.724  11.252  0.41 15.09           C
ATOM     81  OH BTYR A   7      12.921   1.747  12.113  0.41 14.39           O
ATOM     82  OXTBTYR A   7      11.408   3.001   7.529  0.41 17.49           O
TER
HETATM   83  O   HOH A   8      -6.471   5.227   7.124  1.00 22.62           O
HETATM   84  O   HOH A   9      10.431   1.858   3.216  1.00 19.71           O
HETATM   85  O   HOH A  10     -11.286   1.756  -1.468  1.00 17.08           O
HETATM   86  O  AHOH A  11      11.808   4.179   9.970  0.60 23.99           O
HETATM   87  O   HOH A  12      13.605   1.327   9.198  1.00 26.17           O
HETATM   88  O   HOH A  13      -2.749   3.429  10.024  1.00 39.15           O
HETATM   89  O   HOH A  14      -1.500   0.682  10.967  1.00 43.49           O
TER
"""

def prepare_inputs (
    prefix="tst_group_correlated_occupancy",
    create_mtz=False,
    d_min=1.0) :
  pdb_in = "%s_in.pdb" % prefix
  open(pdb_in, "w").write(pdb_raw)
  if (create_mtz) :
    args = [
      pdb_in,
      "high_resolution=%g" % d_min,
      "type=real",
      "label=F",
      "add_sigmas=True",
      "r_free_flags_fraction=0.1",
      "random_seed=12345",
      "output.file_name=%s.mtz" % prefix,
    ]
    fmodel.run(args=args, log=null_out())
  pdb_file = file_reader.any_file(pdb_in)
  hierarchy = pdb_file.file_object.construct_hierarchy()
  xrs = pdb_file.file_object.xray_structure_simple()
  for atom in hierarchy.atoms() :
    atom.b = 5
    if (atom.occ < 1.0) :
      atom.occ = 0.5
  open("%s_start.pdb" % prefix, "w").write(
    hierarchy.as_pdb_string(crystal_symmetry=xrs))

def exercise () :
  prepare_inputs()
  # File #1 (with homogenized occupancies) should work
  # File #2 should fail due to inconsistent occupancies
  pdb_files = [
    "tst_group_correlated_occupancy_start.pdb",
    "tst_group_correlated_occupancy_in.pdb",
  ]
  for i_file, pdb_file in enumerate(pdb_files) :
    processed_pdb_file = pdb_interpretation.run([pdb_file], log=null_out())
    xrs = processed_pdb_file.xray_structure()
    hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
    constraint_groups = occupancy_selections(
      all_chain_proxies=processed_pdb_file.all_chain_proxies,
      xray_structure=xrs)
    try :
      new_groups = occupancies.assemble_constraint_groups_3d(
        xray_structure=xrs,
        pdb_atoms=hierarchy.atoms(),
        constraint_groups=constraint_groups,
        log=null_out())
    except Sorry, s :
      if (i_file == 0) :
        raise
      else :
        assert ("Inconsistent occupancies" in str(s)), str(s)
    else :
      if (i_file == 1) :
        raise Exception_expected
      else :
        assert (len(new_groups) == 1)

if (__name__ == "__main__") :
  exercise()
  print "OK"
