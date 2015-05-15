from __future__ import division
import sys
import time

from mmtbx.conformation_dependent_library.rdl_database import rdl_database
from mmtbx.conformation_dependent_library import rotamers as rdl

from iotbx import pdb

pdbs = {
  "3sgs" : """
CRYST1    4.821   19.500   21.004  90.00  94.23  90.00 P 1 21 1      2
ATOM      1  N   GLY A   1       2.186  -1.809  10.234  1.00 10.29           N
ATOM      2  CA  GLY A   1       2.320  -0.815   9.141  1.00  9.80           C
ATOM      3  C   GLY A   1       1.418  -1.163   7.979  1.00  9.53           C
ATOM      4  O   GLY A   1       0.509  -1.993   8.106  1.00 10.48           O
ATOM      5  N   ASP A   2       1.661  -0.515   6.852  1.00  8.39           N
ATOM      6  CA  ASP A   2       0.811  -0.643   5.673  1.00  7.85           C
ATOM      7  C   ASP A   2       1.470  -1.495   4.602  1.00  6.99           C
ATOM      8  O   ASP A   2       2.696  -1.595   4.534  1.00  6.68           O
ATOM      9  CB  ASP A   2       0.516   0.739   5.087  1.00  7.92           C
ATOM     10  CG  ASP A   2      -0.103   1.696   6.094  1.00  8.88           C
ATOM     11  OD1 ASP A   2      -0.972   1.260   6.879  1.00  9.26           O
ATOM     12  OD2 ASP A   2       0.269   2.895   6.083  1.00  9.54           O
ATOM     13  N   VAL A   3       0.639  -2.089   3.747  1.00  5.97           N
ATOM     14  CA  VAL A   3       1.107  -2.730   2.539  1.00  5.50           C
ATOM     15  C   VAL A   3       0.524  -1.988   1.346  1.00  5.36           C
ATOM     16  O   VAL A   3      -0.705  -1.905   1.195  1.00  4.65           O
ATOM     17  CB  VAL A   3       0.686  -4.205   2.485  1.00  5.06           C
ATOM     18  CG1 VAL A   3       1.198  -4.839   1.211  1.00  5.22           C
ATOM     19  CG2 VAL A   3       1.197  -4.950   3.727  1.00  5.82           C
ATOM     20  N   ILE A   4       1.395  -1.423   0.522  1.00  5.14           N
ATOM     21  CA  ILE A   4       0.955  -0.615  -0.606  1.00  5.38           C
ATOM     22  C   ILE A   4       1.465  -1.216  -1.907  1.00  5.88           C
ATOM     23  O   ILE A   4       2.673  -1.309  -2.125  1.00  5.24           O
ATOM     24  CB  ILE A   4       1.414   0.853  -0.480  1.00  4.86           C
ATOM     25  CG1 ILE A   4       0.807   1.508   0.765  1.00  6.05           C
ATOM     26  CG2 ILE A   4       0.992   1.652  -1.697  1.00  3.89           C
ATOM     27  CD1 ILE A   4       1.343   2.897   1.039  1.00  7.01           C
ATOM     28  N   GLU A   5       0.522  -1.599  -2.768  1.00  7.64           N
ATOM     29  CA  GLU A   5       0.797  -2.207  -4.071  1.00  9.57           C
ATOM     30  C   GLU A   5       0.245  -1.323  -5.182  1.00 11.31           C
ATOM     31  O   GLU A   5      -0.967  -1.167  -5.324  1.00 11.92           O
ATOM     32  CB  GLU A   5       0.155  -3.597  -4.151  1.00  9.45           C
ATOM     33  CG  GLU A   5       1.035  -4.713  -3.663  1.00 10.16           C
ATOM     34  CD  GLU A   5       2.009  -5.243  -4.706  1.00 10.01           C
ATOM     35  OE1 GLU A   5       2.127  -4.680  -5.822  1.00 10.67           O
ATOM     36  OE2 GLU A   5       2.653  -6.260  -4.403  1.00 10.98           O
ATOM     37  N   VAL A   6       1.129  -0.742  -5.983  1.00 13.27           N
ATOM     38  CA  VAL A   6       0.694   0.201  -7.016  1.00 15.08           C
ATOM     39  C   VAL A   6       1.239  -0.152  -8.414  1.00 16.11           C
ATOM     40  O   VAL A   6       2.388   0.136  -8.755  1.00 17.57           O
ATOM     41  CB  VAL A   6       1.089   1.654  -6.625  1.00 15.51           C
ATOM     42  CG1 VAL A   6       0.230   2.140  -5.472  1.00 15.79           C
ATOM     43  CG2 VAL A   6       2.539   1.722  -6.216  1.00 16.27           C
ATOM     44  OXT VAL A   6       0.534  -0.718  -9.323  1.00 17.23           O
END
""",
  }

def get_geometry_restraints_manager(pdb_filename,
                                    #pdb_inp,
                                    #pdb_hierarchy,
                                    ):
  t0=time.time()
  from mmtbx.monomer_library import server
  from mmtbx.monomer_library import pdb_interpretation
  #lines = pdb_hierarchy.as_pdb_string(
  #  crystal_symmetry=pdb_inp.crystal_symmetry(),
  #  )
  mon_lib_srv = server.server()
  ener_lib = server.ener_lib()
  processed_pdb = pdb_interpretation.process(
    mon_lib_srv,
    ener_lib,
    #raw_records=lines,
    file_name=pdb_filename,
    )
  geometry_restraints_manager = processed_pdb.geometry_restraints_manager()
  print 'time',time.time()-t0
  return geometry_restraints_manager

def run(filename):
  print filename
  if 0:
    print "RDL"
    for aa in sorted(rdl_database):
      print "  %s" % aa
      for key, value in rdl_database[aa].items():
        print "    %s" % key
        for names, values in rdl_database[aa][key].items():
          print "      %s : %s" % (names, values)
    assert 0
  #
  pdb_inp = pdb.input(filename)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  geometry_restraints_manager = get_geometry_restraints_manager(filename)
  pdb_hierarchy.reset_i_seq_if_necessary()
  rdl.update_restraints(pdb_hierarchy,
                        geometry_restraints_manager,
                        verbose=True,
    )
  print "OK"

if __name__=="__main__":
  if len(sys.argv)==1:
    f=file("3sgs.pdb", "wb")
    f.write(pdbs["3sgs"])
    f.close()
    run("3sgs.pdb")
  else:
    run(sys.argv[1])
