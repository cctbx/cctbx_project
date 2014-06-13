
from __future__ import division
from mmtbx.command_line import rna_validate
import mmtbx.validation.rna_validate
from iotbx import file_reader
import iotbx.pdb.hierarchy
from libtbx.easy_pickle import loads, dumps
from libtbx.utils import null_out
import libtbx.load_env
from cStringIO import StringIO
import sys, os

def exercise_rna_validate():
  # fragment from 3g8t
  pdb_raw = """\
ATOM   7975  P     G Q 140      10.347 137.422  73.792  1.00118.69           P
ATOM   7976  OP1   G Q 140       9.348 138.439  74.195  1.00118.63           O
ATOM   7977  OP2   G Q 140      11.208 137.681  72.617  1.00118.60           O
ATOM   7978  O5'   G Q 140      11.286 137.083  75.051  1.00119.34           O
ATOM   7979  C5'   G Q 140      11.015 137.599  76.358  1.00120.43           C
ATOM   7980  C4'   G Q 140      11.617 136.733  77.454  1.00121.27           C
ATOM   7981  O4'   G Q 140      11.723 135.361  77.003  1.00121.51           O
ATOM   7982  C3'   G Q 140      13.038 137.080  77.885  1.00121.72           C
ATOM   7983  O3'   G Q 140      13.050 138.143  78.836  1.00122.07           O
ATOM   7984  C2'   G Q 140      13.517 135.770  78.511  1.00121.87           C
ATOM   7985  O2'   G Q 140      13.160 135.640  79.877  1.00121.91           O
ATOM   7986  C1'   G Q 140      12.793 134.717  77.671  1.00121.84           C
ATOM   7987  N9    G Q 140      13.651 134.055  76.684  1.00121.95           N
ATOM   7988  C8    G Q 140      14.226 134.620  75.568  1.00121.96           C
ATOM   7989  N7    G Q 140      14.944 133.787  74.870  1.00121.97           N
ATOM   7990  C5    G Q 140      14.843 132.588  75.563  1.00121.99           C
ATOM   7991  C6    G Q 140      15.412 131.320  75.281  1.00122.03           C
ATOM   7992  O6    G Q 140      16.142 131.004  74.331  1.00121.99           O
ATOM   7993  N1    G Q 140      15.064 130.361  76.235  1.00122.12           N
ATOM   7994  C2    G Q 140      14.264 130.602  77.331  1.00122.14           C
ATOM   7995  N2    G Q 140      14.033 129.557  78.139  1.00122.06           N
ATOM   7996  N3    G Q 140      13.721 131.787  77.604  1.00122.09           N
ATOM   7997  C4    G Q 140      14.050 132.734  76.684  1.00122.00           C
ATOM   7998  P     A Q 141      14.015 139.407  78.645  1.00122.32           P
ATOM   7999  OP1   A Q 141      13.141 140.571  78.383  1.00122.35           O
ATOM   8000  OP2   A Q 141      15.073 139.069  77.665  1.00122.37           O
ATOM   8001  O5'   A Q 141      14.690 139.587  80.088  1.00122.53           O
ATOM   8002  C5'   A Q 141      15.914 138.934  80.438  1.00122.90           C
ATOM   8003  C4'   A Q 141      15.673 137.727  81.331  1.00123.23           C
ATOM   8004  O4'   A Q 141      16.899 137.378  82.025  1.00123.39           O
ATOM   8005  C3'   A Q 141      14.605 137.909  82.412  1.00123.36           C
ATOM   8006  O3'   A Q 141      13.551 136.969  82.215  1.00123.31           O
ATOM   8007  C2'   A Q 141      15.347 137.686  83.732  1.00123.49           C
ATOM   8008  O2'   A Q 141      14.582 136.976  84.689  1.00123.46           O
ATOM   8009  C1'   A Q 141      16.551 136.857  83.291  1.00123.61           C
ATOM   8010  N9    A Q 141      17.714 136.912  84.186  1.00123.79           N
ATOM   8011  C8    A Q 141      17.797 137.531  85.406  1.00123.80           C
ATOM   8012  N7    A Q 141      18.964 137.415  85.992  1.00123.83           N
ATOM   8013  C5    A Q 141      19.706 136.662  85.099  1.00123.90           C
ATOM   8014  C6    A Q 141      21.037 136.187  85.130  1.00123.95           C
ATOM   8015  N6    A Q 141      21.875 136.421  86.146  1.00123.97           N
ATOM   8016  N1    A Q 141      21.476 135.460  84.076  1.00123.93           N
ATOM   8017  C2    A Q 141      20.634 135.225  83.057  1.00123.90           C
ATOM   8018  N3    A Q 141      19.364 135.620  82.917  1.00123.84           N
ATOM   8019  C4    A Q 141      18.953 136.340  83.979  1.00123.84           C
TER    8020        A Q 141
"""
  open("tst_rna_validate_1.pdb", "w").write(pdb_raw)
  rv = rna_validate.run(args=["tst_rna_validate_1.pdb"], out=null_out())
  assert len(rv.puckers.results) == 1
  pickle_unpickle(rv)
  # derived from 2goz
  regression_pdb = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/pdb2goz_refmac_tls.ent",
    test=os.path.isfile)
  if (regression_pdb is None):
    print "Skipping exercise_regression(): input pdb (pdb2goz_refmac_tls.ent) not available"
    return
  rv = rna_validate.run(args=[regression_pdb], out=null_out())
  assert len(rv.puckers.results) == 2
  assert len(rv.bonds.results) == 2
  assert len(rv.angles.results) == 14
  assert len(rv.suites.results) == 4
  pickle_unpickle(rv)
  pdb_in = file_reader.any_file(regression_pdb)
  result = mmtbx.validation.rna_validate.rna_validation(
    pdb_hierarchy=pdb_in.file_object.construct_hierarchy(),
    geometry_restraints_manager=None,
    params=None)
  pickle_unpickle(result)

def pickle_unpickle (result) :
  result2 = loads(dumps(result))
  out1 = StringIO()
  out2 = StringIO()
  result.show(out=out1)
  result2.show(out=out2)
  assert (out1.getvalue() == out2.getvalue())

def run():
  verbose = "--verbose" in sys.argv[1:]
  if (not libtbx.env.has_module(name="suitename")):
    print \
      "Skipping exercise_rna_validate():" \
      " phenix not available"
  else:
    exercise_rna_validate()
    print "OK"

if (__name__ == "__main__"):
  run()
