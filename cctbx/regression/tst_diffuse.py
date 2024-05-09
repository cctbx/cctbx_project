from __future__ import absolute_import, division, print_function
from libtbx import easy_run
import libtbx.load_env
import iotbx.mtz
from libtbx.test_utils import approx_equal
import os, time

pdb_str = """
CRYST1   50.840   42.770   28.950  90.00  90.00  90.00 P 21 21 21
MODEL        1
ATOM      1  N   MET A   1      27.340  24.430   2.614  1.00  0.00           N
ATOM      2  CA  MET A   1      26.266  25.413   2.842  1.00  0.00           C
ATOM      3  C   MET A   1      26.913  26.639   3.531  1.00  0.00           C
ATOM      4  O   MET A   1      27.886  26.463   4.263  1.00  0.00           O
ATOM      5  CB  MET A   1      25.112  24.880   3.649  1.00  0.00           C
ATOM      6  CG  MET A   1      25.353  24.860   5.134  1.00  0.00           C
ATOM      7  SD  MET A   1      23.930  23.959   5.904  1.00  0.00           S
ATOM      8  CE  MET A   1      24.447  23.984   7.620  1.00  0.00           C
ATOM      9  N   GLN A   2      26.335  27.770   3.258  1.00  0.00           N
ATOM     10  CA  GLN A   2      26.850  29.021   3.898  1.00  0.00           C
ATOM     11  C   GLN A   2      26.100  29.253   5.202  1.00  0.00           C
ATOM     12  O   GLN A   2      24.865  29.024   5.330  1.00  0.00           O
ATOM     13  CB  GLN A   2      26.733  30.148   2.905  1.00  0.00           C
ATOM     14  CG  GLN A   2      26.882  31.546   3.409  1.00  0.00           C
ATOM     15  CD  GLN A   2      26.786  32.562   2.270  1.00  0.00           C
ATOM     16  OE1 GLN A   2      27.783  33.160   1.870  1.00  0.00           O
ATOM     17  NE2 GLN A   2      25.562  32.733   1.806  1.00  0.00           N
ATOM     18  N   ILE A   3      26.849  29.656   6.217  1.00  0.00           N
ATOM     19  CA  ILE A   3      26.235  30.058   7.497  1.00  0.00           C
ATOM     20  C   ILE A   3      26.882  31.428   7.862  1.00  0.00           C
ATOM     21  O   ILE A   3      27.906  31.711   7.264  1.00  0.00           O
ATOM     22  CB  ILE A   3      26.344  29.050   8.645  1.00  0.00           C
ATOM     23  CG1 ILE A   3      27.810  28.748   8.999  1.00  0.00           C
ATOM     24  CG2 ILE A   3      25.491  27.771   8.287  1.00  0.00           C
ATOM     25  CD1 ILE A   3      27.967  28.087  10.417  1.00  0.00           C
ATOM     26  N   PHE A   4      26.214  32.097   8.771  1.00  0.00           N
ATOM     27  CA  PHE A   4      26.772  33.436   9.197  1.00  0.00           C
ATOM     28  C   PHE A   4      27.151  33.362  10.650  1.00  0.00           C
ATOM     29  O   PHE A   4      26.350  32.778  11.395  1.00  0.00           O
ATOM     30  CB  PHE A   4      25.695  34.498   8.946  1.00  0.00           C
ATOM     31  CG  PHE A   4      25.288  34.609   7.499  1.00  0.00           C
ATOM     32  CD1 PHE A   4      24.147  33.966   7.038  1.00  0.00           C
ATOM     33  CD2 PHE A   4      26.136  35.346   6.640  1.00  0.00           C
ATOM     34  CE1 PHE A   4      23.812  34.031   5.677  1.00  0.00           C
ATOM     35  CE2 PHE A   4      25.810  35.392   5.267  1.00  0.00           C
ATOM     36  CZ  PHE A   4      24.620  34.778   4.853  1.00  0.00           C
ATOM     37  N   VAL A   5      28.260  33.943  11.096  1.00  0.00           N
ATOM     38  CA  VAL A   5      28.605  33.965  12.503  1.00  0.00           C
ATOM     39  C   VAL A   5      28.638  35.461  12.900  1.00  0.00           C
ATOM     40  O   VAL A   5      29.522  36.103  12.320  1.00  0.00           O
ATOM     41  CB  VAL A   5      29.963  33.317  12.814  1.00  0.00           C
ATOM     42  CG1 VAL A   5      30.211  33.394  14.304  1.00  0.00           C
ATOM     43  CG2 VAL A   5      29.957  31.838  12.352  1.00  0.00           C
TER
ENDMDL
MODEL        2
ATOM      1  N   MET A   1      27.840  24.430   2.614  1.00  0.00           N
ATOM      2  CA  MET A   1      26.766  25.413   2.842  1.00  0.00           C
ATOM      3  C   MET A   1      27.413  26.639   3.531  1.00  0.00           C
ATOM      4  O   MET A   1      28.386  26.463   4.263  1.00  0.00           O
ATOM      5  CB  MET A   1      25.612  24.880   3.649  1.00  0.00           C
ATOM      6  CG  MET A   1      25.853  24.860   5.134  1.00  0.00           C
ATOM      7  SD  MET A   1      24.430  23.959   5.904  1.00  0.00           S
ATOM      8  CE  MET A   1      24.947  23.984   7.620  1.00  0.00           C
ATOM      9  N   GLN A   2      26.835  27.770   3.258  1.00  0.00           N
ATOM     10  CA  GLN A   2      27.350  29.021   3.898  1.00  0.00           C
ATOM     11  C   GLN A   2      26.600  29.253   5.202  1.00  0.00           C
ATOM     12  O   GLN A   2      25.365  29.024   5.330  1.00  0.00           O
ATOM     13  CB  GLN A   2      27.233  30.148   2.905  1.00  0.00           C
ATOM     14  CG  GLN A   2      27.382  31.546   3.409  1.00  0.00           C
ATOM     15  CD  GLN A   2      27.286  32.562   2.270  1.00  0.00           C
ATOM     16  OE1 GLN A   2      28.283  33.160   1.870  1.00  0.00           O
ATOM     17  NE2 GLN A   2      26.062  32.733   1.806  1.00  0.00           N
ATOM     18  N   ILE A   3      27.349  29.656   6.217  1.00  0.00           N
ATOM     19  CA  ILE A   3      26.735  30.058   7.497  1.00  0.00           C
ATOM     20  C   ILE A   3      27.382  31.428   7.862  1.00  0.00           C
ATOM     21  O   ILE A   3      28.406  31.711   7.264  1.00  0.00           O
ATOM     22  CB  ILE A   3      26.844  29.050   8.645  1.00  0.00           C
ATOM     23  CG1 ILE A   3      28.310  28.748   8.999  1.00  0.00           C
ATOM     24  CG2 ILE A   3      25.991  27.771   8.287  1.00  0.00           C
ATOM     25  CD1 ILE A   3      28.467  28.087  10.417  1.00  0.00           C
ATOM     26  N   PHE A   4      26.714  32.097   8.771  1.00  0.00           N
ATOM     27  CA  PHE A   4      27.272  33.436   9.197  1.00  0.00           C
ATOM     28  C   PHE A   4      27.651  33.362  10.650  1.00  0.00           C
ATOM     29  O   PHE A   4      26.850  32.778  11.395  1.00  0.00           O
ATOM     30  CB  PHE A   4      26.195  34.498   8.946  1.00  0.00           C
ATOM     31  CG  PHE A   4      25.788  34.609   7.499  1.00  0.00           C
ATOM     32  CD1 PHE A   4      24.647  33.966   7.038  1.00  0.00           C
ATOM     33  CD2 PHE A   4      26.636  35.346   6.640  1.00  0.00           C
ATOM     34  CE1 PHE A   4      24.312  34.031   5.677  1.00  0.00           C
ATOM     35  CE2 PHE A   4      26.310  35.392   5.267  1.00  0.00           C
ATOM     36  CZ  PHE A   4      25.120  34.778   4.853  1.00  0.00           C
ATOM     37  N   VAL A   5      28.760  33.943  11.096  1.00  0.00           N
ATOM     38  CA  VAL A   5      29.105  33.965  12.503  1.00  0.00           C
ATOM     39  C   VAL A   5      29.138  35.461  12.900  1.00  0.00           C
ATOM     40  O   VAL A   5      30.022  36.103  12.320  1.00  0.00           O
ATOM     41  CB  VAL A   5      30.463  33.317  12.814  1.00  0.00           C
ATOM     42  CG1 VAL A   5      30.711  33.394  14.304  1.00  0.00           C
ATOM     43  CG2 VAL A   5      30.457  31.838  12.352  1.00  0.00           C
TER
ENDMDL
END
"""
#phenix.diffuse pdb=m.pdb probabilities=0.5,0.5 resolution=4.0 prefix=tst
def exercise():
  fo = open("tst_diffuse.pdb","w")
  print(pdb_str, file=fo)
  fo.close()
  cmd_list = [
    "pdb=tst_diffuse.pdb",
    "probabilities=0.5,0.5",
    "resolution=4.0",
    "prefix=tst_diffuse"]
  file_location = os.path.join(abs(libtbx.env.bin_path), 'phenix.diffuse')
  if os.path.isfile(file_location):
    cmd_list.insert(0, 'phenix.diffuse')
  else:
    import cctbx
    file_location = os.path.join(
      os.path.dirname(cctbx.__file__),
      'command_line',
      'diffuse.py')
    cmd_list.insert(0, 'python %s ' % file_location)
  cmd = " ".join(cmd_list)
  if 0: print(cmd)
  assert not easy_run.call(cmd)
  mas = iotbx.mtz.object(file_name="tst_diffuse.mtz").as_miller_arrays()
  assert len(mas) == 1
  ma = mas[0]
  assert ma.data().size() == 937, ma.data().size()
  assert approx_equal(ma.d_max_min(), (32.728, 4.0), 1.e-3)
  assert approx_equal((ma.data()/1000.).min_max_mean().as_tuple(),
    (0.001, 6.118, 0.35456), 1.e-3)
  assert ma.info().labels == ['I', 'SIGI']
  assert approx_equal(ma.unit_cell().parameters(),
    (50.84,42.77,28.95,90.0,90.0,90.0), 1.e-2)

if __name__ == '__main__':
  t0 = time.time()
  exercise()
  print("Time: %5.2f"%(time.time()-t0))
  print("OK")
