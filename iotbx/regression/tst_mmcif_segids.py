from __future__ import absolute_import, division, print_function
import iotbx.pdb
from six.moves import cStringIO as StringIO

pdb_str = """\
ATOM      1  C5'   U     5      97.750 -47.885 -63.217  1.00 83.53      A16S C
ATOM      2  O5'   U     5      98.491 -48.524 -62.188  1.00 79.28      A16S O
ATOM      3  C4'   U     5      96.316 -48.337 -63.206  1.00 86.16      A16S C
ATOM      4  O4'   U     5      95.819 -48.280 -61.845  1.00 86.45      A16S O
ATOM      5  C3'   U     5      95.333 -47.485 -63.995  1.00 87.46      A16S C
ATOM      6  O3'   U     5      95.360 -47.739 -65.396  1.00 87.28      A16S O
ATOM      7  C2'   U     5      94.014 -47.773 -63.286  1.00 84.71      A16S C
ATOM      8  O2'   U     5      93.532 -49.058 -63.634  1.00 85.00      A16S O
ATOM      9  C1'   U     5      94.476 -47.853 -61.835  1.00 84.81      A16S C
ATOM     10  N1    U     5      94.415 -46.549 -61.136  1.00 89.08      A16S N
ATOM     11  C2    U     5      94.346 -46.619 -59.764  1.00 90.22      A16S C
ATOM     12  O2    U     5      94.352 -47.676 -59.164  1.00 86.48      A16S O
ATOM     13  N3    U     5      94.297 -45.413 -59.117  1.00 91.22      A16S N
ATOM     14  C4    U     5      94.296 -44.164 -59.695  1.00 88.42      A16S C
ATOM     15  O4    U     5      94.242 -43.169 -58.977  1.00 87.55      A16S O
ATOM     16  C5    U     5      94.358 -44.166 -61.121  1.00 88.94      A16S C
ATOM     17  C6    U     5      94.408 -45.329 -61.773  1.00 89.00      A16S C
ATOM  99999  P     C   172      78.677  39.674 -73.035  1.00136.18      B23S P
ATOM  99999  C5'   C   172      78.936  39.931 -75.628  1.00129.50      B23S C
ATOM  99999  O5'   C   172      78.085  39.757 -74.507  1.00131.05      B23S O
ATOM  99999  C4'   C   172      78.137  39.999 -76.900  1.00127.67      B23S C
ATOM  99999  O4'   C   172      77.604  38.690 -77.227  1.00128.58      B23S O
ATOM  99999  C3'   C   172      76.907  40.883 -76.846  1.00128.68      B23S C
ATOM  99999  O3'   C   172      77.201  42.258 -76.974  1.00129.72      B23S O
ATOM  99999  C2'   C   172      76.054  40.323 -77.972  1.00130.33      B23S C
ATOM  99999  O2'   C   172      76.544  40.754 -79.232  1.00129.46      B23S O
ATOM  99999  C1'   C   172      76.335  38.830 -77.835  1.00131.16      B23S C
ATOM  99999  N1    C   172      75.328  38.158 -76.984  1.00132.38      B23S N
ATOM  99999  C2    C   172      74.039  37.935 -77.483  1.00130.56      B23S C
ATOM  99999  O2    C   172      73.763  38.302 -78.634  1.00128.13      B23S O
ATOM  99999  N3    C   172      73.131  37.313 -76.697  1.00130.47      B23S N
ATOM  99999  C4    C   172      73.464  36.923 -75.469  1.00131.73      B23S C
ATOM  99999  N4    C   172      72.544  36.315 -74.719  1.00129.68      B23S N
ATOM  99999  C5    C   172      74.765  37.140 -74.941  1.00132.57      B23S C
ATOM  99999  C6    C   172      75.656  37.757 -75.722  1.00131.42      B23S C
ATOM  99999  OP1   C   172      79.530  40.865 -72.814  1.00135.64      B23S O
ATOM  99999  OP2   C   172      77.560  39.387 -72.106  1.00132.48      B23S O1-
ATOM  99999  N   PRO    44      38.238  25.445 -61.632  1.00 58.89      BL28 N
ATOM  99999  CA  PRO    44      38.138  24.574 -62.799  1.00 58.30      BL28 C
ATOM  99999  C   PRO    44      38.531  23.151 -62.439  1.00 58.13      BL28 C
ATOM  99999  O   PRO    44      38.794  22.820 -61.286  1.00 57.98      BL28 O
ATOM  99999  CB  PRO    44      39.136  25.191 -63.776  1.00 56.49      BL28 C
ATOM  99999  CG  PRO    44      40.143  25.804 -62.902  1.00 55.93      BL28 C
ATOM  99999  CD  PRO    44      39.395  26.345 -61.727  1.00 56.81      BL28 C

"""

def exercise(prefix="iotbx_tst_mmcif_segids"):
  """
  Properly output segids when chain ids are empty. Case for Ribosome community.
  After making it work, it should be moved to cctbx_project/iotbx/regression.
  """
  out = StringIO()
  pdb_inp = iotbx.pdb.input(lines=pdb_str, source_info=None)
  h = pdb_inp.construct_hierarchy()
  print("Original:")
  print(h.as_pdb_string())
  assert len(h.only_model().chains()) == 3
  assert [x.id for x in h.only_model().chains()] == [' ', ' ', ' ']
  assert [x.segid for x in h.atoms()] == ['A16S', 'A16S', 'A16S', 'A16S', 'A16S',
      'A16S', 'A16S', 'A16S', 'A16S', 'A16S', 'A16S', 'A16S', 'A16S', 'A16S',
      'A16S', 'A16S', 'A16S', 'B23S', 'B23S', 'B23S', 'B23S', 'B23S', 'B23S',
      'B23S', 'B23S', 'B23S', 'B23S', 'B23S', 'B23S', 'B23S', 'B23S', 'B23S',
      'B23S', 'B23S', 'B23S', 'B23S', 'B23S', 'BL28', 'BL28', 'BL28', 'BL28',
      'BL28', 'BL28', 'BL28']
  cif_block = h.as_cif_block()
  cif = iotbx.cif.model.cif()
  cif['test'] = cif_block
  cif.show(out=out, align_columns=True)
  lines = out.getvalue()
  # print(lines)
  assert "  ATOM  20  OP2  .  C    B23S  172  ?  77.56000   39.38700  -72.10600  1.000  132.48000  O  -1  B  ?  .  OP2  1" in lines

if __name__ == "__main__":
  exercise()
