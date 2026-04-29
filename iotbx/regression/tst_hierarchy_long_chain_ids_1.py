from __future__ import absolute_import, division, print_function
import time
import iotbx.pdb
from libtbx.test_utils import assert_lines_in_text
import mmtbx.model
from mmtbx import monomer_library

# ------------------------------------------------------------------------------

# from https://files.rcsb.org/download/4ZOR-assembly1.cif.gz
mmcif_str = '''
data_XXXX
#
loop_
_audit_author.name
_audit_author.pdbx_ordinal
"Asensio, M.A."     1
"Sankaran, B."      2
"Zwart, P.H."       3
"Tullman-Ercek, D." 4

loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.Cartn_x_esd
_atom_site.Cartn_y_esd
_atom_site.Cartn_z_esd
_atom_site.occupancy_esd
_atom_site.B_iso_or_equiv_esd
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM   1      N N    . SER A    1 2   ? 157.548 157.308 62.978  1.00 68.85  ? ? ? ? ? ? 2   SER A    N    1
ATOM   2      C CA   . SER A    1 2   ? 156.320 156.678 63.447  1.00 78.79  ? ? ? ? ? ? 2   SER A    CA   1
ATOM   3      C C    . SER A    1 2   ? 155.227 156.764 62.393  1.00 67.20  ? ? ? ? ? ? 2   SER A    C    1
ATOM   4      O O    . SER A    1 2   ? 155.070 157.792 61.737  1.00 65.10  ? ? ? ? ? ? 2   SER A    O    1
ATOM   5      C CB   . SER A    1 2   ? 155.840 157.336 64.744  1.00 85.23  ? ? ? ? ? ? 2   SER A    CB   1
ATOM   6      O OG   . SER A    1 2   ? 154.643 156.734 65.207  1.00 91.40  ? ? ? ? ? ? 2   SER A    OG   1
ATOM   7      H HA   . SER A    1 2   ? 156.493 155.741 63.628  1.00 94.55  ? ? ? ? ? ? 2   SER A    HA   1
ATOM   8      H HB2  . SER A    1 2   ? 156.527 157.235 65.422  1.00 102.28 ? ? ? ? ? ? 2   SER A    HB2  1
ATOM   9      H HB3  . SER A    1 2   ? 155.675 158.277 64.578  1.00 102.28 ? ? ? ? ? ? 2   SER A    HB3  1
ATOM   10     H HG   . SER A    1 2   ? 154.392 157.104 65.918  1.00 109.68 ? ? ? ? ? ? 2   SER A    HG   1
ATOM   11     N N    . ASN A    1 3   ? 154.459 155.689 62.251  1.00 63.77  ? ? ? ? ? ? 3   ASN A    N    1
ATOM   12     C CA   . ASN A    1 3   ? 153.352 155.668 61.304  1.00 65.83  ? ? ? ? ? ? 3   ASN A    CA   1
ATOM   13     C C    . ASN A    1 3   ? 152.088 156.238 61.943  1.00 55.45  ? ? ? ? ? ? 3   ASN A    C    1
ATOM   14     O O    . ASN A    1 3   ? 151.014 156.237 61.341  1.00 50.39  ? ? ? ? ? ? 3   ASN A    O    1
ATOM   15     C CB   . ASN A    1 3   ? 153.102 154.243 60.804  1.00 82.56  ? ? ? ? ? ? 3   ASN A    CB   1
ATOM   16     C CG   . ASN A    1 3   ? 152.881 153.251 61.938  1.00 97.37  ? ? ? ? ? ? 3   ASN A    CG   1
ATOM   17     O OD1  . ASN A    1 3   ? 152.730 153.636 63.100  1.00 93.43  ? ? ? ? ? ? 3   ASN A    OD1  1
ATOM   18     N ND2  . ASN A    1 3   ? 152.861 151.967 61.602  1.00 102.91 ? ? ? ? ? ? 3   ASN A    ND2  1
ATOM   19     H H    . ASN A    1 3   ? 154.559 154.958 62.693  1.00 76.52  ? ? ? ? ? ? 3   ASN A    H    1
ATOM   20     H HA   . ASN A    1 3   ? 153.578 156.220 60.540  1.00 78.99  ? ? ? ? ? ? 3   ASN A    HA   1
ATOM   21     H HB2  . ASN A    1 3   ? 152.310 154.239 60.243  1.00 99.07  ? ? ? ? ? ? 3   ASN A    HB2  1
ATOM   22     H HB3  . ASN A    1 3   ? 153.872 153.947 60.293  1.00 99.07  ? ? ? ? ? ? 3   ASN A    HB3  1
ATOM   23     H HD21 . ASN A    1 3   ? 152.739 151.366 62.205  1.00 123.50 ? ? ? ? ? ? 3   ASN A    HD21 1
ATOM   24     H HD22 . ASN A    1 3   ? 152.969 151.735 60.781  1.00 123.50 ? ? ? ? ? ? 3   ASN A    HD22 1
ATOM   25     N N    . PHE A    1 4   ? 152.244 156.739 63.165  1.00 40.42  ? ? ? ? ? ? 4   PHE A    N    1
ATOM   26     C CA   . PHE A    1 4   ? 151.163 157.355 63.923  1.00 39.55  ? ? ? ? ? ? 4   PHE A    CA   1
ATOM   27     C C    . PHE A    1 4   ? 151.025 158.825 63.522  1.00 39.94  ? ? ? ? ? ? 4   PHE A    C    1
ATOM   28     O O    . PHE A    1 4   ? 151.514 159.719 64.217  1.00 40.63  ? ? ? ? ? ? 4   PHE A    O    1
ATOM   29     C CB   . PHE A    1 4   ? 151.452 157.193 65.421  1.00 39.64  ? ? ? ? ? ? 4   PHE A    CB   1
ATOM   30     C CG   . PHE A    1 4   ? 150.306 157.558 66.329  1.00 38.70  ? ? ? ? ? ? 4   PHE A    CG   1
ATOM   31     C CD1  . PHE A    1 4   ? 149.020 157.744 65.846  1.00 37.96  ? ? ? ? ? ? 4   PHE A    CD1  1
ATOM   32     C CD2  . PHE A    1 4   ? 150.526 157.687 67.689  1.00 38.91  ? ? ? ? ? ? 4   PHE A    CD2  1
ATOM   33     C CE1  . PHE A    1 4   ? 147.985 158.072 66.705  1.00 37.01  ? ? ? ? ? ? 4   PHE A    CE1  1
ATOM   34     C CE2  . PHE A    1 4   ? 149.499 158.013 68.549  1.00 38.16  ? ? ? ? ? ? 4   PHE A    CE2  1
ATOM   35     C CZ   . PHE A    1 4   ? 148.226 158.201 68.059  1.00 48.00  ? ? ? ? ? ? 4   PHE A    CZ   1
ATOM   36     H H    . PHE A    1 4   ? 152.992 156.733 63.589  1.00 48.51  ? ? ? ? ? ? 4   PHE A    H    1
ATOM   37     H HA   . PHE A    1 4   ? 150.329 156.902 63.721  1.00 47.46  ? ? ? ? ? ? 4   PHE A    HA   1
ATOM   38     H HB2  . PHE A    1 4   ? 151.679 156.266 65.594  1.00 47.57  ? ? ? ? ? ? 4   PHE A    HB2  1
ATOM   39     H HB3  . PHE A    1 4   ? 152.203 157.761 65.654  1.00 47.57  ? ? ? ? ? ? 4   PHE A    HB3  1
ATOM   40     H HD1  . PHE A    1 4   ? 148.853 157.658 64.935  1.00 45.56  ? ? ? ? ? ? 4   PHE A    HD1  1
ATOM   41     H HD2  . PHE A    1 4   ? 151.383 157.561 68.028  1.00 46.69  ? ? ? ? ? ? 4   PHE A    HD2  1
ATOM   42     H HE1  . PHE A    1 4   ? 147.126 158.199 66.372  1.00 44.41  ? ? ? ? ? ? 4   PHE A    HE1  1
ATOM   43     H HE2  . PHE A    1 4   ? 149.664 158.100 69.460  1.00 45.79  ? ? ? ? ? ? 4   PHE A    HE2  1
ATOM   44     H HZ   . PHE A    1 4   ? 147.532 158.422 68.637  1.00 57.59  ? ? ? ? ? ? 4   PHE A    HZ   1
ATOM   45     N N    . THR A    1 5   ? 150.358 159.064 62.394  1.00 39.85  ? ? ? ? ? ? 5   THR A    N    1
ATOM   46     C CA   . THR A    1 5   ? 150.317 160.394 61.787  1.00 40.09  ? ? ? ? ? ? 5   THR A    CA   1
ATOM   47     C C    . THR A    1 5   ? 148.902 160.860 61.432  1.00 39.18  ? ? ? ? ? ? 5   THR A    C    1
ATOM   48     O O    . THR A    1 5   ? 147.972 160.058 61.325  1.00 38.69  ? ? ? ? ? ? 5   THR A    O    1
ATOM   49     C CB   . THR A    1 5   ? 151.196 160.436 60.514  1.00 64.25  ? ? ? ? ? ? 5   THR A    CB   1
ATOM   50     O OG1  . THR A    1 5   ? 151.389 161.794 60.102  1.00 64.23  ? ? ? ? ? ? 5   THR A    OG1  1
ATOM   51     C CG2  . THR A    1 5   ? 150.563 159.643 59.376  1.00 40.46  ? ? ? ? ? ? 5   THR A    CG2  1
ATOM   52     H H    . THR A    1 5   ? 149.917 158.468 61.958  1.00 47.81  ? ? ? ? ? ? 5   THR A    H    1
ATOM   53     H HA   . THR A    1 5   ? 150.685 161.031 62.418  1.00 48.11  ? ? ? ? ? ? 5   THR A    HA   1
ATOM   54     H HB   . THR A    1 5   ? 152.059 160.040 60.712  1.00 77.10  ? ? ? ? ? ? 5   THR A    HB   1
ATOM   55     H HG1  . THR A    1 5   ? 151.866 161.820 59.411  1.00 77.08  ? ? ? ? ? ? 5   THR A    HG1  1
ATOM   56     H HG21 . THR A    1 5   ? 151.129 159.683 58.590  1.00 48.55  ? ? ? ? ? ? 5   THR A    HG21 1
ATOM   57     H HG22 . THR A    1 5   ? 150.453 158.716 59.639  1.00 48.55  ? ? ? ? ? ? 5   THR A    HG22 1
ATOM   58     H HG23 . THR A    1 5   ? 149.694 160.014 59.158  1.00 48.55  ? ? ? ? ? ? 5   THR A    HG23 1
ATOM   9413   N N    . SER A-2  1 2   ? 70.577  70.817  62.978  1.00 68.85  ? ? ? ? ? ? 2   SER A-2  N    1
ATOM   9414   C CA   . SER A-2  1 2   ? 71.805  71.447  63.447  1.00 78.79  ? ? ? ? ? ? 2   SER A-2  CA   1
ATOM   9415   C C    . SER A-2  1 2   ? 72.898  71.361  62.393  1.00 67.20  ? ? ? ? ? ? 2   SER A-2  C    1
ATOM   9416   O O    . SER A-2  1 2   ? 73.055  70.333  61.737  1.00 65.10  ? ? ? ? ? ? 2   SER A-2  O    1
ATOM   9417   C CB   . SER A-2  1 2   ? 72.285  70.789  64.744  1.00 85.23  ? ? ? ? ? ? 2   SER A-2  CB   1
ATOM   9418   O OG   . SER A-2  1 2   ? 73.482  71.391  65.207  1.00 91.40  ? ? ? ? ? ? 2   SER A-2  OG   1
ATOM   9419   H HA   . SER A-2  1 2   ? 71.632  72.384  63.628  1.00 94.55  ? ? ? ? ? ? 2   SER A-2  HA   1
ATOM   9420   H HB2  . SER A-2  1 2   ? 71.598  70.890  65.422  1.00 102.28 ? ? ? ? ? ? 2   SER A-2  HB2  1
ATOM   9421   H HB3  . SER A-2  1 2   ? 72.450  69.848  64.578  1.00 102.28 ? ? ? ? ? ? 2   SER A-2  HB3  1
ATOM   9422   H HG   . SER A-2  1 2   ? 73.733  71.021  65.918  1.00 109.68 ? ? ? ? ? ? 2   SER A-2  HG   1
ATOM   9423   N N    . ASN A-2  1 3   ? 73.666  72.436  62.251  1.00 63.77  ? ? ? ? ? ? 3   ASN A-2  N    1
ATOM   9424   C CA   . ASN A-2  1 3   ? 74.773  72.457  61.304  1.00 65.83  ? ? ? ? ? ? 3   ASN A-2  CA   1
ATOM   9425   C C    . ASN A-2  1 3   ? 76.037  71.887  61.943  1.00 55.45  ? ? ? ? ? ? 3   ASN A-2  C    1
ATOM   9426   O O    . ASN A-2  1 3   ? 77.111  71.888  61.341  1.00 50.39  ? ? ? ? ? ? 3   ASN A-2  O    1
ATOM   9427   C CB   . ASN A-2  1 3   ? 75.023  73.882  60.804  1.00 82.56  ? ? ? ? ? ? 3   ASN A-2  CB   1
ATOM   9428   C CG   . ASN A-2  1 3   ? 75.244  74.874  61.938  1.00 97.37  ? ? ? ? ? ? 3   ASN A-2  CG   1
ATOM   9429   O OD1  . ASN A-2  1 3   ? 75.395  74.489  63.100  1.00 93.43  ? ? ? ? ? ? 3   ASN A-2  OD1  1
ATOM   9430   N ND2  . ASN A-2  1 3   ? 75.264  76.158  61.602  1.00 102.91 ? ? ? ? ? ? 3   ASN A-2  ND2  1
ATOM   9431   H H    . ASN A-2  1 3   ? 73.566  73.167  62.693  1.00 76.52  ? ? ? ? ? ? 3   ASN A-2  H    1
ATOM   9432   H HA   . ASN A-2  1 3   ? 74.547  71.905  60.540  1.00 78.99  ? ? ? ? ? ? 3   ASN A-2  HA   1
ATOM   9433   H HB2  . ASN A-2  1 3   ? 75.815  73.886  60.243  1.00 99.07  ? ? ? ? ? ? 3   ASN A-2  HB2  1
ATOM   9434   H HB3  . ASN A-2  1 3   ? 74.253  74.178  60.293  1.00 99.07  ? ? ? ? ? ? 3   ASN A-2  HB3  1
ATOM   9435   H HD21 . ASN A-2  1 3   ? 75.386  76.759  62.205  1.00 123.50 ? ? ? ? ? ? 3   ASN A-2  HD21 1
ATOM   9436   H HD22 . ASN A-2  1 3   ? 75.156  76.390  60.781  1.00 123.50 ? ? ? ? ? ? 3   ASN A-2  HD22 1
ATOM   9437   N N    . PHE A-2  1 4   ? 75.881  71.386  63.165  1.00 40.42  ? ? ? ? ? ? 4   PHE A-2  N    1
ATOM   9438   C CA   . PHE A-2  1 4   ? 76.962  70.770  63.923  1.00 39.55  ? ? ? ? ? ? 4   PHE A-2  CA   1
ATOM   9439   C C    . PHE A-2  1 4   ? 77.100  69.300  63.522  1.00 39.94  ? ? ? ? ? ? 4   PHE A-2  C    1
ATOM   9440   O O    . PHE A-2  1 4   ? 76.611  68.406  64.217  1.00 40.63  ? ? ? ? ? ? 4   PHE A-2  O    1
ATOM   9441   C CB   . PHE A-2  1 4   ? 76.673  70.932  65.421  1.00 39.64  ? ? ? ? ? ? 4   PHE A-2  CB   1
ATOM   9442   C CG   . PHE A-2  1 4   ? 77.819  70.567  66.329  1.00 38.70  ? ? ? ? ? ? 4   PHE A-2  CG   1
ATOM   9443   C CD1  . PHE A-2  1 4   ? 79.105  70.381  65.846  1.00 37.96  ? ? ? ? ? ? 4   PHE A-2  CD1  1
ATOM   9444   C CD2  . PHE A-2  1 4   ? 77.599  70.438  67.689  1.00 38.91  ? ? ? ? ? ? 4   PHE A-2  CD2  1
ATOM   9445   C CE1  . PHE A-2  1 4   ? 80.140  70.053  66.705  1.00 37.01  ? ? ? ? ? ? 4   PHE A-2  CE1  1
ATOM   9446   C CE2  . PHE A-2  1 4   ? 78.626  70.112  68.549  1.00 38.16  ? ? ? ? ? ? 4   PHE A-2  CE2  1
ATOM   9447   C CZ   . PHE A-2  1 4   ? 79.899  69.924  68.059  1.00 48.00  ? ? ? ? ? ? 4   PHE A-2  CZ   1
ATOM   9448   H H    . PHE A-2  1 4   ? 75.133  71.392  63.589  1.00 48.51  ? ? ? ? ? ? 4   PHE A-2  H    1
ATOM   9449   H HA   . PHE A-2  1 4   ? 77.796  71.223  63.721  1.00 47.46  ? ? ? ? ? ? 4   PHE A-2  HA   1
ATOM   9450   H HB2  . PHE A-2  1 4   ? 76.446  71.859  65.594  1.00 47.57  ? ? ? ? ? ? 4   PHE A-2  HB2  1
ATOM   9451   H HB3  . PHE A-2  1 4   ? 75.922  70.364  65.654  1.00 47.57  ? ? ? ? ? ? 4   PHE A-2  HB3  1
ATOM   9452   H HD1  . PHE A-2  1 4   ? 79.272  70.467  64.935  1.00 45.56  ? ? ? ? ? ? 4   PHE A-2  HD1  1
ATOM   9453   H HD2  . PHE A-2  1 4   ? 76.742  70.564  68.028  1.00 46.69  ? ? ? ? ? ? 4   PHE A-2  HD2  1
ATOM   9454   H HE1  . PHE A-2  1 4   ? 80.999  69.926  66.372  1.00 44.41  ? ? ? ? ? ? 4   PHE A-2  HE1  1
ATOM   9455   H HE2  . PHE A-2  1 4   ? 78.461  70.025  69.460  1.00 45.79  ? ? ? ? ? ? 4   PHE A-2  HE2  1
ATOM   9456   H HZ   . PHE A-2  1 4   ? 80.593  69.703  68.637  1.00 57.59  ? ? ? ? ? ? 4   PHE A-2  HZ   1
ATOM   9457   N N    . THR A-2  1 5   ? 77.767  69.061  62.394  1.00 39.85  ? ? ? ? ? ? 5   THR A-2  N    1
ATOM   9458   C CA   . THR A-2  1 5   ? 77.808  67.731  61.787  1.00 40.09  ? ? ? ? ? ? 5   THR A-2  CA   1
ATOM   9459   C C    . THR A-2  1 5   ? 79.223  67.265  61.432  1.00 39.18  ? ? ? ? ? ? 5   THR A-2  C    1
ATOM   9460   O O    . THR A-2  1 5   ? 80.153  68.067  61.325  1.00 38.69  ? ? ? ? ? ? 5   THR A-2  O    1
ATOM   9461   C CB   . THR A-2  1 5   ? 76.929  67.689  60.514  1.00 64.25  ? ? ? ? ? ? 5   THR A-2  CB   1
ATOM   9462   O OG1  . THR A-2  1 5   ? 76.736  66.331  60.102  1.00 64.23  ? ? ? ? ? ? 5   THR A-2  OG1  1
ATOM   9463   C CG2  . THR A-2  1 5   ? 77.562  68.482  59.376  1.00 40.46  ? ? ? ? ? ? 5   THR A-2  CG2  1
ATOM   9464   H H    . THR A-2  1 5   ? 78.208  69.657  61.958  1.00 47.81  ? ? ? ? ? ? 5   THR A-2  H    1
ATOM   9465   H HA   . THR A-2  1 5   ? 77.440  67.094  62.418  1.00 48.11  ? ? ? ? ? ? 5   THR A-2  HA   1
ATOM   9466   H HB   . THR A-2  1 5   ? 76.066  68.085  60.712  1.00 77.10  ? ? ? ? ? ? 5   THR A-2  HB   1
ATOM   9467   H HG1  . THR A-2  1 5   ? 76.259  66.305  59.411  1.00 77.08  ? ? ? ? ? ? 5   THR A-2  HG1  1
ATOM   9468   H HG21 . THR A-2  1 5   ? 76.996  68.442  58.590  1.00 48.55  ? ? ? ? ? ? 5   THR A-2  HG21 1
ATOM   9469   H HG22 . THR A-2  1 5   ? 77.672  69.409  59.639  1.00 48.55  ? ? ? ? ? ? 5   THR A-2  HG22 1
ATOM   9470   H HG23 . THR A-2  1 5   ? 78.431  68.111  59.158  1.00 48.55  ? ? ? ? ? ? 5   THR A-2  HG23 1
ATOM   37649  N N    . SER A-5  1 2   ? 62.978  157.548 157.308 1.00 68.85  ? ? ? ? ? ? 2   SER A-5  N    1
ATOM   37650  C CA   . SER A-5  1 2   ? 63.447  156.320 156.678 1.00 78.79  ? ? ? ? ? ? 2   SER A-5  CA   1
ATOM   37651  C C    . SER A-5  1 2   ? 62.393  155.227 156.764 1.00 67.20  ? ? ? ? ? ? 2   SER A-5  C    1
ATOM   37652  O O    . SER A-5  1 2   ? 61.737  155.070 157.792 1.00 65.10  ? ? ? ? ? ? 2   SER A-5  O    1
ATOM   37653  C CB   . SER A-5  1 2   ? 64.744  155.840 157.336 1.00 85.23  ? ? ? ? ? ? 2   SER A-5  CB   1
ATOM   37654  O OG   . SER A-5  1 2   ? 65.207  154.643 156.734 1.00 91.40  ? ? ? ? ? ? 2   SER A-5  OG   1
ATOM   37655  H HA   . SER A-5  1 2   ? 63.628  156.493 155.741 1.00 94.55  ? ? ? ? ? ? 2   SER A-5  HA   1
ATOM   37656  H HB2  . SER A-5  1 2   ? 65.422  156.527 157.235 1.00 102.28 ? ? ? ? ? ? 2   SER A-5  HB2  1
ATOM   37657  H HB3  . SER A-5  1 2   ? 64.578  155.675 158.277 1.00 102.28 ? ? ? ? ? ? 2   SER A-5  HB3  1
ATOM   37658  H HG   . SER A-5  1 2   ? 65.918  154.392 157.104 1.00 109.68 ? ? ? ? ? ? 2   SER A-5  HG   1
ATOM   37659  N N    . ASN A-5  1 3   ? 62.251  154.459 155.689 1.00 63.77  ? ? ? ? ? ? 3   ASN A-5  N    1
ATOM   37660  C CA   . ASN A-5  1 3   ? 61.304  153.352 155.668 1.00 65.83  ? ? ? ? ? ? 3   ASN A-5  CA   1
ATOM   37661  C C    . ASN A-5  1 3   ? 61.943  152.088 156.238 1.00 55.45  ? ? ? ? ? ? 3   ASN A-5  C    1
ATOM   37662  O O    . ASN A-5  1 3   ? 61.341  151.014 156.237 1.00 50.39  ? ? ? ? ? ? 3   ASN A-5  O    1
ATOM   37663  C CB   . ASN A-5  1 3   ? 60.804  153.102 154.243 1.00 82.56  ? ? ? ? ? ? 3   ASN A-5  CB   1
ATOM   37664  C CG   . ASN A-5  1 3   ? 61.938  152.881 153.251 1.00 97.37  ? ? ? ? ? ? 3   ASN A-5  CG   1
ATOM   37665  O OD1  . ASN A-5  1 3   ? 63.100  152.730 153.636 1.00 93.43  ? ? ? ? ? ? 3   ASN A-5  OD1  1
ATOM   37666  N ND2  . ASN A-5  1 3   ? 61.602  152.861 151.967 1.00 102.91 ? ? ? ? ? ? 3   ASN A-5  ND2  1
ATOM   37667  H H    . ASN A-5  1 3   ? 62.693  154.559 154.958 1.00 76.52  ? ? ? ? ? ? 3   ASN A-5  H    1
ATOM   37668  H HA   . ASN A-5  1 3   ? 60.540  153.578 156.220 1.00 78.99  ? ? ? ? ? ? 3   ASN A-5  HA   1
ATOM   37669  H HB2  . ASN A-5  1 3   ? 60.243  152.310 154.239 1.00 99.07  ? ? ? ? ? ? 3   ASN A-5  HB2  1
ATOM   37670  H HB3  . ASN A-5  1 3   ? 60.293  153.872 153.947 1.00 99.07  ? ? ? ? ? ? 3   ASN A-5  HB3  1
ATOM   37671  H HD21 . ASN A-5  1 3   ? 62.205  152.739 151.366 1.00 123.50 ? ? ? ? ? ? 3   ASN A-5  HD21 1
ATOM   37672  H HD22 . ASN A-5  1 3   ? 60.781  152.969 151.735 1.00 123.50 ? ? ? ? ? ? 3   ASN A-5  HD22 1
ATOM   37673  N N    . PHE A-5  1 4   ? 63.165  152.244 156.739 1.00 40.42  ? ? ? ? ? ? 4   PHE A-5  N    1
ATOM   37674  C CA   . PHE A-5  1 4   ? 63.923  151.163 157.355 1.00 39.55  ? ? ? ? ? ? 4   PHE A-5  CA   1
ATOM   37675  C C    . PHE A-5  1 4   ? 63.522  151.025 158.825 1.00 39.94  ? ? ? ? ? ? 4   PHE A-5  C    1
ATOM   37676  O O    . PHE A-5  1 4   ? 64.217  151.514 159.719 1.00 40.63  ? ? ? ? ? ? 4   PHE A-5  O    1
ATOM   37677  C CB   . PHE A-5  1 4   ? 65.421  151.452 157.193 1.00 39.64  ? ? ? ? ? ? 4   PHE A-5  CB   1
ATOM   37678  C CG   . PHE A-5  1 4   ? 66.329  150.306 157.558 1.00 38.70  ? ? ? ? ? ? 4   PHE A-5  CG   1
ATOM   37679  C CD1  . PHE A-5  1 4   ? 65.846  149.020 157.744 1.00 37.96  ? ? ? ? ? ? 4   PHE A-5  CD1  1
ATOM   37680  C CD2  . PHE A-5  1 4   ? 67.689  150.526 157.687 1.00 38.91  ? ? ? ? ? ? 4   PHE A-5  CD2  1
ATOM   37681  C CE1  . PHE A-5  1 4   ? 66.705  147.985 158.072 1.00 37.01  ? ? ? ? ? ? 4   PHE A-5  CE1  1
ATOM   37682  C CE2  . PHE A-5  1 4   ? 68.549  149.499 158.013 1.00 38.16  ? ? ? ? ? ? 4   PHE A-5  CE2  1
ATOM   37683  C CZ   . PHE A-5  1 4   ? 68.059  148.226 158.201 1.00 48.00  ? ? ? ? ? ? 4   PHE A-5  CZ   1
ATOM   37684  H H    . PHE A-5  1 4   ? 63.589  152.992 156.733 1.00 48.51  ? ? ? ? ? ? 4   PHE A-5  H    1
ATOM   37685  H HA   . PHE A-5  1 4   ? 63.721  150.329 156.902 1.00 47.46  ? ? ? ? ? ? 4   PHE A-5  HA   1
ATOM   37686  H HB2  . PHE A-5  1 4   ? 65.594  151.679 156.266 1.00 47.57  ? ? ? ? ? ? 4   PHE A-5  HB2  1
ATOM   37687  H HB3  . PHE A-5  1 4   ? 65.654  152.203 157.761 1.00 47.57  ? ? ? ? ? ? 4   PHE A-5  HB3  1
ATOM   37688  H HD1  . PHE A-5  1 4   ? 64.935  148.853 157.658 1.00 45.56  ? ? ? ? ? ? 4   PHE A-5  HD1  1
ATOM   37689  H HD2  . PHE A-5  1 4   ? 68.028  151.383 157.561 1.00 46.69  ? ? ? ? ? ? 4   PHE A-5  HD2  1
ATOM   37690  H HE1  . PHE A-5  1 4   ? 66.372  147.126 158.199 1.00 44.41  ? ? ? ? ? ? 4   PHE A-5  HE1  1
ATOM   37691  H HE2  . PHE A-5  1 4   ? 69.460  149.664 158.100 1.00 45.79  ? ? ? ? ? ? 4   PHE A-5  HE2  1
ATOM   37692  H HZ   . PHE A-5  1 4   ? 68.637  147.532 158.422 1.00 57.59  ? ? ? ? ? ? 4   PHE A-5  HZ   1
ATOM   37693  N N    . THR A-5  1 5   ? 62.394  150.358 159.064 1.00 39.85  ? ? ? ? ? ? 5   THR A-5  N    1
ATOM   37694  C CA   . THR A-5  1 5   ? 61.787  150.317 160.394 1.00 40.09  ? ? ? ? ? ? 5   THR A-5  CA   1
ATOM   37695  C C    . THR A-5  1 5   ? 61.432  148.902 160.860 1.00 39.18  ? ? ? ? ? ? 5   THR A-5  C    1
ATOM   37696  O O    . THR A-5  1 5   ? 61.325  147.972 160.058 1.00 38.69  ? ? ? ? ? ? 5   THR A-5  O    1
ATOM   37697  C CB   . THR A-5  1 5   ? 60.514  151.196 160.436 1.00 64.25  ? ? ? ? ? ? 5   THR A-5  CB   1
ATOM   37698  O OG1  . THR A-5  1 5   ? 60.102  151.389 161.794 1.00 64.23  ? ? ? ? ? ? 5   THR A-5  OG1  1
ATOM   37699  C CG2  . THR A-5  1 5   ? 59.376  150.563 159.643 1.00 40.46  ? ? ? ? ? ? 5   THR A-5  CG2  1
ATOM   37700  H H    . THR A-5  1 5   ? 61.958  149.917 158.468 1.00 47.81  ? ? ? ? ? ? 5   THR A-5  H    1
ATOM   37701  H HA   . THR A-5  1 5   ? 62.418  150.685 161.031 1.00 48.11  ? ? ? ? ? ? 5   THR A-5  HA   1
ATOM   37702  H HB   . THR A-5  1 5   ? 60.712  152.059 160.040 1.00 77.10  ? ? ? ? ? ? 5   THR A-5  HB   1
ATOM   37703  H HG1  . THR A-5  1 5   ? 59.411  151.866 161.820 1.00 77.08  ? ? ? ? ? ? 5   THR A-5  HG1  1
ATOM   37704  H HG21 . THR A-5  1 5   ? 58.590  151.129 159.683 1.00 48.55  ? ? ? ? ? ? 5   THR A-5  HG21 1
ATOM   37705  H HG22 . THR A-5  1 5   ? 59.639  150.453 158.716 1.00 48.55  ? ? ? ? ? ? 5   THR A-5  HG22 1
ATOM   37706  H HG23 . THR A-5  1 5   ? 59.158  149.694 160.014 1.00 48.55  ? ? ? ? ? ? 5   THR A-5  HG23 1
'''

norm_str = """
ATOM      1  N   GLY A   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      3  C   GLY A   1      -8.015   3.140   4.419  1.00 16.16           C
ATOM      4  O   GLY A   1      -7.523   2.521   5.381  1.00 16.78           O
ATOM      5  N   ASN A   2      -7.656   2.923   3.155  1.00 15.02           N
ATOM      6  CA  ASN A   2      -6.522   2.038   2.831  1.00 14.10           C
ATOM      7  C   ASN A   2      -5.241   2.537   3.427  1.00 13.13           C
ATOM      8  O   ASN A   2      -4.978   3.742   3.426  1.00 11.91           O
ATOM      9  CB  ASN A   2      -6.346   1.881   1.341  1.00 15.38           C
ATOM     10  CG  ASN A   2      -7.584   1.342   0.692  1.00 14.08           C
ATOM     11  OD1 ASN A   2      -8.025   0.227   1.016  1.00 17.46           O
ATOM     12  ND2 ASN A   2      -8.204   2.155  -0.169  1.00 11.72           N
ATOM     13  N   ASN A   3      -4.438   1.590   3.905  1.00 12.26           N
ATOM     14  CA  ASN A   3      -3.193   1.904   4.589  1.00 11.74           C
ATOM     15  C   ASN A   3      -1.955   1.332   3.895  1.00 11.10           C
ATOM     16  O   ASN A   3      -1.872   0.119   3.648  1.00 10.42           O
ATOM     17  CB  ASN A   3      -3.259   1.378   6.042  1.00 12.15           C
ATOM     18  CG  ASN A   3      -2.006   1.739   6.861  1.00 12.82           C
ATOM     19  OD1 ASN A   3      -1.702   2.925   7.072  1.00 15.05           O
"""

def test1():
  """
  Test correct reading of long chain ids from mmCIF file into hierarchy.
  Testing output to pdb and cif.
  """
  inp = iotbx.pdb.input(lines=mmcif_str.split("\n"), source_info=None)
  h = inp.construct_hierarchy()
  chains = list(h.overall_counts().chain_ids.keys())
  chains.sort()
  # print(chains)
  answer = ['A', 'A-2', 'A-5']
  assert (chains == answer), '%s %s' % (chains, answer)

  o_pdb_str = h.as_pdb_string()
  # Note here incorrect/trimmed chain id
  # There's no way to correctly output chain ids longer than 2 char in PDB format
  print(o_pdb_str)
  assert o_pdb_str==""
#   assert_lines_in_text(o_pdb_str, """\
# ATOM     61  C   SERA-2   2      72.898  71.361  62.393  1.00 67.20           C
# ATOM     62  O   SERA-2   2      73.055  70.333  61.737  1.00 65.10           O
#     """)

  o_cif_str = "%s" % h.as_cif_block()
  print(o_cif_str)
  assert_lines_in_text(o_cif_str, """\
  ATOM  116  HG23  .  THR  A-2  5  ?   78.43100   68.11100   59.15800  1.000   48.55000  H  ?  B  ?  4  HG23  1
  ATOM  117  N     .  SER  A-5  2  ?   62.97800  157.54800  157.30800  1.000   68.85000  N  ?  C  ?  1  N     1
    """)

def test2():
  """ Overall counts duplicate atom labels should be 0 if long chain ids
  are processing correctly
  """
  inp = iotbx.pdb.input(lines=mmcif_str.split("\n"), source_info=None)
  h = inp.construct_hierarchy()
  oc = h.overall_counts()
  oc.show()
  assert oc.errors() == []

def test3():
  """ Construct restraints and output .geo file.
  Make sure long chain ids are properly outputted.
  """
  inp = iotbx.pdb.input(lines=mmcif_str.split("\n"), source_info=None)
  m = mmtbx.model.manager(model_input = inp)
  m.process(make_restraints=True)
  geo = m.restraints_as_geo()
  # print (geo)
  for l in [
    'bond pdb=" CA  PHEA-5   4 "',
    'pdb=" CB  PHEA-5   4 "',
    'angle pdb=" C   THRA-5   5 "',
    'dihedral pdb=" N   SERA-5   2 "',
    'chirality pdb=" CA  ASNA-5   3 "',
    'nonbonded pdb=" O   PHEA-5   4 "',]:
    assert_lines_in_text(geo, l)

def test4():
  """
  Test selections
  """
  inp = iotbx.pdb.input(lines=mmcif_str.split("\n"), source_info=None)
  h = inp.construct_hierarchy()
  asc = h.atom_selection_cache()
  s = asc.iselection("chain A")
  assert list(s) == list(range(0,58)), list(s)
  s = asc.iselection("chain A-2")
  assert list(s) == list(range(58,116)), list(s)
  s = asc.iselection("chain A-5")
  assert list(s) == list(range(116,174)), list(s)


if (__name__ == "__main__"):
  t0 = time.time()
  test1()
  test2()
  mon_lib_srv = None
  try:
    mon_lib_srv = monomer_library.server.server()
  except: # intentional
    print("Can not initialize monomer_library, skipping test.")
  if mon_lib_srv is not None:
    test3()
  test4()
  print("OK. Time: %8.3f"%(time.time()-t0))
