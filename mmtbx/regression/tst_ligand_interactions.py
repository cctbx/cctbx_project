from __future__ import absolute_import, division, print_function
'''
Tests for mmtbx.validation.ligand_interactions.

Model: ligand EDO A 1 and other EDO copies placed around it (ideal geostd
geometry, H at electron-cloud X-H lengths after processing):
  EDO A 2   O1 accepts an H-bond from the ligand's O1-HO1 (H...A 2.01 A, 180 deg)
  EDO A 3   H12 clashes with the ligand's H11 (1.75 A, overlap 0.69 A)
  EDO A 4   C1 and H12 in van der Waals contact with the ligand's O2 and HO2
  EDO A 0   15 A away, first in the file: outside validate_ligands' 3 A region,
            so its records' i_seqs must be mapped to the full model
The other copies are environment: probe2 is run by the ligand's selection, not
its resname.
'''
import json
from six.moves import cStringIO as StringIO
import iotbx.cif
import iotbx.pdb
import mmtbx.model
from libtbx import easy_run
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out, Sorry
from scitbx.array_family import flex
from cctbx import sgtbx
import cctbx.geometry_restraints.process_nonbonded_proxies as pnp
from mmtbx.validation import ligand_interactions as LI

LIG_SEL = 'chain A and resseq 1 and resname EDO'

model_str = '''
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
HETATM    1  C1  EDO A   0      30.177  10.059  15.375  1.00 20.00           C
HETATM    2  O1  EDO A   0      31.075  10.413  14.347  1.00 20.00           O
HETATM    3  C2  EDO A   0      28.776  10.464  14.992  1.00 20.00           C
HETATM    4  O2  EDO A   0      28.340   9.696  13.892  1.00 20.00           O
HETATM    5  H11 EDO A   0      30.198   8.980  15.580  1.00 20.00           H
HETATM    6  H12 EDO A   0      30.418  10.564  16.319  1.00 20.00           H
HETATM    7  HO1 EDO A   0      31.953  10.118  14.602  1.00 20.00           H
HETATM    8  H21 EDO A   0      28.139  10.312  15.873  1.00 20.00           H
HETATM    9  H22 EDO A   0      28.758  11.539  14.771  1.00 20.00           H
HETATM   10  HO2 EDO A   0      27.463  10.001  13.644  1.00 20.00           H
HETATM    1  C1  EDO A   1      15.177  10.059  15.375  1.00 20.00           C
HETATM    2  O1  EDO A   1      16.075  10.413  14.347  1.00 20.00           O
HETATM    3  C2  EDO A   1      13.776  10.464  14.992  1.00 20.00           C
HETATM    4  O2  EDO A   1      13.340   9.696  13.892  1.00 20.00           O
HETATM    5  H11 EDO A   1      15.198   8.980  15.580  1.00 20.00           H
HETATM    6  H12 EDO A   1      15.418  10.564  16.319  1.00 20.00           H
HETATM    7  HO1 EDO A   1      16.953  10.118  14.602  1.00 20.00           H
HETATM    8  H21 EDO A   1      13.139  10.312  15.873  1.00 20.00           H
HETATM    9  H22 EDO A   1      13.758  11.539  14.771  1.00 20.00           H
HETATM   10  HO2 EDO A   1      12.463  10.001  13.644  1.00 20.00           H
HETATM   11  C1  EDO A   2      19.869   9.698  14.350  1.00 20.00           C
HETATM   12  O1  EDO A   2      18.690   9.535  15.106  1.00 20.00           O
HETATM   13  C2  EDO A   2      19.563  10.454  13.082  1.00 20.00           C
HETATM   14  O2  EDO A   2      19.185  11.776  13.394  1.00 20.00           O
HETATM   15  H11 EDO A   2      20.643  10.236  14.914  1.00 20.00           H
HETATM   16  H12 EDO A   2      20.303   8.733  14.059  1.00 20.00           H
HETATM   17  HO1 EDO A   2      18.917   9.078  15.920  1.00 20.00           H
HETATM   18  H21 EDO A   2      20.466  10.425  12.458  1.00 20.00           H
HETATM   19  H22 EDO A   2      18.778   9.926  12.526  1.00 20.00           H
HETATM   20  HO2 EDO A   2      18.957  12.223  12.574  1.00 20.00           H
HETATM   21  C1  EDO A   3      15.673   6.534  16.101  1.00 20.00           C
HETATM   22  O1  EDO A   3      16.471   6.047  15.046  1.00 20.00           O
HETATM   23  C2  EDO A   3      14.542   5.574  16.373  1.00 20.00           C
HETATM   24  O2  EDO A   3      15.055   4.361  16.877  1.00 20.00           O
HETATM   25  H11 EDO A   3      16.256   6.673  17.022  1.00 20.00           H
HETATM   26  H12 EDO A   3      15.226   7.507  15.860  1.00 20.00           H
HETATM   27  HO1 EDO A   3      17.199   6.660  14.913  1.00 20.00           H
HETATM   28  H21 EDO A   3      13.865   6.060  17.088  1.00 20.00           H
HETATM   29  H22 EDO A   3      13.968   5.421  15.449  1.00 20.00           H
HETATM   30  HO2 EDO A   3      14.320   3.755  17.002  1.00 20.00           H
HETATM   31  C1  EDO A   4      10.915   8.639  11.757  1.00 20.00           C
HETATM   32  O1  EDO A   4      11.786   9.042  10.725  1.00 20.00           O
HETATM   33  C2  EDO A   4       9.559   8.312  11.184  1.00 20.00           C
HETATM   34  O2  EDO A   4       9.650   7.170  10.362  1.00 20.00           O
HETATM   35  H11 EDO A   4      11.299   7.763  12.296  1.00 20.00           H
HETATM   36  H12 EDO A   4      10.772   9.431  12.504  1.00 20.00           H
HETATM   37  HO1 EDO A   4      12.649   9.215  11.112  1.00 20.00           H
HETATM   38  H21 EDO A   4       8.876   8.153  12.028  1.00 20.00           H
HETATM   39  H22 EDO A   4       9.181   9.182  10.631  1.00 20.00           H
HETATM   40  HO2 EDO A   4       8.785   7.008   9.978  1.00 20.00           H
END
'''


# The ligand alone in a 6 A cell along x: its x+1 and x-1 copies make O-H...O
# H-bonds (pnp, symmetry); EDO A 2 is 15 A away along y as probe2's target.
sym_model_str = '''
CRYST1    6.000   40.000   40.000  90.00  90.00  90.00 P 1
HETATM    1  C1  EDO A   1      15.177  10.059  15.375  1.00 20.00           C
HETATM    2  O1  EDO A   1      16.075  10.413  14.347  1.00 20.00           O
HETATM    3  C2  EDO A   1      13.776  10.464  14.992  1.00 20.00           C
HETATM    4  O2  EDO A   1      13.340   9.696  13.892  1.00 20.00           O
HETATM    5  H11 EDO A   1      15.198   8.980  15.580  1.00 20.00           H
HETATM    6  H12 EDO A   1      15.418  10.564  16.319  1.00 20.00           H
HETATM    7  HO1 EDO A   1      16.953  10.118  14.602  1.00 20.00           H
HETATM    8  H21 EDO A   1      13.139  10.312  15.873  1.00 20.00           H
HETATM    9  H22 EDO A   1      13.758  11.539  14.771  1.00 20.00           H
HETATM   10  HO2 EDO A   1      12.463  10.001  13.644  1.00 20.00           H
HETATM   11  C1  EDO A   2      15.177  25.059  15.375  1.00 20.00           C
HETATM   12  O1  EDO A   2      16.075  25.413  14.347  1.00 20.00           O
HETATM   13  C2  EDO A   2      13.776  25.464  14.992  1.00 20.00           C
HETATM   14  O2  EDO A   2      13.340  24.696  13.892  1.00 20.00           O
HETATM   15  H11 EDO A   2      15.198  23.980  15.580  1.00 20.00           H
HETATM   16  H12 EDO A   2      15.418  25.564  16.319  1.00 20.00           H
HETATM   17  HO1 EDO A   2      16.953  25.118  14.602  1.00 20.00           H
HETATM   18  H21 EDO A   2      13.139  25.312  15.873  1.00 20.00           H
HETATM   19  H22 EDO A   2      13.758  26.539  14.771  1.00 20.00           H
HETATM   20  HO2 EDO A   2      12.463  25.001  13.644  1.00 20.00           H
END
'''

# EDO A 1 folded (O1-C1-C2-O2 52 deg, both hydroxyl H turned inwards): HO1...HO2
# 1.61 A, an intramolecular clash (pnp overlap 0.49 A); EDO A 2, a copy 15 A away
# along y, is probe2's target.
fold_model_str = '''
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
HETATM    1  C1  EDO A   1      15.177  10.059  15.375  1.00 20.00           C
HETATM    2  O1  EDO A   1      16.075  10.413  14.347  1.00 20.00           O
HETATM    3  C2  EDO A   1      13.776  10.464  14.992  1.00 20.00           C
HETATM    4  O2  EDO A   1      13.467   9.962  13.710  1.00 20.00           O
HETATM    5  H11 EDO A   1      15.198   8.980  15.580  1.00 20.00           H
HETATM    6  H12 EDO A   1      15.418  10.564  16.319  1.00 20.00           H
HETATM    7  HO1 EDO A   1      15.708  10.111  13.512  1.00 20.00           H
HETATM    8  H21 EDO A   1      13.099  10.066  15.759  1.00 20.00           H
HETATM    9  H22 EDO A   1      13.692  11.558  15.033  1.00 20.00           H
HETATM   10  HO2 EDO A   1      14.288   9.874  13.218  1.00 20.00           H
HETATM   11  C1  EDO A   2      15.177  25.059  15.375  1.00 20.00           C
HETATM   12  O1  EDO A   2      16.075  25.413  14.347  1.00 20.00           O
HETATM   13  C2  EDO A   2      13.776  25.464  14.992  1.00 20.00           C
HETATM   14  O2  EDO A   2      13.467  24.962  13.710  1.00 20.00           O
HETATM   15  H11 EDO A   2      15.198  23.980  15.580  1.00 20.00           H
HETATM   16  H12 EDO A   2      15.418  25.564  16.319  1.00 20.00           H
HETATM   17  HO1 EDO A   2      15.708  25.111  13.512  1.00 20.00           H
HETATM   18  H21 EDO A   2      13.099  25.066  15.759  1.00 20.00           H
HETATM   19  H22 EDO A   2      13.692  26.558  15.033  1.00 20.00           H
HETATM   20  HO2 EDO A   2      14.288  24.874  13.218  1.00 20.00           H
END
'''

# Inline clash (modelled on 5XH3's C12-H10...His208 NE2): EDO A 1 C1-H11 points
# at the NE2 of a neutral His (HD1 only), H11...NE2 1.78 A, C1...NE2 2.63 A.
inline_model_str = '''
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
HETATM    1  C1  EDO A   1      15.177  10.059  15.375  1.00 20.00           C
HETATM    2  O1  EDO A   1      16.075  10.413  14.347  1.00 20.00           O
HETATM    3  C2  EDO A   1      13.776  10.464  14.992  1.00 20.00           C
HETATM    4  O2  EDO A   1      13.340   9.696  13.892  1.00 20.00           O
HETATM    5  H11 EDO A   1      15.198   8.980  15.580  1.00 20.00           H
HETATM    6  H12 EDO A   1      15.418  10.564  16.319  1.00 20.00           H
HETATM    7  HO1 EDO A   1      16.953  10.118  14.602  1.00 20.00           H
HETATM    8  H21 EDO A   1      13.139  10.312  15.873  1.00 20.00           H
HETATM    9  H22 EDO A   1      13.758  11.539  14.771  1.00 20.00           H
HETATM   10  HO2 EDO A   1      12.463  10.001  13.644  1.00 20.00           H
ATOM     11  N   HIS B  10      16.597   5.849  21.926  1.00 20.00           N
ATOM     12  CA  HIS B  10      16.646   5.791  20.487  1.00 20.00           C
ATOM     13  C   HIS B  10      17.472   4.603  20.014  1.00 20.00           C
ATOM     14  O   HIS B  10      18.460   4.195  20.615  1.00 20.00           O
ATOM     15  CB  HIS B  10      17.213   7.089  19.900  1.00 20.00           C
ATOM     16  CG  HIS B  10      16.818   7.307  18.459  1.00 20.00           C
ATOM     17  ND1 HIS B  10      17.593   6.832  17.459  1.00 20.00           N
ATOM     18  CD2 HIS B  10      15.753   7.933  17.946  1.00 20.00           C
ATOM     19  CE1 HIS B  10      17.014   7.160  16.300  1.00 20.00           C
ATOM     20  NE2 HIS B  10      15.892   7.831  16.583  1.00 20.00           N
ATOM     21  H   HIS B  10      15.833   5.404  22.407  1.00 20.00           H
ATOM     22  HA  HIS B  10      15.618   5.617  20.152  1.00 20.00           H
ATOM     23  HB2 HIS B  10      16.869   7.955  20.482  1.00 20.00           H
ATOM     24  HB3 HIS B  10      18.308   7.099  19.989  1.00 20.00           H
ATOM     25  HD1 HIS B  10      18.463   6.315  17.552  1.00 20.00           H
ATOM     26  HD2 HIS B  10      14.890   8.449  18.319  1.00 20.00           H
ATOM     27  HE1 HIS B  10      17.384   6.926  15.315  1.00 20.00           H
END
'''

# EDO A 1 O1-HO1...EDO A 2 O1 at H...O 2.66 A (electron-cloud H): a pnp H-bond
# without overlap in probe2 (cc by default, wh with weak H-bonds).
weak_model_str = '''
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
HETATM    1  C1  EDO A   1      15.177  10.059  15.375  1.00 20.00           C
HETATM    2  O1  EDO A   1      16.075  10.413  14.347  1.00 20.00           O
HETATM    3  C2  EDO A   1      13.776  10.464  14.992  1.00 20.00           C
HETATM    4  O2  EDO A   1      13.340   9.696  13.892  1.00 20.00           O
HETATM    5  H11 EDO A   1      15.198   8.980  15.580  1.00 20.00           H
HETATM    6  H12 EDO A   1      15.418  10.564  16.319  1.00 20.00           H
HETATM    7  HO1 EDO A   1      16.953  10.118  14.602  1.00 20.00           H
HETATM    8  H21 EDO A   1      13.139  10.312  15.873  1.00 20.00           H
HETATM    9  H22 EDO A   1      13.758  11.539  14.771  1.00 20.00           H
HETATM   10  HO2 EDO A   1      12.463  10.001  13.644  1.00 20.00           H
HETATM   11  C1  EDO A   2      20.463   9.498  14.523  1.00 20.00           C
HETATM   12  O1  EDO A   2      19.284   9.335  15.278  1.00 20.00           O
HETATM   13  C2  EDO A   2      20.158  10.255  13.254  1.00 20.00           C
HETATM   14  O2  EDO A   2      19.779  11.577  13.566  1.00 20.00           O
HETATM   15  H11 EDO A   2      21.237  10.036  15.086  1.00 20.00           H
HETATM   16  H12 EDO A   2      20.897   8.533  14.231  1.00 20.00           H
HETATM   17  HO1 EDO A   2      19.511   8.878  16.092  1.00 20.00           H
HETATM   18  H21 EDO A   2      21.060  10.226  12.631  1.00 20.00           H
HETATM   19  H22 EDO A   2      19.372   9.726  12.698  1.00 20.00           H
HETATM   20  HO2 EDO A   2      19.551  12.023  12.747  1.00 20.00           H
END
'''

# Salt bridges (CCD/GeoStd ideal geometry, partners placed by a clash-free
# orientation search): ACT A 1 with Lys B 10 (NZ...OXT 2.85 A, HZ1 towards OXT)
# and Arg C 20 (NH1...O 2.90 A); ACT A 2 with His D 30 (HD1 and HE2; NE2...OXT
# 2.80 A), Asp E 40 (OD1...O 3.20 A, same charge) and Lys F 50 (NZ...O 5.00 A);
# NH4 A 3 with Asp G 60 (OD1...N 2.85 A).
salt_model_str = '''
CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1
HETATM    1  C   ACT A   1     -40.072  16.000 -30.000  1.00 20.00           C
HETATM    2  O   ACT A   1     -40.682  17.056 -30.000  1.00 20.00           O
HETATM    3  OXT ACT A   1     -40.682  14.944 -30.000  1.00 20.00           O
HETATM    4  CH3 ACT A   1     -38.565  16.000 -30.000  1.00 20.00           C
HETATM    5  H1  ACT A   1     -38.201  16.000 -28.972  1.00 20.00           H
HETATM    6  H2  ACT A   1     -38.201  15.110 -30.514  1.00 20.00           H
HETATM    7  H3  ACT A   1     -38.201  16.890 -30.514  1.00 20.00           H
HETATM    8  C   ACT A   2     -10.072  16.000 -30.000  1.00 20.00           C
HETATM    9  O   ACT A   2     -10.682  17.056 -30.000  1.00 20.00           O
HETATM   10  OXT ACT A   2     -10.682  14.944 -30.000  1.00 20.00           O
HETATM   11  CH3 ACT A   2      -8.565  16.000 -30.000  1.00 20.00           C
HETATM   12  H1  ACT A   2      -8.201  16.000 -28.972  1.00 20.00           H
HETATM   13  H2  ACT A   2      -8.201  15.110 -30.514  1.00 20.00           H
HETATM   14  H3  ACT A   2      -8.201  16.890 -30.514  1.00 20.00           H
HETATM   15  N   NH4 A   3      -9.982  30.026   0.000  1.00 20.00           N
HETATM   16  HN1 NH4 A   3     -10.367  29.481  -0.771  1.00 20.00           H
HETATM   17  HN2 NH4 A   3      -8.962  30.026   0.000  1.00 20.00           H
HETATM   18  HN3 NH4 A   3     -10.322  30.988   0.000  1.00 20.00           H
HETATM   19  HN4 NH4 A   3     -10.367  29.481   0.771  1.00 20.00           H
ATOM     20  N   LYS B  10     -38.213  10.732 -24.880  1.00 20.00           N
ATOM     21  CA  LYS B  10     -38.443   9.616 -25.807  1.00 20.00           C
ATOM     22  C   LYS B  10     -37.177   8.811 -25.947  1.00 20.00           C
ATOM     23  O   LYS B  10     -36.113   9.301 -25.652  1.00 20.00           O
ATOM     24  CB  LYS B  10     -38.851  10.165 -27.175  1.00 20.00           C
ATOM     25  CG  LYS B  10     -40.201  10.877 -27.056  1.00 20.00           C
ATOM     26  CD  LYS B  10     -40.610  11.427 -28.425  1.00 20.00           C
ATOM     27  CE  LYS B  10     -41.958  12.139 -28.306  1.00 20.00           C
ATOM     28  NZ  LYS B  10     -42.352  12.666 -29.619  1.00 20.00           N
ATOM     29  H   LYS B  10     -37.474  11.292 -25.278  1.00 20.00           H
ATOM     30  HA  LYS B  10     -39.237   8.979 -25.419  1.00 20.00           H
ATOM     31  HB2 LYS B  10     -38.098  10.871 -27.524  1.00 20.00           H
ATOM     32  HB3 LYS B  10     -38.936   9.344 -27.887  1.00 20.00           H
ATOM     33  HG2 LYS B  10     -40.954  10.171 -26.708  1.00 20.00           H
ATOM     34  HG3 LYS B  10     -40.117  11.700 -26.345  1.00 20.00           H
ATOM     35  HD2 LYS B  10     -39.856  12.133 -28.773  1.00 20.00           H
ATOM     36  HD3 LYS B  10     -40.694  10.605 -29.135  1.00 20.00           H
ATOM     37  HE2 LYS B  10     -42.713  11.433 -27.958  1.00 20.00           H
ATOM     38  HE3 LYS B  10     -41.875  12.960 -27.594  1.00 20.00           H
ATOM     39  HZ1 LYS B  10     -41.653  13.320 -29.941  1.00 20.00           H
ATOM     40  HZ2 LYS B  10     -42.429  11.906 -30.277  1.00 20.00           H
ATOM     41  HZ3 LYS B  10     -43.241  13.136 -29.541  1.00 20.00           H
ATOM     42  N   ARG C  20     -40.754  18.537 -36.245  1.00 20.00           N
ATOM     43  CA  ARG C  20     -41.671  17.454 -35.896  1.00 20.00           C
ATOM     44  C   ARG C  20     -42.705  17.336 -37.001  1.00 20.00           C
ATOM     45  O   ARG C  20     -42.793  18.075 -37.973  1.00 20.00           O
ATOM     46  CB  ARG C  20     -42.321  17.686 -34.523  1.00 20.00           C
ATOM     47  CG  ARG C  20     -43.220  18.931 -34.444  1.00 20.00           C
ATOM     48  CD  ARG C  20     -43.886  19.083 -33.079  1.00 20.00           C
ATOM     49  NE  ARG C  20     -42.899  19.207 -32.032  1.00 20.00           N
ATOM     50  CZ  ARG C  20     -43.236  19.361 -30.676  1.00 20.00           C
ATOM     51  NH1 ARG C  20     -42.251  19.477 -29.702  1.00 20.00           N
ATOM     52  NH2 ARG C  20     -44.574  19.397 -30.296  1.00 20.00           N
ATOM     53  H   ARG C  20     -39.937  18.656 -35.686  1.00 20.00           H
ATOM     54  HA  ARG C  20     -41.092  16.524 -35.875  1.00 20.00           H
ATOM     55  HB2 ARG C  20     -41.532  17.765 -33.765  1.00 20.00           H
ATOM     56  HB3 ARG C  20     -42.916  16.802 -34.259  1.00 20.00           H
ATOM     57  HG2 ARG C  20     -44.010  18.862 -35.201  1.00 20.00           H
ATOM     58  HG3 ARG C  20     -42.635  19.831 -34.669  1.00 20.00           H
ATOM     59  HD2 ARG C  20     -44.512  19.982 -33.076  1.00 20.00           H
ATOM     60  HD3 ARG C  20     -44.537  18.228 -32.870  1.00 20.00           H
ATOM     61  HE  ARG C  20     -41.897  19.190 -32.257  1.00 20.00           H
ATOM     62 HH11 ARG C  20     -42.497  19.587 -28.720  1.00 20.00           H
ATOM     63 HH12 ARG C  20     -41.261  19.458 -29.933  1.00 20.00           H
ATOM     64 HH21 ARG C  20     -45.319  19.314 -30.983  1.00 20.00           H
ATOM     65 HH22 ARG C  20     -44.841  19.506 -29.320  1.00 20.00           H
ATOM     66  N   HIS D  30     -15.596   8.339 -28.659  1.00 20.00           N
ATOM     67  CA  HIS D  30     -15.301   9.562 -29.361  1.00 20.00           C
ATOM     68  C   HIS D  30     -16.410   9.913 -30.343  1.00 20.00           C
ATOM     69  O   HIS D  30     -17.050   9.066 -30.957  1.00 20.00           O
ATOM     70  CB  HIS D  30     -13.960   9.466 -30.100  1.00 20.00           C
ATOM     71  CG  HIS D  30     -13.354  10.816 -30.401  1.00 20.00           C
ATOM     72  ND1 HIS D  30     -13.642  11.451 -31.559  1.00 20.00           N
ATOM     73  CD2 HIS D  30     -12.519  11.564 -29.671  1.00 20.00           C
ATOM     74  CE1 HIS D  30     -12.981  12.613 -31.565  1.00 20.00           C
ATOM     75  NE2 HIS D  30     -12.294  12.694 -30.421  1.00 20.00           N
ATOM     76  H   HIS D  30     -16.068   8.388 -27.771  1.00 20.00           H
ATOM     77  HA  HIS D  30     -15.284  10.355 -28.604  1.00 20.00           H
ATOM     78  HB2 HIS D  30     -13.236   8.893 -29.505  1.00 20.00           H
ATOM     79  HB3 HIS D  30     -14.080   8.899 -31.032  1.00 20.00           H
ATOM     80  HD1 HIS D  30     -14.249  11.120 -32.304  1.00 20.00           H
ATOM     81  HD2 HIS D  30     -12.028  11.495 -28.720  1.00 20.00           H
ATOM     82  HE1 HIS D  30     -13.000  13.351 -32.350  1.00 20.00           H
ATOM     83  HE2 HIS D  30     -11.700  13.472 -30.150  1.00 20.00           H
ATOM     84  N   ASP E  40     -12.256  16.952 -25.815  1.00 20.00           N
ATOM     85  CA  ASP E  40     -11.206  17.762 -26.447  1.00 20.00           C
ATOM     86  C   ASP E  40     -10.025  17.868 -25.518  1.00 20.00           C
ATOM     87  O   ASP E  40     -10.165  17.656 -24.336  1.00 20.00           O
ATOM     88  CB  ASP E  40     -11.750  19.161 -26.743  1.00 20.00           C
ATOM     89  CG  ASP E  40     -12.852  19.066 -27.767  1.00 20.00           C
ATOM     90  OD1 ASP E  40     -13.168  17.990 -28.215  1.00 20.00           O
ATOM     91  OD2 ASP E  40     -13.480  20.177 -28.182  1.00 20.00           O
ATOM     92  H   ASP E  40     -11.942  16.004 -25.669  1.00 20.00           H
ATOM     93  HA  ASP E  40     -10.893  17.289 -27.379  1.00 20.00           H
ATOM     94  HB2 ASP E  40     -12.144  19.600 -25.826  1.00 20.00           H
ATOM     95  HB3 ASP E  40     -10.947  19.789 -27.130  1.00 20.00           H
ATOM     96  N   LYS F  50     -21.140  20.358 -30.638  1.00 20.00           N
ATOM     97  CA  LYS F  50     -21.346  19.059 -31.293  1.00 20.00           C
ATOM     98  C   LYS F  50     -22.586  18.408 -30.737  1.00 20.00           C
ATOM     99  O   LYS F  50     -23.026  18.759 -29.669  1.00 20.00           O
ATOM    100  CB  LYS F  50     -20.137  18.159 -31.030  1.00 20.00           C
ATOM    101  CG  LYS F  50     -18.900  18.760 -31.702  1.00 20.00           C
ATOM    102  CD  LYS F  50     -17.689  17.860 -31.440  1.00 20.00           C
ATOM    103  CE  LYS F  50     -16.453  18.460 -32.111  1.00 20.00           C
ATOM    104  NZ  LYS F  50     -15.291  17.596 -31.861  1.00 20.00           N
ATOM    105  H   LYS F  50     -21.037  20.172 -29.652  1.00 20.00           H
ATOM    106  HA  LYS F  50     -21.463  19.208 -32.365  1.00 20.00           H
ATOM    107  HB2 LYS F  50     -19.966  18.082 -29.957  1.00 20.00           H
ATOM    108  HB3 LYS F  50     -20.326  17.167 -31.440  1.00 20.00           H
ATOM    109  HG2 LYS F  50     -19.070  18.836 -32.775  1.00 20.00           H
ATOM    110  HG3 LYS F  50     -18.709  19.752 -31.293  1.00 20.00           H
ATOM    111  HD2 LYS F  50     -17.519  17.783 -30.366  1.00 20.00           H
ATOM    112  HD3 LYS F  50     -17.880  16.868 -31.849  1.00 20.00           H
ATOM    113  HE2 LYS F  50     -16.623  18.537 -33.185  1.00 20.00           H
ATOM    114  HE3 LYS F  50     -16.263  19.452 -31.702  1.00 20.00           H
ATOM    115  HZ1 LYS F  50     -15.134  17.526 -30.866  1.00 20.00           H
ATOM    116  HZ2 LYS F  50     -15.468  16.678 -32.239  1.00 20.00           H
ATOM    117  HZ3 LYS F  50     -14.475  17.992 -32.303  1.00 20.00           H
ATOM    118  N   ASP G  60      -7.514  34.467   2.148  1.00 20.00           N
ATOM    119  CA  ASP G  60      -8.545  34.210   3.162  1.00 20.00           C
ATOM    120  C   ASP G  60      -8.453  35.251   4.247  1.00 20.00           C
ATOM    121  O   ASP G  60      -7.433  35.883   4.393  1.00 20.00           O
ATOM    122  CB  ASP G  60      -8.331  32.821   3.767  1.00 20.00           C
ATOM    123  CG  ASP G  60      -8.544  31.772   2.706  1.00 20.00           C
ATOM    124  OD1 ASP G  60      -8.838  32.102   1.582  1.00 20.00           O
ATOM    125  OD2 ASP G  60      -8.408  30.473   3.011  1.00 20.00           O
ATOM    126  H   ASP G  60      -7.672  35.351   1.688  1.00 20.00           H
ATOM    127  HA  ASP G  60      -9.531  34.256   2.698  1.00 20.00           H
ATOM    128  HB2 ASP G  60      -7.314  32.746   4.153  1.00 20.00           H
ATOM    129  HB3 ASP G  60      -9.040  32.666   4.580  1.00 20.00           H
END
'''

# The neutral forms at the same places: acetic acid (ACY, custom restraints
# below, HXT on OXT turned away from Lys) and His D 30 with HD1 only.
neutral_model_str = '''
CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1
HETATM    1  C   ACY A   1     -40.072  16.000 -30.000  1.00 20.00           C
HETATM    2  O   ACY A   1     -40.682  17.056 -30.000  1.00 20.00           O
HETATM    3  OXT ACY A   1     -40.682  14.944 -30.000  1.00 20.00           O
HETATM    4  CH3 ACY A   1     -38.565  16.000 -30.000  1.00 20.00           C
HETATM    5  HXT ACY A   1     -40.627  14.301 -30.783  1.00 20.00           H
HETATM    6  H1  ACY A   1     -38.201  16.000 -28.972  1.00 20.00           H
HETATM    7  H2  ACY A   1     -38.201  15.110 -30.514  1.00 20.00           H
HETATM    8  H3  ACY A   1     -38.201  16.890 -30.514  1.00 20.00           H
HETATM    9  C   ACT A   2     -10.072  16.000 -30.000  1.00 20.00           C
HETATM   10  O   ACT A   2     -10.682  17.056 -30.000  1.00 20.00           O
HETATM   11  OXT ACT A   2     -10.682  14.944 -30.000  1.00 20.00           O
HETATM   12  CH3 ACT A   2      -8.565  16.000 -30.000  1.00 20.00           C
HETATM   13  H1  ACT A   2      -8.201  16.000 -28.972  1.00 20.00           H
HETATM   14  H2  ACT A   2      -8.201  15.110 -30.514  1.00 20.00           H
HETATM   15  H3  ACT A   2      -8.201  16.890 -30.514  1.00 20.00           H
ATOM     16  N   LYS B  10     -38.213  10.732 -24.880  1.00 20.00           N
ATOM     17  CA  LYS B  10     -38.443   9.616 -25.807  1.00 20.00           C
ATOM     18  C   LYS B  10     -37.177   8.811 -25.947  1.00 20.00           C
ATOM     19  O   LYS B  10     -36.113   9.301 -25.652  1.00 20.00           O
ATOM     20  CB  LYS B  10     -38.851  10.165 -27.175  1.00 20.00           C
ATOM     21  CG  LYS B  10     -40.201  10.877 -27.056  1.00 20.00           C
ATOM     22  CD  LYS B  10     -40.610  11.427 -28.425  1.00 20.00           C
ATOM     23  CE  LYS B  10     -41.958  12.139 -28.306  1.00 20.00           C
ATOM     24  NZ  LYS B  10     -42.352  12.666 -29.619  1.00 20.00           N
ATOM     25  H   LYS B  10     -37.474  11.292 -25.278  1.00 20.00           H
ATOM     26  HA  LYS B  10     -39.237   8.979 -25.419  1.00 20.00           H
ATOM     27  HB2 LYS B  10     -38.098  10.871 -27.524  1.00 20.00           H
ATOM     28  HB3 LYS B  10     -38.936   9.344 -27.887  1.00 20.00           H
ATOM     29  HG2 LYS B  10     -40.954  10.171 -26.708  1.00 20.00           H
ATOM     30  HG3 LYS B  10     -40.117  11.700 -26.345  1.00 20.00           H
ATOM     31  HD2 LYS B  10     -39.856  12.133 -28.773  1.00 20.00           H
ATOM     32  HD3 LYS B  10     -40.694  10.605 -29.135  1.00 20.00           H
ATOM     33  HE2 LYS B  10     -42.713  11.433 -27.958  1.00 20.00           H
ATOM     34  HE3 LYS B  10     -41.875  12.960 -27.594  1.00 20.00           H
ATOM     35  HZ1 LYS B  10     -41.653  13.320 -29.941  1.00 20.00           H
ATOM     36  HZ2 LYS B  10     -42.429  11.906 -30.277  1.00 20.00           H
ATOM     37  HZ3 LYS B  10     -43.241  13.136 -29.541  1.00 20.00           H
ATOM     38  N   HIS D  30     -15.596   8.339 -28.659  1.00 20.00           N
ATOM     39  CA  HIS D  30     -15.301   9.562 -29.361  1.00 20.00           C
ATOM     40  C   HIS D  30     -16.410   9.913 -30.343  1.00 20.00           C
ATOM     41  O   HIS D  30     -17.050   9.066 -30.957  1.00 20.00           O
ATOM     42  CB  HIS D  30     -13.960   9.466 -30.100  1.00 20.00           C
ATOM     43  CG  HIS D  30     -13.354  10.816 -30.401  1.00 20.00           C
ATOM     44  ND1 HIS D  30     -13.642  11.451 -31.559  1.00 20.00           N
ATOM     45  CD2 HIS D  30     -12.519  11.564 -29.671  1.00 20.00           C
ATOM     46  CE1 HIS D  30     -12.981  12.613 -31.565  1.00 20.00           C
ATOM     47  NE2 HIS D  30     -12.294  12.694 -30.421  1.00 20.00           N
ATOM     48  H   HIS D  30     -16.068   8.388 -27.771  1.00 20.00           H
ATOM     49  HA  HIS D  30     -15.284  10.355 -28.604  1.00 20.00           H
ATOM     50  HB2 HIS D  30     -13.236   8.893 -29.505  1.00 20.00           H
ATOM     51  HB3 HIS D  30     -14.080   8.899 -31.032  1.00 20.00           H
ATOM     52  HD1 HIS D  30     -14.249  11.120 -32.304  1.00 20.00           H
ATOM     53  HD2 HIS D  30     -12.028  11.495 -28.720  1.00 20.00           H
ATOM     54  HE1 HIS D  30     -13.000  13.351 -32.350  1.00 20.00           H
END
'''

# ACT A 1 and Lys B 10 from salt_model_str, Lys moved by -100 A along x (the cell
# edge): the salt bridge is to its x+1,y,z copy.
salt_sym_model_str = '''
CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1
HETATM    1  C   ACT A   1     -40.072  16.000 -30.000  1.00 20.00           C
HETATM    2  O   ACT A   1     -40.682  17.056 -30.000  1.00 20.00           O
HETATM    3  OXT ACT A   1     -40.682  14.944 -30.000  1.00 20.00           O
HETATM    4  CH3 ACT A   1     -38.565  16.000 -30.000  1.00 20.00           C
HETATM    5  H1  ACT A   1     -38.201  16.000 -28.972  1.00 20.00           H
HETATM    6  H2  ACT A   1     -38.201  15.110 -30.514  1.00 20.00           H
HETATM    7  H3  ACT A   1     -38.201  16.890 -30.514  1.00 20.00           H
ATOM     20  N   LYS B  10    -138.213  10.732 -24.880  1.00 20.00           N
ATOM     21  CA  LYS B  10    -138.443   9.616 -25.807  1.00 20.00           C
ATOM     22  C   LYS B  10    -137.177   8.811 -25.947  1.00 20.00           C
ATOM     23  O   LYS B  10    -136.113   9.301 -25.652  1.00 20.00           O
ATOM     24  CB  LYS B  10    -138.851  10.165 -27.175  1.00 20.00           C
ATOM     25  CG  LYS B  10    -140.201  10.877 -27.056  1.00 20.00           C
ATOM     26  CD  LYS B  10    -140.610  11.427 -28.425  1.00 20.00           C
ATOM     27  CE  LYS B  10    -141.958  12.139 -28.306  1.00 20.00           C
ATOM     28  NZ  LYS B  10    -142.352  12.666 -29.619  1.00 20.00           N
ATOM     29  H   LYS B  10    -137.474  11.292 -25.278  1.00 20.00           H
ATOM     30  HA  LYS B  10    -139.237   8.979 -25.419  1.00 20.00           H
ATOM     31  HB2 LYS B  10    -138.098  10.871 -27.524  1.00 20.00           H
ATOM     32  HB3 LYS B  10    -138.936   9.344 -27.887  1.00 20.00           H
ATOM     33  HG2 LYS B  10    -140.954  10.171 -26.708  1.00 20.00           H
ATOM     34  HG3 LYS B  10    -140.117  11.700 -26.345  1.00 20.00           H
ATOM     35  HD2 LYS B  10    -139.856  12.133 -28.773  1.00 20.00           H
ATOM     36  HD3 LYS B  10    -140.694  10.605 -29.135  1.00 20.00           H
ATOM     37  HE2 LYS B  10    -142.713  11.433 -27.958  1.00 20.00           H
ATOM     38  HE3 LYS B  10    -141.875  12.960 -27.594  1.00 20.00           H
ATOM     39  HZ1 LYS B  10    -141.653  13.320 -29.941  1.00 20.00           H
ATOM     40  HZ2 LYS B  10    -142.429  11.906 -30.277  1.00 20.00           H
ATOM     41  HZ3 LYS B  10    -143.241  13.136 -29.541  1.00 20.00           H
END
'''

# Standard residues by template (CCD ideal geometry, 20 A apart; no OXT): Asp B 10
# with HD2, His B 20 with HD1 only, His B 30 with HD1 and HE2, Lys B 40 with two
# H on NZ, Lys B 50 without any H, Arg B 60 complete.
templates_model_str = '''
CRYST1  200.000  200.000  200.000  90.00  90.00  90.00 P 1
ATOM      1  N   ASP B  10      -0.317   1.688   0.066  1.00 20.00           N
ATOM      2  CA  ASP B  10      -0.470   0.286  -0.344  1.00 20.00           C
ATOM      3  C   ASP B  10      -1.868  -0.180  -0.029  1.00 20.00           C
ATOM      4  O   ASP B  10      -2.534   0.415   0.786  1.00 20.00           O
ATOM      5  CB  ASP B  10       0.539  -0.580   0.413  1.00 20.00           C
ATOM      6  CG  ASP B  10       1.938  -0.195   0.004  1.00 20.00           C
ATOM      7  OD1 ASP B  10       2.109   0.681  -0.810  1.00 20.00           O
ATOM      8  OD2 ASP B  10       2.992  -0.826   0.543  1.00 20.00           O
ATOM      9  H   ASP B  10      -0.928   2.289  -0.467  1.00 20.00           H
ATOM     10  HA  ASP B  10      -0.292   0.199  -1.416  1.00 20.00           H
ATOM     11  HB2 ASP B  10       0.419  -0.425   1.485  1.00 20.00           H
ATOM     12  HB3 ASP B  10       0.367  -1.630   0.176  1.00 20.00           H
ATOM     13  HD2 ASP B  10       3.869  -0.545   0.250  1.00 20.00           H
ATOM     14  N   HIS B  20      19.960  -1.210   0.053  1.00 20.00           N
ATOM     15  CA  HIS B  20      21.172  -1.709   0.652  1.00 20.00           C
ATOM     16  C   HIS B  20      21.083  -3.207   0.905  1.00 20.00           C
ATOM     17  O   HIS B  20      20.040  -3.770   1.222  1.00 20.00           O
ATOM     18  CB  HIS B  20      21.484  -0.975   1.962  1.00 20.00           C
ATOM     19  CG  HIS B  20      22.940  -1.060   2.353  1.00 20.00           C
ATOM     20  ND1 HIS B  20      23.380  -2.075   3.129  1.00 20.00           N
ATOM     21  CD2 HIS B  20      23.960  -0.251   2.046  1.00 20.00           C
ATOM     22  CE1 HIS B  20      24.693  -1.908   3.317  1.00 20.00           C
ATOM     23  NE2 HIS B  20      25.058  -0.801   2.662  1.00 20.00           N
ATOM     24  H   HIS B  20      19.898  -1.155  -0.950  1.00 20.00           H
ATOM     25  HA  HIS B  20      21.965  -1.558  -0.089  1.00 20.00           H
ATOM     26  HB2 HIS B  20      21.215   0.087   1.879  1.00 20.00           H
ATOM     27  HB3 HIS B  20      20.859  -1.368   2.775  1.00 20.00           H
ATOM     28  HD1 HIS B  20      22.828  -2.838   3.511  1.00 20.00           H
ATOM     29  HD2 HIS B  20      24.108   0.647   1.479  1.00 20.00           H
ATOM     30  HE1 HIS B  20      25.340  -2.550   3.892  1.00 20.00           H
ATOM     31  N   HIS B  30      39.960  -1.210   0.053  1.00 20.00           N
ATOM     32  CA  HIS B  30      41.172  -1.709   0.652  1.00 20.00           C
ATOM     33  C   HIS B  30      41.083  -3.207   0.905  1.00 20.00           C
ATOM     34  O   HIS B  30      40.040  -3.770   1.222  1.00 20.00           O
ATOM     35  CB  HIS B  30      41.484  -0.975   1.962  1.00 20.00           C
ATOM     36  CG  HIS B  30      42.940  -1.060   2.353  1.00 20.00           C
ATOM     37  ND1 HIS B  30      43.380  -2.075   3.129  1.00 20.00           N
ATOM     38  CD2 HIS B  30      43.960  -0.251   2.046  1.00 20.00           C
ATOM     39  CE1 HIS B  30      44.693  -1.908   3.317  1.00 20.00           C
ATOM     40  NE2 HIS B  30      45.058  -0.801   2.662  1.00 20.00           N
ATOM     41  H   HIS B  30      39.898  -1.155  -0.950  1.00 20.00           H
ATOM     42  HA  HIS B  30      41.965  -1.558  -0.089  1.00 20.00           H
ATOM     43  HB2 HIS B  30      41.215   0.087   1.879  1.00 20.00           H
ATOM     44  HB3 HIS B  30      40.859  -1.368   2.775  1.00 20.00           H
ATOM     45  HD1 HIS B  30      42.828  -2.838   3.511  1.00 20.00           H
ATOM     46  HD2 HIS B  30      44.108   0.647   1.479  1.00 20.00           H
ATOM     47  HE1 HIS B  30      45.340  -2.550   3.892  1.00 20.00           H
ATOM     48  HE2 HIS B  30      46.002  -0.428   2.627  1.00 20.00           H
ATOM     49  N   LYS B  40      61.422   1.796   0.198  1.00 20.00           N
ATOM     50  CA  LYS B  40      61.394   0.355   0.484  1.00 20.00           C
ATOM     51  C   LYS B  40      62.657  -0.284  -0.032  1.00 20.00           C
ATOM     52  O   LYS B  40      63.316   0.275  -0.876  1.00 20.00           O
ATOM     53  CB  LYS B  40      60.184  -0.278  -0.206  1.00 20.00           C
ATOM     54  CG  LYS B  40      58.898   0.282   0.407  1.00 20.00           C
ATOM     55  CD  LYS B  40      57.687  -0.351  -0.283  1.00 20.00           C
ATOM     56  CE  LYS B  40      56.402   0.208   0.329  1.00 20.00           C
ATOM     57  NZ  LYS B  40      55.239  -0.400  -0.332  1.00 20.00           N
ATOM     58  H   LYS B  40      61.489   1.891  -0.804  1.00 20.00           H
ATOM     59  HA  LYS B  40      61.322   0.200   1.560  1.00 20.00           H
ATOM     60  HB2 LYS B  40      60.210  -0.047  -1.270  1.00 20.00           H
ATOM     61  HB3 LYS B  40      60.211  -1.359  -0.068  1.00 20.00           H
ATOM     62  HG2 LYS B  40      58.872   0.050   1.471  1.00 20.00           H
ATOM     63  HG3 LYS B  40      58.870   1.363   0.269  1.00 20.00           H
ATOM     64  HD2 LYS B  40      57.713  -0.120  -1.348  1.00 20.00           H
ATOM     65  HD3 LYS B  40      57.715  -1.432  -0.145  1.00 20.00           H
ATOM     66  HE2 LYS B  40      56.375  -0.023   1.394  1.00 20.00           H
ATOM     67  HE3 LYS B  40      56.374   1.289   0.192  1.00 20.00           H
ATOM     68  HZ1 LYS B  40      55.264  -0.185  -1.318  1.00 20.00           H
ATOM     69  HZ2 LYS B  40      55.265  -1.400  -0.205  1.00 20.00           H
ATOM     70  N   LYS B  50      81.422   1.796   0.198  1.00 20.00           N
ATOM     71  CA  LYS B  50      81.394   0.355   0.484  1.00 20.00           C
ATOM     72  C   LYS B  50      82.657  -0.284  -0.032  1.00 20.00           C
ATOM     73  O   LYS B  50      83.316   0.275  -0.876  1.00 20.00           O
ATOM     74  CB  LYS B  50      80.184  -0.278  -0.206  1.00 20.00           C
ATOM     75  CG  LYS B  50      78.898   0.282   0.407  1.00 20.00           C
ATOM     76  CD  LYS B  50      77.687  -0.351  -0.283  1.00 20.00           C
ATOM     77  CE  LYS B  50      76.402   0.208   0.329  1.00 20.00           C
ATOM     78  NZ  LYS B  50      75.239  -0.400  -0.332  1.00 20.00           N
ATOM     79  N   ARG B  60      99.531   1.110  -0.993  1.00 20.00           N
ATOM     80  CA  ARG B  60     100.004   2.294  -1.708  1.00 20.00           C
ATOM     81  C   ARG B  60      99.093   2.521  -2.901  1.00 20.00           C
ATOM     82  O   ARG B  60      98.173   1.789  -3.242  1.00 20.00           O
ATOM     83  CB  ARG B  60     101.475   2.150  -2.127  1.00 20.00           C
ATOM     84  CG  ARG B  60     101.745   1.017  -3.130  1.00 20.00           C
ATOM     85  CD  ARG B  60     103.210   0.954  -3.557  1.00 20.00           C
ATOM     86  NE  ARG B  60     104.071   0.726  -2.421  1.00 20.00           N
ATOM     87  CZ  ARG B  60     105.469   0.624  -2.528  1.00 20.00           C
ATOM     88  NH1 ARG B  60     106.259   0.404  -1.405  1.00 20.00           N
ATOM     89  NH2 ARG B  60     106.078   0.744  -3.773  1.00 20.00           N
ATOM     90  H   ARG B  60      99.942   0.903  -0.109  1.00 20.00           H
ATOM     91  HA  ARG B  60      99.897   3.152  -1.034  1.00 20.00           H
ATOM     92  HB2 ARG B  60     102.086   1.988  -1.230  1.00 20.00           H
ATOM     93  HB3 ARG B  60     101.814   3.099  -2.563  1.00 20.00           H
ATOM     94  HG2 ARG B  60     101.136   1.170  -4.029  1.00 20.00           H
ATOM     95  HG3 ARG B  60     101.447   0.054  -2.698  1.00 20.00           H
ATOM     96  HD2 ARG B  60     103.348   0.133  -4.269  1.00 20.00           H
ATOM     97  HD3 ARG B  60     103.505   1.880  -4.062  1.00 20.00           H
ATOM     98  HE  ARG B  60     103.674   0.627  -1.479  1.00 20.00           H
ATOM     99 HH11 ARG B  60     107.271   0.331  -1.484  1.00 20.00           H
ATOM    100 HH12 ARG B  60     105.858   0.307  -0.476  1.00 20.00           H
ATOM    101 HH21 ARG B  60     105.530   0.906  -4.614  1.00 20.00           H
ATOM    102 HH22 ARG B  60     107.088   0.675  -3.874  1.00 20.00           H
END
'''

# Microheterogeneity at B 10: shared blank-altloc backbone (written as ASP), side
# chains Asp (altloc A, no HD2) and Asn (altloc B).
microhet_model_str = '''
CRYST1  200.000  200.000  200.000  90.00  90.00  90.00 P 1
ATOM      1  N   ASP B  10      -0.317   1.688   0.066  1.00 20.00           N
ATOM      2  CA  ASP B  10      -0.470   0.286  -0.344  1.00 20.00           C
ATOM      3  C   ASP B  10      -1.868  -0.180  -0.029  1.00 20.00           C
ATOM      4  O   ASP B  10      -2.534   0.415   0.786  1.00 20.00           O
ATOM      5  H   ASP B  10      -0.928   2.289  -0.467  1.00 20.00           H
ATOM      6  HA  ASP B  10      -0.292   0.199  -1.416  1.00 20.00           H
ATOM      7  CB AASP B  10       0.539  -0.580   0.413  0.50 20.00           C
ATOM      8  CG AASP B  10       1.938  -0.195   0.004  0.50 20.00           C
ATOM      9  OD1AASP B  10       2.109   0.681  -0.810  0.50 20.00           O
ATOM     10  OD2AASP B  10       2.992  -0.826   0.543  0.50 20.00           O
ATOM     11  HB2AASP B  10       0.419  -0.425   1.485  0.50 20.00           H
ATOM     12  HB3AASP B  10       0.367  -1.630   0.176  0.50 20.00           H
ATOM     13  CB BASN B  10       0.862  -0.588   0.401  0.50 20.00           C
ATOM     14  CG BASN B  10       2.260  -0.197  -0.002  0.50 20.00           C
ATOM     15  OD1BASN B  10       2.432   0.697  -0.804  0.50 20.00           O
ATOM     16  ND2BASN B  10       3.319  -0.841   0.527  0.50 20.00           N
ATOM     17  HB2BASN B  10       0.742  -0.451   1.476  0.50 20.00           H
ATOM     18  HB3BASN B  10       0.689  -1.633   0.146  0.50 20.00           H
ATOM     19 HD21BASN B  10       3.181  -1.556   1.168  0.50 20.00           H
ATOM     20 HD22BASN B  10       4.219  -0.590   0.268  0.50 20.00           H
END
'''

# NH4 A 3 and Asp G 60 from salt_model_str; Asp split into A and B (B moved by
# 0.4 A).
altloc_asp_model_str = '''
CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1
HETATM    1  N   NH4 A   3      -9.982  30.026   0.000  1.00 20.00           N
HETATM    2  HN1 NH4 A   3     -10.367  29.481  -0.771  1.00 20.00           H
HETATM    3  HN2 NH4 A   3      -8.962  30.026   0.000  1.00 20.00           H
HETATM    4  HN3 NH4 A   3     -10.322  30.988   0.000  1.00 20.00           H
HETATM    5  HN4 NH4 A   3     -10.367  29.481   0.771  1.00 20.00           H
ATOM      6  N  AASP G  60      -7.514  34.467   2.148  0.50 20.00           N
ATOM      7  CA AASP G  60      -8.545  34.210   3.162  0.50 20.00           C
ATOM      8  C  AASP G  60      -8.453  35.251   4.247  0.50 20.00           C
ATOM      9  O  AASP G  60      -7.433  35.883   4.393  0.50 20.00           O
ATOM     10  CB AASP G  60      -8.331  32.821   3.767  0.50 20.00           C
ATOM     11  CG AASP G  60      -8.544  31.772   2.706  0.50 20.00           C
ATOM     12  OD1AASP G  60      -8.838  32.102   1.582  0.50 20.00           O
ATOM     13  OD2AASP G  60      -8.408  30.473   3.011  0.50 20.00           O
ATOM     14  H  AASP G  60      -7.672  35.351   1.688  0.50 20.00           H
ATOM     15  HA AASP G  60      -9.531  34.256   2.698  0.50 20.00           H
ATOM     16  HB2AASP G  60      -7.314  32.746   4.153  0.50 20.00           H
ATOM     17  HB3AASP G  60      -9.040  32.666   4.580  0.50 20.00           H
ATOM     18  N  BASP G  60      -7.514  34.467   2.548  0.50 20.00           N
ATOM     19  CA BASP G  60      -8.545  34.210   3.562  0.50 20.00           C
ATOM     20  C  BASP G  60      -8.453  35.251   4.647  0.50 20.00           C
ATOM     21  O  BASP G  60      -7.433  35.883   4.793  0.50 20.00           O
ATOM     22  CB BASP G  60      -8.331  32.821   4.167  0.50 20.00           C
ATOM     23  CG BASP G  60      -8.544  31.772   3.106  0.50 20.00           C
ATOM     24  OD1BASP G  60      -8.838  32.102   1.982  0.50 20.00           O
ATOM     25  OD2BASP G  60      -8.408  30.473   3.411  0.50 20.00           O
ATOM     26  H  BASP G  60      -7.672  35.351   2.088  0.50 20.00           H
ATOM     27  HA BASP G  60      -9.531  34.256   3.098  0.50 20.00           H
ATOM     28  HB2BASP G  60      -7.314  32.746   4.553  0.50 20.00           H
ATOM     29  HB3BASP G  60      -9.040  32.666   4.980  0.50 20.00           H
END
'''

# NH4 A 3 and Asp G 60 both split into A and B (the B conformers moved together by
# 0.8 A, so NH4 A...Asp B is also within 4 A).
altloc_both_model_str = '''
CRYST1  100.000  100.000  100.000  90.00  90.00  90.00 P 1
HETATM    1  N  ANH4 A   3      -9.982  30.026   0.000  0.50 20.00           N
HETATM    2  HN1ANH4 A   3     -10.367  29.481  -0.771  0.50 20.00           H
HETATM    3  HN2ANH4 A   3      -8.962  30.026   0.000  0.50 20.00           H
HETATM    4  HN3ANH4 A   3     -10.322  30.988   0.000  0.50 20.00           H
HETATM    5  HN4ANH4 A   3     -10.367  29.481   0.771  0.50 20.00           H
HETATM    6  N  BNH4 A   3      -9.982  30.826   0.000  0.50 20.00           N
HETATM    7  HN1BNH4 A   3     -10.367  30.281  -0.771  0.50 20.00           H
HETATM    8  HN2BNH4 A   3      -8.962  30.826   0.000  0.50 20.00           H
HETATM    9  HN3BNH4 A   3     -10.322  31.788   0.000  0.50 20.00           H
HETATM   10  HN4BNH4 A   3     -10.367  30.281   0.771  0.50 20.00           H
ATOM     11  N  AASP G  60      -7.514  34.467   2.148  0.50 20.00           N
ATOM     12  CA AASP G  60      -8.545  34.210   3.162  0.50 20.00           C
ATOM     13  C  AASP G  60      -8.453  35.251   4.247  0.50 20.00           C
ATOM     14  O  AASP G  60      -7.433  35.883   4.393  0.50 20.00           O
ATOM     15  CB AASP G  60      -8.331  32.821   3.767  0.50 20.00           C
ATOM     16  CG AASP G  60      -8.544  31.772   2.706  0.50 20.00           C
ATOM     17  OD1AASP G  60      -8.838  32.102   1.582  0.50 20.00           O
ATOM     18  OD2AASP G  60      -8.408  30.473   3.011  0.50 20.00           O
ATOM     19  H  AASP G  60      -7.672  35.351   1.688  0.50 20.00           H
ATOM     20  HA AASP G  60      -9.531  34.256   2.698  0.50 20.00           H
ATOM     21  HB2AASP G  60      -7.314  32.746   4.153  0.50 20.00           H
ATOM     22  HB3AASP G  60      -9.040  32.666   4.580  0.50 20.00           H
ATOM     23  N  BASP G  60      -7.514  35.267   2.148  0.50 20.00           N
ATOM     24  CA BASP G  60      -8.545  35.010   3.162  0.50 20.00           C
ATOM     25  C  BASP G  60      -8.453  36.051   4.247  0.50 20.00           C
ATOM     26  O  BASP G  60      -7.433  36.683   4.393  0.50 20.00           O
ATOM     27  CB BASP G  60      -8.331  33.621   3.767  0.50 20.00           C
ATOM     28  CG BASP G  60      -8.544  32.572   2.706  0.50 20.00           C
ATOM     29  OD1BASP G  60      -8.838  32.902   1.582  0.50 20.00           O
ATOM     30  OD2BASP G  60      -8.408  31.273   3.011  0.50 20.00           O
ATOM     31  H  BASP G  60      -7.672  36.151   1.688  0.50 20.00           H
ATOM     32  HA BASP G  60      -9.531  35.056   2.698  0.50 20.00           H
ATOM     33  HB2BASP G  60      -7.314  33.546   4.153  0.50 20.00           H
ATOM     34  HB3BASP G  60      -9.040  33.466   4.580  0.50 20.00           H
END
'''

# Acetic acid (neutral; HXT on OXT) restraints for neutral_model_str.
acy_neutral_cif = ('acy_neutral.cif', '''
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 ACY  ACY  'acetic acid' ligand 8 4 .

data_comp_ACY
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
 ACY  C    C  C     0
 ACY  O    O  O     0
 ACY  OXT  O  OH1   0
 ACY  CH3  C  CH3   0
 ACY  HXT  H  HOH1  0
 ACY  H1   H  HCH3  0
 ACY  H2   H  HCH3  0
 ACY  H3   H  HCH3  0

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
_chem_comp_bond.value_dist_neutron
 ACY  C    O    double  1.210  0.020  1.210
 ACY  C    OXT  single  1.310  0.020  1.310
 ACY  C    CH3  single  1.500  0.020  1.500
 ACY  OXT  HXT  single  0.850  0.020  0.980
 ACY  CH3  H1   single  0.970  0.020  1.090
 ACY  CH3  H2   single  0.970  0.020  1.090
 ACY  CH3  H3   single  0.970  0.020  1.090

loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
 ACY  CH3  C    OXT  113.00  3.000
 ACY  CH3  C    O    124.00  3.000
 ACY  OXT  C    O    123.00  3.000
 ACY  C    OXT  HXT  106.00  3.000
 ACY  H3   CH3  H2   109.50  3.000
 ACY  H3   CH3  H1   109.50  3.000
 ACY  H2   CH3  H1   109.50  3.000
 ACY  H3   CH3  C    109.50  3.000
 ACY  H2   CH3  C    109.50  3.000
 ACY  H1   CH3  C    109.50  3.000

loop_
_chem_comp_plane_atom.comp_id
_chem_comp_plane_atom.plane_id
_chem_comp_plane_atom.atom_id
_chem_comp_plane_atom.dist_esd
 ACY  plan-1  C    0.020
 ACY  plan-1  O    0.020
 ACY  plan-1  OXT  0.020
 ACY  plan-1  CH3  0.020
''')

# GeoStd ACT (acetate) with OXT's formal charge set to 0: conflicts with the
# carboxylate perceived from the model.
act_conflict_cif = ('act_conflict.cif', '''
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
_chem_comp.initial_date
_chem_comp.modified_date
_chem_comp.source
 ACT  ACT  'acetate                  '  ligand  7  4  .  2016-12-20  2022-03-11
;
Directly from eLBOW using geometry from QM method PBEh-3c with CPCM solvent model
Validated by Mogul as PERFECT
;

data_comp_ACT
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
 ACT  C    C  C      0   0.385  50.6323  -7.1466  41.5339
 ACT  O    O  O      0  -0.623  50.8984  -7.3786  42.7319
 ACT  OXT  O  OC     0  -0.618  51.4443  -7.0878  40.5890
 ACT  CH3  C  CH3    0  -0.780  49.1492  -6.8839  41.2121
 ACT  H1   H  HCH3   0   0.201  48.9411  -6.8660  40.1424
 ACT  H2   H  HCH3   0   0.213  48.5069  -7.6339  41.6777
 ACT  H3   H  HCH3   0   0.221  48.8449  -5.9188  41.6260

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
_chem_comp_bond.value_dist_neutron
 ACT  C    O    deloc   1.249  0.020  1.249
 ACT  C    OXT  deloc   1.247  0.020  1.247
 ACT  C    CH3  single  1.540  0.020  1.540
 ACT  CH3  H1   single  0.970  0.020  1.090
 ACT  CH3  H2   single  0.970  0.020  1.090
 ACT  CH3  H3   single  0.970  0.020  1.090

loop_
_chem_comp_angle.comp_id
_chem_comp_angle.atom_id_1
_chem_comp_angle.atom_id_2
_chem_comp_angle.atom_id_3
_chem_comp_angle.value_angle
_chem_comp_angle.value_angle_esd
 ACT  CH3  C    OXT  117.43  3.000
 ACT  CH3  C    O    115.93  3.000
 ACT  OXT  C    O    126.64  3.000
 ACT  H3   CH3  H2   106.33  3.000
 ACT  H3   CH3  H1   107.69  3.000
 ACT  H2   CH3  H1   108.52  3.000
 ACT  H3   CH3  C    109.84  3.000
 ACT  H2   CH3  C    111.12  3.000
 ACT  H1   CH3  C    113.06  3.000

loop_
_chem_comp_tor.comp_id
_chem_comp_tor.id
_chem_comp_tor.atom_id_1
_chem_comp_tor.atom_id_2
_chem_comp_tor.atom_id_3
_chem_comp_tor.atom_id_4
_chem_comp_tor.value_angle
_chem_comp_tor.value_angle_esd
_chem_comp_tor.period
 ACT  Var_01  H1  CH3  C  O  -169.61  30.0  3

loop_
_chem_comp_plane_atom.comp_id
_chem_comp_plane_atom.plane_id
_chem_comp_plane_atom.atom_id
_chem_comp_plane_atom.dist_esd
 ACT  plan-1  C    0.020
 ACT  plan-1  O    0.020
 ACT  plan-1  OXT  0.020
 ACT  plan-1  CH3  0.020

''')

# A probe2-only H-bond (D-H...A 113 deg, below pnp's 120): EDO A 1 O1-HO1 donates
# to EDO A 2 O1. EDO A 1 has conformers A and B, B shifted 0.45 A along O1->HO1,
# so conformer B's O1 is the heavy atom nearest to conformer A's HO1.
alt_model_str = '''
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
HETATM    1  C1 AEDO A   1      15.177  10.059  15.375  0.50 20.00           C
HETATM    2  O1 AEDO A   1      16.075  10.413  14.347  0.50 20.00           O
HETATM    3  C2 AEDO A   1      13.776  10.464  14.992  0.50 20.00           C
HETATM    4  O2 AEDO A   1      13.340   9.696  13.892  0.50 20.00           O
HETATM    5  H11AEDO A   1      15.198   8.980  15.580  0.50 20.00           H
HETATM    6  H12AEDO A   1      15.418  10.564  16.319  0.50 20.00           H
HETATM    7  HO1AEDO A   1      16.953  10.118  14.602  0.50 20.00           H
HETATM    8  H21AEDO A   1      13.139  10.312  15.873  0.50 20.00           H
HETATM    9  H22AEDO A   1      13.758  11.539  14.771  0.50 20.00           H
HETATM   10  HO2AEDO A   1      12.463  10.001  13.644  0.50 20.00           H
HETATM   11  C1 BEDO A   1      15.588   9.921  15.494  0.50 20.00           C
HETATM   12  O1 BEDO A   1      16.486  10.275  14.466  0.50 20.00           O
HETATM   13  C2 BEDO A   1      14.187  10.326  15.111  0.50 20.00           C
HETATM   14  O2 BEDO A   1      13.751   9.558  14.011  0.50 20.00           O
HETATM   15  H11BEDO A   1      15.609   8.842  15.699  0.50 20.00           H
HETATM   16  H12BEDO A   1      15.829  10.426  16.438  0.50 20.00           H
HETATM   17  HO1BEDO A   1      17.364   9.980  14.721  0.50 20.00           H
HETATM   18  H21BEDO A   1      13.550  10.174  15.992  0.50 20.00           H
HETATM   19  H22BEDO A   1      14.169  11.401  14.890  0.50 20.00           H
HETATM   20  HO2BEDO A   1      12.874   9.863  13.763  0.50 20.00           H
HETATM   21  C1  EDO A   2      19.271  12.379  15.630  1.00 20.00           C
HETATM   22  O1  EDO A   2      18.382  11.848  14.673  1.00 20.00           O
HETATM   23  C2  EDO A   2      20.053  11.264  16.277  1.00 20.00           C
HETATM   24  O2  EDO A   2      20.908  10.666  15.327  1.00 20.00           O
HETATM   25  H11 EDO A   2      19.970  13.098  15.181  1.00 20.00           H
HETATM   26  H12 EDO A   2      18.740  12.913  16.428  1.00 20.00           H
HETATM   27  HO1 EDO A   2      17.920  12.579  14.256  1.00 20.00           H
HETATM   28  H21 EDO A   2      20.615  11.700  17.113  1.00 20.00           H
HETATM   29  H22 EDO A   2      19.355  10.537  16.711  1.00 20.00           H
HETATM   30  HO2 EDO A   2      21.363   9.935  15.753  1.00 20.00           H
END
'''

# Names probe2's raw atom field runs together: a 4-character atom name with an
# altloc, a 4-digit residue number with a one-letter chain, an insertion code.
names_str = '''
ATOM      1 HD21AASN A1001B      1.000   2.000   3.000  1.00 20.00           H
ATOM      2 HD21BASN A1001B      1.300   2.000   3.000  1.00 20.00           H
ATOM      3  ND2 ASN A1001B      1.500   2.500   3.000  1.00 20.00           N
ATOM      4  CA  GLY B  52A      5.000   2.000   3.000  1.00 20.00           C
END
'''

# Gly-Asp-Lys-Arg-His-Gly, beta strand from the sequence (ss_idealization), H by
# reduce2 (N-terminal H1-H3, His HD1 and HE2; no OXT)
gdkrhg_model_str = '''
CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1
ATOM      1  N   GLY A   1      27.961   0.504   1.988  1.00  0.00           N
ATOM      2  CA  GLY A   1      29.153   0.205   2.773  1.00  0.00           C
ATOM      3  C   GLY A   1      30.420   0.562   2.003  1.00  0.00           C
ATOM      4  O   GLY A   1      30.753  -0.077   1.005  1.00  0.00           O
ATOM      5  H1  GLY A   1      27.977   0.034   1.233  1.00  0.00           H
ATOM      6  H2  GLY A   1      27.236   0.288   2.456  1.00  0.00           H
ATOM      7  H3  GLY A   1      27.942   1.373   1.796  1.00  0.00           H
ATOM      8  HA2 GLY A   1      29.136   0.714   3.598  1.00  0.00           H
ATOM      9  HA3 GLY A   1      29.174  -0.741   2.986  1.00  0.00           H
ATOM     10  N   ASP A   2      31.123   1.587   2.474  1.00  0.00           N
ATOM     11  CA  ASP A   2      32.355   2.031   1.832  1.00  0.00           C
ATOM     12  C   ASP A   2      33.552   1.851   2.758  1.00  0.00           C
ATOM     13  O   ASP A   2      33.675   2.539   3.772  1.00  0.00           O
ATOM     14  CB  ASP A   2      32.236   3.500   1.416  1.00  0.00           C
ATOM     15  CG  ASP A   2      33.431   3.977   0.614  1.00  0.00           C
ATOM     16  OD1 ASP A   2      33.540   3.602  -0.573  1.00  0.00           O
ATOM     17  OD2 ASP A   2      34.261   4.727   1.168  1.00  0.00           O
ATOM     18  H   ASP A   2      30.906   2.045   3.169  1.00  0.00           H
ATOM     19  HA  ASP A   2      32.507   1.498   1.036  1.00  0.00           H
ATOM     20  HB3 ASP A   2      32.169   4.050   2.212  1.00  0.00           H
ATOM     21  HB2 ASP A   2      31.443   3.611   0.869  1.00  0.00           H
ATOM     22  N   LYS A   3      34.434   0.922   2.403  1.00  0.00           N
ATOM     23  CA  LYS A   3      35.624   0.650   3.202  1.00  0.00           C
ATOM     24  C   LYS A   3      36.893   0.975   2.422  1.00  0.00           C
ATOM     25  O   LYS A   3      37.227   0.298   1.450  1.00  0.00           O
ATOM     26  CB  LYS A   3      35.644  -0.817   3.641  1.00  0.00           C
ATOM     27  CG  LYS A   3      34.550  -1.184   4.630  1.00  0.00           C
ATOM     28  CD  LYS A   3      34.642  -2.643   5.046  1.00  0.00           C
ATOM     29  CE  LYS A   3      33.535  -3.015   6.019  1.00  0.00           C
ATOM     30  NZ  LYS A   3      33.606  -4.445   6.426  1.00  0.00           N
ATOM     31  H   LYS A   3      34.366   0.432   1.699  1.00  0.00           H
ATOM     32  HA  LYS A   3      35.607   1.207   3.996  1.00  0.00           H
ATOM     33  HB3 LYS A   3      36.497  -1.003   4.063  1.00  0.00           H
ATOM     34  HB2 LYS A   3      35.532  -1.377   2.857  1.00  0.00           H
ATOM     35  HG2 LYS A   3      34.637  -0.634   5.424  1.00  0.00           H
ATOM     36  HG3 LYS A   3      33.683  -1.038   4.219  1.00  0.00           H
ATOM     37  HD3 LYS A   3      35.495  -2.800   5.480  1.00  0.00           H
ATOM     38  HD2 LYS A   3      34.560  -3.207   4.261  1.00  0.00           H
ATOM     39  HE2 LYS A   3      33.615  -2.469   6.817  1.00  0.00           H
ATOM     40  HE3 LYS A   3      32.675  -2.863   5.597  1.00  0.00           H
ATOM     41  HZ1 LYS A   3      33.529  -4.969   5.711  1.00  0.00           H
ATOM     42  HZ3 LYS A   3      32.947  -4.633   6.994  1.00  0.00           H
ATOM     43  HZ2 LYS A   3      34.385  -4.610   6.823  1.00  0.00           H
ATOM     44  N   ARG A   4      37.598   2.015   2.856  1.00  0.00           N
ATOM     45  CA  ARG A   4      38.832   2.432   2.201  1.00  0.00           C
ATOM     46  C   ARG A   4      40.027   2.284   3.136  1.00  0.00           C
ATOM     47  O   ARG A   4      40.148   3.009   4.124  1.00  0.00           O
ATOM     48  CB  ARG A   4      38.719   3.886   1.729  1.00 10.00           C
ATOM     49  CG  ARG A   4      37.617   4.134   0.703  1.00 10.00           C
ATOM     50  CD  ARG A   4      37.883   3.436  -0.629  1.00 10.00           C
ATOM     51  NE  ARG A   4      39.078   3.951  -1.302  1.00 10.00           N
ATOM     52  CZ  ARG A   4      39.739   3.351  -2.295  1.00 10.00           C
ATOM     53  NH1 ARG A   4      39.359   2.176  -2.796  1.00 10.00           N
ATOM     54  NH2 ARG A   4      40.810   3.945  -2.804  1.00 10.00           N
ATOM     55  H   ARG A   4      37.381   2.499   3.533  1.00  0.00           H
ATOM     56  HA  ARG A   4      38.983   1.869   1.426  1.00  0.00           H
ATOM     57  HB3 ARG A   4      39.561   4.145   1.323  1.00 10.00           H
ATOM     58  HB2 ARG A   4      38.535   4.447   2.498  1.00 10.00           H
ATOM     59  HG3 ARG A   4      37.550   5.087   0.534  1.00 10.00           H
ATOM     60  HG2 ARG A   4      36.777   3.800   1.054  1.00 10.00           H
ATOM     61  HD2 ARG A   4      38.012   2.488  -0.470  1.00 10.00           H
ATOM     62  HD3 ARG A   4      37.124   3.574  -1.217  1.00 10.00           H
ATOM     63  HE  ARG A   4      39.382   4.709  -1.033  1.00 10.00           H
ATOM     64 HH11 ARG A   4      39.807   1.816  -3.436  1.00 10.00           H
ATOM     65 HH12 ARG A   4      38.666   1.777  -2.479  1.00 10.00           H
ATOM     66 HH21 ARG A   4      41.246   3.571  -3.444  1.00 10.00           H
ATOM     67 HH22 ARG A   4      41.069   4.704  -2.494  1.00 10.00           H
ATOM     68  N   HIS A   5      40.907   1.340   2.818  1.00  0.00           N
ATOM     69  CA  HIS A   5      42.094   1.095   3.630  1.00  0.00           C
ATOM     70  C   HIS A   5      43.367   1.388   2.841  1.00  0.00           C
ATOM     71  O   HIS A   5      43.701   0.675   1.896  1.00  0.00           O
ATOM     72  CB  HIS A   5      42.110  -0.354   4.123  1.00  0.00           C
ATOM     73  CG  HIS A   5      40.982  -0.690   5.047  1.00  0.00           C
ATOM     74  ND1 HIS A   5      40.946  -0.276   6.360  1.00  0.00           N
ATOM     75  CD2 HIS A   5      39.847  -1.402   4.846  1.00  0.00           C
ATOM     76  CE1 HIS A   5      39.840  -0.718   6.930  1.00  0.00           C
ATOM     77  NE2 HIS A   5      39.154  -1.404   6.033  1.00  0.00           N
ATOM     78  H   HIS A   5      40.839   0.825   2.133  1.00  0.00           H
ATOM     79  HA  HIS A   5      42.076   1.681   4.403  1.00  0.00           H
ATOM     80  HB3 HIS A   5      42.940  -0.512   4.600  1.00  0.00           H
ATOM     81  HB2 HIS A   5      42.050  -0.945   3.356  1.00  0.00           H
ATOM     82  HD1 HIS A   5      41.550   0.197   6.748  1.00  0.00           H
ATOM     83  HD2 HIS A   5      39.586  -1.812   4.053  1.00  0.00           H
ATOM     84  HE1 HIS A   5      39.587  -0.571   7.813  1.00  0.00           H
ATOM     85  HE2 HIS A   5      38.397  -1.789   6.170  1.00  0.00           H
ATOM     86  N   GLY A   6      44.072   2.443   3.238  1.00  0.00           N
ATOM     87  CA  GLY A   6      45.309   2.832   2.570  1.00  0.00           C
ATOM     88  C   GLY A   6      46.501   2.716   3.514  1.00  0.00           C
ATOM     89  O   GLY A   6      46.622   3.477   4.474  1.00  0.00           O
ATOM     90  H   GLY A   6      43.854   2.953   3.895  1.00  0.00           H
ATOM     91  HA3 GLY A   6      45.240   3.751   2.266  1.00  0.00           H
ATOM     92  HA2 GLY A   6      45.461   2.256   1.804  1.00  0.00           H
END
'''

# DC 1 - DG 2 from tst_add_hydrogen_16 (DNA control), H by reduce2
dna_model_str = '''
CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1
ATOM      1  O5'  DC A   1      19.545  18.136  17.917  1.00  3.07           O
ATOM      2  C5'  DC A   1      19.769  17.119  18.884  1.00  2.46           C
ATOM      3  C4'  DC A   1      18.610  16.148  19.001  1.00  2.19           C
ATOM      4  O4'  DC A   1      17.462  16.852  19.514  1.00  2.62           O
ATOM      5  C3'  DC A   1      18.161  15.506  17.674  1.00  2.26           C
ATOM      6  O3'  DC A   1      17.782  14.139  17.875  1.00  2.47           O
ATOM      7  C2'  DC A   1      16.906  16.282  17.315  1.00  2.41           C
ATOM      8  C1'  DC A   1      16.340  16.624  18.692  1.00  2.34           C
ATOM      9  N1   DC A   1      15.516  17.837  18.704  1.00  2.30           N
ATOM     10  C2   DC A   1      14.145  17.720  18.492  1.00  2.34           C
ATOM     11  O2   DC A   1      13.658  16.581  18.329  1.00  3.07           O
ATOM     12  N3   DC A   1      13.385  18.831  18.454  1.00  2.35           N
ATOM     13  C4   DC A   1      13.943  20.039  18.611  1.00  2.37           C
ATOM     14  N4   DC A   1      13.161  21.108  18.580  1.00  2.79           N
ATOM     15  C5   DC A   1      15.357  20.189  18.812  1.00  2.71           C
ATOM     16  C6   DC A   1      16.103  19.061  18.841  1.00  2.68           C
ATOM     17  H5'  DC A   1      19.913  17.537  19.747  1.00  2.46           H
ATOM     18 H5''  DC A   1      20.566  16.624  18.636  1.00  2.46           H
ATOM     19  H4'  DC A   1      18.857  15.438  19.613  1.00  2.19           H
ATOM     20  H3'  DC A   1      18.861  15.595  17.008  1.00  2.26           H
ATOM     21 H2''  DC A   1      17.117  17.083  16.809  1.00  2.41           H
ATOM     22  H2'  DC A   1      16.287  15.738  16.802  1.00  2.41           H
ATOM     23  H1'  DC A   1      15.829  15.876  19.039  1.00  2.34           H
ATOM     24  H41  DC A   1      12.314  21.018  18.461  1.00  2.79           H
ATOM     25  H42  DC A   1      13.500  21.892  18.679  1.00  2.79           H
ATOM     26  H5   DC A   1      15.744  21.028  18.918  1.00  2.71           H
ATOM     27  H6   DC A   1      17.024  19.120  18.955  1.00  2.68           H
ATOM     28  P    DG A   2      18.825  12.942  17.684  1.00  2.51           P
ATOM     29  OP1  DG A   2      19.788  13.206  16.573  1.00  3.24           O
ATOM     30  OP2  DG A   2      17.976  11.710  17.621  1.00  3.25           O
ATOM     31  O5'  DG A   2      19.719  12.937  19.002  1.00  2.66           O
ATOM     32  C5'  DG A   2      19.103  12.733  20.284  1.00  2.85           C
ATOM     33  C4'  DG A   2      20.140  13.045  21.335  1.00  2.73           C
ATOM     34  O4'  DG A   2      20.546  14.399  21.207  1.00  2.76           O
ATOM     35  C3'  DG A   2      19.598  12.910  22.753  1.00  2.85           C
ATOM     36  O3'  DG A   2      19.812  11.563  23.230  1.00  3.51           O
ATOM     37  C2'  DG A   2      20.430  13.919  23.526  1.00  3.30           C
ATOM     38  C1'  DG A   2      20.834  14.964  22.481  1.00  2.85           C
ATOM     39  N9   DG A   2      20.140  16.222  22.572  1.00  2.76           N
ATOM     40  C8   DG A   2      20.744  17.451  22.654  1.00  3.56           C
ATOM     41  N7   DG A   2      19.903  18.441  22.626  1.00  3.75           N
ATOM     42  C5   DG A   2      18.658  17.840  22.510  1.00  2.68           C
ATOM     43  C6   DG A   2      17.359  18.399  22.397  1.00  2.62           C
ATOM     44  O6   DG A   2      17.067  19.612  22.382  1.00  3.52           O
ATOM     45  N1   DG A   2      16.371  17.435  22.313  1.00  2.19           N
ATOM     46  C2   DG A   2      16.610  16.084  22.263  1.00  1.99           C
ATOM     47  N2   DG A   2      15.540  15.294  22.108  1.00  2.48           N
ATOM     48  N3   DG A   2      17.814  15.534  22.356  1.00  2.10           N
ATOM     49  C4   DG A   2      18.782  16.466  22.477  1.00  2.26           C
ATOM     50 H5''  DG A   2      18.338  13.321  20.384  1.00  2.85           H
ATOM     51  H5'  DG A   2      18.807  11.813  20.370  1.00  2.85           H
ATOM     52  H4'  DG A   2      20.914  12.473  21.215  1.00  2.73           H
ATOM     53  H3'  DG A   2      18.657  13.137  22.813  1.00  2.85           H
ATOM     54 H2''  DG A   2      19.908  14.325  24.236  1.00  3.30           H
ATOM     55  H2'  DG A   2      21.213  13.497  23.912  1.00  3.30           H
ATOM     56  H1'  DG A   2      21.789  15.119  22.557  1.00  2.85           H
ATOM     57  H8   DG A   2      21.664  17.564  22.723  1.00  3.56           H
ATOM     58  H1   DG A   2      15.553  17.699  22.291  1.00  2.19           H
ATOM     59  H21  DG A   2      14.755  15.640  22.047  1.00  2.48           H
ATOM     60  H22  DG A   2      15.638  14.440  22.070  1.00  2.48           H
END
'''

# N-methylpyridinium and nitromethane (RDKit embedding; formal charges, explicit orders)
zpy_cif = ('zpy.cif', '''
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 ZPY  ZPY  'ZPY' ligand 15 7 .

data_comp_ZPY
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
 ZPY  CM   C  C     0   0.000  -2.1240  -0.2930  -0.3397
 ZPY  N1   N  N     1   0.000  -0.6851  -0.0938  -0.0998
 ZPY  C2   C  C     0   0.000  -0.1762   1.1597  -0.0920
 ZPY  C3   C  C     0   0.000   1.1812   1.3698   0.1225
 ZPY  C4   C  C     0   0.000   2.0158   0.2782   0.3244
 ZPY  C5   C  C     0   0.000   1.4792  -1.0026   0.3054
 ZPY  C6   C  C     0   0.000   0.1162  -1.1675   0.0874
 ZPY  HCM1 H  H     0   0.000  -2.4531  -1.1950   0.1827
 ZPY  HCM2 H  H     0   0.000  -2.6746   0.5688   0.0462
 ZPY  HCM3 H  H     0   0.000  -2.2648  -0.3943  -1.4186
 ZPY  HC21 H  H     0   0.000  -0.8620   1.9853  -0.2597
 ZPY  HC31 H  H     0   0.000   1.5872   2.3792   0.1305
 ZPY  HC41 H  H     0   0.000   3.0819   0.4252   0.4935
 ZPY  HC51 H  H     0   0.000   2.1208  -1.8680   0.4580
 ZPY  HC61 H  H     0   0.000  -0.3423  -2.1520   0.0593

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
 ZPY  CM   N1   single   1.472  0.020
 ZPY  N1   C2   aromatic  1.353  0.020
 ZPY  C2   C3   aromatic  1.390  0.020
 ZPY  C3   C4   aromatic  1.389  0.020
 ZPY  C4   C5   aromatic  1.389  0.020
 ZPY  C5   C6   aromatic  1.390  0.020
 ZPY  C6   N1   aromatic  1.353  0.020
 ZPY  CM   HCM1 single   1.093  0.020
 ZPY  CM   HCM2 single   1.093  0.020
 ZPY  CM   HCM3 single   1.093  0.020
 ZPY  C2   HC21 single   1.086  0.020
 ZPY  C3   HC31 single   1.088  0.020
 ZPY  C4   HC41 single   1.089  0.020
 ZPY  C5   HC51 single   1.088  0.020
 ZPY  C6   HC61 single   1.086  0.020
''')

znx_cif = ('znx.cif', '''
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 ZNX  ZNX  'ZNX' ligand 7 4 .

data_comp_ZNX
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
 ZNX  C1   C  C     0   0.000  -0.6519  -0.0458   0.0124
 ZNX  N1   N  N     1   0.000   0.8315   0.0577  -0.0234
 ZNX  O1   O  O     0   0.000   1.4709  -0.9962   0.0763
 ZNX  O2   O  O    -1   0.000   1.3133   1.1917  -0.1300
 ZNX  HC11 H  H     0   0.000  -0.9554   0.0310   1.0581
 ZNX  HC12 H  H     0   0.000  -0.9400  -1.0097  -0.4127
 ZNX  HC13 H  H     0   0.000  -1.0683   0.7713  -0.5807

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
 ZNX  C1   N1   single   1.487  0.020
 ZNX  N1   O1   double   1.237  0.020
 ZNX  N1   O2   single   1.237  0.020
 ZNX  C1   HC11 single   1.092  0.020
 ZNX  C1   HC12 single   1.092  0.020
 ZNX  C1   HC13 single   1.092  0.020
''')

# GeoStd BEN (neutral benzamidine) superposed with C, N1, N2 on ACT 2's C, O, OXT
# of salt_model_str (least squares)
ben_at_act2_lines = '''
HETATM  900  C1  BEN A   2      -8.557  15.920 -29.956  1.00 20.00           C
HETATM  901  C2  BEN A   2      -7.932  14.966 -29.155  1.00 20.00           C
HETATM  902  C3  BEN A   2      -6.546  14.904 -29.093  1.00 20.00           C
HETATM  903  C4  BEN A   2      -5.775  15.783 -29.840  1.00 20.00           C
HETATM  904  C5  BEN A   2      -6.393  16.732 -30.645  1.00 20.00           C
HETATM  905  C6  BEN A   2      -7.777  16.804 -30.697  1.00 20.00           C
HETATM  906  C   BEN A   2     -10.043  16.029 -30.000  1.00 20.00           C
HETATM  907  N1  BEN A   2     -10.688  17.132 -30.000  1.00 20.00           N
HETATM  908  N2  BEN A   2     -10.705  14.840 -30.000  1.00 20.00           N
HETATM  909  H2  BEN A   2      -8.524  14.278 -28.565  1.00 20.00           H
HETATM  910  H3  BEN A   2      -6.071  14.168 -28.458  1.00 20.00           H
HETATM  911  H4  BEN A   2      -4.695  15.731 -29.797  1.00 20.00           H
HETATM  912  H5  BEN A   2      -5.798  17.416 -31.236  1.00 20.00           H
HETATM  913  H6  BEN A   2      -8.253  17.542 -31.331  1.00 20.00           H
HETATM  914  HN1 BEN A   2     -10.061  17.919 -29.891  1.00 20.00           H
HETATM  915 HN21 BEN A   2     -11.703  14.867 -30.130  1.00 20.00           H
HETATM  916 HN22 BEN A   2     -10.246  14.018 -30.353  1.00 20.00           H
'''

# methyl phosphate, neutral (HOP2, HOP3; RDKit embedding, explicit orders, formal
# charges 0), superposed with P, O1P, O2P on ACT 1's C, O, OXT of salt_model_str
zmp_cif = ('zmp.cif', '''
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 ZMP  ZMP  'ZMP' ligand 11 6 .

data_comp_ZMP
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
 ZMP  C1   C  C      0  0.000   1.6955  -0.2505   0.2106
 ZMP  O1   O  O      0  0.000   0.4623  -0.8150  -0.1845
 ZMP  P    P  P      0  0.000  -0.7439   0.2030  -0.4358
 ZMP  O1P  O  O      0  0.000  -0.4609   1.3912  -1.2861
 ZMP  O2P  O  O      0  0.000  -1.9356  -0.7349  -0.9257
 ZMP  O3P  O  O      0  0.000  -1.2561   0.6007   1.0230
 ZMP  HC11 H  H      0  0.000   1.5670   0.3465   1.1176
 ZMP  HC12 H  H      0  0.000   2.4022  -1.0585   0.4161
 ZMP  HC13 H  H      0  0.000   2.0962   0.3740  -0.5924
 ZMP  HOP2 H  H      0  0.000  -2.0560  -1.4850  -0.3172
 ZMP  HOP3 H  H      0  0.000  -1.7707   1.4284   0.9742

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
 ZMP  C1   O1   single    1.413  0.020
 ZMP  O1   P    single    1.598  0.020
 ZMP  P    O1P  double    1.488  0.020
 ZMP  P    O2P  single    1.594  0.020
 ZMP  P    O3P  single    1.596  0.020
 ZMP  C1   HC11 single    1.093  0.020
 ZMP  C1   HC12 single    1.093  0.020
 ZMP  C1   HC13 single    1.093  0.020
 ZMP  O2P  HOP2 single    0.973  0.020
 ZMP  O3P  HOP3 single    0.976  0.020
''')

zmp_at_act1_lines = '''
HETATM  900  C1  ZMP A   1     -38.040  16.940 -31.473  1.00 20.00           C
HETATM  901  O1  ZMP A   1     -38.888  15.846 -31.194  1.00 20.00           O
HETATM  902  P   ZMP A   1     -39.934  16.037 -30.000  1.00 20.00           P
HETATM  903  O1P ZMP A   1     -40.740  17.288 -30.000  1.00 20.00           O
HETATM  904  O2P ZMP A   1     -40.762  14.675 -30.000  1.00 20.00           O
HETATM  905  O3P ZMP A   1     -39.072  15.863 -28.668  1.00 20.00           O
HETATM  906 HC11 ZMP A   1     -37.478  17.224 -30.579  1.00 20.00           H
HETATM  907 HC12 ZMP A   1     -37.335  16.648 -32.256  1.00 20.00           H
HETATM  908 HC13 ZMP A   1     -38.629  17.790 -31.828  1.00 20.00           H
HETATM  909 HOP2 ZMP A   1     -40.163  13.908 -29.999  1.00 20.00           H
HETATM  910 HOP3 ZMP A   1     -39.551  16.257 -27.915  1.00 20.00           H
'''

# paracetamol (amide N, phenol) and 4-aminophenol (aniline, phenol)
zpa_cif = ('zpa.cif', '''
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 ZPA  ZPA  'ZPA' ligand 20 11 .

data_comp_ZPA
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
 ZPA  C1   C  C      0  0.000   3.7026   0.2029   0.1412
 ZPA  C2   C  C      0  0.000   2.3493  -0.4473  -0.0075
 ZPA  O1   O  O      0  0.000   2.2538  -1.6469  -0.2404
 ZPA  N1   N  N      0  0.000   1.3157   0.4628   0.1208
 ZPA  C3   C  C      0  0.000  -0.0688   0.2085   0.0358
 ZPA  C4   C  C      0  0.000  -0.6190  -1.0561  -0.1885
 ZPA  C5   C  C      0  0.000  -2.0066  -1.2270  -0.2591
 ZPA  C6   C  C      0  0.000  -2.8430  -0.1285  -0.1040
 ZPA  O2   O  O      0  0.000  -4.1986  -0.2512  -0.1656
 ZPA  C7   C  C      0  0.000  -2.3155   1.1352   0.1202
 ZPA  C8   C  C      0  0.000  -0.9303   1.3035   0.1902
 ZPA  HC11 H  H      0  0.000   3.6804   0.9818   0.9089
 ZPA  HC12 H  H      0  0.000   4.4358  -0.5501   0.4444
 ZPA  HC13 H  H      0  0.000   4.0026   0.6376  -0.8158
 ZPA  HN11 H  H      0  0.000   1.5765   1.4263   0.2876
 ZPA  HC41 H  H      0  0.000   0.0055  -1.9352  -0.3131
 ZPA  HC51 H  H      0  0.000  -2.4054  -2.2212  -0.4346
 ZPA  HO21 H  H      0  0.000  -4.4164  -1.1842  -0.3271
 ZPA  HC71 H  H      0  0.000  -2.9786   1.9871   0.2402
 ZPA  HC81 H  H      0  0.000  -0.5400   2.3020   0.3665

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
 ZPA  C1   C2   single    1.509  0.020
 ZPA  C2   O1   double    1.226  0.020
 ZPA  C2   N1   single    1.383  0.020
 ZPA  N1   C3   single    1.410  0.020
 ZPA  C3   C4   aromatic  1.397  0.020
 ZPA  C4   C5   aromatic  1.400  0.020
 ZPA  C5   C6   aromatic  1.389  0.020
 ZPA  C6   O2   single    1.363  0.020
 ZPA  C6   C7   aromatic  1.388  0.020
 ZPA  C7   C8   aromatic  1.397  0.020
 ZPA  C8   C3   aromatic  1.402  0.020
 ZPA  C1   HC11 single    1.094  0.020
 ZPA  C1   HC12 single    1.094  0.020
 ZPA  C1   HC13 single    1.093  0.020
 ZPA  N1   HN11 single    1.012  0.020
 ZPA  C4   HC41 single    1.086  0.020
 ZPA  C5   HC51 single    1.085  0.020
 ZPA  O2   HO21 single    0.972  0.020
 ZPA  C7   HC71 single    1.086  0.020
 ZPA  C8   HC81 single    1.086  0.020
''')

zap_cif = ('zap.cif', '''
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 ZAP  ZAP  'ZAP' ligand 15 8 .

data_comp_ZAP
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
 ZAP  N1   N  N      0  0.000   2.5622   0.2304  -0.1204
 ZAP  C1   C  C      0  0.000   1.1699   0.1532   0.0007
 ZAP  C2   C  C      0  0.000   0.5933  -0.4850   1.1025
 ZAP  C3   C  C      0  0.000  -0.7894  -0.6704   1.1750
 ZAP  C4   C  C      0  0.000  -1.5940  -0.2563   0.1220
 ZAP  O1   O  O      0  0.000  -2.9378  -0.4576   0.2296
 ZAP  C5   C  C      0  0.000  -1.0357   0.3269  -1.0093
 ZAP  C6   C  C      0  0.000   0.3490   0.5103  -1.0737
 ZAP  HN11 H  H      0  0.000   3.0489   0.2524   0.7702
 ZAP  HN12 H  H      0  0.000   2.8797   0.9593  -0.7512
 ZAP  HC21 H  H      0  0.000   1.2160  -0.8436   1.9173
 ZAP  HC31 H  H      0  0.000  -1.2322  -1.1493   2.0432
 ZAP  HO11 H  H      0  0.000  -3.3637  -0.1369  -0.5820
 ZAP  HC51 H  H      0  0.000  -1.6503   0.6304  -1.8502
 ZAP  HC61 H  H      0  0.000   0.7841   0.9362  -1.9738

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
 ZAP  N1   C1   single    1.400  0.020
 ZAP  C1   C2   aromatic  1.398  0.020
 ZAP  C2   C3   aromatic  1.397  0.020
 ZAP  C3   C4   aromatic  1.388  0.020
 ZAP  C4   O1   single    1.363  0.020
 ZAP  C4   C5   aromatic  1.390  0.020
 ZAP  C5   C6   aromatic  1.398  0.020
 ZAP  C6   C1   aromatic  1.399  0.020
 ZAP  N1   HN11 single    1.015  0.020
 ZAP  N1   HN12 single    1.015  0.020
 ZAP  C2   HC21 single    1.086  0.020
 ZAP  C3   HC31 single    1.086  0.020
 ZAP  O1   HO11 single    0.971  0.020
 ZAP  C5   HC51 single    1.085  0.020
 ZAP  C6   HC61 single    1.087  0.020
''')

# arginine zwitterion with eLBOW-style energy types (N NT3, NH1/NH2 NC2, NE NC1,
# O/OXT OC)
zar_cif = ('zar.cif', '''
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 ZAR  ZAR  'ZAR' ligand 27 12 .

data_comp_ZAR
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
 ZAR  N    N  NT3    1  0.000  -3.1231   0.1902   1.4525
 ZAR  CA   C  C      0  0.000  -1.9573  -0.0227   0.4936
 ZAR  CB   C  C      0  0.000  -0.8550   0.9922   0.8122
 ZAR  CG   C  C      0  0.000  -0.1736   1.6637  -0.3997
 ZAR  CD   C  C      0  0.000   0.3905   0.7596  -1.4947
 ZAR  NE   N  NC1    0  0.000   1.1320  -0.3687  -0.9509
 ZAR  CZ   C  C      0  0.000   2.4255  -0.6406  -0.9079
 ZAR  NH1  N  NC2    0  0.000   2.8163  -1.7886  -0.3504
 ZAR  NH2  N  NC2    1  0.000   3.3474   0.1926  -1.3916
 ZAR  C    C  C      0  0.000  -1.5107  -1.5248   0.6290
 ZAR  O    O  OC     0  0.000  -0.3802  -1.8453   0.1776
 ZAR  OXT  O  OC    -1  0.000  -2.4029  -2.1844   1.2423
 ZAR  HN1  H  H      0  0.000  -3.4611  -0.7997   1.5810
 ZAR  HN2  H  H      0  0.000  -3.8992   0.7494   1.0961
 ZAR  HN3  H  H      0  0.000  -2.8351   0.4615   2.3967
 ZAR  HCA1 H  H      0  0.000  -2.3826   0.1273  -0.5039
 ZAR  HCB1 H  H      0  0.000  -0.0789   0.5184   1.4267
 ZAR  HCB2 H  H      0  0.000  -1.2601   1.8077   1.4261
 ZAR  HCG1 H  H      0  0.000   0.6483   2.2775  -0.0090
 ZAR  HCG2 H  H      0  0.000  -0.8825   2.3618  -0.8619
 ZAR  HCD1 H  H      0  0.000   1.0175   1.3289  -2.1874
 ZAR  HCD2 H  H      0  0.000  -0.4297   0.3438  -2.0903
 ZAR  HNE1 H  H      0  0.000   0.5550  -1.1445  -0.5736
 ZAR  HH11 H  H      0  0.000   3.7813  -2.0738  -0.2862
 ZAR  HH12 H  H      0  0.000   2.1049  -2.4126   0.0327
 ZAR  HH21 H  H      0  0.000   3.0828   1.0678  -1.8186
 ZAR  HH22 H  H      0  0.000   4.3304  -0.0368  -1.3544

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
 ZAR  N    CA   single    1.524  0.020
 ZAR  CA   CB   single    1.532  0.020
 ZAR  CB   CG   single    1.544  0.020
 ZAR  CG   CD   single    1.528  0.020
 ZAR  CD   NE   single    1.456  0.020
 ZAR  NE   CZ   single    1.322  0.020
 ZAR  CZ   NH1  single    1.335  0.020
 ZAR  CZ   NH2  double    1.333  0.020
 ZAR  CA   C    single    1.573  0.020
 ZAR  C    O    double    1.259  0.020
 ZAR  C    OXT  single    1.268  0.020
 ZAR  N    HN1  single    1.054  0.020
 ZAR  N    HN2  single    1.021  0.020
 ZAR  N    HN3  single    1.024  0.020
 ZAR  CA   HCA1 single    1.095  0.020
 ZAR  CB   HCB1 single    1.097  0.020
 ZAR  CB   HCB2 single    1.098  0.020
 ZAR  CG   HCG1 single    1.098  0.020
 ZAR  CG   HCG2 single    1.097  0.020
 ZAR  CD   HCD1 single    1.094  0.020
 ZAR  CD   HCD2 single    1.096  0.020
 ZAR  NE   HNE1 single    1.038  0.020
 ZAR  NH1  HH11 single    1.008  0.020
 ZAR  NH1  HH12 single    1.021  0.020
 ZAR  NH2  HH21 single    1.009  0.020
 ZAR  NH2  HH22 single    1.010  0.020
''')

# neutral methylguanidine, 5-methyl-1H-tetrazole, methanesulfonic acid,
# N,N-dimethylglycine N-methylamide (tertiary amine and amide N)
zgn_cif = ('zgn.cif', '''
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 ZGN  ZGN  'ZGN' ligand 12 5 .

data_comp_ZGN
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
 ZGN  C1   C  C      0  0.000  -1.7473  -0.2767  -0.1645
 ZGN  N1   N  N      0  0.000  -0.3641  -0.6448   0.0622
 ZGN  C2   C  C      0  0.000   0.6058   0.2979   0.0462
 ZGN  N2   N  N      0  0.000   1.8307  -0.2590   0.1177
 ZGN  N3   N  N      0  0.000   0.4057   1.5613  -0.0157
 ZGN  HC11 H  H      0  0.000  -2.3830  -1.1608  -0.0573
 ZGN  HC12 H  H      0  0.000  -1.8831   0.1223  -1.1747
 ZGN  HC13 H  H      0  0.000  -2.0813   0.4677   0.5652
 ZGN  HN11 H  H      0  0.000  -0.0970  -1.4818  -0.4416
 ZGN  HN21 H  H      0  0.000   2.5683   0.3875   0.3656
 ZGN  HN22 H  H      0  0.000   1.8462  -1.0665   0.7302
 ZGN  HN31 H  H      0  0.000   1.2990   2.0529  -0.0334

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
 ZGN  C1   N1   single    1.449  0.020
 ZGN  N1   C2   single    1.353  0.020
 ZGN  C2   N2   single    1.347  0.020
 ZGN  C2   N3   double    1.281  0.020
 ZGN  C1   HC11 single    1.094  0.020
 ZGN  C1   HC12 single    1.095  0.020
 ZGN  C1   HC13 single    1.095  0.020
 ZGN  N1   HN11 single    1.013  0.020
 ZGN  N2   HN21 single    1.012  0.020
 ZGN  N2   HN22 single    1.014  0.020
 ZGN  N3   HN31 single    1.020  0.020
''')

ztz_cif = ('ztz.cif', '''
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 ZTZ  ZTZ  'ZTZ' ligand 10 6 .

data_comp_ZTZ
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
 ZTZ  C1   C  C      0  0.000  -1.4097  -0.0368   0.0429
 ZTZ  C2   C  C      0  0.000   0.0443  -0.2591  -0.0411
 ZTZ  N1   N  N      0  0.000   0.6896  -1.3889  -0.2361
 ZTZ  N2   N  N      0  0.000   2.0265  -1.0692  -0.2332
 ZTZ  N3   N  N      0  0.000   2.2041   0.2303  -0.0406
 ZTZ  N4   N  N      0  0.000   0.9708   0.7362   0.0792
 ZTZ  HC11 H  H      0  0.000  -1.9523  -0.9788  -0.0825
 ZTZ  HC12 H  H      0  0.000  -1.7415   0.6529  -0.7393
 ZTZ  HC13 H  H      0  0.000  -1.6812   0.3846   1.0157
 ZTZ  HN41 H  H      0  0.000   0.8494   1.7287   0.2351

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
 ZTZ  C1   C2   single    1.473  0.020
 ZTZ  C2   N1   aromatic  1.316  0.020
 ZTZ  N1   N2   aromatic  1.375  0.020
 ZTZ  N2   N3   aromatic  1.326  0.020
 ZTZ  N3   N4   aromatic  1.338  0.020
 ZTZ  N4   C2   aromatic  1.365  0.020
 ZTZ  C1   HC11 single    1.094  0.020
 ZTZ  C1   HC12 single    1.094  0.020
 ZTZ  C1   HC13 single    1.094  0.020
 ZTZ  N4   HN41 single    1.012  0.020
''')

zsa_cif = ('zsa.cif', '''
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 ZSA  ZSA  'ZSA' ligand 9 5 .

data_comp_ZSA
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
 ZSA  C1   C  C      0  0.000  -1.1285  -0.0423  -0.0372
 ZSA  S1   S  S      0  0.000   0.5990  -0.3692   0.2262
 ZSA  O1   O  O      0  0.000   1.1211  -1.1321  -0.8844
 ZSA  O2   O  O      0  0.000   0.8185  -0.7449   1.6024
 ZSA  O3   O  O      0  0.000   1.1458   1.1355   0.0207
 ZSA  HC11 H  H      0  0.000  -1.6681  -0.9877   0.0494
 ZSA  HC12 H  H      0  0.000  -1.4748   0.6575   0.7255
 ZSA  HC13 H  H      0  0.000  -1.2606   0.3738  -1.0376
 ZSA  HO31 H  H      0  0.000   1.8476   1.1094  -0.6652

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
 ZSA  C1   S1   single    1.778  0.020
 ZSA  S1   O1   double    1.445  0.020
 ZSA  S1   O2   double    1.443  0.020
 ZSA  S1   O3   single    1.614  0.020
 ZSA  C1   HC11 single    1.092  0.020
 ZSA  C1   HC12 single    1.092  0.020
 ZSA  C1   HC13 single    1.092  0.020
 ZSA  O3   HO31 single    0.982  0.020
''')

zam_cif = ('zam.cif', '''
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 ZAM  ZAM  'ZAM' ligand 20 8 .

data_comp_ZAM
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
 ZAM  C1   C  C      0  0.000  -2.6433  -0.8280   0.2400
 ZAM  N1   N  N      0  0.000  -1.4282  -0.0146   0.2777
 ZAM  C2   C  C      0  0.000  -1.3787   0.8755  -0.8812
 ZAM  C3   C  C      0  0.000  -0.2393  -0.8881   0.3629
 ZAM  C4   C  C      0  0.000   1.0389  -0.0872   0.6476
 ZAM  O1   O  O      0  0.000   1.3281   0.3403   1.7625
 ZAM  N2   N  N      0  0.000   1.8397   0.1122  -0.4629
 ZAM  C5   C  C      0  0.000   3.1035   0.7939  -0.3533
 ZAM  HC11 H  H      0  0.000  -2.7171  -1.4540   1.1363
 ZAM  HC12 H  H      0  0.000  -3.5347  -0.1908   0.2316
 ZAM  HC13 H  H      0  0.000  -2.6766  -1.4799  -0.6405
 ZAM  HC21 H  H      0  0.000  -1.2894   0.3224  -1.8230
 ZAM  HC22 H  H      0  0.000  -2.2838   1.4918  -0.9316
 ZAM  HC23 H  H      0  0.000  -0.5412   1.5767  -0.8076
 ZAM  HC31 H  H      0  0.000  -0.3429  -1.5931   1.1976
 ZAM  HC32 H  H      0  0.000  -0.1174  -1.4883  -0.5474
 ZAM  HN21 H  H      0  0.000   1.6160  -0.3640  -1.3258
 ZAM  HC51 H  H      0  0.000   3.8417   0.0979   0.0538
 ZAM  HC52 H  H      0  0.000   3.4113   1.1186  -1.3498
 ZAM  HC53 H  H      0  0.000   3.0136   1.6588   0.3094

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
 ZAM  C1   N1   single    1.463  0.020
 ZAM  N1   C2   single    1.462  0.020
 ZAM  N1   C3   single    1.478  0.020
 ZAM  C3   C4   single    1.535  0.020
 ZAM  C4   O1   double    1.229  0.020
 ZAM  C4   N2   single    1.384  0.020
 ZAM  N2   C5   single    1.440  0.020
 ZAM  C1   HC11 single    1.096  0.020
 ZAM  C1   HC12 single    1.096  0.020
 ZAM  C1   HC13 single    1.096  0.020
 ZAM  C2   HC21 single    1.096  0.020
 ZAM  C2   HC22 single    1.096  0.020
 ZAM  C2   HC23 single    1.095  0.020
 ZAM  C3   HC31 single    1.098  0.020
 ZAM  C3   HC32 single    1.097  0.020
 ZAM  N2   HN21 single    1.011  0.020
 ZAM  C5   HC51 single    1.093  0.020
 ZAM  C5   HC52 single    1.092  0.020
 ZAM  C5   HC53 single    1.093  0.020
''')

# acetohydrazide, N-methylhydroxylamine, dimethylcyanamide, 1-methyltetrazole
zhz_cif = ('zhz.cif', '''
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 ZHZ  ZHZ  'ZHZ' ligand 11 5 .

data_comp_ZHZ
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
 ZHZ  C1   C  C      0  0.000  -1.6095   0.0569   0.3409
 ZHZ  C2   C  C      0  0.000  -0.1439   0.3905   0.2836
 ZHZ  O1   O  O      0  0.000   0.3275   1.3922   0.8171
 ZHZ  N1   N  N      0  0.000   0.6164  -0.5225  -0.4325
 ZHZ  N2   N  N      0  0.000   1.9895  -0.2661  -0.6381
 ZHZ  HC11 H  H      0  0.000  -1.9524   0.1310   1.3766
 ZHZ  HC12 H  H      0  0.000  -2.1618   0.7642  -0.2831
 ZHZ  HC13 H  H      0  0.000  -1.8067  -0.9587  -0.0134
 ZHZ  HN11 H  H      0  0.000   0.1686  -1.2465  -0.9907
 ZHZ  HN21 H  H      0  0.000   2.4751  -0.4969   0.2350
 ZHZ  HN22 H  H      0  0.000   2.0971   0.7559  -0.6954

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
 ZHZ  C1   C2   single    1.504  0.020
 ZHZ  C2   O1   double    1.229  0.020
 ZHZ  C2   N1   single    1.387  0.020
 ZHZ  N1   N2   single    1.412  0.020
 ZHZ  C1   HC11 single    1.093  0.020
 ZHZ  C1   HC12 single    1.093  0.020
 ZHZ  C1   HC13 single    1.093  0.020
 ZHZ  N1   HN11 single    1.018  0.020
 ZHZ  N2   HN21 single    1.025  0.020
 ZHZ  N2   HN22 single    1.029  0.020
''')

zha_cif = ('zha.cif', '''
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 ZHA  ZHA  'ZHA' ligand 8 3 .

data_comp_ZHA
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
 ZHA  C1   C  C      0  0.000  -0.4928  -0.0423  -0.6962
 ZHA  N1   N  N      0  0.000  -0.0967   0.2973   0.6659
 ZHA  O1   O  O      0  0.000   1.2119  -0.3142   0.8381
 ZHA  HC11 H  H      0  0.000  -1.4437   0.4439  -0.9333
 ZHA  HC12 H  H      0  0.000  -0.6325  -1.1229  -0.8001
 ZHA  HC13 H  H      0  0.000   0.2526   0.2929  -1.4248
 ZHA  HN11 H  H      0  0.000   0.1253   1.2958   0.7067
 ZHA  HO11 H  H      0  0.000   1.0759  -0.8505   1.6438

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
 ZHA  C1   N1   single    1.459  0.020
 ZHA  N1   O1   single    1.455  0.020
 ZHA  C1   HC11 single    1.094  0.020
 ZHA  C1   HC12 single    1.094  0.020
 ZHA  C1   HC13 single    1.095  0.020
 ZHA  N1   HN11 single    1.024  0.020
 ZHA  O1   HO11 single    0.977  0.020
''')

zcy_cif = ('zcy.cif', '''
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 ZCY  ZCY  'ZCY' ligand 11 5 .

data_comp_ZCY
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
 ZCY  C1   C  C      0  0.000  -1.1558   0.6061  -0.1188
 ZCY  N1   N  N      0  0.000   0.2818   0.2467  -0.1800
 ZCY  C2   C  C      0  0.000   0.4950  -1.2130  -0.0272
 ZCY  C3   C  C      0  0.000   1.2364   1.1299   0.1572
 ZCY  N2   N  N      0  0.000   2.0562   1.8890   0.4589
 ZCY  HC11 H  H      0  0.000  -1.7388  -0.0232  -0.7982
 ZCY  HC12 H  H      0  0.000  -1.5330   0.4787   0.9007
 ZCY  HC13 H  H      0  0.000  -1.3007   1.6486  -0.4203
 ZCY  HC21 H  H      0  0.000   0.2804  -1.5197   1.0013
 ZCY  HC22 H  H      0  0.000  -0.1536  -1.7700  -0.7103
 ZCY  HC23 H  H      0  0.000   1.5321  -1.4731  -0.2632

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
 ZCY  C1   N1   single    1.483  0.020
 ZCY  N1   C2   single    1.483  0.020
 ZCY  N1   C3   single    1.344  0.020
 ZCY  C3   N2   triple    1.157  0.020
 ZCY  C1   HC11 single    1.094  0.020
 ZCY  C1   HC12 single    1.094  0.020
 ZCY  C1   HC13 single    1.095  0.020
 ZCY  C2   HC21 single    1.094  0.020
 ZCY  C2   HC22 single    1.094  0.020
 ZCY  C2   HC23 single    1.095  0.020
''')

zmt_cif = ('zmt.cif', '''
data_comp_list
loop_
_chem_comp.id
_chem_comp.three_letter_code
_chem_comp.name
_chem_comp.group
_chem_comp.number_atoms_all
_chem_comp.number_atoms_nh
_chem_comp.desc_level
 ZMT  ZMT  'ZMT' ligand 10 6 .

data_comp_ZMT
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.charge
_chem_comp_atom.partial_charge
_chem_comp_atom.x
_chem_comp_atom.y
_chem_comp_atom.z
 ZMT  C1   C  C      0  0.000  -1.3779  -0.0728  -0.0198
 ZMT  N1   N  N      0  0.000   0.0507  -0.1820   0.0317
 ZMT  C2   C  C      0  0.000   0.9863   0.7893  -0.1094
 ZMT  N2   N  N      0  0.000   2.1730   0.2395   0.0103
 ZMT  N3   N  N      0  0.000   1.9452  -1.0999   0.2296
 ZMT  N4   N  N      0  0.000   0.6450  -1.3645   0.2438
 ZMT  HC11 H  H      0  0.000  -1.7468  -0.7027  -0.8331
 ZMT  HC12 H  H      0  0.000  -1.6576   0.9680  -0.2008
 ZMT  HC13 H  H      0  0.000  -1.7879  -0.4059   0.9369
 ZMT  HC21 H  H      0  0.000   0.7700   1.8312  -0.2891

loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
 ZMT  C1   N1   single    1.434  0.020
 ZMT  N1   C2   aromatic  1.356  0.020
 ZMT  C2   N2   aromatic  1.313  0.020
 ZMT  N2   N3   aromatic  1.376  0.020
 ZMT  N3   N4   aromatic  1.327  0.020
 ZMT  N4   N1   aromatic  1.340  0.020
 ZMT  C1   HC11 single    1.093  0.020
 ZMT  C1   HC12 single    1.093  0.020
 ZMT  C1   HC13 single    1.093  0.020
 ZMT  C2   HC21 single    1.079  0.020
''')

# ------------------------------------------------------------------------------

def get_model(lines=None, cifs=()):
  '''cifs: (file name, restraint cif text) supplied as restraint objects.'''
  if lines is None:
    lines = model_str.split('\n')
  model = mmtbx.model.manager(
    model_input=iotbx.pdb.input(lines=lines, source_info=None),
    restraint_objects=[(n, iotbx.cif.reader(input_string=t).model()) for n, t in cifs] or None,
    log=null_out())
  model.process(make_restraints=True)
  model.set_hydrogen_bond_length(use_neutron_distances=False)
  return model

def get_manager(model, sel=LIG_SEL, **kwargs):
  params = LI.master_params().extract().ligand_interactions
  for k, v in kwargs.items():
    setattr(params, k, v)
  isel = model.selection(sel).iselection()
  return LI.manager(model, isel, sel, params=params).run()

def short(label):
  '''"A EDO 3 H12" -> "EDO 3 H12"'''
  return None if label is None else label.split(None, 1)[1]

def entries_by_type(m):
  result = {}
  for e in m.entries:
    result.setdefault(e['type'], []).append(e)
  return result

# ------------------------------------------------------------------------------

def exercise_entries_and_counts(model):
  '''
  One H-bond and one clash (pnp and probe2 agree), vdW contacts with probe's
  class as subtype; counts; criteria recorded; as_dict is JSON-able; show() runs.
  '''
  m = get_manager(model)
  by = entries_by_type(m)
  assert len(by['hbond']) == 1, by['hbond']
  hb = by['hbond'][0]
  assert [short(l) for l in hb['labels']] == ['EDO 1 O1', 'EDO 1 HO1', 'EDO 2 O1']
  assert hb['sources'] == ['pnp', 'probe2']
  assert hb['cross_check'] == 'pnp and probe2'
  assert hb['residue'] == 'A EDO 2' and hb['symop'] is None
  assert approx_equal(hb['geometry']['pnp']['d_HA'], 2.011, eps=0.01)
  assert hb['geometry']['probe2']['dots']['hb'] > 0
  assert len(by['clash']) == 1, by['clash']
  cl = by['clash'][0]
  assert [short(l) for l in cl['labels']] == ['EDO 1 H11', 'EDO 3 H12']
  assert cl['sources'] == ['pnp', 'probe2']
  assert approx_equal(cl['geometry']['pnp']['overlap'], 0.69, eps=0.01)
  assert approx_equal(cl['geometry']['probe2']['min_gap'], -0.69, eps=0.01)
  vdw = dict([(tuple([short(l) for l in e['labels']]), e) for e in by['vdw']])
  assert vdw[('EDO 1 HO2', 'EDO 4 H12')]['subtype'] == 'so'
  assert vdw[('EDO 1 HO2', 'EDO 4 C1')]['subtype'] == 'cc'
  assert vdw[('EDO 1 O2', 'EDO 4 C1')]['subtype'] == 'wc'
  for e in by['vdw']:
    assert e['sources'] == ['probe2'] and e['cross_check'] is None
    assert e['geometry']['probe2']['min_gap'] is not None
  for e in m.entries:
    assert e['model_support'] is None
  assert m.disagreements == [] and m.internal == []
  assert m.probe_unmapped == set()
  c = m.counts()
  assert c.per_type['hbond'] == 1 and c.per_type['clash'] == 1
  assert c.per_type['vdw:so'] == 1 and c.per_type['vdw:cc'] == 1
  assert sum(c.per_type.values()) == len(m.entries)
  assert c.per_residue['A EDO 3'] == {'clash': 1}
  assert c.per_ligand_atom['A EDO 1 H11'] == {'clash': 1}
  d = json.loads(json.dumps(m.as_dict()))
  assert d['pair_class_order'] == ['bo', 'hb', 'so', 'cc', 'wc']
  assert d['probe']['density'] == 16
  assert d['hbond_criteria']['d_HA_cutoff'] == [1.4, 2.8]
  assert d['hbond_criteria']['d_DA_cutoff'] == [2.4, 4.1]
  assert d['hbond_criteria']['a_DHA_cutoff'] == 120
  assert d['hbond_criteria']['min_bonds_H_A'] == 5
  assert d['clash_criteria']['min_overlap'] == 0.4
  assert len(d['entries']) == len(m.entries)
  log = StringIO()
  m.show(log=log)
  assert 'contact patches' in log.getvalue()
  assert 'd_HA_cutoff=[1.4, 2.8]' in log.getvalue()
  return m

def exercise_pair_class_order(model):
  '''
  Plumbing: pair_class_order changes the label (the clash pair has so and bo
  dots; with so first it becomes a vdW contact and pnp's clash a disagreement)
  and invalid orders are refused. Not a recommended setting.
  '''
  m = get_manager(model, pair_class_order=['so', 'bo', 'hb', 'cc', 'wc'])
  by = entries_by_type(m)
  cl = by['clash'][0]
  assert cl['sources'] == ['pnp'] and cl['cross_check'] == 'pnp only', cl
  pair = [e for e in by['vdw'] if [short(l) for l in e['labels']] ==
    ['EDO 1 H11', 'EDO 3 H12']]
  assert len(pair) == 1 and pair[0]['subtype'] == 'so'
  assert pair[0]['geometry']['probe2']['dots']['bo'] > 0
  assert len(m.disagreements) == 1
  d = m.disagreements[0]
  assert d['type'] == 'clash' and d['missing'] == ['probe2']
  assert d['probe_class'] == 'so'
  # each class once
  for bad in (['bo', 'hb', 'so', 'cc'], ['bo', 'hb', 'so', 'cc', 'wc', 'wc'],
              ['bo', 'hb', 'so', 'cc', 'xx']):
    try:
      get_manager(model, pair_class_order=bad)
    except Sorry:
      pass
    else:
      raise AssertionError('pair_class_order %s accepted' % bad)

def exercise_internal():
  '''
  An intramolecular clash (folded EDO, HO1...HO2): listed as internal, not an
  entry, not a disagreement.
  '''
  m = get_manager(get_model(fold_model_str.split('\n')), sel='chain A and resseq 1')
  assert len(m.overlaps.clash_records) == 1
  assert [(d['type'], [short(l) for l in d['labels']]) for d in m.internal] == \
    [('clash', ['EDO 1 HO1', 'EDO 1 HO2'])], m.internal
  assert approx_equal(m.internal[0]['geometry']['pnp']['overlap'], 0.49, eps=0.01)
  assert m.entries == [] and m.disagreements == []
  assert len(json.loads(json.dumps(m.as_dict()))['internal']) == 1

def exercise_symmetry():
  '''
  H-bonds between the ligand and its symmetry-related copies: kept with the
  operator, "symmetry, probe2 not applicable", not disagreements, not internal.
  '''
  model = get_model(sym_model_str.split('\n'))
  m = get_manager(model, sel='chain A and resseq 1')
  by = entries_by_type(m)
  assert len(by['hbond']) == 2, m.entries
  # pnp's operator belongs to its pair: x+1,y,z for both records
  for r in m.overlaps.hbond_records:
    assert r['symop'] == 'x+1,y,z'
  uc = model.crystal_symmetry().unit_cell()
  atoms = model.get_hierarchy().atoms()
  def site(i, op):
    return uc.orthogonalize(sgtbx.rt_mx(op) * uc.fractionalize(atoms[i].xyz))
  ops = []
  for e in by['hbond']:
    assert e['cross_check'] == 'symmetry, probe2 not applicable'
    assert e['sources'] == ['pnp']
    d, h, a = e['atoms']
    # the ligand donor stays, the acceptor (its own copy) moves
    assert e['operators'][:2] == ['x,y,z', 'x,y,z'] and e['operators'][2] == e['symop']
    d_HA = e['geometry']['pnp']['d_HA']
    assert approx_equal(d_HA, 2.62, eps=0.01)
    assert approx_equal(atoms[h].distance(site(a, e['operators'][2])), d_HA, eps=1.e-6)
    # pnp's own operator on the acceptor gives the other copy
    assert approx_equal(atoms[h].distance(site(a, 'x+1,y,z')) if e['symop'] != 'x+1,y,z'
      else atoms[h].distance(site(a, 'x-1,y,z')), 9.55, eps=0.01)
    assert e['residue'] == 'A EDO 1 (%s)' % e['symop']
    assert e['labels'][2].endswith(' (%s)' % e['symop'])
    assert e['ligand_atoms'] == [e['labels'][0], e['labels'][1]]
    ops.append(e['symop'])
  assert sorted(ops) == ['x+1,y,z', 'x-1,y,z'], ops
  assert m.disagreements == [] and m.internal == []
  per_residue = m.counts().per_residue
  assert sorted(per_residue) == ['A EDO 1 (x+1,y,z)', 'A EDO 1 (x-1,y,z)'], per_residue
  assert '-' not in per_residue and None not in per_residue

def exercise_donor_conformers():
  '''
  A probe2-only H-bond from a donor with alternate conformers. probe2's hb class
  is overlap-based and ignores the angle, which is why only probe2 reports this
  H-bond (D-H...A 113 deg, below pnp's 120; recorded from the model). The donor
  is the H's bonded O1 of the same conformer (restraints' connectivity), not the
  nearest heavy atom (conformer B's O1); the partner residue is right with the
  donor, the acceptor or one donor conformer selected as the ligand.
  '''
  model = get_model(alt_model_str.split('\n'))
  atoms = model.get_hierarchy().atoms()
  by_name = dict([((a.parent().parent().resseq.strip(), a.name.strip(),
    a.parent().altloc), a) for a in atoms])
  h_a = by_name[('1', 'HO1', 'A')]
  o_a, o_b = by_name[('1', 'O1', 'A')], by_name[('1', 'O1', 'B')]
  assert h_a.distance(o_b) < h_a.distance(o_a) < 1.3  # the old lookup's choice
  for sel, residue in (('chain A and resseq 2', 'A EDO 1'),
                       ('chain A and resseq 1', 'A EDO 2'),
                       ('chain A and resseq 1 and altloc A', 'A EDO 2')):
    m = get_manager(model, sel=sel)
    assert m.overlaps.hbond_records == [] and m.unresolved_donors == []
    hb = [e for e in m.entries if e['type'] == 'hbond' and e['atoms'][1] == h_a.i_seq]
    assert len(hb) == 1, (sel, m.entries)
    e = hb[0]
    assert e['atoms'][0] == o_a.i_seq, (sel, e['labels'])
    assert e['labels'] == ['A EDO 1 O1 alt A', 'A EDO 1 HO1 alt A', 'A EDO 2 O1'], e
    assert e['residue'] == residue, (sel, e['residue'])
    assert e['cross_check'] == 'probe2 only'
    assert e['operators'] == ['x,y,z'] * 3 and e['symop'] is None
    g = e['geometry']['probe2']
    assert approx_equal(g['a_DHA'], h_a.angle(by_name[('2', 'O1', '')], o_a, deg=True),
      eps=1.e-6)
    assert 110 < g['a_DHA'] < 120, g['a_DHA']
    assert approx_equal(g['d_HA'], 2.28, eps=0.02), g['d_HA']
    n_b = len([x for x in m.entries if x['type'] == 'hbond'])
    assert n_b == (1 if 'altloc A' in sel else 2), n_b
  # an H without a bonded heavy atom: reported, donor unset
  m._fsc0 = [[] for a in atoms]
  assert m._donor_of(h_a.i_seq) is None
  assert m.unresolved_donors == [dict(h='A EDO 1 HO1 alt A', bonded_heavy=[])]

def exercise_probe_names():
  '''
  probe2's atom field (its own format) for names that run together; the key
  parsed by columns equals atom_key, the table maps each field to its atom.
  '''
  h = iotbx.pdb.input(lines=names_str.split('\n'), source_info=None).construct_hierarchy()
  atoms = h.atoms()
  fields = [LI.probe_atom_field(a) for a in atoms]
  by_name = dict([((a.name.strip(), a.parent().altloc), f) for a, f in zip(atoms, fields)])
  assert by_name[('HD21', 'A')] == ' A1001BASN HD21A', by_name
  assert by_name[('HD21', 'B')] == ' A1001BASN HD21B'
  assert by_name[('CA', '')] == ' B  52AGLY  CA  '
  assert by_name[('HD21', 'A')].split() == ['A1001BASN', 'HD21A']  # why not split
  for a, f in zip(atoms, fields):
    assert LI.probe_key(f) == LI.atom_key(a), (f, LI.probe_key(f), LI.atom_key(a))
  assert LI.probe_key(by_name[('HD21', 'A')]) == ('A', '1001B', 'ASN', 'HD21', 'A')
  assert LI.probe_key(by_name[('CA', '')]) == ('B', '52A', 'GLY', 'CA', '')
  table = LI.probe_atom_table(atoms)
  assert [table[f] for f in fields] == [a.i_seq for a in atoms]

def exercise_probe_names_in_run():
  '''
  The same on a probe2 run: ligand A1001, H-bond partner 52A (insertion code),
  vdW partner EDO 4 in two altlocs; every field maps to the atom it names.
  '''
  lines = []
  for l in model_str.split('\n'):
    if l[17:26] == 'EDO A   1':
      l = l[:22] + '1001' + l[26:]
    elif l[17:26] == 'EDO A   2':
      l = l[:22] + '  52A' + l[27:]
    elif l[17:26] == 'EDO A   4':
      lb = l[:16] + 'B' + l[17:46] + '%8.3f' % (float(l[46:54]) - 0.3) + l[54:]
      l = l[:16] + 'A' + l[17:54] + '  0.50' + l[60:]
      lb = lb[:54] + '  0.50' + lb[60:]
      lines.append(l)
      l = lb
    lines.append(l)
  model = get_model(lines)
  m = get_manager(model, sel='chain A and resseq 1001')
  atoms = model.get_hierarchy().atoms()
  assert m.probe_unmapped == set()
  table = LI.probe_atom_table(atoms)
  fields = set()
  for line in m.probe_output.splitlines():
    f = line.split(':')
    if len(f) > 4 and f[2] in LI.probe_classes:
      fields.update([f[3], f[4]])
  for f in fields:
    assert LI.probe_key(f) == LI.atom_key(atoms[table[f]]), f
  assert any([f.startswith(' A1001 EDO') for f in fields])
  assert any([f.startswith(' A  52AEDO') for f in fields])
  alts = set([f[-1] for f in fields if f.startswith(' A   4 EDO')])
  assert alts == set(['A', 'B']), alts
  hb = entries_by_type(m)['hbond']
  assert [short(l) for l in hb[0]['labels']] == ['EDO 1001 O1', 'EDO 1001 HO1',
    'EDO 52A O1'], hb
  assert hb[0]['cross_check'] == 'pnp and probe2'
  vdw_alts = set([l.split()[-1] for e in m.entries if e['type'] == 'vdw'
    for l in e['labels'] if ' EDO 4 ' in l])
  assert vdw_alts == set(['A', 'B']), vdw_alts

def exercise_dot_patches():
  '''Two clusters 2 A apart; a chain of points 0.4 A apart links them at 0.5 A only.'''
  a = [(x * 0.25, y * 0.25, 0) for x in range(4) for y in range(4)]
  b = [(p[0] + 2.75, p[1], p[2]) for p in a]
  groups = LI.dot_patches(flex.vec3_double(a + b), 0.5)
  assert [len(g) for g in groups] == [16, 16]
  assert sorted(groups[0] + groups[1]) == list(range(32))
  chain = [(0.75 + 0.4 * k, 0, 0) for k in range(1, 5)]
  groups = LI.dot_patches(flex.vec3_double(a + b + chain), 0.5)
  assert [len(g) for g in groups] == [36]
  groups = LI.dot_patches(flex.vec3_double(a + b + chain), 0.35)
  assert [len(g) for g in groups] == [16, 16, 1, 1, 1, 1]
  assert LI.dot_patches(flex.vec3_double(), 0.5) == []

def exercise_patches(m):
  '''
  Patches on the model: each patch sees one environment copy here. Dots are
  at their location on the ligand atom's surface (vdW radius from the atom), not
  at the spike end, which overlap dots push into the neighbour.
  '''
  lig = set(m.ligand_isel)
  out = [d for d in m.probe_dots if d['source'] in lig and d['target'] not in lig]
  assert sum([p['n_dots'] for p in m.patches]) == len(out)
  residues = [p['residues'] for p in m.patches]
  assert ['A EDO 3'] in residues and ['A EDO 4'] in residues
  big = [p for p in m.patches if p['residues'] == ['A EDO 3']][0]
  assert big['ligand_atoms'] == ['A EDO 1 H11']
  assert big['dots']['bo'] > 0 and big['dots']['so'] > 0
  assert sum(big['dots'].values()) == big['n_dots']
  assert sum([q['dots'] for q in big['pairs']]) == big['n_dots']
  atoms = m.model.get_hierarchy().atoms()
  h11 = [d for d in out if atoms[d['source']].name.strip() == 'H11']
  r_loc = [atoms[d['source']].distance(d['loc']) for d in h11]
  r_spike = [atoms[d['source']].distance(d['spike']) for d in h11]
  assert max(r_loc) - min(r_loc) < 0.01, (min(r_loc), max(r_loc))
  assert max(r_spike) - min(r_spike) > 0.1
  xyz = flex.vec3_double([d['loc'] for d in out])
  groups = LI.dot_patches(xyz, 0.5)
  assert [len(g) for g in groups] == [p['n_dots'] for p in m.patches]
  for g, p in zip(groups, m.patches):
    assert approx_equal(xyz.select(flex.size_t(g)).mean(), p['center'], eps=1.e-9)

def exercise_library_vs_command_line(model):
  '''
  probe2 as a library call gives the same dots as mmtbx.probe2 on the file (both
  from the file: the PDB format rounds coordinates to 0.001 A).
  '''
  src, tgt = '(%s)' % LIG_SEL, 'not (%s)' % LIG_SEL
  fn = 'tst_ligand_interactions.pdb'
  with open(fn, 'w') as f:
    f.write(model.model_as_pdb())
  from_file = mmtbx.model.manager(model_input=iotbx.pdb.input(file_name=fn),
    log=null_out())
  from_file.process(make_restraints=True)
  text = LI.run_probe2(from_file, src, tgt)
  out = 'tst_ligand_interactions_probe2.txt'
  cmd = ('mmtbx.probe2 %s approach=both "source_selection=%s" "target_selection=%s"'
    ' output.format=raw output.filename=%s --overwrite' % (fn, src, tgt, out))
  r = easy_run.fully_buffered(cmd)
  assert r.return_code == 0, '\n'.join(r.stderr_lines)
  with open(out) as f:
    cli = f.read()
  lib_lines = [l for l in text.splitlines() if l.strip()]
  cli_lines = [l for l in cli.splitlines() if l.strip()]
  assert len(lib_lines) > 100
  assert lib_lines == cli_lines

def exercise_ligand_overlaps(model):
  '''
  ligand_overlaps is what validate_ligands' get_overlaps reports; its records
  point to the right atoms of the full model. EDO 0 comes first in the file but
  lies outside the 3 A region, so model_within's atom numbers differ from the
  full model's by 10; the assertions fail if the mapping back to the full model
  is missing.
  '''
  from mmtbx.validation import validate_ligands
  isel = model.selection(LIG_SEL).iselection()
  ov = LI.ligand_overlaps(model, LIG_SEL)
  params = validate_ligands.master_params().extract().validate_ligands
  lr = validate_ligands.ligand_result(model=model, fmodel=None, map_manager=None,
    ligand_isel=isel, sel_str=LIG_SEL, params=params)
  vl = lr.get_overlaps()
  for k in ('n_clashes', 'clashscore', 'n_clashes_sym', 'clashes_str', 'n_hbonds',
            'clash_records', 'hbond_records'):
    assert getattr(vl, k) == getattr(ov, k), k
  atoms = model.get_hierarchy().atoms()
  assert ov.n_clashes == len(ov.clash_records) == 1
  assert ov.n_hbonds == len(ov.hbond_records) == 1
  r = ov.clash_records[0]
  assert r['symop'] == ''
  names = [atoms[r['i_seq']].id_str(), atoms[r['j_seq']].id_str()]
  assert 'H11 EDO A   1' in names[0] and 'H12 EDO A   3' in names[1], names
  assert approx_equal(atoms[r['i_seq']].distance(atoms[r['j_seq']]), r['distance'],
    eps=1.e-6)
  assert '%.2f' % r['distance'] in ov.clashes_str
  assert '%.2f' % r['overlap'] in ov.clashes_str
  for field in ('H11 EDO A   1', 'H12 EDO A   3'):
    assert field in ov.clashes_str
  h = ov.hbond_records[0]
  assert [atoms[h[k]].name.strip() for k in ('d_seq', 'h_seq', 'a_seq')] == \
    ['O1', 'HO1', 'O1']
  assert [atoms[h[k]].parent().parent().resseq.strip() for k in
    ('d_seq', 'h_seq', 'a_seq')] == ['1', '1', '2']
  assert approx_equal(atoms[h['h_seq']].distance(atoms[h['a_seq']]), h['d_HA'],
    eps=1.e-6)

def exercise_imports():
  '''
  A fresh process that imports ligand_interactions and runs it loads neither
  validate_ligands nor mmtbx.nci.hbond through the module itself (pnp imports
  mmtbx.nci.hbond), nor anything from nci_analysis.
  '''
  code = '\n'.join([
    'from __future__ import print_function',
    'import sys',
    'import mmtbx.validation.ligand_interactions',
    'from mmtbx.regression import tst_ligand_interactions as T',
    'T.get_manager(T.get_model())',
    'print(" ".join(sorted([k for k in sys.modules if',
    '  k == "mmtbx.validation.validate_ligands" or',
    '  k.split(".")[0] == "nci_analysis"])))'])
  fn = 'tst_ligand_interactions_imports.py'
  with open(fn, 'w') as f:
    f.write(code + '\n')
  r = easy_run.fully_buffered('libtbx.python %s' % fn)
  assert r.return_code == 0, '\n'.join(r.stderr_lines)
  assert [l for l in r.stdout_lines if l.strip()] == [], r.stdout_lines
  src = open(LI.__file__.replace('.pyc', '.py')).read()
  for name in ('validate_ligands', 'mmtbx.nci', 'nci_analysis'):
    assert ('import %s' % name) not in src and ('from %s' % name) not in src, name

def exercise_inline_clash():
  '''
  An atom clashing with two bonded atoms in line with it is one clash (pnp's rule
  in _process_clashes): EDO C1-H11 points at His NE2 (H11...NE2 1.78 A, C1...NE2
  2.63 A). probe2 has both pairs as bad overlaps, pnp keeps H11...NE2 only. One
  entry lists both pairs and marks the one pnp kept; it counts once. With pnp's
  record removed (a probe2-only inline case), the shorter pair represents it.
  '''
  model = get_model(inline_model_str.split('\n'))
  m = get_manager(model, sel='chain A and resseq 1')
  cl = [e for e in m.entries if e['type'] == 'clash']
  assert len(cl) == 1, cl
  e = cl[0]
  assert [short(l) for l in e['labels']] == ['EDO 1 H11', 'HIS 10 NE2']
  assert e['cross_check'] == 'pnp and probe2' and m.disagreements == []
  pairs = [([short(l) for l in q['labels']], q['sources'], q['pnp_kept'],
    q['representative']) for q in e['pairs']]
  assert pairs == [(['EDO 1 H11', 'HIS 10 NE2'], ['pnp', 'probe2'], True, True),
                   (['EDO 1 C1', 'HIS 10 NE2'], ['probe2'], False, False)], pairs
  assert approx_equal(e['pairs'][1]['distance'], 2.63, eps=0.01)
  assert m.overlaps.clash_records[0]['distance'] < e['pairs'][1]['distance']
  c = m.counts()
  assert c.per_type['clash'] == 1 and c.per_residue['B HIS 10']['clash'] == 1
  assert 'clash' in c.per_ligand_atom['A EDO 1 H11']
  assert 'clash' not in c.per_ligand_atom.get('A EDO 1 C1', {})
  assert [v for k, v in m.probe_pairs.items()
    if [short(atom_label_of(model, x)) for x in k] == ['EDO 1 C1', 'HIS 10 NE2']
    ][0]['pair_class'] == 'bo'
  r = m.clash_criteria['inline_merge']
  assert r['cos_min'] == 0.707 and r['function'] == 'pnp.cos_vec'
  # probe2 only
  m.overlaps.clash_records = []
  m.entries, m.disagreements = [], []
  m._build_clashes()
  assert len(m.entries) == 1
  e = m.entries[0]
  assert e['cross_check'] == 'probe2 only' and len(m.disagreements) == 1
  assert [short(l) for l in e['labels']] == ['EDO 1 H11', 'HIS 10 NE2']
  assert [(q['pnp_kept'], q['representative']) for q in e['pairs']] == \
    [(False, True), (False, False)]

def atom_label_of(model, i):
  return LI.atom_label(model.get_hierarchy().atoms()[i])

def exercise_probe_classes(model):
  '''
  Every probe2 class is parsed: wh (allow_weak_hydrogen_bonds), wo
  (separate_worse_clashes); a class not in pair_class_order or unknown to probe2
  raises Sorry; the order must include wh and wo when they can occur. wo and bo
  are clashes, hb and wh H-bonds, wh-only pairs with subtype "weak (probe2)".
  With the default options nothing changes.
  '''
  line = ':1->2:%s: A   1 EDO  H11 : A   3 EDO  H12 :-0.694:-0.393:0:0:0:0:0:C:C:0:0:0:1:1'
  text = '\n'.join([line % c for c in ('wo', 'bo', 'bo', 'wh')])
  p = LI.parse_probe_raw(text, LI.extended_pair_class_order)
  (v,) = p.pairs.values()
  assert v['dots'] == dict(wo=1, bo=2, hb=0, wh=1, so=0, cc=0, wc=0), v['dots']
  assert v['pair_class'] == 'wo'
  for order, cls in ((('bo', 'hb', 'so', 'cc', 'wc'), 'wh'),
                     (LI.extended_pair_class_order, 'zz')):
    try:
      LI.parse_probe_raw(line % cls, order)
    except Sorry as e:
      assert cls in str(e)
    else:
      raise AssertionError('class %s accepted' % cls)
  # the order must include the classes that can occur
  for weak, worse in ((True, False), (False, True)):
    params = LI.master_params().extract().ligand_interactions
    params.probe.allow_weak_hydrogen_bonds = weak
    params.separate_worse_clashes = worse
    try:
      LI.manager(model, model.selection(LIG_SEL).iselection(), LIG_SEL, params=params)
    except Sorry as e:
      assert 'wo bo hb wh so cc wc' in str(e), str(e)
    else:
      raise AssertionError('default order accepted with wh/wo')
  default = get_manager(model)
  both = get_manager(model, pair_class_order=list(LI.extended_pair_class_order),
    separate_worse_clashes=True, probe=_probe(allow_weak_hydrogen_bonds=True))
  assert both.counts().per_type == default.counts().per_type
  hb = [e for e in both.entries if e['type'] == 'hbond'][0]
  assert hb['geometry']['probe2']['dots']['wh'] > 0 and hb['subtype'] is None
  cl = [e for e in both.entries if e['type'] == 'clash'][0]
  assert cl['geometry']['probe2']['dots']['wo'] > 0
  assert set(default.patches[0]['dots']) == set(LI.probe_classes)
  # a wh-only pair: EDO 1 O1-HO1...EDO 2 O1 with H...O 2.66 A (pnp H-bond; probe2
  # has no overlap: cc by default, wh with weak H-bonds)
  wm = get_model(weak_model_str.split('\n'))
  w0 = get_manager(wm, sel='chain A and resseq 1')
  e0 = [e for e in w0.entries if e['type'] == 'hbond'][0]
  assert e0['cross_check'] == 'pnp only' and e0['subtype'] is None
  assert w0.disagreements[0]['probe_class'] == 'cc'
  w1 = get_manager(wm, sel='chain A and resseq 1',
    pair_class_order=list(LI.extended_pair_class_order),
    probe=_probe(allow_weak_hydrogen_bonds=True))
  e1 = [e for e in w1.entries if e['type'] == 'hbond'][0]
  assert e1['cross_check'] == 'pnp and probe2' and e1['subtype'] == 'weak (probe2)'
  assert e1['geometry']['probe2']['dots']['wh'] > 0
  assert e1['geometry']['probe2']['dots']['hb'] == 0
  assert w1.counts().per_type['hbond:weak (probe2)'] == 1

def _probe(**kwargs):
  p = LI.master_params().extract().ligand_interactions.probe
  for k, v in kwargs.items():
    setattr(p, k, v)
  return p

def exercise_hbond_phil():
  '''
  pnp's H-bond criteria from PHIL: the defaults are pnp.h_bond()'s; at
  a_DHA_cutoff 110 pnp finds the fixture's 113-degree H-bond, which it misses at
  the default 120; validate_ligands keeps the defaults.
  '''
  p = LI.master_params().extract().ligand_interactions.hbond
  h, d = LI.h_bond_params(p), pnp.h_bond()
  for k in ('Hs', 'As', 'Ds', 'd_HA_cutoff', 'd_DA_cutoff', 'a_DHA_cutoff',
            'a_YAH_cutoff', 'min_bonds_H_A'):
    assert list(getattr(h, k)) == list(getattr(d, k)) if isinstance(getattr(d, k), list) \
      else getattr(h, k) == getattr(d, k), k
  model = get_model(alt_model_str.split('\n'))
  sel = 'chain A and resseq 1'
  m0 = get_manager(model, sel=sel)
  assert m0.overlaps.hbond_records == []
  hp = LI.master_params().extract().ligand_interactions.hbond
  hp.a_DHA_cutoff = 110
  m1 = get_manager(model, sel=sel, hbond=hp)
  atoms = model.get_hierarchy().atoms()
  found = [(atom_label_of(model, r['h_seq']), atom_label_of(model, r['a_seq']),
    r['a_DHA']) for r in m1.overlaps.hbond_records]
  assert [(f[0], f[1]) for f in found if f[2] < 120] == [('A EDO 1 HO1 alt A', 'A EDO 2 O1')], found
  assert m1.as_dict()['hbond_criteria']['a_DHA_cutoff'] == 110
  e = [x for x in m1.entries if x['type'] == 'hbond' and x['labels'][1] == 'A EDO 1 HO1 alt A'][0]
  assert e['cross_check'] == 'pnp and probe2'
  assert approx_equal(e['geometry']['pnp']['a_DHA'], e['geometry']['probe2']['a_DHA'], eps=1.e-6)
  from mmtbx.validation import validate_ligands
  params = validate_ligands.master_params().extract().validate_ligands
  lr = validate_ligands.ligand_result(model=model, fmodel=None, map_manager=None,
    ligand_isel=model.selection(sel).iselection(), sel_str=sel, params=params)
  assert lr.get_overlaps().n_hbonds == 0

def exercise_contact_area(m):
  '''
  Contact area = dots / probe2's dot density, per side (ligand atom's surface:
  ligand -> environment dots; environment atom's surface: the reverse), per class
  and in total; the two sides add up to the pair's dots. Patches: ligand surface.
  '''
  density = m.params.probe.density
  n = 0
  for e in m.entries:
    g = e['geometry'].get('probe2')
    if g is None:
      continue
    n += 1
    if e['type'] == 'hbond':
      h, a = e['atoms'][1:]
      i, j = (h, a) if h in m._lig else (a, h)
    else:
      i, j = e['atoms']
    for side, (s_, t_) in (('ligand', (i, j)), ('environment', (j, i))):
      recount = {}
      for d in m.probe_dots:
        if d['source'] == s_ and d['target'] == t_:
          recount[d['cls']] = recount.get(d['cls'], 0) + 1
      for c, k in g['dots_%s' % side].items():
        assert k == recount.get(c, 0), (side, c, k, recount)
        assert approx_equal(g['area_%s' % side][c], k / density, eps=1.e-12)
      assert approx_equal(g['area_%s_total' % side],
        sum(g['dots_%s' % side].values()) / density, eps=1.e-12)
    for c in g['dots']:
      assert g['dots'][c] == g['dots_ligand'][c] + g['dots_environment'][c]
    assert 'area' not in g and 'area_total' not in g
  assert n >= 5
  sides = [e['geometry']['probe2'] for e in m.entries if e['type'] == 'vdw']
  assert [g for g in sides if g['area_ligand_total'] != g['area_environment_total']]
  for p in m.patches:
    assert approx_equal(p['area_total'], p['n_dots'] / density, eps=1.e-12)
    for c, k in p['dots'].items():
      assert approx_equal(p['area'][c], k / density, eps=1.e-12)

def salt_bridges(m):
  return [e for e in m.entries if e['type'] == 'salt_bridge']

def exercise_salt_bridges():
  '''
  Salt bridges on inline models (ideal CCD/GeoStd geometry; ligand groups from the
  builder, amino acids by template):
    ACT 1 with Lys 10 (NZ...OXT 2.85 A, NZ-HZ1...OXT H-bond) and Arg 20;
    ACT 2 with His 30 (HD1 and HE2), Asp 40 (same charge), Lys 50 (5.0 A);
    NH4 3 with Asp 60; the neutral model: acetic acid (ACY, custom restraints)
    at ACT 1's place, His 30 with HD1 only.
  '''
  model = get_model(salt_model_str.split('\n'))
  m1 = get_manager(model, sel='chain A and resseq 1')
  assert m1.charged_groups == [dict(kind='carboxylate', charge=-1, usual_charge=None,
    state='modelled', source='builder', altloc='', atoms=['A ACT 1 O', 'A ACT 1 OXT'],
    residue='A ACT 1', certain=True, metal_bound=False, charge_source='formal charges',
    hydrogens='model', notes=[])], m1.charged_groups
  sb = dict([(e['residue'], e) for e in salt_bridges(m1)])
  assert sorted(sb) == ['B LYS 10', 'C ARG 20'], sorted(sb)
  lys = sb['B LYS 10']
  g = lys['geometry']['charged_groups']
  assert lys['subtype'] == 'salt bridge (K&N)', lys['subtype']
  assert approx_equal(g['min_atom_distance'], 2.85, eps=0.01)
  assert g['closest_pair'] == ['A ACT 1 OXT', 'B LYS 10 NZ']
  assert g['partner_group']['kind'] == 'ammonium' and g['partner_group']['charge'] == 1
  assert lys['ligand_atoms'] == ['A ACT 1 O', 'A ACT 1 OXT']
  assert [h['labels'] for h in lys['hbonds']] == [['B LYS 10 NZ', 'B LYS 10 HZ1', 'A ACT 1 OXT']]
  assert m1.entries[lys['hbonds'][0]['index']]['type'] == 'hbond'
  arg = sb['C ARG 20']
  assert arg['geometry']['charged_groups']['partner_group']['kind'] == 'guanidinium'
  assert arg['subtype'] == 'N-O bridge (K&N)'
  assert arg['geometry']['charged_groups']['charge_centre_distance'] > 4.0
  c = m1.counts()
  assert c.per_residue['B LYS 10']['hbond'] == 1
  assert c.per_residue['B LYS 10']['salt_bridge:salt bridge (K&N)'] == 1
  m2 = get_manager(model, sel='chain A and resseq 2')
  sb = salt_bridges(m2)
  assert [e['residue'] for e in sb] == ['D HIS 30'], sb
  assert sb[0]['geometry']['charged_groups']['partner_group']['kind'] == 'imidazolium'
  m2c = get_manager(model, sel='chain A and resseq 2',
    salt_bridge=_salt(criterion='charge_centre'))
  sbc = dict([(e['residue'], e['subtype']) for e in salt_bridges(m2c)])
  assert sbc == {'D HIS 30': 'N-O bridge (K&N)', 'F LYS 50': 'longer-range (K&N)'}, sbc
  m3 = get_manager(model, sel='chain A and resseq 3')
  sb = salt_bridges(m3)
  assert [(e['residue'], e['geometry']['charged_groups']['ligand_group']['kind'])
    for e in sb] == [('G ASP 60', 'ammonium')], sb
  assert m1.formal_charge_conflicts() == [] and m3.formal_charge_conflicts() == []
  assert [(x['source'], x['status']) for x in m1.formal_charges
    if x['residue'] == 'A ACT 1'] == [('restraints', 'agrees'), ('CCD', 'agrees')]
  # neutral: protonated carboxylate, singly protonated His
  neutral = get_model(neutral_model_str.split('\n'), cifs=(acy_neutral_cif,))
  n1 = get_manager(neutral, sel='chain A and resseq 1')
  assert n1.charged_groups == [] and salt_bridges(n1) == []
  assert [e['labels'] for e in n1.entries if e['type'] == 'hbond'] == \
    [['B LYS 10 NZ', 'B LYS 10 HZ1', 'A ACY 1 OXT']]
  n2 = get_manager(neutral, sel='chain A and resseq 2')
  assert salt_bridges(n2) == []
  # the CCD's HIS (+1 on ND1) is the doubly protonated form: with HD1 only the
  # protonation differs (no conflict); with HD1 and HE2 the charges agree
  assert n2.formal_charge_conflicts() == []
  assert [(x['group'], x['perceived_charge'], x['status']) for x in n2.formal_charges
    if x['residue'] == 'D HIS 30' and x['group'] == 'imidazolium'] == \
    [('imidazolium', 0, 'protonation differs')]
  assert [(x['group'], x['status']) for x in m2.formal_charges
    if x['residue'] == 'D HIS 30' and x['group'] == 'imidazolium'] == [('imidazolium', 'agrees')]

def _salt(**kwargs):
  p = LI.master_params().extract().ligand_interactions.salt_bridge
  for k, v in kwargs.items():
    setattr(p, k, v)
  return p

def exercise_salt_bridge_symmetry():
  '''
  Lys 10 moved by -a (100 A): its x+1,y,z copy makes the salt bridge; the ligand
  stays, the operator moves the Lys atoms and reproduces the recorded distance.
  '''
  model = get_model(salt_sym_model_str.split('\n'))
  m = get_manager(model, sel='chain A and resseq 1')
  (e,) = salt_bridges(m)
  assert e['residue'] == 'B LYS 10 (x+1,y,z)' and e['symop'] == 'x+1,y,z'
  assert e['operators'] == ['x,y,z', 'x,y,z', 'x+1,y,z']
  assert e['labels'] == ['A ACT 1 O', 'A ACT 1 OXT', 'B LYS 10 NZ (x+1,y,z)']
  uc = model.crystal_symmetry().unit_cell()
  atoms = model.get_hierarchy().atoms()
  def site(i, op):
    return uc.orthogonalize(sgtbx.rt_mx(op) * uc.fractionalize(atoms[i].xyz))
  d = min([atoms[e['atoms'][k]].distance(site(e['atoms'][2], e['operators'][2]))
    for k in (0, 1)])
  g = e['geometry']['charged_groups']
  assert approx_equal(d, g['min_atom_distance'], eps=1.e-6)
  assert approx_equal(d, 2.85, eps=0.01)
  assert atoms[e['atoms'][1]].distance(atoms[e['atoms'][2]]) > 90
  assert [h['labels'][1] for h in e['hbonds']] == ['B LYS 10 HZ1 (x+1,y,z)']

def exercise_formal_charge_conflict():
  '''
  ACT restraints with OXT's formal charge set to 0 (no H on the carboxylate, as in
  the model): the builder takes the file's charges, DetermineBondOrders fails at
  total 0, so ACT has no groups (listed), no salt bridge, and its formal charges
  are not compared. The conflict status itself: a dictionary charge differing from
  the perception with the same H.
  '''
  model = get_model(salt_sym_model_str.split('\n'), cifs=(act_conflict_cif,))
  m = get_manager(model, sel='chain A and resseq 1')
  assert m.charged_groups == [] and len(salt_bridges(m)) == 0
  (f,) = m.charged_group_failures
  assert f['residue'] == 'A ACT 1' and f['reason'].startswith(
    'DetermineBondOrders fails for ACT with the formal total 0:'), f
  act = [x for x in m.formal_charges if x['residue'] == 'A ACT 1']
  assert [(x['source'], x['file'], x['status'].startswith(
    'not compared (residue_molecule failed: ')) for x in act] == [
    ('restraints', 'act_conflict.cif', True), ('CCD', None, True)], act
  assert m.formal_charge_conflicts() == []
  d = dict(atoms=dict(C=('C', 0), O=('O', 0), OXT=('O', 0)), h=dict(C=set(), O=set(),
    OXT=set()))
  c = m._charge_check('A ACT 1', 'restraints', 'x.cif', 'carboxylate', ['C', 'O', 'OXT'],
    ['C', 'O', 'OXT'], -1, d, dict(C=set(), O=set(), OXT=set()))
  assert (c['status'], c['dictionary_charge'], c['perceived_charge']) == ('conflict', 0, -1)
  m.formal_charges.append(c)
  log = StringIO()
  m.show(log=log)
  assert 'formal-charge conflicts' in log.getvalue()
  assert 'residues without charged groups (residue_molecule failed):' in log.getvalue()

def builder_rows(model, found):
  atoms = model.get_hierarchy().atoms()
  return sorted([(LI.residue_label(atoms[g['center']], g['resname']), g['kind'], g['charge'],
    atoms[g['center']].name.strip(), tuple(sorted([atoms[i].name.strip() for i in g['charged']])),
    g['certain'], g['metal_bound']) for g in found.groups if g['source'] == 'builder'])

def exercise_builder_groups():
  '''
  Non-standard residues: groups from rdkit_utils.residue_molecule (fixtures of
  tst_rdkit_utils_molecule). Taurine zwitterion ZZW: ammonium and sulfonate;
  methyl triphosphate ZTP: -1, -1, -2 (-4); N-methylpyridinium: "other"; nitro:
  no group (balanced by the bonded O-); partial charges only (ZPC): uncertain;
  carboxylate O on a Zn of another residue: metal_bound; ZZW's N linked to an
  acetyl C: dropped (capped); ACT ester-linked to Ser OG: no carboxylate; ZC5: the
  builder fails, listed, no groups.
  '''
  from mmtbx.regression import tst_rdkit_utils_molecule as M
  def found_for(model):
    return LI.find_charged_groups(model, all_atoms(model))
  model = M.get_model(M.pdb_from_cif('ZZW', M.zzw_cif), cifs=(('ZZW', M.zzw_cif),))
  f = found_for(model)
  assert builder_rows(model, f) == [
    ('A ZZW 1', 'ammonium', 1, 'N1', ('N1',), True, False),
    ('A ZZW 1', 'sulfonate', -1, 'S1', ('O1', 'O2', 'O3'), True, False)], builder_rows(model, f)
  g = f.groups[0]
  assert (g['source'], g['charge_source'], g['hydrogens'], g['state'], g['notes'],
    g['usual_charge']) == ('builder', 'formal charges', 'model', 'modelled', [], None), g
  # an N-H missing in the model: completed from the file, the ammonium uncertain
  model = M.get_model(M.pdb_from_cif('ZZW', M.zzw_cif, drop=('HN13',)),
    cifs=(('ZZW', M.zzw_cif),))
  f = found_for(model)
  assert [(g['kind'], g['certain'], g['notes']) for g in f.groups] == [
    ('ammonium', False, ['H completed from the restraint file on N1']),
    ('sulfonate', True, [])], f.groups
  # without H: assumed, the file's H not held against the groups
  hs = [l.split()[1] for l in M.zzw_cif.split('\n') if l.startswith(' ZZW ') and
    len(l.split()) == 9 and l.split()[2] == 'H']
  model = M.get_model(M.pdb_from_cif('ZZW', M.zzw_cif, drop=hs), cifs=(('ZZW', M.zzw_cif),))
  f = found_for(model)
  assert [(g['kind'], g['certain'], g['state']) for g in f.groups] == [
    ('ammonium', True, 'assumed (no H)'), ('sulfonate', True, 'assumed (no H)')], f.groups
  model = M.get_model(M.pdb_from_cif('ZTP', M.ztp_cif), cifs=(('ZTP', M.ztp_cif),))
  rows = sorted(builder_rows(model, found_for(model)), key=lambda r: r[3])
  assert [(r[1], r[2], r[3], r[4]) for r in rows] == [
    ('phosphate', -1, 'P1', ('O11', 'O12')), ('phosphate', -1, 'P2', ('O21', 'O22')),
    ('phosphate', -2, 'P3', ('O31', 'O32', 'O33'))], rows
  model = get_model(pdb_from_cif_text('ZPY', zpy_cif[1]), cifs=(zpy_cif,))
  assert builder_rows(model, found_for(model)) == [
    ('A ZPY 1', 'other', 1, 'N1', ('N1',), True, False)]
  model = get_model(pdb_from_cif_text('ZNX', znx_cif[1]), cifs=(znx_cif,))
  f = found_for(model)
  assert f.groups == [] and f.failures == [] and f.dropped == [dict(residue='A ZNX 1',
    altloc='', kind='charge-separated', charge=0, atoms=['N1', 'O2'],
    reason='charge-separated, net 0')], f.dropped
  # partial charges only: uncertain, not paired with the ZZW ammonium 3 A away
  model = M.get_model(M.pdb_from_cif('ZPC', M.zpc_cif), cifs=(('ZPC', M.zpc_cif),))
  f = found_for(model)
  (g,) = f.groups
  assert (g['kind'], g['charge'], g['certain'], g['charge_source']) == ('carboxylate', -1,
    False, 'search'), g
  assert g['notes'] == ['total charge by search (no formal charges)']
  lines = [l for l in M.pdb_from_cif('ZPC', M.zpc_cif).split('\n') if l.startswith('HETATM')]
  o1 = [float(x) for x in [l for l in lines if l[12:16] == ' O1 '][0][30:54].split()]
  zzw = [l for l in M.pdb_from_cif('ZZW', M.zzw_cif).split('\n') if l.startswith('HETATM')]
  n1 = [float(x) for x in [l for l in zzw if l[12:16] == ' N1 '][0][30:54].split()]
  shift = [o1[0] - n1[0] + 3.0, o1[1] - n1[1], o1[2] - n1[2]]
  moved = [l[:21] + 'B' + l[22:30] + ''.join(['%8.3f' % (float(l[30 + 8 * k:38 + 8 * k]) +
    shift[k]) for k in range(3)]) + l[54:] for l in zzw]
  model = M.get_model('\n'.join(['CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1'] +
    lines + moved + ['END']), cifs=(('ZPC', M.zpc_cif), ('ZZW', M.zzw_cif)))
  f = found_for(model)
  assert [(g['resname'], g['kind'], g['certain']) for g in f.groups] == [
    ('ZPC', 'carboxylate', False), ('ZZW', 'ammonium', True), ('ZZW', 'sulfonate', True)], \
    [(g['resname'], g['kind'], g['certain']) for g in f.groups]
  def with_zpc(groups):
    return [p for p in LI.charged_group_pairs(model, groups, 4.0) if 0 in p[:2]]
  assert with_zpc(f.groups) == []
  assert len(with_zpc([dict(g, certain=True) for g in f.groups])) == 2
  # carboxylate O on a Zn
  b = iotbx.cif.reader(input_string=M.zac_cif).model()['comp_ZAC']
  xyz = dict([(n, [float(b['_chem_comp_atom.%s' % c][k]) + 20 for c in 'xyz'])
    for k, n in enumerate(b['_chem_comp_atom.atom_id'])])
  o, c = xyz['O2'], xyz['C1']
  d = [o[k] - c[k] for k in range(3)]
  n = sum([v * v for v in d]) ** 0.5
  zn = 'HETATM   99 ZN    ZN Z   1    %8.3f%8.3f%8.3f  1.00 20.00          ZN' % tuple(
    [o[k] + 2.0 * d[k] / n for k in range(3)])
  edits = M.zinc_edits.replace('name S1', 'name O2').replace('2.30', '2.00')
  model = M.get_model(M.pdb_from_cif('ZAC', M.zac_cif, extra=(zn,)), cifs=(('ZAC', M.zac_cif),),
    edits=edits)
  f = found_for(model)
  assert builder_rows(model, f) == [('A ZAC 1', 'carboxylate', -1, 'C1', ('O1', 'O2'), True,
    True)], builder_rows(model, f)
  assert f.groups[0]['notes'] == ['metal-bound: O2']
  # ZZW N1 bonded to an acetyl C (one N-H less): the ammonium's N is capped
  a = M.pdb_from_cif('ZZW', M.zzw_cif, drop=('HN13',)).split('\n')
  acetyl = [l[:17] + 'ZAC B' + l[22:] for l in M.pdb_from_cif('ZAC', M.zac_cif,
    drop=('O2',)).split('\n') if l.startswith('HETATM')]
  edits = M.zinc_edits.replace('chain A and resseq 1 and name S1',
    'chain A and resseq 1 and name N1').replace('chain Z and resseq 1 and name ZN',
    'chain B and resseq 1 and name C1').replace('2.30', '1.33')
  model = M.get_model('\n'.join(a[:-1] + acetyl + ['END']), cifs=(('ZZW', M.zzw_cif),
    ('ZAC', M.zac_cif)), edits=edits)
  f = found_for(model)
  assert [(g['resname'], g['kind']) for g in f.groups] == [('ZZW', 'sulfonate')], f.groups
  assert f.dropped == [dict(residue='A ZZW 1', altloc='', kind='ammonium', charge=1,
    atoms=['N1'], reason='capped: N1 linked to B ZAC 1 C1')], f.dropped
  # ACT B 2 ester-linked to Ser OG: no carboxylate
  model = M.get_model(M.ester_pdb, edits=M.ester_edits)
  f = found_for(model)
  assert [g for g in f.groups if g['resname'] == 'ACT'] == [] and f.failures == []
  # builder failure
  model = M.get_model(M.pdb_from_cif('ZC5', M.zc5_cif), cifs=(('ZC5', M.zc5_cif),))
  f = found_for(model)
  assert f.groups == [] and [x['residue'] for x in f.failures] == ['A ZC5 1']
  assert f.failures[0]['reason'].startswith('DetermineBondOrders fails for ZC5')
  # the manager: ZZW 6 A away gives probe2 polar H
  zzw = [l[:21] + 'B' + l[22:30] + '%8.3f' % (float(l[30:38]) + 6.0) + l[38:]
    for l in M.pdb_from_cif('ZZW', M.zzw_cif).split('\n') if l.startswith('HETATM')]
  zc5 = M.pdb_from_cif('ZC5', M.zc5_cif).split('\n')
  model = M.get_model('\n'.join(zc5[:-1] + zzw + ['END']), cifs=(('ZC5', M.zc5_cif),
    ('ZZW', M.zzw_cif)))
  m = get_manager(model, sel='resname ZC5')
  assert m.charged_groups == [] and m.charged_group_failures == f.failures
  assert m.as_dict()['charged_group_failures'] == f.failures
  assert [(x['source'], x['status'].startswith('not compared (residue_molecule failed: '))
    for x in m.formal_charges if x['residue'] == 'A ZC5 1'] == [('restraints', True)], \
    m.formal_charges
  log = StringIO()
  m.show(log=log)
  assert 'residues without charged groups (residue_molecule failed):' in log.getvalue()
  assert '    A ZC5 1: DetermineBondOrders fails for ZC5' in log.getvalue()

def pdb_from_cif_text(code, text):
  from mmtbx.regression import tst_rdkit_utils_molecule as M
  return M.pdb_from_cif(code, text).split('\n')

def geostd_text(code):
  from mmtbx.monomer_library import server
  import os
  return open(os.path.join(server.server().geostd_path, code[0].lower(),
    'data_%s.cif' % code)).read()

def exercise_charge_separated():
  '''
  Charged atoms bonded to an opposite charge form clusters. GeoStd nitrate NO3
  (N +1, O2 and O3 -1): a group of -1 on the three O (centre N); azide ion AZI
  (N1 -1, N2 +1, N3 -1): -1 on N1 and N3; trimethylamine N-oxide TMO and the ZNX
  nitro: net 0, no group, listed as dropped. Nitrate at ACT 1's place in
  salt_model_str: salt bridges with Lys 10 and Arg 20, its formal charges agreeing.
  '''
  from mmtbx.regression import tst_rdkit_utils_molecule as M
  rows = {}
  for code in ('NO3', 'AZI', 'TMO'):
    model = M.get_model(M.pdb_from_cif(code, geostd_text(code)))
    atoms = model.get_hierarchy().atoms()
    f = LI.find_charged_groups(model, all_atoms(model))
    rows[code] = ([(g['kind'], g['charge'], atoms[g['center']].name.strip(),
      sorted([atoms[i].name.strip() for i in g['charged']]), g['certain']) for g in f.groups],
      [(d['atoms'], d['charge'], d['reason']) for d in f.dropped])
  assert rows['NO3'] == ([('other', -1, 'N', ['O1', 'O2', 'O3'], True)], []), rows['NO3']
  assert rows['AZI'] == ([('other', -1, 'N2', ['N1', 'N3'], True)], []), rows['AZI']
  assert rows['TMO'] == ([], [(['NAC', 'OAE'], 0, 'charge-separated, net 0')]), rows['TMO']
  # nitrate replacing ACT 1: N at C, O1 at O, O2 at OXT, O3 1.25 A from N away from CH3
  lines = []
  for l in salt_model_str.split('\n'):
    if l[17:26] == 'ACT A   1':
      n = l[12:16].strip()
      if n in ('C', 'O', 'OXT'):
        lines.append(l[:12] + {'C': ' N  ', 'O': ' O1 ', 'OXT': ' O2 '}[n] + ' NO3' + l[20:76] +
          (' N' if n == 'C' else ' O'))
      if n == 'C':
        x = float(l[30:38]) + 1.25
        lines.append(l[:12] + ' O3  NO3' + l[20:30] + '%8.3f' % x + l[38:76] + ' O')
      continue
    lines.append(l)
  model = get_model(lines)
  m = get_manager(model, sel='chain A and resseq 1')
  assert [(g['kind'], g['charge'], g['atoms']) for g in m.charged_groups] == [('other', -1,
    ['A NO3 1 O1', 'A NO3 1 O2', 'A NO3 1 O3'])], m.charged_groups
  assert sorted([e['residue'] for e in salt_bridges(m)]) == ['B LYS 10', 'C ARG 20'], \
    [e['residue'] for e in salt_bridges(m)]
  no3 = [(x['source'], x['atoms'], x['perceived_charge'], x['dictionary_charge'],
    x['status']) for x in m.formal_charges if x['residue'] == 'A NO3 1']
  assert ('restraints', ['N'], 1, 1, 'agrees') in no3 and \
    ('restraints', ['O1', 'O2', 'O3'], -2, -2, 'agrees') in no3, no3
  assert m.formal_charge_conflicts() == []

def exercise_atomic_charge_comparison():
  '''
  Builder residues are compared atom by atom (resonance-equivalent atoms summed):
  the ZNX nitro (N1 +1, O2 -1 in the file and the builder) agrees, no conflict
  (before, N1 and O2 were compared with 0 as atoms outside the groups; the CCD's
  ZNX is another compound: protonation differs). A conflict
  is still reported: ACT against a stand-in CCD with OXT's charge 0 and the same H
  (O, OXT: builder -1, stand-in 0).
  '''
  from mmtbx.regression import tst_rdkit_utils_molecule as M
  zzw = [l[:21] + 'B' + l[22:30] + '%8.3f' % (float(l[30:38]) + 6.0) + l[38:]
    for l in M.pdb_from_cif('ZZW', M.zzw_cif).split('\n') if l.startswith('HETATM')]
  znx = pdb_from_cif_text('ZNX', znx_cif[1])
  model = get_model(znx[:-1] + zzw + ['END'], cifs=(znx_cif, ('zzw.cif', M.zzw_cif)))
  m = get_manager(model, sel='resname ZNX')
  rows = sorted([(x['source'], x['atoms'], x['perceived_charge'], x['dictionary_charge'],
    x['status']) for x in m.formal_charges if x['residue'] == 'A ZNX 1' and
    x['source'] == 'restraints'])
  assert rows == [('restraints', ['N1'], 1, 1, 'agrees'),
    ('restraints', ['O1', 'O2'], -1, -1, 'agrees')], rows
  assert m.formal_charge_conflicts() == []
  assert [d['reason'] for d in m.charged_groups_dropped if d['residue'] == 'A ZNX 1'] == [
    'charge-separated, net 0']
  original = LI.ccd_formal_charges
  def stand_in(resname):
    d = original(resname)
    if resname == 'ACT':
      d['atoms']['OXT'] = ('O', 0)
    return d
  LI.ccd_formal_charges = stand_in
  try:
    model = get_model(salt_sym_model_str.split('\n'))
    m = get_manager(model, sel='chain A and resseq 1')
  finally:
    LI.ccd_formal_charges = original
  (c,) = m.formal_charge_conflicts()
  assert (c['residue'], c['source'], c['atoms'], c['perceived_charge'],
    c['dictionary_charge']) == ('A ACT 1', 'CCD', ['O', 'OXT'], -1, 0), c
  assert [x['status'] for x in m.formal_charges if x['residue'] == 'A ACT 1' and
    x['source'] == 'restraints'] == ['agrees']

def possible(m):
  return [(e['residue'], e['subtype'], e['charges'], e['reasons']) for e in
    m.possible_salt_bridges]

def add_h(lines, residue, heavy, parent, name, alt=None):
  '''lines with an H on heavy (0.97 A from it, away from parent) after it.'''
  o = [l for l in lines if l[17:26] == residue and l[12:16].strip() == heavy and
    (alt is None or l[16] == alt)][0]
  c = [l for l in lines if l[17:26] == residue and l[12:16].strip() == parent and
    l[16] == o[16]][0]
  xo = [float(o[30 + 8 * k:38 + 8 * k]) for k in range(3)]
  xc = [float(c[30 + 8 * k:38 + 8 * k]) for k in range(3)]
  d = [xo[k] - xc[k] for k in range(3)]
  n = sum([v * v for v in d]) ** 0.5
  h = o[:12] + ' %-3s' % name + o[16:30] + '%8.3f%8.3f%8.3f' % tuple(
    [xo[k] + 0.97 * d[k] / n for k in range(3)]) + o[54:76] + ' H'
  k = lines.index(o)
  return lines[:k + 1] + [h] + lines[k + 1:]

def exercise_possible_salt_bridges():
  '''
  Possible salt bridges: opposite charges, each group's modelled charge if certain
  and charged, else its usual charge; at least one not certainly charged; the
  salt-bridge criterion; listed apart, not counted. Neutral benzamidine (GeoStd
  BEN, at ACT 2's place) next to Asp 40; NH4 next to Asp 60 with HD2; ACT 1 next
  to Lys 10 with two H on NZ (Arg 20 stays a salt bridge); a neutral methyl
  phosphate (HOP2, HOP3; usual -2) at ACT 1's place next to Lys 10 and Arg 20; the
  partial-charges-only ZPC (uncertain) next to ZZW's ammonium; the Lys case moved
  by -a (symmetry); Asp 60 with HD2 in both conformers next to NH4 A/B (altlocs).
  '''
  from mmtbx.regression import tst_rdkit_utils_molecule as M
  salt = salt_model_str.split('\n')
  # neutral amidine
  lines = [l for l in salt if l[17:26] != 'ACT A   2' and l != 'END']
  model = get_model(lines + ben_at_act2_lines.strip('\n').split('\n') + ['END'])
  m = get_manager(model, sel='chain A and resseq 2')
  assert salt_bridges(m) == []
  assert possible(m) == [('E ASP 40', 'N-O bridge (K&N)', [1, -1],
    ['A BEN 2 amidine neutral as modelled (restraint file)'])], possible(m)
  e = m.possible_salt_bridges[0]
  g = e['geometry']['charged_groups']
  assert (g['ligand_group']['kind'], g['ligand_group']['charge'],
    g['ligand_group']['usual_charge'], g['ligand_group']['atoms']) == ('amidine', 0, 1,
    ['A BEN 2 N1', 'A BEN 2 N2']), g['ligand_group']
  assert (g['partner_group']['charge'], g['partner_group']['usual_charge'],
    g['partner_group']['state']) == (-1, -1, 'modelled')
  assert e['type'] == 'possible_salt_bridge' and e not in m.entries
  c = m.counts()
  assert not [t for t in c.per_type if 'salt' in t], c.per_type
  assert not [t for r in c.per_residue.values() for t in r if 'salt' in t], c.per_residue
  assert not [t for a in c.per_ligand_atom.values() for t in a if 'salt' in t]
  assert [(g['kind'], g['usual_charge']) for g in m.possible_groups] == [('amidine', 1)]
  d = m.as_dict()
  assert d['possible_salt_bridges'] == m.possible_salt_bridges
  assert [x['kind'] for x in d['possible_groups']] == ['amidine']
  json.dumps(d, default=str)
  log = StringIO()
  m.show(log=log)
  assert 'possible salt bridges (not counted):' in log.getvalue()
  assert 'A BEN 2 amidine neutral as modelled (restraint file)' in log.getvalue()
  # Asp 60 protonated (HD2) next to NH4
  model = get_model(add_h(salt, 'ASP G  60', 'OD2', 'CG', 'HD2'))
  m = get_manager(model, sel='chain A and resseq 3')
  assert salt_bridges(m) == []
  assert possible(m) == [('G ASP 60', 'salt bridge (K&N)', [1, -1],
    ['G ASP 60 carboxylate protonated (HD2)'])], possible(m)
  # Lys 10 with two H on NZ next to ACT 1
  lines = [l for l in salt if not (l[17:26] == 'LYS B  10' and l[12:16].strip() == 'HZ3')]
  m = get_manager(get_model(lines), sel='chain A and resseq 1')
  assert [e['residue'] for e in salt_bridges(m)] == ['C ARG 20']
  assert possible(m) == [('B LYS 10', 'salt bridge (K&N)', [-1, 1],
    ['B LYS 10 ammonium neutral (NZ with 2 H)'])], possible(m)
  assert [h['labels'] for h in m.possible_salt_bridges[0]['hbonds']] == [
    ['B LYS 10 NZ', 'B LYS 10 HZ1', 'A ACT 1 OXT']]
  # neutral phosphoric acid monoester (usual -2)
  lines = [l for l in salt if l[17:26] != 'ACT A   1' and l != 'END']
  model = get_model(lines + zmp_at_act1_lines.strip('\n').split('\n') + ['END'],
    cifs=(zmp_cif,))
  m = get_manager(model, sel='chain A and resseq 1')
  assert salt_bridges(m) == []
  reason = ['A ZMP 1 phosphoric acid neutral as modelled (restraint file)']
  assert sorted(possible(m)) == [('B LYS 10', 'salt bridge (K&N)', [-2, 1], reason),
    ('C ARG 20', 'N-O bridge (K&N)', [-2, 1], reason)], possible(m)
  assert [a.split()[-1] for a in m.possible_groups[0]['atoms']] == ['O1P', 'O2P', 'O3P']
  # an uncertain group (partial charges only) next to a certain ammonium
  lines = [l for l in M.pdb_from_cif('ZPC', M.zpc_cif).split('\n') if l.startswith('HETATM')]
  o1 = [float(x) for x in [l for l in lines if l[12:16] == ' O1 '][0][30:54].split()]
  zzw = [l for l in M.pdb_from_cif('ZZW', M.zzw_cif).split('\n') if l.startswith('HETATM')]
  n1 = [float(x) for x in [l for l in zzw if l[12:16] == ' N1 '][0][30:54].split()]
  shift = [o1[0] - n1[0] + 3.0, o1[1] - n1[1], o1[2] - n1[2]]
  moved = [l[:21] + 'B' + l[22:30] + ''.join(['%8.3f' % (float(l[30 + 8 * k:38 + 8 * k]) +
    shift[k]) for k in range(3)]) + l[54:] for l in zzw]
  model = get_model(['CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1'] + lines +
    moved + ['END'], cifs=(('zpc.cif', M.zpc_cif), ('zzw.cif', M.zzw_cif)))
  m = get_manager(model, sel='resname ZPC')
  assert salt_bridges(m) == []
  assert [(r, q, why) for r, s, q, why in possible(m)] == [('B ZZW 1', [-1, 1],
    ['A ZPC 1 carboxylate uncertain (total charge by search (no formal charges))'])], \
    possible(m)
  # symmetry: Lys 10 moved by -a, two H on NZ
  lines = [l for l in salt_sym_model_str.split('\n') if not (l[17:20] == 'LYS' and
    l[12:16].strip() == 'HZ3')]
  m = get_manager(get_model(lines), sel='chain A and resseq 1')
  (e,) = m.possible_salt_bridges
  assert (e['residue'], e['symop']) == ('B LYS 10 (x+1,y,z)', 'x+1,y,z'), e['residue']
  assert e['labels'] == ['A ACT 1 O', 'A ACT 1 OXT', 'B LYS 10 NZ (x+1,y,z)']
  assert approx_equal(e['geometry']['charged_groups']['min_atom_distance'], 2.85, eps=0.01)
  # altlocs: A with A, B with B
  lines = altloc_both_model_str.split('\n')
  for alt in 'AB':
    lines = add_h(lines, 'ASP G  60', 'OD2', 'CG', 'HD2', alt=alt)
  m = get_manager(get_model(lines), sel='chain A and resseq 3')
  assert salt_bridges(m) == []
  assert sorted([(e['geometry']['charged_groups']['ligand_group']['altloc'],
    e['geometry']['charged_groups']['partner_group']['altloc'])
    for e in m.possible_salt_bridges]) == [('A', 'A'), ('B', 'B')]

def exercise_possible_kinds():
  '''
  Neutral groups usually charged at pH 7 by SMARTS on the builder's molecule:
  guanidine (+1, the three N), 1H-tetrazole (-1, the NH N and its resonance
  partners), sulfonic acid (-1, the three O), a tertiary amine (+1; the amide N of
  the same molecule is not a group); SEP in a chain: phosphoric acid, the capped
  backbone N skipped silently (not listed).
  '''
  from mmtbx.regression import tst_rdkit_utils_molecule as M
  rows = []
  for code, cif in (('ZGN', zgn_cif), ('ZTZ', ztz_cif), ('ZSA', zsa_cif), ('ZAM', zam_cif)):
    model = M.get_model(M.pdb_from_cif(code, cif[1]), cifs=((code, cif[1]),))
    atoms = model.get_hierarchy().atoms()
    f = LI.find_charged_groups(model, all_atoms(model))
    assert f.groups == [] and f.failures == []
    rows += [(code, g['kind'], g['charge'], g['usual_charge'], atoms[g['center']].name.strip(),
      sorted([atoms[i].name.strip() for i in g['charged']]), g['notes']) for g in f.possible]
  why = ['neutral as modelled (restraint file)']
  assert rows == [('ZGN', 'guanidine', 0, 1, 'C2', ['N1', 'N2', 'N3'], why),
    ('ZTZ', 'tetrazole', 0, -1, 'C2', ['N1', 'N2', 'N4'], why),
    ('ZSA', 'sulfonic acid', 0, -1, 'S1', ['O1', 'O2', 'O3'], why),
    ('ZAM', 'amine', 0, 1, 'N1', ['N1'], why)], rows
  # SEP between two Gly (GeoStd SEP, HOP2 and HOP3 modelled): the phosphoric acid
  # (usual -2); the backbone N, linked to Gly 1 (capped), dropped
  model = M.get_model(M.chain_pdb)
  atoms = model.get_hierarchy().atoms()
  f = LI.find_charged_groups(model, all_atoms(model))
  assert [(g['resname'], g['kind'], g['usual_charge'], sorted([atoms[i].name.strip()
    for i in g['charged']])) for g in f.possible] == [('SEP', 'phosphoric acid', -2,
    ['O1P', 'O2P', 'O3P'])], f.possible
  assert f.dropped == [] and not hasattr(f, 'possible_dropped')

def placed(code, cif, atom, site, chain='A', resseq=9):
  '''pdb lines of the cif's residue translated so that atom lies on site.'''
  from mmtbx.regression import tst_rdkit_utils_molecule as M
  lines = [l for l in M.pdb_from_cif(code, cif[1]).split('\n') if l.startswith('HETATM')]
  a = [float(x) for x in [l for l in lines if l[12:16].strip() == atom][0][30:54].split()]
  return [l[:21] + chain + '%4d' % resseq + l[26:30] + ''.join(['%8.3f' % (
    float(l[30 + 8 * k:38 + 8 * k]) - a[k] + site[k]) for k in range(3)]) + l[54:]
    for l in lines]

def superposed(code, cif, pick, sites, chain='A', resseq=9):
  '''pdb lines of the cif's residue superposed so that the atoms pick lie on sites.'''
  from mmtbx.regression import tst_rdkit_utils_molecule as M
  from scitbx.math import superpose
  lines = [l for l in M.pdb_from_cif(code, cif[1]).split('\n') if l.startswith('HETATM')]
  xyz = dict([(l[12:16].strip(), [float(l[30 + 8 * k:38 + 8 * k]) for k in range(3)])
    for l in lines])
  f = superpose.least_squares_fit(reference_sites=flex.vec3_double(sites),
    other_sites=flex.vec3_double([xyz[n] for n in pick]))
  r, t = f.r.elems, f.t.elems
  out = []
  for l in lines:
    x = xyz[l[12:16].strip()]
    y = [sum([r[3 * i + j] * x[j] for j in range(3)]) + t[i] for i in range(3)]
    out.append(l[:21] + chain + '%4d' % resseq + l[26:30] + '%8.3f%8.3f%8.3f' % tuple(y) +
      l[54:])
  return out

def exercise_amine_and_tetrazole_rules():
  '''
  No possible salt bridge for an N bonded to N or O, or to a nitrile C: acetohydrazide
  (terminal N2) and N-methylhydroxylamine (N1) at NH4 3's N, next to Asp 60;
  dimethylcyanamide: no group. A tetrazole only with an H on a ring N:
  1-methyltetrazole at ACT 1's place, next to Lys 10: nothing; 5-methyl-1H-
  tetrazole (N4-H) there: possible salt bridges with Lys 10 and Arg 20 (both rings
  superposed with C5 and its two ring N on ACT 1's C, O, OXT).
  '''
  from mmtbx.regression import tst_rdkit_utils_molecule as M
  salt = [l for l in salt_model_str.split('\n') if l != 'END']
  nh4 = [float(x) for x in [l for l in salt if l[17:26] == 'NH4 A   3'][0][30:54].split()]
  act = dict([(l[12:16].strip(), [float(x) for x in l[30:54].split()]) for l in salt
    if l[17:26] == 'ACT A   1'])
  sites = [act['C'], act['O'], act['OXT']]
  for code, cif, atom in (('ZHZ', zhz_cif, 'N2'), ('ZHA', zha_cif, 'N1')):
    lines = [l for l in salt if l[17:26] != 'NH4 A   3'] + placed(code, cif, atom, nh4) + ['END']
    m = get_manager(get_model(lines, cifs=(cif,)), sel='resname %s' % code)
    assert m.possible_salt_bridges == [] and m.possible_groups == [], (code, possible(m))
    assert min([a.distance(b) for a in m.model.get_hierarchy().atoms() for b in
      m.model.get_hierarchy().atoms() if a.name.strip() == atom and a.parent().resname == code
      and b.name.strip() in ('OD1', 'OD2') and b.parent().parent().resseq_as_int() == 60]) < 4
  model = M.get_model(M.pdb_from_cif('ZCY', zcy_cif[1]), cifs=(('ZCY', zcy_cif[1]),))
  f = LI.find_charged_groups(model, all_atoms(model))
  assert f.groups == [] and f.possible == []
  for code, cif, pick, partners in (('ZMT', zmt_cif, ['C2', 'N2', 'N1'], []),
      ('ZTZ', ztz_cif, ['C2', 'N1', 'N4'], ['B LYS 10', 'C ARG 20'])):
    lines = [l for l in salt if l[17:26] != 'ACT A   1'] + superposed(code, cif, pick,
      sites) + ['END']
    m = get_manager(get_model(lines, cifs=(cif,)), sel='resname %s' % code)
    assert m.charged_group_failures == [], m.charged_group_failures
    assert sorted([e['residue'] for e in m.possible_salt_bridges]) == partners, (code,
      possible(m))

def exercise_not_possible():
  '''
  No possible salt bridge: His with HD1 only (usual 0) next to ACT 2; paracetamol
  (amide N, phenol) and 4-aminophenol (aniline, phenol): no possible groups; Lys 50
  with two H on NZ, 5 A from ACT 2 (beyond the 4 A atom-pair cutoff); ACT 1 with
  Lys 10 and Arg 20, all charged: salt bridges only.
  '''
  from mmtbx.regression import tst_rdkit_utils_molecule as M
  neutral = get_model(neutral_model_str.split('\n'), cifs=(acy_neutral_cif,))
  m = get_manager(neutral, sel='chain A and resseq 2')
  assert m.possible_salt_bridges == [] and salt_bridges(m) == []
  for code, cif in (('ZPA', zpa_cif), ('ZAP', zap_cif)):
    model = M.get_model(M.pdb_from_cif(code, cif[1]), cifs=((code, cif[1]),))
    f = LI.find_charged_groups(model, all_atoms(model))
    assert f.groups == [] and f.possible == [] and f.failures == [], (code, f.possible)
  salt = salt_model_str.split('\n')
  lines = [l for l in salt if not (l[17:26] == 'LYS F  50' and l[12:16].strip() == 'HZ3')]
  m = get_manager(get_model(lines), sel='chain A and resseq 2')
  assert m.possible_salt_bridges == [], possible(m)
  m = get_manager(get_model(salt), sel='chain A and resseq 1')
  assert sorted([e['residue'] for e in salt_bridges(m)]) == ['B LYS 10', 'C ARG 20']
  assert m.possible_salt_bridges == [] and m.possible_groups == []

def exercise_charge_types_and_ccd_match():
  '''
  Formal-charge records carry the restraint file's type_energy (ZAR: NT3, OC, NC1
  NC2; None for the CCD). A CCD entry whose names, elements or heavy-atom bonds do
  not match the residue (ZNX, ZAR: other compounds with these codes) gives "CCD
  entry does not match" instead of a comparison.
  '''
  from mmtbx.regression import tst_rdkit_utils_molecule as M
  zzw = [l[:21] + 'B' + l[22:30] + '%8.3f' % (float(l[30:38]) + 8.0) + l[38:]
    for l in M.pdb_from_cif('ZZW', M.zzw_cif).split('\n') if l.startswith('HETATM')]
  rows = {}
  for code, cif in (('ZAR', zar_cif), ('ZNX', znx_cif)):
    model = get_model(pdb_from_cif_text(code, cif[1])[:-1] + zzw + ['END'],
      cifs=(cif, ('zzw.cif', M.zzw_cif)))
    m = get_manager(model, sel='resname %s' % code)
    rows[code] = [(x['source'], x['atoms'], x['status'], x['type_energy'])
      for x in m.formal_charges if x['residue'] == 'A %s 1' % code]
    assert m.formal_charge_conflicts() == []
  assert rows['ZAR'] == [
    ('restraints', ['N'], 'agrees', ['NT3']),
    ('restraints', ['O', 'OXT'], 'agrees', ['OC', 'OC']),
    ('restraints', ['NE', 'NH1', 'NH2'], 'agrees', ['NC1', 'NC2', 'NC2']),
    ('CCD', ['C', 'CA', 'CB', 'CD', 'CG', 'CZ', 'N', 'NE', 'NH1', 'NH2', 'O', 'OXT'],
      'CCD entry does not match', None)], rows['ZAR']
  assert rows['ZNX'] == [
    ('restraints', ['N1'], 'agrees', ['N']),
    ('restraints', ['O1', 'O2'], 'agrees', ['O', 'O']),
    ('CCD', ['C1', 'N1', 'O1', 'O2'], 'CCD entry does not match', None)], rows['ZNX']
  # a matching CCD entry: compared (ACT, salt_sym_model_str)
  m = get_manager(get_model(salt_sym_model_str.split('\n')), sel='chain A and resseq 1')
  act = [(x['source'], x['status'], x['type_energy']) for x in m.formal_charges
    if x['residue'] == 'A ACT 1']
  assert act == [('restraints', 'agrees', ['O', 'OC']), ('CCD', 'agrees', None)], act

def exercise_builder_chain():
  '''
  SEP between two Gly with a SEP restraint file without HOP2/HOP3 (GeoStd's SEP is
  the neutral acid) and the model without them: the phosphate from the builder,
  -2, uncertain (no formal charges: total by search); the backbone is not a group
  (caps on N and C). With GeoStd's SEP and HOP2/HOP3 modelled: no group.
  '''
  from mmtbx.regression import tst_rdkit_utils_molecule as M
  from mmtbx.monomer_library import server
  import os
  model = M.get_model(M.chain_pdb)
  f = LI.find_charged_groups(model, all_atoms(model))
  assert [g for g in f.groups if g['resname'] == 'SEP'] == [] and f.failures == []
  text = open(os.path.join(server.server().geostd_path, 's', 'data_SEP.cif')).read()
  dianion = '\n'.join([l for l in text.split('\n') if 'HOP2' not in l and 'HOP3' not in l])
  pdb = '\n'.join([l for l in M.chain_pdb.split('\n') if l[12:16].strip() not in ('HOP2',
    'HOP3')])
  model = M.get_model(pdb, cifs=(('SEP', dianion),))
  f = LI.find_charged_groups(model, all_atoms(model))
  rows = builder_rows(model, f)
  assert rows == [('A SEP 2', 'phosphate', -2, 'P', ('O1P', 'O2P', 'O3P'), False, False)], rows
  (g,) = [g for g in f.groups if g['resname'] == 'SEP']
  assert g['charge_source'] == 'search' and g['hydrogens'] == 'model'

def exercise_templates_vs_builder():
  '''
  Gly-Asp-Lys-Arg-His-Gly, H complete (reduce2): the templates and the builder
  (use_templates=False) give the same side-chain groups; the builder's are
  uncertain (the amino-acid restraint files have no formal charges: total by
  search). The N-terminal Gly with H1-H3 fails in the builder (the in-chain GLY
  file has one H on N); the C-terminal Gly without OXT: no group either way (the
  builder caps its C for the absent residue).
  '''
  model = get_model(gdkrhg_model_str.split('\n'))
  atoms = model.get_hierarchy().atoms()
  t = LI.find_charged_groups(model, all_atoms(model))
  b = LI.find_charged_groups(model, all_atoms(model), use_templates=False)
  def core(groups):
    return sorted([(LI.residue_label(atoms[g['center']]), g['kind'], g['charge'],
      tuple(sorted([atoms[i].name.strip() for i in g['charged']]))) for g in groups
      if g['charge']])
  side = [r for r in core(t.groups) if r[0] != 'A GLY 1']
  assert side == [('A ARG 4', 'guanidinium', 1, ('NE', 'NH1', 'NH2')),
    ('A ASP 2', 'carboxylate', -1, ('OD1', 'OD2')), ('A HIS 5', 'imidazolium', 1, ('ND1', 'NE2')),
    ('A LYS 3', 'ammonium', 1, ('NZ',))], side
  assert core(b.groups) == side, core(b.groups)
  assert [r for r in core(t.groups) if r[0] == 'A GLY 1'] == [('A GLY 1', 'ammonium', 1,
    ('N',))]
  assert set([(g['source'], g['certain'], g['charge_source']) for g in b.groups]) == set(
    [('builder', False, 'search')])
  assert set([(g['source'], g['certain']) for g in t.groups]) == set([('template', True)])
  assert [(x['residue'], x['reason']) for x in b.failures] == [('A GLY 1',
    'H not in the restraint file: N: 3 H in the model, 1 in the restraint file')]
  assert t.failures == [] and t.missing == []
  # the C-terminal Gly without OXT: the absent following residue capped (not an
  # acylium cation), neutral
  from mmtbx.ligands import rdkit_utils
  rg = [x for x in model.get_hierarchy().residue_groups() if x.resseq_as_int() == 6][0]
  r = rdkit_utils.residue_molecule(model, rg)
  assert r.ok and r.total_charge == 0, r.reason
  assert sorted([(c['kind'], atoms[c['on']].name.strip(), c['partner'] if c['kind'] ==
    'missing' else atoms[c['partner']].name.strip()) for c in r.caps]) == [
    ('linked', 'N', 'C'), ('missing', 'C', '(no following residue)')], r.caps

def exercise_nucleotides():
  '''
  Nucleotides by name: DG 2 in the chain -1 (OP1, OP2); DG with OP3 (5' phosphate)
  -2; with HOP3 on OP3 -1; without H the charge at pH 7 ("assumed (no H)"); DC 1
  without P: no group.
  '''
  lines = [l for l in dna_model_str.split('\n') if l.startswith('ATOM')]
  def rows(lines):
    model = get_model(['CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1'] + lines +
      ['END'])
    return group_rows(model, LI.find_charged_groups(model, all_atoms(model)).groups)
  assert rows(lines) == [('A DG 2', 'phosphate', -1, -1, 'modelled', '', 'template',
    ('OP1', 'OP2'))], rows(lines)
  o3 = [l for l in lines if l[12:16] == " O3'" and l[22:26] == '   1'][0]
  dg = [l for l in lines if l[22:26] == '   2']
  op3 = o3[:12] + ' OP3  DG A   2' + o3[26:]
  assert rows(dg + [op3]) == [('A DG 2', 'phosphate', -2, -2, 'modelled', '', 'template',
    ('OP1', 'OP2', 'OP3'))]
  hop3 = op3[:12] + 'HOP3' + op3[16:30] + '%8.3f' % (float(op3[30:38]) + 0.96) + \
    op3[38:76] + ' H'
  assert rows(dg + [op3, hop3])[0][2:5] == (-1, -2, 'modelled')
  no_h = [l for l in lines if l[76:78].strip() != 'H']
  assert rows(no_h) == [('A DG 2', 'phosphate', -1, -1, 'assumed (no H)', '', 'template',
    ('OP1', 'OP2'))]
  assert rows([l for l in dg if l[76:78].strip() != 'H'] + [op3])[0][2:5] == (-2, -2,
    'assumed (no H)')

def exercise_termini():
  '''
  Without H, the N-terminus is assumed +1 only for the first residue of its chain:
  LYS 50 (no H) of templates_model_str as residue 1 and 5 of a chain: residue 1
  only. Neither has OXT: chain breaks, no C-terminal group and no missing OXT.
  '''
  lys = [l for l in templates_model_str.split('\n') if l[17:26] == 'LYS B  50']
  first = [l[:21] + 'X   1' + l[26:] for l in lys]
  fifth = [l[:21] + 'X   5' + l[26:30] + '%8.3f' % (float(l[30:38]) + 20) + l[38:] for l in lys]
  model = get_model(['CRYST1  200.000  200.000  200.000  90.00  90.00  90.00 P 1'] + first +
    fifth + ['END'])
  f = LI.find_charged_groups(model, all_atoms(model))
  assert group_rows(model, f.groups) == [
    ('X LYS 1', 'ammonium', 1, 1, 'assumed (no H)', '', 'template', ('N',)),
    ('X LYS 1', 'ammonium', 1, 1, 'assumed (no H)', '', 'template', ('NZ',)),
    ('X LYS 5', 'ammonium', 1, 1, 'assumed (no H)', '', 'template', ('NZ',))], \
    group_rows(model, f.groups)
  assert f.missing == []

def all_atoms(model):
  return flex.bool(model.get_number_of_atoms(), True)

def group_rows(model, groups, charged_only=False):
  atoms = model.get_hierarchy().atoms()
  return sorted([(LI.residue_label(atoms[g['center']]), g['kind'], g['charge'],
    g['usual_charge'], g['state'], g['altloc'], g['source'],
    tuple(sorted([atoms[i].name.strip() for i in g['charged']])))
    for g in groups if g['charge'] or not charged_only])

def exercise_charged_groups():
  '''
  Every group kind on the salt fixtures (charged groups): the standard residues by
  template, ACT and NH4 from the builder (templates vs builder:
  exercise_templates_vs_builder).
  '''
  model = get_model(salt_model_str.split('\n'))
  atoms = model.get_hierarchy().atoms()
  found = LI.find_charged_groups(model, all_atoms(model))
  kinds = sorted([(g['kind'], g['charge'], g['resname']) for g in found.groups if g['charge']])
  assert kinds == [('ammonium', 1, 'LYS'), ('ammonium', 1, 'LYS'), ('ammonium', 1, 'NH4'),
    ('carboxylate', -1, 'ACT'), ('carboxylate', -1, 'ACT'), ('carboxylate', -1, 'ASP'),
    ('carboxylate', -1, 'ASP'), ('guanidinium', 1, 'ARG'), ('imidazolium', 1, 'HIS')], kinds
  assert set([g['source'] for g in found.groups if g['resname'] in ('LYS', 'ARG', 'HIS',
    'ASP')]) == set(['template'])
  builder = [g for g in found.groups if g['resname'] in ('ACT', 'NH4')]
  assert sorted(set([(g['resname'], g['source'], g['certain'], g['usual_charge'],
    g['charge_source']) for g in builder])) == [('ACT', 'builder', True, None,
    'formal charges'), ('NH4', 'builder', True, None, 'restraint file')], builder
  assert found.failures == [] and found.dropped == []

def exercise_amino_acid_templates():
  '''
  Charge as modelled from the H on the group atoms, the charge at pH 7 recorded:
  Asp with HD2 neutral, His with HD1 only neutral and with HD1 and HE2 charged,
  Lys with two H on NZ neutral, a residue without H assumed at pH 7 and flagged,
  Arg complete charged; the N-terminal N (two H here) neutral; without H (Lys 50,
  not the chain's first residue) no N-terminal group; no OXT: chain breaks, no
  C-terminal group and nothing reported.
  '''
  model = get_model(templates_model_str.split('\n'))
  found = LI.find_charged_groups(model, all_atoms(model))
  rows = [r for r in group_rows(model, found.groups) if r[1] != 'ammonium' or
    r[7] != ('N',)]
  assert rows == [
    ('B ARG 60', 'guanidinium', 1, 1, 'modelled', '', 'template', ('NE', 'NH1', 'NH2')),
    ('B ASP 10', 'carboxylate', 0, -1, 'modelled', '', 'template', ('OD1', 'OD2')),
    ('B HIS 20', 'imidazolium', 0, 0, 'modelled', '', 'template', ('ND1', 'NE2')),
    ('B HIS 30', 'imidazolium', 1, 0, 'modelled', '', 'template', ('ND1', 'NE2')),
    ('B LYS 40', 'ammonium', 0, 1, 'modelled', '', 'template', ('NZ',)),
    ('B LYS 50', 'ammonium', 1, 1, 'assumed (no H)', '', 'template', ('NZ',))], rows
  nterm = [r for r in group_rows(model, found.groups) if r[7] == ('N',)]
  assert [(r[0], r[2], r[4]) for r in nterm] == [('B ARG 60', 0, 'modelled'),
    ('B ASP 10', 0, 'modelled'), ('B HIS 20', 0, 'modelled'), ('B HIS 30', 0, 'modelled'),
    ('B LYS 40', 0, 'modelled')], nterm
  assert found.missing == []

def exercise_conformers():
  '''
  Per conformer: Asp/Asn microheterogeneity gives the carboxylate for Asp (A)
  only (the residue split into two residue groups: the N-terminal N of each
  conformer); an Asp split into A and B next to a blank-altloc ligand ammonium makes one
  salt bridge per conformer; with the ligand split too, A pairs only with A and B
  with B (all four distances within 4 A).
  '''
  model = get_model(microhet_model_str.split('\n'))
  found = LI.find_charged_groups(model, all_atoms(model))
  rows = [r for r in group_rows(model, found.groups) if r[7] != ('N',)]
  assert rows == [('B ASP 10', 'carboxylate', -1, -1, 'modelled', 'A', 'template',
    ('OD1', 'OD2'))], rows
  assert [g['resname'] for g in found.groups if g['altloc'] == 'B'] == ['ASN']
  assert found.missing == []
  atoms = model.get_hierarchy().atoms()
  assert [LI.residue_label(atoms[g['center']], g['resname']) for g in found.groups
    if g['altloc'] == 'B'] == ['B ASN 10']
  model = get_model(altloc_asp_model_str.split('\n'))
  m = get_manager(model, sel='chain A and resseq 3')
  assert [(g['kind'], g['altloc']) for g in m.charged_groups] == [('ammonium', '')]
  sb = [(e['geometry']['charged_groups']['ligand_group']['altloc'],
    e['geometry']['charged_groups']['partner_group']['altloc'], e['labels'][1])
    for e in salt_bridges(m)]
  assert sb == [('', 'A', 'G ASP 60 OD1 alt A'), ('', 'B', 'G ASP 60 OD1 alt B')], sb
  model = get_model(altloc_both_model_str.split('\n'))
  atoms = model.get_hierarchy().atoms()
  n = dict([(a.parent().altloc, a) for a in atoms if a.name.strip() == 'N' and
    a.parent().resname == 'NH4'])
  o = [a for a in atoms if a.parent().resname == 'ASP' and a.name.strip() in ('OD1', 'OD2')]
  assert max([n['A'].distance(x) for x in o if x.parent().altloc == 'B']) < 4.0
  m = get_manager(model, sel='chain A and resseq 3')
  sb = sorted([(e['geometry']['charged_groups']['ligand_group']['altloc'],
    e['geometry']['charged_groups']['partner_group']['altloc']) for e in salt_bridges(m)])
  assert sb == [('A', 'A'), ('B', 'B')], sb

def exercise_group_scope():
  '''
  Groups are examined on the ligand and the residues within the search distance
  only: for NH4 A 3 the charged residues of the other two sites (> 20 A) are not
  examined.
  '''
  model = get_model(salt_model_str.split('\n'))
  m = get_manager(model, sel='chain A and resseq 3')
  residues = sorted(set([g['residue'] for g in m.examined_groups]))
  assert residues == ['A NH4 3', 'G ASP 60'], residues
  atoms = model.get_hierarchy().atoms()
  lys = [a for a in atoms if a.parent().resname == 'LYS' and a.name.strip() == 'NZ']
  nh4 = [a for a in atoms if a.parent().resname == 'NH4' and a.name.strip() == 'N'][0]
  assert min([nh4.distance(a) for a in lys]) > 20

def exercise_errors(model):
  isel = model.selection(LIG_SEL).iselection()
  try:
    LI.manager(model, isel, 'chain A and resseq 2').run()
  except Sorry:
    pass
  else:
    raise AssertionError('mismatched sel_str accepted')
  no_h = model.select(~model.get_hierarchy().atom_selection_cache().selection(
    'element H'))
  try:
    LI.manager(no_h, no_h.selection(LIG_SEL).iselection(), LIG_SEL)
  except Sorry:
    pass
  else:
    raise AssertionError('model without H accepted')
  # one molecule: two separate EDO residues refused
  sel = 'chain A and resname EDO and (resseq 1 or resseq 3)'
  try:
    LI.manager(model, model.selection(sel).iselection(), sel).run()
  except Sorry as e:
    assert 'not one molecule' in str(e) and '2 fragments' in str(e), str(e)
  else:
    raise AssertionError('two molecules accepted')

# ------------------------------------------------------------------------------

def run():
  model = get_model()
  exercise_dot_patches()
  exercise_probe_names()
  m = exercise_entries_and_counts(model)
  exercise_patches(m)
  exercise_contact_area(m)
  exercise_probe_classes(model)
  exercise_inline_clash()
  exercise_hbond_phil()
  exercise_charged_groups()
  exercise_amino_acid_templates()
  exercise_conformers()
  exercise_builder_groups()
  exercise_charge_separated()
  exercise_atomic_charge_comparison()
  exercise_possible_salt_bridges()
  exercise_possible_kinds()
  exercise_amine_and_tetrazole_rules()
  exercise_not_possible()
  exercise_charge_types_and_ccd_match()
  exercise_builder_chain()
  exercise_templates_vs_builder()
  exercise_nucleotides()
  exercise_termini()
  exercise_group_scope()
  exercise_salt_bridges()
  exercise_salt_bridge_symmetry()
  exercise_formal_charge_conflict()
  exercise_pair_class_order(model)
  exercise_internal()
  exercise_symmetry()
  exercise_donor_conformers()
  exercise_probe_names_in_run()
  exercise_library_vs_command_line(model)
  exercise_ligand_overlaps(model)
  exercise_imports()
  exercise_errors(model)
  print('OK')

if __name__ == '__main__':
  run()
