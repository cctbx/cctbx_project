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
  Salt bridges on inline models (ideal CCD/GeoStd geometry; charges from the H):
    ACT 1 with Lys 10 (NZ...OXT 2.85 A, NZ-HZ1...OXT H-bond) and Arg 20;
    ACT 2 with His 30 (HD1 and HE2), Asp 40 (same charge), Lys 50 (5.0 A);
    NH4 3 with Asp 60; the neutral model: acetic acid (ACY, custom restraints)
    at ACT 1's place, His 30 with HD1 only.
  '''
  model = get_model(salt_model_str.split('\n'))
  m1 = get_manager(model, sel='chain A and resseq 1')
  assert m1.charged_groups == [dict(kind='carboxylate', charge=-1, usual_charge=None,
    state='modelled', source='rules', altloc='', atoms=['A ACT 1 O', 'A ACT 1 OXT'],
    residue='A ACT 1')], m1.charged_groups
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
  the model): the perceived -1 conflicts with the dictionary's 0 and is reported.
  '''
  model = get_model(salt_sym_model_str.split('\n'), cifs=(act_conflict_cif,))
  m = get_manager(model, sel='chain A and resseq 1')
  (c,) = m.formal_charge_conflicts()
  assert (c['residue'], c['source'], c['group'], c['atoms'], c['dictionary_charge'],
    c['perceived_charge'], c['file']) == ('A ACT 1', 'restraints', 'carboxylate',
    ['C', 'O', 'OXT'], 0, -1, 'act_conflict.cif'), c
  assert [x['status'] for x in m.formal_charges if x['residue'] == 'A ACT 1' and
    x['source'] == 'CCD'] == ['agrees']
  assert len(salt_bridges(m)) == 1
  log = StringIO()
  m.show(log=log)
  assert 'formal-charge conflicts' in log.getvalue()

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
  Every group kind on the salt fixtures (charged groups); for H-complete,
  altloc-free standard residues the templates and the rules give the same groups.
  '''
  model = get_model(salt_model_str.split('\n'))
  atoms = model.get_hierarchy().atoms()
  found = LI.find_charged_groups(model, all_atoms(model))
  kinds = sorted([(g['kind'], g['charge'], g['resname']) for g in found.groups if g['charge']])
  assert kinds == [('ammonium', 1, 'LYS'), ('ammonium', 1, 'LYS'), ('ammonium', 1, 'NH4'),
    ('carboxylate', -1, 'ACT'), ('carboxylate', -1, 'ACT'), ('carboxylate', -1, 'ASP'),
    ('carboxylate', -1, 'ASP'), ('guanidinium', 1, 'ARG'), ('imidazolium', 1, 'HIS')], kinds
  rules = LI.find_charged_groups(model, all_atoms(model), use_templates=False)
  def core(groups):
    return sorted([(g['kind'], g['charge'], g['center'], tuple(sorted(g['charged'])))
      for g in groups if g['charge']])
  assert core(found.groups) == core(rules.groups)
  assert set([g['source'] for g in found.groups if g['resname'] in ('LYS', 'ARG', 'HIS',
    'ASP')]) == set(['template'])
  assert set([g['source'] for g in rules.groups]) == set(['rules'])
  assert set([g['usual_charge'] for g in rules.groups]) == set([None])

def exercise_amino_acid_templates():
  '''
  Charge as modelled from the H on the group atoms, the charge at pH 7 recorded:
  Asp with HD2 neutral, His with HD1 only neutral and with HD1 and HE2 charged,
  Lys with two H on NZ neutral, a residue without H assumed at pH 7 and flagged,
  Arg complete charged; the N-terminal N (two H here) neutral; a missing OXT is
  reported and the C-terminal group skipped.
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
    ('B LYS 40', 0, 'modelled'), ('B LYS 50', 1, 'assumed (no H)')], nterm
  assert sorted([(x['residue'], x['kind'], x['atoms']) for x in found.missing]) == [
    (r, 'carboxylate', ['OXT']) for r in ('B ARG 60', 'B ASP 10', 'B HIS 20', 'B HIS 30',
    'B LYS 40', 'B LYS 50')]

def exercise_conformers():
  '''
  Per conformer: Asp/Asn microheterogeneity gives the carboxylate for Asp (A)
  only; an Asp split into A and B next to a blank-altloc ligand ammonium makes one
  salt bridge per conformer; with the ligand split too, A pairs only with A and B
  with B (all four distances within 4 A).
  '''
  model = get_model(microhet_model_str.split('\n'))
  found = LI.find_charged_groups(model, all_atoms(model))
  rows = [r for r in group_rows(model, found.groups) if r[7] != ('N',)]
  assert rows == [('B ASP 10', 'carboxylate', -1, -1, 'modelled', 'A', 'template',
    ('OD1', 'OD2'))], rows
  assert [g['resname'] for g in found.groups if g['altloc'] == 'B'] == ['ASN']
  assert sorted([(x['residue'], x['altloc']) for x in found.missing]) == [
    ('B ASN 10', 'B'), ('B ASP 10', 'A')]
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
