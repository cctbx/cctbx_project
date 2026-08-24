"""Regression test: ``mmtbx.development.endo_exo`` produces a stable
QM region around two contrasting Fe sites:

* **1BQ8** (Pyrococcus rubredoxin, P 21 21 21).  All four Cys ligands
  live inside the ASU, so symmetry plays no role.  Exercises the
  baseline BFS + capping pipeline.

* **2C2U** (Halobacterium ferritin core, P 2 3).  The Fe sits on a
  3-fold special position; three symmetry-related copies of Asp 93 and
  HOH 2154 from the ASU map onto the metal.  Exercises:

  - symmetry-aware adjacency (the BFS picks up the Asp 93 / HOH 2154
    symmetry images even though the parent atoms are far from the
    metal in the ASU);
  - special-position deduplication (three sym_ops put Fe at the same
    physical point; only one Fe atom survives in the materialized
    region);
  - per-sym_op chain IDs (the identity copy stays in chain "A"; the
    two extra images land in single-character chain IDs).
"""

from __future__ import absolute_import, division, print_function

import io
import random
from collections import defaultdict, deque

import iotbx.pdb
import libtbx.phil
from libtbx.test_utils import approx_equal
from libtbx.utils import Sorry, format_cpu_times

import mmtbx.model
from iotbx.data_manager import DataManager
from mmtbx.programs.endo_exo import Program as EndoexoProgram
from scitbx import matrix
from cctbx import geometry_restraints, sgtbx

from mmtbx.geometry_restraints.endo_exo.util import _canon_op
from mmtbx.geometry_restraints.endo_exo.capping import HydrogenCapper
from mmtbx.geometry_restraints.endo_exo.cutting import (
  BondCutDetector, PREFERRED_CUTS)
from mmtbx.geometry_restraints.endo_exo.graph import AtomGraphBuilder
from mmtbx.geometry_restraints.endo_exo.grow import QMRegionGrower


# 8 A sphere around the Fe of 1BQ8 (29 residues, 154 atoms).  Slightly
# larger than the QM region endo_exo extracts (72 atoms) so the BFS has
# scaffold to cap into.
_1BQ8_FE_SPHERE_PDB = """\
CRYST1   33.823   34.705   43.204  90.00  90.00  90.00 P 21 21 21
ATOM      1  N   VAL A   5      16.394   5.495   0.266  1.00  3.20           N
ATOM      2  CA  VAL A   5      15.980   6.892   0.493  1.00  3.49           C
ATOM      3  C   VAL A   5      14.480   7.032   0.278  1.00  2.98           C
ATOM      4  O   VAL A   5      13.922   6.570  -0.728  1.00  3.30           O
ATOM      5  CB  VAL A   5      16.763   7.863  -0.366  1.00  3.68           C
ATOM      6  CG1 VAL A   5      16.514   7.724  -1.870  1.00  4.62           C
ATOM      7  CG2 VAL A   5      16.523   9.338   0.046  1.00  5.32           C
ATOM      8  H   VAL A   5      16.393   5.260  -0.561  1.00  3.20           H
ATOM      9  HA  VAL A   5      16.173   7.127   1.414  1.00  3.49           H
ATOM     10  HB  VAL A   5      17.689   7.625  -0.205  1.00  3.68           H
ATOM     11 HG11 VAL A   5      15.573   7.877  -2.050  1.00  4.62           H
ATOM     12 HG12 VAL A   5      16.763   6.829  -2.150  1.00  4.62           H
ATOM     13 HG13 VAL A   5      17.052   8.380  -2.341  1.00  4.62           H
ATOM     14 HG21 VAL A   5      16.798   9.457   0.968  1.00  5.32           H
ATOM     15 HG22 VAL A   5      15.579   9.542  -0.046  1.00  5.32           H
ATOM     16 HG23 VAL A   5      17.045   9.916  -0.532  1.00  5.32           H
ATOM     17  N   CYS A   6      13.876   7.840   1.115  1.00  2.56           N
ATOM     18  CA  CYS A   6      12.515   8.342   0.910  1.00  2.77           C
ATOM     19  C   CYS A   6      12.644   9.458  -0.142  1.00  2.56           C
ATOM     20  O   CYS A   6      13.271  10.473   0.122  1.00  3.23           O
ATOM     21  CB  CYS A   6      11.937   8.897   2.187  1.00  2.83           C
ATOM     22  SG  CYS A   6      10.256   9.599   1.876  1.00  3.04           S
ATOM     23  H   CYS A   6      14.236   8.131   1.840  1.00  2.56           H
ATOM     24  HA  CYS A   6      11.906   7.642   0.627  1.00  2.77           H
ATOM     25  HB2 CYS A   6      12.511   9.602   2.524  1.00  2.83           H
ATOM     26  HB3 CYS A   6      11.861   8.188   2.845  1.00  2.83           H
ATOM     28  N   LYS A   7      12.120   9.227  -1.346  1.00  2.98           N
ATOM     29  CA  LYS A   7      12.260  10.201  -2.421  1.00  3.51           C
ATOM     30  C   LYS A   7      11.541  11.504  -2.103  1.00  3.88           C
ATOM     31  O   LYS A   7      11.857  12.540  -2.672  1.00  5.65           O
ATOM     32  CB  LYS A   7      11.742   9.650  -3.758  1.00  5.11           C
ATOM     33  CG  LYS A   7      12.524   8.451  -4.265  1.00  6.96           C
ATOM     34  CD  LYS A   7      11.938   7.989  -5.583  1.00  8.84           C
ATOM     35  CE  LYS A   7      10.500   7.483  -5.529  1.00  9.73           C
ATOM     36  NZ  LYS A   7      10.079   6.961  -6.833  1.00 11.95           N
ATOM     37  H   LYS A   7      11.683   8.518  -1.561  1.00  2.98           H
ATOM     38  HA  LYS A   7      13.207  10.388  -2.518  1.00  3.51           H
ATOM     39  HB2 LYS A   7      11.802  10.349  -4.428  1.00  5.11           H
ATOM     40  HB3 LYS A   7      10.818   9.376  -3.646  1.00  5.11           H
ATOM     41  HG2 LYS A   7      13.452   8.698  -4.404  1.00  6.96           H
ATOM     42  HG3 LYS A   7      12.466   7.725  -3.624  1.00  6.96           H
ATOM     43  HD2 LYS A   7      11.956   8.735  -6.203  1.00  8.84           H
ATOM     44  HD3 LYS A   7      12.485   7.263  -5.922  1.00  8.84           H
ATOM     45  HE2 LYS A   7       9.909   8.213  -5.285  1.00  9.73           H
ATOM     46  HE3 LYS A   7      10.433   6.769  -4.876  1.00  9.73           H
ATOM     47  HZ1 LYS A   7      10.606   6.286  -7.075  1.00 11.95           H
ATOM     48  HZ2 LYS A   7      10.129   7.602  -7.448  1.00 11.95           H
ATOM     49  HZ3 LYS A   7       9.239   6.669  -6.788  1.00 11.95           H
ATOM     50  N   ILE A   8      10.572  11.487  -1.162  1.00  3.53           N
ATOM     51  CA  ILE A   8       9.858  12.717  -0.799  1.00  3.44           C
ATOM     52  C   ILE A   8      10.661  13.609   0.123  1.00  3.26           C
ATOM     53  O   ILE A   8      10.797  14.826  -0.111  1.00  5.13           O
ATOM     54  CB  ILE A   8       8.487  12.375  -0.162  1.00  5.11           C
ATOM     55  CG1 ILE A   8       7.671  11.360  -0.956  1.00  5.47           C
ATOM     56  CG2 ILE A   8       7.719  13.652   0.102  1.00  5.65           C
ATOM     57  CD1 ILE A   8       7.349  11.801  -2.358  1.00  8.43           C
ATOM     58  H   ILE A   8      10.320  10.786  -0.732  1.00  3.53           H
ATOM     59  HA  ILE A   8       9.719  13.222  -1.616  1.00  3.44           H
ATOM     60  HB  ILE A   8       8.666  11.930   0.681  1.00  5.11           H
ATOM     61 HG12 ILE A   8       8.174  10.533  -1.015  1.00  5.47           H
ATOM     62 HG13 ILE A   8       6.832  11.206  -0.495  1.00  5.47           H
ATOM     63 HG21 ILE A   8       7.582  14.117  -0.738  1.00  5.65           H
ATOM     64 HG22 ILE A   8       8.231  14.209   0.709  1.00  5.65           H
ATOM     65 HG23 ILE A   8       6.863  13.429   0.500  1.00  5.65           H
ATOM     66 HD11 ILE A   8       8.177  11.950  -2.840  1.00  8.43           H
ATOM     67 HD12 ILE A   8       6.835  12.623  -2.320  1.00  8.43           H
ATOM     68 HD13 ILE A   8       6.832  11.107  -2.796  1.00  8.43           H
ATOM     69  N   CYS A   9      11.210  13.074   1.223  1.00  3.35           N
ATOM     70  CA  CYS A   9      11.769  13.851   2.272  1.00  3.49           C
ATOM     71  C   CYS A   9      13.219  13.660   2.581  1.00  3.53           C
ATOM     72  O   CYS A   9      13.796  14.406   3.390  1.00  3.48           O
ATOM     73  CB  CYS A   9      10.948  13.686   3.564  1.00  3.16           C
ATOM     74  SG  CYS A   9      11.311  12.173   4.470  1.00  3.34           S
ATOM     75  H   CYS A   9      11.261  12.228   1.367  1.00  3.35           H
ATOM     76  HA  CYS A   9      11.723  14.760   1.938  1.00  3.49           H
ATOM     77  HB2 CYS A   9      11.137  14.435   4.151  1.00  3.16           H
ATOM     78  HB3 CYS A   9      10.006  13.674   3.335  1.00  3.16           H
ATOM     80  N   GLY A  10      13.879  12.653   2.014  1.00  2.95           N
ATOM     81  CA  GLY A  10      15.282  12.416   2.247  1.00  2.72           C
ATOM     82  C   GLY A  10      15.614  11.546   3.426  1.00  2.73           C
ATOM     83  O   GLY A  10      16.820  11.254   3.626  1.00  3.00           O
ATOM     84  H   GLY A  10      13.520  12.084   1.479  1.00  2.95           H
ATOM     85  HA2 GLY A  10      15.719  13.271   2.384  1.00  2.72           H
ATOM     86  HA3 GLY A  10      15.657  11.993   1.459  1.00  2.72           H
ATOM     87  N   TYR A  11      14.667  11.081   4.202  1.00  2.81           N
ATOM     88  CA  TYR A  11      14.915  10.112   5.274  1.00  2.91           C
ATOM     89  C   TYR A  11      15.553   8.845   4.649  1.00  2.90           C
ATOM     90  O   TYR A  11      15.159   8.433   3.551  1.00  3.10           O
ATOM     91  CB  TYR A  11      13.580   9.729   5.930  1.00  2.88           C
ATOM     92  CG  TYR A  11      13.703   8.536   6.859  1.00  2.57           C
ATOM     93  CD1 TYR A  11      14.232   8.652   8.139  1.00  3.41           C
ATOM     94  CD2 TYR A  11      13.305   7.275   6.410  1.00  3.49           C
ATOM     95  CE1 TYR A  11      14.336   7.547   8.970  1.00  4.29           C
ATOM     96  CE2 TYR A  11      13.427   6.178   7.235  1.00  3.09           C
ATOM     97  CZ  TYR A  11      13.946   6.307   8.493  1.00  3.73           C
ATOM     98  OH  TYR A  11      14.053   5.156   9.266  1.00  4.47           O
ATOM     99  H   TYR A  11      13.841  11.310   4.137  1.00  2.81           H
ATOM    100  HA  TYR A  11      15.502  10.482   5.952  1.00  2.91           H
ATOM    101  HB2 TYR A  11      13.256  10.482   6.449  1.00  2.88           H
ATOM    102  HB3 TYR A  11      12.940   9.505   5.236  1.00  2.88           H
ATOM    103  HD1 TYR A  11      14.520   9.483   8.442  1.00  3.41           H
ATOM    104  HD2 TYR A  11      12.957   7.175   5.553  1.00  3.49           H
ATOM    105  HE1 TYR A  11      14.663   7.637   9.836  1.00  4.29           H
ATOM    106  HE2 TYR A  11      13.154   5.342   6.933  1.00  3.09           H
ATOM    107  HH  TYR A  11      13.383   5.083   9.767  1.00  4.47           H
ATOM    108  N   ILE A  12      16.530   8.325   5.360  1.00  2.50           N
ATOM    109  CA  ILE A  12      17.176   7.071   4.961  1.00  2.77           C
ATOM    110  C   ILE A  12      16.631   5.931   5.778  1.00  2.84           C
ATOM    111  O   ILE A  12      16.815   5.923   6.989  1.00  3.43           O
ATOM    112  CB  ILE A  12      18.702   7.193   5.103  1.00  3.75           C
ATOM    113  CG1 ILE A  12      19.290   8.421   4.373  1.00  4.12           C
ATOM    114  CG2 ILE A  12      19.388   5.922   4.674  1.00  4.75           C
ATOM    115  CD1 ILE A  12      19.194   8.339   2.859  1.00  6.27           C
ATOM    116  H   ILE A  12      16.846   8.672   6.081  1.00  2.50           H
ATOM    117  HA  ILE A  12      16.980   6.882   4.030  1.00  2.77           H
ATOM    118  HB  ILE A  12      18.873   7.336   6.047  1.00  3.75           H
ATOM    119 HG12 ILE A  12      18.808   9.213   4.657  1.00  4.12           H
ATOM    120 HG13 ILE A  12      20.228   8.502   4.606  1.00  4.12           H
ATOM    121 HG21 ILE A  12      19.171   5.744   3.745  1.00  4.75           H
ATOM    122 HG22 ILE A  12      19.077   5.192   5.232  1.00  4.75           H
ATOM    123 HG23 ILE A  12      20.347   6.031   4.775  1.00  4.75           H
ATOM    124 HD11 ILE A  12      18.261   8.268   2.605  1.00  6.27           H
ATOM    125 HD12 ILE A  12      19.680   7.557   2.554  1.00  6.27           H
ATOM    126 HD13 ILE A  12      19.581   9.141   2.474  1.00  6.27           H
ATOM    127  N   TRP A  37       3.798  -0.020   6.436  1.00  4.44           N
ATOM    128  CA  TRP A  37       4.860   0.808   7.049  1.00  3.68           C
ATOM    129  C   TRP A  37       5.032   2.052   6.206  1.00  3.09           C
ATOM    130  O   TRP A  37       4.931   2.011   5.003  1.00  4.07           O
ATOM    131  CB  TRP A  37       6.150  -0.001   7.149  1.00  3.69           C
ATOM    132  CG  TRP A  37       7.317   0.766   7.681  1.00  3.57           C
ATOM    133  CD1 TRP A  37       7.768   0.734   8.982  1.00  3.65           C
ATOM    134  CD2 TRP A  37       8.215   1.617   6.994  1.00  3.73           C
ATOM    135  NE1 TRP A  37       8.896   1.532   9.118  1.00  3.98           N
ATOM    136  CE2 TRP A  37       9.187   2.061   7.896  1.00  3.73           C
ATOM    137  CE3 TRP A  37       8.319   2.027   5.639  1.00  3.78           C
ATOM    138  CZ2 TRP A  37      10.240   2.931   7.536  1.00  3.88           C
ATOM    139  CZ3 TRP A  37       9.358   2.853   5.292  1.00  4.00           C
ATOM    140  CH2 TRP A  37      10.308   3.297   6.238  1.00  4.03           C
ATOM    141  H   TRP A  37       3.998  -0.854   6.372  1.00  4.44           H
ATOM    142  HA  TRP A  37       4.637   1.078   7.954  1.00  3.68           H
ATOM    143  HB2 TRP A  37       6.386  -0.318   6.263  1.00  3.69           H
ATOM    144  HB3 TRP A  37       5.999  -0.754   7.742  1.00  3.69           H
ATOM    145  HD1 TRP A  37       7.374   0.248   9.670  1.00  3.65           H
ATOM    146  HE1 TRP A  37       9.334   1.669   9.845  1.00  3.98           H
ATOM    147  HE3 TRP A  37       7.701   1.744   5.005  1.00  3.78           H
ATOM    148  HZ2 TRP A  37      10.856   3.237   8.162  1.00  3.88           H
ATOM    149  HZ3 TRP A  37       9.439   3.129   4.408  1.00  4.00           H
ATOM    150  HH2 TRP A  37      10.998   3.855   5.960  1.00  4.03           H
ATOM    151  N   VAL A  38       5.270   3.185   6.882  1.00  2.88           N
ATOM    152  CA  VAL A  38       5.473   4.448   6.200  1.00  3.57           C
ATOM    153  C   VAL A  38       6.750   5.130   6.648  1.00  3.48           C
ATOM    154  O   VAL A  38       7.242   4.891   7.767  1.00  3.30           O
ATOM    155  CB  VAL A  38       4.286   5.394   6.306  1.00  4.43           C
ATOM    156  CG1 VAL A  38       3.017   4.799   5.695  1.00  5.23           C
ATOM    157  CG2 VAL A  38       4.031   5.795   7.745  1.00  5.51           C
ATOM    158  H   VAL A  38       5.317   3.238   7.739  1.00  2.88           H
ATOM    159  HA  VAL A  38       5.567   4.231   5.259  1.00  3.57           H
ATOM    160  HB  VAL A  38       4.511   6.189   5.798  1.00  4.43           H
ATOM    161 HG11 VAL A  38       2.796   3.978   6.163  1.00  5.23           H
ATOM    162 HG12 VAL A  38       3.176   4.612   4.756  1.00  5.23           H
ATOM    163 HG13 VAL A  38       2.293   5.438   5.787  1.00  5.23           H
ATOM    164 HG21 VAL A  38       4.819   6.241   8.093  1.00  5.51           H
ATOM    165 HG22 VAL A  38       3.844   4.999   8.266  1.00  5.51           H
ATOM    166 HG23 VAL A  38       3.270   6.396   7.775  1.00  5.51           H
ATOM    167  N   CYS A  39       7.214   6.059   5.822  1.00  3.08           N
ATOM    168  CA  CYS A  39       8.351   6.911   6.228  1.00  3.01           C
ATOM    169  C   CYS A  39       8.015   7.538   7.573  1.00  2.69           C
ATOM    170  O   CYS A  39       6.957   8.180   7.696  1.00  2.90           O
ATOM    171  CB  CYS A  39       8.518   8.006   5.167  1.00  2.63           C
ATOM    172  SG  CYS A  39       9.909   9.101   5.614  1.00  3.02           S
ATOM    173  H   CYS A  39       6.902   6.220   5.037  1.00  3.08           H
ATOM    174  HA  CYS A  39       9.180   6.414   6.305  1.00  3.01           H
ATOM    175  HB2 CYS A  39       7.707   8.536   5.115  1.00  2.63           H
ATOM    176  HB3 CYS A  39       8.705   7.600   4.306  1.00  2.63           H
ATOM    178  N   PRO A  40       8.872   7.421   8.563  1.00  3.42           N
ATOM    179  CA  PRO A  40       8.561   7.954   9.916  1.00  3.79           C
ATOM    180  C   PRO A  40       8.613   9.480   9.954  1.00  3.62           C
ATOM    181  O   PRO A  40       8.148  10.060  10.929  1.00  5.67           O
ATOM    182  CB  PRO A  40       9.655   7.382  10.800  1.00  4.40           C
ATOM    183  CG  PRO A  40      10.750   7.051   9.885  1.00  5.10           C
ATOM    184  CD  PRO A  40      10.139   6.613   8.599  1.00  3.49           C
ATOM    185  HA  PRO A  40       7.697   7.631  10.217  1.00  3.79           H
ATOM    186  HB2 PRO A  40       9.935   8.046  11.450  1.00  4.40           H
ATOM    187  HB3 PRO A  40       9.331   6.589  11.254  1.00  4.40           H
ATOM    188  HG2 PRO A  40      11.301   7.837   9.747  1.00  5.10           H
ATOM    189  HG3 PRO A  40      11.282   6.336  10.267  1.00  5.10           H
ATOM    190  HD2 PRO A  40      10.713   6.829   7.847  1.00  3.49           H
ATOM    191  HD3 PRO A  40       9.948   5.662   8.608  1.00  3.49           H
ATOM    192  N   ILE A  41       9.158  10.104   8.922  1.00  3.58           N
ATOM    193  CA  ILE A  41       9.257  11.549   8.867  1.00  3.55           C
ATOM    194  C   ILE A  41       8.095  12.177   8.150  1.00  3.59           C
ATOM    195  O   ILE A  41       7.399  13.021   8.678  1.00  4.62           O
ATOM    196  CB  ILE A  41      10.614  11.965   8.241  1.00  3.97           C
ATOM    197  CG1 ILE A  41      11.762  11.284   8.954  1.00  4.60           C
ATOM    198  CG2 ILE A  41      10.737  13.460   8.152  1.00  5.47           C
ATOM    199  CD1 ILE A  41      11.735  11.381  10.447  1.00  7.33           C
ATOM    200  H   ILE A  41       9.482   9.706   8.232  1.00  3.58           H
ATOM    201  HA  ILE A  41       9.216  11.890   9.774  1.00  3.55           H
ATOM    202  HB  ILE A  41      10.650  11.651   7.324  1.00  3.97           H
ATOM    203 HG12 ILE A  41      12.591  11.687   8.653  1.00  4.60           H
ATOM    204 HG13 ILE A  41      11.749  10.342   8.725  1.00  4.60           H
ATOM    205 HG21 ILE A  41      10.676  13.836   9.044  1.00  5.47           H
ATOM    206 HG22 ILE A  41      10.018  13.804   7.598  1.00  5.47           H
ATOM    207 HG23 ILE A  41      11.594  13.683   7.757  1.00  5.47           H
ATOM    208 HD11 ILE A  41      10.919  10.971  10.775  1.00  7.33           H
ATOM    209 HD12 ILE A  41      11.762  12.316  10.703  1.00  7.33           H
ATOM    210 HD13 ILE A  41      12.507  10.917  10.808  1.00  7.33           H
ATOM    211  N   CYS A  42       7.817  11.737   6.928  1.00  3.39           N
ATOM    212  CA  CYS A  42       6.813  12.362   6.092  1.00  2.32           C
ATOM    213  C   CYS A  42       5.546  11.599   5.836  1.00  2.80           C
ATOM    214  O   CYS A  42       4.605  12.109   5.233  1.00  3.24           O
ATOM    215  CB  CYS A  42       7.416  12.822   4.772  1.00  3.36           C
ATOM    216  SG  CYS A  42       7.709  11.447   3.585  1.00  3.47           S
ATOM    217  H   CYS A  42       8.206  11.065   6.557  1.00  3.39           H
ATOM    218  HA  CYS A  42       6.515  13.122   6.616  1.00  2.32           H
ATOM    219  HB2 CYS A  42       8.269  13.249   4.949  1.00  3.36           H
ATOM    220  HB3 CYS A  42       6.809  13.453   4.353  1.00  3.36           H
ATOM    222  N   GLY A  43       5.529  10.292   6.224  1.00  2.45           N
ATOM    223  CA  GLY A  43       4.397   9.433   6.041  1.00  3.04           C
ATOM    224  C   GLY A  43       4.241   8.817   4.697  1.00  3.00           C
ATOM    225  O   GLY A  43       3.225   8.146   4.401  1.00  4.16           O
ATOM    226  H   GLY A  43       6.193   9.898   6.603  1.00  2.45           H
ATOM    227  HA2 GLY A  43       4.457   8.709   6.684  1.00  3.04           H
ATOM    228  HA3 GLY A  43       3.593   9.947   6.217  1.00  3.04           H
ATOM    229  N   ALA A  44       5.237   8.944   3.824  1.00  3.45           N
ATOM    230  CA  ALA A  44       5.165   8.344   2.481  1.00  3.40           C
ATOM    231  C   ALA A  44       5.111   6.836   2.552  1.00  3.56           C
ATOM    232  O   ALA A  44       5.790   6.203   3.365  1.00  3.89           O
ATOM    233  CB  ALA A  44       6.378   8.747   1.654  1.00  3.81           C
ATOM    234  H   ALA A  44       5.967   9.371   3.979  1.00  3.45           H
ATOM    235  HA  ALA A  44       4.358   8.669   2.051  1.00  3.40           H
ATOM    236  HB1 ALA A  44       7.182   8.437   2.100  1.00  3.81           H
ATOM    237  HB2 ALA A  44       6.397   9.713   1.572  1.00  3.81           H
ATOM    238  HB3 ALA A  44       6.310   8.342   0.775  1.00  3.81           H
ATOM    239  N   PRO A  45       4.346   6.215   1.643  1.00  3.69           N
ATOM    240  CA  PRO A  45       4.298   4.754   1.585  1.00  4.02           C
ATOM    241  C   PRO A  45       5.576   4.171   1.101  1.00  3.00           C
ATOM    242  O   PRO A  45       6.468   4.880   0.549  1.00  3.12           O
ATOM    243  CB  PRO A  45       3.135   4.457   0.624  1.00  5.15           C
ATOM    244  CG  PRO A  45       2.917   5.686  -0.103  1.00  6.41           C
ATOM    245  CD  PRO A  45       3.406   6.842   0.655  1.00  4.15           C
ATOM    246  HA  PRO A  45       4.079   4.377   2.452  1.00  4.02           H
ATOM    247  HB2 PRO A  45       3.382   3.738   0.022  1.00  5.15           H
ATOM    248  HB3 PRO A  45       2.345   4.210   1.130  1.00  5.15           H
ATOM    249  HG2 PRO A  45       3.390   5.629  -0.948  1.00  6.41           H
ATOM    250  HG3 PRO A  45       1.965   5.782  -0.262  1.00  6.41           H
ATOM    251  HD2 PRO A  45       3.874   7.464   0.076  1.00  4.15           H
ATOM    252  HD3 PRO A  45       2.677   7.294   1.108  1.00  4.15           H
ATOM    253  N   GLU A  48       7.481   4.933  -2.272  1.00  2.79           N
ATOM    254  CA  GLU A  48       8.221   6.151  -2.488  1.00  3.20           C
ATOM    255  C   GLU A  48       9.676   6.047  -2.095  1.00  2.89           C
ATOM    256  O   GLU A  48      10.378   7.082  -2.151  1.00  3.92           O
ATOM    257  CB  GLU A  48       7.533   7.364  -1.863  1.00  3.74           C
ATOM    258  CG  GLU A  48       6.066   7.505  -2.237  1.00  3.93           C
ATOM    259  CD  GLU A  48       5.870   7.509  -3.737  1.00  3.85           C
ATOM    260  OE1 GLU A  48       6.135   8.569  -4.357  1.00  5.88           O
ATOM    261  OE2 GLU A  48       5.533   6.437  -4.316  1.00  4.39           O
ATOM    262  H   GLU A  48       7.138   4.867  -1.486  1.00  2.79           H
ATOM    263  HA  GLU A  48       8.228   6.313  -3.444  1.00  3.20           H
ATOM    264  HB2 GLU A  48       7.991   8.167  -2.157  1.00  3.74           H
ATOM    265  HB3 GLU A  48       7.585   7.287  -0.897  1.00  3.74           H
ATOM    266  HG2 GLU A  48       5.724   8.341  -1.883  1.00  3.93           H
ATOM    267  HG3 GLU A  48       5.568   6.759  -1.867  1.00  3.93           H
ATOM    268  N   PHE A  49      10.156   4.875  -1.748  1.00  2.54           N
ATOM    269  CA  PHE A  49      11.543   4.634  -1.421  1.00  3.16           C
ATOM    270  C   PHE A  49      12.287   4.035  -2.598  1.00  3.23           C
ATOM    271  O   PHE A  49      11.739   3.260  -3.399  1.00  4.57           O
ATOM    272  CB  PHE A  49      11.623   3.646  -0.221  1.00  3.07           C
ATOM    273  CG  PHE A  49      11.333   4.340   1.105  1.00  2.65           C
ATOM    274  CD1 PHE A  49      10.072   4.713   1.442  1.00  2.91           C
ATOM    275  CD2 PHE A  49      12.383   4.683   1.940  1.00  2.86           C
ATOM    276  CE1 PHE A  49       9.830   5.384   2.660  1.00  3.94           C
ATOM    277  CE2 PHE A  49      12.152   5.386   3.132  1.00  3.21           C
ATOM    278  CZ  PHE A  49      10.853   5.721   3.494  1.00  3.61           C
ATOM    279  H   PHE A  49       9.675   4.164  -1.691  1.00  2.54           H
ATOM    280  HA  PHE A  49      11.967   5.474  -1.186  1.00  3.16           H
ATOM    281  HB2 PHE A  49      12.515   3.268  -0.178  1.00  3.07           H
ATOM    282  HB3 PHE A  49      10.969   2.940  -0.344  1.00  3.07           H
ATOM    283  HD1 PHE A  49       9.365   4.526   0.868  1.00  2.91           H
ATOM    284  HD2 PHE A  49      13.251   4.445   1.708  1.00  2.86           H
ATOM    285  HE1 PHE A  49       8.957   5.600   2.897  1.00  3.94           H
ATOM    286  HE2 PHE A  49      12.866   5.627   3.678  1.00  3.21           H
ATOM    287  HZ  PHE A  49      10.685   6.168   4.292  1.00  3.61           H
HETATM  288 FE    FE A  55       9.794  10.608   3.873  1.00  2.90          FE
HETATM  289  O   HOH A 107       5.079   8.927   9.538  1.00 12.57           O
HETATM  290  O   HOH A 123       9.120  16.310   5.902  1.00 10.75           O
HETATM  291  O   HOH A 124      12.913  16.679   4.634  1.00  5.39           O
HETATM  292  O   HOH A 148      13.863  16.664   0.689  1.00 19.89           O
HETATM  293  O   HOH A 158      13.763  14.830  -1.239  1.00 16.61           O
HETATM  294  O   HOH A 164       8.536  15.979   3.044  1.00  8.33           O
HETATM  295  O   HOH A 183      11.838  17.652   2.302  1.00 20.40           O
HETATM  296  O   HOH A 307       9.914  17.434   0.789  0.66 11.76           O
HETATM  297  O   HOH A 422       9.389  18.060   1.780  0.34  9.73           O
END
"""


def _run_endo_exo_on_string(pdb_str, radius=None, include=None, selection=None):
  """Drive ``mmtbx.programs.endo_exo.Program`` in-memory on a PDB string
  with default settings (metal scan, radius=5.0, depth=3).  Parses the
  string with ``iotbx.pdb`` directly -- no disk roundtrip.  Returns
  the single result dict produced for the seed.  *radius* overrides
  ``params.buffer.radius`` when given; *include*, when given, is a
  ``(selection, scope, proximity)`` tuple configuring
  ``residues_to_include``; *selection*, when given, seeds from that CCTBX
  selection string instead of scanning for metals."""
  pdb_in = iotbx.pdb.input(source_info=None, lines=pdb_str.split("\n"))
  model = mmtbx.model.manager(model_input=pdb_in)
  dm = DataManager(["model"])
  dm.add_model("model", model)
  dm.set_default_model("model")

  master = libtbx.phil.parse(EndoexoProgram.master_phil_str)
  params = master.extract()
  params.write_files = False
  if radius is not None:
    params.buffer.radius = radius
  if selection is not None:
    params.selection = [selection]
  if include is not None:
    (params.residues_to_include.selection,
     params.residues_to_include.scope,
     params.residues_to_include.proximity) = include

  prog = EndoexoProgram(dm, params, master_phil=master, logger=io.StringIO())
  prog.validate()
  prog.run()
  results = prog.get_results()
  assert len(results) == 1, (
    f"expected 1 submodel (1 seed atom in input); got {len(results)}")
  return results[0]


def _residue_atom_names(result, resseq):
  """Return the sorted set of atom names for residue *resseq* in the
  materialized submodel of *result* (empty set if the residue is absent)."""
  names = set()
  for rg in result["model"].get_hierarchy().residue_groups():
    if rg.resseq.strip() == resseq:
      for ag in rg.atom_groups():
        for a in ag.atoms():
          names.add(a.name.strip())
  return names


def exercise_residues_to_include():
  """residues_to_include pulls a residue into the region whole, exempt from
  the sidechain cut rules, gated by the per_seed proximity sphere.

  Target: Lys 7 of 1BQ8.  Ordinary growth reaches its backbone but cuts the
  sidechain back at CB, and its nearest atom sits 5.88 A from the Fe, so it
  straddles the proximity threshold -- ideal for exercising both the gate
  and the include path.  What the include path changes is the sidechain
  beyond CB, so that is what is compared."""
  full_lys7 = {"N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ"}
  # What ordinary growth leaves: the backbone and CB, the sidechain cut
  # there.  Compared exactly, so creep past CB fails and so does losing CB.
  cut_back_lys7 = {"N", "CA", "C", "O", "CB"}

  def heavy_lys7(result):
    return {name for name in _residue_atom_names(result, "7")
            if not name.startswith("H")}

  # Baseline: the sidechain is cut back at CB.
  base = _run_endo_exo_on_string(_1BQ8_FE_SPHERE_PDB)
  assert base["model"].get_number_of_atoms() == 170
  assert heavy_lys7(base) == cut_back_lys7, (
    f"expected the sidechain cut at CB, got {sorted(heavy_lys7(base))}")

  # per_seed, proximity below the 5.88 A closest approach -> still excluded.
  near_excl = _run_endo_exo_on_string(
    _1BQ8_FE_SPHERE_PDB, include=("resseq 7", "per_seed", 5.0))
  assert near_excl["model"].get_number_of_atoms() == 170
  assert heavy_lys7(near_excl) == cut_back_lys7

  # per_seed, proximity above 5.88 A -> included whole (all 9 heavy atoms),
  # regardless of the sidechain cut rules.
  near_incl = _run_endo_exo_on_string(
    _1BQ8_FE_SPHERE_PDB, include=("resseq 7", "per_seed", 7.0))
  assert near_incl["model"].get_number_of_atoms() == 178
  assert heavy_lys7(near_incl) == full_lys7

  # global ignores proximity: included even with a tiny sphere.
  glob = _run_endo_exo_on_string(
    _1BQ8_FE_SPHERE_PDB, include=("resseq 7", "global", 1.0))
  assert glob["model"].get_number_of_atoms() == 178
  assert heavy_lys7(glob) == full_lys7

  # Whole-residue expansion: a single-atom selection still pulls the
  # complete residue.
  expand = _run_endo_exo_on_string(
    _1BQ8_FE_SPHERE_PDB, include=("resseq 7 and name NZ", "global", 1.0))
  assert heavy_lys7(expand) == full_lys7


def exercise_submodel_shape():
  """Lock in the structural shape of the submodel: total atom count,
  element distribution, seed and cap iseq counts."""
  result = _run_endo_exo_on_string(_1BQ8_FE_SPHERE_PDB)

  atoms = list(result["model"].get_hierarchy().atoms())
  assert len(atoms) == 170, (
    f"submodel atom count drifted: expected 170, got {len(atoms)}")

  elements = {}
  for a in atoms:
    el = a.element.strip().upper()
    elements[el] = elements.get(el, 0) + 1
  expected_elements = {"C": 56, "FE": 1, "H": 84, "N": 12, "O": 13, "S": 4}
  assert elements == expected_elements, (
    f"submodel element distribution drifted:\n"
    f"  expected: {expected_elements}\n"
    f"  got     : {elements}")

  assert len(result["seed_iseqs"]) == 1
  assert len(result["cap_iseqs"]) == 8
  seed_iseq = result["seed_iseqs"][0]
  assert atoms[seed_iseq].element.strip().upper() == "FE"

  # Each cap's heavy-atom anchor is recorded (caps sit 1.1 A from their
  # anchor): anchors are heavy atoms, and every cap has one within bonding
  # distance.
  anchors = result["cap_anchor_iseqs"]
  assert anchors, "no cap anchors recorded"
  for ai in anchors:
    assert atoms[ai].element.strip().upper() != "H", (
      f"cap anchor {ai} is not a heavy atom")
  for ci in result["cap_iseqs"]:
    cap_xyz = matrix.col(atoms[ci].xyz)
    assert any((matrix.col(atoms[ai].xyz) - cap_xyz).length() < 1.2
               for ai in anchors), (
      f"cap {ci} has no recorded anchor within bonding distance")


def exercise_cys_coordination():
  """Verify the chemistry: the four Cys side chains coordinate the Fe
  with Sg atoms within 3 A of the seed."""
  result = _run_endo_exo_on_string(_1BQ8_FE_SPHERE_PDB)
  hier = result["model"].get_hierarchy()
  atoms = list(hier.atoms())

  cys_residues = [
    ag for ag in hier.atom_groups()
    if ag.resname.strip().upper() == "CYS"]
  assert len(cys_residues) == 4, (
    f"expected 4 Cys residues in submodel; got {len(cys_residues)}")

  seed_iseq = result["seed_iseqs"][0]
  fe_xyz = matrix.col(atoms[seed_iseq].xyz)
  sg_distances = []
  for ag in cys_residues:
    sg = next((a for a in ag.atoms()
               if a.name.strip().upper() == "SG"), None)
    assert sg is not None, (
      f"Cys residue {ag.parent().resseq.strip()} has no SG atom")
    d = (matrix.col(sg.xyz) - fe_xyz).length()
    sg_distances.append(d)
    assert d < 3.0, (
      f"Cys-SG -> Fe distance {d:.2f} A is outside coordination range")

  mean_d = sum(sg_distances) / len(sg_distances)
  assert approx_equal(mean_d, 2.27, eps=0.1)


# Asp 93 / Arg 89 / Arg 92 / Phe 90 + Fe 1210 + three Fe-shell waters
# from 2C2U (HsFt-like ferritin core, P 2 3, a=b=c=90.369).  Fe sits on
# a 3-fold special position: three symmetry operators (1_555, 6_566 and
# 12_665) map ASU Asp 93 / HOH 2154 onto the metal at 2.5-2.6 A, so the
# materialized QM region must contain three Asp 93 copies and three HOH
# 2154 copies plus a single deduplicated Fe.  Fe occupancy is 0.5 in
# 2C2U; bumped to 1.00 here so the BFS metal-detection threshold (which
# uses crystallographic occupancy) is not the thing being exercised.
_2C2U_FE_SPHERE_PDB = """\
CRYST1   90.369   90.369   90.369  90.00  90.00  90.00 P 2 3
ATOM      1  N   ARG A  89      29.832  67.632  19.953  1.00  6.64           N
ATOM      2  CA  ARG A  89      29.359  67.682  21.326  1.00  7.99           C
ATOM      3  C   ARG A  89      29.123  66.317  21.969  1.00  8.16           C
ATOM      4  O   ARG A  89      28.978  66.252  23.181  1.00  8.85           O
ATOM      5  CB  ARG A  89      28.123  68.590  21.438  1.00  8.08           C
ATOM      6  CG  ARG A  89      26.860  67.986  20.809  1.00  9.06           C
ATOM      7  CD  ARG A  89      25.899  69.081  20.404  1.00  9.62           C
ATOM      8  NE  ARG A  89      24.608  68.586  19.936  1.00  9.83           N
ATOM      9  CZ  ARG A  89      24.381  68.082  18.723  1.00  9.19           C
ATOM     10  NH1 ARG A  89      25.372  68.013  17.839  1.00 10.52           N
ATOM     11  NH2 ARG A  89      23.190  67.671  18.374  1.00 11.31           N
ATOM     12  H   ARG A  89      29.235  67.839  19.370  1.00  6.64           H
ATOM     13  HA  ARG A  89      30.070  68.072  21.858  1.00  7.99           H
ATOM     14  HB2 ARG A  89      28.309  69.428  20.986  1.00  8.08           H
ATOM     15  HB3 ARG A  89      27.939  68.753  22.376  1.00  8.08           H
ATOM     16  HG2 ARG A  89      27.101  67.478  20.019  1.00  9.06           H
ATOM     17  HG3 ARG A  89      26.419  67.410  21.453  1.00  9.06           H
ATOM     18  HD2 ARG A  89      26.297  69.596  19.684  1.00  9.62           H
ATOM     19  HD3 ARG A  89      25.736  69.653  21.170  1.00  9.62           H
ATOM     20  HE  ARG A  89      23.945  68.622  20.483  1.00  9.83           H
ATOM     21 HH11 ARG A  89      25.228  67.688  17.056  1.00 10.52           H
ATOM     22 HH12 ARG A  89      26.157  68.293  18.050  1.00 10.52           H
ATOM     23 HH21 ARG A  89      23.060  67.348  17.588  1.00 11.31           H
ATOM     24 HH22 ARG A  89      22.537  67.723  18.931  1.00 11.31           H
ATOM     25  N   PHE A  90      29.143  65.230  21.179  1.00  7.41           N
ATOM     26  CA  PHE A  90      29.067  63.858  21.670  1.00  8.07           C
ATOM     27  C   PHE A  90      30.440  63.196  21.800  1.00  6.75           C
ATOM     28  O   PHE A  90      30.528  62.012  22.041  1.00  6.99           O
ATOM     29  CB  PHE A  90      28.107  63.005  20.854  1.00  9.51           C
ATOM     30  CG  PHE A  90      26.747  63.660  20.727  1.00  9.88           C
ATOM     31  CD1 PHE A  90      25.917  63.757  21.804  1.00 12.35           C
ATOM     32  CD2 PHE A  90      26.356  64.247  19.523  1.00 11.28           C
ATOM     33  CE1 PHE A  90      24.652  64.408  21.670  1.00 13.72           C
ATOM     34  CE2 PHE A  90      25.134  64.895  19.415  1.00 13.59           C
ATOM     35  CZ  PHE A  90      24.329  64.988  20.488  1.00 12.95           C
ATOM     36  H   PHE A  90      29.203  65.270  20.322  1.00  7.41           H
ATOM     37  HA  PHE A  90      28.700  63.898  22.567  1.00  8.07           H
ATOM     38  HB2 PHE A  90      28.469  62.879  19.963  1.00  9.51           H
ATOM     39  HB3 PHE A  90      27.992  62.146  21.290  1.00  9.51           H
ATOM     40  HD1 PHE A  90      26.173  63.400  22.624  1.00 12.35           H
ATOM     41  HD2 PHE A  90      26.920  64.203  18.785  1.00 11.28           H
ATOM     42  HE1 PHE A  90      24.058  64.433  22.385  1.00 13.72           H
ATOM     43  HE2 PHE A  90      24.872  65.264  18.602  1.00 13.59           H
ATOM     44  HZ  PHE A  90      23.532  65.461  20.418  1.00 12.95           H
ATOM     45  N   ARG A  92      32.664  63.114  24.285  1.00  6.21           N
ATOM     46  CA  ARG A  92      32.906  62.301  25.482  1.00  8.24           C
ATOM     47  C   ARG A  92      32.082  61.069  25.480  1.00  7.94           C
ATOM     48  O   ARG A  92      32.550  59.985  25.866  1.00 10.38           O
ATOM     49  CB AARG A  92      32.838  63.133  26.771  0.50  8.24           C
ATOM     50  CG AARG A  92      33.908  64.231  26.908  0.50  8.00           C
ATOM     51  CD AARG A  92      35.345  63.768  27.104  0.50  8.06           C
ATOM     52  NE AARG A  92      35.471  63.103  28.390  0.50  7.84           N
ATOM     53  CZ AARG A  92      36.586  62.621  28.904  0.50  9.02           C
ATOM     54  NH1AARG A  92      37.706  62.705  28.233  0.50 14.66           N
ATOM     55  NH2AARG A  92      36.555  62.048  30.101  0.50 11.64           N
ATOM     56  H  AARG A  92      33.127  63.838  24.254  1.00  6.21           H
ATOM     57  HA AARG A  92      33.821  61.979  25.477  1.00  8.24           H
ATOM     58  HB2AARG A  92      32.941  62.534  27.527  0.50  8.24           H
ATOM     59  HB3AARG A  92      31.972  63.568  26.808  0.50  8.24           H
ATOM     60  HG2AARG A  92      33.681  64.778  27.676  0.50  8.00           H
ATOM     61  HG3AARG A  92      33.892  64.769  26.101  0.50  8.00           H
ATOM     62  HD2AARG A  92      35.942  64.533  27.087  0.50  8.06           H
ATOM     63  HD3AARG A  92      35.587  63.143  26.403  0.50  8.06           H
ATOM     64  HE AARG A  92      34.755  63.016  28.858  0.50  7.84           H
ATOM     65 HH11AARG A  92      38.433  62.391  28.568  0.50 14.66           H
ATOM     66 HH12AARG A  92      37.715  63.075  27.457  0.50 14.66           H
ATOM     67 HH21AARG A  92      37.277  61.731  30.445  0.50 11.64           H
ATOM     68 HH22AARG A  92      35.813  61.994  30.533  0.50 11.64           H
ATOM     69  CB BARG A  92      32.507  63.088  26.731  0.50  8.15           C
ATOM     70  CG BARG A  92      33.428  64.191  27.114  0.50  7.60           C
ATOM     71  CD BARG A  92      34.674  63.668  27.784  0.50  8.27           C
ATOM     72  NE BARG A  92      34.344  63.138  29.101  0.50  7.94           N
ATOM     73  CZ BARG A  92      35.190  62.540  29.954  0.50  8.28           C
ATOM     74  NH1BARG A  92      36.484  62.428  29.633  0.50  7.37           N
ATOM     75  NH2BARG A  92      34.767  62.156  31.175  0.50  6.81           N
ATOM     76  H  BARG A  92      33.127  63.838  24.254  1.00  6.21           H
ATOM     77  HA BARG A  92      33.847  62.068  25.500  1.00  8.24           H
ATOM     78  HB2BARG A  92      32.467  62.472  27.480  0.50  8.15           H
ATOM     79  HB3BARG A  92      31.634  63.481  26.578  0.50  8.15           H
ATOM     80  HG2BARG A  92      32.976  64.786  27.733  0.50  7.60           H
ATOM     81  HG3BARG A  92      33.692  64.679  26.318  0.50  7.60           H
ATOM     82  HD2BARG A  92      35.316  64.388  27.889  0.50  8.27           H
ATOM     83  HD3BARG A  92      35.059  62.956  27.249  0.50  8.27           H
ATOM     84  HE BARG A  92      33.527  63.216  29.357  0.50  7.94           H
ATOM     85 HH11BARG A  92      37.030  62.046  30.176  0.50  7.37           H
ATOM     86 HH12BARG A  92      36.770  62.738  28.884  0.50  7.37           H
ATOM     87 HH21BARG A  92      35.313  61.773  31.718  0.50  6.81           H
ATOM     88 HH22BARG A  92      33.951  62.294  31.411  0.50  6.81           H
ATOM     89  N   ASP A  93      30.815  61.167  25.066  1.00  7.08           N
ATOM     90  CA  ASP A  93      29.930  59.993  25.026  1.00  7.52           C
ATOM     91  C   ASP A  93      30.579  58.852  24.237  1.00  6.59           C
ATOM     92  O   ASP A  93      30.664  57.719  24.694  1.00  8.25           O
ATOM     93  CB  ASP A  93      28.578  60.370  24.419  1.00  9.20           C
ATOM     94  CG  ASP A  93      27.912  61.536  25.111  1.00 11.51           C
ATOM     95  OD1 ASP A  93      27.283  61.299  26.128  1.00 18.61           O
ATOM     96  OD2 ASP A  93      28.068  62.711  24.724  1.00 16.53           O
ATOM     97  H   ASP A  93      30.444  61.897  24.804  1.00  7.08           H
ATOM     98  HA  ASP A  93      29.770  59.675  25.928  1.00  7.52           H
ATOM     99  HB2 ASP A  93      28.709  60.612  23.489  1.00  9.20           H
ATOM    100  HB3 ASP A  93      27.982  59.607  24.483  1.00  9.20           H
HETATM  101 FE   FE  A1210      26.664  63.705  26.664  1.00  8.65          FE
HETATM  102  O   HOH A2153      30.815  65.185  24.984  1.00 16.78           O
HETATM  103  O   HOH A2154      26.856  65.362  24.718  1.00 13.80           O
HETATM  104  O   HOH A2162      28.848  59.920  28.597  1.00 33.42           O
END
"""


# A Zn on a 2-fold special position (P 2 2 2; the 2-fold along a is (x,-y,-z))
# with a coordinating water ON the axis (its O maps onto itself while its two H
# sit in a general orientation, so the fixing op places two more H) plus a
# second water OFF the axis (a genuine two-copy symmetry pair).
_SPECIAL_POS_WATER_PDB = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 2 2 2
HETATM    1 ZN    ZN A 262      10.000   0.000   0.000  1.00 10.00          ZN
HETATM    2  O   HOH A 300      12.200   0.000   0.000  1.00 10.00           O
ATOM      3  H1  HOH A 300      12.780   0.620   0.420  1.00 10.00           H
ATOM      4  H2  HOH A 300      12.780   0.420   0.620  1.00 10.00           H
HETATM    5  O   HOH A 301      11.000   1.500   1.500  1.00 10.00           O
ATOM      6  H1  HOH A 301      11.500   2.100   1.900  1.00 10.00           H
ATOM      7  H2  HOH A 301      11.500   1.900   2.100  1.00 10.00           H
END
"""


def exercise_2c2u_symmetry_materialization():
  """Lock in the symmetry-expanded shape of the 2C2U Fe region.

  Asserts:

  * exactly one Fe in the materialized model (the three sym_ops put the
    metal at the same point; only one survives deduplication);
  * three Asp 93 atom groups, one per coordinating sym_op;
  * three HOH 2154 oxygens, one per coordinating sym_op;
  * chain IDs are all single character (PDB-compatible) and include
    the parent "A" plus at least two additional chains for the
    symmetry images.
  """
  result = _run_endo_exo_on_string(_2C2U_FE_SPHERE_PDB)
  hier = result["model"].get_hierarchy()

  fe_atoms = 0
  asp93_groups = 0
  hoh2154_groups = 0
  chain_ids = set()
  for ch in hier.chains():
    chain_ids.add(ch.id.strip())
    assert len(ch.id.strip()) == 1, (
      f"chain id {ch.id!r} is not a single character; this would "
      f"overflow the PDB chain field and break downstream restraint "
      f"reconstruction")
    for rg in ch.residue_groups():
      for ag in rg.atom_groups():
        rn = ag.resname.strip().upper()
        if rn == "FE":
          fe_atoms += len(list(ag.atoms()))
        elif rn == "ASP" and rg.resseq.strip() == "93":
          asp93_groups += 1
        elif rn == "HOH" and rg.resseq.strip() == "2154":
          hoh2154_groups += 1

  assert fe_atoms == 1, (
    f"expected 1 Fe atom after special-position dedup; got {fe_atoms}")
  assert asp93_groups == 3, (
    f"expected 3 Asp 93 atom groups (one per sym_op); got {asp93_groups}")
  assert hoh2154_groups == 3, (
    f"expected 3 HOH 2154 atom groups (one per sym_op); got "
    f"{hoh2154_groups}")

  assert "A" in chain_ids, (
    f"parent chain 'A' missing from materialized region: {chain_ids}")
  assert len(chain_ids) >= 3, (
    f"expected >=3 chain ids (identity + two sym images); got "
    f"{sorted(chain_ids)}")


def exercise_2c2u_fe_coordination_distances():
  """The three Asp 93 OD1/OD2 and three HOH 2154 oxygens that survive
  materialization should all be within Fe coordination range
  (<= 3.0 A from the deduplicated Fe)."""
  result = _run_endo_exo_on_string(_2C2U_FE_SPHERE_PDB)
  hier = result["model"].get_hierarchy()
  atoms = list(hier.atoms())

  seed_iseq = result["seed_iseqs"][0]
  fe_xyz = matrix.col(atoms[seed_iseq].xyz)
  assert atoms[seed_iseq].element.strip().upper() == "FE"

  asp_o_distances = []
  hoh_o_distances = []
  for a in atoms:
    rn = a.parent().resname.strip().upper()
    rg = a.parent().parent()
    name = a.name.strip().upper()
    if rn == "ASP" and rg.resseq.strip() == "93" and name in ("OD1", "OD2"):
      asp_o_distances.append((matrix.col(a.xyz) - fe_xyz).length())
    elif rn == "HOH" and rg.resseq.strip() == "2154":
      hoh_o_distances.append((matrix.col(a.xyz) - fe_xyz).length())

  # 3 Asp93 * (OD1 + OD2) = 6 carboxylate oxygens; at least 3 of them
  # (one per Asp93 image) coordinate the Fe.
  close_asp_o = [d for d in asp_o_distances if d < 3.0]
  assert len(close_asp_o) >= 3, (
    f"expected >=3 Asp 93 oxygens within 3 A of Fe; got "
    f"{sorted(asp_o_distances)}")

  assert len(hoh_o_distances) == 3
  for d in hoh_o_distances:
    assert d < 3.0, (
      f"HOH 2154 O-Fe distance {d:.2f} A outside coordination range")


def exercise_2c2u_symmetry_truncation_consistency():
  """At radius=6 the three symmetry copies of Asp 93 must be truncated
  identically.

  Asp 93's CA (5.2 A) and CB (4.5 A) both sit inside a 6 A sphere around
  the Fe.  The radius search is symmetry-aware, so every copy's CA/CB are
  seeded and protected -- none is cut at its preferred CA-CB site.  Before
  the symmetry-aware seeding the identity copy kept CA/CB (its atoms were
  ASU seeds) while the symmetry images, reached only by BFS, were cut at
  CA-CB; this pins that asymmetry shut."""
  result = _run_endo_exo_on_string(_2C2U_FE_SPHERE_PDB, radius=6.0)
  hier = result["model"].get_hierarchy()

  asp_copies = []
  for ch in hier.chains():
    for rg in ch.residue_groups():
      for ag in rg.atom_groups():
        if ag.resname.strip().upper() == "ASP" and rg.resseq.strip() == "93":
          # real (non-cap) heavy-atom names; a capped atom carries element H
          real_heavy = {a.name.strip() for a in ag.atoms()
                        if a.element.strip().upper() != "H"}
          asp_copies.append((ch.id.strip(), real_heavy))

  assert len(asp_copies) == 3, (
    f"expected 3 Asp 93 copies (one per 3-fold sym_op); got "
    f"{len(asp_copies)}")
  for chain_id, real_heavy in asp_copies:
    assert {"CA", "CB"} <= real_heavy, (
      f"Asp 93 copy in chain {chain_id} is missing real CA/CB (cut at "
      f"CA-CB instead of keeping the in-radius backbone): "
      f"{sorted(real_heavy)}")
  # All three copies must keep the *same* heavy-atom set.
  atom_sets = {frozenset(real_heavy) for _id, real_heavy in asp_copies}
  assert len(atom_sets) == 1, (
    f"Asp 93 symmetry copies truncated differently: "
    f"{[(cid, sorted(s)) for cid, s in asp_copies]}")


def exercise_residue_composition():
  """The default buffer (radius=5, depth=3) pulls in a stable scaffold
  around the four coordinating Cys.  Pin the residue type counts."""
  result = _run_endo_exo_on_string(_1BQ8_FE_SPHERE_PDB)
  hier = result["model"].get_hierarchy()

  resname_counts = {}
  for ag in hier.atom_groups():
    rn = ag.resname.strip().upper()
    resname_counts[rn] = resname_counts.get(rn, 0) + 1

  expected = {
    "ALA": 1, "CYS": 4, "FE": 1, "GLY": 2, "ILE": 2,
    "LYS": 1, "PHE": 1, "PRO": 1, "TYR": 1, "VAL": 2,
  }
  assert resname_counts == expected, (
    f"residue composition drifted:\n"
    f"  expected: {expected}\n"
    f"  got     : {resname_counts}")


# mmtbx.regression.model_1yjp with hydrogens, since the program requires
# them and the shared fixture has none.  A local copy rather than an edit
# to the shared one, which other tests use.  TYR 7 came back with CD1/CD2
# and CE1/CE2 exchanged; the ring is two-fold symmetric, so only the
# naming differs.
_1YJP_PDB = """\
CRYST1   21.937    4.866   23.477  90.00 107.08  90.00 P 1 21 1
ATOM      1  N   GLY A   1      -9.009   4.612   6.102  1.00 16.77           N
ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      3  C   GLY A   1      -8.015   3.140   4.419  1.00 16.16           C
ATOM      4  O   GLY A   1      -7.523   2.521   5.381  1.00 16.78           O
ATOM      5  H1  GLY A   1      -9.618   5.243   6.252  1.00 16.77           H
ATOM      6  H2  GLY A   1      -9.183   3.904   6.613  1.00 16.77           H
ATOM      7  H3  GLY A   1      -8.200   4.928   6.298  1.00 16.77           H
ATOM      8  HA2 GLY A   1      -8.861   4.970   4.084  1.00 16.57           H
ATOM      9  HA3 GLY A   1      -9.929   3.858   4.426  1.00 16.57           H
ATOM     10  N   ASN A   2      -7.656   2.923   3.155  1.00 15.02           N
ATOM     11  CA  ASN A   2      -6.522   2.038   2.831  1.00 14.10           C
ATOM     12  C   ASN A   2      -5.241   2.537   3.427  1.00 13.13           C
ATOM     13  O   ASN A   2      -4.978   3.742   3.426  1.00 11.91           O
ATOM     14  CB  ASN A   2      -6.346   1.881   1.341  1.00 15.38           C
ATOM     15  CG  ASN A   2      -7.584   1.342   0.692  1.00 14.08           C
ATOM     16  OD1 ASN A   2      -8.025   0.227   1.016  1.00 17.46           O
ATOM     17  ND2 ASN A   2      -8.204   2.155  -0.169  1.00 11.72           N
ATOM     18  H   ASN A   2      -8.045   3.269   2.470  1.00 15.02           H
ATOM     19  HA  ASN A   2      -6.719   1.165   3.206  1.00 14.10           H
ATOM     20  HB2 ASN A   2      -6.148   2.746   0.949  1.00 15.38           H
ATOM     21  HB3 ASN A   2      -5.618   1.264   1.168  1.00 15.38           H
ATOM     22 HD21 ASN A   2      -8.919   1.893  -0.569  1.00 11.72           H
ATOM     23 HD22 ASN A   2      -7.888   2.940  -0.323  1.00 11.72           H
ATOM     24  N   ASN A   3      -4.438   1.590   3.905  1.00 12.26           N
ATOM     25  CA  ASN A   3      -3.193   1.904   4.589  1.00 11.74           C
ATOM     26  C   ASN A   3      -1.955   1.332   3.895  1.00 11.10           C
ATOM     27  O   ASN A   3      -1.872   0.119   3.648  1.00 10.42           O
ATOM     28  CB  ASN A   3      -3.259   1.378   6.042  1.00 12.15           C
ATOM     29  CG  ASN A   3      -2.006   1.739   6.861  1.00 12.82           C
ATOM     30  OD1 ASN A   3      -1.702   2.925   7.072  1.00 15.05           O
ATOM     31  ND2 ASN A   3      -1.271   0.715   7.306  1.00 13.48           N
ATOM     32  H   ASN A   3      -4.596   0.747   3.844  1.00 12.26           H
ATOM     33  HA  ASN A   3      -3.087   2.868   4.586  1.00 11.74           H
ATOM     34  HB2 ASN A   3      -3.339   0.411   6.025  1.00 12.15           H
ATOM     35  HB3 ASN A   3      -4.030   1.767   6.485  1.00 12.15           H
ATOM     36 HD21 ASN A   3      -0.561   0.864   7.768  1.00 13.48           H
ATOM     37 HD22 ASN A   3      -1.508  -0.093   7.130  1.00 13.48           H
ATOM     38  N   GLN A   4      -1.005   2.228   3.598  1.00 10.29           N
ATOM     39  CA  GLN A   4       0.384   1.888   3.199  1.00 10.53           C
ATOM     40  C   GLN A   4       1.435   2.606   4.088  1.00 10.24           C
ATOM     41  O   GLN A   4       1.547   3.843   4.115  1.00  8.86           O
ATOM     42  CB  GLN A   4       0.656   2.148   1.711  1.00  9.80           C
ATOM     43  CG  GLN A   4       1.944   1.458   1.213  1.00 10.25           C
ATOM     44  CD  GLN A   4       2.504   2.044  -0.089  1.00 12.43           C
ATOM     45  OE1 GLN A   4       2.744   3.268  -0.190  1.00 14.62           O
ATOM     46  NE2 GLN A   4       2.750   1.161  -1.091  1.00  9.05           N
ATOM     47  H   GLN A   4      -1.138   3.077   3.619  1.00 10.29           H
ATOM     48  HA  GLN A   4       0.492   0.933   3.329  1.00 10.53           H
ATOM     49  HB2 GLN A   4       0.752   3.103   1.569  1.00  9.80           H
ATOM     50  HB3 GLN A   4      -0.088   1.808   1.190  1.00  9.80           H
ATOM     51  HG2 GLN A   4       2.628   1.548   1.895  1.00 10.25           H
ATOM     52  HG3 GLN A   4       1.753   0.520   1.056  1.00 10.25           H
ATOM     53 HE21 GLN A   4       3.065   1.439  -1.842  1.00  9.05           H
ATOM     54 HE22 GLN A   4       2.591   0.324  -0.975  1.00  9.05           H
ATOM     55  N   GLN A   5       2.154   1.821   4.871  1.00 10.38           N
ATOM     56  CA  GLN A   5       3.270   2.361   5.640  1.00 11.39           C
ATOM     57  C   GLN A   5       4.594   1.768   5.172  1.00 11.52           C
ATOM     58  O   GLN A   5       4.768   0.546   5.054  1.00 12.05           O
ATOM     59  CB  GLN A   5       3.056   2.183   7.147  1.00 11.96           C
ATOM     60  CG  GLN A   5       1.829   2.950   7.647  1.00 10.81           C
ATOM     61  CD  GLN A   5       1.344   2.414   8.954  1.00 13.10           C
ATOM     62  OE1 GLN A   5       0.774   1.325   9.002  1.00 10.65           O
ATOM     63  NE2 GLN A   5       1.549   3.187  10.039  1.00 12.30           N
ATOM     64  H   GLN A   5       2.020   0.978   4.976  1.00 10.38           H
ATOM     65  HA  GLN A   5       3.321   3.318   5.493  1.00 11.39           H
ATOM     66  HB2 GLN A   5       3.835   2.514   7.621  1.00 11.96           H
ATOM     67  HB3 GLN A   5       2.925   1.242   7.341  1.00 11.96           H
ATOM     68  HG2 GLN A   5       2.061   3.884   7.767  1.00 10.81           H
ATOM     69  HG3 GLN A   5       1.113   2.867   6.998  1.00 10.81           H
ATOM     70 HE21 GLN A   5       1.287   2.920  10.813  1.00 12.30           H
ATOM     71 HE22 GLN A   5       1.942   3.948   9.957  1.00 12.30           H
ATOM     72  N   ASN A   6       5.514   2.664   4.856  1.00 11.99           N
ATOM     73  CA  ASN A   6       6.831   2.310   4.318  1.00 12.30           C
ATOM     74  C   ASN A   6       7.854   2.761   5.324  1.00 13.40           C
ATOM     75  O   ASN A   6       8.219   3.943   5.374  1.00 13.92           O
ATOM     76  CB  ASN A   6       7.065   3.016   2.993  1.00 12.13           C
ATOM     77  CG  ASN A   6       5.961   2.735   2.003  1.00 12.77           C
ATOM     78  OD1 ASN A   6       5.798   1.604   1.551  1.00 14.27           O
ATOM     79  ND2 ASN A   6       5.195   3.747   1.679  1.00 10.07           N
ATOM     80  H   ASN A   6       5.403   3.512   4.944  1.00 11.99           H
ATOM     81  HA  ASN A   6       6.906   1.357   4.152  1.00 12.30           H
ATOM     82  HB2 ASN A   6       7.101   3.974   3.143  1.00 12.13           H
ATOM     83  HB3 ASN A   6       7.902   2.708   2.610  1.00 12.13           H
ATOM     84 HD21 ASN A   6       4.552   3.637   1.119  1.00 10.07           H
ATOM     85 HD22 ASN A   6       5.335   4.521   2.027  1.00 10.07           H
ATOM     86  N   TYR A   7       8.292   1.817   6.147  1.00 14.70           N
ATOM     87  CA  TYR A   7       9.159   2.144   7.299  1.00 15.18           C
ATOM     88  C   TYR A   7      10.603   2.331   6.885  1.00 15.91           C
ATOM     89  O   TYR A   7      11.041   1.811   5.855  1.00 15.76           O
ATOM     90  CB  TYR A   7       9.061   1.065   8.369  1.00 15.35           C
ATOM     91  CG  TYR A   7       7.665   0.929   8.902  1.00 14.45           C
ATOM     92  CD1 TYR A   7       7.210   1.756   9.920  1.00 14.80           C
ATOM     93  CD2 TYR A   7       6.771   0.021   8.327  1.00 15.68           C
ATOM     94  CE1 TYR A   7       5.904   1.649  10.416  1.00 14.33           C
ATOM     95  CE2 TYR A   7       5.480  -0.094   8.796  1.00 13.46           C
ATOM     96  CZ  TYR A   7       5.047   0.729   9.831  1.00 15.09           C
ATOM     97  OH  TYR A   7       3.766   0.589  10.291  1.00 14.39           O
ATOM     98  OXT TYR A   7      11.358   2.999   7.612  1.00 17.49           O
ATOM     99  H   TYR A   7       8.108   0.980   6.070  1.00 14.70           H
ATOM    100  HA  TYR A   7       8.853   2.982   7.680  1.00 15.18           H
ATOM    101  HB2 TYR A   7       9.647   1.294   9.107  1.00 15.35           H
ATOM    102  HB3 TYR A   7       9.325   0.213   7.988  1.00 15.35           H
ATOM    103  HD1 TYR A   7       7.783   2.394  10.280  1.00 14.80           H
ATOM    104  HD2 TYR A   7       7.051  -0.513   7.619  1.00 15.68           H
ATOM    105  HE1 TYR A   7       5.618   2.183  11.122  1.00 14.33           H
ATOM    106  HE2 TYR A   7       4.901  -0.719   8.422  1.00 13.46           H
ATOM    107  HH  TYR A   7       3.396   1.341  10.344  1.00 14.39           H
END
"""


def exercise_selection_seed_terminates():
  """A non-metal ``selection`` seed on a symmetric crystal grows a bounded,
  covalent-only region containing the seeded residue with its side chain
  kept whole."""
  result = _run_endo_exo_on_string(_1YJP_PDB, selection="resid 3")

  # Finite region, not an unbounded lattice walk.
  n = result["model"].get_number_of_atoms()
  assert 0 < n < 400, f"resid 3 region atom count looks unbounded: {n}"

  # The seeded ASN 3 is present with its side chain kept whole.
  asn3 = _residue_atom_names(result, "3")
  assert {"CB", "CG", "OD1", "ND2"} <= asn3, (
    f"seeded ASN 3 side chain incomplete: {sorted(asn3)}")


# ===========================================================================
# Engine unit tests
#
# These drive the pure engine pieces directly (no metal scan, no full
# Program pipeline) and cover the subtle logic -- symmetry-op
# canonicalisation, cap placement, bond-cut heuristics, and adjacency
# construction -- that the count-based integration tests above do not.
# ===========================================================================

def _atoms_by_name(pdb_str):
  """Parse *pdb_str* and return ``{atom_name: atom}`` with i_seqs reset to
  positional order."""
  pdb_in = iotbx.pdb.input(source_info=None, lines=pdb_str.split("\n"))
  hier = pdb_in.construct_hierarchy()
  atoms = hier.atoms()
  atoms.reset_i_seq()
  return {a.name.strip(): a for a in atoms}


def exercise_canon_op():
  """``_canon_op`` must give a canonical, hash-stable representative: two
  rt_mx values that compare equal (same xyz) must hash equal and satisfy
  set membership after canonicalisation, regardless of how they were
  built.  This is the contract the symmetry-aware BFS relies on when it
  keys ``visited`` / ``cap_candidates`` on ``(i_seq, rt_mx)`` nodes."""
  op = sgtbx.rt_mx("x,y,z+1")
  # Same operation reached via composition with the identity -- this can
  # leave a different internal denominator / representation.
  op_composed = op.multiply(sgtbx.rt_mx())
  assert op.as_xyz() == op_composed.as_xyz()

  c1 = _canon_op(op)
  c2 = _canon_op(op_composed)

  # canonical form preserves the operation ...
  assert c1.as_xyz() == op.as_xyz()
  # ... is idempotent ...
  assert _canon_op(c1).as_xyz() == c1.as_xyz()
  # ... and is hash-stable: equal ops hash equal and live as one set member.
  assert c1 == c2
  assert hash(c1) == hash(c2)
  assert c2 in {c1}

  # The inverse of a pure translation negates it; canonicalisation must
  # survive that too (the reverse-edge case in build_adjacency).
  inv = _canon_op(op.inverse())
  assert inv.as_xyz() == sgtbx.rt_mx("x,y,z-1").as_xyz()


_CAP_PDB = """\
ATOM      1  C1  XXX A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  C2  XXX A   1       1.500   0.000   0.000  1.00  0.00           C
"""


def exercise_hydrogen_capper():
  """``cap_atom`` retypes the cap to hydrogen and moves it to 1.1 A along
  the anchor->cap direction; ``None`` arguments are a no-op."""
  capper = HydrogenCapper(log=io.StringIO())

  by_name = _atoms_by_name(_CAP_PDB)
  anchor, cap = by_name["C1"], by_name["C2"]
  capper.cap_atom(anchor, cap)

  assert cap.element.strip().upper() == "H"
  d = (matrix.col(cap.xyz) - matrix.col(anchor.xyz)).length()
  assert approx_equal(d, 1.1, eps=1e-6), d
  # cap lies on the original +x direction from the anchor
  assert approx_equal(cap.xyz, (1.1, 0.0, 0.0), eps=1e-6), cap.xyz

  # None on either side is a no-op (must not raise / must not mutate).
  fresh = _atoms_by_name(_CAP_PDB)
  before = fresh["C2"].xyz
  capper.cap_atom(None, fresh["C2"])
  capper.cap_atom(fresh["C1"], None)
  assert fresh["C2"].xyz == before
  assert fresh["C2"].element.strip().upper() == "C"


# Lysine fragment: backbone N-CA-C plus the CD-CE sidechain bond that is a
# PREFERRED_CUTS site for LYS.  Geometry is only nominal; the preferred /
# backbone checks are name + adjacency based.
_LYS_PDB = """\
ATOM      1  N   LYS A   1       0.000   0.000   0.000  1.00  0.00           N
ATOM      2  CA  LYS A   1       1.450   0.000   0.000  1.00  0.00           C
ATOM      3  C   LYS A   1       2.000   1.400   0.000  1.00  0.00           C
ATOM      4  CD  LYS A   1       3.000   0.000   0.000  1.00  0.00           C
ATOM      5  CE  LYS A   1       4.500   0.000   0.000  1.00  0.00           C
"""

# Two sp3 carbons 1.54 A apart for the geometric C-C heuristic branch.
_CC_PDB = """\
ATOM      1  CA  XXX A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  CB  XXX A   1       1.540   0.000   0.000  1.00  0.00           C
"""


def _edge(adj, i, j):
  """Add an undirected bare edge (the sym_op is irrelevant to the cut
  checks, which drop it via _neighbour_iseqs)."""
  adj.setdefault(i, set()).add((j, None))
  adj.setdefault(j, set()).add((i, None))


def exercise_bond_cut_preferred():
  """With ``use_preferred_cuts=True`` the PREFERRED_CUTS table decides:
  LYS CD-CE is a cut site, LYS CA-CD is not."""
  by_name = _atoms_by_name(_LYS_PDB)
  det = BondCutDetector(use_preferred_cuts=True, log=io.StringIO())

  adj = {}
  _edge(adj, by_name["CD"].i_seq, by_name["CE"].i_seq)
  _edge(adj, by_name["CA"].i_seq, by_name["CD"].i_seq)

  assert det.is_cc_single_sp3_bond(
    "LYS", by_name["CD"], by_name["CE"], adj) is True
  # CA is not in LYS's preferred {CD, CE} set.
  assert det.is_cc_single_sp3_bond(
    "LYS", by_name["CA"], by_name["CD"], adj) is False


def exercise_bond_cut_heuristic():
  """For a residue absent from PREFERRED_CUTS the geometric heuristic
  applies: two degree-4 sp3 carbons 1.42-1.68 A apart are a cut site;
  shifting the distance out of range disqualifies the bond."""
  by_name = _atoms_by_name(_CC_PDB)
  ca, cb = by_name["CA"], by_name["CB"]
  det = BondCutDetector(use_preferred_cuts=True, log=io.StringIO())

  # Degree-4 carbons: the real C-C bond plus three dummy neighbours each.
  adj = {}
  _edge(adj, ca.i_seq, cb.i_seq)
  for d in (101, 102, 103):
    _edge(adj, ca.i_seq, d)
  for d in (201, 202, 203):
    _edge(adj, cb.i_seq, d)

  # 'XXX' is not in PREFERRED_CUTS -> heuristic branch; 1.54 A is in range.
  assert det.is_cc_single_sp3_bond("XXX", ca, cb, adj) is True

  # Move CB out to 2.0 A: now outside the 1.42-1.68 A window.
  cb.set_xyz((2.0, 0.0, 0.0))
  assert det.is_cc_single_sp3_bond("XXX", ca, cb, adj) is False


def exercise_proline_ring_is_never_cut():
  """PRO's empty PREFERRED_CUTS entry vetoes its sidechain outright.

  The ring closes back onto the residue's own N, so one cut detaches
  nothing and two strand the atoms between them.

  The region half seeds mid-ring on purpose. Growing from the metal, BFS
  reaches PRO at CA or N and demotion heals whatever it cuts, so the ring
  comes out whole either way and pins nothing; entering at CG cuts both
  ways out of it at once and nothing is left to trigger a repair."""
  by_name = _atoms_by_name(_CC_PDB)
  ca, cb = by_name["CA"], by_name["CB"]
  det = BondCutDetector(use_preferred_cuts=True, log=io.StringIO())

  adj = {}
  _edge(adj, ca.i_seq, cb.i_seq)
  for d in (101, 102, 103):
    _edge(adj, ca.i_seq, d)
  for d in (201, 202, 203):
    _edge(adj, cb.i_seq, d)

  # Same atoms, same graph: cuttable for a residue absent from the table,
  # refused for PRO. An entry of None rather than frozenset() fails here.
  assert det.is_cc_single_sp3_bond("XXX", ca, cb, adj) is True
  assert det.is_cc_single_sp3_bond("PRO", ca, cb, adj) is False

  # The veto is the entry, not the knob: turning preferred cuts off puts
  # PRO back on the geometric path like any untabled residue.
  off = BondCutDetector(use_preferred_cuts=False, log=io.StringIO())
  assert off.is_cc_single_sp3_bond("PRO", ca, cb, adj) is True

  ring = {"N", "CA", "CB", "CG", "CD"}
  for seed in ("resseq 40 and name CG", "resseq 45 and name CG"):
    result = _run_endo_exo_params(
      _1BQ8_FE_SPHERE_PDB, buffer__skip_search=True, selection=[seed])[0]
    hierarchy = result["model"].get_hierarchy()
    caps = set(result["cap_iseqs"])
    present = {
      atom.name.strip()
      for i_seq, atom in enumerate(hierarchy.atoms())
      if atom.parent().resname.strip() == "PRO" and i_seq not in caps
      and atom.name.strip() in ring
    }
    assert present == ring, (
      f"PRO seeded at '{seed}' came out in pieces: {sorted(present)}")


# What each PREFERRED_CUTS row must leave attached to the sidechain, stated
# as the chemistry rather than as a restatement of the table: cutting the
# named bond reduces the residue to the minimal model compound of its
# functional group.
_PREFERRED_CUT_FRAGMENTS = {
  "ALA": {"CB"},                                              # methane
  "ARG": {"CD", "NE", "CZ", "NH1", "NH2"},                    # methylguanidine
  "ASN": {"CB", "CG", "OD1", "ND2"},                          # acetamide
  "ASP": {"CB", "CG", "OD1", "OD2"},                          # acetate
  "CYS": {"CB", "SG"},                                        # methanethiol
  "GLN": {"CG", "CD", "OE1", "NE2"},                          # acetamide
  "GLU": {"CG", "CD", "OE1", "OE2"},                          # acetate
  "HIS": {"CB", "CG", "ND1", "CD2", "CE1", "NE2"},            # methylimidazole
  "ILE": {"CB", "CG1", "CG2", "CD1"},                         # butane
  "LEU": {"CG", "CD1", "CD2"},                                # propane
  "LYS": {"CE", "NZ"},                                        # methylamine
  "MET": {"CG", "SD", "CE"},                                  # dimethyl sulfide
  "PHE": {"CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"},      # toluene
  "SER": {"CB", "OG"},                                        # methanol
  "THR": {"CB", "OG1", "CG2"},                                # propan-2-ol
  "TRP": {"CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3",
          "CZ2", "CZ3", "CH2"},                               # methylindole
  "TYR": {"CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"},  # cresol
  "VAL": {"CB", "CG1", "CG2"},                                # propane
}


_MONOMER_SERVER = []


def _monomer_server():
  """Return a shared monomer library server.

  Constructing one re-reads the library's CIF files, which costs far more
  than every lookup made of it here put together, and it holds no cache of
  its own."""
  from mmtbx.monomer_library import server

  if not _MONOMER_SERVER:
    _MONOMER_SERVER.append(server.server())
  return _MONOMER_SERVER[0]


def _ideal_heavy_bonds(resname):
  """Return ``{atom name: set of bonded atom names}`` for the ideal monomer.

  Read from the monomer library rather than from a fixture, so the check
  does not depend on which residues a test structure happens to contain."""
  comp = _monomer_server().get_comp_comp_id_direct(resname)
  adjacency = {}
  for bond in comp.bond_list:
    name1 = bond.atom_id_1.strip()
    name2 = bond.atom_id_2.strip()
    if name1.startswith("H") or name2.startswith("H"):
      continue
    adjacency.setdefault(name1, set()).add(name2)
    adjacency.setdefault(name2, set()).add(name1)
  return adjacency


def exercise_preferred_cut_rows():
  """Each PREFERRED_CUTS row names a real bond of the ideal monomer, and
  cutting it severs exactly the fragment the row exists to retain.

  Asserting the fragment, not merely that the named atoms are bonded, is
  what makes this a pin: LEU ``{CA, CB}`` is a genuine bond and still
  retains CD1/CD2, but it drags an extra carbon into the region."""
  covered = {r for r, e in PREFERRED_CUTS.items() if e}
  assert covered == set(_PREFERRED_CUT_FRAGMENTS), (
    "PREFERRED_CUTS and the expected fragments have drifted apart: "
    f"{sorted(covered ^ set(_PREFERRED_CUT_FRAGMENTS))}")

  for resname, expected in sorted(_PREFERRED_CUT_FRAGMENTS.items()):
    entry = PREFERRED_CUTS[resname]
    assert len(entry) == 2, f"{resname} entry names {len(entry)} atoms"
    name1, name2 = sorted(entry)
    adjacency = _ideal_heavy_bonds(resname)
    assert name2 in adjacency.get(name1, ()), (
      f"{resname} entry {name1}-{name2} is not a bond of the ideal monomer")

    # Sever the bond, then collect what is no longer reachable from CA.
    adjacency[name1].discard(name2)
    adjacency[name2].discard(name1)
    reachable = {"CA"}
    queue = deque(["CA"])
    while queue:
      current = queue.popleft()
      for neighbour in adjacency.get(current, ()):
        if neighbour not in reachable:
          reachable.add(neighbour)
          queue.append(neighbour)

    fragment = set(adjacency) - reachable
    assert fragment == expected, (
      f"{resname} cut at {name1}-{name2} retains {sorted(fragment)}, "
      f"expected {sorted(expected)}")


def exercise_preferred_cuts_without_a_cut():
  """GLY and PRO are the two rows that name no bond, and they mean
  different things: GLY is absent from the table (geometric rule applies),
  PRO is present but empty (its sidechain is never cut)."""
  assert PREFERRED_CUTS["GLY"] is None
  assert PREFERRED_CUTS["PRO"] == frozenset()


def exercise_cut_refused_when_it_would_strand_the_anchor():
  """A residue entered at a terminal sidechain atom is not cut away from
  itself.

  ALA 44's CB has one heavy neighbour, CA, and the table names precisely
  that bond.  Reached at CB -- which happens when the buffer sphere
  catches a methyl hydrogen whose own carbon lies outside it -- the cut
  would leave CB with nothing holding it in the region, and it would come
  out as a free methane.  The guard reads the covalent graph, so it does
  not depend on where BFS arrived from.

  The composition alone proves nothing: with the guard gone the cut is
  made and then undone by :meth:`_reattach_stranded`, so CA comes back
  either way.  What separates them is that the guard needs no repair, so
  the repairs are counted too."""
  from mmtbx.geometry_restraints.endo_exo.grow import QMRegionGrower

  repairs = {"n": 0}
  original = QMRegionGrower._demote_cap_candidate

  def counting(self, *args, **kwargs):
    repairs["n"] += 1
    return original(self, *args, **kwargs)

  QMRegionGrower._demote_cap_candidate = counting
  try:
    result = _run_endo_exo_params(
      _1BQ8_FE_SPHERE_PDB, buffer__skip_search=True,
      selection=["resseq 44 and name CB"])[0]
  finally:
    QMRegionGrower._demote_cap_candidate = original

  hierarchy = result["model"].get_hierarchy()
  caps = set(result["cap_iseqs"])
  kept = {
    atom.name.strip()
    for i_seq, atom in enumerate(hierarchy.atoms())
    if atom.parent().parent().resseq.strip() == "44"
    and i_seq not in caps and not atom.element_is_hydrogen()
  }
  assert "CB" in kept, "the seeded atom itself left the region"
  assert "CA" in kept, (
    f"ALA 44 CB was cut from its only heavy neighbour and stranded: "
    f"kept {sorted(kept)}")
  assert repairs["n"] == 0, (
    f"the cut was made and then repaired ({repairs['n']} repair(s)); the "
    f"guard should have refused it in the first place")


def exercise_no_atom_is_left_unattached():
  """No heavy atom comes out of a region bonded only to caps.

  GLU 48 at radius 7 is the reachable case: CB sits in the buffer sphere
  while both of its heavy bonds, CA-CB and CB-CG, are cuttable, so the
  cuts leave it as a free methane.  Which of the two bonds should survive
  cannot be settled while BFS runs, since a bond it has yet to traverse
  looks like one it is about to cut, so it is settled afterwards."""
  # Caps are hydrogens by the time a region is materialised, so an atom cut
  # away from everything shows up as a heavy atom with no heavy neighbour.
  # A metal ion or a water oxygen is legitimately alone; a carbon is not.
  for radius in (5.0, 6.0, 7.0, 8.0):
    for result in _run_endo_exo_params(
        _1BQ8_FE_SPHERE_PDB, buffer__radius=radius):
      hierarchy = result["model"].get_hierarchy()
      atoms = list(hierarchy.atoms())
      heavy = [
        (i_seq, atom) for i_seq, atom in enumerate(atoms)
        if not atom.element_is_hydrogen()
      ]
      for i_seq, atom in heavy:
        if atom.element.strip().upper() != "C":
          continue
        assert any(
          j_seq != i_seq and atom.distance(other) <= 1.95
          for j_seq, other in heavy
        ), (
          f"at radius {radius}, {atom.parent().resname.strip()} "
          f"{atom.parent().parent().resseq.strip()}/{atom.name.strip()} "
          f"was cut away from every heavy neighbour")


def exercise_bond_cut_backbone():
  """``is_ca_c_bond`` / ``is_ca_n_bond`` fire only on a genuine, bonded
  CA->C and CA->N pair (direction and adjacency both matter)."""
  by_name = _atoms_by_name(_LYS_PDB)
  det = BondCutDetector(use_preferred_cuts=True, log=io.StringIO())

  adj = {}
  _edge(adj, by_name["CA"].i_seq, by_name["C"].i_seq)
  _edge(adj, by_name["CA"].i_seq, by_name["N"].i_seq)

  assert det.is_ca_c_bond(by_name["CA"], by_name["C"], adj) is True
  assert det.is_ca_n_bond(by_name["CA"], by_name["N"], adj) is True
  # Wrong direction / wrong names.
  assert det.is_ca_c_bond(by_name["C"], by_name["CA"], adj) is False
  assert det.is_ca_n_bond(by_name["CA"], by_name["C"], adj) is False
  # Right names but not adjacent -> not a bond.
  assert det.is_ca_c_bond(by_name["CA"], by_name["C"], {}) is False


# Real 1RYO chain-A LYS 296 (heavy atoms + the CD/CG hydrogens).  CD and
# NZ are 4.8 / 4.3 A from the Fe in the full structure; CE is 5.0 A out.
# The CD/CG hydrogens give those carbons covalent degree 4 so the
# geometric C-C heuristic (used by the fallback) accepts the CD-CG bond;
# CD-CG is 1.52 A apart, inside the 1.42-1.68 A sp3 window.
_LYS_296_PDB = """\
CRYST1   44.092   57.252  135.988  90.00  90.00  90.00 P 21 21 21
ATOM      1  N   LYS A 296      46.335  42.452  48.863  1.00  4.54           N
ATOM      2  CA  LYS A 296      47.390  42.927  47.971  1.00  4.20           C
ATOM      3  C   LYS A 296      48.088  41.772  47.253  1.00  4.42           C
ATOM      4  O   LYS A 296      48.337  40.715  47.847  1.00  4.42           O
ATOM      5  CB  LYS A 296      48.426  43.737  48.755  1.00  5.02           C
ATOM      6  CG  LYS A 296      47.999  45.167  49.009  1.00  4.18           C
ATOM      7  CD  LYS A 296      48.823  45.848  50.096  1.00  4.88           C
ATOM      8  CE  LYS A 296      50.305  45.946  49.760  1.00  5.79           C
ATOM      9  NZ  LYS A 296      50.995  46.848  50.737  1.00  7.00           N
ATOM     10  HD2 LYS A 296      48.728  45.278  51.020  1.00  4.88           H
ATOM     11  HD3 LYS A 296      48.447  46.860  50.242  1.00  4.88           H
ATOM     12  HG2 LYS A 296      48.116  45.741  48.090  1.00  4.18           H
ATOM     13  HG3 LYS A 296      46.955  45.175  49.324  1.00  4.18           H
"""


def exercise_consumed_preferred_cut():
  """Reproduce the 1RYO Fe / LYS A 296 case: the radius search seeds CD
  and NZ (4.8 / 4.3 A from the Fe) while CE sits just outside (5.0 A), so
  the LYS preferred cut CD-CE has both endpoints inside the region and can
  no longer be made.

  The geometric heuristic then takes the next sp3 C-C bond, CD-CG, leaving
  the {CD, CE, NZ} tip capped at CG rather than walking inward through the
  sidechain and dragging the backbone in."""
  pdb_in = iotbx.pdb.input(
    source_info=None, lines=_LYS_296_PDB.split("\n"))
  model = mmtbx.model.manager(model_input=pdb_in)
  atoms = model.get_hierarchy().atoms()
  atoms.reset_i_seq()
  iseq = {a.name.strip(): a.i_seq for a in atoms}
  name_of = {a.i_seq: a.name.strip() for a in atoms}
  elem_of = {a.i_seq: a.element.strip().upper() for a in atoms}

  # Tagged adjacency must carry a real rt_mx on each edge (the BFS
  # composes ops); the identity stands in for an intra-ASU bond.
  identity = _canon_op(sgtbx.rt_mx())
  adjacency = defaultdict(set)
  def _bond(a, b):
    adjacency[iseq[a]].add((iseq[b], identity))
    adjacency[iseq[b]].add((iseq[a], identity))
  for a, b in [("N", "CA"), ("CA", "C"), ("C", "O"), ("CA", "CB"),
               ("CB", "CG"), ("CG", "CD"), ("CD", "CE"), ("CE", "NZ"),
               ("CD", "HD2"), ("CD", "HD3"), ("CG", "HG2"), ("CG", "HG3")]:
    _bond(a, b)

  # CD and NZ are the atoms the radius search would seed.
  seeds = {iseq["CD"], iseq["NZ"]}

  det = BondCutDetector(use_preferred_cuts=True, log=io.StringIO())
  grower = QMRegionGrower(det, log=io.StringIO())

  def _split(visited, caps):
    cap_iseqs = {c[0] for c in caps}
    interior = {name_of[i] for (i, _op) in visited
                if i not in cap_iseqs and elem_of[i] not in ("H", "D")}
    return interior, {name_of[i] for i in cap_iseqs}

  for label, (visited, caps) in (
      ("grow_by_depth", grower.grow_by_depth(seeds, adjacency, model)),
      ("grow_region", grower.grow_region(seeds, adjacency, model))):
    interior, cap_names = _split(visited, caps)
    assert interior == {"CD", "CE", "NZ"}, (
      f"{label}: expected the {{CD, CE, NZ}} tip, got {sorted(interior)}")
    assert cap_names == {"CG"}, (
      f"{label}: expected CG as the cap, got {sorted(cap_names)}")
    assert not ({"N", "CA", "CB"} & interior), (
      f"{label}: the backbone must not be pulled in; "
      f"interior={sorted(interior)}")


# Real 1RYO chain-A LEU 62 - ASP 63 dipeptide (peptide bond LEU62 C - ASP63
# N).  In the full structure ASP 63 coordinates the Fe (it is radius-seeded)
# while LEU 62 is >6 A away -- it only enters the region as backbone
# overgrowth off ASP 63.
_LEU62_ASP63_PDB = """\
CRYST1   44.092   57.252  135.988  90.00  90.00  90.00 P 21 21 21
ATOM    867  N   LEU A  62      42.072  46.636  54.794  1.00  4.17           N
ATOM    868  CA  LEU A  62      42.362  48.055  54.653  1.00  4.15           C
ATOM    869  C   LEU A  62      43.824  48.433  54.764  1.00  4.47           C
ATOM    870  O   LEU A  62      44.579  47.852  55.543  1.00  4.55           O
ATOM    871  CB  LEU A  62      41.613  48.852  55.726  1.00  5.28           C
ATOM    872  CG  LEU A  62      40.088  48.813  55.774  1.00  5.85           C
ATOM    873  CD1 LEU A  62      39.590  49.620  56.961  1.00  6.23           C
ATOM    874  CD2 LEU A  62      39.523  49.367  54.480  1.00  5.86           C
ATOM    875  H   LEU A  62      42.395  46.287  55.502  1.00  4.17           H
ATOM    876  HA  LEU A  62      42.047  48.357  53.787  1.00  4.15           H
ATOM    877  HB2 LEU A  62      41.924  48.540  56.590  1.00  5.28           H
ATOM    878  HB3 LEU A  62      41.862  49.783  55.624  1.00  5.28           H
ATOM    886  N   ASP A  63      44.202  49.436  53.978  1.00  4.55           N
ATOM    887  CA  ASP A  63      45.538  50.008  54.033  1.00  4.24           C
ATOM    888  C   ASP A  63      45.620  50.561  55.466  1.00  4.40           C
ATOM    889  O   ASP A  63      44.602  50.957  56.042  1.00  4.34           O
ATOM    890  CB  ASP A  63      45.644  51.143  53.015  1.00  4.42           C
ATOM    891  CG  ASP A  63      46.806  52.071  53.293  1.00  4.13           C
ATOM    892  OD1 ASP A  63      47.954  51.770  52.885  1.00  4.01           O
ATOM    893  OD2 ASP A  63      46.556  53.112  53.936  1.00  4.67           O
ATOM    895  HA  ASP A  63      46.214  49.346  53.894  1.00  4.24           H
ATOM    896  HB2 ASP A  63      45.767  50.764  52.131  1.00  4.42           H
ATOM    897  HB3 ASP A  63      44.829  51.668  53.041  1.00  4.42           H
"""


def exercise_cut_depends_on_where_bfs_arrives():
  """Where BFS meets a residue decides where it is cut.

  ASP 63 is seeded through its sidechain, so BFS reaches its configured
  CA-CB bond from the tip and the carboxylate survives.  LEU 62 enters
  only via the peptide bond; its configured CB-CG bond lies further along
  the sidechain than BFS would otherwise go, so the geometric heuristic
  takes CA-CB first and CG, CD1 and CD2 are dropped."""
  pdb_in = iotbx.pdb.input(
    source_info=None, lines=_LEU62_ASP63_PDB.split("\n"))
  model = mmtbx.model.manager(model_input=pdb_in)
  atoms = model.get_hierarchy().atoms()
  atoms.reset_i_seq()

  by_res = {}
  for a in atoms:
    resseq = a.parent().parent().resseq.strip()
    by_res.setdefault(resseq, {})[a.name.strip()] = a.i_seq
  leu, asp = by_res["62"], by_res["63"]
  name_of = {a.i_seq: a.name.strip() for a in atoms}
  resseq_of = {a.i_seq: a.parent().parent().resseq.strip() for a in atoms}
  elem_of = {a.i_seq: a.element.strip().upper() for a in atoms}

  identity = _canon_op(sgtbx.rt_mx())
  adjacency = defaultdict(set)
  def _bond(i, j):
    adjacency[i].add((j, identity))
    adjacency[j].add((i, identity))
  for a, b in [("N", "CA"), ("CA", "C"), ("C", "O"), ("CA", "CB"),
               ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"), ("N", "H"),
               ("CA", "HA"), ("CB", "HB2"), ("CB", "HB3")]:
    _bond(leu[a], leu[b])
  for a, b in [("N", "CA"), ("CA", "C"), ("C", "O"), ("CA", "CB"),
               ("CB", "CG"), ("CG", "OD1"), ("CG", "OD2"),
               ("CA", "HA"), ("CB", "HB2"), ("CB", "HB3")]:
    _bond(asp[a], asp[b])
  _bond(leu["C"], asp["N"])  # peptide bond

  # ASP 63 is the radius-seeded ligand; LEU 62 carries no seed atom.
  seeds = set(asp.values())

  det = BondCutDetector(use_preferred_cuts=True, log=io.StringIO())
  grower = QMRegionGrower(det, log=io.StringIO())

  def _leu_atoms(visited, caps):
    cap_iseqs = {c[0] for c in caps}
    present = {name_of[i] for (i, _o) in visited if resseq_of[i] == "62"}
    caps_leu = {name_of[i] for i in cap_iseqs if resseq_of[i] == "62"}
    interior = {n for n in present
                if n not in caps_leu and elem_of[leu[n]] not in ("H", "D")}
    return interior, caps_leu, present

  # LEU is a spectator, so it is cut at CA-CB and its sidechain beyond CB
  # is dropped.  Both entry points agree.
  for label, (visited, caps) in (
      ("grow_by_depth", grower.grow_by_depth(seeds, adjacency, model)),
      ("grow_region", grower.grow_region(seeds, adjacency, model))):
    _interior, caps_leu, present = _leu_atoms(visited, caps)
    assert "CB" in caps_leu, (
      f"{label}: spectator LEU should cap CB (CA-CB cut); "
      f"caps={sorted(caps_leu)}")
    assert not ({"CG", "CD1", "CD2"} & present), (
      f"{label}: spectator LEU sidechain beyond CB should be dropped; "
      f"present={sorted(present)}")

  # ASP coordinates through its sidechain, so it keeps its preferred cut
  # rather than being trimmed at CA-CB.
  visited, caps = grower.grow_by_depth(seeds, adjacency, model)
  cap_iseqs = {c[0] for c in caps}
  asp_present = {name_of[i] for (i, _o) in visited if resseq_of[i] == "63"}
  asp_caps = {name_of[i] for i in cap_iseqs if resseq_of[i] == "63"}
  assert "CB" not in asp_caps, (
    f"coordinating ASP must not be cut at CA-CB; caps={sorted(asp_caps)}")
  assert {"CB", "CG", "OD1", "OD2"} <= asp_present, (
    f"coordinating ASP should keep its sidechain; "
    f"present={sorted(asp_present)}")

  # Backbone-only seeding is not coordination: seed ASP's mainchain alone
  # and its sidechain is trimmed like any other spectator.
  backbone_seeds = {asp[n] for n in ("N", "CA", "C", "O")}
  visited, caps = grower.grow_by_depth(backbone_seeds, adjacency, model)
  cap_iseqs = {c[0] for c in caps}
  asp_caps = {name_of[i] for i in cap_iseqs if resseq_of[i] == "63"}
  asp_present = {name_of[i] for (i, _o) in visited if resseq_of[i] == "63"}
  assert "CB" in asp_caps, (
    f"backbone-seeded ASP should be cut at CA-CB; caps={sorted(asp_caps)}")
  assert not ({"CG", "OD1", "OD2"} & asp_present), (
    f"backbone-seeded ASP sidechain should be dropped; "
    f"present={sorted(asp_present)}")


def _simple_bond_proxies(pdb_str):
  pdb_in = iotbx.pdb.input(source_info=None, lines=pdb_str.split("\n"))
  model = mmtbx.model.manager(model_input=pdb_in)
  model.process(
    pdb_interpretation_params=model.get_current_pdb_interpretation_params(),
    make_restraints=True)
  grm = model.get_restraints_manager().geometry
  simple, _asu = grm.get_all_bond_proxies(sites_cart=model.get_sites_cart())
  return simple


def _asu_bond_proxies(pdb_str, distance_cutoff=3.2):
  pdb_in = iotbx.pdb.input(source_info=None, lines=pdb_str.split("\n"))
  model = mmtbx.model.manager(model_input=pdb_in)
  pair_asu_table = model.get_xray_structure().pair_asu_table(
    distance_cutoff=distance_cutoff)
  sorted_asu_proxies = geometry_restraints.bond_sorted_asu_proxies(
    pair_asu_table=pair_asu_table)
  return list(sorted_asu_proxies.asu), pair_asu_table.asu_mappings()


# Four zinc ions at the corners of a regular tetrahedron, with waters placed
# to pin the hull rather than merely a bounding sphere.  Neither test
# structure puts a water inside a region's hull at any radius, so this needs a
# fixture built for it.
#
# Distances from the centroid (30, 30, 30), whose insphere is 2.309 A:
#   HOH 5  2.00 A, 0.31 A INSIDE  a facet, and carries hydrogens
#   HOH 6  3.46 A, 1.16 A inside  (out towards a vertex)
#   HOH 7  3.00 A, 0.69 A OUTSIDE a facet
#   HOH 8  far outside
# HOH 6 sits FURTHER from the centre than HOH 7 while being inside, so no
# sphere about the centroid can take both waters and reject HOH 7; the plane
# test has to be a real one.  The 0.31/0.69 A margins leave no room to shift
# the plane tolerance either.
_HULL_WATER_PDB = """\
CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1
HETATM    1 ZN    ZN A   1      34.000  34.000  34.000  1.00 10.00          ZN
HETATM    2 ZN    ZN A   2      34.000  26.000  26.000  1.00 10.00          ZN
HETATM    3 ZN    ZN A   3      26.000  34.000  26.000  1.00 10.00          ZN
HETATM    4 ZN    ZN A   4      26.000  26.000  34.000  1.00 10.00          ZN
HETATM    5  O   HOH A   5      28.845  28.845  28.845  1.00 10.00           O
HETATM    6  H1  HOH A   5      29.805  28.845  28.845  1.00 10.00           H
HETATM    7  H2  HOH A   5      28.605  29.775  28.845  1.00 10.00           H
HETATM    8  O   HOH A   6      32.000  32.000  32.000  1.00 10.00           O
HETATM    9  O   HOH A   7      28.268  28.268  28.268  1.00 10.00           O
HETATM   10  O   HOH A   8      52.000  52.000  52.000  1.00 10.00           O
END
"""

# A trigonal cell, where the rotations -y,x-y,z and -x+y,-x,z are NOT
# symmetric matrices.  Both other fixtures are P 1 and P 1 2 1, whose
# rotations are diagonal, so R and its transpose are identical there and a
# search that used the wrong one would still find every water.  The hull sits
# at fractional (0.30, 0.20, 0.25) and the water at (0.90, 0.70, 0.25), which
# -y,x-y,z carries onto the hull centre.
_HULL_WATER_TRIGONAL_PDB = """\
CRYST1   40.000   40.000   40.000  90.00  90.00 120.00 P 3
HETATM    1 ZN    ZN A   1       9.000   8.660  12.000  1.00 10.00          ZN
HETATM    2 ZN    ZN A   2      11.000   5.196   8.000  1.00 10.00          ZN
HETATM    3 ZN    ZN A   3       5.000   8.660   8.000  1.00 10.00          ZN
HETATM    4 ZN    ZN A   4       7.000   5.196  12.000  1.00 10.00          ZN
HETATM    5  O   HOH A   5      22.000  24.249  10.000  1.00 10.00           O
END
"""

# The same idea under P 1 2 1, where -x,y,-z carries the hull centre
# (20, 30, 20) onto (40, 30, 40).  The water sits there, so its own position
# is nowhere near the hull and only its symmetry image falls inside.
_HULL_WATER_SYMMETRY_PDB = """\
CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1 2 1
HETATM    1 ZN    ZN A   1      24.000  34.000  24.000  1.00 10.00          ZN
HETATM    2 ZN    ZN A   2      24.000  26.000  16.000  1.00 10.00          ZN
HETATM    3 ZN    ZN A   3      16.000  34.000  16.000  1.00 10.00          ZN
HETATM    4 ZN    ZN A   4      16.000  26.000  24.000  1.00 10.00          ZN
HETATM    5  O   HOH A   5      40.000  30.000  40.000  1.00 10.00           O
HETATM    6  O   HOH A   6      30.000  30.000  30.000  1.00 10.00           O
END
"""


def _hull_water_nodes(pdb_str, enabled=True):
  """Run ``_add_hull_waters`` over a hull made of the fixture's zinc ions.

  The method reads only the crystal symmetry, the sites and the hierarchy,
  so the fixture needs no restraints.  The builder is constructed normally
  rather than through ``__new__``: the method reaches for state with
  ``getattr(self, ..., default)``, so a half-built object would answer
  differently from the one the program uses and the difference would be
  silent."""
  from mmtbx.geometry_restraints.endo_exo.builder import QMRegionBuilder

  pdb_in = iotbx.pdb.input(source_info=None, lines=pdb_str.split("\n"))
  model = mmtbx.model.manager(model_input=pdb_in)
  atoms = list(model.get_hierarchy().atoms())
  identity = _canon_op(sgtbx.rt_mx())

  master = libtbx.phil.parse(EndoexoProgram.master_phil_str)
  params = master.extract()
  params.include_waters_in_convex_hull = enabled
  builder = QMRegionBuilder(params, logger=io.StringIO())

  visited = {(i_seq, identity) for i_seq, atom in enumerate(atoms)
             if atom.element.strip().upper() == "ZN"}
  assert len(visited) == 4, f"fixture changed: {len(visited)} hull vertices"
  return model, atoms, visited, builder._add_hull_waters(model, visited)


def exercise_hbond_partner_symmetry_composition():
  """A partner is placed under the operation that reaches it, composed with
  the one its near atom already carries.

  Every hydrogen bond on the test structures is intra-ASU, so the symmetry
  half of this is unreachable through the pipeline and several distinct
  ways of getting it wrong pass unnoticed.  Driving the method with records
  carrying a three-fold makes them visible: that operation is not its own
  inverse, so composing the wrong way round, or reading the same side of
  the bond in both branches, lands the atom on a different image.

  ``hbond.find`` applies its operation to the acceptor when the pair's
  first atom is the hydrogen and to the donor otherwise, which is why the
  method has to consult ``atom_H.index``."""
  from libtbx import group_args
  from mmtbx.geometry_restraints.endo_exo.builder import QMRegionBuilder

  pdb_in = iotbx.pdb.input(
    source_info=None, lines=_HULL_WATER_TRIGONAL_PDB.split("\n"))
  model = mmtbx.model.manager(model_input=pdb_in)
  master = libtbx.phil.parse(EndoexoProgram.master_phil_str)
  params = master.extract()
  params.include_hbond_partners = True
  builder = QMRegionBuilder(params, logger=io.StringIO())

  turn = sgtbx.rt_mx("-y,x-y,z")
  seed_op = _canon_op(sgtbx.rt_mx("y,x,-z"))
  assert turn.as_xyz() != turn.inverse().as_xyz(), (
    "the operation must not be its own inverse or this proves nothing")
  # The two must not commute either, or composing them the wrong way round
  # gives the same answer.  Two rotations about one axis always commute, so
  # the pair here is a three-fold and a two-fold across it.
  assert (seed_op.multiply(turn).as_xyz()
          != turn.multiply(seed_op).as_xyz()), (
    "the operations commute, so the order of composition is not tested")
  # The hydrogen is one of the pair; the donor recorded beside it is the
  # heavy atom carrying it, a third atom.  Conflating the two makes both
  # records take the same branch and the exercise proves nothing.
  donor, hydrogen, acceptor = 0, 1, 4

  def record(hydrogen_is_first):
    """One hbond.find record, naming which of the pair is the hydrogen."""
    first, second = ((hydrogen, acceptor) if hydrogen_is_first
                     else (acceptor, hydrogen))
    return group_args(
      i=first, j=second, symop=turn,
      atom_D=group_args(index=donor),
      atom_A=group_args(index=acceptor),
      atom_H=group_args(index=hydrogen))

  def placed(hydrogen_is_first, seeded):
    builder._hbond_cache = [record(hydrogen_is_first)]
    found = builder._hbond_partner_nodes(model, {(seeded, seed_op)})
    return {op.as_xyz() for _i_seq, op in found}

  forward = _canon_op(seed_op.multiply(turn)).as_xyz()
  backward = _canon_op(seed_op.multiply(turn.inverse())).as_xyz()
  assert forward != backward

  # With the hydrogen first the operation carries the acceptor, so seeding
  # the donor is the forward direction; with it second the roles swap, and
  # reading one side in both branches gets exactly this the wrong way round.
  for hydrogen_is_first, near, far in ((True, donor, acceptor),
                                       (False, acceptor, donor)):
    assert placed(hydrogen_is_first, near) == {forward}, (
      f"hydrogen first={hydrogen_is_first}: seeding {near} placed the "
      f"partner at {placed(hydrogen_is_first, near)}, expected {forward}")
    assert placed(hydrogen_is_first, far) == {backward}, (
      f"hydrogen first={hydrogen_is_first}: seeding {far} placed the "
      f"partner at {placed(hydrogen_is_first, far)}, expected {backward}")


def exercise_hull_waters():
  """Waters inside the region's convex hull are pulled in whole, waters
  outside are not, and the knob turns the whole thing off.

  Worth pinning because it is reached on real structures -- 8FUM adds 23
  water atoms this way around Fe at radius 8 -- while no test structure
  here triggers it at any radius.

  The fixture is shaped so the shape of the test matters: HOH 6 is further
  from the centroid than HOH 7 yet inside, so no sphere about the centre
  can take the one and reject the other."""
  _model, atoms, visited, out = _hull_water_nodes(_HULL_WATER_PDB)
  added = {atoms[i_seq].parent().parent().resseq.strip()
           for i_seq, _op in out - visited}
  assert added == {"5", "6"}, (
    f"expected the two waters inside the hull, got {sorted(added)}")

  # A residue group enters whole, not just the atom that was tested.
  taken = {(atoms[i_seq].parent().parent().resseq.strip(),
            atoms[i_seq].name.strip()) for i_seq, _op in out - visited}
  assert taken >= {("5", "O"), ("5", "H1"), ("5", "H2")}, (
    f"HOH 5 came in without its hydrogens: {sorted(taken)}")

  _model_off, _atoms_off, visited_off, out_off = _hull_water_nodes(
    _HULL_WATER_PDB, enabled=False)
  assert out_off == visited_off, (
    "include_waters_in_convex_hull=False still added waters")


def exercise_hull_waters_use_symmetry_images():
  """A water is taken under the operation that brings it into the hull.

  Its own position is far outside, so a search reading only the deposited
  coordinates would find nothing here.  The operation recorded has to be
  the space-group one composed with the lattice shift carrying the water to
  its nearest image, since that pair is what places the atom when the
  region is materialised."""
  model, atoms, visited, out = _hull_water_nodes(_HULL_WATER_SYMMETRY_PDB)
  assert (model.crystal_symmetry().space_group().type().lookup_symbol()
          == "P 1 2 1")
  new = out - visited
  assert len(new) == 1, f"expected one water image, got {len(new)}"

  (i_seq, op), = new
  assert atoms[i_seq].parent().resname.strip() == "HOH"
  assert op.as_xyz() == "-x+1,y,-z+1", (
    f"expected the 2-fold composed with its lattice shift, got {op.as_xyz()}")

  # A rotation that is not its own transpose. Both fixtures above are
  # diagonal, so either convention finds their waters and the distinction
  # is invisible; here it is not.
  model, atoms, visited, out = _hull_water_nodes(_HULL_WATER_TRIGONAL_PDB)
  assert (model.crystal_symmetry().space_group().type().lookup_symbol()
          == "P 3")
  new = out - visited
  assert len(new) == 1, f"expected one water image, got {len(new)}"
  (i_seq, op), = new
  assert atoms[i_seq].parent().resname.strip() == "HOH"
  assert op.as_xyz() == "-y+1,x-y,z", (
    f"expected the 3-fold composed with its lattice shift, got {op.as_xyz()}")


def exercise_hull_spans_symmetry_images_already_in_the_region():
  """The hull is built from where the region's atoms actually sit, so a
  symmetry image already in it enlarges the bounding volume.

  With only the deposited copy of the zincs the hull is one small
  tetrahedron and HOH 6 lies 9.2 A outside it.  Adding the image the region
  also holds stretches the hull across both, and the water falls inside.
  Building the hull from deposited coordinates instead would find
  nothing."""
  from mmtbx.geometry_restraints.endo_exo.builder import QMRegionBuilder

  pdb_in = iotbx.pdb.input(
    source_info=None, lines=_HULL_WATER_SYMMETRY_PDB.split("\n"))
  model = mmtbx.model.manager(model_input=pdb_in)
  atoms = list(model.get_hierarchy().atoms())
  identity = _canon_op(sgtbx.rt_mx())
  image = _canon_op(sgtbx.rt_mx("-x+1,y,-z+1"))

  master = libtbx.phil.parse(EndoexoProgram.master_phil_str)
  params = master.extract()
  params.include_waters_in_convex_hull = True
  builder = QMRegionBuilder(params, logger=io.StringIO())

  zinc = [i_seq for i_seq, atom in enumerate(atoms)
          if atom.element.strip().upper() == "ZN"]
  assert len(zinc) == 4

  def waters_found(ops):
    visited = {(i_seq, op) for i_seq in zinc for op in ops}
    out = builder._add_hull_waters(model, visited)
    return {atoms[i_seq].parent().parent().resseq.strip()
            for i_seq, _op in out - visited}

  assert "6" not in waters_found([identity]), (
    "HOH 6 is outside the hull of the deposited zincs alone")
  assert "6" in waters_found([identity, image]), (
    "HOH 6 lies inside the hull spanning both copies of the zincs, so the "
    "hull is not being built from the region's actual positions")


def exercise_hbond_partners():
  """``include_hbond_partners`` brings in the far side of a hydrogen bond
  the cuts ran through, and is off unless asked for.

  The bond here is LYS 7 N-H donating to GLU 48's carbonyl oxygen (2.03 A,
  170 degrees).  At radius 7 the region holds the lysine but caps GLU 48's
  carbonyl carbon away, so the acceptor is missing.  Seeding it brings the
  oxygen in, and growth carries on over the peptide bond into PHE 49, which
  is why atoms of a second residue appear as well."""
  master = libtbx.phil.parse(EndoexoProgram.master_phil_str)
  assert master.extract().include_hbond_partners is False, (
    "include_hbond_partners must stay off by default: it grows a region by "
    "over a tenth on a real structure")

  # Nothing here can be decided without hydrogens, so a model that has
  # none is refused rather than truncated on wrong evidence.
  stripped = "\n".join(line for line in _1BQ8_FE_SPHERE_PDB.split("\n")
                       if line[76:78].strip() != "H")
  assert stripped != _1BQ8_FE_SPHERE_PDB, "the fixture carries no hydrogens"
  try:
    _run_endo_exo_params(stripped, buffer__radius=7.0)
  except Sorry as e:
    assert "hydrogens" in str(e), f"unexpected refusal: {e}"
  else:
    raise AssertionError("a model with no hydrogens was accepted")

  def composition(enabled):
    kept, caps = set(), set()
    for result in _run_endo_exo_params(
        _1BQ8_FE_SPHERE_PDB, buffer__radius=7.0,
        include_hbond_partners=enabled):
      hierarchy = result["model"].get_hierarchy()
      cap_iseqs = set(result["cap_iseqs"])
      for i_seq, atom in enumerate(hierarchy.atoms()):
        if atom.element_is_hydrogen() and i_seq not in cap_iseqs:
          continue
        key = (atom.parent().parent().resseq.strip(), atom.name.strip())
        (caps if i_seq in cap_iseqs else kept).add(key)
    return kept, caps

  kept_off, caps_off = composition(False)
  kept_on, caps_on = composition(True)

  assert kept_on - kept_off, (
    "the knob added nothing; this fixture and radius no longer cut through "
    "a hydrogen bond")
  # Atoms move between interior and cap as the cuts land further out, so
  # compare presence rather than status: nothing may leave the region.
  assert not (kept_off | caps_off) - (kept_on | caps_on), (
    f"atoms left the region: "
    f"{sorted((kept_off | caps_off) - (kept_on | caps_on))}")

  # The acceptor itself, and the carbonyl carrying it, are now interior.
  for residue, name in (("48", "O"), ("48", "C")):
    assert (residue, name) in kept_on, (
      f"{residue}/{name} carries the accepting oxygen and should be in the "
      f"region")
  assert ("48", "C") in caps_off, (
    "the fixture no longer caps GLU 48's carbonyl away, so the acceptor is "
    "present either way and this proves nothing")

  # Growth continues over the peptide bond, so the next residue follows.
  for residue, name in (("49", "N"), ("49", "CA")):
    assert (residue, name) in kept_on, (
      f"{residue}/{name} should follow the acceptor over the peptide bond")


def exercise_build_adjacency():
  """Intra-ASU bonds carry the identity op on both directed edges;
  symmetry-crossing bonds carry the rt_mx forward and its inverse on the
  reverse edge."""
  builder = AtomGraphBuilder()
  identity = _canon_op(sgtbx.rt_mx())

  simple = _simple_bond_proxies(_1BQ8_FE_SPHERE_PDB)
  assert len(simple) > 0
  adj = builder.build_adjacency(simple, [], None)
  for proxy in simple:
    i_seq, j_seq = proxy.i_seqs
    assert (j_seq, identity) in adj[i_seq]
    assert (i_seq, identity) in adj[j_seq]

  asu, asu_mappings = _asu_bond_proxies(_2C2U_FE_SPHERE_PDB)
  assert len(asu) > 0
  adj = builder.build_adjacency([], asu, asu_mappings)
  n_non_identity = 0
  for proxy in asu:
    op = _canon_op(asu_mappings.get_rt_mx_ji(proxy))
    if op != identity:
      n_non_identity += 1
    assert (proxy.j_seq, op) in adj[proxy.i_seq]
    assert (proxy.i_seq, _canon_op(op.inverse())) in adj[proxy.j_seq]
  assert n_non_identity > 0


def exercise_2c2u_sym_image_provenance():
  """Symmetry-image donors of the 2C2U Fe site carry provenance back to their
  ASU parent's selection identity (chain, resseq, resname, name, altloc) + the
  operator, so a metal-ligand bond to a symmetry image can be restrained via
  geometry_restraints.edits (atom_1 vs symmetry_operation * atom_2). The identity
  is captured, not the atom object, so the hand-off does not depend on the parent
  model staying alive."""
  result = _run_endo_exo_on_string(_2C2U_FE_SPHERE_PDB)
  prov = result["sym_image_provenance"]
  atoms = list(result["model"].get_hierarchy().atoms())
  identity = sgtbx.rt_mx().as_xyz()

  assert prov, "expected symmetry-image provenance for the special-position Fe"
  for iseq, (ident, op_xyz) in prov.items():
    image = atoms[iseq]
    assert image.parent().parent().parent().id.strip() != "A"   # renamed image
    assert op_xyz != identity                                    # a real sym op
    sgtbx.rt_mx(op_xyz)                                          # parseable operator
    chain_id, resseq, resname, name, altloc = ident
    # the ASU parent is on an original chain and shares the image's identity
    assert chain_id == "A"
    assert name == image.name.strip()
    assert resseq == image.parent().parent().resseq.strip()
  # the coordinating Asp 93 / HOH 2154 images are represented
  resnames = {atoms[i].parent().resname.strip().upper() for i in prov}
  assert "ASP" in resnames and "HOH" in resnames


def exercise_special_position_water_dedup():
  """A water ON a special position must not keep the symmetry-flipped H as extra
  protons (an H4O): the redundant image is dropped whole, so the on-axis water
  survives with exactly two H. An OFF-axis water is a genuine symmetry pair and
  keeps both two-H copies."""
  result = _run_endo_exo_on_string(_SPECIAL_POS_WATER_PDB)
  waters = []
  for ch in result["model"].get_hierarchy().chains():
    for rg in ch.residue_groups():
      for ag in rg.atom_groups():
        if ag.resname.strip().upper() != "HOH":
          continue
        n_h = sum(1 for a in ag.atoms() if a.element_is_hydrogen())
        n_o = sum(1 for a in ag.atoms() if not a.element_is_hydrogen())
        waters.append((rg.resseq.strip(), n_o, n_h))

  # every surviving water is a normal 2-H water: no orphan H (O deduped away) and
  # no H4O (flipped H kept on a shared O).
  for resseq, n_o, n_h in waters:
    assert (n_o, n_h) == (1, 2), (
      f"HOH {resseq}: {n_o} O, {n_h} H (expected 1 O, 2 H)")
  counts = {}
  for resseq, _o, _h in waters:
    counts[resseq] = counts.get(resseq, 0) + 1
  assert counts.get("300") == 1, (
    f"on-axis water 300 should survive once; got {counts}")
  assert counts.get("301") == 2, (
    f"off-axis water 301 should keep both symmetry copies; got {counts}")


_ALTLOC_PDB = """\
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
ATOM      1  N   ALA A   1       5.000   5.000   5.000  1.00 10.00           N
ATOM      2  CA  ALA A   1       6.500   5.000   5.000  1.00 10.00           C
ATOM      3  C   ALA A   1       7.000   6.400   5.000  1.00 10.00           C
ATOM      4  O   ALA A   1       6.300   7.400   5.000  1.00 10.00           O
ATOM      5  H1  ALA A   1       4.703   4.672   4.228  1.00 10.00           H
ATOM      6  H2  ALA A   1       4.703   4.495   5.670  1.00 10.00           H
ATOM      7  H3  ALA A   1       4.703   5.833   5.102  1.00 10.00           H
ATOM      8  CB AALA A   1       7.000   4.200   6.200  0.30 10.00           C
ATOM      9  HA AALA A   1       6.854   4.571   4.206  1.00 10.00           H
ATOM     10  HB1AALA A   1       6.670   3.290   6.133  0.30 10.00           H
ATOM     11  HB2AALA A   1       6.670   4.612   7.014  0.30 10.00           H
ATOM     12  HB3AALA A   1       7.970   4.201   6.198  0.30 10.00           H
ATOM     13  CB BALA A   1       7.000   4.200   3.800  0.70 10.00           C
ATOM     14  HA BALA A   1       6.854   4.571   5.794  1.00 10.00           H
ATOM     15  HB1BALA A   1       6.670   4.612   2.986  0.70 10.00           H
ATOM     16  HB2BALA A   1       6.670   3.290   3.867  0.70 10.00           H
ATOM     17  HB3BALA A   1       7.970   4.201   3.802  0.70 10.00           H
END
"""


def _run_endo_exo_params(pdb_str, **overrides):
  """Drive the program with dotted PHIL overrides (``__`` for ``.``).

  Returns the full result list rather than a single region, so it can drive
  the knobs that change how many regions are produced, or that raise."""
  pdb_in = iotbx.pdb.input(source_info=None, lines=pdb_str.split("\n"))
  model = mmtbx.model.manager(model_input=pdb_in)
  dm = DataManager(["model"])
  dm.add_model("model", model)
  dm.set_default_model("model")
  master = libtbx.phil.parse(EndoexoProgram.master_phil_str)
  params = master.extract()
  params.write_files = False
  for dotted, value in overrides.items():
    scope = params
    parts = dotted.split("__")
    for part in parts[:-1]:
      scope = getattr(scope, part)
    setattr(scope, parts[-1], value)
  prog = EndoexoProgram(dm, params, master_phil=master, logger=io.StringIO())
  prog.validate()
  prog.run()
  return prog.get_results()


def region_invariant_violations(result, cap_bond=1.1, tol=1.0e-4):
  """Return a list of cap-invariant violations in a capped region *result*.

  Assumes ``capping.enable=True``.  The sub-model carries no restraints
  manager, so the anchor is identified as each cap's nearest heavy atom;
  the assertions are then against exact expected values rather than a
  bond-length threshold.  Checks:

  * every cap carries element H;
  * every cap sits *cap_bond* Angstrom from its nearest heavy atom, which is
    where :class:`HydrogenCapper` puts it;
  * the set of those anchors equals ``cap_anchor_iseqs`` exactly;
  * ``cap_iseqs`` has no duplicates, is in range, and is disjoint from
    ``seed_iseqs``;
  * ``cap_original_elements`` is parallel to ``cap_iseqs``.
  """
  # Keep the hierarchy referenced: atoms taken from a discarded one dangle.
  hierarchy = result["model"].get_hierarchy()
  atoms = list(hierarchy.atoms())
  caps = list(result["cap_iseqs"])
  bad = []

  if len(result["cap_original_elements"]) != len(caps):
    bad.append("cap_original_elements is not parallel to cap_iseqs")
  if len(set(caps)) != len(caps):
    bad.append(f"cap_iseqs has duplicates: {caps}")
  if any(ci < 0 or ci >= len(atoms) for ci in caps):
    bad.append(f"cap_iseqs out of range for {len(atoms)} atoms: {caps}")
  overlap = set(caps) & set(result["seed_iseqs"])
  if overlap:
    bad.append(f"atoms are both cap and seed: {sorted(overlap)}")

  observed_anchors = set()
  for ci in caps:
    if ci < 0 or ci >= len(atoms):
      continue
    if atoms[ci].element.strip().upper() != "H":
      bad.append(f"cap {ci} ({atoms[ci].name.strip()}) element is "
                 f"{atoms[ci].element.strip()!r}, expected H")
    heavy = sorted(((atoms[ci].distance(a), j) for j, a in enumerate(atoms)
                    if j != ci and not a.element_is_hydrogen()))
    if not heavy:
      bad.append(f"cap {ci} has no heavy atom to anchor to")
      continue
    distance, anchor = heavy[0]
    if abs(distance - cap_bond) > tol:
      bad.append(f"cap {ci} ({atoms[ci].name.strip()}) sits {distance:.3f} A "
                 f"from its nearest heavy atom, expected {cap_bond}")
    observed_anchors.add(anchor)

  if observed_anchors != set(result["cap_anchor_iseqs"]):
    bad.append(f"cap_anchor_iseqs {sorted(result['cap_anchor_iseqs'])} does "
               f"not match the anchors actually used "
               f"{sorted(observed_anchors)}")
  return bad


def exercise_capping_disabled():
  """capping.enable gates the hydrogen placement, not the bookkeeping.

  With capping off the boundary atoms stay in the region carrying their
  original element and position, leaving dangling bonds, but they are still
  reported in ``cap_iseqs`` with their anchors, and ``cap_original_elements``
  reports the element actually present."""
  for pdb_str in (_1BQ8_FE_SPHERE_PDB, _2C2U_FE_SPHERE_PDB):
    on = _run_endo_exo_params(pdb_str, capping__enable=True)[0]
    off = _run_endo_exo_params(pdb_str, capping__enable=False)[0]

    assert on["n_atoms"] == off["n_atoms"], (on["n_atoms"], off["n_atoms"])
    assert on["cap_iseqs"] == off["cap_iseqs"], "cap set must not depend on capping"
    assert on["cap_anchor_iseqs"] == off["cap_anchor_iseqs"], (
      "anchors must be recorded whether or not a hydrogen is placed")
    assert len(on["cap_iseqs"]) > 0, "fixture is expected to have caps"

    assert on["cap_original_elements"] == off["cap_original_elements"], (
      on["cap_original_elements"], off["cap_original_elements"])

    on_hier = on["model"].get_hierarchy()
    off_hier = off["model"].get_hierarchy()
    on_atoms = list(on_hier.atoms())
    off_atoms = list(off_hier.atoms())
    for k, ci in enumerate(on["cap_iseqs"]):
      assert on_atoms[ci].element.strip().upper() == "H", (
        f"capping on: cap {ci} should be H")
      assert off_atoms[ci].element.strip() == on["cap_original_elements"][k], (
        f"capping off: cap {ci} should keep its original element")


def exercise_special_position_partial_overlap():
  """A group straddling a symmetry element does not duplicate a nucleus.

  A ligand with one atom ON a rotation axis and the rest off it has an image
  that coincides only in part.  Dropping only groups that coincide ENTIRELY
  keeps such an image whole, and the shared atoms then sit two to a place,
  which no QM code can run.

  Altlocs are the reason the comparison is per parent atom group rather than
  positional: two conformers share their backbone and would otherwise cancel
  each other, so both modes are checked here as well."""
  for radius in (4.0, 5.0, 6.0):
    result = _run_endo_exo_params(
      _SPECIAL_PARTIAL_PDB, buffer__radius=radius, selection=["resname ZN"])
    atoms = list(result[0]["model"].get_hierarchy().atoms())
    assert len(atoms) > 10, f"region collapsed to {len(atoms)} atoms"
    for i, atom in enumerate(atoms):
      for other in atoms[i + 1:]:
        assert atom.distance(other) > 0.5, (
          f"{atom.name.strip()} and {other.name.strip()} coincide at "
          f"{atom.distance(other):.3f} A, r={radius}")

  # The image is dropped whole, not trimmed to the atoms that coincide.
  names = sorted(atom.name.strip() for atom in atoms)
  assert names.count("C") == 1, f"central carbon appears {names.count('C')}x"

  # Disorder is not symmetry: both conformers survive altloc=all.
  for mode, expected in (("auto", {".", "B"}), ("all", {".", "A", "B"})):
    result = _run_endo_exo_params(
      _ALTLOC_PDB, buffer__radius=4.0, altloc=mode,
      selection=["name CA"], buffer__skip_search=True)
    hierarchy = result[0]["model"].get_hierarchy()
    present = {ag.altloc.strip() or "." for ag in hierarchy.atom_groups()}
    assert present == expected, f"altloc={mode} gave {present}"


def exercise_altloc_filter():
  """_apply_altloc_filter keeps one conformer per residue group.

  'auto' takes the highest mean occupancy, a letter takes that letter, an
  absent letter falls back to the highest occupancy, and 'all' disables the
  filter.  The shared blank-altloc backbone is kept in every mode."""
  from mmtbx.geometry_restraints.endo_exo.builder import QMRegionBuilder

  def _filter(altloc):
    pdb_in = iotbx.pdb.input(source_info=None, lines=_ALTLOC_PDB.split("\n"))
    model = mmtbx.model.manager(model_input=pdb_in)
    master = libtbx.phil.parse(EndoexoProgram.master_phil_str)
    params = master.extract()
    params.altloc = altloc
    builder = QMRegionBuilder(params, logger=io.StringIO())
    return builder._apply_altloc_filter(model).get_hierarchy()

  def _cb_altlocs(altloc):
    hierarchy = _filter(altloc)
    return sorted(a.parent().altloc.strip() for a in hierarchy.atoms()
                  if a.name.strip() == "CB")

  # CB is present as altloc A at occupancy 0.30 and B at 0.70.  B is listed
  # second, so picking it distinguishes "highest occupancy" from "first".
  assert _cb_altlocs("auto") == ["B"], _cb_altlocs("auto")
  assert _cb_altlocs("A") == ["A"], _cb_altlocs("A")
  assert _cb_altlocs("B") == ["B"], _cb_altlocs("B")
  assert _cb_altlocs("all") == ["A", "B"], _cb_altlocs("all")
  assert _cb_altlocs("Z") == ["B"], _cb_altlocs("Z")

  for mode in ("auto", "A", "B", "all", "Z"):
    hierarchy = _filter(mode)
    names = {a.name.strip() for a in hierarchy.atoms()}
    assert {"N", "CA", "C", "O"} <= names, (mode, names)


def exercise_element_filter():
  """element_filter restricts which metals seed a region; a miss raises Sorry.

  Driven on a two-metal model so that filtering to one element is
  distinguishable from ignoring the filter."""
  from libtbx.utils import Sorry

  # A Zn well away from the Fe site, so it seeds its own region.
  zn = ("HETATM  700 ZN    ZN A  56      20.000  20.000  20.000"
        "  1.00 10.00          ZN\n")
  two_metals = _1BQ8_FE_SPHERE_PDB.replace("END\n", zn + "END\n") \
    if "END\n" in _1BQ8_FE_SPHERE_PDB else _1BQ8_FE_SPHERE_PDB + zn

  def _seed_elements(results):
    out = []
    for res in results:
      hierarchy = res["model"].get_hierarchy()
      atoms = list(hierarchy.atoms())
      out.append(sorted(atoms[i].element.strip().upper()
                        for i in res["seed_iseqs"]))
    return sorted(out)

  both = _run_endo_exo_params(two_metals)
  assert len(both) == 2, f"expected one region per metal, got {len(both)}"
  assert _seed_elements(both) == [["FE"], ["ZN"]], _seed_elements(both)

  only_fe = _run_endo_exo_params(two_metals, element_filter=["FE"])
  assert len(only_fe) == 1, f"element_filter=FE should give 1 region, got {len(only_fe)}"
  assert _seed_elements(only_fe) == [["FE"]], _seed_elements(only_fe)

  only_zn = _run_endo_exo_params(two_metals, element_filter=["ZN"])
  assert len(only_zn) == 1, f"element_filter=ZN should give 1 region, got {len(only_zn)}"
  assert _seed_elements(only_zn) == [["ZN"]], _seed_elements(only_zn)

  try:
    _run_endo_exo_params(_1BQ8_FE_SPHERE_PDB, element_filter=["ZN"])
  except Sorry:
    pass
  else:
    raise AssertionError("element_filter=ZN on an Fe-only model should Sorry")


def exercise_depth_and_skip_search():
  """max_search_depth does not bound the region, and a metal seed alone
  cannot expand.

  Growth is bounded by the cut rules rather than by a hop count, so every
  ``max_search_depth`` yields the same region.  With ``buffer.skip_search``
  the seed is the metal alone, and metal coordination bonds are absent from
  the graph (``link_metals=False``), so there is nothing to expand through."""
  for pdb_str in (_1BQ8_FE_SPHERE_PDB, _2C2U_FE_SPHERE_PDB):
    sizes = {_run_endo_exo_params(pdb_str, max_search_depth=d)[0]["n_atoms"]
             for d in (0, 1, 2, 3, 10)}
    assert len(sizes) == 1, f"max_search_depth changed the region: {sizes}"

    skipped = _run_endo_exo_params(pdb_str, buffer__skip_search=True)[0]
    # The region is exactly the seed set: nothing was added by the radius
    # search, and a metal has no covalent edge to expand along.
    assert sorted(skipped["seed_iseqs"]) == list(range(skipped["n_atoms"])), (
      skipped["seed_iseqs"], skipped["n_atoms"])
    assert skipped["cap_iseqs"] == [], skipped["cap_iseqs"]


def exercise_write_files_roundtrip():
  """write_files emits a PDB, an mmCIF and a sidecar that re-parses.

  The sidecar indices address the written atoms, and cap_atoms point at
  hydrogens."""
  import os
  import shutil
  import tempfile

  from mmtbx.geometry_restraints.endo_exo.builder import master_sidecar_phil_str

  tmp = tempfile.mkdtemp(prefix="tst_endo_exo_")
  cwd = os.getcwd()
  try:
    os.chdir(tmp)
    result = _run_endo_exo_params(_1BQ8_FE_SPHERE_PDB, write_files=True)[0]
    stem = os.path.basename(result["file_name"])
    for ext in (".pdb", ".cif", ".phil"):
      assert os.path.isfile(stem + ext), f"missing {stem + ext}"

    # Fetch against the schema so the values come back typed, which also
    # checks the sidecar carries nothing the schema does not declare.
    schema = libtbx.phil.parse(master_sidecar_phil_str)
    written_phil = libtbx.phil.parse(file_name=stem + ".phil")
    unused = []
    region = schema.fetch(
      source=written_phil,
      track_unused_definitions=unused).extract().endo_exo_region
    assert not unused, f"sidecar has undeclared parameters: {unused}"
    assert list(region.cap_atoms) == list(result["cap_iseqs"]), (
      region.cap_atoms, result["cap_iseqs"])
    assert list(region.seed_atoms) == list(result["seed_iseqs"]), (
      region.seed_atoms, result["seed_iseqs"])
    assert list(region.cap_original_elements) == list(
      result["cap_original_elements"])
    assert region.selection_string is None

    pdb_in = iotbx.pdb.input(file_name=stem + ".pdb")
    written = pdb_in.construct_hierarchy()
    watoms = list(written.atoms())
    assert len(watoms) == result["n_atoms"], (len(watoms), result["n_atoms"])
    assert pdb_in.crystal_symmetry() is not None, "written PDB lost CRYST1"
    for ci in region.cap_atoms:
      assert watoms[ci].element.strip().upper() == "H", (
        f"sidecar cap index {ci} does not address a hydrogen in the PDB")

    # The mmCIF copy exists so the capped atoms keep their element through a
    # re-parse rather than being re-derived from their (heavy-atom) names.
    cif_in = iotbx.pdb.input(file_name=stem + ".cif")
    cif_hier = cif_in.construct_hierarchy()
    catoms = list(cif_hier.atoms())
    assert len(catoms) == result["n_atoms"], (len(catoms), result["n_atoms"])
    assert ([a.element.strip().upper() for a in catoms]
            == [a.element.strip().upper() for a in watoms]), (
      "mmCIF and PDB copies disagree on elements")
    for ci in region.cap_atoms:
      assert catoms[ci].element.strip().upper() == "H", (
        f"sidecar cap index {ci} does not address a hydrogen in the mmCIF")

    # A region with no caps must still write a sidecar that re-parses: an
    # empty list cannot be expressed in PHIL, so it is written as None.
    capless = _run_endo_exo_params(
      _1BQ8_FE_SPHERE_PDB, write_files=True, buffer__skip_search=True)[0]
    capless_stem = os.path.basename(capless["file_name"])
    assert capless["cap_iseqs"] == [], capless["cap_iseqs"]
    capless_region = schema.fetch(source=libtbx.phil.parse(
      file_name=capless_stem + ".phil")).extract().endo_exo_region
    assert capless_region.cap_atoms is None, capless_region.cap_atoms
    assert capless_region.seed_atoms is not None
  finally:
    os.chdir(cwd)
    shutil.rmtree(tmp, ignore_errors=True)


def exercise_region_invariants():
  """Every capped region satisfies the cap invariants across radii and both
  preferred-cut settings."""
  for name, pdb_str in (("1BQ8", _1BQ8_FE_SPHERE_PDB),
                        ("2C2U", _2C2U_FE_SPHERE_PDB)):
    for radius in (4.0, 5.0, 6.0, 7.0):
      for preferred in (True, False):
        result = _run_endo_exo_params(
          pdb_str, buffer__radius=radius,
          capping__preferred_cuts=preferred)[0]
        bad = region_invariant_violations(result)
        assert not bad, (
          f"{name} radius={radius} preferred_cuts={preferred}: "
          + "; ".join(bad))


# One alanine per asymmetric unit in a 3.6 A cell: the chain leaves the
# unit through a lattice translation, so the residue after it is itself.
_POLYMER_PDB = """\
CRYST1    3.600   20.000   20.000  90.00  90.00  90.00 P 1
ATOM      1 N    ALA A   1      -1.200   0.000   0.000  1.00 20.00           N
ATOM      2 CA   ALA A   1       0.000   0.000   0.000  1.00 20.00           C
ATOM      3 C    ALA A   1       1.000   1.050   0.000  1.00 20.00           C
ATOM      4 O    ALA A   1       0.700   2.220   0.000  1.00 20.00           O
ATOM      5 CB   ALA A   1       0.550  -0.800  -1.200  1.00 20.00           C
END
"""


# Tris on the 2-fold of P 1 2 1, central C on the axis and every other
# heavy atom off it, so the symmetry image coincides only in part.
_SPECIAL_PARTIAL_PDB = """\
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1 2 1
HETATM    1 N    TRS A 401       0.849   5.849   0.849  1.00 20.00           N
HETATM    2 C    TRS A 401       0.000   5.000   0.000  1.00 20.00           C
HETATM    3 C1   TRS A 401       0.883   4.117  -0.883  1.00 20.00           C
HETATM    4 C2   TRS A 401      -0.883   5.883  -0.883  1.00 20.00           C
HETATM    5 C3   TRS A 401      -0.883   4.117   0.883  1.00 20.00           C
HETATM    6 O1   TRS A 401       1.709   3.291  -1.709  1.00 20.00           O
HETATM    7 O2   TRS A 401      -1.709   6.709  -1.709  1.00 20.00           O
HETATM    8 O3   TRS A 401      -1.709   3.291   1.709  1.00 20.00           O
HETATM    9  H11 TRS A 401       1.118   4.980  -0.594  1.00 20.00           H
HETATM   10  H12 TRS A 401      -0.020   3.922  -0.504  1.00 20.00           H
HETATM   11  H21 TRS A 401       0.056   6.024  -0.677  1.00 20.00           H
HETATM   12  H22 TRS A 401      -1.156   5.076  -0.423  1.00 20.00           H
HETATM   13  H31 TRS A 401       0.024   3.922   0.579  1.00 20.00           H
HETATM   14  H32 TRS A 401      -1.124   4.977   0.522  1.00 20.00           H
HETATM   15  HN1 TRS A 401       0.366   6.367   1.388  1.00 20.00           H
HETATM   16  HN2 TRS A 401       1.367   6.388   0.366  1.00 20.00           H
HETATM   17  HN3 TRS A 401       1.388   5.366   1.367  1.00 20.00           H
HETATM   18  HO1 TRS A 401       1.284   2.587  -1.755  1.00 20.00           H
HETATM   19  HO2 TRS A 401      -1.250   7.372  -1.979  1.00 20.00           H
HETATM   20  HO3 TRS A 401      -2.466   3.313   1.364  1.00 20.00           H
HETATM   21 ZN    ZN A 402       2.061   7.061   2.061  1.00 20.00          ZN
END
"""


# Tris on a zinc: a ligand that names an atom C and an atom N, both
# bonded entirely within the residue.  Protonated with mmtbx.reduce2.
_TRS_LIGAND_PDB = """\
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
HETATM    1 N    TRS A 401       0.849   0.849   0.849  1.00 20.00           N
HETATM    2 C    TRS A 401       0.000   0.000   0.000  1.00 20.00           C
HETATM    3 C1   TRS A 401       0.883  -0.883  -0.883  1.00 20.00           C
HETATM    4 C2   TRS A 401      -0.883   0.883  -0.883  1.00 20.00           C
HETATM    5 C3   TRS A 401      -0.883  -0.883   0.883  1.00 20.00           C
HETATM    6 O1   TRS A 401       1.709  -1.709  -1.709  1.00 20.00           O
HETATM    7 O2   TRS A 401      -1.709   1.709  -1.709  1.00 20.00           O
HETATM    8 O3   TRS A 401      -1.709  -1.709   1.709  1.00 20.00           O
HETATM    9  H11 TRS A 401       1.118  -0.020  -0.594  1.00 20.00           H
HETATM   10  H12 TRS A 401      -0.020  -1.078  -0.504  1.00 20.00           H
HETATM   11  H21 TRS A 401       0.056   1.024  -0.677  1.00 20.00           H
HETATM   12  H22 TRS A 401      -1.156   0.076  -0.423  1.00 20.00           H
HETATM   13  H31 TRS A 401       0.024  -1.078   0.579  1.00 20.00           H
HETATM   14  H32 TRS A 401      -1.124  -0.023   0.522  1.00 20.00           H
HETATM   15  HN1 TRS A 401       0.366   1.367   1.388  1.00 20.00           H
HETATM   16  HN2 TRS A 401       1.367   1.388   0.366  1.00 20.00           H
HETATM   17  HN3 TRS A 401       1.388   0.366   1.367  1.00 20.00           H
HETATM   18  HO1 TRS A 401       1.284  -2.413  -1.755  1.00 20.00           H
HETATM   19  HO2 TRS A 401      -1.250   2.372  -1.979  1.00 20.00           H
HETATM   20  HO3 TRS A 401      -2.466  -1.687   1.364  1.00 20.00           H
HETATM   21 ZN    ZN A 402       2.061   2.061   2.061  1.00 20.00          ZN
END
"""


_TRIPEPTIDE_PDB = """\
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00 10.00           N
ATOM      2  CA  ALA A   1      11.458  10.000  10.000  1.00 10.00           C
ATOM      3  C   ALA A   1      12.009  11.420  10.000  1.00 10.00           C
ATOM      4  O   ALA A   1      11.257  12.394  10.000  1.00 10.00           O
ATOM      5  CB  ALA A   1      11.987   9.230  11.207  1.00 10.00           C
ATOM      6  H1  ALA A   1       9.703  10.451   9.293  1.00 10.00           H
ATOM      7  H2  ALA A   1       9.703   9.162   9.963  1.00 10.00           H
ATOM      8  H3  ALA A   1       9.703  10.387  10.745  1.00 10.00           H
ATOM      9  HA  ALA A   1      11.770   9.558   9.195  1.00 10.00           H
ATOM     10  HB1 ALA A   1      11.666   8.316  11.164  1.00 10.00           H
ATOM     11  HB2 ALA A   1      11.666   9.655  12.018  1.00 10.00           H
ATOM     12  HB3 ALA A   1      12.957   9.242  11.189  1.00 10.00           H
ATOM     13  N   CYS A   2      13.339  11.540  10.000  1.00 10.00           N
ATOM     14  CA  CYS A   2      13.998  12.840  10.000  1.00 10.00           C
ATOM     15  C   CYS A   2      15.508  12.660  10.000  1.00 10.00           C
ATOM     16  O   CYS A   2      16.010  11.535  10.000  1.00 10.00           O
ATOM     17  CB  CYS A   2      13.560  13.680  11.200  1.00 10.00           C
ATOM     18  SG  CYS A   2      14.100  15.400  11.200  1.00 10.00           S
ATOM     19  H   CYS A   2      13.884  10.875  10.000  1.00 10.00           H
ATOM     20  HA  CYS A   2      13.743  13.329   9.202  1.00 10.00           H
ATOM     21  HB2 CYS A   2      13.915  13.268  12.004  1.00 10.00           H
ATOM     22  HB3 CYS A   2      12.590  13.684  11.229  1.00 10.00           H
ATOM     23  HG  CYS A   2      15.197  15.465  11.682  1.00 10.00           H
ATOM     24  N   ALA A   3      16.230  13.780  10.000  1.00 10.00           N
ATOM     25  CA  ALA A   3      17.690  13.790  10.000  1.00 10.00           C
ATOM     26  C   ALA A   3      18.220  15.220  10.000  1.00 10.00           C
ATOM     27  O   ALA A   3      17.450  16.180  10.000  1.00 10.00           O
ATOM     28  CB  ALA A   3      18.230  13.010  11.200  1.00 10.00           C
ATOM     29  OXT ALA A   3      19.450  15.380  10.000  1.00 10.00           O
ATOM     30  H   ALA A   3      15.887  14.569  10.000  1.00 10.00           H
ATOM     31  HA  ALA A   3      18.012  13.352   9.197  1.00 10.00           H
ATOM     32  HB1 ALA A   3      17.917  12.093  11.148  1.00 10.00           H
ATOM     33  HB2 ALA A   3      17.908  13.423  12.016  1.00 10.00           H
ATOM     34  HB3 ALA A   3      19.200  13.030  11.179  1.00 10.00           H
END
"""


def _atom_names_by_residue(result):
  """Return ``{resseq: set of atom names}`` for a region."""
  hierarchy = result["model"].get_hierarchy()
  out = {}
  for atom in hierarchy.atoms():
    resseq = atom.parent().parent().resseq.strip()
    out.setdefault(resseq, set()).add(atom.name.strip())
  return out


def _protonated_graph(pdb_str):
  """Return ``(model, atoms, adjacency)`` for *pdb_str*.

  The fixtures carry hydrogens already, so nothing is placed here.  The
  model is returned so callers hold it: atoms taken from a hierarchy
  nothing references dangle."""
  pdb_in = iotbx.pdb.input(source_info=None, lines=pdb_str.split("\n"))
  model = mmtbx.model.manager(model_input=pdb_in)
  model.process(make_restraints=True)

  geometry = model.get_restraints_manager().geometry
  sites_cart = model.get_sites_cart()
  simple, asu = geometry.get_all_bond_proxies(sites_cart=sites_cart)
  asu_mappings = geometry.pair_proxies(
    sites_cart=sites_cart).bond_proxies.asu_mappings()
  adjacency = AtomGraphBuilder().build_adjacency(simple, asu, asu_mappings)
  return model, list(model.get_hierarchy().atoms()), adjacency


def exercise_radius_only_grows_the_region():
  """Raising ``buffer.radius`` adds atoms and never takes any away.

  The knob is swept throughout these tests but nothing pins what it does.
  A cut rule firing on the wider sphere but not the narrower one would
  drop atoms as the radius grew, which a plain count hides and a subset
  check does not."""
  for pdb_str, overrides in (
      (_1BQ8_FE_SPHERE_PDB, {}),):
    previous, previous_interior = None, None
    # 5.0 -> 6.0 is the step that matters: it is where a guard keyed on the
    # seed set drops TRP 37 as the sphere widens.
    for radius in (4.0, 5.0, 6.0, 8.0):
      results = _run_endo_exo_params(
        pdb_str, buffer__radius=radius, **overrides)
      present, interior = set(), set()
      for result in results:
        hierarchy = result["model"].get_hierarchy()
        caps = set(result["cap_iseqs"])
        for i_seq, atom in enumerate(hierarchy.atoms()):
          residue_group = atom.parent().parent()
          key = (residue_group.parent().id.strip(),
                 residue_group.resseq.strip(),
                 atom.parent().resname.strip(),
                 atom.name.strip())
          present.add(key)
          if i_seq not in caps:
            interior.add(key)
      if previous is not None:
        lost = previous - present
        assert not lost, (
          f"radius {radius} dropped {len(lost)} atom(s) kept at the smaller "
          f"radius, e.g. {sorted(lost)[:4]}")
        # Capping leaves the name alone, so an atom demoted from interior
        # to a cap still matches above; that is a loss of chemistry rather
        # than of an atom, and has to be looked for separately.
        demoted = previous_interior - interior
        assert not demoted, (
          f"radius {radius} turned {len(demoted)} interior atom(s) into "
          f"caps, e.g. {sorted(demoted)[:4]}")
      previous, previous_interior = present, interior


def exercise_region_does_not_depend_on_queue_order():
  """The region is the same however BFS orders its queue.

  This has to reorder the queue itself.  Handing ``grow_by_depth`` the
  seeds in a different order proves nothing: it collects them into a set on
  the way in, so the ordering is gone before BFS starts, and a test written
  that way passes against code that is plainly racy.

  Shuffled rather than rotated.  ``visited`` is pre-seeded with the whole
  seed set, hundreds of nodes, so rotating the start point leaves relative
  order intact and barely perturbs anything: three rotations collapse to
  two behaviours.  The shuffles are seeded, so a failure reproduces exactly.

  Radius 8 on 1BQ8 is the configuration that ever varied, and 2C2U is here
  for the composed-operator path rather than for a race of its own.  Each
  region costs seconds to grow, so the sweep is kept to those.

  Compared on the grower's own state.  The materialised region cannot show
  which anchor a cap belongs to -- ``cap_anchor_iseqs`` is deduplicated, so
  it is not parallel to ``cap_iseqs`` -- and a cap swapping anchors leaves
  every atom count untouched."""
  from mmtbx.geometry_restraints.endo_exo.grow import QMRegionGrower

  original = QMRegionGrower._grow_until_exhausted
  grow_by_depth = QMRegionGrower.grow_by_depth
  captured = []
  trial = [0]

  def reordered(self, queue, *args, **kwargs):
    rng = random.Random(trial[0])
    items = list(queue)
    rng.shuffle(items)
    queue.clear()
    queue.extend(items)
    return original(self, queue, *args, **kwargs)

  def capturing(self, seeds, adjacency, model, *args, **kwargs):
    visited, caps = grow_by_depth(
      self, seeds, adjacency, model, *args, **kwargs)
    captured.append((
      tuple(sorted((iseq, op.as_xyz()) for (iseq, op) in visited)),
      tuple(sorted(((cap[0], cap[1].as_xyz()),
                    (anchor[0], anchor[1].as_xyz()))
                   for cap, anchor in caps.items())),
    ))
    return visited, caps

  outcomes = {}
  QMRegionGrower._grow_until_exhausted = reordered
  QMRegionGrower.grow_by_depth = capturing
  try:
    for index in range(4):
      trial[0] = index
      captured.clear()
      for pdb_str, radius in ((_1BQ8_FE_SPHERE_PDB, 8.0),
                              (_2C2U_FE_SPHERE_PDB, 7.0)):
        _run_endo_exo_params(pdb_str, buffer__radius=radius)
      outcomes.setdefault(tuple(captured), []).append(index)
  finally:
    QMRegionGrower._grow_until_exhausted = original
    QMRegionGrower.grow_by_depth = grow_by_depth

  assert outcomes, "no region was grown"
  assert len(outcomes) == 1, (
    f"queue order changed the region: {len(outcomes)} distinct outcomes, "
    f"grouped as {sorted(outcomes.values())}")


def exercise_amide_unit_is_never_split():
  """A retained backbone carbonyl keeps both its oxygen and its amide
  nitrogen.

  This is what makes a truncated backbone come out as CH3CONHCH3. Neither
  C-O nor the peptide C-N is a cut site, so BFS reaching a carbonyl carbon
  always continues through both: the amide is carried whole or left out
  entirely, and the cuts fall on CA-C and CA-N either side of it. A cap
  appearing on the oxygen or the nitrogen would mean an amide carved in
  half, leaving an aldehyde or a bare amine in the QM region."""
  checked = 0
  for pdb_str, overrides in (
      (_1BQ8_FE_SPHERE_PDB, {}),
      (_2C2U_FE_SPHERE_PDB, {})):
    for radius in (5.0, 7.0):
      for result in _run_endo_exo_params(
          pdb_str, buffer__radius=radius, **overrides):
        hierarchy = result["model"].get_hierarchy()
        atoms = list(hierarchy.atoms())
        caps = set(result["cap_iseqs"])
        # Capping sets the element to H but leaves the name alone, so a
        # severed nitrogen still reads as "N" and no longer reads as
        # nitrogen.  Ask what each atom WAS, or the check walks straight
        # past the atoms it exists to catch.
        was = {i_seq: atom.element.strip().upper()
               for i_seq, atom in enumerate(atoms)}
        was.update({
          i_seq: element.strip().upper() for i_seq, element
          in zip(result["cap_iseqs"], result["cap_original_elements"])})
        partners = [i_seq for i_seq, element in was.items()
                    if element in ("O", "N")]

        for i_seq, atom in enumerate(atoms):
          if i_seq in caps or atom.name.strip() != "C":
            continue
          if atom.parent().resname.strip() == "HOH":
            continue
          oxygens, nitrogens = [], []
          for j_seq in partners:
            # A cap sits 1.1 A from its anchor, so it stays within reach
            # of the carbonyl it was cut from and is still found here.
            if j_seq == i_seq or atom.distance(atoms[j_seq]) > 1.45:
              continue
            (oxygens if was[j_seq] == "O" else nitrogens).append(j_seq)
          where = (f"{atom.parent().resname.strip()}"
                   f"{atom.parent().parent().resseq.strip()} at radius "
                   f"{radius}")
          assert any(j not in caps for j in oxygens), (
            f"{where}: carbonyl oxygen capped away, leaving a bare carbon")
          # A C-terminus has no amide nitrogen at all, which is fine; one
          # that exists must not have been replaced by a hydrogen.
          assert not nitrogens or any(j not in caps for j in nitrogens), (
            f"{where}: amide nitrogen capped away, splitting the peptide "
            f"unit")
          checked += len(nitrogens)

  assert checked > 30, (
    f"only {checked} carbonyl-nitrogen pairs examined; the nitrogen half "
    f"of this check is not exercising anything")


# A proline whose preceding residue is absent: its ring nitrogen holds
# CA and CD and no hydrogen at all.  CD is a carbon bearing no oxygen,
# so it must not be mistaken for an acyl group.
_PROLINE_GAP_PDB = """\
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
ATOM      1  N   PRO A   2       8.872   7.421   8.563  1.00  3.42           N
ATOM      2  CA  PRO A   2       8.561   7.954   9.916  1.00  3.79           C
ATOM      3  C   PRO A   2       8.613   9.480   9.954  1.00  3.62           C
ATOM      4  O   PRO A   2       8.148  10.060  10.929  1.00  5.67           O
ATOM      5  CB  PRO A   2       9.655   7.382  10.800  1.00  4.40           C
ATOM      6  CG  PRO A   2      10.750   7.051   9.885  1.00  5.10           C
ATOM      7  CD  PRO A   2      10.139   6.613   8.599  1.00  3.49           C
ATOM      8  HA  PRO A   2       7.697   7.631  10.217  1.00  3.79           H
ATOM      9  HB2 PRO A   2       9.935   8.046  11.450  1.00  4.40           H
ATOM     10  HB3 PRO A   2       9.331   6.589  11.254  1.00  4.40           H
ATOM     11  HG2 PRO A   2      11.301   7.837   9.747  1.00  5.10           H
ATOM     12  HG3 PRO A   2      11.282   6.336  10.267  1.00  5.10           H
ATOM     13  HD2 PRO A   2      10.713   6.829   7.847  1.00  3.49           H
ATOM     14  HD3 PRO A   2       9.948   5.662   8.608  1.00  3.49           H
ATOM     15  N   ALA A   3       9.180  10.100   8.930  1.00  3.50           N
ATOM     16  CA  ALA A   3       9.290  11.550   8.850  1.00  3.60           C
ATOM     17  C   ALA A   3      10.010  11.960   7.570  1.00  3.70           C
ATOM     18  O   ALA A   3      10.400  11.100   6.780  1.00  3.80           O
ATOM     19  CB  ALA A   3       7.910  12.190   8.880  1.00  3.90           C
ATOM     20  OXT ALA A   3      10.190  13.160   7.330  1.00  3.80           O
ATOM     21  H   ALA A   3       9.518   9.694   8.251  1.00  3.50           H
ATOM     22  HA  ALA A   3       9.795  11.869   9.615  1.00  3.60           H
ATOM     23  HB1 ALA A   3       7.468  11.948   9.709  1.00  3.90           H
ATOM     24  HB2 ALA A   3       7.395  11.867   8.124  1.00  3.90           H
ATOM     25  HB3 ALA A   3       8.008  13.154   8.826  1.00  3.90           H
END
"""


# An N-terminal proline, whose secondary nitrogen carries two hydrogens
# rather than the three of a primary amine.
_PROLINE_NTERM_PDB = """\
CRYST1   40.000   40.000   40.000  90.00  90.00  90.00 P 1
ATOM      1  N   PRO A   1       8.872   7.421   8.563  1.00  3.42           N
ATOM      2  CA  PRO A   1       8.561   7.954   9.916  1.00  3.79           C
ATOM      3  C   PRO A   1       8.613   9.480   9.954  1.00  3.62           C
ATOM      4  O   PRO A   1       8.148  10.060  10.929  1.00  5.67           O
ATOM      5  CB  PRO A   1       9.655   7.382  10.800  1.00  4.40           C
ATOM      6  CG  PRO A   1      10.750   7.051   9.885  1.00  5.10           C
ATOM      7  CD  PRO A   1      10.139   6.613   8.599  1.00  3.49           C
ATOM      8  H2  PRO A   1       8.149   6.863   8.269  1.00  3.42           H
ATOM      9  H3  PRO A   1       8.981   8.153   7.952  1.00  3.42           H
ATOM     10  HA  PRO A   1       7.697   7.631  10.217  1.00  3.79           H
ATOM     11  HB2 PRO A   1       9.935   8.046  11.450  1.00  4.40           H
ATOM     12  HB3 PRO A   1       9.331   6.589  11.254  1.00  4.40           H
ATOM     13  HG2 PRO A   1      11.301   7.837   9.747  1.00  5.10           H
ATOM     14  HG3 PRO A   1      11.282   6.336  10.267  1.00  5.10           H
ATOM     15  HD2 PRO A   1      10.713   6.829   7.847  1.00  3.49           H
ATOM     16  HD3 PRO A   1       9.948   5.662   8.608  1.00  3.49           H
ATOM     17  N   ALA A   2       9.180  10.100   8.930  1.00  3.50           N
ATOM     18  CA  ALA A   2       9.290  11.550   8.850  1.00  3.60           C
ATOM     19  C   ALA A   2      10.010  11.960   7.570  1.00  3.70           C
ATOM     20  O   ALA A   2      10.400  11.100   6.780  1.00  3.80           O
ATOM     21  CB  ALA A   2       7.910  12.190   8.880  1.00  3.90           C
ATOM     22  H   ALA A   2       9.518   9.694   8.251  1.00  3.50           H
ATOM     23  HA  ALA A   2       9.795  11.869   9.615  1.00  3.60           H
ATOM     24  HB1 ALA A   2       7.468  11.948   9.709  1.00  3.90           H
ATOM     25  HB2 ALA A   2       7.395  11.867   8.124  1.00  3.90           H
ATOM     26  HB3 ALA A   2       8.008  13.154   8.826  1.00  3.90           H
END
"""


def exercise_terminal_amine_covers_proline():
  """An N-terminal proline is protected too.

  Its nitrogen is secondary, carrying CD as well as CA within its own
  residue, so a free proline terminus holds two hydrogens where a primary
  amine holds three.  A guard demanding three would quietly stop protecting
  every proline terminus, and no other fixture here has one."""
  _model, atoms, adjacency = _protonated_graph(_PROLINE_NTERM_PDB)

  def nitrogen(resseq):
    return next(i_seq for i_seq, atom in enumerate(atoms)
                if atom.parent().parent().resseq.strip() == resseq
                and atom.name.strip() == "N")

  proline_n = nitrogen("1")
  hydrogens = sum(1 for j_seq, _op in adjacency[proline_n]
                  if atoms[j_seq].element_is_hydrogen())
  assert hydrogens == 2, (
    f"fixture changed: the proline terminus carries {hydrogens} hydrogens, "
    f"expected 2, so it no longer pins the threshold")

  assert QMRegionGrower._is_terminal_amine(
    proline_n, adjacency, atoms) is True, (
    "an N-terminal proline was not recognised as a chain terminus")
  assert QMRegionGrower._is_terminal_amine(
    nitrogen("2"), adjacency, atoms) is False, (
    "the mid-chain nitrogen following proline was called a terminus")


def exercise_chain_break_amine_is_not_any_nitrogen():
  """The chain-break guard covers a backbone nitrogen whose predecessor is
  absent, and nothing else.

  The predicate decides a cut on its own, without consulting the growing
  region, so nothing downstream will catch it being too broad.  Three ways
  to get it wrong: a sidechain nitrogen carrying one hydrogen (TRP NE1)
  reads as a chain break unless the CA bond is required; a mid-chain
  nitrogen reads as one unless the preceding carbonyl vetoes it; and a real
  N-terminus reads as one unless the hydrogen count separates them."""
  _model, atoms, adjacency = _protonated_graph(_1BQ8_FE_SPHERE_PDB)

  def nitrogen(resseq, name):
    return next(i_seq for i_seq, atom in enumerate(atoms)
                if atom.parent().parent().resseq.strip() == resseq
                and atom.name.strip() == name)

  # 1BQ8 runs 5-12, 37-45, 48-49: three backbone nitrogens have no
  # preceding residue in the file.
  for resseq in ("5", "37", "48"):
    assert QMRegionGrower._is_chain_break_amine(
      nitrogen(resseq, "N"), adjacency, atoms) is True, (
      f"residue {resseq} starts a chain break and was not recognised")

  # Mid-chain nitrogens carry the preceding carbonyl and are not breaks.
  # 40 is proline, whose N also carries CD within its own residue.
  for resseq in ("6", "7", "40", "45", "49"):
    assert QMRegionGrower._is_chain_break_amine(
      nitrogen(resseq, "N"), adjacency, atoms) is False, (
      f"residue {resseq} is mid-chain, its predecessor is present")

  # Sidechain nitrogens are not backbone nitrogens, whatever their
  # hydrogen count.  TRP NE1 carries exactly one, as an amide does, and is
  # the case the CA bond requirement exists for.
  for resseq, name in (("37", "NE1"), ("7", "NZ")):
    assert QMRegionGrower._is_chain_break_amine(
      nitrogen(resseq, name), adjacency, atoms) is False, (
      f"sidechain nitrogen {resseq}/{name} read as a chain break")

  # Nor is anything that is not a nitrogen.  A backbone carbonyl carbon
  # satisfies every other clause -- bonded to its own CA, carrying no
  # hydrogen, and its one neighbour in another residue is named N rather
  # than C, so the peptide bond does not veto it.
  for resseq in ("6", "48"):
    assert QMRegionGrower._is_chain_break_amine(
      nitrogen(resseq, "C"), adjacency, atoms) is False, (
      f"carbonyl carbon {resseq}/C read as a chain-break nitrogen")

  # An acyl group is not always the previous residue's atom named C.  Read
  # from the ideal monomers so the N-acylated residues are covered whether
  # or not a fixture happens to contain one: FME names its formyl carbon
  # CN, AYA its acetyl CT, SAC its C1A, and all three sit in the residue's
  # own atom group.  Each is a complete amide.
  monomers = _monomer_server()
  identity = _canon_op(sgtbx.rt_mx())
  for code, expected in (("FME", False), ("AYA", False), ("SAC", False),
                         ("MSE", True), ("ALA", True)):
    comp = monomers.get_comp_comp_id_direct(code)
    root = iotbx.pdb.hierarchy.root()
    hier_model = iotbx.pdb.hierarchy.model()
    root.append_model(hier_model)
    chain = iotbx.pdb.hierarchy.chain(id="A")
    hier_model.append_chain(chain)
    residue_group = iotbx.pdb.hierarchy.residue_group(resseq="   1")
    chain.append_residue_group(residue_group)
    atom_group = iotbx.pdb.hierarchy.atom_group(resname=code)
    residue_group.append_atom_group(atom_group)
    index = {}
    for position, entry in enumerate(comp.atom_list):
      atom = iotbx.pdb.hierarchy.atom()
      atom.name = entry.atom_id.strip().ljust(4)
      atom.element = entry.type_symbol.strip().rjust(2)
      atom.xyz = (position * 2.0, 0, 0)
      atom_group.append_atom(atom)
      index[entry.atom_id.strip()] = position
    ideal_atoms = list(root.atoms())
    ideal_adjacency = defaultdict(set)
    for bond in comp.bond_list:
      first, second = bond.atom_id_1.strip(), bond.atom_id_2.strip()
      if first in index and second in index:
        ideal_adjacency[index[first]].add((index[second], identity))
        ideal_adjacency[index[second]].add((index[first], identity))
    for position, atom in enumerate(ideal_atoms):
      if atom.name.strip().upper() != "N":
        continue
      assert QMRegionGrower._is_chain_break_amine(
        position, ideal_adjacency, ideal_atoms) is expected, (
        f"{code} nitrogen read as chain break {not expected}")

  # Hydrogens are counted by element, not by name.  A neutron structure
  # names them D and this is the same terminus, so a name-based count would
  # see none and cut a real amine away.
  deuterated = []
  for line in _TRIPEPTIDE_PDB.split("\n"):
    if line.startswith(("ATOM", "HETATM")) and line[76:78].strip() == "H":
      line = (line[:12] + line[12:16].replace("H", "D", 1)
              + line[16:76] + " D" + line[78:])
    deuterated.append(line)
  _m5, atoms5, adjacency5 = _protonated_graph("\n".join(deuterated))
  heavy_water_n = next(
    i for i, a in enumerate(atoms5)
    if a.parent().parent().resseq.strip() == "1" and a.name.strip() == "N")
  assert QMRegionGrower._is_chain_break_amine(
    heavy_water_n, adjacency5, atoms5) is False, (
    "a deuterated N-terminus read as a chain break")

  # A ring carbon is not an acyl carbon.  A proline at a chain break holds
  # CA and CD and no hydrogen, and CD carries no oxygen, so the nitrogen is
  # genuinely short and the ring must not be read as its amide.
  _m4, atoms4, adjacency4 = _protonated_graph(_PROLINE_GAP_PDB)
  gap_proline_n = next(
    i for i, a in enumerate(atoms4)
    if a.parent().resname.strip() == "PRO" and a.name.strip() == "N")
  assert QMRegionGrower._is_chain_break_amine(
    gap_proline_n, adjacency4, atoms4) is True, (
    "a proline at a chain break was not recognised")

  # A real N-terminus is not a chain break: it carries its own hydrogens.
  _m2, atoms2, adjacency2 = _protonated_graph(_TRIPEPTIDE_PDB)
  first_n = next(i for i, a in enumerate(atoms2)
                 if a.parent().parent().resseq.strip() == "1"
                 and a.name.strip() == "N")
  assert QMRegionGrower._is_chain_break_amine(
    first_n, adjacency2, atoms2) is False, (
    "a real N-terminus read as a chain break")

  # The threshold from the other side.  An N-terminal proline is secondary
  # and holds exactly two hydrogens, so it is the case that separates
  # "fewer than two" from "two or fewer"; the `not _is_terminal_amine`
  # conjunct in _cut_reason hides a wrong answer here, which is why it is
  # asked of the predicate directly.
  _m3, atoms3, adjacency3 = _protonated_graph(_PROLINE_NTERM_PDB)
  proline_n = next(i for i, a in enumerate(atoms3)
                   if a.parent().parent().resseq.strip() == "1"
                   and a.name.strip() == "N")
  assert QMRegionGrower._is_chain_break_amine(
    proline_n, adjacency3, atoms3) is False, (
    "an N-terminal proline read as a chain break")


def exercise_terminal_amine_is_not_any_nitrogen():
  """The amine guard covers the chain terminus and nothing else.

  Two ways to get this wrong, both of which fired while writing it. A
  sidechain nitrogen carries only its own residue's carbons, exactly as a
  terminus does, so LYS NZ and TRP NE1 read as termini unless the test
  asks for the one bonded to CA. And a residue whose predecessor is simply
  absent from the file looks terminal too; a region carved around a metal
  is full of those, and what separates them is that an amide nitrogen
  keeps one hydrogen where a free amino group keeps more."""
  _model, atoms, adjacency = _protonated_graph(_1BQ8_FE_SPHERE_PDB)

  def nitrogen(resseq, name):
    return next(i_seq for i_seq, atom in enumerate(atoms)
                if atom.parent().parent().resseq.strip() == resseq
                and atom.name.strip() == name)

  # Sidechain nitrogens, one of them carrying three hydrogens.
  for resseq, name in (("7", "NZ"), ("37", "NE1")):
    assert QMRegionGrower._is_terminal_amine(
      nitrogen(resseq, name), adjacency, atoms) is False, (
      f"sidechain nitrogen {resseq}/{name} read as a chain terminus")

  # Backbone nitrogens whose preceding residue is outside this extract.
  # They keep the single amide hydrogen, so they are not free amines.
  for resseq in ("5", "37", "48"):
    assert QMRegionGrower._is_terminal_amine(
      nitrogen(resseq, "N"), adjacency, atoms) is False, (
      f"residue {resseq} is mid-chain, its predecessor merely absent")


def exercise_terminal_amine_kept():
  """Cutting CA-N away from the N-terminus would change the region's charge.

  The mirror of :func:`exercise_terminal_carboxylate_kept`: a terminal
  -NH3(+) replaced by a single hydrogen silently takes its positive charge
  with it, so that cut is refused.  A nitrogen mid-chain carries the
  preceding carbonyl carbon and is still cut."""
  seed = ["chain A and resseq 2 and name SG"]
  for radius in (3.0, 5.0):
    result = _run_endo_exo_params(
      _TRIPEPTIDE_PDB, selection=seed, buffer__radius=radius)[0]
    hierarchy = result["model"].get_hierarchy()
    caps = set(result["cap_iseqs"])
    for i_seq, atom in enumerate(hierarchy.atoms()):
      residue_group = atom.parent().parent()
      if residue_group.resseq.strip() != "1" or atom.name.strip() != "N":
        continue
      assert i_seq not in caps, (
        f"at radius {radius} the N-terminal amine was capped away")

  # The mid-chain nitrogen of residue 2 carries C of residue 1, so it is
  # not a terminus and the guard must not protect it.
  _model, atoms, adjacency = _protonated_graph(_TRIPEPTIDE_PDB)

  def nitrogen_of(resseq):
    return next(i_seq for i_seq, atom in enumerate(atoms)
                if atom.parent().parent().resseq.strip() == resseq
                and atom.name.strip() == "N")

  assert QMRegionGrower._is_terminal_amine(
    nitrogen_of("1"), adjacency, atoms) is True
  for resseq in ("2", "3"):
    assert QMRegionGrower._is_terminal_amine(
      nitrogen_of(resseq), adjacency, atoms) is False, (
      f"residue {resseq}'s nitrogen is mid-chain, not a terminus")


def exercise_terminal_carboxylate_kept():
  """Cutting CA-C away from a carboxylate would change the region's charge.

  A terminal -COO(-) replaced by a single hydrogen silently takes its
  negative charge with it, so that cut is refused and the carboxylate is
  retained whether or not the radius search reached it.  A chain break is
  not the same case: only a neutral carbonyl leaves there, so it is still
  cut and capped."""
  seed = ["chain A and resseq 2 and name SG"]
  reached = _run_endo_exo_params(
    _TRIPEPTIDE_PDB, selection=seed, buffer__radius=4.0)[0]
  not_reached = _run_endo_exo_params(
    _TRIPEPTIDE_PDB, selection=seed, buffer__radius=3.0)[0]

  for tag, result in (("radius 4.0", reached), ("radius 3.0", not_reached)):
    names = _atom_names_by_residue(result)
    assert "3" in names, f"{tag}: C-terminal residue absent"
    assert {"C", "O", "OXT"} <= names["3"], (
      f"{tag}: carboxylate not retained, got {sorted(names['3'])}")

  # Whether OXT sat inside the sphere must not decide whether the
  # carboxylate survives.
  assert _atom_names_by_residue(reached)["3"] == \
      _atom_names_by_residue(not_reached)["3"], (
    "C-terminal residue depends on the search radius")

  # Same peptide with OXT removed is a chain break, not a terminus: the
  # cut must still be made, leaving a capped C rather than a carbonyl
  # carbon with an unsatisfied valence.
  # Drop the line by atom name rather than by its whole text, which
  # carries a serial number that moves whenever the fixture changes.
  broken = "\n".join(line for line in _TRIPEPTIDE_PDB.split("\n")
                     if line[12:16].strip() != "OXT")
  assert "OXT" not in broken
  result = _run_endo_exo_params(
    broken, selection=seed, buffer__radius=3.0)[0]
  hierarchy = result["model"].get_hierarchy()
  atoms = list(hierarchy.atoms())
  caps = set(result["cap_iseqs"])
  for i, atom in enumerate(atoms):
    if atom.parent().parent().resseq.strip() != "3":
      continue
    if atom.name.strip() != "C":
      continue
    assert i in caps, (
      "at a chain break the carbonyl C must be capped, not retained "
      "with an unsatisfied valence")
    break
  else:
    raise AssertionError("residue 3 C absent from the region")


def exercise_scaffold_residue_keeps_no_functional_group():
  """A residue reached through its backbone keeps no functional group.

  Lys 7 of 1BQ8 comes within 5.88 A of the Fe through its backbone N,
  while its sidechain lies 7.9-11.3 A out, so BFS arrives from the
  mainchain and never from the tip.  Its configured cut is CD-CE, which
  BFS would only reach after walking the whole sidechain; the geometric
  heuristic takes the first sp3 C-C instead.  The ammonium must not
  survive at any radius, though atoms the sphere itself demanded do."""
  for radius in (5.0, 6.0, 7.0, 8.0, 9.0):
    result = _run_endo_exo_params(
      _1BQ8_FE_SPHERE_PDB, buffer__radius=radius)[0]
    hierarchy = result["model"].get_hierarchy()
    heavy = set()
    for atom in hierarchy.atoms():
      if atom.element_is_hydrogen():
        continue
      if atom.parent().parent().resseq.strip() == "7":
        heavy.add(atom.name.strip())
    if not heavy:
      continue  # outside the region at this radius
    assert {"N", "CA", "C", "O"} <= heavy, (
      f"radius {radius}: Lys 7 lost backbone atoms: {sorted(heavy)}")
    assert not ({"CD", "CE", "NZ"} & heavy), (
      f"radius {radius}: Lys 7 is reached through its backbone, so its "
      f"ammonium must not survive; got {sorted(heavy)}")


def exercise_ligand_keeps_its_coordinating_group():
  """A residue reached from its tip is cut at its configured bond.

  BFS enters the Cys ligands of 1BQ8 at SG, so it meets the CA-CB entry
  first and the thiolate survives."""
  result = _run_endo_exo_params(_1BQ8_FE_SPHERE_PDB)[0]
  hierarchy = result["model"].get_hierarchy()
  by_res = {}
  for atom in hierarchy.atoms():
    if atom.element_is_hydrogen():
      continue
    if atom.parent().resname.strip() != "CYS":
      continue
    rg = atom.parent().parent()
    by_res.setdefault(rg.resseq.strip(), set()).add(atom.name.strip())
  assert by_res, "no Cys ligands in the region"
  for resseq, names in sorted(by_res.items()):
    assert {"CB", "SG"} <= names, (
      f"coordinating Cys {resseq} lost its thiolate: {sorted(names)}")


def exercise_symmetry_images_truncate_alike():
  """Symmetry images of one residue are truncated the same way.

  On the 2C2U fixture the metal sits on a 3-fold, so several images of one
  residue enter the region; each is grown from its own seeds and they must
  agree."""
  result = _run_endo_exo_params(
    _2C2U_FE_SPHERE_PDB, buffer__radius=6.0)[0]
  hierarchy = result["model"].get_hierarchy()
  # One shape per residue COPY: union the atom groups of a residue group,
  # so a shared backbone plus an altloc'd sidechain counts once rather than
  # as two shapes.  Each symmetry image lands in its own chain.
  shapes = {}
  for chain in hierarchy.chains():
    for rg in chain.residue_groups():
      names = set()
      resname = ""
      for ag in rg.atom_groups():
        resname = resname or ag.resname.strip()
        names.update(a.name.strip() for a in ag.atoms()
                     if not a.element_is_hydrogen())
      if names:
        shapes.setdefault((resname, rg.resseq.strip()), set()).add(
          frozenset(names))
  assert shapes, "region is empty"
  for (resname, resseq), variants in sorted(shapes.items()):
    assert len(variants) == 1, (
      f"{resname} {resseq} truncated differently across symmetry images: "
      f"{[sorted(v) for v in variants]}")


def exercise_terminal_carboxylate_detection():
  """A carboxylate terminus is recognised from the graph, not from names.

  Two oxygens on the carbonyl carbon and no nitrogen in another residue.
  So a terminal oxygen named OT1/OT2 counts, and a stray OXT on a
  mid-chain residue does not, because that carbon still carries a peptide
  nitrogen."""
  from mmtbx.geometry_restraints.endo_exo.grow import QMRegionGrower
  from mmtbx.geometry_restraints.endo_exo.graph import AtomGraphBuilder

  def _verdicts(pdb_str):
    pdb_in = iotbx.pdb.input(source_info=None, lines=pdb_str.split("\n"))
    model = mmtbx.model.manager(model_input=pdb_in)
    model.process(make_restraints=True)
    grm = model.get_restraints_manager().geometry
    sites_cart = model.get_sites_cart()
    simple, asu = grm.get_all_bond_proxies(sites_cart=sites_cart)
    mappings = grm.pair_proxies(
      sites_cart=sites_cart).bond_proxies.asu_mappings()
    adjacency = AtomGraphBuilder().build_adjacency(simple, asu, mappings)
    hierarchy = model.get_hierarchy()
    atoms = hierarchy.atoms()
    out = {}
    for i, atom in enumerate(atoms):
      if atom.name.strip() != "C":
        continue
      resseq = atom.parent().parent().resseq.strip()
      out[resseq] = QMRegionGrower._is_terminal_carboxylate(i, adjacency, atoms)
    return out

  # Residue 3 is the real C-terminus; 1 and 2 are mid-chain.
  assert _verdicts(_TRIPEPTIDE_PDB) == {"1": False, "2": False, "3": True}

  # Legacy naming: the terminal oxygen is OT2, not OXT.
  legacy = _TRIPEPTIDE_PDB.replace(" OXT ALA A   3", " OT2 ALA A   3")
  assert _verdicts(legacy) == {"1": False, "2": False, "3": True}, (
    "a terminus named OT2 must still be recognised")

  # A stray OXT on a mid-chain residue is not a terminus: its carbonyl
  # carbon still carries the next residue's nitrogen.
  stray = _TRIPEPTIDE_PDB.replace(
    "END",
    "ATOM     18  OXT ALA A   1      12.700  11.500   9.100  1.00 10.00"
    "           O\nEND")
  assert _verdicts(stray) == {"1": False, "2": False, "3": True}, (
    "a mid-chain residue with a stray OXT must not read as a terminus")


def exercise_polymer_bonded_to_its_own_image_terminates():
  """A polymer with one residue per asymmetric unit stops growing.

  The chain leaves the asymmetric unit through a lattice translation, so the
  residue after this one is this one, under another symmetry image.  The
  CA-N cut is what bounds the walk, and its guard asks whether the next
  residue has been reached; comparing residue groups alone answers "there is
  no next residue" forever and BFS follows the chain out of the crystal.

  Driven on a hand-built adjacency because the point is the graph, and
  bounded by counting neighbour lookups, because the failure is a hang."""
  from mmtbx.geometry_restraints.endo_exo.grow import QMRegionGrower
  from mmtbx.geometry_restraints.endo_exo.util import _canon_op as canon

  class BoundedAdjacency(dict):
    """Raise rather than hang if growth does not settle."""
    def __init__(self, *args, budget=500):
      dict.__init__(self, *args)
      self.budget = budget

    def __getitem__(self, key):
      self.budget -= 1
      if self.budget < 0:
        raise RuntimeError("growth did not terminate")
      return dict.get(self, key, set())

  pdb_in = iotbx.pdb.input(source_info=None,
                           lines=_POLYMER_PDB.split("\n"))
  model = mmtbx.model.manager(model_input=pdb_in)
  atoms = model.get_hierarchy().atoms()
  iseq = {atom.name.strip(): atom.i_seq for atom in atoms}

  identity = canon(sgtbx.rt_mx())
  forward = canon(sgtbx.rt_mx("x+1,y,z"))
  backward = canon(sgtbx.rt_mx("x-1,y,z"))

  adjacency = BoundedAdjacency()
  def bond(first, second, there, back):
    adjacency.setdefault(iseq[first], set()).add((iseq[second], there))
    adjacency.setdefault(iseq[second], set()).add((iseq[first], back))

  bond("N", "CA", identity, identity)
  bond("CA", "C", identity, identity)
  bond("CA", "CB", identity, identity)
  bond("C", "O", identity, identity)
  # The peptide bond crosses a lattice translation, so the chain is endless.
  bond("C", "N", forward, backward)

  grower = QMRegionGrower(
    BondCutDetector(use_preferred_cuts=True, log=io.StringIO()),
    log=io.StringIO())
  visited, caps = grower.grow_by_depth(
    {(iseq["CA"], identity)}, adjacency, model)

  assert len(visited) < 20, f"region grew to {len(visited)} nodes"
  images = {op.as_xyz() for (_iseq, op) in visited}
  assert len(images) <= 2, f"growth ran through {len(images)} symmetry images"


def exercise_growth_converges_in_a_couple_of_rounds():
  """Growth and stranding repair alternate only a round or two.

  They terminate because the visited set only grows, a node is enqueued
  only while absent from it, and a demoted cap stays visited so it cannot
  be capped again -- each round therefore leaves strictly fewer caps to
  demote.  ``max_repair_rounds`` is a backstop for that reasoning being
  wrong, not a working limit, so a region approaching it means one of those
  properties has been broken and the loop is no longer known to end."""
  from mmtbx.geometry_restraints.endo_exo.grow import QMRegionGrower

  original = QMRegionGrower._grow_until_exhausted
  grow_by_depth = QMRegionGrower.grow_by_depth
  passes, current = [], [0]

  def counting_grow(self, *args, **kwargs):
    current[0] += 1
    return original(self, *args, **kwargs)

  # Counted per region rather than on a zero return from the repair: a region
  # that gives up through the max_repair_rounds backstop must still be
  # counted, or the worse the regression the less it records.
  def counting_region(self, *args, **kwargs):
    current[0] = 0
    try:
      return grow_by_depth(self, *args, **kwargs)
    finally:
      passes.append(current[0])

  QMRegionGrower._grow_until_exhausted = counting_grow
  QMRegionGrower.grow_by_depth = counting_region
  try:
    for pdb_str, overrides in (
        (_1BQ8_FE_SPHERE_PDB, {}),
        (_2C2U_FE_SPHERE_PDB, {})):
      for radius in (4.0, 6.0, 8.0):
        _run_endo_exo_params(pdb_str, buffer__radius=radius, **overrides)
  finally:
    QMRegionGrower._grow_until_exhausted = original
    QMRegionGrower.grow_by_depth = grow_by_depth

  assert len(passes) == 6, (
    f"expected one count per region, got {len(passes)}")
  assert max(passes) == 1, (
    f"a region took {max(passes)} growth passes; the loop is converging "
    f"more slowly than the reasoning behind max_repair_rounds allows")


def exercise_open_valences_are_reported():
  """An atom whose neighbour is missing from the model is named in the log.

  Capping covers bonds this code severs.  Where the input stops mid-chain,
  in a disordered loop or at the edge of an extract, there is no bond to
  sever and none to cap, so a nitrogen or carbonyl carbon reaches the
  region a bond short.  Nothing here can place the missing atom, so it is
  reported rather than repaired, and a real terminus -- which is not short
  of anything -- is not reported.

  1BQ8 runs from residue 45 to 48, so PRO 45's carbonyl has no following
  nitrogen once the sphere is wide enough to keep it."""
  def note_for(pdb_str, radius, selection=None):
    pdb_in = iotbx.pdb.input(source_info=None, lines=pdb_str.split("\n"))
    model = mmtbx.model.manager(model_input=pdb_in)
    dm = DataManager(["model"])
    dm.add_model("model", model)
    dm.set_default_model("model")
    master = libtbx.phil.parse(EndoexoProgram.master_phil_str)
    params = master.extract()
    params.write_files = False
    params.buffer.radius = radius
    if selection is not None:
      params.selection = [selection]
    log = io.StringIO()
    program = EndoexoProgram(dm, params, master_phil=master, logger=log)
    program.validate()
    program.run()
    return [line for line in log.getvalue().split("\n")
            if "a bond short" in line]

  # The whole set, not the number of lines: there is only ever one line, so
  # a count of it cannot fail, and a predicate that names extra atoms -- a
  # cap about to become a hydrogen, say -- would go unnoticed.
  def atoms_named(pdb_str, radius, selection=None):
    lines = note_for(pdb_str, radius, selection)
    if not lines:
      return set()
    assert len(lines) == 1, f"expected one note, got {lines}"
    return set(lines[0].rsplit(": ", 1)[1].rstrip(".").split(", "))

  # The extract jumps 45 -> 48, so PRO 45's carbonyl has no following
  # nitrogen.  GLU 48's nitrogen has no preceding one either, but it is far
  # enough out to be cut and capped, and a cap is not short of anything.
  # A nitrogen like it inside the radius floor is a seed, never cut-tested,
  # and is still reported: see ARG 92 below.
  assert atoms_named(_1BQ8_FE_SPHERE_PDB, 8.0) == {"PRO 45 C"}, (
    atoms_named(_1BQ8_FE_SPHERE_PDB, 8.0))

  # Narrower, and the same carbonyl is capped away instead: nothing to say.
  assert not note_for(_1BQ8_FE_SPHERE_PDB, 5.0), (
    "reported an open valence where the atom is not in the region")

  # Symmetry images of one atom are one problem, named once.  The extract
  # skips residue 91 and stops at 93, so two carbonyls have no following
  # nitrogen and one nitrogen has no preceding carbonyl.  ARG 92's nitrogen
  # is a seed, so it is retained rather than capped and stays reportable.
  assert atoms_named(_2C2U_FE_SPHERE_PDB, 7.0) == {
    "ARG 92 N", "ASP 93 C", "PHE 90 C"}, (
      atoms_named(_2C2U_FE_SPHERE_PDB, 7.0))

  # Only a peptide backbone has a neighbouring residue to be missing.  Tris
  # names an atom C and an atom N; the C is quaternary and complete, and
  # there is no residue next to either to be absent.
  assert not atoms_named(_TRS_LIGAND_PDB, 5.0, selection="resname ZN"), (
    atoms_named(_TRS_LIGAND_PDB, 5.0, selection="resname ZN"))

  # The tripeptide is whole, and its own termini are not short of anything.
  assert not note_for(_TRIPEPTIDE_PDB, 8.0,
                      selection="chain A and resseq 2 and name SG"), (
    "a complete molecule should report nothing")


def exercise_surviving_cap_has_only_its_anchor():
  """A cap that survives growth is bonded to its anchor and nothing else
  in the region.

  This is what keeps the repair machinery small.  Two branches it used to
  carry -- invalidating a demoted cap's neighbours, and promoting two
  adjacent caps to interior -- exist only for regions that break this, and
  were removed once measurement showed they never ran.  Their absence is
  safe exactly as long as this holds, so it is checked rather than implied.

  Checked on the grower's own state: once a region is materialised the caps
  are hydrogens 1.1 A from their anchors, which no longer shows which region
  atoms they were bonded to.

  Driven through the program rather than by seeding the grower here, because
  most of a region is symmetry images -- which hand-built seeds within a
  radius of the metal do not produce, leaving the composed-operator path the
  invariant most depends on unexamined."""
  from mmtbx.geometry_restraints.endo_exo.grow import QMRegionGrower
  from mmtbx.geometry_restraints.endo_exo.util import _canon_op as canon

  grow_by_depth = QMRegionGrower.grow_by_depth
  regions = []

  def capturing(self, seeds, adjacency, model, *args, **kwargs):
    visited, caps = grow_by_depth(
      self, seeds, adjacency, model, *args, **kwargs)
    regions.append((visited, dict(caps), adjacency, model))
    return visited, caps

  QMRegionGrower.grow_by_depth = capturing
  try:
    for pdb_str in (_1BQ8_FE_SPHERE_PDB, _2C2U_FE_SPHERE_PDB):
      for radius in (5.0, 7.0):
        _run_endo_exo_params(pdb_str, buffer__radius=radius)
  finally:
    QMRegionGrower.grow_by_depth = grow_by_depth

  assert len(regions) == 4, f"expected four regions, got {len(regions)}"
  checked, under_symmetry = 0, 0
  for (visited, caps, adjacency, model) in regions:
    atoms = model.get_hierarchy().atoms()
    assert set(caps) <= visited, "a cap is recorded outside the region"
    for cap, anchor in caps.items():
      cap_iseq, cap_op = cap
      if cap_op.as_xyz() != "x,y,z":
        under_symmetry += 1
      for (nb_iseq, edge_op) in adjacency[cap_iseq]:
        if atoms[nb_iseq].element_is_hydrogen():
          continue
        neighbour = (nb_iseq, canon(cap_op.multiply(edge_op)))
        if neighbour == anchor or neighbour not in visited:
          continue
        raise AssertionError(
          f"cap {atoms[cap_iseq].name.strip()} is also bonded to "
          f"{atoms[nb_iseq].name.strip()} in the region")
      checked += 1

  assert checked > 20, f"only {checked} caps examined"
  assert under_symmetry, "no cap under a symmetry image was examined"


def exercise_cut_not_made_into_the_region():
  """A cut is refused when the far atom is already bonded into the region.

  Capping there would leave the cap covalently bonded to a region atom
  that is not its anchor.  Refusing produces the same region while leaving
  less for the cap repair to undo afterwards."""
  from mmtbx.geometry_restraints.endo_exo.grow import QMRegionGrower

  demotions = {"n": 0}
  original = QMRegionGrower._demote_cap_candidate

  def counting(self, *args, **kwargs):
    demotions["n"] += 1
    return original(self, *args, **kwargs)

  QMRegionGrower._demote_cap_candidate = counting
  try:
    for pdb_str in (_1BQ8_FE_SPHERE_PDB, _2C2U_FE_SPHERE_PDB):
      for radius in (4.0, 5.0, 6.0, 7.0):
        for preferred in (True, False):
          result = _run_endo_exo_params(
            pdb_str, buffer__radius=radius,
            capping__preferred_cuts=preferred)[0]
          assert not region_invariant_violations(result), (
            f"radius={radius} preferred_cuts={preferred}")
  finally:
    QMRegionGrower._demote_cap_candidate = original

  # Repairs still happen, for the cases no cut-time rule can see: a bond
  # inside a ring, and a neighbour of the cap that BFS reaches only later.
  # The bound is what matters -- more than this means cuts are being made
  # that were knowably wrong when they were made.  Removing the refusal
  # entirely takes these same regions to 36, so the bound sits ON the
  # measured value: slack here is enough to hide a whole guard.  Note this
  # count only sees the far-atom half; the near-atom one is pinned by
  # exercise_cut_refused_when_it_would_strand_the_anchor, since removing it
  # changes no count here.
  assert demotions["n"] <= 17, (
    f"expected at most 17 cap repairs over these 16 regions, "
    f"got {demotions['n']}")


def run():
  exercise_submodel_shape()
  exercise_cys_coordination()
  exercise_residue_composition()
  exercise_residues_to_include()
  exercise_2c2u_symmetry_materialization()
  exercise_2c2u_fe_coordination_distances()
  exercise_2c2u_symmetry_truncation_consistency()
  exercise_2c2u_sym_image_provenance()
  exercise_special_position_water_dedup()
  exercise_selection_seed_terminates()
  exercise_region_invariants()
  # parameter coverage
  exercise_capping_disabled()
  exercise_special_position_partial_overlap()
  exercise_altloc_filter()
  exercise_element_filter()
  exercise_depth_and_skip_search()
  exercise_write_files_roundtrip()
  exercise_scaffold_residue_keeps_no_functional_group()
  exercise_ligand_keeps_its_coordinating_group()
  exercise_symmetry_images_truncate_alike()
  exercise_terminal_carboxylate_kept()
  exercise_terminal_amine_kept()
  exercise_terminal_amine_covers_proline()
  exercise_chain_break_amine_is_not_any_nitrogen()
  exercise_terminal_amine_is_not_any_nitrogen()
  exercise_amide_unit_is_never_split()
  exercise_radius_only_grows_the_region()
  exercise_region_does_not_depend_on_queue_order()
  exercise_terminal_carboxylate_detection()
  exercise_polymer_bonded_to_its_own_image_terminates()
  exercise_growth_converges_in_a_couple_of_rounds()
  exercise_open_valences_are_reported()
  exercise_surviving_cap_has_only_its_anchor()
  exercise_cut_not_made_into_the_region()
  # engine unit tests
  exercise_canon_op()
  exercise_hydrogen_capper()
  exercise_bond_cut_preferred()
  exercise_bond_cut_heuristic()
  exercise_proline_ring_is_never_cut()
  exercise_preferred_cut_rows()
  exercise_preferred_cuts_without_a_cut()
  exercise_cut_refused_when_it_would_strand_the_anchor()
  exercise_no_atom_is_left_unattached()
  exercise_bond_cut_backbone()
  exercise_consumed_preferred_cut()
  exercise_cut_depends_on_where_bfs_arrives()
  exercise_hbond_partners()
  exercise_hbond_partner_symmetry_composition()
  exercise_hull_waters()
  exercise_hull_waters_use_symmetry_images()
  exercise_hull_spans_symmetry_images_already_in_the_region()
  exercise_build_adjacency()
  print(format_cpu_times())
  print("OK")


if __name__ == "__main__":
  run()
