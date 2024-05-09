from __future__ import absolute_import, division, print_function
import sys
from six.moves import cStringIO as StringIO
from six.moves import range

pdbs = [
  # 0
  ["""
ATOM      1  N   ALA A   1      -1.521   3.898   0.713  1.00 20.00      A    N+1
ATOM      2  H   ALA A   1      -1.675   4.606   1.277  1.00 20.00      A    H
ATOM      3  H2  ALA A   1      -2.026   3.996  -0.045  1.00 20.00      A    H
ATOM      4  H3  ALA A   1      -0.632   3.875   0.487  1.00 20.00      A    H
ATOM      5  CA  ALA A   1      -1.879   2.678   1.370  1.00 20.00      A    C
ATOM      6  HA  ALA A   1      -1.349   2.598   2.189  1.00 20.00      A    H
ATOM      7  CB  ALA A   1      -3.346   2.715   1.747  1.00 20.00      A    C
ATOM      8  HB1 ALA A   1      -3.507   3.466   2.360  1.00 20.00      A    H
ATOM      9  HB2 ALA A   1      -3.592   1.873   2.188  1.00 20.00      A    H
ATOM     10  HB3 ALA A   1      -3.886   2.831   0.942  1.00 20.00      A    H
ATOM     11  C   ALA A   1      -1.551   1.486   0.458  1.00 20.00      A    C
ATOM     12  O   ALA A   1      -2.007   1.439  -0.662  1.00 20.00      A    O
ATOM     13  N   ALA A   2      -0.744   0.404   0.952  1.00 20.00      A    N
ATOM     14  H   ALA A   2      -0.359   0.462   1.781  1.00 20.00      A    H
ATOM     15  CA  ALA A   2      -0.493  -0.749   0.119  1.00 20.00      A    C
ATOM     16  HA  ALA A   2      -0.668  -0.439  -0.791  1.00 20.00      A    H
ATOM     17  CB  ALA A   2      -1.476  -1.885   0.331  1.00 20.00      A    C
ATOM     18  HB1 ALA A   2      -1.345  -2.266   1.221  1.00 20.00      A    H
ATOM     19  HB2 ALA A   2      -1.326  -2.579  -0.349  1.00 20.00      A    H
ATOM     20  HB3 ALA A   2      -2.395  -1.542   0.252  1.00 20.00      A    H
ATOM     21  C   ALA A   2       0.972  -1.193   0.100  1.00 20.00      A    C
ATOM     22  O   ALA A   2       1.603  -1.231   1.129  1.00 20.00      A    O
ATOM     23  N   ALA A   3       1.569  -1.642  -1.120  1.00 20.00      A    N
ATOM     24  H   ALA A   3       1.050  -1.766  -1.865  1.00 20.00      A    H
ATOM     25  CA  ALA A   3       2.971  -1.947  -1.163  1.00 20.00      A    C
ATOM     26  HA  ALA A   3       3.432  -1.386  -0.509  1.00 20.00      A    H
ATOM     27  CB  ALA A   3       3.520  -1.630  -2.540  1.00 20.00      A    C
ATOM     28  HB1 ALA A   3       4.474  -1.845  -2.568  1.00 20.00      A    H
ATOM     29  HB2 ALA A   3       3.393  -0.674  -2.732  1.00 20.00      A    H
ATOM     30  HB3 ALA A   3       3.045  -2.163  -3.211  1.00 20.00      A    H
ATOM     31  C   ALA A   3       3.187  -3.395  -0.790  1.00 20.00      A    C
ATOM     32  O   ALA A   3       2.646  -4.306  -1.466  1.00 20.00      A    O
ATOM     33  OXT ALA A   3       3.917  -3.690   0.195  1.00 20.00      A    O-1
""", """
ATOM      1  N   ALA A   1      -1.521   3.898   0.713  1.00 20.00      A    N+1
ATOM      2  CA  ALA A   1      -1.879   2.678   1.370  1.00 20.00      A    C
ATOM      3  C   ALA A   1      -1.551   1.486   0.458  1.00 20.00      A    C
ATOM      4  O   ALA A   1      -2.007   1.439  -0.662  1.00 20.00      A    O
ATOM      5  CB  ALA A   1      -3.346   2.715   1.747  1.00 20.00      A    C
ATOM      6  H   ALA A   1      -1.675   4.606   1.277  1.00 20.00      A    H
ATOM      7  H2  ALA A   1      -2.026   3.996  -0.045  1.00 20.00      A    H
ATOM      8  H3  ALA A   1      -0.632   3.875   0.487  1.00 20.00      A    H
ATOM      9  HA  ALA A   1      -1.349   2.598   2.189  1.00 20.00      A    H
ATOM     10  HB1 ALA A   1      -3.507   3.466   2.360  1.00 20.00      A    H
ATOM     11  HB2 ALA A   1      -3.592   1.873   2.188  1.00 20.00      A    H
ATOM     12  HB3 ALA A   1      -3.886   2.831   0.942  1.00 20.00      A    H
ATOM     13  N   ALA A   2      -0.744   0.404   0.952  1.00 20.00      A    N
ATOM     14  CA  ALA A   2      -0.493  -0.749   0.119  1.00 20.00      A    C
ATOM     15  C   ALA A   2       0.972  -1.193   0.100  1.00 20.00      A    C
ATOM     16  O   ALA A   2       1.603  -1.231   1.129  1.00 20.00      A    O
ATOM     17  CB  ALA A   2      -1.476  -1.885   0.331  1.00 20.00      A    C
ATOM     18  H   ALA A   2      -0.359   0.462   1.781  1.00 20.00      A    H
ATOM     19  HA  ALA A   2      -0.668  -0.439  -0.791  1.00 20.00      A    H
ATOM     20  HB1 ALA A   2      -1.345  -2.266   1.221  1.00 20.00      A    H
ATOM     21  HB2 ALA A   2      -1.326  -2.579  -0.349  1.00 20.00      A    H
ATOM     22  HB3 ALA A   2      -2.395  -1.542   0.252  1.00 20.00      A    H
ATOM     23  N   ALA A   3       1.569  -1.642  -1.120  1.00 20.00      A    N
ATOM     24  CA  ALA A   3       2.971  -1.947  -1.163  1.00 20.00      A    C
ATOM     25  C   ALA A   3       3.187  -3.395  -0.790  1.00 20.00      A    C
ATOM     26  O   ALA A   3       2.646  -4.306  -1.466  1.00 20.00      A    O
ATOM     27  CB  ALA A   3       3.520  -1.630  -2.540  1.00 20.00      A    C
ATOM     28  OXT ALA A   3       3.917  -3.690   0.195  1.00 20.00      A    O-1
ATOM     29  H   ALA A   3       1.050  -1.766  -1.865  1.00 20.00      A    H
ATOM     30  HA  ALA A   3       3.432  -1.386  -0.509  1.00 20.00      A    H
ATOM     31  HB1 ALA A   3       4.474  -1.845  -2.568  1.00 20.00      A    H
ATOM     32  HB2 ALA A   3       3.393  -0.674  -2.732  1.00 20.00      A    H
ATOM     33  HB3 ALA A   3       3.045  -2.163  -3.211  1.00 20.00      A    H
"""],
# no H
  ["""
ATOM      1  N   ALA A   1      -1.521   3.898   0.713  1.00 20.00      A    N+1
ATOM      5  CA  ALA A   1      -1.879   2.678   1.370  1.00 20.00      A    C
ATOM      7  CB  ALA A   1      -3.346   2.715   1.747  1.00 20.00      A    C
ATOM     11  C   ALA A   1      -1.551   1.486   0.458  1.00 20.00      A    C
ATOM     12  O   ALA A   1      -2.007   1.439  -0.662  1.00 20.00      A    O
ATOM     13  N   ALA A   2      -0.744   0.404   0.952  1.00 20.00      A    N
ATOM     15  CA  ALA A   2      -0.493  -0.749   0.119  1.00 20.00      A    C
ATOM     17  CB  ALA A   2      -1.476  -1.885   0.331  1.00 20.00      A    C
ATOM     21  C   ALA A   2       0.972  -1.193   0.100  1.00 20.00      A    C
ATOM     22  O   ALA A   2       1.603  -1.231   1.129  1.00 20.00      A    O
ATOM     23  N   ALA A   3       1.569  -1.642  -1.120  1.00 20.00      A    N
ATOM     25  CA  ALA A   3       2.971  -1.947  -1.163  1.00 20.00      A    C
ATOM     27  CB  ALA A   3       3.520  -1.630  -2.540  1.00 20.00      A    C
ATOM     31  C   ALA A   3       3.187  -3.395  -0.790  1.00 20.00      A    C
ATOM     32  O   ALA A   3       2.646  -4.306  -1.466  1.00 20.00      A    O
ATOM     33  OXT ALA A   3       3.917  -3.690   0.195  1.00 20.00      A    O-1
""", """
ATOM      1  N   ALA A   1      -1.521   3.898   0.713  1.00 20.00      A    N+1
ATOM      2  CA  ALA A   1      -1.879   2.678   1.370  1.00 20.00      A    C
ATOM      3  C   ALA A   1      -1.551   1.486   0.458  1.00 20.00      A    C
ATOM      4  O   ALA A   1      -2.007   1.439  -0.662  1.00 20.00      A    O
ATOM      5  CB  ALA A   1      -3.346   2.715   1.747  1.00 20.00      A    C
ATOM      6  N   ALA A   2      -0.744   0.404   0.952  1.00 20.00      A    N
ATOM      7  CA  ALA A   2      -0.493  -0.749   0.119  1.00 20.00      A    C
ATOM      8  C   ALA A   2       0.972  -1.193   0.100  1.00 20.00      A    C
ATOM      9  O   ALA A   2       1.603  -1.231   1.129  1.00 20.00      A    O
ATOM     10  CB  ALA A   2      -1.476  -1.885   0.331  1.00 20.00      A    C
ATOM     11  N   ALA A   3       1.569  -1.642  -1.120  1.00 20.00      A    N
ATOM     12  CA  ALA A   3       2.971  -1.947  -1.163  1.00 20.00      A    C
ATOM     13  C   ALA A   3       3.187  -3.395  -0.790  1.00 20.00      A    C
ATOM     14  O   ALA A   3       2.646  -4.306  -1.466  1.00 20.00      A    O
ATOM     15  CB  ALA A   3       3.520  -1.630  -2.540  1.00 20.00      A    C
ATOM     16  OXT ALA A   3       3.917  -3.690   0.195  1.00 20.00      A    O-1
"""],
# just CB
  ["""
ATOM      1  N   ALA A   1      -1.521   3.898   0.713  1.00 20.00      A    N+1
ATOM      5  CA  ALA A   1      -1.879   2.678   1.370  1.00 20.00      A    C
ATOM      7  CB  ALA A   1      -3.346   2.715   1.747  1.00 20.00      A    C
ATOM     11  C   ALA A   1      -1.551   1.486   0.458  1.00 20.00      A    C
ATOM     12  O   ALA A   1      -2.007   1.439  -0.662  1.00 20.00      A    O
ATOM     13  N   ALA A   2      -0.744   0.404   0.952  1.00 20.00      A    N
ATOM     15  CA  ALA A   2      -0.493  -0.749   0.119  1.00 20.00      A    C
ATOM     17  CB AALA A   2      -1.476  -1.885   0.331  1.00 20.00      A    C
ATOM     17  CB BALA A   2      -1.476  -1.885   0.331  1.00 20.00      A    C
ATOM     21  C   ALA A   2       0.972  -1.193   0.100  1.00 20.00      A    C
ATOM     22  O   ALA A   2       1.603  -1.231   1.129  1.00 20.00      A    O
ATOM     23  N   ALA A   3       1.569  -1.642  -1.120  1.00 20.00      A    N
ATOM     25  CA  ALA A   3       2.971  -1.947  -1.163  1.00 20.00      A    C
ATOM     27  CB  ALA A   3       3.520  -1.630  -2.540  1.00 20.00      A    C
ATOM     31  C   ALA A   3       3.187  -3.395  -0.790  1.00 20.00      A    C
ATOM     32  O   ALA A   3       2.646  -4.306  -1.466  1.00 20.00      A    O
ATOM     33  OXT ALA A   3       3.917  -3.690   0.195  1.00 20.00      A    O-1
""", """
ATOM      1  N   ALA A   1      -1.521   3.898   0.713  1.00 20.00      A    N+1
ATOM      2  CA  ALA A   1      -1.879   2.678   1.370  1.00 20.00      A    C
ATOM      3  C   ALA A   1      -1.551   1.486   0.458  1.00 20.00      A    C
ATOM      4  O   ALA A   1      -2.007   1.439  -0.662  1.00 20.00      A    O
ATOM      5  CB  ALA A   1      -3.346   2.715   1.747  1.00 20.00      A    C
ATOM      6  N   ALA A   2      -0.744   0.404   0.952  1.00 20.00      A    N
ATOM      7  CA  ALA A   2      -0.493  -0.749   0.119  1.00 20.00      A    C
ATOM      8  C   ALA A   2       0.972  -1.193   0.100  1.00 20.00      A    C
ATOM      9  O   ALA A   2       1.603  -1.231   1.129  1.00 20.00      A    O
ATOM     10  CB AALA A   2      -1.476  -1.885   0.331  0.50 20.00      A    C
ATOM     11  CB BALA A   2      -1.476  -1.885   0.331  0.50 20.00      A    C
ATOM     12  N   ALA A   3       1.569  -1.642  -1.120  1.00 20.00      A    N
ATOM     13  CA  ALA A   3       2.971  -1.947  -1.163  1.00 20.00      A    C
ATOM     14  C   ALA A   3       3.187  -3.395  -0.790  1.00 20.00      A    C
ATOM     15  O   ALA A   3       2.646  -4.306  -1.466  1.00 20.00      A    O
ATOM     16  CB  ALA A   3       3.520  -1.630  -2.540  1.00 20.00      A    C
ATOM     17  OXT ALA A   3       3.917  -3.690   0.195  1.00 20.00      A    O-1
"""],
# just CB diff
  ["""
ATOM      1  N   ALA A   1      -1.521   3.898   0.713  1.00 20.00      A    N+1
ATOM      5  CA  ALA A   1      -1.879   2.678   1.370  1.00 20.00      A    C
ATOM      7  CB  ALA A   1      -3.346   2.715   1.747  1.00 20.00      A    C
ATOM     11  C   ALA A   1      -1.551   1.486   0.458  1.00 20.00      A    C
ATOM     12  O   ALA A   1      -2.007   1.439  -0.662  1.00 20.00      A    O
ATOM     13  N   ALA A   2      -0.744   0.404   0.952  1.00 20.00      A    N
ATOM     15  CA  ALA A   2      -0.493  -0.749   0.119  1.00 20.00      A    C
ATOM     17  CB AALA A   2      -1.976  -1.885   0.331  1.00 20.00      A    C
ATOM     17  CB BALA A   2      -1.476  -1.885   0.331  1.00 20.00      A    C
ATOM     21  C   ALA A   2       0.972  -1.193   0.100  1.00 20.00      A    C
ATOM     22  O   ALA A   2       1.603  -1.231   1.129  1.00 20.00      A    O
ATOM     23  N   ALA A   3       1.569  -1.642  -1.120  1.00 20.00      A    N
ATOM     25  CA  ALA A   3       2.971  -1.947  -1.163  1.00 20.00      A    C
ATOM     27  CB  ALA A   3       3.520  -1.630  -2.540  1.00 20.00      A    C
ATOM     31  C   ALA A   3       3.187  -3.395  -0.790  1.00 20.00      A    C
ATOM     32  O   ALA A   3       2.646  -4.306  -1.466  1.00 20.00      A    O
ATOM     33  OXT ALA A   3       3.917  -3.690   0.195  1.00 20.00      A    O-1
""", """
ATOM      1  N   ALA A   1      -1.521   3.898   0.713  1.00 20.00      A    N+1
ATOM      2  CA  ALA A   1      -1.879   2.678   1.370  1.00 20.00      A    C
ATOM      3  CB  ALA A   1      -3.346   2.715   1.747  1.00 20.00      A    C
ATOM      4  C  AALA A   1      -1.551   1.486   0.458  0.50 20.00      A    C
ATOM      5  O  AALA A   1      -2.007   1.439  -0.662  0.50 20.00      A    O
ATOM      6  C  BALA A   1      -1.551   1.486   0.458  0.50 20.00      A    C
ATOM      7  O  BALA A   1      -2.007   1.439  -0.662  0.50 20.00      A    O
ATOM      8  CB AALA A   2      -1.976  -1.885   0.331  0.50 20.00      A    C
ATOM      9  N  AALA A   2      -0.744   0.404   0.952  0.50 20.00      A    N
ATOM     10  CA AALA A   2      -0.493  -0.749   0.119  0.50 20.00      A    C
ATOM     11  C  AALA A   2       0.972  -1.193   0.100  0.50 20.00      A    C
ATOM     12  O  AALA A   2       1.603  -1.231   1.129  0.50 20.00      A    O
ATOM     13  CB BALA A   2      -1.476  -1.885   0.331  0.50 20.00      A    C
ATOM     14  N  BALA A   2      -0.744   0.404   0.952  0.50 20.00      A    N
ATOM     15  CA BALA A   2      -0.493  -0.749   0.119  0.50 20.00      A    C
ATOM     16  C  BALA A   2       0.972  -1.193   0.100  0.50 20.00      A    C
ATOM     17  O  BALA A   2       1.603  -1.231   1.129  0.50 20.00      A    O
ATOM     18  CA  ALA A   3       2.971  -1.947  -1.163  1.00 20.00      A    C
ATOM     19  C   ALA A   3       3.187  -3.395  -0.790  1.00 20.00      A    C
ATOM     20  O   ALA A   3       2.646  -4.306  -1.466  1.00 20.00      A    O
ATOM     21  CB  ALA A   3       3.520  -1.630  -2.540  1.00 20.00      A    C
ATOM     22  OXT ALA A   3       3.917  -3.690   0.195  1.00 20.00      A    O-1
ATOM     23  N  AALA A   3       1.569  -1.642  -1.120  0.50 20.00      A    N
ATOM     24  N  BALA A   3       1.569  -1.642  -1.120  0.50 20.00      A    N
"""],
# all 2
  ["""
ATOM      1  N   ALA A   1      -1.521   3.898   0.713  1.00 20.00      A    N+1
ATOM      5  CA  ALA A   1      -1.879   2.678   1.370  1.00 20.00      A    C
ATOM      7  CB  ALA A   1      -3.346   2.715   1.747  1.00 20.00      A    C
ATOM     11  C   ALA A   1      -1.551   1.486   0.458  1.00 20.00      A    C
ATOM     12  O   ALA A   1      -2.007   1.439  -0.662  1.00 20.00      A    O
ATOM     13  N  AALA A   2      -0.744   0.404   0.952  1.00 20.00      A    N
ATOM     15  CA AALA A   2      -0.493  -0.749   0.119  1.00 20.00      A    C
ATOM     17  CB AALA A   2      -1.976  -1.885   0.331  1.00 20.00      A    C
ATOM     21  C  AALA A   2       0.972  -1.193   0.100  1.00 20.00      A    C
ATOM     22  O  AALA A   2       1.603  -1.231   1.129  1.00 20.00      A    O
ATOM     13  N  BALA A   2      -0.744   0.404   0.952  1.00 20.00      A    N
ATOM     15  CA BALA A   2      -0.493  -0.749   0.119  1.00 20.00      A    C
ATOM     17  CB BALA A   2      -1.976  -1.885   0.331  1.00 20.00      A    C
ATOM     21  C  BALA A   2       0.972  -1.193   0.100  1.00 20.00      A    C
ATOM     22  O  BALA A   2       1.603  -1.231   1.129  1.00 20.00      A    O
ATOM     23  N   ALA A   3       1.569  -1.642  -1.120  1.00 20.00      A    N
ATOM     25  CA  ALA A   3       2.971  -1.947  -1.163  1.00 20.00      A    C
ATOM     27  CB  ALA A   3       3.520  -1.630  -2.540  1.00 20.00      A    C
ATOM     31  C   ALA A   3       3.187  -3.395  -0.790  1.00 20.00      A    C
ATOM     32  O   ALA A   3       2.646  -4.306  -1.466  1.00 20.00      A    O
ATOM     33  OXT ALA A   3       3.917  -3.690   0.195  1.00 20.00      A    O-1
""", """
ATOM      1  N   ALA A   1      -1.521   3.898   0.713  1.00 20.00      A    N+1
ATOM      2  CA  ALA A   1      -1.879   2.678   1.370  1.00 20.00      A    C
ATOM      3  CB  ALA A   1      -3.346   2.715   1.747  1.00 20.00      A    C
ATOM      4  C  AALA A   1      -1.551   1.486   0.458  0.50 20.00      A    C
ATOM      5  O  AALA A   1      -2.007   1.439  -0.662  0.50 20.00      A    O
ATOM      6  C  BALA A   1      -1.551   1.486   0.458  0.50 20.00      A    C
ATOM      7  O  BALA A   1      -2.007   1.439  -0.662  0.50 20.00      A    O
ATOM      8  N  AALA A   2      -0.744   0.404   0.952  0.50 20.00      A    N
ATOM      9  CA AALA A   2      -0.493  -0.749   0.119  0.50 20.00      A    C
ATOM     10  C  AALA A   2       0.972  -1.193   0.100  0.50 20.00      A    C
ATOM     11  O  AALA A   2       1.603  -1.231   1.129  0.50 20.00      A    O
ATOM     12  CB AALA A   2      -1.976  -1.885   0.331  0.50 20.00      A    C
ATOM     13  N  BALA A   2      -0.744   0.404   0.952  0.50 20.00      A    N
ATOM     14  CA BALA A   2      -0.493  -0.749   0.119  0.50 20.00      A    C
ATOM     15  C  BALA A   2       0.972  -1.193   0.100  0.50 20.00      A    C
ATOM     16  O  BALA A   2       1.603  -1.231   1.129  0.50 20.00      A    O
ATOM     17  CB BALA A   2      -1.976  -1.885   0.331  0.50 20.00      A    C
ATOM     18  CA  ALA A   3       2.971  -1.947  -1.163  1.00 20.00      A    C
ATOM     19  C   ALA A   3       3.187  -3.395  -0.790  1.00 20.00      A    C
ATOM     20  O   ALA A   3       2.646  -4.306  -1.466  1.00 20.00      A    O
ATOM     21  CB  ALA A   3       3.520  -1.630  -2.540  1.00 20.00      A    C
ATOM     22  OXT ALA A   3       3.917  -3.690   0.195  1.00 20.00      A    O-1
ATOM     23  N  AALA A   3       1.569  -1.642  -1.120  0.50 20.00      A    N
ATOM     24  N  BALA A   3       1.569  -1.642  -1.120  0.50 20.00      A    N
"""],
# all 2 interlocked
  ["""
ATOM      1  N   ALA A   1      -1.521   3.898   0.713  1.00 20.00      A    N+1
ATOM      5  CA  ALA A   1      -1.879   2.678   1.370  1.00 20.00      A    C
ATOM      7  CB  ALA A   1      -3.346   2.715   1.747  1.00 20.00      A    C
ATOM     11  C   ALA A   1      -1.551   1.486   0.458  1.00 20.00      A    C
ATOM     12  O   ALA A   1      -2.007   1.439  -0.662  1.00 20.00      A    O
ATOM     13  N  AALA A   2      -0.744   0.404   0.952  1.00 20.00      A    N
ATOM     13  N  BALA A   2      -0.744   0.404   0.952  1.00 20.00      A    N
ATOM     15  CA AALA A   2      -0.493  -0.749   0.119  1.00 20.00      A    C
ATOM     15  CA BALA A   2      -0.493  -0.749   0.119  1.00 20.00      A    C
ATOM     17  CB AALA A   2      -1.976  -1.885   0.331  1.00 20.00      A    C
ATOM     17  CB BALA A   2      -1.976  -1.885   0.331  1.00 20.00      A    C
ATOM     21  C  AALA A   2       0.972  -1.193   0.100  1.00 20.00      A    C
ATOM     21  C  BALA A   2       0.972  -1.193   0.100  1.00 20.00      A    C
ATOM     22  O  AALA A   2       1.603  -1.231   1.129  1.00 20.00      A    O
ATOM     22  O  BALA A   2       1.603  -1.231   1.129  1.00 20.00      A    O
ATOM     23  N   ALA A   3       1.569  -1.642  -1.120  1.00 20.00      A    N
ATOM     25  CA  ALA A   3       2.971  -1.947  -1.163  1.00 20.00      A    C
ATOM     27  CB  ALA A   3       3.520  -1.630  -2.540  1.00 20.00      A    C
ATOM     31  C   ALA A   3       3.187  -3.395  -0.790  1.00 20.00      A    C
ATOM     32  O   ALA A   3       2.646  -4.306  -1.466  1.00 20.00      A    O
ATOM     33  OXT ALA A   3       3.917  -3.690   0.195  1.00 20.00      A    O-1
""", """
ATOM      1  N   ALA A   1      -1.521   3.898   0.713  1.00 20.00      A    N+1
ATOM      2  CA  ALA A   1      -1.879   2.678   1.370  1.00 20.00      A    C
ATOM      3  CB  ALA A   1      -3.346   2.715   1.747  1.00 20.00      A    C
ATOM      4  C  AALA A   1      -1.551   1.486   0.458  0.50 20.00      A    C
ATOM      5  O  AALA A   1      -2.007   1.439  -0.662  0.50 20.00      A    O
ATOM      6  C  BALA A   1      -1.551   1.486   0.458  0.50 20.00      A    C
ATOM      7  O  BALA A   1      -2.007   1.439  -0.662  0.50 20.00      A    O
ATOM      8  N  AALA A   2      -0.744   0.404   0.952  0.50 20.00      A    N
ATOM      9  CA AALA A   2      -0.493  -0.749   0.119  0.50 20.00      A    C
ATOM     10  C  AALA A   2       0.972  -1.193   0.100  0.50 20.00      A    C
ATOM     11  O  AALA A   2       1.603  -1.231   1.129  0.50 20.00      A    O
ATOM     12  CB AALA A   2      -1.976  -1.885   0.331  0.50 20.00      A    C
ATOM     13  N  BALA A   2      -0.744   0.404   0.952  0.50 20.00      A    N
ATOM     14  CA BALA A   2      -0.493  -0.749   0.119  0.50 20.00      A    C
ATOM     15  C  BALA A   2       0.972  -1.193   0.100  0.50 20.00      A    C
ATOM     16  O  BALA A   2       1.603  -1.231   1.129  0.50 20.00      A    O
ATOM     17  CB BALA A   2      -1.976  -1.885   0.331  0.50 20.00      A    C
ATOM     18  CA  ALA A   3       2.971  -1.947  -1.163  1.00 20.00      A    C
ATOM     19  C   ALA A   3       3.187  -3.395  -0.790  1.00 20.00      A    C
ATOM     20  O   ALA A   3       2.646  -4.306  -1.466  1.00 20.00      A    O
ATOM     21  CB  ALA A   3       3.520  -1.630  -2.540  1.00 20.00      A    C
ATOM     22  OXT ALA A   3       3.917  -3.690   0.195  1.00 20.00      A    O-1
ATOM     23  N  AALA A   3       1.569  -1.642  -1.120  0.50 20.00      A    N
ATOM     24  N  BALA A   3       1.569  -1.642  -1.120  0.50 20.00      A    N
"""],
# correcting partials
  ["""
ATOM    757  N   GLY A 103      21.999  10.712   0.961  1.00  4.33           N
ATOM    758  CA  GLY A 103      21.040  11.791   1.055  1.00  4.23           C
ATOM    759  C   GLY A 103      19.704  11.419   0.439  1.00  3.85           C
ATOM    760  O   GLY A 103      19.580  10.539  -0.422  1.00  4.32           O
ATOM    761  N   GLU A 104      18.680  12.148   0.879  1.00  3.87           N
ATOM    762  CA  GLU A 104      17.299  11.909   0.484  1.00  3.91           C
ATOM    763  O   GLU A 104      17.826  13.338  -1.410  1.00  4.76           O
ATOM    764  CB  GLU A 104      16.372  12.308   1.645  1.00  4.08           C
ATOM    765  CG  GLU A 104      16.686  11.548   2.930  1.00  4.59           C
ATOM    766  CD  GLU A 104      17.829  12.093   3.786  1.00  5.00           C
ATOM    767  OE1 GLU A 104      18.152  11.408   4.823  1.00  6.45           O
ATOM    768  OE2 GLU A 104      18.390  13.159   3.470  1.00  5.26           O
ATOM    769  C  AGLU A 104      16.973  12.646  -0.808  0.50  3.96           C
ATOM    770  C  BGLU A 104      16.975  12.684  -0.800  0.50  3.89           C
ATOM    771  N  ALEU A 105      15.717  12.499  -1.259  0.50  4.06           N
ATOM    772  CA ALEU A 105      15.284  13.059  -2.541  0.50  4.88           C
ATOM    773  C  ALEU A 105      15.637  14.444  -2.871  0.50  4.43           C
ATOM    774  O  ALEU A 105      15.735  14.769  -4.050  0.50  6.07           O
ATOM    775  CB ALEU A 105      13.730  13.346  -2.655  0.50  3.54           C
ATOM    776  CG ALEU A 105      13.190  11.968  -3.122  0.50  3.81           C
ATOM    777  CD1ALEU A 105      11.685  12.054  -3.263  0.50  4.54           C
ATOM    778  CD2ALEU A 105      13.804  11.453  -4.425  0.50  4.05           C
ATOM    779  N  BLEU A 105      15.711  12.640  -1.220  0.50  3.57           N
ATOM    780  CA BLEU A 105      15.237  13.396  -2.364  0.50  4.37           C
ATOM    781  C  BLEU A 105      15.537  14.890  -2.046  0.50  3.47           C
ATOM    782  O  BLEU A 105      15.411  15.372  -0.929  0.50  3.89           O
ATOM    783  CB BLEU A 105      13.764  12.747  -2.519  0.50  4.65           C
ATOM    784  CG BLEU A 105      13.014  12.862  -3.842  0.50  6.20           C
ATOM    785  CD1BLEU A 105      13.508  11.889  -4.865  0.50  7.64           C
ATOM    786  CD2BLEU A 105      11.548  12.657  -3.620  0.50  5.53           C
ATOM    787  N  AGLY A 106      15.734  15.299  -1.857  0.50  4.10           N
ATOM    788  CA AGLY A 106      16.043  16.715  -2.048  0.50  5.76           C
ATOM    789  C  AGLY A 106      17.449  17.119  -1.650  0.50  4.81           C
ATOM    790  O  AGLY A 106      17.778  18.310  -1.696  0.50  5.91           O
ATOM    791  N  BGLY A 106      15.975  15.616  -3.081  0.50  3.71           N
ATOM    792  CA BGLY A 106      16.294  17.031  -2.939  0.50  4.59           C
ATOM    793  C  BGLY A 106      17.649  17.368  -2.295  0.50  4.30           C
ATOM    794  O  BGLY A 106      18.056  18.524  -2.339  0.50  4.74           O
ATOM    795  CA  SER A 107      19.612  16.580  -0.930  1.00  4.51           C
ATOM    796  C   SER A 107      20.554  17.027  -2.048  1.00  4.08           C
ATOM    797  O   SER A 107      20.441  16.592  -3.191  1.00  4.85           O
ATOM    798  CB  SER A 107      20.205  15.355  -0.212  1.00  4.55           C
ATOM    799  OG  SER A 107      20.403  14.244  -1.072  1.00  4.35           O
ATOM    800  N  ASER A 107      18.312  16.190  -1.350  0.50  5.31           N
ATOM    801  N  BSER A 107      18.278  16.387  -1.593  0.50  5.42           N
ATOM    802  N   PRO A 108      21.538  17.876  -1.717  1.00  4.04           N
ATOM    803  CA  PRO A 108      22.472  18.325  -2.750  1.00  4.37           C
ATOM    804  C   PRO A 108      23.241  17.187  -3.437  1.00  4.03           C
ATOM    805  O   PRO A 108      23.746  16.255  -2.804  1.00  4.51           O
ATOM    806  CB  PRO A 108      23.413  19.264  -1.990  1.00  4.85           C
ATOM    807  CG  PRO A 108      22.520  19.852  -0.907  1.00  5.00           C
ATOM    808  CD  PRO A 108      21.676  18.668  -0.477  1.00  4.48           C
""", """
ATOM      1  N   GLY A 103      21.999  10.712   0.961  1.00  4.33           N
ATOM      2  CA  GLY A 103      21.040  11.791   1.055  1.00  4.23           C
ATOM      3  C   GLY A 103      19.704  11.419   0.439  1.00  3.85           C
ATOM      4  O   GLY A 103      19.580  10.539  -0.422  1.00  4.32           O
ATOM      5  N   GLU A 104      18.680  12.148   0.879  1.00  3.87           N
ATOM      6  CA  GLU A 104      17.299  11.909   0.484  1.00  3.91           C
ATOM      7  CB  GLU A 104      16.372  12.308   1.645  1.00  4.08           C
ATOM      8  CG  GLU A 104      16.686  11.548   2.930  1.00  4.59           C
ATOM      9  CD  GLU A 104      17.829  12.093   3.786  1.00  5.00           C
ATOM     10  OE1 GLU A 104      18.152  11.408   4.823  1.00  6.45           O
ATOM     11  OE2 GLU A 104      18.390  13.159   3.470  1.00  5.26           O
ATOM     12  C  AGLU A 104      16.973  12.646  -0.808  0.50  3.96           C
ATOM     13  O  AGLU A 104      17.826  13.338  -1.410  0.50  4.76           O
ATOM     14  C  BGLU A 104      16.975  12.684  -0.800  0.50  3.89           C
ATOM     15  O  BGLU A 104      17.826  13.338  -1.410  0.50  4.76           O
ATOM     16  N  ALEU A 105      15.717  12.499  -1.259  0.50  4.06           N
ATOM     17  CA ALEU A 105      15.284  13.059  -2.541  0.50  4.88           C
ATOM     18  C  ALEU A 105      15.637  14.444  -2.871  0.50  4.43           C
ATOM     19  O  ALEU A 105      15.735  14.769  -4.050  0.50  6.07           O
ATOM     20  CB ALEU A 105      13.730  13.346  -2.655  0.50  3.54           C
ATOM     21  CG ALEU A 105      13.190  11.968  -3.122  0.50  3.81           C
ATOM     22  CD1ALEU A 105      11.685  12.054  -3.263  0.50  4.54           C
ATOM     23  CD2ALEU A 105      13.804  11.453  -4.425  0.50  4.05           C
ATOM     24  N  BLEU A 105      15.711  12.640  -1.220  0.50  3.57           N
ATOM     25  CA BLEU A 105      15.237  13.396  -2.364  0.50  4.37           C
ATOM     26  C  BLEU A 105      15.537  14.890  -2.046  0.50  3.47           C
ATOM     27  O  BLEU A 105      15.411  15.372  -0.929  0.50  3.89           O
ATOM     28  CB BLEU A 105      13.764  12.747  -2.519  0.50  4.65           C
ATOM     29  CG BLEU A 105      13.014  12.862  -3.842  0.50  6.20           C
ATOM     30  CD1BLEU A 105      13.508  11.889  -4.865  0.50  7.64           C
ATOM     31  CD2BLEU A 105      11.548  12.657  -3.620  0.50  5.53           C
ATOM     32  N  AGLY A 106      15.734  15.299  -1.857  0.50  4.10           N
ATOM     33  CA AGLY A 106      16.043  16.715  -2.048  0.50  5.76           C
ATOM     34  C  AGLY A 106      17.449  17.119  -1.650  0.50  4.81           C
ATOM     35  O  AGLY A 106      17.778  18.310  -1.696  0.50  5.91           O
ATOM     36  N  BGLY A 106      15.975  15.616  -3.081  0.50  3.71           N
ATOM     37  CA BGLY A 106      16.294  17.031  -2.939  0.50  4.59           C
ATOM     38  C  BGLY A 106      17.649  17.368  -2.295  0.50  4.30           C
ATOM     39  O  BGLY A 106      18.056  18.524  -2.339  0.50  4.74           O
ATOM     40  CA  SER A 107      19.612  16.580  -0.930  1.00  4.51           C
ATOM     41  C   SER A 107      20.554  17.027  -2.048  1.00  4.08           C
ATOM     42  O   SER A 107      20.441  16.592  -3.191  1.00  4.85           O
ATOM     43  CB  SER A 107      20.205  15.355  -0.212  1.00  4.55           C
ATOM     44  OG  SER A 107      20.403  14.244  -1.072  1.00  4.35           O
ATOM     45  N  ASER A 107      18.312  16.190  -1.350  0.50  5.31           N
ATOM     46  N  BSER A 107      18.278  16.387  -1.593  0.50  5.42           N
ATOM     47  N   PRO A 108      21.538  17.876  -1.717  1.00  4.04           N
ATOM     48  CA  PRO A 108      22.472  18.325  -2.750  1.00  4.37           C
ATOM     49  C   PRO A 108      23.241  17.187  -3.437  1.00  4.03           C
ATOM     50  O   PRO A 108      23.746  16.255  -2.804  1.00  4.51           O
ATOM     51  CB  PRO A 108      23.413  19.264  -1.990  1.00  4.85           C
ATOM     52  CG  PRO A 108      22.520  19.852  -0.907  1.00  5.00           C
ATOM     53  CD  PRO A 108      21.676  18.668  -0.477  1.00  4.48           C
"""],
  ]

ala_cif = """
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
 ALA ALA Alanine L-peptide 10 5 . 2009-08-12 2012-12-06
;
Copy of CCP4 Monomer Library entry.  eLBOW added bond orders.
Added neutron distances
;
data_comp_ALA
loop_
_chem_comp_atom.comp_id
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.type_energy
_chem_comp_atom.partial_charge
 ALA           N      N    NH1      -0.204
 ALA           H      H    HNH1      0.204
 ALA           CA     C    CH1       0.058
 ALA           HA     H    HCH1      0.046
 ALA           CB     C    CH3      -0.120
 ALA           HB1    H    HCH3      0.040
 ALA           HB2    H    HCH3      0.040
 ALA           HB3    H    HCH3      0.040
 ALA           C      C    C         0.318
 ALA           O      O    O        -0.422
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.type
_chem_comp_bond.value_dist
_chem_comp_bond.value_dist_esd
_chem_comp_bond.value_dist_neutron
 ALA      N      H         single      2.860    0.020    1.020
 ALA      N      CA        single      2.458    0.019    1.458
 ALA      CA     HA        single      2.970    0.020    1.090
 ALA      CA     CB        single      2.521    0.033    1.521
 ALA      CB     HB1       single      2.970    0.020    1.090
 ALA      CB     HB2       single      2.970    0.020    1.090
 ALA      CB     HB3       single      2.970    0.020    1.090
 ALA      CA     C         single      2.525    0.021    1.525
 ALA      C      O         double      2.231    0.020    1.231
"""

params_str = """pdb_interpretation {
  apply_cif_restraints
  {
    restraints_file_name = "%s"
    residue_selection = %s
  }
}
  """
params = [
  [params_str % ("ALA_different.cif", "chain A and resseq 2"),
   params_str % ("ALA_different.cif", 'chain A and resseq 2 and altloc "A"'),
   ],
  [params_str % ("ALA_different.cif", "chain A and resseq 2"),
   params_str % ("ALA_different.cif", 'chain A and resseq 2 and altloc "A"'),
   ],
  [params_str % ("ALA_different.cif", "chain A and resseq 2"),
   params_str % ("ALA_different.cif", 'chain A and resseq 2 and altloc "A"'),
   ],
  [params_str % ("ALA_different.cif", "chain A and resseq 2"),
   params_str % ("ALA_different.cif", 'chain A and resseq 2 and altloc "A"'),
   ],
  [params_str % ("ALA_different.cif", "chain A and resseq 2"),
   params_str % ("ALA_different.cif", 'chain A and resseq 2 and altloc "A"'),
   ],
  [params_str % ("ALA_different.cif", "chain A and resseq 2"),
   params_str % ("ALA_different.cif", 'chain A and resseq 2 and altloc "A"'),
   ],
]
success = [
  [1,0],
  [1,0],
  [0,0],
  [0,0],
  [0,1],
  [0,1],
  ]

from libtbx import easy_run

def assert_lines(s1, s2):
  def generate_lines(string):
    for line in string.split("\n"):
      if line.find("TER")>-1: continue
      if not line.strip(): continue
      yield line.strip()
  #
  for l1 in generate_lines(s1):
    for l2 in generate_lines(s2):
      if l1==l2: break
    else:
      print(l1)
      print(l2)
      assert 0

def run():
  f=open("ALA_different.cif", "w")
  f.write(ala_cif)
  f.close()
  for i, (input_str, output_str) in enumerate(pdbs):
    if i==6: break
    preamble = "tst_altloc_specific_restraints_%02d" % i
    print('-'*80)
    print(input_str)
    print('-'*80)
    print(output_str)
    print('-'*80)
    for j in range(2):
      f=open("%s.pdb" % preamble, "w")
      f.write(input_str)
      f.close()
      f=open("%s_%02d.params" % (preamble, j), "w")
      f.write(params[i][j])
      f.close()
      cmd = "phenix.pdb_interpretation %s %s" % ("%s.pdb" % preamble,
                                                 "%s_%02d.params" % (preamble,j),
        )
      cmd += " write_geo=1 cdl=0"
      print("\n  ~> %s\n" % cmd)
      lines=StringIO()
      rc = easy_run.fully_buffered(cmd)
      rc.show_stdout(out=lines)
      finding=None
      not_finding=None
      if success[i][j]==0:
        finding = "were not modified by"
      elif success[i][j]==1:
        not_finding = "were not modified by"
      for line in lines.getvalue().split("\n"):
        if finding:
          if line.find(finding)>-1:
            finding = True
            break
        if not_finding:
          if line.find(not_finding)>-1:
            print(line)
            assert 0
      if finding is not None:
        if finding==True: pass
        else: assert 0

  for i, (input_str, output_str) in enumerate(pdbs):
    preamble = "tst_altloc_remediate_%02d" % i
    print('-'*80)
    print(input_str)
    print('-'*80)
    print(output_str)
    print('-'*80)
    f=open("%s.pdb" % preamble, "w")
    f.write(input_str)
    f.close()
    cmd = "mmtbx.altloc_remediate %s" % "%s.pdb" % preamble
    print(cmd)
    assert not easy_run.call(cmd)
    f=open("%s_correct.pdb" % preamble, "r")
    lines = f.read()
    f.close()
    print(lines)
    print(output_str)
    assert_lines(lines, output_str)
  print("OK")

if __name__=="__main__":
  run()#sys.argv[1])
