from __future__ import division
import sys

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
ATOM      3  CB  ALA A   1      -3.346   2.715   1.747  1.00 20.00      A    C
ATOM      4  C   ALA A   1      -1.551   1.486   0.458  1.00 20.00      A    C
ATOM      5  O   ALA A   1      -2.007   1.439  -0.662  1.00 20.00      A    O
ATOM      6  N   ALA A   2      -0.744   0.404   0.952  1.00 20.00      A    N
ATOM      7  CA  ALA A   2      -0.493  -0.749   0.119  1.00 20.00      A    C
ATOM      8  CB  ALA A   2      -1.476  -1.885   0.331  1.00 20.00      A    C
ATOM      9  C   ALA A   2       0.972  -1.193   0.100  1.00 20.00      A    C
ATOM     10  O   ALA A   2       1.603  -1.231   1.129  1.00 20.00      A    O
ATOM     11  N   ALA A   3       1.569  -1.642  -1.120  1.00 20.00      A    N
ATOM     12  CA  ALA A   3       2.971  -1.947  -1.163  1.00 20.00      A    C
ATOM     13  CB  ALA A   3       3.520  -1.630  -2.540  1.00 20.00      A    C
ATOM     14  C   ALA A   3       3.187  -3.395  -0.790  1.00 20.00      A    C
ATOM     15  O   ALA A   3       2.646  -4.306  -1.466  1.00 20.00      A    O
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
ATOM      3  CB  ALA A   1      -3.346   2.715   1.747  1.00 20.00      A    C
ATOM      4  C   ALA A   1      -1.551   1.486   0.458  1.00 20.00      A    C
ATOM      5  O   ALA A   1      -2.007   1.439  -0.662  1.00 20.00      A    O
ATOM      6  N   ALA A   2      -0.744   0.404   0.952  1.00 20.00      A    N
ATOM      7  CA  ALA A   2      -0.493  -0.749   0.119  1.00 20.00      A    C
ATOM      8  C   ALA A   2       0.972  -1.193   0.100  1.00 20.00      A    C
ATOM      9  O   ALA A   2       1.603  -1.231   1.129  1.00 20.00      A    O
ATOM     10  CB AALA A   2      -1.476  -1.885   0.331  0.50 20.00      A    C
ATOM     11  CB BALA A   2      -1.476  -1.885   0.331  0.50 20.00      A    C
ATOM     12  N   ALA A   3       1.569  -1.642  -1.120  1.00 20.00      A    N
ATOM     13  CA  ALA A   3       2.971  -1.947  -1.163  1.00 20.00      A    C
ATOM     14  CB  ALA A   3       3.520  -1.630  -2.540  1.00 20.00      A    C
ATOM     15  C   ALA A   3       3.187  -3.395  -0.790  1.00 20.00      A    C
ATOM     16  O   ALA A   3       2.646  -4.306  -1.466  1.00 20.00      A    O
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
ATOM     19  CB  ALA A   3       3.520  -1.630  -2.540  1.00 20.00      A    C
ATOM     20  C   ALA A   3       3.187  -3.395  -0.790  1.00 20.00      A    C
ATOM     21  O   ALA A   3       2.646  -4.306  -1.466  1.00 20.00      A    O
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
ATOM     10  CB AALA A   2      -1.976  -1.885   0.331  0.50 20.00      A    C
ATOM     11  C  AALA A   2       0.972  -1.193   0.100  0.50 20.00      A    C
ATOM     12  O  AALA A   2       1.603  -1.231   1.129  0.50 20.00      A    O
ATOM     13  N  BALA A   2      -0.744   0.404   0.952  0.50 20.00      A    N
ATOM     14  CA BALA A   2      -0.493  -0.749   0.119  0.50 20.00      A    C
ATOM     15  CB BALA A   2      -1.976  -1.885   0.331  0.50 20.00      A    C
ATOM     16  C  BALA A   2       0.972  -1.193   0.100  0.50 20.00      A    C
ATOM     17  O  BALA A   2       1.603  -1.231   1.129  0.50 20.00      A    O
ATOM     18  CA  ALA A   3       2.971  -1.947  -1.163  1.00 20.00      A    C
ATOM     19  CB  ALA A   3       3.520  -1.630  -2.540  1.00 20.00      A    C
ATOM     20  C   ALA A   3       3.187  -3.395  -0.790  1.00 20.00      A    C
ATOM     21  O   ALA A   3       2.646  -4.306  -1.466  1.00 20.00      A    O
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
ATOM     10  CB AALA A   2      -1.976  -1.885   0.331  0.50 20.00      A    C
ATOM     11  C  AALA A   2       0.972  -1.193   0.100  0.50 20.00      A    C
ATOM     12  O  AALA A   2       1.603  -1.231   1.129  0.50 20.00      A    O
ATOM     13  N  BALA A   2      -0.744   0.404   0.952  0.50 20.00      A    N
ATOM     14  CA BALA A   2      -0.493  -0.749   0.119  0.50 20.00      A    C
ATOM     15  CB BALA A   2      -1.976  -1.885   0.331  0.50 20.00      A    C
ATOM     16  C  BALA A   2       0.972  -1.193   0.100  0.50 20.00      A    C
ATOM     17  O  BALA A   2       1.603  -1.231   1.129  0.50 20.00      A    O
ATOM     18  CA  ALA A   3       2.971  -1.947  -1.163  1.00 20.00      A    C
ATOM     19  CB  ALA A   3       3.520  -1.630  -2.540  1.00 20.00      A    C
ATOM     20  C   ALA A   3       3.187  -3.395  -0.790  1.00 20.00      A    C
ATOM     21  O   ALA A   3       2.646  -4.306  -1.466  1.00 20.00      A    O
ATOM     22  OXT ALA A   3       3.917  -3.690   0.195  1.00 20.00      A    O-1
ATOM     23  N  AALA A   3       1.569  -1.642  -1.120  0.50 20.00      A    N
ATOM     24  N  BALA A   3       1.569  -1.642  -1.120  0.50 20.00      A    N
"""],
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
      print l1
      print l2
      assert 0

def run():
  for i, (input_str, output_str) in enumerate(pdbs):
    preamble = "tst_altloc_remediate_%02d" % i
    print '-'*80
    print input_str
    print '-'*80
    print output_str
    print '-'*80
    f=file("%s.pdb" % preamble, "wb")
    f.write(input_str)
    f.close()
    cmd = "mmtbx.altloc_remediate %s" % "%s.pdb" % preamble
    print cmd
    easy_run.call(cmd)
    f=file("%s_correct.pdb" % preamble, "rb")
    lines = f.read()
    f.close()
    print lines
    print output_str
    assert_lines(lines, output_str)
  print "OK"

if __name__=="__main__":
  run()#sys.argv[1])
