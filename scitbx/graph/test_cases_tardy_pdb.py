from __future__ import absolute_import, division, print_function
import scitbx.graph.utils as graph_utils
from scitbx import matrix
from libtbx.utils import Sorry

class pdb_extract(object):

  def __init__(O,
        tag,
        pdb,
        bonds,
        clusters=None,
        hinge_edges=None,
        loop_edges=None,
        loop_edge_bendings=None,
        find_cluster_loop_repeats=0,
        merge_clusters_with_multiple_connections_passes=1):
    O.index = None
    O.tag = tag
    O.labels, O.sites = [], []
    for line in pdb.splitlines():
      O.labels.append(line[22:26].strip()+"."+line[12:16].strip())
      O.sites.append(matrix.col([float(line[30+i*8:38+i*8]) for i in [0,1,2]]))
    O.bonds = bonds
    O.clusters = clusters
    O.hinge_edges = hinge_edges
    O.loop_edges = loop_edges
    O.loop_edge_bendings = loop_edge_bendings
    O.find_cluster_loop_repeats = find_cluster_loop_repeats
    O.merge_clusters_with_multiple_connections_passes = \
      merge_clusters_with_multiple_connections_passes

  def angles(O):
    return graph_utils.extract_edge_list(
      graph_utils.bond_bending_edge_sets(
        edge_sets=graph_utils.construct_edge_sets(
          n_vertices=len(O.sites),
          edge_list=O.bonds),
        omit_bonds=True))

  def tardy_tree_construct(O, fixed_vertex_lists=[]):
    from scitbx.graph import tardy_tree
    tt = tardy_tree.construct(
      sites=O.sites,
      edge_list=O.bonds,
      fixed_vertex_lists=fixed_vertex_lists)
    if (len(fixed_vertex_lists) != 0):
      return tt
    cm = tt.cluster_manager
    if (O.clusters is None):
      print("tag:", O.tag)
    if (O.clusters is None):
      print("clusters:", cm.clusters)
    else:
      assert cm.clusters == O.clusters
    if (O.hinge_edges is None):
      print("hinge_edges:", cm.hinge_edges)
    else:
      assert cm.hinge_edges == O.hinge_edges
    if (O.loop_edges is None):
      print("loop_edges:", cm.loop_edges)
    else:
      assert cm.loop_edges == O.loop_edges
    if (O.loop_edge_bendings is None):
      print("loop_edge_bendings:", cm.loop_edge_bendings)
    else:
      assert cm.loop_edge_bendings == O.loop_edge_bendings
    assert tt.find_cluster_loop_repeats == O.find_cluster_loop_repeats
    assert cm.merge_clusters_with_multiple_connections_passes \
        ==  O.merge_clusters_with_multiple_connections_passes
    if (O.clusters is None):
      print()
    return tt

test_cases = [

pdb_extract(
  tag="point",
  pdb="""\
ATOM      0  FE  FE  A   1      10.949  12.815  15.189  0.00  0.00          FE
""",
  bonds=[],
  clusters=[[0]],
  hinge_edges=[(-1, 0)],
  loop_edges=[],
  loop_edge_bendings=[]),

pdb_extract(
  tag="gly_no_h",
  pdb="""\
ATOM      0  N   GLY A   1      10.949  12.815  15.189  0.00  0.00           N
ATOM      1  CA  GLY A   1      10.405  13.954  15.917  0.00  0.00           C
ATOM      2  C   GLY A   1      10.779  15.262  15.227  0.00  0.00           C
ATOM      3  O   GLY A   1       9.916  16.090  14.936  0.00  0.00           O
""",
  bonds=[(0,1),(1,2),(2,3)],
  clusters=[[0, 1, 2], [3]],
  hinge_edges=[(-1, 1), (1, 2)],
  loop_edges=[],
  loop_edge_bendings=[]),

pdb_extract(
  tag="gly_with_nh",
  pdb="""\
ATOM      0  N   GLY A   1      10.949  12.815  15.189  0.00  0.00           N
ATOM      1  CA  GLY A   1      10.405  13.954  15.917  0.00  0.00           C
ATOM      2  C   GLY A   1      10.779  15.262  15.227  0.00  0.00           C
ATOM      3  O   GLY A   1       9.916  16.090  14.936  0.00  0.00           O
ATOM      4  H   GLY A   1      11.792  12.691  15.311  0.00  0.00           H
""",
  bonds=[(0,1),(0,4),(1,2),(2,3)],
  clusters=[[0, 1, 4], [2], [3]],
  hinge_edges=[(-1, 0), (0, 1), (1, 2)],
  loop_edges=[],
  loop_edge_bendings=[]),

pdb_extract(
  tag="ala_no_h",
  pdb="""\
ATOM      0  N   ALA A   1      10.949  12.815  15.189  0.00  0.00           N
ATOM      1  CA  ALA A   1      10.405  13.954  15.917  0.00  0.00           C
ATOM      2  C   ALA A   1      10.779  15.262  15.227  0.00  0.00           C
ATOM      3  CB  ALA A   1      10.908  13.950  17.351  0.00  0.00           C
ATOM      4  O   ALA A   1       9.916  16.090  14.936  0.00  0.00           O
""",
  bonds=[(0,1),(1,2),(1,3),(2,4)],
  clusters=[[0, 1, 2, 3], [4]],
  hinge_edges=[(-1, 1), (1, 2)],
  loop_edges=[],
  loop_edge_bendings=[]),

pdb_extract(
  tag="ala_with_h",
  pdb="""\
ATOM      0  N   ALA A   1      10.949  12.815  15.189  0.00  0.00           N
ATOM      1  CA  ALA A   1      10.405  13.954  15.917  0.00  0.00           C
ATOM      2  C   ALA A   1      10.779  15.262  15.227  0.00  0.00           C
ATOM      3  HA  ALA A   1       9.428  13.887  15.936  0.00  0.00           H
ATOM      4  O   ALA A   1       9.916  16.090  14.936  0.00  0.00           O
ATOM      5  H   ALA A   1      11.792  12.691  15.311  0.00  0.00           H
ATOM      6  CB  ALA A   1      10.908  13.950  17.351  0.00  0.00           C
ATOM      7  HB1 ALA A   1      10.627  13.138  17.778  0.00  0.00           H
ATOM      8  HB2 ALA A   1      10.540  14.707  17.813  0.00  0.00           H
ATOM      9  HB3 ALA A   1      11.867  14.004  17.346  0.00  0.00           H
""",
  bonds=[(0,1),(0,5),(1,2),(1,3),(1,6),(2,4),(6,7),(6,8),(6,9)],
  clusters=[[0, 1, 2, 3, 6], [7, 8, 9], [5], [4]],
  hinge_edges=[(-1, 1), (1, 6), (1, 0), (1, 2)],
  loop_edges=[],
  loop_edge_bendings=[]),

pdb_extract(
  tag="tyr_with_h",
  pdb="""\
ATOM      0  CG  TYR A   1      11.007   9.417   9.446  1.00  0.79           C
ATOM      1  CD1 TYR A   1       9.923  10.155   8.940  1.00  1.42           C
ATOM      2  CD2 TYR A   1      10.765   8.288  10.238  1.00  1.41           C
ATOM      3  CE1 TYR A   1       8.612   9.760   9.229  1.00  1.61           C
ATOM      4  CE2 TYR A   1       9.453   7.895  10.525  1.00  1.42           C
ATOM      5  CZ  TYR A   1       8.377   8.631  10.021  1.00  1.11           C
ATOM      6  HD1 TYR A   1      10.092  11.024   8.328  1.00  2.14           H
ATOM      7  HD2 TYR A   1      11.596   7.718  10.630  1.00  2.21           H
ATOM      8  HE1 TYR A   1       7.780  10.329   8.841  1.00  2.44           H
ATOM      9  HE2 TYR A   1       9.270   7.023  11.135  1.00  2.13           H
ATOM     10  OH  TYR A   1       7.083   8.244  10.304  1.00  1.32           O
ATOM     11  HH  TYR A   1       6.494   8.723   9.717  1.00  2.00           H
ATOM     12  CB  TYR A   1      12.440   9.818   9.148  1.00  0.74           C
ATOM     13  HB2 TYR A   1      12.827   9.193   8.358  1.00  0.78           H
ATOM     14  HB3 TYR A   1      13.036   9.677  10.037  1.00  0.78           H
ATOM     15  N   TYR A   1      11.593  12.101   9.550  1.00  0.82           N
ATOM     16  CA  TYR A   1      12.527  11.286   8.721  1.00  0.75           C
ATOM     17  C   TYR A   1      12.160  11.413   7.239  1.00  0.76           C
ATOM     18  HA  TYR A   1      13.536  11.638   8.870  1.00  0.85           H
ATOM     19  O   TYR A   1      12.298  12.462   6.643  1.00  0.83           O
ATOM     20  H   TYR A   1      10.948  12.701   9.122  1.00  0.88           H
""",
  bonds=[
    (0, 1), (0, 2), (0, 12), (1, 3), (1, 6), (2, 4), (2, 7), (3, 5),
    (3, 8), (4, 5), (4, 9), (5, 10), (10, 11), (12, 13), (12, 14),
    (12, 16), (15, 16), (15, 20), (16, 17), (16, 18), (17, 19)],
  clusters=[
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12],
    [13, 14, 16], [15, 17, 18], [11], [20], [19]],
  hinge_edges=[
    (-1, 0), (0, 12), (12, 16), (5, 10), (16, 15), (16, 17)],
  loop_edges=[],
  loop_edge_bendings=[]),

pdb_extract(
  tag="van_fragment", # PDB code 1qd8
  pdb="""\
HETATM   47  C44 VAN     1       0.718   5.411   2.269  1.00  3.66           C
HETATM   48  C47 VAN     1       0.913   3.899   2.010  1.00  3.90           C
HETATM   49  C48 VAN     1       1.937   3.222   2.700  1.00  3.87           C
HETATM   50  C50 VAN     1       2.013   1.800   2.726  1.00  4.31           C
HETATM   51  C52 VAN     1       1.044   1.071   2.004  1.00  4.71           C
HETATM   52  O53 VAN     1       1.077  -0.297   2.083  1.00  6.04           O
HETATM   53  C51 VAN     1       0.066   1.767   1.246  1.00  5.04           C
HETATM   54  C49 VAN     1       0.006   3.153   1.237  1.00  4.22           C
HETATM   55  C45 VAN     1      -0.704   5.507   2.917  1.00  3.80           C
HETATM   56  O46 VAN     1      -1.666   5.899   2.216  1.00  4.55           O
HETATM   57  N54 VAN     1      -0.869   5.098   4.194  1.00  4.02           N
HETATM   58  C55 VAN     1       0.141   4.679   5.156  1.00  4.23           C
HETATM   69  C56 VAN     1       0.045   3.151   5.360  1.00  4.51           C
HETATM   70  O57 VAN     1      -1.022   2.542   5.259  1.00  5.51           O
HETATM   71  N68 VAN     1       1.226   2.525   5.641  1.00  4.67           N
HETATM   72  C69 VAN     1       1.340   1.069   5.484  1.00  5.14           C
HETATM   73  C72 VAN     1       2.754   0.737   4.924  1.00  4.58           C
HETATM   74  C73 VAN     1       3.049   1.131   3.577  1.00  4.34           C
HETATM   75  C75 VAN     1       4.354   0.899   3.099  1.00  4.29           C
HETATM   76  O79 VAN     1       4.736   1.286   1.834  1.00  4.73           O
HETATM   77  C76 VAN     1       5.342   0.269   3.883  1.00  4.57           C
HETATM   78  C77 VAN     1       5.023  -0.107   5.195  1.00  4.56           C
HETATM   79  O78 VAN     1       5.912  -0.707   6.053  1.00  5.51           O
HETATM   80  C74 VAN     1       3.724   0.123   5.713  1.00  4.95           C
HETATM   81  C70 VAN     1       1.069   0.287   6.838  1.00  5.82           C
HETATM   82  O71 VAN     1       1.149   0.912   7.924  1.00  7.21           O
HETATM   83  O80 VAN     1       0.816  -0.957   6.693  1.00  8.10           O
""",
  bonds=[
    (0, 1), (0, 8), (1, 2), (1, 7), (2, 3), (3, 4), (3, 17), (4, 5), (4, 6),
    (6, 7), (8, 9), (8, 10), (10, 11), (11, 12), (12, 13), (12, 14), (14, 15),
    (15, 16), (15, 24), (16, 17), (16, 23), (17, 18), (18, 19), (18, 20),
    (20, 21), (21, 22), (21, 23), (24, 25), (24, 26)],
  clusters=[
    [3, 15, 16, 17, 18, 19, 20, 21, 22, 23],
    [0, 1, 2, 4, 5, 6, 7],
    [14, 24], [25, 26], [8], [9, 10], [11], [12], [13]],
  hinge_edges=[
    (-1, 16), (17, 3), (16, 15), (15, 24), (1, 0), (0, 8), (8, 10),
    (10, 11), (11, 12)],
  loop_edges=[(12, 14)],
  loop_edge_bendings=[(11, 14), (12, 15), (13, 14)]),

pdb_extract(
  tag="ZINC00000015", # C[C@@H](C(=O)[O-])[NH+](CCCl)CCCl
  pdb="""\
ATOM      1  C00 LIG A   1       2.209   1.054   0.461  1.00 20.00      A    C
ATOM      2  C01 LIG A   1       0.708   1.054   0.461  1.00 20.00      A    C
ATOM      3  C02 LIG A   1       0.094   1.054   1.832  1.00 20.00      A    C
ATOM      4  O03 LIG A   1      -0.670   1.994   2.178  1.00 20.00      A    O
ATOM      5  O04 LIG A   1       0.433   0.179   2.672  1.00 20.00      A    O-1
ATOM      6  N05 LIG A   1       0.121   0.117  -0.463  1.00 20.00      A    N+1
ATOM      7  C06 LIG A   1       0.808  -1.157  -0.363  1.00 20.00      A    C
ATOM      8  C07 LIG A   1       1.154  -1.925  -1.637  1.00 20.00      A    C
ATOM      9 CL08 LIG A   1       1.969  -3.508  -1.515  1.00 20.00      A   CL
ATOM     10  C09 LIG A   1      -1.281  -0.063  -0.144  1.00 20.00      A    C
ATOM     11  C10 LIG A   1      -2.233  -0.595  -1.214  1.00 20.00      A    C
ATOM     12 CL11 LIG A   1      -3.974  -0.726  -0.854  1.00 20.00      A   CL
ATOM     13  H12 LIG A   1       0.449   2.022   0.057  1.00 20.00      A    H
ATOM     14  H13 LIG A   1       0.214   0.498  -1.470  1.00 20.00      A    H
""",
  bonds=[
    (0, 1), (1, 2), (1, 5), (1, 12), (2, 3), (2, 4), (5, 6), (5, 9), (5, 13),
    (6, 7), (7, 8), (9, 10), (10, 11)],
  clusters=[[0, 1, 2, 5, 12], [6, 9, 13], [3, 4], [7], [8], [10], [11]],
  hinge_edges=[(-1, 1), (1, 5), (1, 2), (5, 6), (6, 7), (5, 9), (9, 10)],
  loop_edges=[],
  loop_edge_bendings=[]),

pdb_extract(
  tag="C1CCCCCCC1",
  pdb="""\
ATOM      1  C01 LIG A   1      -0.602   0.044  -1.758  1.00 20.00      A    C
ATOM      2  C02 LIG A   1       0.946   0.044  -1.758  1.00 20.00      A    C
ATOM      3  C03 LIG A   1       1.549   0.044  -0.333  1.00 20.00      A    C
ATOM      4  C04 LIG A   1       0.985   1.155   0.557  1.00 20.00      A    C
ATOM      5  C05 LIG A   1       0.430   0.619   1.870  1.00 20.00      A    C
ATOM      6  C06 LIG A   1      -0.536  -0.585   1.672  1.00 20.00      A    C
ATOM      7  C07 LIG A   1      -1.590  -0.345   0.547  1.00 20.00      A    C
ATOM      8  C08 LIG A   1      -1.184  -0.969  -0.799  1.00 20.00      A    C
ATOM      9 H011 LIG A   1      -0.949  -0.183  -2.755  1.00 20.00      A    H
ATOM     10 H012 LIG A   1      -0.952   1.028  -1.483  1.00 20.00      A    H
ATOM     11 H021 LIG A   1       1.293  -0.834  -2.283  1.00 20.00      A    H
ATOM     12 H022 LIG A   1       1.293   0.921  -2.283  1.00 20.00      A    H
ATOM     13 H031 LIG A   1       2.618   0.172  -0.411  1.00 20.00      A    H
ATOM     14 H032 LIG A   1       1.349  -0.911   0.130  1.00 20.00      A    H
ATOM     15 H041 LIG A   1       1.773   1.861   0.775  1.00 20.00      A    H
ATOM     16 H042 LIG A   1       0.196   1.666   0.025  1.00 20.00      A    H
ATOM     17 H051 LIG A   1       1.254   0.304   2.494  1.00 20.00      A    H
ATOM     18 H052 LIG A   1      -0.102   1.413   2.373  1.00 20.00      A    H
ATOM     19 H061 LIG A   1      -1.057  -0.766   2.600  1.00 20.00      A    H
ATOM     20 H062 LIG A   1       0.045  -1.461   1.425  1.00 20.00      A    H
ATOM     21 H071 LIG A   1      -2.531  -0.775   0.857  1.00 20.00      A    H
ATOM     22 H072 LIG A   1      -1.722   0.718   0.413  1.00 20.00      A    H
ATOM     23 H081 LIG A   1      -0.452  -1.743  -0.620  1.00 20.00      A    H
ATOM     24 H082 LIG A   1      -2.056  -1.415  -1.253  1.00 20.00      A    H
""",
  bonds=[
    (0, 1), (0, 7), (0, 8), (0, 9), (1, 2), (1, 10), (1, 11),
    (2, 3), (2, 12), (2, 13), (3, 4), (3, 14), (3, 15),
    (4, 5), (4, 16), (4, 17), (5, 6), (5, 18), (5, 19),
    (6, 7), (6, 20), (6, 21), (7, 22), (7, 23)],
  clusters=[
    [0, 1, 7, 8, 9],
    [2, 10, 11], [3, 12, 13], [4, 14, 15], [5, 16, 17], [6, 18, 19],
    [20, 21], [22, 23]],
  hinge_edges=[
    (-1, 0), (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (0, 7)],
  loop_edges=[(6, 7)],
  loop_edge_bendings=[(0, 6), (5, 7), (6, 22), (6, 23), (7, 20), (7, 21)]),

pdb_extract(
  tag="c1cc2ccc1CCc3ccc(cc3)CC2",
  pdb="""\
ATOM      1  C00 LIG A   1      -0.703  -1.478  -1.316  1.00 20.00      A    C
ATOM      2  C01 LIG A   1       0.705  -1.478  -1.316  1.00 20.00      A    C
ATOM      3  C02 LIG A   1       1.403  -1.478  -0.101  1.00 20.00      A    C
ATOM      4  C03 LIG A   1       0.695  -1.478   1.107  1.00 20.00      A    C
ATOM      5  C04 LIG A   1      -0.713  -1.478   1.095  1.00 20.00      A    C
ATOM      6  C05 LIG A   1      -1.400  -1.478  -0.101  1.00 20.00      A    C
ATOM      7  C06 LIG A   1      -2.732  -0.755  -0.050  1.00 20.00      A    C
ATOM      8  C07 LIG A   1      -2.724   0.768   0.057  1.00 20.00      A    C
ATOM      9  C08 LIG A   1      -1.389   1.474   0.106  1.00 20.00      A    C
ATOM     10  C09 LIG A   1      -0.698   1.663  -1.087  1.00 20.00      A    C
ATOM     11  C10 LIG A   1       0.700   1.663  -1.087  1.00 20.00      A    C
ATOM     12  C11 LIG A   1       1.392   1.474   0.106  1.00 20.00      A    C
ATOM     13  C12 LIG A   1       0.685   1.286   1.290  1.00 20.00      A    C
ATOM     14  C13 LIG A   1      -0.682   1.286   1.290  1.00 20.00      A    C
ATOM     15  C14 LIG A   1       2.727   0.768   0.057  1.00 20.00      A    C
ATOM     16  C15 LIG A   1       2.734  -0.755  -0.050  1.00 20.00      A    C
""",
  bonds=[
    (0, 1), (0, 5), (1, 2), (2, 3), (2, 15), (3, 4), (4, 5), (5, 6),
    (6, 7), (7, 8), (8, 9), (8, 13), (9, 10), (10, 11),
    (11, 12), (11, 14), (12, 13), (14, 15)],
  clusters=[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]],
  hinge_edges=[(-1, 0)],
  loop_edges=[],
  loop_edge_bendings=[],
  find_cluster_loop_repeats=1,
  merge_clusters_with_multiple_connections_passes=2),

pdb_extract(
  tag="ZINC04656480", # C[NH+]1CC2C[NH+](CC(C1)NC2=O)C
  pdb="""\
ATOM      1  C01 LIG A   1       3.508  -0.708  -0.412  1.00 20.00      A    C
ATOM      2  N02 LIG A   1       2.081  -0.708  -0.412  1.00 20.00      A    N+1
ATOM      3  C03 LIG A   1       1.365  -0.708   0.929  1.00 20.00      A    C
ATOM      4  C04 LIG A   1       0.082   0.043   1.236  1.00 20.00      A    C
ATOM      5  C05 LIG A   1      -1.306  -0.550   0.994  1.00 20.00      A    C
ATOM      6  N06 LIG A   1      -2.192  -0.518  -0.222  1.00 20.00      A    N+1
ATOM      7  C07 LIG A   1      -1.365   0.265  -1.189  1.00 20.00      A    C
ATOM      8  C08 LIG A   1       0.016   0.896  -1.052  1.00 20.00      A    C
ATOM      9  C09 LIG A   1       1.306   0.176  -1.376  1.00 20.00      A    C
ATOM     10  N10 LIG A   1       0.128   1.826  -0.060  1.00 20.00      A    N
ATOM     11  C11 LIG A   1       0.183   1.401   1.073  1.00 20.00      A    C
ATOM     12  O12 LIG A   1       0.029   2.140   1.997  1.00 20.00      A    O
ATOM     13  C13 LIG A   1      -3.570  -0.330  -0.116  1.00 20.00      A    C
ATOM     14 H021 LIG A   1       1.870  -1.698  -0.789  1.00 20.00      A    H
ATOM     15 H061 LIG A   1      -2.133  -1.527  -0.603  1.00 20.00      A    H
""",
  bonds=[
    (0, 1), (1, 2), (1, 8), (1, 13), (2, 3), (3, 4), (3, 10), (4, 5),
    (5, 6), (5, 12), (5, 14), (6, 7), (7, 8), (7, 9), (9, 10), (10, 11)],
  clusters=[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]],
  hinge_edges=[(-1, 1)],
  loop_edges=[],
  loop_edge_bendings=[]),

pdb_extract(
  tag="ZINC00196949", # C1C[NH+]2C[NH+]3CC[NH+](C2)C[NH+]1C3
  pdb="""\
ATOM      1  C01 LIG A   1      -0.794  -0.210  -2.068  1.00 20.00      A    C
ATOM      2  C02 LIG A   1       0.777  -0.210  -2.068  1.00 20.00      A    C
ATOM      3  N03 LIG A   1       1.469  -0.210  -0.777  1.00 20.00      A    N+1
ATOM      4  C04 LIG A   1       1.272   1.121   0.048  1.00 20.00      A    C
ATOM      5  N05 LIG A   1      -0.055   1.295   0.938  1.00 20.00      A    N+1
ATOM      6  C06 LIG A   1      -0.049   0.433   2.224  1.00 20.00      A    C
ATOM      7  C07 LIG A   1      -0.000  -0.931   1.970  1.00 20.00      A    C
ATOM      8  N08 LIG A   1       0.028  -1.145   0.699  1.00 20.00      A    N+1
ATOM      9  C09 LIG A   1       1.169  -1.328   0.066  1.00 20.00      A    C
ATOM     10  C10 LIG A   1      -1.086  -1.335   0.014  1.00 20.00      A    C
ATOM     11  N11 LIG A   1      -1.435  -0.222  -0.800  1.00 20.00      A    N+1
ATOM     12  C12 LIG A   1      -1.284   1.103   0.012  1.00 20.00      A    C
ATOM     13 H031 LIG A   1       2.523  -0.274  -1.003  1.00 20.00      A    H
ATOM     14 H051 LIG A   1      -0.075   2.324   1.264  1.00 20.00      A    H
ATOM     15 H081 LIG A   1       0.030  -0.086   0.488  1.00 20.00      A    H
ATOM     16 H111 LIG A   1      -2.490  -0.325  -1.006  1.00 20.00      A    H
""",
  bonds=[
    (0, 1), (0, 10), (1, 2), (2, 3), (2, 8), (2, 12), (3, 4),
    (4, 5), (4, 11), (4, 13), (5, 6), (6, 7), (7, 8), (7, 9), (7, 14),
    (9, 10), (10, 11), (10, 15)],
  clusters=[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]],
  hinge_edges=[(-1, 0)],
  loop_edges=[],
  loop_edge_bendings=[]),

pdb_extract(
  tag="ZINC03847120", # c1ccc2ccccc(c1)C2(C#N)Br
  pdb="""\
ATOM      1  C01 LIG A   1      -0.432   0.012  -2.178  1.00 20.00      A    C
ATOM      2  C02 LIG A   1       0.935   0.041  -2.063  1.00 20.00      A    C
ATOM      3  C03 LIG A   1       1.687   0.100  -0.911  1.00 20.00      A    C
ATOM      4  C04 LIG A   1       1.243   0.139   0.405  1.00 20.00      A    C
ATOM      5  C05 LIG A   1       1.591  -1.206   0.440  1.00 20.00      A    C
ATOM      6  C06 LIG A   1       0.763  -2.304   0.362  1.00 20.00      A    C
ATOM      7  C07 LIG A   1      -0.604  -2.323   0.236  1.00 20.00      A    C
ATOM      8  C08 LIG A   1      -1.462  -1.249   0.162  1.00 20.00      A    C
ATOM      9  C09 LIG A   1      -1.152   0.105   0.202  1.00 20.00      A    C
ATOM     10  C10 LIG A   1      -1.368   0.038  -1.168  1.00 20.00      A    C
ATOM     11  C11 LIG A   1       0.001   0.809   0.965  1.00 20.00      A    C
ATOM     12  C12 LIG A   1      -0.253   2.219   0.533  1.00 20.00      A    C
ATOM     13  N13 LIG A   1      -0.445   3.287   0.205  1.00 20.00      A    N
ATOM     14 BR14 LIG A   1      -0.504   0.332   2.810  1.00 20.00      A   BR
""",
  bonds=[
    (0, 1), (0, 9), (1, 2), (2, 3), (3, 4), (3, 10), (4, 5), (5, 6),
    (6, 7), (7, 8), (8, 9), (8, 10), (10, 11), (10, 13), (11, 12)],
  clusters=[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]],
  hinge_edges=[(-1, 0)],
  loop_edges=[],
  loop_edge_bendings=[]),

pdb_extract(
  tag="-ZINC12362269", # C/C(=N\NC(=O)[C@H]1CC1(C2CC2)C3CC3)/c4cccc(c4)Br
                       # atoms 16-23 removed
  pdb="""\
ATOM      1  C01 LIG A   1      -1.077  -0.498  -3.089  1.00 20.00      A    C
ATOM      2  C02 LIG A   1      -1.796  -0.766  -1.768  1.00 20.00      A    C
ATOM      3  N03 LIG A   1      -1.180  -0.649  -0.687  1.00 20.00      A    N
ATOM      4  N04 LIG A   1       0.179  -0.269  -0.687  1.00 20.00      A    N
ATOM      5  C05 LIG A   1       0.896  -0.141   0.569  1.00 20.00      A    C
ATOM      6  O06 LIG A   1       0.341  -0.355   1.593  1.00 20.00      A    O
ATOM      7  C07 LIG A   1       2.307   0.249   0.569  1.00 20.00      A    C
ATOM      8  C08 LIG A   1       2.760  -0.479   1.787  1.00 20.00      A    C
ATOM      9  C09 LIG A   1       3.212   0.911   1.550  1.00 20.00      A    C
ATOM     10  C10 LIG A   1       3.211   2.174   2.408  1.00 20.00      A    C
ATOM     11  C11 LIG A   1       3.039   3.693   2.491  1.00 20.00      A    C
ATOM     12  C12 LIG A   1       2.935   2.841   3.758  1.00 20.00      A    C
ATOM     13  C13 LIG A   1       4.738   0.961   1.496  1.00 20.00      A    C
ATOM     14  C14 LIG A   1       6.113   0.494   1.979  1.00 20.00      A    C
ATOM     15  C15 LIG A   1       6.006   1.517   0.846  1.00 20.00      A    C
""",
  bonds=[
    (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (4, 6), (6, 7), (6, 8),
    (7, 8), (8, 9), (8, 12), (9, 10), (9, 11), (10, 11),
    (12, 13), (12, 14), (13, 14)],
  clusters=[[4, 6, 7, 8, 9, 12], [10, 11], [13, 14], [3, 5], [2], [1], [0]],
  hinge_edges=[(-1, 6), (8, 9), (8, 12), (6, 4), (4, 3), (3, 2), (2, 1)],
  loop_edges=[],
  loop_edge_bendings=[]),

pdb_extract(
  tag="ZINC01638508", # C1CC1C(C2CC2)(C3CC3)O
  pdb="""\
ATOM      1  C01 LIG A   1      -0.942   1.244  -2.155  1.00 20.00      A    C
ATOM      2  C02 LIG A   1       0.586   1.277  -2.221  1.00 20.00      A    C
ATOM      3  C03 LIG A   1      -0.120   1.216  -0.865  1.00 20.00      A    C
ATOM      4  C04 LIG A   1      -0.056  -0.069  -0.042  1.00 20.00      A    C
ATOM      5  C05 LIG A   1      -0.281   0.259   1.432  1.00 20.00      A    C
ATOM      6  C06 LIG A   1      -1.658   0.020   2.055  1.00 20.00      A    C
ATOM      7  C07 LIG A   1      -1.166   1.448   1.812  1.00 20.00      A    C
ATOM      8  C08 LIG A   1       1.316  -0.717  -0.216  1.00 20.00      A    C
ATOM      9  C09 LIG A   1       1.940  -1.483   0.953  1.00 20.00      A    C
ATOM     10  C10 LIG A   1       1.427  -2.242  -0.273  1.00 20.00      A    C
ATOM     11  O11 LIG A   1      -1.048  -0.953  -0.480  1.00 20.00      A    O
""",
  bonds=[
    (0, 1), (0, 2), (1, 2), (2, 3), (3, 4), (3, 7), (3, 10),
    (4, 5), (4, 6), (5, 6), (7, 8), (7, 9), (8, 9)],
  clusters=[[2, 3, 4, 7, 10], [0, 1], [5, 6], [8, 9]],
  hinge_edges=[(-1, 3), (3, 2), (3, 4), (3, 7)],
  loop_edges=[],
  loop_edge_bendings=[]),

pdb_extract(
  tag="ZINC01638509", # C1CC(OC1)(C2CC2)C3CC3
  pdb="""\
ATOM      1  C01 LIG A   1      -1.976   0.015  -1.638  1.00 20.00      A    C
ATOM      2  C02 LIG A   1      -0.474   0.091  -1.532  1.00 20.00      A    C
ATOM      3  C03 LIG A   1      -0.172   0.077  -0.038  1.00 20.00      A    C
ATOM      4  O04 LIG A   1      -1.402   0.464   0.622  1.00 20.00      A    O
ATOM      5  C05 LIG A   1      -2.482  -0.041  -0.201  1.00 20.00      A    C
ATOM      6  C06 LIG A   1       0.849   1.194   0.167  1.00 20.00      A    C
ATOM      7  C07 LIG A   1       2.063   2.028  -0.248  1.00 20.00      A    C
ATOM      8  C08 LIG A   1       1.374   2.295   1.091  1.00 20.00      A    C
ATOM      9  C09 LIG A   1       0.546  -1.212   0.356  1.00 20.00      A    C
ATOM     10  C10 LIG A   1       1.137  -2.583   0.018  1.00 20.00      A    C
ATOM     11  C11 LIG A   1       0.535  -2.329   1.402  1.00 20.00      A    C
""",
  bonds=[
    (0, 1), (0, 4), (1, 2), (2, 3), (2, 5), (2, 8), (3, 4),
    (5, 6), (5, 7), (6, 7), (8, 9), (8, 10), (9, 10)],
  clusters=[[0, 1, 2, 3, 4, 5, 8], [6, 7], [9, 10]],
  hinge_edges=[(-1, 0), (2, 5), (2, 8)],
  loop_edges=[],
  loop_edge_bendings=[]),

pdb_extract(
  tag="ZINC00057527", # COc1cc2c(c(c1OC)OC)-c3ccc(c(=O)cc3[C@H](CC2)[NH3+])O
  pdb="""\
ATOM      1  C01 LIG A   1      -3.707  -1.073  -3.995  1.00 20.00      A    C
ATOM      2  O02 LIG A   1      -2.733  -0.091  -3.791  1.00 20.00      A    O
ATOM      3  C03 LIG A   1      -1.979  -0.103  -2.612  1.00 20.00      A    C
ATOM      4  C04 LIG A   1      -2.587  -0.442  -1.394  1.00 20.00      A    C
ATOM      5  C05 LIG A   1      -1.817  -0.515  -0.228  1.00 20.00      A    C
ATOM      6  C06 LIG A   1      -0.427  -0.384  -0.301  1.00 20.00      A    C
ATOM      7  C07 LIG A   1       0.172   0.035  -1.496  1.00 20.00      A    C
ATOM      8  C08 LIG A   1      -0.601   0.164  -2.656  1.00 20.00      A    C
ATOM      9  O09 LIG A   1      -0.012   0.608  -3.845  1.00 20.00      A    O
ATOM     10  C10 LIG A   1      -0.529   1.736  -4.489  1.00 20.00      A    C
ATOM     11  O11 LIG A   1       1.495   0.490  -1.488  1.00 20.00      A    O
ATOM     12  C12 LIG A   1       2.426  -0.085  -2.359  1.00 20.00      A    C
ATOM     13  C13 LIG A   1       0.459   0.004   0.693  1.00 20.00      A    C
ATOM     14  C14 LIG A   1       1.662  -0.498   0.151  1.00 20.00      A    C
ATOM     15  C15 LIG A   1       2.869  -0.758   0.803  1.00 20.00      A    C
ATOM     16  C16 LIG A   1       3.180  -0.590   2.138  1.00 20.00      A    C
ATOM     17  C17 LIG A   1       2.375  -0.122   3.157  1.00 20.00      A    C
ATOM     18  O18 LIG A   1       2.907   0.021   4.205  1.00 20.00      A    O
ATOM     19  C19 LIG A   1       1.045   0.298   3.092  1.00 20.00      A    C
ATOM     20  C20 LIG A   1       0.191   0.348   1.967  1.00 20.00      A    C
ATOM     21  C21 LIG A   1      -1.251   0.702   2.179  1.00 20.00      A    C
ATOM     22  C22 LIG A   1      -2.144  -0.498   2.330  1.00 20.00      A    C
ATOM     23  C23 LIG A   1      -2.468  -1.194   0.991  1.00 20.00      A    C
ATOM     24  N24 LIG A   1      -1.457   1.705   3.192  1.00 20.00      A    N+1
ATOM     25  O25 LIG A   1       4.488  -0.926   2.503  1.00 20.00      A    O
ATOM     26 H211 LIG A   1      -1.563   1.168   1.254  1.00 20.00      A    H
""",
  bonds=[
    (0, 1), (1, 2), (2, 3), (2, 7), (3, 4), (4, 5), (4, 22), (5, 6), (5, 12),
    (6, 7), (6, 10), (7, 8), (8, 9), (10, 11), (12, 13), (12, 19), (13, 14),
    (14, 15), (15, 16), (15, 24), (16, 17), (16, 18), (18, 19), (19, 20),
    (20, 21), (20, 23), (20, 25), (21, 22)],
  clusters=[
    [1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 13, 14, 15, 16, 17, 18, 19,
     20, 21, 22, 23, 24, 25], [0], [9], [11]],
  hinge_edges=[(-1, 2), (2, 1), (7, 8), (6, 10)],
  loop_edges=[],
  loop_edge_bendings=[],
  find_cluster_loop_repeats=1,
  merge_clusters_with_multiple_connections_passes=2),

pdb_extract(
  # scitbx/rigid_body/essence/tst_tardy.py, exercise_near_singular_hinges()
  tag="collinear",
  pdb="""\
ATOM      1  C1  COL     1      -1.334  -0.770   0.000  1.00  0.00
ATOM      2  C2  COL     1      -1.334   0.770   0.000  1.00  0.00
ATOM      3  C3  COL     1       0.000   0.000   0.000  1.00  0.00
ATOM      4  C4  COL     1       1.540   0.000   0.000  1.00  0.00
ATOM      5  C5  COL     1       3.080   0.000   0.000  1.00  0.00
ATOM      6  C6  COL     1       4.620   0.000   0.000  1.00  0.00
ATOM      7  C7  COL     1       5.954  -0.770   0.000  1.00  0.00
""",
  bonds=[(0, 1), (0, 2), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6)],
  clusters=[[0, 1, 2, 3, 4, 5], [6]],
  hinge_edges=[(-1, 0), (4, 5)],
  loop_edges=[],
  loop_edge_bendings=[]),

]

def _set_test_case_indices():
  for i,test_case in enumerate(test_cases):
    test_case.index = i
_set_test_case_indices()

def select_test_cases(tags_or_indices):
  result = []
  test_case_dict = dict([(test_case.tag,i)
    for i,test_case in enumerate(test_cases)])
  for tag_or_index in tags_or_indices:
    i = test_case_dict.get(tag_or_index)
    if (i is None):
      try: i = int(tag_or_index)
      except ValueError: i = -1
      if (i < 0 or i >= len(test_cases)):
        msg = "No such test case: %s" % tag_or_index
        print(msg)
        print("List of available test cases:")
        for i,test_case in enumerate(test_cases):
          print("  %2d: %s" % (i, test_case.tag))
        raise Sorry(msg)
    result.append(test_cases[i])
  if (len(result) == 0):
    result = test_cases
  return result
