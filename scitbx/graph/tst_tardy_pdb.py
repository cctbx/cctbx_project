from scitbx import matrix

class pdb_extract(object):

  def __init__(O,
        tag,
        pdb,
        bonds,
        clusters=None,
        hinge_edges=None,
        loop_edges=None,
        loop_edge_bendings=None):
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

  def tardy_tree_construct(O, rigid_loop_size_max=8):
    from scitbx.graph import tardy_tree
    tt = tardy_tree.construct(
      n_vertices=len(O.sites),
      edge_list=O.bonds,
      rigid_loop_size_max=rigid_loop_size_max).finalize()
    cm = tt.cluster_manager
    if (O.clusters is None):
      print "tag:", O.tag
      print "clusters:", cm.clusters
    else:
      assert cm.clusters == O.clusters
    if (O.hinge_edges is None):
      print "hinge_edges:", cm.hinge_edges
    else:
      assert cm.hinge_edges == O.hinge_edges
    if (O.loop_edges is None):
      print "loop_edges:", cm.loop_edges
    else:
      assert cm.loop_edges == O.loop_edges
    if (O.loop_edge_bendings is None):
      print "loop_edge_bendings:", cm.loop_edge_bendings
    else:
      assert cm.loop_edge_bendings == O.loop_edge_bendings
    if (O.clusters is None):
      print
    return tt

test_cases = [

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

]
