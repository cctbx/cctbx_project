from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
from libtbx.utils import search_for
from libtbx.test_utils import block_show_diff
import libtbx.load_env
from cStringIO import StringIO
import os

def exercise_pdb_string(mon_lib_srv, ener_lib):
  raw_records = """\
CRYST1   50.066   67.126   47.862  90.00  92.41  90.00 P 1 21 1
ATOM      0  N   MET     0      18.670  12.527  40.988  1.00 52.89           N
ATOM      1  CA  MET     0      17.631  13.191  40.117  1.00 52.89           C
ATOM      2  CB  MET     0      18.110  13.132  38.643  1.00 52.89           C
ATOM      3  CG  MET     0      18.621  11.738  38.198  1.00 52.60           C
ATOM      4  SD  MET     0      19.279  11.706  36.483  1.00 52.89           S
ATOM      5  CE  MET     0      20.263  13.253  36.435  1.00 52.89           C
ATOM      6  C   MET     0      16.228  12.530  40.276  1.00 50.97           C
ATOM      7  O   MET     0      16.082  11.507  40.963  1.00 52.89           O
ATOM      8  HA  MET     0      17.568  14.136  40.375  1.00 52.89           D
ATOM      9  HB1 MET     0      18.846  13.755  38.537  1.00 45.86           D
ATOM     10  HB2 MET     0      17.375  13.382  38.060  1.00 52.35           D
ATOM     11  HG1 MET     0      17.893  11.102  38.251  1.00 47.31           D
ATOM     12  HG2 MET     0      19.339  11.465  38.791  1.00 52.89           D
ATOM     13  HE1 MET     0      20.554  13.413  35.532  1.00 52.89           D
ATOM     14  HE2 MET     0      21.028  13.143  37.008  1.00 52.89           D
ATOM     15  HE3 MET     0      19.722  13.985  36.740  1.00 52.89           D
ATOM     16  H1  MET     0      19.528  12.921  40.813  1.00 52.89           D
ATOM     17  H2  MET     0      18.706  11.586  40.792  1.00 52.89           D
ATOM     18  H3  MET     0      18.445  12.648  41.914  1.00 52.89           D
END
""".splitlines()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=None,
    raw_records=raw_records,
    force_symmetry=True)
  xray_structure = processed_pdb_file.xray_structure()
  assert xray_structure.scattering_type_registry().type_count_dict() \
      == {"S": 1, "N": 1, "C": 5, "O": 1, "D": 11}
  raw_records = [line[:66] for line in raw_records]
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=None,
    raw_records=raw_records,
    force_symmetry=True)
  xray_structure = processed_pdb_file.xray_structure()
  assert xray_structure.scattering_type_registry().type_count_dict() \
      == {"S": 1, "N": 1, "C": 5, "O": 1, "H": 11}
  #
  raw_records = """\
HETATM    1  N   NH3     1       0.000   0.000   0.000  1.00  0.00           N
HETATM    2  N   NH3     2       0.000   0.000   0.000  1.00  0.00
HETATM    3 NH3  NH3     3       0.000   0.000   0.000  1.00  0.00
HETATM    4  X   CH4     4       0.000   0.000   0.000  1.00  0.00           C
HETATM    5  C   CH4     5       0.000   0.000   0.000  1.00  0.00
HETATM    6  CH4 CH4     6       0.000   0.000   0.000  1.00  0.00
END
""".splitlines()
  log = StringIO()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=None,
    raw_records=raw_records,
    log=log)
  looking_for = "Ad-hoc single atom residues: "
  for line in log.getvalue().splitlines():
    line = line.strip()
    if (not line.startswith(looking_for)): continue
    counts = eval(line[len(looking_for):])
    assert counts == {"CH4": 3, "NH3": 3}
    break
  else:
    raise RuntimeError('Expected string not found in output: "%s"'
      % looking_for)
  #
  raw_records = """\
CRYST1  106.820   62.340  114.190  90.00  90.00  90.00 P 21 21 21
ATOM      7  N   SER A   4      64.059  32.579  18.554  1.00 14.95
ATOM      8  CA ASER A   4      64.561  31.418  17.802  1.00  9.04
ATOM      9  CA BSER A   4      64.804  31.635  18.406  1.00 12.07
ATOM     10  C   SER A   4      65.282  31.515  17.187  1.00 10.69
ATOM     38  N   ALA C   3     109.043  27.391  28.663  1.00 15.05
ATOM     39  CA  ALA C   3     109.073  26.531  28.433  1.00  3.14
ATOM     40  C   ALA C   3     108.930  26.867  26.637  1.00 15.22
ATOM     41  N   SER C   4     109.311  25.012  25.749  1.00  6.18
ATOM     42  CA BSER C   4     109.271  25.165  25.154  1.00 11.32
ATOM     43  CA ASER C   4     109.508  25.337  25.273  1.00  4.97
ATOM     44  C   SER C   4     109.045  24.408  24.187  1.00 15.28
END
""".splitlines()
  log = StringIO()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=None,
    raw_records=raw_records,
    log=log)
  lines = search_for(
    pattern=" +bond proxies already assigned to first conformer: 3$",
    mode="re.match",
    lines=log.getvalue().splitlines())
  assert len(lines) == 1
  #
  raw_records = """\
CRYST1   53.910   23.100   23.100  90.00 110.40  90.00 C 2
ATOM    561  N   TYR X  39      11.814  -3.041  22.005  1.00 15.52           N
ATOM    562  CA  TYR X  39      13.044  -2.727  21.259  1.00 17.08           C
ATOM    563  C   TYR X  39      13.804  -1.593  21.948  1.00 19.46           C
ATOM    564  O   TYR X  39      14.993  -1.524  21.704  1.00 28.79           O
ATOM    582  N   VAL X  40      13.093  -0.787  22.711  1.00 23.32           N
ATOM    583  CA  VAL X  40      13.567   0.356  23.491  1.00 37.02           C
ATOM    584  C   VAL X  40      13.408   0.141  24.979  1.00 47.42           C
ATOM    585  O   VAL X  40      13.919   0.926  25.809  1.00 58.31           O
ATOM    589  OXT VAL X  40      12.720  -0.842  25.372  1.00 48.85           O
HETATM  600  C1 AEOH X 200      12.974   5.558  13.017  0.54 29.75           C
HETATM  601  C1 BEOH X 200      14.446   4.322  12.408  0.46 32.58           C
HETATM  602  C2 AEOH X 200      12.259   6.535  13.707  0.54 26.57           C
HETATM  603  C2 BEOH X 200      13.905   4.879  13.576  0.46 32.29           C
END
""".splitlines()
  log = StringIO()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=None,
    raw_records=raw_records,
    log=log)
  lines = search_for(
    pattern="""\
 +Inner-chain residues flagged as termini: \\['pdbres="VAL X  40 "'\\]$""",
    mode="re.match",
    lines=log.getvalue().splitlines())
  assert len(lines) == 2
  #
  raw_records = """\
CRYST1   53.910   23.100   23.100  90.00 110.40  90.00 C 2
ATOM    264  N   GLU    21       3.822  -2.789  -0.368  1.00 11.74           N
ATOM    265  CA  GLU    21       5.032  -3.595  -0.485  1.00 11.86           C
ATOM    266  C   GLU    21       6.226  -2.698  -0.821  1.00 10.70           C
ATOM    267  O   GLU    21       6.118  -1.718  -1.512  1.00 11.93           O
ATOM    268  CB AGLU    21       4.980  -4.649  -1.566  0.32 13.49           C
ATOM    269  CB BGLU    21       4.769  -4.484  -1.766  0.68 14.37           C
ATOM    270  CG AGLU    21       3.942  -5.699  -1.498  0.32 13.60           C
ATOM    271  CG BGLU    21       5.737  -5.573  -1.952  0.68 19.43           C
ATOM    272  CD AGLU    21       4.064  -6.621  -2.739  0.32 15.27           C
ATOM    273  CD BGLU    21       5.251  -6.534  -3.036  0.68 19.91           C
ATOM    274  OE1AGLU    21       5.043  -6.411  -3.478  0.32 16.46           O
ATOM    276  OE2AGLU    21       3.193  -7.474  -2.870  0.32 19.03           O
ATOM    277  OE2BGLU    21       6.037  -6.819  -3.934  0.68 23.06           O
END
""".splitlines()
  for i_pass in [0,1]:
    if (i_pass == 0):
      pattern = " +2.60 -     3.06: 2$"
      ci_counts = {0: 4, 1: 5, 2: 4}
    else:
      raw_records = raw_records[:1] + raw_records[5:]
      pattern = " +2.60 -     2.80: 1$"
      ci_counts = {1: 5, 2: 4}
    log = StringIO()
    processed_pdb_file = monomer_library.pdb_interpretation.process(
      mon_lib_srv=mon_lib_srv,
      ener_lib=ener_lib,
      file_name=None,
      raw_records=raw_records,
      log=log)
    grm = processed_pdb_file.geometry_restraints_manager()
    assert dict(grm.conformer_indices.counts()) == ci_counts
    lines = search_for(
      pattern=pattern, mode="re.match", lines=log.getvalue().splitlines())
    assert len(lines) == 1
  #
  raw_records = """\
REMARK    HYDROLASE                               05-SEP-07   2VB1
CRYST1   27.070   31.250   33.760  87.98 108.00 112.11 P 1           1
ATOM      1  N   LYS A   1       1.984   5.113  14.226  1.00  6.14           N
ATOM      2  CA  LYS A   1       1.811   6.069  13.092  1.00  5.23           C
ATOM      3  C   LYS A   1       2.502   7.339  13.503  1.00  4.84           C
ATOM      4  O   LYS A   1       2.414   7.746  14.659  1.00  5.74           O
ATOM      5  CB  LYS A   1       0.325   6.305  12.891  1.00  5.48           C
ATOM      6  CG  LYS A   1       0.013   7.371  11.851  1.00  5.12           C
ATOM      7  CD  LYS A   1      -1.494   7.455  11.617  1.00  5.56           C
ATOM      8  HA  LYS A   1       2.215   5.709  12.275  1.00  6.28           H
ATOM      9  HB2 LYS A   1      -0.091   5.472  12.620  1.00  6.57           H
ATOM     10  HB3 LYS A   1      -0.068   6.569  13.738  1.00  6.57           H
ATOM     11  HG2 LYS A   1       0.344   8.231  12.157  1.00  6.14           H
ATOM     12  HG3 LYS A   1       0.461   7.155  11.018  1.00  6.14           H
ATOM     13  CE ALYS A   1      -1.966   8.606  10.745  0.69  5.10           C
ATOM     14  NZ ALYS A   1      -1.548   8.473   9.287  0.69  4.56           N
ATOM     15  HD2ALYS A   1      -1.786   6.625  11.210  0.50  6.67           H
ATOM     16  HD3ALYS A   1      -1.933   7.526  12.479  0.50  6.67           H
ATOM     17  HE2ALYS A   1      -2.934   8.658  10.791  0.69  6.12           H
ATOM     18  HE3ALYS A   1      -1.610   9.436  11.100  0.69  6.12           H
ATOM     19  HZ1ALYS A   1      -1.884   7.721   8.949  0.69  6.84           H
ATOM     20  HZ2ALYS A   1      -1.854   9.169   8.825  0.69  6.84           H
ATOM     21  HZ3ALYS A   1      -0.659   8.450   9.234  0.69  6.84           H
ATOM     22  CE BLYS A   1      -1.791   8.575  10.680  0.31  5.79           C
ATOM     23  NZ BLYS A   1      -3.148   8.376  10.084  0.31  6.86           N
ATOM     24  H1 BLYS A   1       2.852   4.975  14.368  0.50  9.20           H
ATOM     25  H2 BLYS A   1       1.614   5.454  14.960  0.50  9.20           H
ATOM     26  H3 BLYS A   1       1.589   4.341  14.026  0.50  9.20           H
ATOM     27  HE2BLYS A   1      -1.763   9.419  11.157  0.31  6.95           H
ATOM     28  HE3BLYS A   1      -1.124   8.601   9.977  0.31  6.95           H
ATOM     29  HZ1BLYS A   1      -3.755   8.339  10.733  0.31 10.29           H
ATOM     30  HZ2BLYS A   1      -3.335   9.055   9.540  0.31 10.29           H
ATOM     31  HZ3BLYS A   1      -3.160   7.615   9.622  0.31 10.29           H
END
""".splitlines()
  log = StringIO()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=None,
    raw_records=raw_records,
    log=log)
  lines = search_for(
    pattern="""\
  Number of resolved chirality restraint conflicts: 1""",
    mode="==",
    lines=log.getvalue().splitlines())
  assert len(lines) == 1

def exercise_rna(
      mon_lib_srv,
      ener_lib,
      file_name,
      expected_block,
      expected_block_last_startswith=True,
      expected_modifications_used=None):
  file_path = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/"+file_name,
    test=os.path.isfile)
  if (file_path is None):
    print 'Skipping exercise_rna("%s"): input file not available' % file_name
    return
  log = StringIO()
  processed_pdb_file = monomer_library.pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name=file_path,
    log=log)
  lines = []
  lines_modifications_used = []
  for line in log.getvalue().splitlines():
    if (line.startswith("""\
          Modifications used: {""")):
      lines_modifications_used.append(line)
    else:
      lines.append(line)
  assert not block_show_diff(
    lines, expected_block, last_startswith=expected_block_last_startswith)
  if (expected_modifications_used is None):
    assert len(lines_modifications_used) == 0
  else:
    assert len(lines_modifications_used) == 1
    modifications_used = eval(lines_modifications_used[0].split(":", 1)[1])
    assert modifications_used == expected_modifications_used

def exercise_cns_rna(mon_lib_srv, ener_lib):
  exercise_rna(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name="cns_rna.pdb",
    expected_block= """\
  Total number of atoms: 646
  Number of models: 1
  Model: ""
    Number of chains: 1
    Chain: " "
      Number of atoms: 646
      Number of conformers: 1
      Conformer: ""
        Number of residues, atoms: 20, 646
          Classifications: {'RNA': 20}
          Link IDs: {'p': 19}
  Time building chain proxies: """,
    expected_modifications_used = {
      'rnaEsd': 19, 'p5*END': 1, 'rnaC3': 19, '3*END': 1})

def exercise_rna_3p_2p(mon_lib_srv, ener_lib):
  exercise_rna(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name="rna_3p_2p.pdb",
    expected_block= """\
  Total number of atoms: 63
  Number of models: 1
  Model: ""
    Number of chains: 1
    Chain: "A"
      Number of atoms: 63
      Number of conformers: 1
      Conformer: ""
        Number of residues, atoms: 3, 63
          Classifications: {'RNA': 3}
          Link IDs: {'p': 2}
  Time building chain proxies: """,
    expected_modifications_used = {'rnaEsd': 2, 'rnaC3': 1, 'rnaC2': 1})
  #
  exercise_rna(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    file_name="coords_rna_gap.pdb",
    expected_block= """\
  Total number of atoms: 90
  Number of models: 1
  Model: ""
    Number of chains: 1
    Chain: "A"
      Number of atoms: 90
      Number of conformers: 1
      Conformer: ""
        Number of residues, atoms: 5, 90
          Classifications: {'RNA': 5}
          Link IDs: {'p': 4}
          Chain breaks: 1
""",
    expected_block_last_startswith=False,
    expected_modifications_used={'rnaEsd': 2, 'rnaC3': 2})

def exercise():
  mon_lib_srv = monomer_library.server.server()
  ener_lib = monomer_library.server.ener_lib()
  exercise_pdb_string(mon_lib_srv, ener_lib)
  exercise_cns_rna(mon_lib_srv, ener_lib)
  exercise_rna_3p_2p(mon_lib_srv, ener_lib)
  print "OK"

if (__name__ == "__main__"):
  exercise()
