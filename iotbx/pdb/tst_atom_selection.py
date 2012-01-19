from iotbx import pdb
from cctbx.array_family import flex
from libtbx.test_utils import Exception_expected, show_diff

def exercise_selection():
  hierarchy = pdb.input(source_info=None, lines=flex.split_lines("""\
CRYST1   50.800   50.800  155.300  90.00  90.00  90.00 P 43 21 2     8
MODEL        1
ATOM      4  N   SER     1       8.753  29.755  61.685  1.00 49.13
ATOM      5  CA  SER     1       9.242  30.200  62.974  1.00 46.62
ANISOU    5  CA  SER     1    343    490   2719    -45   -169    617
ATOM      6  C   SER     1      10.453  29.500  63.579  1.00 41.99
ATOM      7  O   SER     1      10.593  29.607  64.814  1.00 43.24
ANISOU    7  O   SER     1    343    490   2719    -45   -169    617
ATOM      8  CB  SER     1       8.052  30.189  63.974  1.00 53.00
ATOM      9  OG  SER     1       7.294  31.409  63.930  1.00 57.79
ATOM     10  N   ARG     2      11.360  28.819  62.827  1.00 36.48
ATOM     11  CA  ARG     2      12.548  28.316  63.532  1.00 30.20
ATOM     12  C   ARG     2      13.502  29.501  63.500  1.00 25.54
ATOM     13  O   ARG     2      13.730  30.037  62.407  1.00 23.86
ATOM     14  CB  ARG     2      13.241  27.119  62.861  1.00 27.44
ATOM     15  CG  ARG     2      12.412  25.849  62.964  1.00 23.66
ATOM     16  CD  ARG     2      13.267  24.651  63.266  1.00 23.98
ATOM     17  NE  ARG     2      13.948  24.115  62.135  1.00 22.71
ATOM     18  CZ  ARG     2      15.114  23.487  62.201  1.00 21.38
ATOM     19  NH1 ARG     2      15.845  23.331  63.301  1.00 19.34
ATOM     20  NH2 ARG     2      15.575  23.030  61.051  1.00 26.66
ATOM     21  N   PRO     3J     13.947  29.997  64.680  1.00 22.94
ATOM     22  CA  PRO     3J     14.902  31.100  64.827  1.00 20.19
ATOM     23  C   PRO     3J     16.195  30.718  64.086  1.00 18.44
ATOM     24  O   PRO     3J     16.545  29.521  64.086  1.00 19.76
ATOM     25  CB  PRO     3J     15.133  31.218  66.313  1.00 19.17
ATOM     26  CG  PRO     3J     14.065  30.364  66.951  1.00 15.12
ATOM     27  CD  PRO     3J     13.816  29.289  65.966  1.00 19.56
ATOM     28  N  AILE     4      16.953  31.648  63.512  1.00 15.29
ATOM     29  CA AILE     4      18.243  31.372  62.859  1.00 14.32
ATOM     30  C  AILE     4      19.233  32.112  63.743  1.00 13.54
ATOM     31  O  AILE     4      19.105  33.315  64.009  1.00 11.84
ATOM     32  CB AILE     4      18.298  31.951  61.406  1.00 13.62
ATOM     33  CG1AILE     4      17.157  31.300  60.620  1.00 18.39
ATOM     34  CG2AILE     4      19.661  31.747  60.743  1.00 13.64
ATOM     35  CD1AILE     4      16.879  32.102  59.355  1.00 16.69
ATOM     28  N  BILE     4      16.953  31.648  63.512  1.00 15.29
ATOM     29  CA BILE     4      18.243  31.372  62.859  1.00 14.32
ATOM     30  C  BILE     4      19.233  32.112  63.743  1.00 13.54
ATOM     31  O  BILE     4      19.105  33.315  64.009  1.00 11.84
ATOM     32  CB BILE     4      18.298  31.951  61.406  1.00 13.62
ATOM     33  CG1BILE     4      17.157  31.300  60.620  1.00 18.39
ATOM     34  CG2BILE     4      19.661  31.747  60.743  1.00 13.64
ATOM1200035  CD1BILE     4      16.879  32.102  59.355  1.00 16.69
TER      36      ILE     4
ENDMDL
MODEL        2
HETATM 1451  PA  5GP H 187      29.875  44.488  69.823  1.00 19.62
HETATM 1452  O1A 5GP H 187      28.526  44.888  69.143  1.00 19.86
HETATM 1453  O2A 5GP H 187      30.764  44.617  68.702  1.00 23.42
HETATM 1454  O3A 5GP H 187      30.319  45.004  71.073  1.00 20.20
HETATM 1455  O5* 5GP H 187      29.683  43.016  70.027  1.00 20.32
HETATM 1456  C5* 5GP H 187      30.740  42.297  70.837  1.00 21.47
HETATM 1457  C4* 5GP H 187      30.677  40.747  70.770  1.00 21.56
HETATM 1458  O4* 5GP H 187      29.608  40.160  71.599  1.00 20.50
HETATM 1459  C3* 5GP H 187      30.547  40.121  69.352  1.00 20.18
HETATM 1460  O3* 5GP H 187      31.228  38.864  69.416  1.00 23.65
HETATM 1461  C2* 5GP H 187      29.031  39.871  69.248  1.00 18.78
HETATM 1462  O2* 5GP H 187      28.685  38.690  68.496  1.00 20.45
HETATM 1463  C1* 5GP H 187      28.634  39.641  70.688  1.00 17.09
HETATM 1464  N9  5GP H 187      27.238  39.525  71.076  1.00 15.35
HETATM 1465  C8  5GP H 187      26.330  40.535  70.852  1.00 12.57
HETATM 1466  N7' 5GP H 187      25.175  40.314  71.417  1.00 12.88
HETATM 1467  C5  5GP H 187      25.278  39.082  72.070  1.00 10.75
HETATM 1468  C6  5GP H 187      24.326  38.354  72.827  1.00  9.77
HETATM 1469  O6  5GP H 187      23.169  38.678  73.029  1.00  8.66
HETATM 1470  N1' 5GP H 187      24.836  37.190  73.270  1.00  9.67
HETATM 1471  C2  5GP H 187      26.075  36.701  73.001  1.00  9.84
HETATM 1472  N2  5GP H 187      26.361  35.490  73.520  1.00  9.77
HETATM 1473  N3  5GP H 187      27.005  37.353  72.310  1.00 10.31
HETATM 1474  C4  5GP H 187      26.583  38.559  71.844  1.00 12.50
ENDMDL
MODEL        3
HETATM 1475  S   SO4 S 188      31.424  42.923  60.396  1.00 55.69           S4+
HETATM 1476  O1  SO4 S 188      31.631  41.513  60.336  1.00 59.84           o1-
HETATM 1477  O2  SO4 S 188      32.533  43.699  59.932  1.00 49.98           O1-
HETATM 1478  O3  SO4 S 188      31.128  43.217  61.738  1.00 59.44           O1-
HETATM 1479  O4  SO4 S 188      30.353  43.201  59.539  1.00 60.54           O1-
HETATM 1480  O   HOH W 200      29.478  23.354  61.364  1.00  8.67      WATE
ATOM   2000  A1  AAA X   1       8.753  29.755  61.685  1.00 49.13
ATOM   2001  A2  AAA X   1       9.242  30.200  62.974  1.00 46.62
ATOM   2002  A1  BBB X   2      11.360  28.819  62.827  1.00 36.48
ATOM   2003  A2  BBB X   2      12.548  28.316  63.532  1.00 30.20
ATOM   2004  A1  AAA Y   1       8.753  29.755  61.685  1.00 49.13
ATOM   2005  A2  AAA Y   1       9.242  30.200  62.974  1.00 46.62
ATOM   2006  A1  CCC Y   5       9.242  30.200  62.974  1.00 46.62
ATOM   2007  A2  BBB Y   2      12.548  28.316  63.532  1.00 30.20
ATOM   2008  A1  AAA Z   1K      8.753  29.755  61.685  1.00 49.13
ATOM   2009  A1  BBB Z   2      11.360  28.819  62.827  1.00 36.48
ATOM   2010  A2  BBB Z   2      12.548  28.316  63.532  1.00 30.20
ATOM   2011  A1  AAAZZ   1K      8.753  29.755  61.685  1.00 49.13
ATOM   2012  A1  BBBZZ   2      11.360  28.819  62.827  1.00 36.48
ATOM   2013  A1  CCCZZ   5       9.242  30.200  62.974  1.00 46.62
ATOM   2014  A1  CCCZZA001       9.242  30.200  62.974  1.00 46.62
ATOM   2015  A1  CCCZZA002       9.242  30.200  62.974  1.00 46.62
ATOM   2016  A1  CCCZZA003       9.242  30.200  62.974  1.00 46.62
ATOM   2017  A1  AAAUU  1K       8.753  29.755  61.685  1.00 49.13
ENDMDL
END
""")).construct_hierarchy()
  sel_cache = hierarchy.atom_selection_cache()
  assert sel_cache.n_seq == hierarchy.atoms_size()
  isel = sel_cache.iselection
  assert isel("").size() == 0
  assert isel("all").size() == sel_cache.n_seq
  assert isel("none").size() == 0
  assert isel("optional none", optional=True).size() == 0
  assert isel("optional none", optional=False) is None
  assert isel("not all").size() == 0
  assert isel("not none").size() == sel_cache.n_seq
  assert list(isel(r"name c?\*")) == [45,46,48,50,52]
  assert list(isel(r"name 'C?\*'")) == []
  assert list(isel(r"name ' C?\*'")) == [45,46,48,50,52]
  assert list(isel(r"name ' c?\*'")) == [45,46,48,50,52]
  assert list(isel(r"name n?'")) == [55, 59]
  for conj in ["and ", ""]:
    assert list(isel(r"altloc a %sname n" % conj)) == [24]
    assert list(isel(r"altloc b %sname n" % conj)) == [32]
    assert list(isel(r"altloc ' ' %sname n" % conj)) == [0,6,17]
    assert list(isel(r"altid ' ' %sname n" % conj)) == [0,6,17]
  assert list(isel(r"resname hoh")) == [69]
  assert list(isel(r"resname SO4")) == [64,65,66,67,68]
  assert list(isel(r"resname so4")) == [64,65,66,67,68]
  assert list(isel(r"resname So4")) == [64,65,66,67,68]
  assert list(isel(r"resname S?4")) == [64,65,66,67,68]
  assert list(isel(r"resname pro and name cg")) == [22]
  assert list(isel(r"resname pro and (name cg or name ca)")) == [18,22]
  assert list(isel(r"resname pro AND (name cg or name ca)")) == [18,22]
  assert list(isel(r"resname pro and (name cg OR name ca)")) == [18,22]
  assert list(isel(r"resname pro AND (name cg OR name ca)")) == [18,22]
  assert list(isel(r"not resname pro and (name cg or name ca)")
              ) == [1,7,11,25,33]
  assert list(isel(r"chain h and name o*")) == [41,42,43,44,47,49,51,58]
  assert list(isel(r"(chain h or chain s) and name o[2-46]")) == [58,66,67,68]
  assert list(isel(r"resseq 188")) == [64,65,66,67,68]
  assert list(isel(r"resseq 188")) == [64,65,66,67,68]
  assert list(isel(r"resseq 1:1")) == [0,1,2,3,4,5,70,71,74,75,78,81]
  assert list(isel(r"resseq 2:2")) == range(6,17) + [72,73,77,79,80,82]
  assert list(isel(r"resseq 5:5")) == [76,83]
  assert list(isel(r"resseq 1:5")) == range(40)+range(70,84)
  assert list(isel(r"resseq 2:3")) == range(6,24)+[72,73,77,79,80,82]
  assert list(isel(r"resseq 188:188")) == [64,65,66,67,68]
  assert list(isel(r"resseq 200:200")) == [69]
  assert list(isel(r"resseq 188:200")) == [64,65,66,67,68,69]
  assert list(isel(r"resseq 9999:A002")) == [84,85]
  assert list(isel(r"resseq A002:A003")) == [85,86]
  assert list(isel(r"resseq :")) == range(88)
  assert list(isel(r"resseq :2 and name n*")) == [0,6,13,15,16]
  assert list(isel(r"resseq 2: and name cb")) == [10,21,28,36]
  assert list(isel(r"resseq 1:2 and name n*")) == [0,6,13,15,16]
  assert list(isel(r"resseq 2:4 and name cb")) == [10,21,28,36]
  assert list(isel(r"model 1 and name cb")) == [4,10,21,28,36]
  assert list(isel(r"model 2:3 and name o1*")) == [41,65]
  assert list(isel(r"icode j and name c?")) == [18,21,22,23]
  assert list(isel(r"resid 188")) == [64,65,66,67,68]
  assert list(isel(r"resid 3J")) == [17,18,19,20,21,22,23]
  assert list(isel(r"resid 1K")) == [78,81,87]
  assert list(isel(r"resid '   1K'")) == [78,81]
  assert list(isel(r"resid '  1K '")) == [87]
  assert list(isel(r"resi '  1K'")) == []
  assert list(isel(r"resid 1:2")) \
      == range(17) + [70,71,72,73,74,75,77,78,79,80,81,82]
  expected = range(6,17) + [72,73,77,78,79,80,81,82]
  assert list(isel(r"resid 1K:2")) == expected
  assert list(isel(r"resid '   1K:2'")) == expected
  assert list(isel(r"resid '  1K:2'")) == expected
  expected = range(6,40) + [72,73,76,77,78,79,80,81,82,83,87]
  assert list(isel(r"resi '  1K:  1K '")) == expected
  #
  expected = [7,18,25,33]
  assert list(isel(r"resseq 2:4 and name ca")) == expected
  assert list(isel(r"resseq 2 : 4 and name ca")) == expected
  assert list(isel(r"resseq 2: 4 and name ca")) == expected
  assert list(isel(r"resseq 2 :4 and name ca")) == expected
  expected = [1,7,18]
  assert list(isel(r"resseq :3 and name ca")) == expected
  assert list(isel(r"resseq : 3 and name ca")) == expected
  assert list(isel(r"name ca and resseq :3")) == expected
  assert list(isel(r"name ca and resseq : 3")) == expected
  expected = [18,25,33]
  assert list(isel(r"resseq 3: and name ca")) == expected
  assert list(isel(r"resseq 3 : and name ca")) == expected
  assert list(isel(r"name ca and resseq 3:")) == expected
  assert list(isel(r"name ca and resseq 3 :")) == expected
  assert list(isel(r"name ca and resseq 3")) == [18]
  expected = [1,7,18,25,33]
  assert list(isel(r"resseq : and name ca")) == expected
  assert list(isel(r"name ca and resseq :")) == expected
  #
  assert list(isel(r"segid wate")) == [69]
  assert list(isel(r"element o")) == [65,66,67,68]
  assert list(isel(r"charge 4+")) == [64]
  assert list(isel(r"anisou")) == [1, 3]
  assert list(isel(r"pepnames")) == range(40)
  assert list(isel(r"single_atom_residue")) == [
    69, 76, 77, 78, 81, 82, 83, 84, 85, 86, 87]
  #
  try: isel(r"resseq")
  except pdb.atom_selection.AtomSelectionError, e:
    assert str(e).find(
      "Missing argument for resseq.") >= 0
  else: raise Exception_expected
  try: isel(r"resseq 3:2")
  except pdb.atom_selection.AtomSelectionError, e:
    assert str(e).find(
      "range with first index > last index: resseq 3:2") >= 0
  else: raise Exception_expected
  try: isel(r"resid ' 1K :2'")
  except pdb.atom_selection.AtomSelectionError, e:
    assert str(e).find(
      "range with first index > last index: resid  1K :2") >= 0
  else: raise Exception_expected
  for s in ["altloc a and and name n",
            "altloc a and or name n",
            "altloc a or and name n",
            "altloc a or or name n",
            "and name n",
            "or name n",
            "not",
            "not not",
            "altloc a optional",
            "optional optional altloc a"]:
    try: isel(string=s)
    except pdb.atom_selection.AtomSelectionError, e:
      assert str(e).endswith("""\
Atom selection string leading to error:
  %s""" % s)
    else: raise Exception_expected
  #
  sel = sel_cache.get_labels(name=" CA ")
  assert len(sel) == 1
  assert list(sel[0]) == [1,7,18,25,33]
  atoms = hierarchy.atoms()
  for i_seq in sel[0]:
    assert atoms[i_seq].name == " CA "
  sel = sel_cache.get_labels(resseq="   5")
  assert len(sel) == 1
  assert list(sel[0]) == [76, 83]
  #
  link_records = [
    pdb.records.link(pdb_str=pdb_str)
      for pdb_str in """\
LINK         S   SO4 S 188                 O1  SO4 S 188
LINK         S   SO4 S 188                 O2  SO4 S 188
LINK         NZ  LYS A 680        1.260    C4A PLP D   1                LYS-PLP
""".splitlines()]
  expected_results = [
    [[64], [65]],
    [[64], [66]],
    [[], []]]
  for link_record,expected in zip(link_records, expected_results):
    assert [list(sel) for sel in sel_cache.link_iselections(link_record)] \
        == expected
  #
  hierarchy = pdb.input(source_info=None, lines=flex.split_lines("""\
CRYST1   21.937    4.866   23.477  90.00 107.08  90.00 P 1 21 1      2
ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      6  CA  ASN A   2      -6.522   2.038   2.831  1.00 14.10           C
ATOM     14  CA  Asn A   3      -3.193   1.904   4.589  1.00 11.74           C
ATOM     22  CA  GLN a   4       0.384   1.888   3.199  1.00 10.53           C
ATOM     31  CA  GLN a   5       3.270   2.361   5.640  1.00 11.39           C
ATOM     40  CA  ASN a   6       6.831   2.310   4.318  1.00 12.30           C
END
""")).construct_hierarchy()
  sel_cache = hierarchy.atom_selection_cache()
  isel = sel_cache.iselection
  assert list(isel("chain A")) == [0,1,2]
  assert list(isel("chain a")) == [3,4,5]
  assert list(isel("name ca")) == range(6)
  assert list(isel("resname asn")) == [1,2,5]
  assert list(isel("resname ASN")) == [1,2,5]
  assert list(isel("resname Asn")) == [1,2,5]
  hierarchy = pdb.input(source_info=None, lines=flex.split_lines("""\
CRYST1   21.937    4.866   23.477  90.00 107.08  90.00 P 1 21 1      2
ATOM      2  CA  GLY A   1      -9.052   4.207   4.651  1.00 16.57           C
ATOM      6  CA  ASN A   2      -6.522   2.038   2.831  1.00 14.10           C
ATOM     14  CA  ASN A   3      -3.193   1.904   4.589  1.00 11.74           C
ATOM     22  CA  GLN b   4       0.384   1.888   3.199  1.00 10.53           C
ATOM     31  CA  GLN b   5       3.270   2.361   5.640  1.00 11.39           C
ATOM     40  CA  ASN b   6       6.831   2.310   4.318  1.00 12.30           C
END
""")).construct_hierarchy()
  sel_cache = hierarchy.atom_selection_cache()
  isel = sel_cache.iselection
  assert list(isel("resname asn")) == [1,2,5]
  assert list(isel("resname ASN")) == [1,2,5]
  assert list(isel("resname Asn")) == [1,2,5]
  assert list(isel("chain A")) == [0,1,2]
  assert list(isel("chain a")) == []
  assert list(isel("chain B")) == []
  assert list(isel("chain b")) == [3,4,5]
  #
  hierarchy = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM      5  CA  SER     1       9.242  30.200  62.974  1.00 46.62
ATOM     11  CA  ARG     2      12.548  28.316  63.532  1.00 30.20
ATOM     21  N   NON     3J     13.947  29.997  64.680  1.00 22.94
ATOM     22  CA  NON     3J     14.902  31.100  64.827  1.00 20.19
ATOM     24  O   NON     3J     16.545  29.521  64.086  1.00 19.76
ATOM     28  N  AILE     4      16.953  31.648  63.512  1.00 15.29
ATOM     29  CA AILE     4      18.243  31.372  62.859  1.00 14.32
ATOM     30  C  AILE     4      19.233  32.112  63.743  1.00 13.54
ATOM     31  O  AILE     4      19.105  33.315  64.009  1.00 11.84
ATOM     41  N   PRO     5      13.947  29.997  64.680  1.00 22.94
ATOM     42  CA  PRO     5      14.902  31.100  64.827  1.00 20.19
ATOM     44  O   PRO     5      16.545  29.521  64.086  1.00 19.76
ATOM     45  CA  CA      6      16.545  29.521  64.086  1.00 19.76
""")).construct_hierarchy()
  sel_cache = hierarchy.atom_selection_cache()
  isel = sel_cache.iselection
  assert list(isel("pepnames")) == [0,1,5,6,7,8,9,10,11]
  #
  for s in ["peptide", "protein"]:
    try:
      isel(s)
    except pdb.atom_selection.AtomSelectionError, e:
      assert not show_diff(str(e), """\
Sorry: "%s" atom selection keyword not available:
  Please try using "pepnames" instead.
Atom selection string leading to error:
  %s""" % (s, s))
    else: raise Exception_expected
  #
  try:
    isel("chain A or (peptyde and name ca)")
  except pdb.atom_selection.AtomSelectionError, e:
    assert not show_diff(str(e), """\
RuntimeError: Atom selection syntax error at word "peptyde".
Atom selection string leading to error:
  chain A or (peptyde and name ca)""")
  else: raise Exception_expected
  #
  hierarchy = pdb.input(source_info=None, lines=flex.split_lines("""\
ATOM    459  CA  SER A  58
ATOM    463  OG ASER A  58
ATOM    464  OG BSER A  58
""")).construct_hierarchy()
  sel_cache = hierarchy.atom_selection_cache()
  assert sel_cache.iselection("single_atom_residue").size() == 0
  hierarchy = pdb.input(source_info=None, lines=flex.split_lines("""\
CRYST1   21.937    4.866   23.477  90.00 107.08  90.00 P 1 21 1      2
ATOM      2  CA  GLY A 666      -9.052   4.207   4.651  1.00 16.57           C
ATOM      6  CA  ASN A 777      -6.522   2.038   2.831  1.00 14.10           C
ATOM     14  CA  ASN A   1      -3.193   1.904   4.589  1.00 11.74           C
ATOM     22  CA  GLN A   2       0.384   1.888   3.199  1.00 10.53           C
ATOM     31  CA  GLN A   3       3.270   2.361   5.640  1.00 11.39           C
ATOM     40  CA  ASN A   4       6.831   2.310   4.318  1.00 12.30           C
TER
ATOM      6  CA  ASN B 777      -6.522   2.038   2.831  1.00 14.10           C
ATOM     14  CA  ASN B   1      -3.193   1.904   4.589  1.00 11.74           C
ATOM     22  CA  GLN B   2       0.384   1.888   3.199  1.00 10.53           C
END
""")).construct_hierarchy()
  sel_cache = hierarchy.atom_selection_cache()
  sele = sel_cache.iselection("resid 777 through 3 and chain A")
  assert (sele.size() == 4)
  sele = sel_cache.iselection("resid 777 through 3")
  assert (sele.size() == 7)
  try :
    sele = sel_cache.iselection("resid 777 through chain A")
  except pdb.atom_selection.AtomSelectionError, e :
    pass
  else : raise Exception_expected
  hierarchy = pdb.input(source_info=None, lines=flex.split_lines("""\
CRYST1   21.937    4.866   23.477  90.00 107.08  90.00 P 1 21 1      2
ATOM      6  CA  ASN B 777      -6.522   2.038   2.831  1.00 14.10           C
ATOM     14  CA  ASN B   1      -3.193   1.904   4.589  1.00 11.74           C
ATOM     22  CA  GLN B   2       0.384   1.888   3.199  1.00 10.53           C
TER
ATOM      2  CA  GLY A 666      -9.052   4.207   4.651  1.00 16.57           C
ATOM      6  CA  ASN A 777      -6.522   2.038   2.831  1.00 14.10           C
ATOM     14  CA  ASN A   1      -3.193   1.904   4.589  1.00 11.74           C
ATOM     22  CA  GLN A   2       0.384   1.888   3.199  1.00 10.53           C
ATOM     31  CA  GLN A   3       3.270   2.361   5.640  1.00 11.39           C
ATOM     40  CA  ASN A   4       6.831   2.310   4.318  1.00 12.30           C
TER
END
""")).construct_hierarchy()
  sel_cache = hierarchy.atom_selection_cache()
  sele = sel_cache.iselection("resid 2 through 4")
  assert (list(sele) == [2,6,7,8])
  sele = sel_cache.iselection("chain A and resid 2 through 4")
  assert (list(sele) == [6,7,8])

def run():
  exercise_selection()
  print "OK"

if (__name__ == "__main__"):
  run()
