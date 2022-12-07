from __future__ import absolute_import, division, print_function
import iotbx.pdb

# All H ------------------------------------------------------------------------

pdb_str1_H = """
CRYST1   13.022   15.843   16.570  90.00  90.00  90.00 P 1
ATOM      1  N   HIS             8.022   8.797   9.927  1.00  0.00           N
ATOM      2  CA  HIS             7.020   9.705   9.381  1.00  0.00           C
ATOM      3  C   HIS             6.734  10.843  10.354  1.00  0.00           C
ATOM      4  O   HIS             6.781  10.659  11.570  1.00  0.00           O
ATOM      5  CB  HIS             5.729   8.948   9.062  1.00  0.00           C
ATOM      6  CG  HIS             5.905   7.850   8.060  1.00  0.00           C
ATOM      7  ND1 HIS             5.768   8.050   6.704  1.00  0.00           N
ATOM      8  CD2 HIS             6.207   6.539   8.218  1.00  0.00           C
ATOM      9  CE1 HIS             5.978   6.911   6.069  1.00  0.00           C
ATOM     10  NE2 HIS             6.246   5.978   6.965  1.00  0.00           N
ATOM     11  H   HIS             7.644   8.098  10.566  1.00  0.00           H
ATOM     12  HA  HIS             7.398  10.137   8.454  1.00  0.00           H
ATOM     13  HB2 HIS             5.000   9.652   8.660  1.00  0.00           H
ATOM     14  HB3 HIS             5.346   8.504   9.981  1.00  0.00           H
ATOM     15  HD1 HIS             5.540   8.939   6.259  1.00  0.00           H
ATOM     16  HD2 HIS             6.384   6.030   9.154  1.00  0.00           H
ATOM     17  HE1 HIS             5.937   6.766   5.000  1.00  0.00           H
ATOM     18  HE2 HIS             6.450   5.000   6.760  1.00  0.00           H
TER
END
"""

pdb_str2_H = """
CRYST1   13.022   15.843   16.570  90.00  90.00  90.00 P 1
ATOM      1  N   HIS             8.022   8.797   9.927  1.00  0.00           N
ATOM      2  CA  HIS             7.020   9.705   9.381  1.00  0.00           C
ATOM      3  C   HIS             6.734  10.843  10.354  1.00  0.00           C
ATOM      4  O   HIS             6.781  10.659  11.570  1.00  0.00           O
ATOM      5  CB  HIS             5.729   8.948   9.062  1.00  0.00           C
ATOM      6  CG  HIS             5.905   7.850   8.060  1.00  0.00           C
ATOM      7  ND1 HIS             5.768   8.050   6.704  1.00  0.00           N
ATOM      8  CD2 HIS             6.207   6.539   8.218  1.00  0.00           C
ATOM      9  CE1 HIS             5.978   6.911   6.069  1.00  0.00           C
ATOM     10  NE2 HIS             6.246   5.978   6.965  1.00  0.00           N
ATOM     11  H   HIS             7.644   8.098  10.566  1.00  0.00           H
ATOM     12  HA  HIS             7.398  10.137   8.454  1.00  0.00           H
ATOM     13  HB2 HIS             5.000   9.652   8.660  1.00  0.00           H
ATOM     14  HB3 HIS             5.346   8.504   9.981  1.00  0.00           H
ATOM     15  HD1 HIS             5.540   8.939   6.259  1.00  0.00           H
ATOM     16  HD2 HIS             6.384   6.030   9.154  1.00  0.00           H
ATOM     17  HE1 HIS             5.937   6.766   5.000  1.00  0.00           H
TER
END
"""

pdb_str3_H = """
CRYST1   13.022   15.843   16.570  90.00  90.00  90.00 P 1
ATOM      1  N   HIS             8.022   8.797   9.927  1.00  0.00           N
ATOM      2  CA  HIS             7.020   9.705   9.381  1.00  0.00           C
ATOM      3  C   HIS             6.734  10.843  10.354  1.00  0.00           C
ATOM      4  O   HIS             6.781  10.659  11.570  1.00  0.00           O
ATOM      5  CB  HIS             5.729   8.948   9.062  1.00  0.00           C
ATOM      6  CG  HIS             5.905   7.850   8.060  1.00  0.00           C
ATOM      7  ND1 HIS             5.768   8.050   6.704  1.00  0.00           N
ATOM      8  CD2 HIS             6.207   6.539   8.218  1.00  0.00           C
ATOM      9  CE1 HIS             5.978   6.911   6.069  1.00  0.00           C
ATOM     10  NE2 HIS             6.246   5.978   6.965  1.00  0.00           N
ATOM     11  H   HIS             7.644   8.098  10.566  1.00  0.00           H
ATOM     12  HA  HIS             7.398  10.137   8.454  1.00  0.00           H
ATOM     13  HB2 HIS             5.000   9.652   8.660  1.00  0.00           H
ATOM     14  HB3 HIS             5.346   8.504   9.981  1.00  0.00           H
ATOM     16  HD2 HIS             6.384   6.030   9.154  1.00  0.00           H
ATOM     17  HE1 HIS             5.937   6.766   5.000  1.00  0.00           H
ATOM     18  HE2 HIS             6.450   5.000   6.760  1.00  0.00           H
TER
END
"""

# All D ------------------------------------------------------------------------

pdb_str1_D = """
CRYST1   13.022   15.843   16.570  90.00  90.00  90.00 P 1
ATOM      1  N   HIS             8.022   8.797   9.927  1.00  0.00           N
ATOM      2  CA  HIS             7.020   9.705   9.381  1.00  0.00           C
ATOM      3  C   HIS             6.734  10.843  10.354  1.00  0.00           C
ATOM      4  O   HIS             6.781  10.659  11.570  1.00  0.00           O
ATOM      5  CB  HIS             5.729   8.948   9.062  1.00  0.00           C
ATOM      6  CG  HIS             5.905   7.850   8.060  1.00  0.00           C
ATOM      7  ND1 HIS             5.768   8.050   6.704  1.00  0.00           N
ATOM      8  CD2 HIS             6.207   6.539   8.218  1.00  0.00           C
ATOM      9  CE1 HIS             5.978   6.911   6.069  1.00  0.00           C
ATOM     10  NE2 HIS             6.246   5.978   6.965  1.00  0.00           N
ATOM     11  D   HIS             7.644   8.098  10.566  1.00  0.00           D
ATOM     12  DA  HIS             7.398  10.137   8.454  1.00  0.00           D
ATOM     13  DB2 HIS             5.000   9.652   8.660  1.00  0.00           D
ATOM     14  DB3 HIS             5.346   8.504   9.981  1.00  0.00           D
ATOM     15  DD1 HIS             5.540   8.939   6.259  1.00  0.00           D
ATOM     16  DD2 HIS             6.384   6.030   9.154  1.00  0.00           D
ATOM     17  DE1 HIS             5.937   6.766   5.000  1.00  0.00           D
ATOM     18  DE2 HIS             6.450   5.000   6.760  1.00  0.00           D
TER
END
"""

pdb_str2_D = """
CRYST1   13.022   15.843   16.570  90.00  90.00  90.00 P 1
ATOM      1  N   HIS             8.022   8.797   9.927  1.00  0.00           N
ATOM      2  CA  HIS             7.020   9.705   9.381  1.00  0.00           C
ATOM      3  C   HIS             6.734  10.843  10.354  1.00  0.00           C
ATOM      4  O   HIS             6.781  10.659  11.570  1.00  0.00           O
ATOM      5  CB  HIS             5.729   8.948   9.062  1.00  0.00           C
ATOM      6  CG  HIS             5.905   7.850   8.060  1.00  0.00           C
ATOM      7  ND1 HIS             5.768   8.050   6.704  1.00  0.00           N
ATOM      8  CD2 HIS             6.207   6.539   8.218  1.00  0.00           C
ATOM      9  CE1 HIS             5.978   6.911   6.069  1.00  0.00           C
ATOM     10  NE2 HIS             6.246   5.978   6.965  1.00  0.00           N
ATOM     11  D   HIS             7.644   8.098  10.566  1.00  0.00           D
ATOM     12  DA  HIS             7.398  10.137   8.454  1.00  0.00           D
ATOM     13  DB2 HIS             5.000   9.652   8.660  1.00  0.00           D
ATOM     14  DB3 HIS             5.346   8.504   9.981  1.00  0.00           D
ATOM     15  DD1 HIS             5.540   8.939   6.259  1.00  0.00           D
ATOM     16  DD2 HIS             6.384   6.030   9.154  1.00  0.00           D
ATOM     17  DE1 HIS             5.937   6.766   5.000  1.00  0.00           D
TER
END
"""

pdb_str3_D = """
CRYST1   13.022   15.843   16.570  90.00  90.00  90.00 P 1
ATOM      1  N   HIS             8.022   8.797   9.927  1.00  0.00           N
ATOM      2  CA  HIS             7.020   9.705   9.381  1.00  0.00           C
ATOM      3  C   HIS             6.734  10.843  10.354  1.00  0.00           C
ATOM      4  O   HIS             6.781  10.659  11.570  1.00  0.00           O
ATOM      5  CB  HIS             5.729   8.948   9.062  1.00  0.00           C
ATOM      6  CG  HIS             5.905   7.850   8.060  1.00  0.00           C
ATOM      7  ND1 HIS             5.768   8.050   6.704  1.00  0.00           N
ATOM      8  CD2 HIS             6.207   6.539   8.218  1.00  0.00           C
ATOM      9  CE1 HIS             5.978   6.911   6.069  1.00  0.00           C
ATOM     10  NE2 HIS             6.246   5.978   6.965  1.00  0.00           N
ATOM     11  D   HIS             7.644   8.098  10.566  1.00  0.00           D
ATOM     12  DA  HIS             7.398  10.137   8.454  1.00  0.00           D
ATOM     13  DB2 HIS             5.000   9.652   8.660  1.00  0.00           D
ATOM     14  DB3 HIS             5.346   8.504   9.981  1.00  0.00           D
ATOM     16  DD2 HIS             6.384   6.030   9.154  1.00  0.00           D
ATOM     17  DE1 HIS             5.937   6.766   5.000  1.00  0.00           D
ATOM     18  DE2 HIS             6.450   5.000   6.760  1.00  0.00           D
TER
END
"""

# All H/D ----------------------------------------------------------------------

pdb_str1_HD = """
CRYST1   13.022   15.843   16.570  90.00  90.00  90.00 P 1
ATOM      1  N   HIS             8.022   8.797   9.927  1.00  0.00           N
ATOM      2  CA  HIS             7.020   9.705   9.381  1.00  0.00           C
ATOM      3  C   HIS             6.734  10.843  10.354  1.00  0.00           C
ATOM      4  O   HIS             6.781  10.659  11.570  1.00  0.00           O
ATOM      5  CB  HIS             5.729   8.948   9.062  1.00  0.00           C
ATOM      6  CG  HIS             5.905   7.850   8.060  1.00  0.00           C
ATOM      7  ND1 HIS             5.768   8.050   6.704  1.00  0.00           N
ATOM      8  CD2 HIS             6.207   6.539   8.218  1.00  0.00           C
ATOM      9  CE1 HIS             5.978   6.911   6.069  1.00  0.00           C
ATOM     10  NE2 HIS             6.246   5.978   6.965  1.00  0.00           N
ATOM     11  H  AHIS             7.644   8.098  10.566  0.50  0.00           H
ATOM     11  D  BHIS             7.644   8.098  10.566  0.50  0.00           D
ATOM     12  HA  HIS             7.398  10.137   8.454  1.00  0.00           H
ATOM     13  HB2 HIS             5.000   9.652   8.660  1.00  0.00           H
ATOM     14  HB3 HIS             5.346   8.504   9.981  1.00  0.00           H
ATOM     15  HD1AHIS             5.540   8.939   6.259  0.50  0.00           H
ATOM     15  DD1BHIS             5.540   8.939   6.259  0.50  0.00           D
ATOM     16  HD2 HIS             6.384   6.030   9.154  1.00  0.00           H
ATOM     17  HE1 HIS             5.937   6.766   5.000  1.00  0.00           H
ATOM     18  HE2AHIS             6.450   5.000   6.760  0.50  0.00           H
ATOM     18  DE2BHIS             6.450   5.000   6.760  0.50  0.00           D
TER
END
"""

pdb_str2_HD = """
CRYST1   13.022   15.843   16.570  90.00  90.00  90.00 P 1
ATOM      1  N   HIS             8.022   8.797   9.927  1.00  0.00           N
ATOM      2  CA  HIS             7.020   9.705   9.381  1.00  0.00           C
ATOM      3  C   HIS             6.734  10.843  10.354  1.00  0.00           C
ATOM      4  O   HIS             6.781  10.659  11.570  1.00  0.00           O
ATOM      5  CB  HIS             5.729   8.948   9.062  1.00  0.00           C
ATOM      6  CG  HIS             5.905   7.850   8.060  1.00  0.00           C
ATOM      7  ND1 HIS             5.768   8.050   6.704  1.00  0.00           N
ATOM      8  CD2 HIS             6.207   6.539   8.218  1.00  0.00           C
ATOM      9  CE1 HIS             5.978   6.911   6.069  1.00  0.00           C
ATOM     10  NE2 HIS             6.246   5.978   6.965  1.00  0.00           N
ATOM     11  H  AHIS             7.644   8.098  10.566  0.50  0.00           H
ATOM     11  D  BHIS             7.644   8.098  10.566  0.50  0.00           D
ATOM     12  HA  HIS             7.398  10.137   8.454  1.00  0.00           H
ATOM     13  HB2 HIS             5.000   9.652   8.660  1.00  0.00           H
ATOM     14  HB3 HIS             5.346   8.504   9.981  1.00  0.00           H
ATOM     15  HD1AHIS             5.540   8.939   6.259  0.50  0.00           H
ATOM     15  DD1BHIS             5.540   8.939   6.259  0.50  0.00           D
ATOM     16  HD2 HIS             6.384   6.030   9.154  1.00  0.00           H
ATOM     17  HE1 HIS             5.937   6.766   5.000  1.00  0.00           H
TER
END
"""

pdb_str3_HD = """
CRYST1   13.022   15.843   16.570  90.00  90.00  90.00 P 1
ATOM      1  N   HIS             8.022   8.797   9.927  1.00  0.00           N
ATOM      2  CA  HIS             7.020   9.705   9.381  1.00  0.00           C
ATOM      3  C   HIS             6.734  10.843  10.354  1.00  0.00           C
ATOM      4  O   HIS             6.781  10.659  11.570  1.00  0.00           O
ATOM      5  CB  HIS             5.729   8.948   9.062  1.00  0.00           C
ATOM      6  CG  HIS             5.905   7.850   8.060  1.00  0.00           C
ATOM      7  ND1 HIS             5.768   8.050   6.704  1.00  0.00           N
ATOM      8  CD2 HIS             6.207   6.539   8.218  1.00  0.00           C
ATOM      9  CE1 HIS             5.978   6.911   6.069  1.00  0.00           C
ATOM     10  NE2 HIS             6.246   5.978   6.965  1.00  0.00           N
ATOM     11  H  AHIS             7.644   8.098  10.566  0.50  0.00           H
ATOM     11  D  BHIS             7.644   8.098  10.566  0.50  0.00           D
ATOM     12  HA  HIS             7.398  10.137   8.454  1.00  0.00           H
ATOM     13  HB2 HIS             5.000   9.652   8.660  1.00  0.00           H
ATOM     14  HB3 HIS             5.346   8.504   9.981  1.00  0.00           H
ATOM     16  HD2 HIS             6.384   6.030   9.154  1.00  0.00           H
ATOM     17  HE1 HIS             5.937   6.766   5.000  1.00  0.00           H
ATOM     18  HE2AHIS             6.450   5.000   6.760  0.50  0.00           H
ATOM     18  DE2BHIS             6.450   5.000   6.760  0.50  0.00           D
TER
END
"""

class load(object):
  def __init__(self):
    self.r = {}
    all_his = [
      ["h1",  pdb_str1_H],
      ["h2",  pdb_str2_H],
      ["h3",  pdb_str3_H],
      ["d1",  pdb_str1_D],
      ["d2",  pdb_str2_D],
      ["d3",  pdb_str3_D],
      ["hd1", pdb_str1_HD],
      ["hd2", pdb_str2_HD],
      ["hd3", pdb_str3_HD]
    ]
    for his in all_his:
      self.r[his[0]] = iotbx.pdb.input(source_info = None, lines = his[1]
        ).construct_hierarchy().only_residue_group().detached_copy()

  def main_three(self):
    result = []
    for ps in ["h1", "h2", "h3"]:
      result.append(self.r[ps])
    return result
