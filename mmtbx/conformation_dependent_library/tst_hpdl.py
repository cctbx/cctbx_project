from __future__ import absolute_import, division, print_function
import sys

from mmtbx.conformation_dependent_library.hpdl_database import get_hpdl_database
from libtbx import easy_run
try:
  from elbow.formats import refine_geo_parser
except ImportError as e:
  refine_geo_parser = None

pdbs = {
  "his_double" : """
HETATM    1  N   HIS A   1      -1.844  -0.558   1.760  1.00 20.00      A    N
HETATM    2  CA  HIS A   1      -0.434  -0.568   1.417  1.00 20.00      A    C
HETATM    3  C   HIS A   1       0.399  -0.690   2.691  1.00 20.00      A    C
HETATM    4  O   HIS A   1       0.256   0.148   3.620  1.00 20.00      A    O
HETATM    5  CB  HIS A   1      -0.076   0.731   0.698  1.00 20.00      A    C
HETATM    6  CG  HIS A   1       0.278   0.426  -0.756  1.00 20.00      A    C
HETATM    7  ND1 HIS A   1       1.450   0.029  -1.199  1.00 20.00      A    N+1
HETATM    8  CD2 HIS A   1      -0.547   0.574  -1.828  1.00 20.00      A    C
HETATM    9  CE1 HIS A   1       1.375  -0.076  -2.534  1.00 20.00      A    C
HETATM   10  NE2 HIS A   1       0.133   0.264  -2.909  1.00 20.00      A    N
HETATM   11  OXT HIS A   1       1.247  -1.615   2.803  1.00 20.00      A    O-1
HETATM   12  H   HIS A   1      -2.107   0.359   2.068  1.00 20.00      A    H
HETATM   13  H2  HIS A   1      -2.386  -0.808   0.954  1.00 20.00      A    H
HETATM   14  HA  HIS A   1      -0.226  -1.410   0.768  1.00 20.00      A    H
HETATM   15  HB2 HIS A   1       0.774   1.193   1.188  1.00 20.00      A    H
HETATM   16  HB3 HIS A   1      -0.923   1.407   0.730  1.00 20.00      A    H
HETATM   17  HD1 HIS A   1       2.256  -0.162  -0.638  1.00 20.00      A    H
HETATM   18  HD2 HIS A   1      -1.591   0.861  -1.793  1.00 20.00      A    H
HETATM   19  HE1 HIS A   1       2.177  -0.383  -3.194  1.00 20.00      A    H
HETATM   20  HE2 HIS A   1      -0.215   0.280  -3.847  1.00 20.00      A    H
""",
  "his_nd1" : """
HETATM    1  N   HIS A   1      -1.844  -0.558   1.760  1.00 20.00      A    N
HETATM    2  CA  HIS A   1      -0.434  -0.568   1.417  1.00 20.00      A    C
HETATM    3  C   HIS A   1       0.399  -0.690   2.691  1.00 20.00      A    C
HETATM    4  O   HIS A   1       0.256   0.148   3.620  1.00 20.00      A    O
HETATM    5  CB  HIS A   1      -0.076   0.731   0.698  1.00 20.00      A    C
HETATM    6  CG  HIS A   1       0.278   0.426  -0.756  1.00 20.00      A    C
HETATM    7  ND1 HIS A   1       1.450   0.029  -1.199  1.00 20.00      A    N+1
HETATM    8  CD2 HIS A   1      -0.547   0.574  -1.828  1.00 20.00      A    C
HETATM    9  CE1 HIS A   1       1.375  -0.076  -2.534  1.00 20.00      A    C
HETATM   10  NE2 HIS A   1       0.133   0.264  -2.909  1.00 20.00      A    N
HETATM   11  OXT HIS A   1       1.247  -1.615   2.803  1.00 20.00      A    O-1
HETATM   12  H   HIS A   1      -2.107   0.359   2.068  1.00 20.00      A    H
HETATM   13  H2  HIS A   1      -2.386  -0.808   0.954  1.00 20.00      A    H
HETATM   14  HA  HIS A   1      -0.226  -1.410   0.768  1.00 20.00      A    H
HETATM   15  HB2 HIS A   1       0.774   1.193   1.188  1.00 20.00      A    H
HETATM   16  HB3 HIS A   1      -0.923   1.407   0.730  1.00 20.00      A    H
HETATM   17  HD1 HIS A   1       2.256  -0.162  -0.638  1.00 20.00      A    H
HETATM   18  HD2 HIS A   1      -1.591   0.861  -1.793  1.00 20.00      A    H
HETATM   19  HE1 HIS A   1       2.177  -0.383  -3.194  1.00 20.00      A    H
""",
  "his_ne2" : """
HETATM    1  N   HIS A   1      -1.844  -0.558   1.760  1.00 20.00      A    N
HETATM    2  CA  HIS A   1      -0.434  -0.568   1.417  1.00 20.00      A    C
HETATM    3  C   HIS A   1       0.399  -0.690   2.691  1.00 20.00      A    C
HETATM    4  O   HIS A   1       0.256   0.148   3.620  1.00 20.00      A    O
HETATM    5  CB  HIS A   1      -0.076   0.731   0.698  1.00 20.00      A    C
HETATM    6  CG  HIS A   1       0.278   0.426  -0.756  1.00 20.00      A    C
HETATM    7  ND1 HIS A   1       1.450   0.029  -1.199  1.00 20.00      A    N+1
HETATM    8  CD2 HIS A   1      -0.547   0.574  -1.828  1.00 20.00      A    C
HETATM    9  CE1 HIS A   1       1.375  -0.076  -2.534  1.00 20.00      A    C
HETATM   10  NE2 HIS A   1       0.133   0.264  -2.909  1.00 20.00      A    N
HETATM   11  OXT HIS A   1       1.247  -1.615   2.803  1.00 20.00      A    O-1
HETATM   12  H   HIS A   1      -2.107   0.359   2.068  1.00 20.00      A    H
HETATM   13  H2  HIS A   1      -2.386  -0.808   0.954  1.00 20.00      A    H
HETATM   14  HA  HIS A   1      -0.226  -1.410   0.768  1.00 20.00      A    H
HETATM   15  HB2 HIS A   1       0.774   1.193   1.188  1.00 20.00      A    H
HETATM   16  HB3 HIS A   1      -0.923   1.407   0.730  1.00 20.00      A    H
HETATM   18  HD2 HIS A   1      -1.591   0.861  -1.793  1.00 20.00      A    H
HETATM   19  HE1 HIS A   1       2.177  -0.383  -3.194  1.00 20.00      A    H
HETATM   20  HE2 HIS A   1      -0.215   0.280  -3.847  1.00 20.00      A    H
""",
  }

def check_ideals(geo_filename, restraints):
  rc = refine_geo_parser.run(geo_filename,
                             selection1="HIS",
                             return_internals=True)
  bonds, angles, dihedrals, chirals, planes = rc
  bond_remove=[]
  angle_remove=[]
  for restraint, values in restraints.items():
    if len(restraint)==2:
      for i, bond in enumerate(bonds):
        if ((bond[0].find(" %s " % restraint[0])>-1 and
             bond[1].find(" %s " % restraint[1])>-1) or
            (bond[0].find(" %s " % restraint[1])>-1 and
             bond[1].find(" %s " % restraint[0])>-1)):
          assert float(bond[2][0])==values[0]
          bond_remove.append(i)
          break
    elif len(restraint)==3:
      for i, angle in enumerate(angles):
        if angle[1].find(" %s " % restraint[1])==-1: continue
        if ((angle[0].find(" %s " % restraint[0])>-1 and
             angle[2].find(" %s " % restraint[2])>-1) or
            (angle[0].find(" %s " % restraint[2])>-1 and
             angle[2].find(" %s " % restraint[0])>-1)):
          assert float(angle[3][0])==values[0], "mismatch %s %s to %s" % (
            restraint,
            values,
            angle,
            )
          angle_remove.append(i)
          break
  assert len(bond_remove)==6
  assert len(angle_remove)>=13
  bond_remove.sort()
  bond_remove.reverse()
  for i in bond_remove: del bonds[i]
  print("OK")

def run(filename=None):
  assert filename is None
  hpdl_database = get_hpdl_database()
  if 0:
    print("HPDL")
    print(hpdl_database)
    for aa in sorted(hpdl_database):
      print("  %s" % aa)
      for key, value in sorted(hpdl_database[aa].items()):
        print("    %s : %s" % (key, value))
    assert 0
  #
  for pdb, lines in pdbs.items():
    print(pdb)
    filename="hpdl_%s.pdb" % pdb
    f=open(filename, "w")
    f.write(lines)
    f.close()
    cmd = "phenix.pdb_interpretation %s write_geo=1 hpdl=%s" % (filename, True)
    print(cmd)
    assert not easy_run.call(cmd)
    if filename.find("his_double")>-1: key="ND1 and NE2 protonated"
    elif filename.find("his_nd1")>-1: key="Only ND1 protonated"
    elif filename.find("his_ne2")>-1: key="Only NE2 protonated"
    if refine_geo_parser:
      check_ideals("%s.geo" % filename, hpdl_database[key])
  print("OK")

if __name__=="__main__":
  run(*tuple(sys.argv[1:]))
