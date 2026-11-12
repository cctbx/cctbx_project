from __future__ import absolute_import, division, print_function
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
from mmtbx.conformation_dependent_library.hpdl_database import get_hpdl_database

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

# HPDL key; angles: 7 heavy-atom + 8 ring-H, less the 2 of an absent HD1/HE2
protonations = {
  "his_double" : ("ND1 and NE2 protonated", 15),
  "his_nd1"    : ("Only ND1 protonated", 13),
  "his_ne2"    : ("Only NE2 protonated", 13),
  }

def check_ideals(lines, restraints, n_angles_expected):
  """Process one His with hpdl=True and compare the bond and angle ideals of
  the geometry restraints manager with the HPDL table for its protonation."""
  model = mmtbx.model.manager(
    model_input=iotbx.pdb.input(source_info=None, lines=lines), log=null_out())
  params = mmtbx.model.manager.get_default_pdb_interpretation_params()
  params.pdb_interpretation.restraints_library.hpdl = True
  model.process(make_restraints=True, pdb_interpretation_params=params)
  grm = model.get_restraints_manager().geometry
  names = [a.name.strip() for a in model.get_hierarchy().atoms()]
  bonds = {}
  for p in grm.get_all_bond_proxies(sites_cart=model.get_sites_cart())[0]:
    bonds[frozenset(names[i] for i in p.i_seqs)] = p.distance_ideal
  angles = {}
  for p in grm.angle_proxies:
    i, j, k = [names[i_seq] for i_seq in p.i_seqs]
    angles[(j, frozenset((i, k)))] = p.angle_ideal
  n_bonds = 0
  n_angles = 0
  for restraint, values in restraints.items():
    if len(restraint) == 2:
      ideal = bonds.get(frozenset(restraint))
      if ideal is None: continue
      assert ideal == values[0], (restraint, ideal, values)
      n_bonds += 1
    elif len(restraint) == 3:
      i, j, k = restraint
      ideal = angles.get((j, frozenset((i, k))))
      if ideal is None: continue
      assert ideal == values[0], (restraint, ideal, values)
      n_angles += 1
  assert n_bonds == 6, n_bonds
  assert n_angles == n_angles_expected, (n_angles, n_angles_expected)

def run():
  hpdl_database = get_hpdl_database()
  for name, lines in pdbs.items():
    print(name)
    protonation, n_angles = protonations[name]
    check_ideals(lines, hpdl_database[protonation], n_angles)
  print("OK")

if __name__=="__main__":
  run()
