from __future__ import absolute_import, division, print_function
import random
from scitbx.array_family import flex
import iotbx.map_model_manager
from mmtbx.building import ligands
import libtbx.load_env
import os

if (1):
  random.seed(5198425)
  flex.set_random_seed(5198425)

pdb_str = """
REMARK Idealized model of ATP from a library
HETATM    1  C1' ATP     1      25.100  19.996  20.709  1.00 30.00           C
HETATM    2  C2  ATP     1      26.706  23.954  19.757  1.00 30.00           C
HETATM    3  C2' ATP     1      26.458  19.303  20.846  1.00 30.00           C
HETATM    4  C3' ATP     1      26.049  17.830  20.809  1.00 30.00           C
HETATM    5  C4  ATP     1      25.589  22.474  20.949  1.00 30.00           C
HETATM    6  C4' ATP     1      24.723  17.840  21.563  1.00 30.00           C
HETATM    7  C5  ATP     1      25.243  23.419  21.902  1.00 30.00           C
HETATM    8  C5' ATP     1      24.831  17.625  23.052  1.00 30.00           C
HETATM    9  C6  ATP     1      25.709  24.731  21.693  1.00 30.00           C
HETATM   10  C8  ATP     1      24.377  21.591  22.542  1.00 30.00           C
HETATM   11  N1  ATP     1      26.452  24.967  20.589  1.00 30.00           N
HETATM   12  N3  ATP     1      26.321  22.679  19.849  1.00 30.00           N
HETATM   13  N6  ATP     1      25.459  25.751  22.522  1.00 30.00           N
HETATM   14  N7  ATP     1      24.475  22.846  22.905  1.00 30.00           N
HETATM   15  N9  ATP     1      25.025  21.294  21.370  1.00 30.00           N
HETATM   16  O1A ATP     1      22.686  19.051  25.764  1.00 30.00           O
HETATM   17  O1B ATP     1      20.672  16.689  23.539  1.00 30.00           O
HETATM   18  O1G ATP     1      20.680  12.823  22.913  1.00 30.00           O
HETATM   19  O2' ATP     1      27.306  19.645  19.755  1.00 30.00           O
HETATM   20  O2A ATP     1      24.704  17.520  25.870  1.00 30.00           O
HETATM   21  O2B ATP     1      20.013  16.346  25.973  1.00 30.00           O
HETATM   22  O2G ATP     1      20.809  12.141  25.347  1.00 30.00           O
HETATM   23  O3' ATP     1      25.876  17.362  19.476  1.00 30.00           O
HETATM   24  O3A ATP     1      22.442  16.557  25.361  1.00 30.00           O
HETATM   25  O3B ATP     1      21.206  14.511  24.739  1.00 30.00           O
HETATM   26  O3G ATP     1      18.918  13.574  24.568  1.00 30.00           O
HETATM   27  O4' ATP     1      24.141  19.137  21.303  1.00 30.00           O
HETATM   28  O5' ATP     1      23.571  17.969  23.664  1.00 30.00           O
HETATM   29  PA  ATP     1      23.358  17.852  25.225  1.00 30.00           P
HETATM   30  PB  ATP     1      21.011  16.087  24.844  1.00 30.00           P
HETATM   31  PG  ATP     1      20.385  13.201  24.324  1.00 30.00           P
TER
"""

def run(prefix="tst_00_mmtbx_building_ligands"):
  # Ligand file
  with open("%s.pdb"%prefix,"w") as fo:
    fo.write(pdb_str)
  # Read map and model
  from iotbx.data_manager import DataManager
  dm = DataManager()
  map_file = libtbx.env.find_in_repositories(
    relative_path="mmtbx/building/ligands/tst_00_mmtbx_building_ligands.map",
    test=os.path.isfile)
  mm = dm.get_real_map(map_file)
  model = dm.get_model("%s.pdb"%prefix)
  model.set_crystal_symmetry(mm.crystal_symmetry())
  model.process(make_restraints=True)
  # Create map_model_manager
  mmm = iotbx.map_model_manager.map_model_manager(map_manager=mm, model=model)
  # Build ligand
  o = ligands.lifi.run(map_model_manager = mmm, d_min = 2.5)

if(__name__ == "__main__"):
  run()
