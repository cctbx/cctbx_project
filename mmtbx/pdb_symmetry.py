
import libtbx.load_env
from libtbx import easy_pickle
from libtbx import group_args
from math import sqrt
import os

def parse_database (file_name) :
  from cctbx import crystal
  from cctbx import sgtbx
  from cctbx import uctbx
  db = []
  lines = open(file_name).readlines()
  process = False
  for line in lines :
    line = line.strip()
    if (process) :
      fields = line.split()
      pdb_id = fields[0]
      try :
        uc = uctbx.unit_cell([ float(x.replace(",","")) for x in fields[2:8] ])
      except RuntimeError, e :
        print "Unit cell error:"
        print line
        continue
      try :
        sg = sgtbx.space_group_info(" ".join(fields[8:-1]))
      except RuntimeError, e :
        print "Unrecognized space group:"
        print line
        continue
      try :
        symm = crystal.symmetry(unit_cell=uc, space_group_info=sg)
      except AssertionError :
        print "Incompatible unit cell parameters:"
        print line
        continue
      niggli_symm = symm.niggli_cell()
      db.append(group_args(
        pdb_id=pdb_id,
        crystal_symmetry=symm,
        niggli_cell=niggli_symm))
    elif (line.startswith("-----")) :
      process = True
  return db

db_ = None
def load_db (save_pickle=False) :
  global db_
  if (db_ is not None) :
    return db_
  pkl_file_name = libtbx.env.find_in_repositories(
    relative_path="chem_data/pdb/crystal.pkl",
    test=os.path.isfile)
  txt_file_name = libtbx.env.find_in_repositories(
    relative_path="chem_data/pdb/crystal.idx",
    test=os.path.isfile)
  if (pkl_file_name is None) or (save_pickle) :
    assert (txt_file_name is not None)
    db_ = parse_database(txt_file_name)
    if (save_pickle) :
      easy_pickle.dump(txt_file_name[:-4] + ".pkl", db_)
  else :
    db_ = easy_pickle.load(pkl_file_name)
  return db_

def rms_difference (params1, params2) :
  assert (len(params1) == len(params2) != 0)
  sum = 0
  n = 0
  for x1, x2 in zip(params1, params2) :
    sum += (x1-x2)**2
  return sqrt(sum / len(params1))

def symmetry_search (
    symmetry_db,
    crystal_symmetry,
    max_rmsd=None) :
  scores = []
  niggli_cell = crystal_symmetry.niggli_cell()
  input_volume = niggli_cell.unit_cell().volume()
  niggli_uc = niggli_cell.unit_cell().parameters()
  for entry in symmetry_db :
    entry_volume = entry.niggli_cell.unit_cell().volume()
    entry_uc = entry.niggli_cell.unit_cell().parameters()
    # XXX suggested by James Holton - sum separate RMSDs for edge lengths and
    # angles
    rmsd = rms_difference(niggli_uc[0:3], entry_uc[0:3]) + \
           rms_difference(niggli_uc[3:6], entry_uc[3:6])
    if (max_rmsd is not None) and (rmsd > max_rmsd) :
      continue
    scores.append(group_args(entry=entry,
      rmsd=rmsd,
      volume_ratio=entry_volume/input_volume))
  scores.sort(lambda x,y: cmp(x.rmsd, y.rmsd))
  return scores

if (__name__ == "__main__") :
  load_db(save_pickle=True)
