from __future__ import division
from __future__ import print_function
import os, sys
from cctbx.array_family import flex
from libtbx import easy_pickle
import libtbx.load_env
from libtbx.test_utils import approx_equal

def load_dict():
  file = libtbx.env.find_in_repositories(relative_path=
    "chem_data/polygon_data/all_mvd.pickle", test=os.path.isfile)
  database_dict = easy_pickle.load(file)
  # Python 3 pickle fix
  # =========================================================================
  if sys.version_info.major == 3:
    database_dict = easy_pickle.fix_py2_pickle(database_dict)
  # =========================================================================

  return database_dict

def exercise_00():
  database_dict = load_dict()
  #
  i_1l2h = list(database_dict["pdb_code"]).index("1l2h")
  found = 0
  for key in database_dict.keys():
    #print key, ": ",database_dict[key][i_1l2h]
    #
    value = database_dict[key][i_1l2h]
    if(key == "unit_cell"):
      assert "(53.89, 53.89, 77.36, 90, 90, 90)" == value.strip()
      found += 1
    if(key == "high_resolution"):
      assert approx_equal(1.54, float(value))
      found += 1
    if(key == "pdb_header_tls"):
      assert "false (n_groups: 0)" == value.strip()
      found += 1
    if(key == "test_set_size"):
      assert approx_equal(0.0476, float(value),.0002)
      found += 1
    if(key == "twinned"):
      assert "h,-k,-l" == value.strip()
      found += 1
  #
  assert found == 5

def exercise_01():
  database_dict = load_dict()
  #
  twinned = database_dict["twinned"]
  sel  = twinned != "none"
  sel &= twinned != "false"
  n_twinned = sel.count(True)
  print ("TWINNED: %d (percent: %6.2f), TOTAL: %d" % (n_twinned,
    n_twinned*100./sel.size(), sel.size()))
  r_work_pdb = database_dict["pdb_header_r_work"]
  r_work_cutoff = database_dict["r_work_cutoffs"]
  r_work_re_computed = database_dict["r_work"]
  name = database_dict["pdb_code"]
  #
  sel &= r_work_cutoff != "none"
  sel &= r_work_pdb != "none"
  #
  r_work_pdb         = r_work_pdb.select(sel)
  r_work_cutoff      = r_work_cutoff.select(sel)
  r_work_re_computed = r_work_re_computed.select(sel)
  name               = name.select(sel)
  twinned            = twinned.select(sel)
  #
  def str_to_float(x):
    tmp = flex.double()
    for x_ in x:
      tmp.append(float(x_))
    return tmp
  #
  r_work_cutoff = str_to_float(r_work_cutoff)
  r_work_re_computed = str_to_float(r_work_re_computed)
  r_work_pdb = str_to_float(r_work_pdb)
  #
  delta = (r_work_cutoff - r_work_pdb)*100.
  #
  sp = flex.sort_permutation(delta)
  name               = name         .select(sp)
  delta              = delta        .select(sp)
  r_work_cutoff      = r_work_cutoff.select(sp)
  r_work_pdb         = r_work_pdb   .select(sp)
  r_work_re_computed = r_work_re_computed.select(sp)
  twinned            = twinned.select(sp)
  #
  for n,d,rwc,rwp,rw,t in zip(name,delta,r_work_cutoff,r_work_pdb, r_work_re_computed, twinned):
    print ("%s diff=%6.2f rw_c=%6.4f rw_p=%6.4f rw_ad=%6.4f %s" % (n,d,rwc,rwp,rw, t))

def exercise_02():
  database_dict = load_dict()
  #
  r_work_pdb = database_dict["pdb_header_r_work"]
  r_work_cutoff = database_dict["r_work_cutoffs"]
  r_work_re_computed = database_dict["r_work"]
  name = database_dict["pdb_code"]
  print ("size: ",name.size())
  #
  sel = r_work_pdb != "none"
  print ("NO R in PDB:", sel.count(False))
  #
  r_work_pdb = r_work_pdb.select(sel)
  r_work_cutoff = r_work_cutoff.select(sel)
  r_work_re_computed = r_work_re_computed.select(sel)
  name = name.select(sel)
  #
  def str_to_float(x):
    tmp = flex.double()
    for x_ in x:
      tmp.append(float(x_))
    return tmp
  #
  r_work_cutoff = str_to_float(r_work_cutoff)
  r_work_re_computed = str_to_float(r_work_re_computed)
  r_work_pdb = str_to_float(r_work_pdb)
  #
  delta1 = (r_work_cutoff - r_work_pdb)*100.
  delta2 = (r_work_re_computed - r_work_pdb)*100.
  print ("Number with better R: ", (delta1 <=0.).count(True), (delta2 <=0.).count(True))
  print ("Number with worse  R: ", (delta1 >0.).count(True), (delta2 >0.).count(True))
  #
  #sp = flex.sort_permutation(delta)
  #name          = name         .select(sp)
  #delta         = delta        .select(sp)
  #r_work_cutoff = r_work_cutoff.select(sp)
  #r_work_pdb    = r_work_pdb   .select(sp)
  #r_work_re_computed = r_work_re_computed.select(sp)
  #for n,d,rwc,rwp,rw,t in zip(name,delta,r_work_cutoff,r_work_pdb, r_work_re_computed, twinned):
  #  print "%s %8.4f %8.4f %8.4f %8.4f %s" % (n,d,rwc,rwp,rw, t)
  #print len(twinned)


if (__name__ == "__main__"):
  exercise_00()
  exercise_01()
  exercise_02()
  print("OK")
