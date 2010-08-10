# LIBTBX_SET_DISPATCHER_NAME phenix.r_factor_statistics

import os,sys
from cctbx.array_family import flex
from libtbx import easy_pickle
import libtbx.load_env
from libtbx.str_utils import format_value

def show_histogram(data, n_slots):
  hm = flex.histogram(data = data, n_slots = n_slots)
  lc_1 = hm.data_min()
  s_1 = enumerate(hm.slots())
  for (i_1,n_1) in s_1:
    hc_1 = hm.data_min() + hm.slot_width() * (i_1+1)
    print "%10.3f - %-10.3f : %d" % (lc_1, hc_1, n_1)
    lc_1 = hc_1

help_message = """\
phenix.r_factor_statistics: a tool to show Rwork, Rfree and Rfree-Rwork
                            statistics for PDB structures.

Usage:
  - running phenix.r_factor_statistics without arguments will consider all PDB
    structures;
  - runing with single numeric argument
      phenix.r_factor_statistics 2.5
    will consider structures in resolution range [2.5-0.1, 2.5+0.1]
  - one can specify left and right offsets as well as the number of bins
      phenix.r_factor_statistics 2.5 left_offset=0.1 right_offset=0.5 n_bins=10
"""

def run(args, left_offset=0.1, right_offset=0.1, n_bins=10):
  need_help_msg = False
  if(len(args)==0): need_help_msg = True
  else:
    for arg in args:
      if(str(arg).lower() in ["help","-h","--help","--h","h"]):
        need_help_msg = True
        break
  if(need_help_msg):
    print help_message
  #
  def get_value(x):
    result = None
    try:
      result = float(x[x.index("=")+1:])
    except IndexError: pass
    except ValueError: pass
    return result
  d_min = None
  for arg in args:
    arg = str(arg).lower()
    try: d_min = float(arg)
    except ValueError: pass
    if(arg.count("left_offset")):
      x = get_value(arg)
      if(x is not None): left_offset=x
    if(arg.count("right_offset")):
      x = get_value(arg)
      if(x is not None): right_offset=x
    if(arg.count("n_bins")):
      x = int(get_value(arg))
      if(x is not None): n_bins=x
  dl,dr = [0,1.e+6]
  if(d_min is not None):
    cmd = "phenix.r_factor_statistics %s left_offset=%s right_offset=%s n_bins=%s"%(
      format_value("%6.3f",d_min).strip(),
      format_value("%6.3f",left_offset).strip(),
      format_value("%6.3f",right_offset).strip(),
      format_value("%d",n_bins).strip())
    print "Command used:\n\n%s\n\n"%cmd
    dl = d_min-left_offset
    dr = d_min+right_offset
  file = libtbx.env.find_in_repositories(relative_path=
    "chem_data/polygon_data/all_mvd.pickle",
    test=os.path.isfile)
  database_dict = easy_pickle.load(file)
  r_work_pdb = database_dict["pdb_header_r_work"]
  r_free_pdb = database_dict["pdb_header_r_free"]
  d_min = database_dict["high_resolution"]
  sel = r_work_pdb != "none"
  sel &= r_free_pdb != "none"
  sel &= d_min != "none"
  #
  r_work_pdb = r_work_pdb.select(sel)
  r_free_pdb = r_free_pdb.select(sel)
  d_min      = d_min.select(sel)
  #
  def str_to_float(x):
    tmp = flex.double()
    for x_ in x:
      tmp.append(float(x_))
    return tmp
  #
  d_min = str_to_float(d_min)
  r_work_pdb = str_to_float(r_work_pdb)
  r_free_pdb = str_to_float(r_free_pdb)
  diff = r_free_pdb - r_work_pdb
  #
  sel  = d_min > dl
  sel &= d_min < dr
  sel &= diff > 0.
  sel &= diff < 0.1
  #
  r_work_pdb = r_work_pdb.select(sel)
  r_free_pdb = r_free_pdb.select(sel)
  d_min      = d_min.select(sel)
  diff = diff.select(sel)
  #
  print "Histogram of Rwork for models in PDB at resolution %4.2f-%4.2f A:"%(dl,dr)
  show_histogram(data = r_work_pdb, n_slots=n_bins)
  print "Histogram of Rfree for models in PDB at resolution %4.2f-%4.2f A:"%(dl,dr)
  show_histogram(data = r_free_pdb, n_slots=n_bins)
  print "Histogram of Rfree-Rwork for all model in PDB at resolution %4.2f-%4.2f A:"%(dl,dr)
  show_histogram(data = diff, n_slots=n_bins)
  print "Number of structures considered:", diff.size()

if (__name__ == "__main__"):
  run(args=sys.argv[1:])
