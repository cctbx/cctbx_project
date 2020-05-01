from __future__ import absolute_import, division, print_function
import time
from mmtbx import monomer_library
import mmtbx.refinement.real_space.fit_residue
import iotbx.pdb
from mmtbx.rotamer.rotamer_eval import RotamerEval
import mmtbx.utils
import sys
from libtbx import easy_pickle
import libtbx.load_env
from six.moves import range
from libtbx import easy_mp
from libtbx import group_args
import math

mon_lib_srv = monomer_library.server.server()
rotamer_eval = RotamerEval()

pdb_files = [
  "ala.pdb",
  "asn.pdb",
  "asp.pdb",
  "cys.pdb",
  "gln.pdb",
  "glu.pdb",
  "gly.pdb",
  "his.pdb",
  "ile.pdb",
  "leu.pdb",
  "met.pdb",
  "mse.pdb", # is ignored with rotamer named None
  "phe.pdb",
  "pro.pdb", # BAD all-rotamers files
  "ser.pdb",
  "thr.pdb",
  "trp.pdb",
  "tyr.pdb",
  "val.pdb",
  "arg.pdb",
  "lys.pdb"
]

def get_nested_loop(n, fine, start=0, end=360):
  assert n >= 1 and n<=4
  result = []
  #
  if(fine):
    if   (n in [1,2]): step=1
    elif (n == 3):     step=3
  else:
    step = 10
  #
  if(n==1):
    for a1 in range(start,end+step,step):
      result.append([a1])
  elif(n==2):
    for a1 in range(start,end+step,step):
      for a2 in range(start,end+step,step):
        result.append([a1, a2])
  elif(n==3):
    for a1 in range(start,end+step,step):
      for a2 in range(start,end+step,step):
        for a3 in range(start,end+step,step):
          result.append([a1, a2, a3])
  elif(n==4):
    if(fine):
      for a1 in range(start,end+7,7):
        for a2 in range(start,end+8,8):
          for a3 in range(start,end+9,9):
            for a4 in range(start,end+10,10):
              result.append([a1, a2, a3, a4])
    else:
      for a1 in range(start,end+step,step):
        for a2 in range(start,end+step,step):
          for a3 in range(start,end+step,step):
            for a4 in range(start,end+step,step):
              result.append([a1, a2, a3, a4])
  return result

def get_clusters_and_angles(file_name, fine):
  path=libtbx.env.find_in_repositories("mmtbx/idealized_aa_residues/data")
  pdb_inp = iotbx.pdb.input(file_name=path+"/"+file_name)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  xrs = pdb_inp.xray_structure_simple()
  residue = pdb_hierarchy.only_residue()
  clusters = mmtbx.refinement.real_space.aa_residue_axes_and_clusters(
    residue         = residue,
    mon_lib_srv     = mon_lib_srv,
    backbone_sample = False).clusters
  if(len(clusters)==0): return None,None
  nested_loop = get_nested_loop(n=len(clusters), fine=fine)
  return clusters, nested_loop

def chunker(x, dim):
  return (x[i::dim] for i in range(dim))

def run_one(args):
  clusters, chunk, file_name, include, collect_states = args
  #
  path=libtbx.env.find_in_repositories("mmtbx/idealized_aa_residues/data")
  pdb_inp = iotbx.pdb.input(file_name=path+"/"+file_name)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  xrs = pdb_inp.xray_structure_simple()
  residue = pdb_hierarchy.only_residue()
  ri = mmtbx.refinement.real_space.fit_residue.get_rotamer_iterator(
    mon_lib_srv = mon_lib_srv,
    residue     = residue)
  if(len(clusters)==0): return
  for rotamer, rotamer_sites_cart in ri:
    residue.atoms().set_xyz(rotamer_sites_cart)
    xrs = xrs.replace_sites_cart(rotamer_sites_cart)
    if(collect_states):
      states = mmtbx.utils.states(
        xray_structure=xrs, pdb_hierarchy=pdb_hierarchy)
    else:
      states = None # Collecting states with multiprocessing won't work!
    good_angles = mmtbx.refinement.real_space.generate_angles_nested(
      clusters     = clusters,
      residue      = residue,
      rotamer_eval = rotamer_eval,
      nested_loop  = chunk,
      include      = include,
      states       = states)
    break
    #
  good_angles_ = []
  for chi in good_angles:
    chi_ = []
    for ch in chi:
      chi_.append(ch*math.pi/180)
    good_angles_.append(chi_)
  good_angles = good_angles_
  #
  return group_args(good_angles = good_angles, states = states)

def exercise(file_name, include, NPROC=96):
  fine = False
  if(len(include)==2): fine = True
  suffix = "_".join([s.lower() for s in include])
  clusters, nested_loop = get_clusters_and_angles(file_name=file_name, fine=fine)
  if(clusters is None): return
  chunks = list(chunker(nested_loop, NPROC))
  tt = 0
  if(NPROC>1):
    argss = []
    for chunk in chunks:
      argss.append([clusters, chunk, file_name, include, False])
    stdout_and_results = easy_mp.pool_map(
      processes    = NPROC,
      fixed_func   = run_one,
      args         = argss,
      func_wrapper = "buffer_stdout_stderr")
    good_angles = []
    for result in stdout_and_results:
      good_angles.extend(result[1].good_angles)
    states = None # Undefined if multiprocessing is used
  else:
    t0 = time.time()
    args = [clusters, nested_loop, file_name, include, True]
    result = run_one(args)
    good_angles = result.good_angles
    states      = result.states
    tt = time.time()-t0
  if(states is not None):
    states.write(file_name="%s_%s.pdb"%(file_name[:-4],suffix))
  print("file_name, n_clusters, n_good_angles, total:", file_name, \
    len(clusters), len(good_angles), len(nested_loop), tt)
  sys.stdout.flush()
  easy_pickle.dump(
    file_name="%s_%s.pkl"%(file_name[:-4], suffix),
    obj=good_angles)

if(__name__ == "__main__"):
  for it in [["FAVORED","ALLOWED"], ["FAVORED"]]:
    print (it, "-"*20)
    if(len(sys.argv[1:])==1):
      exercise(file_name=sys.argv[1:][0], include=it)
    else:
      for fn in pdb_files:
        exercise(file_name=fn, include=it)
