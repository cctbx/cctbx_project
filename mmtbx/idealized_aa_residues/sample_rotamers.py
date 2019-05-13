from __future__ import division
from __future__ import print_function
import time
from mmtbx import monomer_library
import mmtbx.refinement.real_space.fit_residue
import iotbx.pdb
from mmtbx.rotamer.rotamer_eval import RotamerEval
import mmtbx.utils
from scitbx.matrix import rotate_point_around_axis
import sys
from libtbx import easy_pickle
import libtbx.load_env

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

def torsion_search_nested(residue, clusters, rotamer_eval, states):
  n_angles = len(clusters)
  rid = rotamer_eval.evaluate_residue(residue = residue)
  xyz_moved_dc = residue.atoms().extract_xyz().deep_copy()
  nested_loop = []
  if(n_angles==1):
    for a1 in range(0,363,3):
      nested_loop.append([a1])
  elif(n_angles==2):
    for a1 in range(0,370,10):#range(0,363,3):
      for a2 in range(0,370,10):#range(0,363,3):
        nested_loop.append([a1, a2])
  elif(n_angles==3):
    for a1 in range(0,370,10):#range(0,365,5):
      for a2 in range(0,370,10):#range(0,365,5):
        for a3 in range(0,370,10):#range(0,365,5):
          nested_loop.append([a1, a2, a3])
  #lif(n_angles==4):
  # for a1 in range(0,365,5):
  #   for a2 in range(0,365,5):
  #     for a3 in range(0,365,5):
  #       for a4 in range(0,365,5):
  #         nested_loop.append([a1, a2, a3, a4])
  elif(n_angles==4):
    for a1 in range(0,370,10):
      for a2 in range(0,370,10):
        for a3 in range(0,370,10):
          for a4 in range(0,370,10):
            nested_loop.append([a1, a2, a3, a4])
  #elif(n_angles==4):
  #  for a1 in range(0,365,5):
  #    for a2 in range(0,367,7):
  #      for a3 in range(0,369,9):
  #        for a4 in range(0,371,11):
  #          nested_loop.append([a1, a2, a3, a4])
  else: assert 0
  good = []
  good_sites = []
  for ia, angles in enumerate(nested_loop):
    xyz_moved = xyz_moved_dc.deep_copy()
    for i, angle in enumerate(angles):
      cl = clusters[i]
      for atom in cl.atoms_to_rotate:
        new_xyz = rotate_point_around_axis(
          axis_point_1 = xyz_moved[cl.axis[0]],
          axis_point_2 = xyz_moved[cl.axis[1]],
          point        = xyz_moved[atom],
          angle        = angle, deg=True)
        xyz_moved[atom] = new_xyz
    residue.atoms().set_xyz(xyz_moved)
    fl = rotamer_eval.evaluate_residue_2(residue = residue)
    #if(fl != "OUTLIER" and str(fl).upper() != "NONE"):
    if(fl == "Favored"):
      states.add(sites_cart=xyz_moved)
      good.append(angles)
  return states, good, nested_loop

def exercise(file_name):
  path=libtbx.env.find_in_repositories("mmtbx/idealized_aa_residues/data")
  pdb_inp = iotbx.pdb.input(file_name=path+"/"+file_name)
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  xrs = pdb_inp.xray_structure_simple()
  residue = pdb_hierarchy.only_residue()
  clusters = mmtbx.refinement.real_space.aa_residue_axes_and_clusters(
    residue         = residue,
    mon_lib_srv     = mon_lib_srv,
    backbone_sample = False).clusters
  ri = mmtbx.refinement.real_space.fit_residue.get_rotamer_iterator(
    mon_lib_srv = mon_lib_srv,
    residue     = residue)
  if(len(clusters)==0): return
  for rotamer, rotamer_sites_cart in ri:
    residue.atoms().set_xyz(rotamer_sites_cart)
    xrs= xrs.replace_sites_cart(rotamer_sites_cart)
    states = mmtbx.utils.states(xray_structure=xrs, pdb_hierarchy=pdb_hierarchy)
    t0 = time.time()
    states, good_angles, nested_loop = torsion_search_nested(
      residue      = residue,
      clusters     = clusters,
      rotamer_eval = rotamer_eval,
      states       = states)
    tt = time.time()-t0
    states.write(file_name="%s_all-coarse_step10.pdb"%file_name[:-4])
    break
  print("file_name, n_clusters, n_good_angles, total:", file_name, \
    len(clusters), len(good_angles), len(nested_loop), tt)
  easy_pickle.dump(
    file_name="%s-coarse_step10_favored.pickle"%file_name[:-4],
    obj=good_angles)

if(__name__ == "__main__"):
  if(len(sys.argv[1:])==1):
    exercise(file_name=sys.argv[1:][0])
  else:
    for fn in pdb_files:
      exercise(file_name=fn)
