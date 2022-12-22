from __future__ import absolute_import, division, print_function
import time
from mmtbx import monomer_library
import mmtbx.refinement.real_space.fit_residue
import iotbx.pdb
from mmtbx.rotamer.rotamer_eval import RotamerEval
import mmtbx.utils
from scitbx.matrix import rotate_point_around_axis
import libtbx.load_env

mon_lib_srv = monomer_library.server.server()
rotamer_eval = RotamerEval()
rotamer_manager = mmtbx.idealized_aa_residues.rotamer_manager.load(
  rotamers="favored")

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
  "ser.pdb",
  "thr.pdb",
  "trp.pdb",
  "tyr.pdb",
  "val.pdb",
  "arg.pdb",
  "lys.pdb"
]

def sample(residue, clusters, states):
  xyz_moved_dc = residue.atoms().extract_xyz().deep_copy()
  chi_angles = rotamer_manager.get_chi_angles(resname=residue.resname)
  print(residue.resname, len(chi_angles))
  if(residue.resname in ["ARG","LYS"]):
    return None # skip: too long to run
  if(residue.resname in ["ALA","GLY"]):
    assert len(chi_angles) == 0
  else:
    len(chi_angles) > 0
  for ia, angles in enumerate(chi_angles):
    xyz_moved = xyz_moved_dc.deep_copy()
    for i, angle in enumerate(angles):
      cl = clusters[i]
      for atom in cl.atoms_to_rotate:
        new_xyz = rotate_point_around_axis(
          axis_point_1 = xyz_moved[cl.axis[0]],
          axis_point_2 = xyz_moved[cl.axis[1]],
          point        = xyz_moved[atom],
          angle        = angle, deg=False)
        xyz_moved[atom] = new_xyz
    residue.atoms().set_xyz(xyz_moved)
    fl = rotamer_eval.evaluate_residue(residue = residue)
    assert fl != "OUTLIER"
    if(states is not None): states.add(sites_cart=xyz_moved)
  return states

def exercise(file_name, write_pdb_file=False):
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
    states = None
    if(write_pdb_file):
      states = mmtbx.utils.states(xray_structure=xrs,
        pdb_hierarchy=pdb_hierarchy)
    states_result = sample(
      residue      = residue,
      clusters     = clusters,
      states       = states)
    if(write_pdb_file):
      states.write(file_name="%s_all-coarse_step10.pdb"%file_name[:-4])
    break # By convention, use first rotamer!

if(__name__ == "__main__"):
  t0 = time.time()
  for fn in pdb_files:
    exercise(file_name=fn)
  print("Time: %8.4f"%(time.time()-t0))
  print("OK")
