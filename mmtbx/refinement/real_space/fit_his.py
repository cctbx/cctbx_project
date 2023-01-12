from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from cctbx import maptbx
import mmtbx.idealized_aa_residues.rotamer_manager
from libtbx import group_args

from mmtbx.refinement.real_space import fit_residue

import boost_adaptbx.boost.python as bp
ext = bp.import_ext("mmtbx_rotamer_fit_ext")

from mmtbx.idealized_aa_residues import all_his
from scitbx.math import superpose

import mmtbx.monomer_library.server
mon_lib_srv = mmtbx.monomer_library.server.server()

non_h_ring_atom_names = ["CG", "ND1", "CD2", "CE1", "NE2"]

def manageable_HIS(rg):
  result = group_args(
    is_true_altloc = False,
    has_h          = False,
    has_d          = False)
  resnames=[]
  elements = []
  atom_groups = rg.atom_groups()
  coords = []
  if(len(atom_groups)>3): return None # True altlocs
  for i, ag in enumerate(rg.atom_groups()):
    if(ag.resname.strip().upper()!="HIS"): return None # Not HIS
    elements.extend( [e.strip().upper() for e in ag.atoms().extract_element()] )
    for a in ag.atoms():
      if(a.name.strip().upper()=="CE1"):
        coords.append(a.xyz)
  if(len(coords)!=1): return None
  elements = list(set(elements))
  print(elements)
  result.has_h = "H" in elements
  result.has_d = "D" in elements
  return result

def match_residue_groups(rg_working, rg_to_insert):
  # Rigid-body superpose rings
  atoms_fixed  = []
  atoms_moving = []
  for n in non_h_ring_atom_names:
    for atom_to_insert in rg_to_insert.atoms():
      if(atom_to_insert.name.strip().upper() == n):
        atoms_moving.append(atom_to_insert)
    for atom_working in rg_working.atoms():
      if(atom_working.name.strip().upper() == n):
        atoms_fixed.append(atom_working)
  sites_fixed  = flex.vec3_double([a.xyz for a in atoms_fixed])
  sites_moving = flex.vec3_double([a.xyz for a in atoms_moving])
  lsq_fit_obj = superpose.least_squares_fit(
    reference_sites = sites_fixed,
    other_sites     = sites_moving)
  sites_moved = lsq_fit_obj.r.elems*rg_to_insert.atoms().extract_xyz()+lsq_fit_obj.t.elems
  # Overlay ring
  for i, atom_moving in enumerate(rg_to_insert.atoms()):
    atom_moving.set_xyz(sites_moved[i])
  # Then match the rest
  for atom_to_insert in rg_to_insert.atoms():
    for atom_working in rg_working.atoms():
      if(atom_to_insert.name.strip().upper() !=
           atom_working.name.strip().upper()): continue
      atom_to_insert.set_xyz(atom_working.xyz)
      atom_to_insert.set_b(atom_working.b)
      atom_to_insert.set_charge(atom_working.charge)
      atom_to_insert.set_fdp(atom_working.fdp)
      atom_to_insert.set_fp(atom_working.fp)
      atom_to_insert.set_occ(atom_working.occ)
      atom_to_insert.set_segid(atom_working.segid)
      atom_to_insert.set_serial(atom_working.serial)
      atom_to_insert.set_uij(atom_working.uij)
      if(not atom_to_insert.name.strip().upper() in non_h_ring_atom_names):
        atom_to_insert.set_xyz(atom_working.xyz)

def run(hierarchy, map_data, mon_lib_srv, rotamer_manager, unit_cell):
  allhis = all_his.load()

  def compute_target(rg):
     sites_cart  = rg.atoms().extract_xyz()
     t = maptbx.real_space_target_simple(
        unit_cell   = unit_cell,
        density_map = map_data,
        sites_cart  = sites_cart)
     return t/sites_cart.size()

  def per_atom(rg):
    #hs = ["CG", "ND1", "CE1", "NE2", "CD2", "HD1", "HE1", "HE2", "HD2"]
    hs = ["HD1", "HE1", "HE2", "HD2",   "DD1", "DE1", "DE2", "DD2"]
    vals = []
    for a in rg.atoms():
      if a.name.strip().upper() in hs:
        t = maptbx.real_space_target_simple(
          unit_cell   = unit_cell,
          density_map = map_data,
          sites_cart  = flex.vec3_double([a.xyz,]))
        vals.append("%s %5.2f"%(a.name.strip().upper(),t))
    return " ".join(vals)


  for c in hierarchy.chains():
    for i, rg in enumerate(c.residue_groups()):
      o = manageable_HIS(rg=rg)
      if(o is None): continue

      print("chain %s resid %s"%(c.id, rg.resid()),"-"*20)

      rg_dc = rg.detached_copy()
      rg_w  = rg
      print("  t orig:", compute_target(rg), per_atom(rg))
      for i, rg_to_insert in enumerate(allhis.main_three()):
        print("  protonation", i)

        rg_to_insert = allhis.r["d1"]

        rg_to_insert.resseq = rg_w.resseq
        match_residue_groups(rg_working=rg_w, rg_to_insert=rg_to_insert)
        print("    t start:", compute_target(rg_to_insert), per_atom(rg_to_insert))

        o = fit_residue.tune_up(
          target_map           = map_data,
          residue              = rg_to_insert,
          mon_lib_srv          = mon_lib_srv,
          rotamer_manager      = None,
          unit_cell            = unit_cell,
          exclude_hd           = False,
          torsion_search_start = -5,
          torsion_search_stop  = 5,
          torsion_search_step  = 0.1)
        print("    t final:", compute_target(o.residue), per_atom(o.residue))

        rg_to_insert = o.residue.detached_copy()
        c.remove_residue_group(rg_w)
        rg_to_insert.link_to_previous=True
        c.insert_residue_group(i=i, residue_group=rg_to_insert)
        rg_w  = rg_to_insert
        break

        if 0:
          c.remove_residue_group(rg_w)
          rg_to_insert.link_to_previous=True
          c.insert_residue_group(i=i, residue_group=rg_to_insert)
          rg_w  = rg_to_insert

  #
  hierarchy.write_pdb_file(file_name = "inserted.pdb")
