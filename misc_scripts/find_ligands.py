"""
Utility script to process bulk (model,data) pairs
"""

from __future__ import absolute_import, division, print_function
import os, sys
from scitbx.array_family import flex
import mmtbx.model
from libtbx import easy_mp
from collections import OrderedDict
import iotbx.cif
import traceback
import os
from iotbx import reflection_file_reader
import mmtbx.model
import mmtbx.monomer_library.server
mon_lib_srv = mmtbx.monomer_library.server.server()
import cctbx.miller
from cctbx import maptbx, crystal
from cctbx.crystal import super_cell

cif_files = os.getenv('PDB_MIRROR_MMCIF')
pdbmtz = os.getenv('PDBMTZ')
SUPPORTED = ["H","C","N","O","S","SE","B","SI","P","AS","F","CL","BR","I"]


def get_model_files_dict(path):
  index_file = "/".join([path,"INDEX"])
  ifn = open(index_file,"r")
  result = {}
  for l in ifn.readlines():
    l = l.strip()
    file_name = "/".join([path,l])
    #assert os.path.isfile(file_name) # Terribly runtime expensive
    code = l[-11:-7]
    result[code] = file_name
  return result

def get_mtz_dict(path):
  result     = OrderedDict()
  codes      = flex.std_string()
  file_names = flex.std_string()
  sizes      = flex.double()
  cntr = 0
  for l in os.listdir(path):
    l = l.strip()
    if not l.endswith(".mtz"): continue
    cntr += 1
    code = l[:-4]
    #
    # FOR DEBUGGING
    #
    #if(cntr==1000): break
    #
    #if code != "6cg3": continue
    #
    file_name = "/".join([path,l])
    #assert os.path.isfile(file_name) # Terribly runtime expensive
    codes.append(code)
    file_names.append(file_name)
    sizes.append(os.path.getsize(file_name))
  sel = flex.sort_permutation(sizes)
  codes      = codes     .select(sel)
  file_names = file_names.select(sel)
  sizes      = sizes     .select(sel)
  for c,f,s in zip(codes, file_names, sizes):
    result[c] = [f,s]
  return result

def has_ligands(h):
  get_class = iotbx.pdb.common_residue_names_get_class
  skip_classes = {
      "common_water",
      "common_amino_acid",
      "modified_amino_acid",
      "d_amino_acid",
      "common_rna_dna"
  }
  ligand_selections = []
  if len(list(h.models())) > 1: return []
  for m in h.models():
    for c in m.chains():
      for rg in c.residue_groups():
        # Alt-Conf Check: If there is more than 1 atom_group, it has altlocs
        if len(rg.atom_groups()) > 1: continue
        ag = rg.only_atom_group()
        # Failsafe: Skip if it's a pure alt-conf with no main conformer
        if ag.altloc.strip() != '': continue
        if get_class(name=ag.resname) in skip_classes: continue
        atoms = ag.atoms()
        nat = atoms.size()
        if nat < 10: continue
        elements = [a.element.strip().upper() for a in atoms]
        if any(e not in SUPPORTED for e in elements): continue
        occ = atoms.extract_occ()
        if not occ.all_eq(1.0): ontinue
        r = rg.conformers()[0].only_residue()
        ma = r.missing_atoms(mon_lib_srv=mon_lib_srv)
        if ma is None: continue
        ma = [x for x in ma if x != "OXT"]
        if len(ma) > 0: continue
        print(ag.resname, nat, list(occ), occ.all_eq(1.0))
        ligand_selections.append(atoms.extract_i_seq())
  return ligand_selections

def check_maps(miller_arrays, ligand_selections, code, sites_frac, hd_selection):
  def _get_real_map(mc):
    fft_map = cctbx.miller.fft_map(
      crystal_gridding     = crystal_gridding,
      fourier_coefficients = mc)
    fft_map.apply_sigma_scaling()
    return fft_map.real_map_unpadded()
  mc_tfofc, mc_fofc = None, None
  for ma in miller_arrays:
    lbls = ma.info().label_string()
    if "2FOFCWT,PHI2FOFCWT" in lbls: mc_tfofc = ma.deep_copy()
    if "FOFCWT,PHIFOFCWT" in lbls:   mc_fofc  = ma.deep_copy()
  assert mc_tfofc is not None and mc_fofc is not None
  crystal_gridding = mc_tfofc.crystal_gridding(
    d_min             = mc_tfofc.d_min(),
    symmetry_flags    = maptbx.use_space_group_symmetry,
    resolution_factor = 1./3)
  tfofc = _get_real_map(mc = mc_tfofc)
  fofc  = _get_real_map(mc = mc_fofc)
  result = []
  for ls in ligand_selections:
    good = True
    for i in ls:
      if hd_selection[i]: continue
      site_frac = sites_frac[i]
      mv_tfofc = tfofc.tricubic_interpolation(site_frac)
      mv_fofc  = fofc .tricubic_interpolation(site_frac)
      if mv_tfofc < 0.7 or abs(mv_fofc) > 2.5:
        good = False
        break
    if good: result.append(ls)
  return result

def get_complete_hierarchy(h, sel_within):
  keep = flex.size_t()
  for m in h.models():
    for c in m.chains():
      for rg in c.residue_groups():
        if len(rg.conformers())>1: continue
        iseqs = rg.atoms().extract_i_seq()
        for i in iseqs:
          if i in sel_within:
            keep.extend(iseqs)
            break
  return keep

def run_one(args):
  cif, hkl, code = args
  #
  try:
    #
    # DATA
    #
    miller_arrays = reflection_file_reader.any_reflection_file(file_name =
      hkl).as_miller_arrays()
    d_min = miller_arrays[0].d_min()
    if d_min > 2.5: return None
    #
    # MODEL
    #
    model = mmtbx.model.manager(
      model_input         = iotbx.pdb.input(file_name=cif),
      stop_for_unknowns   = True,
      skip_ss_annotations = True,
      process_biomt       = False)
    model = model.remove_hydrogens()
    h = model.get_hierarchy()
    sp = model.selection(string="protein")
    if sp.count(True)*100./sp.size() < 80.: return None # not protein

    if len(h.models()) != 1: return None
    ligand_selections = has_ligands(h)
    if ligand_selections is None or len(ligand_selections)==0: return None
    sites_frac = model.get_sites_frac()
    ligand_selections = check_maps(
      miller_arrays     = miller_arrays,
      ligand_selections = ligand_selections,
      code              = code,
      sites_frac        = sites_frac,
      hd_selection      = model.get_xray_structure().hd_selection())
    if len(ligand_selections)==0: return None
    #
    suse = super_cell.manager(
      pdb_hierarchy        = model.get_hierarchy(),
      crystal_symmetry     = model.crystal_symmetry(),
      select_within_radius = 4.0,
      #box                  = True
      box_super_sphere=True
      )
    suseh = suse.super_sphere_hierarchy
    size = suseh.atoms().size()
    asu_mappings=suse.cs_super_sphere.special_position_settings().asu_mappings(
        buffer_thickness=4.0,
        sites_cart=suseh.atoms().extract_xyz())
    #suseh.write_pdb_file("suseh.pdb")
    #
    atoms = suseh.atoms()
    for i, ls in enumerate(ligand_selections):
      sel_within = crystal.neighbors_fast_pair_generator(
        asu_mappings=asu_mappings,
        distance_cutoff=4.0).neighbors_of(
          primary_selection=flex.bool(size, ls)).iselection()

      #diff = flex.size_t([x for x in sel_within if x not in ls])
      #suseh.select(sel_within).write_pdb_file("zz.pdb")
      #suseh.select(ls).write_pdb_file("zz_ls.pdb")
      #suseh.select(diff).write_pdb_file("zz_diff.pdb")

      keep = get_complete_hierarchy(h=suseh, sel_within=sel_within)
      resname = atoms[ls[0]].parent().resname.strip()
      resseq  = atoms[ls[0]].parent().parent().resseq.strip()
      cid     = atoms[ls[0]].parent().parent().parent().id.strip()
      rid = "%s_%s_%s"%(resname, cid, resseq)

      suseh_keep = suseh.select(keep)

      elements  = [a.element.strip().upper() for a in suseh_keep.atoms()]
      found_unsupported = False
      for e in elements:
        if not e in SUPPORTED: found_unsupported = True
      if found_unsupported: continue

      suseh_keep.write_pdb_or_mmcif_file(
        target_filename  = "%s_%s.cif" % (rid, code),
        target_format    = "mmcif",
        crystal_symmetry = suse.cs_super_sphere)

  #
  except Exception as e:
    of = open("%s.log"%code, "w")
    traceback.print_exc(file=of)
    of.close()

def run(cmdargs, NPROC=128):
  #
  processed = []
  for f in os.listdir("."):
    if(f.endswith(".cif")):
      code = f.rsplit("_", 1)[-1].split(".")[0]
      processed.append(code)
  #
  cifs = get_model_files_dict(path=cif_files)
  hkls = get_mtz_dict(path=pdbmtz)
  print("Number of cifs:", len(cifs.keys()))
  print("Number of hkls:", len(hkls.keys()))
  argss = []
  for code in hkls.keys():
    if(code in processed): continue
    hkl,size = hkls[code]
    try:
      cif = cifs[code]
      argss.append([cif, hkl, code])
    except: pass # intentional
  print("NEW FILES TO PROCESS:", len(argss))
  if(NPROC>1):
    stdout_and_results = easy_mp.pool_map(
      processes    = NPROC,
      fixed_func   = run_one,
      args         = argss,
      func_wrapper = "buffer_stdout_stderr")
  else:
    for args in argss:
      run_one(args)

if (__name__ == "__main__"):
  run(sys.argv[1:])
