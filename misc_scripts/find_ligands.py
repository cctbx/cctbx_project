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
from scitbx.array_family import flex
from cctbx import crystal, sgtbx
import iotbx.pdb
from cctbx import uctbx
from libtbx.utils import Sorry

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
    #if(cntr==10000): break
    #
    #if code != "6p1y": continue
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

def has_ligands(h, code):
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
        if not occ.all_eq(1.0): continue
        r = rg.conformers()[0].only_residue()
        ma = r.missing_atoms(mon_lib_srv=mon_lib_srv)
        if ma is None: continue
        ma = [x for x in ma if x != "OXT"]
        if len(ma) > 0: continue
        print(code, ag.resname, nat, occ.all_eq(1.0))
        ligand_selections.append(atoms.extract_i_seq())
  return ligand_selections

def check_maps(tfofc, fofc, ligand_selections, sites_frac):
  result = []
  for ls in ligand_selections:
    good = True
    for i in ls:
      site_frac = sites_frac[i]
      mv_tfofc = tfofc.tricubic_interpolation(site_frac)
      mv_fofc  = fofc .tricubic_interpolation(site_frac)
      if mv_tfofc < 0.7 or abs(mv_fofc) > 2.5:
        good = False
        break
    if good: result.append(ls)
  return result

def check_maps2(tfofc, fofc, sites_cart, unit_cell):
  m1 = flex.double()
  m2 = flex.double()
  for site_cart in sites_cart:
    site_frac = unit_cell.fractionalize(site_cart)
    mv_tfofc = tfofc.tricubic_interpolation(site_frac)
    mv_fofc  = fofc .tricubic_interpolation(site_frac)
    m1.append(mv_tfofc)
    m2.append(mv_fofc)
  print("  min tfofc, max |fofc|:", flex.min(m1), flex.max(flex.abs(m2)))
  if flex.min(m1) < 0.5 or flex.max(flex.abs(m2)) > 3.0:
    return False
  return True



from scitbx.array_family import flex
from cctbx import crystal, sgtbx
import iotbx.pdb
import iotbx.pdb.utils

def extract_ligand_environment(hierarchy, ligand_iselection, crystal_symmetry, distmax):
    """
    Extracts the ligand and all complete surrounding residues within distmax.
    Explicitly transforms symmetry-related residues to map them next to the ligand,
    assigning them unique chain IDs to prevent duplicate atom label errors.

    Returns None if any extracted residue has an alternative conformation.
    """
    sites_cart = hierarchy.atoms().extract_xyz()
    n_atoms = sites_cart.size()

    # 1. Create a boolean selection array for the ligand
    ligand_sel = flex.bool(n_atoms, False)
    ligand_sel.set_selected(ligand_iselection, True)

    unit_cell = crystal_symmetry.unit_cell()
    sps = crystal_symmetry.special_position_settings()

    # 2. Use special_position_settings to account for crystal symmetry
    asu_mappings = sps.asu_mappings(buffer_thickness=distmax)
    asu_mappings.process_sites_cart(
        original_sites=sites_cart,
        site_symmetry_table=sps.site_symmetry_table(sites_cart)
    )

    pair_generator = crystal.neighbors_fast_pair_generator(
        asu_mappings=asu_mappings,
        distance_cutoff=distmax
    )

    # 3. Map each i_seq to its parent residue_group and chain ID for easy retrieval
    atom_to_rg = {}
    rg_to_chain_id = {}
    for model in hierarchy.models():
        for chain in model.chains():
            for rg in chain.residue_groups():
                for atom in rg.atoms():
                    atom_to_rg[atom.i_seq] = rg
                rg_to_chain_id[rg] = chain.id

    # 4. Collect required residue groups and their specific symmetry operations
    env_rgs = {}

    identity_op = sgtbx.rt_mx("x,y,z")
    for i_seq in ligand_iselection:
        rg = atom_to_rg[i_seq]
        env_rgs[(rg, str(identity_op))] = (rg, identity_op, rg_to_chain_id[rg])

    for pair in pair_generator:
        i_in_lig = ligand_sel[pair.i_seq]
        j_in_lig = ligand_sel[pair.j_seq]

        if not i_in_lig and not j_in_lig:
            continue

        if i_in_lig:
            rt_mx_i = asu_mappings.get_rt_mx_i(pair)
            rt_mx_j = asu_mappings.get_rt_mx_j(pair)
            rt_mx_ji = rt_mx_i.inverse().multiply(rt_mx_j)

            rg_j = atom_to_rg[pair.j_seq]
            env_rgs[(rg_j, str(rt_mx_ji))] = (rg_j, rt_mx_ji, rg_to_chain_id[rg_j])

        if j_in_lig:
            rt_mx_i = asu_mappings.get_rt_mx_i(pair)
            rt_mx_j = asu_mappings.get_rt_mx_j(pair)
            rt_mx_ij = rt_mx_j.inverse().multiply(rt_mx_i)

            rg_i = atom_to_rg[pair.i_seq]
            env_rgs[(rg_i, str(rt_mx_ij))] = (rg_i, rt_mx_ij, rg_to_chain_id[rg_i])

    # 5. Build the new hierarchy
    new_hierarchy = iotbx.pdb.hierarchy.root()
    new_model = iotbx.pdb.hierarchy.model()
    new_hierarchy.append_model(new_model)

    chains_dict = {}

    # Pre-generate all valid PDB chain IDs ("A", "B", ... "AA", "AB"...)
    available_chain_ids = iotbx.pdb.utils.all_chain_ids()

    # Remove chain IDs that are already used in the original hierarchy to be safe
    for cid in hierarchy.chain_ids():
        if cid in available_chain_ids:
            available_chain_ids.remove(cid)

    for (rg, rt_mx, chain_id) in env_rgs.values():

        # Reject if the residue has an alternative conformation
        if len(rg.atom_groups()) > 1 or rg.atom_groups()[0].altloc.strip() != "":
            return None

        rt_str = str(rt_mx)

        # Determine the chain ID to use for this specific (chain, symmetry_operation) pair
        if (chain_id, rt_str) not in chains_dict:

            # If it is the original ASU position, keep the original chain ID
            if rt_mx.is_unit_mx():
                new_chain_id = chain_id
            else:
                # If it is a symmetry mate, pop a completely new chain ID
                new_chain_id = available_chain_ids.pop(0)

            new_chain = iotbx.pdb.hierarchy.chain(id=new_chain_id)
            chains_dict[(chain_id, rt_str)] = new_chain
            new_model.append_chain(new_chain)

        rg_copy = rg.detached_copy()

        # Apply the exact symmetry operation so the residue maps right next to the ligand
        if not rt_mx.is_unit_mx():
            for atom in rg_copy.atoms():
                site_frac = unit_cell.fractionalize(atom.xyz)
                new_site_frac = rt_mx * site_frac
                atom.xyz = unit_cell.orthogonalize(new_site_frac)

        chains_dict[(chain_id, rt_str)].append_residue_group(rg_copy)

    new_hierarchy.atoms().reset_i_seq()

    return new_hierarchy

def get_maps(miller_arrays):
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
  return tfofc, fofc

def center_hierarchy_in_box(hierarchy, default_buffer=5.0):
  """
  Creates a new P1 crystal symmetry box around the hierarchy,
  centers the molecule exactly in the middle of the box,
  and updates the coordinates in-place.

  Returns the modified hierarchy and the new P1 crystal_symmetry object.
  """
  # Extract the current Cartesian coordinates
  sites_cart = hierarchy.atoms().extract_xyz()
  # Calculate the new box and the required centering shift
  box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
      sites_cart=sites_cart,
      buffer_layer=default_buffer
  )
  # Apply the shifted, centered coordinates back to the hierarchy in-place
  hierarchy.atoms().set_xyz(new_xyz=box.sites_cart)
  # Extract the newly generated P1 crystal symmetry
  new_crystal_symmetry = box.crystal_symmetry()
  return hierarchy, new_crystal_symmetry

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
    cs = model.crystal_symmetry()
    model = model.remove_hydrogens()
    h = model.get_hierarchy()
    atoms = h.atoms()
    atoms.reset_i_seq()
    sp = model.selection(string="protein")
    if sp.count(True)*100./sp.size() < 80.: return None # not protein

    if len(h.models()) != 1: return None
    ligand_selections = has_ligands(h, code)
    if ligand_selections is None or len(ligand_selections)==0: return None
    #
    # Compute maps
    #
    tfofc, fofc = get_maps(miller_arrays = miller_arrays)
    #
    # Check ligands against maps
    #
    sites_frac = model.get_sites_frac()
    ligand_selections = check_maps(
      tfofc             = tfofc,
      fofc              = fofc,
      ligand_selections = ligand_selections,
      sites_frac        = sites_frac)
    if len(ligand_selections)==0:
      print("poor ligand map")
      return None
    for i, ls in enumerate(ligand_selections):
      around_ligand = extract_ligand_environment(
        hierarchy         = h,
        ligand_iselection = ls,
        crystal_symmetry  = cs,
        distmax           = 4.0)
      if around_ligand is None:
        print(i, "poor env (altlocs)")
        continue
      # Check for unsupported
      elements  = [a.element.strip().upper() for a in around_ligand.atoms()]
      found_unsupported = False
      for e in elements:
        if not e in SUPPORTED: found_unsupported = True
      if found_unsupported:
        print(i, "poor env (unsupp)")
        continue
      # Check whole thing against maps
      mcheck = check_maps2(
        tfofc      = tfofc,
        fofc       = fofc,
        sites_cart = around_ligand.atoms().extract_xyz(),
        unit_cell  = cs.unit_cell())
      if mcheck is False:
        print(i, "poor env (map)")
        continue
      # Make box around new hierarchy
      box_h, box_cs = center_hierarchy_in_box(hierarchy = around_ligand)
      # Create model object
      try:
        model_box = mmtbx.model.manager(
          pdb_hierarchy    = box_h,
          crystal_symmetry = box_cs)
        model_box.process(make_restraints = True)
      except Sorry as e:
        if "Fatal problems interpreting model file" in str(e): continue
      except KeyboardInterrupt: raise
      # Check geometry
      #g = model_box.geometry_statistics().result(slim = True)
      #
      #print("     bond ", g.bond.mean,  g.bond.max)
      #print("     angle", g.angle.mean, g.angle.max)
      #print("     nb", g.nonbonded.min)

      # Format output file name prefix
      resname = atoms[ls[0]].parent().resname.strip()
      resseq  = atoms[ls[0]].parent().parent().resseq.strip()
      cid     = atoms[ls[0]].parent().parent().parent().id.strip()
      rid = "%s_%s_%s"%(resname, cid, resseq)
      # Write as mmCIF
      file_name = "%s_%s.cif" % (rid, code)
      assert not os.path.isfile(file_name)
      with open(file_name, "w") as fo:
        fo.write(model_box.model_as_mmcif())
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
  print("processed (will be skept):", len(processed))
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
