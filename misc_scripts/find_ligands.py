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
    #if code != "5nmw": continue
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

import iotbx.pdb
from cctbx import crystal, sgtbx
from scitbx.array_family import flex

def extract_ligand_environment(pdb_hierarchy, ligand_iselection, crystal_symmetry, distmax):
  """
  Extracts a ligand and its surrounding environment, explicitly applying
  crystal symmetry to map interacting symmetry-related residues next to the ligand.

  This function creates a new PDB hierarchy containing the specified ligand
  and any complete residue that has at least one atom within `distmax` of
  any ligand atom. If an interacting residue comes from a symmetry-related
  molecule, the required rigid-body transformation is explicitly applied to
  its coordinates.

  To avoid duplicate atom-label errors, symmetry-transformed chains are
  assigned new, unique chain IDs starting with the prefix 'S' (e.g., chain 'A'
  becomes 'SA').

  Importantly, if consecutive residues from the same original chain are mapped
  by the *same* symmetry operation, they are kept together as consecutive
  residues within the *same* single chain ID (e.g., 'SA'). Numerical suffixes
  (e.g., 'SA1', 'SA2') are only appended if the *same original chain* interacts
  with the ligand via *multiple, distinct* symmetry operations. This isolates
  distinct symmetry mappings of the same chain into separate new chains to
  prevent unphysical spatial jumps and duplicate residue labels.

  If any of the extracted residues contain alternative conformations (altlocs),
  the extraction is aborted.
  """
  # Ensure ligand_iselection is a standard set of integers (i_seqs)
  if isinstance(ligand_iselection, flex.bool):
      ligand_iselection = ligand_iselection.iselection()
  ligand_set = set(ligand_iselection)

  # Map each i_seq to its corresponding residue_group and chain
  iseq_to_rg = {}
  rg_to_chain = {}
  for model in pdb_hierarchy.models():
      for chain in model.chains():
          for rg in chain.residue_groups():
              for atom in rg.atoms():
                  iseq_to_rg[atom.i_seq] = rg
              rg_to_chain[rg] = chain

  # Helper function to check if a residue group has alternative conformations
  def has_altlocs(rg):
      if len(rg.atom_groups()) > 1:
          return True
      for ag in rg.atom_groups():
          if ag.altloc.strip() != '':
              return True
      return False

  sites_cart = pdb_hierarchy.atoms().extract_xyz()

  # Setup symmetry and map atoms across unit cell boundaries
  sst = crystal_symmetry.special_position_settings().site_symmetry_table(sites_cart=sites_cart)
  conn_asu_mappings = crystal_symmetry.special_position_settings().asu_mappings(buffer_thickness=distmax)
  conn_asu_mappings.process_sites_cart(
      original_sites=sites_cart,
      site_symmetry_table=sst
  )

  pair_generator = crystal.neighbors_fast_pair_generator(
      conn_asu_mappings,
      distance_cutoff=distmax)

  # Dictionary to store required RGs grouped by transformation and chain:
  # extracted[op_str][chain_id] = set(rgs)
  extracted = {}
  identity_op = sgtbx.rt_mx("x,y,z").as_xyz()

  # Step 1. Always include the full ligand in its original position
  for iseq in ligand_set:
      rg = iseq_to_rg[iseq]
      chain_id = rg_to_chain[rg].id
      extracted.setdefault(identity_op, {}).setdefault(chain_id, set()).add(rg)

  # Step 2. Collect all interacting residues, mapped back to the ligand's environment
  for pair in pair_generator:
      i_seq = pair.i_seq
      j_seq = pair.j_seq

      rt_mx_i = conn_asu_mappings.get_rt_mx_i(pair)
      rt_mx_j = conn_asu_mappings.get_rt_mx_j(pair)

      if i_seq in ligand_set:
          op = rt_mx_i.inverse().multiply(rt_mx_j)
          rg = iseq_to_rg[j_seq]
          chain_id = rg_to_chain[rg].id
          extracted.setdefault(op.as_xyz(), {}).setdefault(chain_id, set()).add(rg)

      if j_seq in ligand_set:
          op = rt_mx_j.inverse().multiply(rt_mx_i)
          rg = iseq_to_rg[i_seq]
          chain_id = rg_to_chain[rg].id
          extracted.setdefault(op.as_xyz(), {}).setdefault(chain_id, set()).add(rg)

  # Step 3. Validate that none of the extracted residues have alternative conformations
  for op_str, chain_dict in extracted.items():
      for chain_id, rgs in chain_dict.items():
          for rg in rgs:
              if has_altlocs(rg):
                  return None

  # Step 4. Construct the new hierarchy
  new_hierarchy = iotbx.pdb.hierarchy.root()
  new_model = iotbx.pdb.hierarchy.model()
  new_hierarchy.append_model(new_model)

  used_chain_ids = set()

  # 4a. Add residues from the identity operation first (original chains stay together)
  if identity_op in extracted:
      chain_dict = extracted[identity_op]
      for chain_id in sorted(chain_dict.keys()):
          new_chain = iotbx.pdb.hierarchy.chain(id=chain_id)
          used_chain_ids.add(chain_id)

          # Sort RGs by residue sequence & insert code to maintain consecutive order
          rgs = sorted(list(chain_dict[chain_id]), key=lambda x: (x.resseq_as_int(), x.icode))
          for rg in rgs:
              new_chain.append_residue_group(rg.detached_copy())
          new_model.append_chain(new_chain)

  # 4b. Add residues from explicit symmetry operations
  unit_cell = crystal_symmetry.unit_cell()

  for op_str, chain_dict in sorted(extracted.items()):
      if op_str == identity_op:
          continue

      rt_mx = sgtbx.rt_mx(op_str)

      for chain_id in sorted(chain_dict.keys()):
          # Create a unique chain ID starting with 'S'
          base_new_id = "S" + chain_id.strip()
          new_id = base_new_id
          suffix = 1
          while new_id in used_chain_ids:
              new_id = base_new_id + str(suffix)
              suffix += 1
          used_chain_ids.add(new_id)

          new_chain = iotbx.pdb.hierarchy.chain(id=new_id)

          # Sorted iteration guarantees that sequential pairs remain adjacent
          rgs = sorted(list(chain_dict[chain_id]), key=lambda x: (x.resseq_as_int(), x.icode))

          for rg in rgs:
              new_rg = rg.detached_copy()

              # Apply the symmetry operator to rotate/translate the rigid-body residue
              for atom in new_rg.atoms():
                  site_frac = unit_cell.fractionalize(atom.xyz)
                  new_site_frac = rt_mx * site_frac
                  atom.xyz = unit_cell.orthogonalize(new_site_frac)

              new_chain.append_residue_group(new_rg)
          new_model.append_chain(new_chain)

  # Clean indices
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
        pdb_hierarchy     = h,
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
