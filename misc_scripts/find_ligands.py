from __future__ import absolute_import, division, print_function
import os, sys
import traceback
from io import StringIO
from collections import OrderedDict
from scitbx.array_family import flex
from libtbx import easy_mp
from libtbx.utils import Sorry, null_out
from cctbx import crystal, maptbx, uctbx
import cctbx.miller
import iotbx.pdb
from iotbx import reflection_file_reader
import mmtbx.model
import mmtbx.monomer_library.server
mon_lib_srv = mmtbx.monomer_library.server.server()

cif_files = os.getenv('PDB_MIRROR_MMCIF')
pdbmtz = os.getenv('PDBMTZ')
SUPPORTED = ["H","C","N","O","S","SE","B","SI","P","AS","F","CL","BR","I"]


def get_model_files_dict(path):
  """
  PDB code -> model file, from the INDEX file of the mirror. The code is the
  file name without directories and extensions (ab/1abc.cif.gz -> 1abc).
  """
  result = {}
  with open(os.path.join(path, "INDEX"), "r") as ifn:
    for l in ifn.readlines():
      l = l.strip()
      if len(l) == 0: continue
      code = os.path.basename(l).split(".")[0].lower()
      result[code] = os.path.join(path, l)
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
    file_name = os.path.join(path, l)
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

def has_ligands(h, code, log):
  get_class = iotbx.pdb.common_residue_names_get_class
  skip_classes = {
      "common_water",
      "common_amino_acid",
      "modified_amino_acid",
      "d_amino_acid",
      "common_rna_dna",
      "modified_rna_dna",
      "ccp4_mon_lib_rna_dna"
  }
  ligand_selections = []
  for rg in h.residue_groups():
    # Skip residues with alternative conformations
    if len(rg.atom_groups()) != 1 or rg.have_conformers(): continue
    ag = rg.only_atom_group()
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
    print(code, ag.resname, nat, file=log)
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

def check_maps2(tfofc, fofc, sites_cart, unit_cell, log):
  m1 = flex.double()
  m2 = flex.double()
  for site_cart in sites_cart:
    site_frac = unit_cell.fractionalize(site_cart)
    mv_tfofc = tfofc.tricubic_interpolation(site_frac)
    mv_fofc  = fofc .tricubic_interpolation(site_frac)
    m1.append(mv_tfofc)
    m2.append(mv_fofc)
  print("  min tfofc, max |fofc|:", flex.min(m1), flex.max(flex.abs(m2)),
    file=log)
  if flex.min(m1) < 0.5 or flex.max(flex.abs(m2)) > 3.0:
    return False
  return True

def append_chains_merged(hierarchy, model, new_chain_id):
  """
  Append the residue groups of hierarchy (one model) to model, one new chain
  per chain ID: segments of a chain (split at TER in the input) are merged
  and their residue groups sorted by residue number and insertion code.
  The new chain is named new_chain_id(chain_id).
  """
  for chain_id in hierarchy.chain_ids(unique_only=True):
    rgs = []
    for chain in hierarchy.only_model().chains():
      if chain.id == chain_id:
        rgs.extend(chain.residue_groups())
    rgs.sort(key=lambda rg: (rg.resseq_as_int(), rg.icode))
    new_chain = iotbx.pdb.hierarchy.chain(id=new_chain_id(chain_id))
    for rg in rgs:
      new_chain.append_residue_group(rg.detached_copy())
    model.append_chain(new_chain)

def extract_ligand_environment(pdb_hierarchy, ligand_iselection,
                               crystal_symmetry, distmax):
  """
  Extracts a ligand and its surrounding environment: a new hierarchy with
  the ligand and every complete residue that has at least one atom within
  `distmax` of a ligand atom, in the asymmetric unit or in a symmetry
  mate.

  Residues of the asymmetric unit keep their chain IDs, one chain per ID,
  sorted by residue number.  Symmetry mates
  (iotbx.pdb.hierarchy.root.sym_equiv_residue_groups) are added as new
  chains, one per original chain and symmetry operation, named with the
  prefix 'S' (chain 'A' becomes 'SA'); a numerical suffix ('SA1', 'SA2') is
  appended only if the same original chain interacts with the ligand
  through several distinct symmetry operations.

  Returns None if any of the extracted residues has alternative
  conformations.
  """
  if isinstance(ligand_iselection, flex.bool):
    ligand_iselection = ligand_iselection.iselection()
  atoms = pdb_hierarchy.atoms()
  sites_cart = atoms.extract_xyz()
  # Atoms of the asymmetric unit within distmax of the ligand, without
  # symmetry mates, plus the ligand itself
  keep = flex.bool(atoms.size(), False)
  keep.set_selected(ligand_iselection, True)
  for i, j, rt_mx, distance in crystal.pairs_within(
      sites_cart=sites_cart, distance_cutoff=distmax, i_seqs=ligand_iselection):
    keep[j] = True
  # Complete the residue groups, by position in atoms()
  first = 0
  for rg in pdb_hierarchy.residue_groups():
    i_seqs = flex.size_t_range(first, first + rg.atoms_size())
    if keep.select(i_seqs).count(True) > 0:
      if rg.have_conformers(): return None
      keep.set_selected(i_seqs, True)
    first += rg.atoms_size()
  # Symmetry mates of residues within distmax of the ligand
  copies = pdb_hierarchy.sym_equiv_residue_groups(
    crystal_symmetry = crystal_symmetry,
    radius           = distmax,
    selection        = ligand_iselection)
  for copy in copies:
    for rg in copy.hierarchy.residue_groups():
      if rg.have_conformers(): return None
  # Construct the new hierarchy
  new_hierarchy = iotbx.pdb.hierarchy.root()
  new_model = iotbx.pdb.hierarchy.model()
  new_hierarchy.append_model(new_model)
  asu = pdb_hierarchy.select(keep, copy_atoms=True)
  append_chains_merged(asu, new_model, lambda chain_id: chain_id)
  used_chain_ids = set(new_hierarchy.chain_ids())
  def sym_chain_id(chain_id):
    base_new_id = "S" + chain_id.strip()
    new_id = base_new_id
    suffix = 1
    while new_id in used_chain_ids:
      new_id = base_new_id + str(suffix)
      suffix += 1
    used_chain_ids.add(new_id)
    return new_id
  for copy in copies:
    append_chains_merged(copy.hierarchy, new_model, sym_chain_id)
  new_hierarchy.atoms().reset_i_seq()
  return new_hierarchy

def get_map_coefficients(miller_arrays):
  mc_tfofc, mc_fofc = None, None
  for ma in miller_arrays:
    lbls = ma.info().label_string()
    if lbls == "2FOFCWT,PHI2FOFCWT": mc_tfofc = ma.deep_copy()
    if lbls == "FOFCWT,PHIFOFCWT":   mc_fofc  = ma.deep_copy()
  assert mc_tfofc is not None and mc_fofc is not None
  return mc_tfofc, mc_fofc

def get_maps(mc_tfofc, mc_fofc):
  def _get_real_map(mc):
    fft_map = cctbx.miller.fft_map(
      crystal_gridding     = crystal_gridding,
      fourier_coefficients = mc)
    fft_map.apply_sigma_scaling()
    return fft_map.real_map_unpadded()
  crystal_gridding = mc_tfofc.crystal_gridding(
    d_min             = mc_tfofc.d_min(),
    symmetry_flags    = maptbx.use_space_group_symmetry,
    resolution_factor = 1./3)
  tfofc = _get_real_map(mc = mc_tfofc)
  fofc  = _get_real_map(mc = mc_fofc)
  return tfofc, fofc

def write_skip(code, msg):
  file_name = "skip_%s.log"%code
  with open(file_name, "w") as fo:
    fo.write(msg)

def write_info(code, log):
  """Diagnostics of one entry, written only if there are any."""
  msg = log.getvalue()
  if len(msg) == 0: return
  file_name = "info_%s.log"%code
  with open(file_name, "w") as fo:
    fo.write(msg)

def run_one(args):
  cif, hkl, code = args
  log = StringIO()
  try:
    _run_one(cif=cif, hkl=hkl, code=code, log=log)
  except Exception as e:
    msg = "Other:\n %s"%traceback.format_exc()
    write_skip(code=code, msg=msg)
  write_info(code=code, log=log)

def _run_one(cif, hkl, code, log):
  #
  # DATA
  #
  miller_arrays = reflection_file_reader.any_reflection_file(file_name =
    hkl).as_miller_arrays()
  mc_tfofc, mc_fofc = get_map_coefficients(miller_arrays = miller_arrays)
  if mc_tfofc.d_min() > 2.5:
    write_skip(code=code, msg="Bad resolution")
    return
  #
  # MODEL
  #
  model = mmtbx.model.manager(
    model_input         = iotbx.pdb.input(file_name=cif),
    stop_for_unknowns   = False,
    skip_ss_annotations = True,
    process_biomt       = False)
  cs = model.crystal_symmetry()
  model = model.remove_hydrogens()
  h = model.get_hierarchy()
  atoms = h.atoms()
  atoms.reset_i_seq()
  sp = model.selection(string="protein")
  if sp.count(True)*100./sp.size() < 80.:
    write_skip(code=code, msg="Not a protein")
    return # not protein
  if len(h.models()) != 1:
    write_skip(code=code, msg="Multi-model")
    return
  ligand_selections = has_ligands(h, code, log)
  if len(ligand_selections)==0:
    write_skip(code=code, msg="No ligand(s) selected")
    return
  #
  # Compute maps
  #
  tfofc, fofc = get_maps(mc_tfofc = mc_tfofc, mc_fofc = mc_fofc)
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
    write_skip(code=code, msg="No ligands with good map")
    return
  file_was_created = False
  for i, ls in enumerate(ligand_selections):
    # A failure for one ligand does not skip the entry: the files of the
    # other ligands are still written
    try:
      file_was_created |= process_ligand(
        h=h, ls=ls, i=i, cs=cs, tfofc=tfofc, fofc=fofc, code=code, log=log)
    except KeyboardInterrupt:
      raise
    except Exception as e:
      print(i, "failed:\n%s" % traceback.format_exc(), file=log)
  if not file_was_created:
    write_skip(code=code, msg="File was not created")

def process_ligand(h, ls, i, cs, tfofc, fofc, code, log):
  """
  Extract, check and write the environment of one ligand. Returns True if
  the file was written.
  """
  around_ligand = extract_ligand_environment(
    pdb_hierarchy     = h,
    ligand_iselection = ls,
    crystal_symmetry  = cs,
    distmax           = 4.0)
  if around_ligand is None:
    print(i, "poor env (altlocs)", file=log)
    return False
  # Check for unsupported
  elements  = [a.element.strip().upper() for a in around_ligand.atoms()]
  found_unsupported = False
  for e in elements:
    if not e in SUPPORTED: found_unsupported = True
  if found_unsupported:
    print(i, "poor env (unsupp)", file=log)
    return False
  # Check whole thing against maps
  mcheck = check_maps2(
    tfofc      = tfofc,
    fofc       = fofc,
    sites_cart = around_ligand.atoms().extract_xyz(),
    unit_cell  = cs.unit_cell(),
    log        = log)
  if mcheck is False:
    print(i, "poor env (map)", file=log)
    return False
  # Move the environment to the center of a P1 box
  box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
    sites_cart   = around_ligand.atoms().extract_xyz(),
    buffer_layer = 5.0)
  around_ligand.atoms().set_xyz(box.sites_cart)
  # Create model object
  try:
    model_box = mmtbx.model.manager(
      pdb_hierarchy    = around_ligand,
      crystal_symmetry = box.crystal_symmetry(),
      log              = null_out())
    model_box.process(make_restraints=True)
  except Sorry as e:
    print(i, "model processing failed:", str(e), file=log)
    return False
  # Check geometry
  #g = model_box.geometry_statistics().result(slim = True)
  #
  #print("     bond ", g.bond.mean,  g.bond.max)
  #print("     angle", g.angle.mean, g.angle.max)
  #print("     nb", g.nonbonded.min)

  # Format output file name prefix
  atoms = h.atoms()
  resname = atoms[ls[0]].parent().resname.strip()
  resseq  = atoms[ls[0]].parent().parent().resseq.strip()
  cid     = atoms[ls[0]].parent().parent().parent().id.strip()
  rid = "%s_%s_%s"%(resname, cid, resseq)
  # Write as mmCIF
  file_name = "%s_%s.cif" % (rid, code)
  assert not os.path.isfile(file_name)
  with open(file_name, "w") as fo:
    fo.write(model_box.model_as_mmcif())
  return True

def run(cmdargs, NPROC=128):
  #
  processed = []
  skept     = []
  for f in os.listdir("."):
    code = f.rsplit("_", 1)[-1].split(".")[0]
    if(f.endswith(".cif")):
      processed.append(code)
    if(f.endswith(".log") and f.startswith("skip_")):
      skept.append(code)
  ignore = processed + skept
  #
  print("processed (will be skept):", len(processed))
  print("skept     (will be skept):", len(skept))
  print("ignored   (will be skept):", len(ignore))
  #
  cifs = get_model_files_dict(path=cif_files)
  hkls = get_mtz_dict(path=pdbmtz)
  print("Number of cifs:", len(cifs.keys()))
  print("Number of hkls:", len(hkls.keys()))
  argss = []
  for code in hkls.keys():
    if(code in ignore): continue
    hkl,size = hkls[code]
    if code not in cifs: continue
    argss.append([cifs[code], hkl, code])
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
