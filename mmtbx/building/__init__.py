
from __future__ import absolute_import, division, print_function
from libtbx import group_args
import sys
import mmtbx.refinement.real_space
from six.moves import zip
from six.moves import range

#-----------------------------------------------------------------------
# MODEL UTILITIES
#-----------------------------------------------------------------------

def reprocess_pdb(pdb_hierarchy, crystal_symmetry, cif_objects, out):
  from mmtbx.monomer_library import pdb_interpretation
  from mmtbx.monomer_library import server
  mon_lib_srv = server.server()
  ener_lib = server.ener_lib()
  for cif_object in cif_objects :
    for srv in [mon_lib_srv, ener_lib]:
      srv.process_cif_object(cif_object=cif_object)
  return pdb_interpretation.process(
    mon_lib_srv=mon_lib_srv,
    ener_lib=ener_lib,
    raw_records=pdb_hierarchy.as_pdb_or_mmcif_string(
      target_format='pdb',crystal_symmetry=crystal_symmetry),
    crystal_symmetry=crystal_symmetry,
    log=out)

def get_nearby_water_selection(
    pdb_hierarchy,
    xray_structure,
    selection,
    selection_buffer_radius=5):
  sel_cache = pdb_hierarchy.atom_selection_cache()
  water_selection = sel_cache.selection("resname HOH")
  selection_within = xray_structure.selection_within(
    radius=selection_buffer_radius,
    selection=selection)
  return (selection_within & water_selection)

def iter_residue_groups(hierarchy):
  for chain in hierarchy.only_model().chains():
    for residue_group in chain.residue_groups():
      yield residue_group

def extract_iselection(pdb_objects):
  from scitbx.array_family import flex
  isel = flex.size_t()
  for pdb_obj in pdb_objects :
    isel.extend(pdb_obj.atoms().extract_i_seq())
  assert not isel.all_eq(0)
  return isel

def remove_sidechain_atoms(pdb_objects):
  for pdb_obj in pdb_objects :
    atoms = pdb_obj.atoms()
    for atom in atoms :
      if (not atom.name.strip() in ["CA","C","N","O","CB"]):
        atom.parent().remove_atom(atom)

def is_stub_residue(atom_group):
  if (atom_group.resname in ["GLY","ALA"]):
    return True
  has_sidechain_atoms = False
  for atom in atom_group.atoms():
    if (not atom.name.strip() in ["H","HA","C","CA","CB","N","O",]):
      has_sidechain_atoms = True
  return not has_sidechain_atoms

def get_non_hydrogen_atom_indices(pdb_object):
  from scitbx.array_family import flex
  i_seqs_no_hd = flex.size_t()
  for atom in pdb_object.atoms():
    if (atom.element.strip() not in ["H","D"]):
      i_seqs_no_hd.append(atom.i_seq)
  return i_seqs_no_hd

def get_window_around_residue(residue, window_size=2):
  residue_group = residue.parent()
  chain = residue_group.parent()
  all_residues = chain.residue_groups()
  sel_residues = []
  for i_res, other_rg in enumerate(all_residues):
    if (other_rg == residue_group):
      for j_res in range(-window_size, window_size+1):
        k_res = i_res + j_res
        if (k_res >= 0) and (k_res < len(all_residues)):
          sel_residues.append(all_residues[k_res])
  assert (len(sel_residues) > 0)
  return sel_residues

def atom_group_as_hierarchy(atom_group):
  """
  Convert an iotbx.pdb.hierarchy.atom_group object to a separate hierarchy
  (leaving the original hierarchy untouched), which allows us to use the
  pickling capabilities of iotbx.pdb.hierarchy.root.
  """
  import iotbx.pdb.hierarchy
  old_rg = atom_group.parent()
  assert (old_rg is not None)
  old_chain = old_rg.parent()
  root = iotbx.pdb.hierarchy.root()
  model = iotbx.pdb.hierarchy.model()
  chain = iotbx.pdb.hierarchy.chain(id=old_chain.id)
  residue_group = iotbx.pdb.hierarchy.residue_group(
    resseq=old_rg.resseq,
    icode=old_rg.icode)
  residue_group.append_atom_group(atom_group.detached_copy())
  chain.append_residue_group(residue_group)
  model.append_chain(chain)
  root.append_model(model)
  atoms = root.atoms()
  atoms.reset_i_seq()
  atoms.reset_serial()
  return root

def residues_are_adjacent(residue1, residue2, max_sep=2.5):
  if (type(residue1).__name__ == "residue_group"):
    residue1 = residue1.atom_groups()[0]
  if (type(residue2).__name__ == "residue_group"):
    residue2 = residue2.atom_groups()[0]
  c_pep = n_pep = None
  for atom in residue1.atoms():
    if (atom.name.strip() == "C"):
      c_pep = atom
      break
  for atom in residue2.atoms():
    if (atom.name.strip() == "N"):
      n_pep = atom
      break
  if (None in [c_pep, n_pep]):
    return False
  return c_pep.distance(n_pep) < max_sep

#-----------------------------------------------------------------------
# MAP-RELATED
#-----------------------------------------------------------------------

def residue_density_quality(
    atom_group,
    unit_cell,
    two_fofc_map,
    fofc_map):
  from scitbx.array_family import flex
  sites_cart = flex.vec3_double()
  atom_names = flex.std_string()
  for atom in atom_group.atoms():
    if (not atom.element.strip() in ["H","D"]):
      sites_cart.append(atom.xyz)
      atom_names.append(atom.name)
  if (len(sites_cart) == 0):
    raise RuntimeError("No non-hydrogen atoms in %s" % atom_group.id_str())
  density = local_density_quality(
    fofc_map=fofc_map,
    two_fofc_map=two_fofc_map,
    sites_cart=sites_cart,
    atom_names=atom_names,
    unit_cell=unit_cell)
  return density

def get_difference_maps(fmodel, resolution_factor=0.25):
  fofc_coeffs = fmodel.map_coefficients(map_type="mFo-DFc",
    exclude_free_r_reflections=True)
  fofc_fft = fofc_coeffs.fft_map(
    resolution_factor=resolution_factor)
  fofc_map = fofc_fft.apply_sigma_scaling().real_map_unpadded()
  two_fofc_coeffs = fmodel.map_coefficients(map_type="2mFo-DFc",
    exclude_free_r_reflections=True)
  two_fofc_fft = two_fofc_coeffs.fft_map(resolution_factor=resolution_factor)
  two_fofc_map = two_fofc_fft.apply_sigma_scaling().real_map_unpadded()
  return two_fofc_map, fofc_map

# XXX this should probably be merged with something else, e.g. the
# local_density_quality class
def get_model_map_stats(
    selection,
    target_map,
    model_map,
    unit_cell,
    sites_cart,
    pdb_atoms,
    local_sampling=False):
  """
  Collect basic statistics for a model map and some target map (usually an
  mFo-DFc map), including CC, mean, and minimum density at the atomic
  positions.
  """
  assert (len(target_map) == len(model_map))
  iselection = selection
  if (type(selection).__name__ == 'bool'):
    iselection = selection.iselection()
  from scitbx.array_family import flex
  sites_cart_refined = sites_cart.select(selection)
  sites_selected = flex.vec3_double()
  map1 = flex.double()
  map2 = flex.double()
  min_density = sys.maxsize
  sum_density = n_sites = 0
  worst_atom = None
  # XXX I'm not sure the strict density cutoff is a good idea here
  for i_seq, xyz in zip(iselection, sites_cart_refined):
    if (pdb_atoms[i_seq].element.strip() != "H"):
      sites_selected.append(xyz)
      site_frac = unit_cell.fractionalize(site_cart=xyz)
      target_value = target_map.tricubic_interpolation(site_frac)
      if (target_value < min_density):
        min_density = target_value
        worst_atom = pdb_atoms[i_seq]
      sum_density += target_value
      n_sites += 1
      if (not local_sampling):
        map1.append(target_value)
        map2.append(model_map.tricubic_interpolation(site_frac))
  assert (n_sites > 0)
  if (local_sampling):
    from cctbx import maptbx
    map_sel = maptbx.grid_indices_around_sites(
      unit_cell=unit_cell,
      fft_n_real=target_map.focus(),
      fft_m_real=target_map.all(),
      sites_cart=sites_selected,
      site_radii=flex.double(sites_selected.size(), 1.0))
    map1 = target_map.select(map_sel)
    map2 = model_map.select(map_sel)
  assert (len(map1) > 0) and (len(map1) == len(map2))
  cc = flex.linear_correlation(x=map1, y=map2).coefficient()
  return group_args(
    cc=cc,
    min=min_density,
    mean=sum_density/n_sites)

#-----------------------------------------------------------------------
# SAMPLING UTILITIES
#-----------------------------------------------------------------------
def generate_sidechain_clusters(residue, mon_lib_srv):
  """
  Extract Chi angle indices (including rotation axis) from the atom_group
  """
  from mmtbx.utils import rotatable_bonds
  from scitbx.array_family import flex
  atoms = residue.atoms()
  axes_and_atoms_aa_specific = \
      rotatable_bonds.axes_and_atoms_aa_specific(residue = residue,
        mon_lib_srv = mon_lib_srv)
  result = []
  if(axes_and_atoms_aa_specific is not None):
    for i_aa, aa in enumerate(axes_and_atoms_aa_specific):
      n_heavy = 0
      #for i_seq in aa[1] :
      #  if (atoms[i_seq].element.strip() != "H"):
      #    n_heavy += 1
      #if (n_heavy == 0) : continue
      if(i_aa == len(axes_and_atoms_aa_specific)-1):
        result.append(mmtbx.refinement.real_space.cluster(
          axis=aa[0],
          atoms_to_rotate=aa[1],
          selection=flex.size_t(aa[1]),
          vector=None)) # XXX
      else:
        result.append(mmtbx.refinement.real_space.cluster(
          axis=aa[0],
          atoms_to_rotate=aa[1],
          selection=flex.size_t([aa[1][0]]),
          vector=None)) # XXX
  return result
