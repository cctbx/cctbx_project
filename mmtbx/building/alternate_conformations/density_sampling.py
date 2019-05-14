
from __future__ import absolute_import, division, print_function
import libtbx.phil
import time
import sys

sidechain_density_params = """
  min_2fofc = 0.8
    .type = float
    .optional = False
  min_fofc = 2.0
    .type = float
    .optional = False
  min_fofc_only = 4.0
    .type = float
    .optional = False
  mean_2fofc = 1.0
    .type = float
    .optional = False
  mean_fofc = 3.0
    .type = float
    .optional = False
  min_cbeta_2fofc = 2.0
    .type = float
"""

master_params = """
set_partial_occupancy = None
  .type = float
rmsd_min = 1.0
  .type = float
rigid_body_refine = False
  .type = bool
translation = 0.5 #None
  .type = float
  .help = Specifies the radius for the optional translation search, in \
    Angstroms
translation_sampling = 10
  .type = int
  .optional = False
  .help = Number of distances within the translation search to sample.
torsion_search {
  include scope mmtbx.building.alternate_conformations.conformer_generation.torsion_search_params
}
sidechain_scoring {
  %s
}
""" % sidechain_density_params

def master_phil():
  return libtbx.phil.parse(master_params, process_includes=True)

def screen_residue(
    residue,
    prev_residue,
    next_residue,
    next_next_residue,
    sites_cart,
    fmodel,
    mon_lib_srv,
    params,
    fofc_map=None,
    two_fofc_map=None,
    map_file_name=None,
    backrub=False,
    shear=False,
    verbose=False,
    out=None):
  from mmtbx.building.alternate_conformations import conformer_generation
  if (out is None):
    out = sys.stdout
  if (params is None):
    params = master_phil().extract()
  xrs = fmodel.xray_structure
  occ_start = xrs.scatterers().extract_occupancies()
  i_seqs = residue.atoms().extract_i_seq()
  assert (not i_seqs.all_eq(0))
  if (params.set_partial_occupancy is not None):
    occ = occ_start.deep_copy()
    occ.set_selected(i_seqs, params.set_partial_occupancy)
    xrs.set_occupancies(occ)
    fmodel.update_xray_structure(xrs, update_f_calc=True)
  if (fofc_map is None) or (params.set_partial_occupancy is not None):
    two_fofc_coeffs = fmodel.map_coefficients(map_type="2mFo-DFc")
    fofc_coeffs = fmodel.map_coefficients(map_type="mFo-DFc")
    two_fofc_fft = two_fofc_coeffs.fft_map(resolution_factor=1/4.)
    fofc_fft = fofc_coeffs.fft_map(resolution_factor=1/4.)
    two_fofc_map = two_fofc_fft.apply_sigma_scaling().real_map_unpadded()
    fofc_map = fofc_fft.apply_sigma_scaling().real_map_unpadded()
    if (map_file_name is not None):
      import iotbx.map_tools
      iotbx.map_tools.write_map_coeffs(two_fofc_coeffs, fofc_coeffs,
        map_file_name)
  unit_cell = xrs.unit_cell()
  sites_start = sites_cart.deep_copy()
  good_confs = []
  n_confs = 0
  t_start = time.time()
  for conf in conformer_generation.generate_single_residue_confs(
      atom_group=residue,
      sites_cart=sites_cart.deep_copy(),
      mon_lib_srv=mon_lib_srv,
      params=params.torsion_search,
      prev_residue=prev_residue,
      next_residue=next_residue,
      next_next_residue=next_next_residue,
      backrub=backrub,
      shear=shear):
    n_confs += 1
    # basically a very crude rigid-body grid search
    translations = []
    for sites_new, translation in translation_search(
        sites=conf.sites_selected(),
        translation=params.translation,
        translation_sampling=params.translation_sampling):
      conf_trans = conf.translate_sites(
        sites_new=sites_new,
        translation=translation).update_rmsd(sites_start)
      good_maps = score_density(
        self=conf_trans,
        two_fofc_map=two_fofc_map,
        fofc_map=fofc_map,
        unit_cell=unit_cell,
        params=params.sidechain_scoring)
      translations.append(conf_trans)
      if (good_maps) : #and (conf_trans.rmsd > params.rmsd_min):
        translations.append(conf_trans)
    if (len(translations) > 0):
      translations.sort(lambda a,b: cmp(b.mean_fofc, a.mean_fofc))
      good_confs.append(translations[0])
  if (n_confs > 1) and (len(good_confs) != 0):
    for conf in good_confs :
      conf.show_summary(out=out)
  xrs.set_occupancies(occ_start)
  t_end = time.time()
  if (verbose):
    print("time to sample residue: %.2fs" % (t_end - t_start), file=out)
  return good_confs

def score_density(self, two_fofc_map, fofc_map, unit_cell, params,
    sidechain_only=True):
  if (two_fofc_map is None):
    return None
  backbone_atoms = ["N","CA","C","O","H"]
  if (sidechain_only):
    sum_fofc = 0
    sum_2fofc = 0
    min_fofc = sys.maxsize
    min_2fofc = sys.maxsize
    min_cbeta_2fofc = sys.maxsize
    min_2fofc_atom_name = None
    n_atoms = 0
    for atom in self.residue.atoms():
      name = atom.name.strip()
      if (name in backbone_atoms) or (atom.element.strip() == "H"):
        continue
      site_frac = unit_cell.fractionalize(site_cart=self.sites_cart[atom.i_seq])
      value_fofc = fofc_map.eight_point_interpolation(site_frac)
      value_2fofc = two_fofc_map.eight_point_interpolation(site_frac)
      if (name == "CB"):
        if (value_2fofc < min_cbeta_2fofc):
          min_cbeta_2fofc = value_2fofc
      sum_fofc += value_fofc
      sum_2fofc += value_2fofc
      if (value_fofc < min_fofc):
        min_fofc = value_fofc
      if (value_2fofc < min_2fofc):
        min_2fofc = value_2fofc
        min_2fofc_atom_name = atom.name
      n_atoms += 1
    if (n_atoms != 0):
      self.mean_fofc = sum_fofc / n_atoms
      self.mean_2fofc = sum_2fofc / n_atoms
  return ((self.mean_fofc >= params.mean_fofc) and
          (((self.mean_2fofc >= params.mean_2fofc) and
            (min_2fofc >= params.min_2fofc) and
            (min_fofc >= params.min_fofc)) or
           (min_fofc >= params.min_fofc_only)))# and
          #(min_cbeta_2fofc >= params.min_cbeta_2fofc))

# XXX from mmtbx/refinement/real_space/fit_residue.py
def rigid_body_refine(residue, target_map, unit_cell, max_iterations = 250):
  assert 0, "Never tested, now definitely broken."
  import mmtbx.geometry_restraints
  import mmtbx.refinement.real_space.rigid_body
  import cctbx.geometry_restraints.manager
  from scitbx.array_family import flex
  import scitbx.lbfgs
  reference_sites = flex.vec3_double()
  reference_selection = flex.size_t()
  cntr = 0
  for atom in residue.atoms():
    if(atom.name.strip().upper() in ["N", "C"]):
      reference_sites.append(atom.xyz)
      reference_selection.append(cntr)
      cntr += 1
  generic_restraints_manager = mmtbx.geometry_restraints.manager()
  restraints_manager = cctbx.geometry_restraints.manager.manager(
    generic_restraints_manager = generic_restraints_manager)
  restraints_manager.generic_restraints_manager.reference_manager.\
    add_coordinate_restraints(
      sites_cart = reference_sites,
      selection  = reference_selection,
      sigma      = 0.5)
  flags = cctbx.geometry_restraints.flags.flags(generic_restraints=True)
  lbfgs_termination_params=scitbx.lbfgs.termination_parameters(
    max_iterations = max_iterations)
  minimized = mmtbx.refinement.real_space.rigid_body.refine(
    residue                     = residue,
    density_map                 = target_map,
    geometry_restraints_manager = restraints_manager,
    real_space_target_weight    = 1,
    real_space_gradients_delta  = 2.0*0.25,
    lbfgs_termination_params    = lbfgs_termination_params,
    unit_cell                   = unit_cell,
    cctbx_geometry_restraints_flags = flags)
  return minimized.sites_cart_residue

#def translation_search(sites, translation=None, translation_sampling=1):
def translation_search(sites, translation, translation_sampling):
  from scitbx.array_family import flex
  results = []
  if (translation is None):
    results.append((sites, None))
  else :
    assert (translation_sampling >= 1)
    results.append((sites.deep_copy(), None))
    increment = translation / translation_sampling
    for i in range(3):
      for j in [-1, 1] :
        for k in range(translation_sampling):
          t = [ 0 ] * 3
          t[i] = j * (k+1) * increment
          tvec = flex.vec3_double(sites.size(), t)
          trans_sites = sites.deep_copy() + tvec
          assert (len(trans_sites) == len(sites))
          results.append((trans_sites, t))
  return results
