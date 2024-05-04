from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from cctbx import geometry_restraints
from libtbx.utils import Sorry
import boost_adaptbx.boost.python as bp
ext = bp.import_ext("mmtbx_reference_coordinate_ext")
from mmtbx.rotamer.sidechain_angles import collect_residue_torsion_angles


def generate_torsion_restraints(
      pdb_hierarchy,
      sites_cart,
      selection=None,
      sigma=2.5,
      limit=15.0,
      chi_angles_only=False,
      top_out_potential=False,
      origin_id=None):
  """
  selection ties up hierarchy with sites cart. Basically, if applied
    to hierarchy, remaining atoms should correspond to sites_cart.
    Therefore it is necessary that
    len(sites_cart) == len(actual_selection) below.
    This seems to be done this way to make it possible to use sites_cart
    from another source (reference model) which is not necessary
    of the same size as hierarchy.
  """
  pdb_hierarchy.atoms().reset_i_seq()
  torsion_proxies = geometry_restraints.shared_dihedral_proxy()
  if pdb_hierarchy.atoms_size() < 4:
    return torsion_proxies
  assert not pdb_hierarchy.atoms().extract_i_seq().all_eq(0)
  bool_pdbh_selection = flex.bool(pdb_hierarchy.atoms_size(), False)
  if (selection is not None):
    if (isinstance(selection, flex.bool)):
      assert len(selection) == pdb_hierarchy.atoms_size()
      bool_pdbh_selection = selection
    elif (isinstance(selection, flex.size_t)):
      bool_pdbh_selection.set_selected(selection, True)
  if selection is None:
    bool_pdbh_selection = flex.bool(pdb_hierarchy.atoms_size(), True)
  actual_selection = bool_pdbh_selection.iselection()

  assert len(sites_cart) == len(actual_selection)

  if abs(sigma) < 1e-6:
    raise Sorry("Please set non-zero sigma for reference model restraints.")
  weight = 1.0 / (sigma**2)
  selection_to_sites_map = get_selection_to_sites_map(
                             sites_cart=sites_cart,
                             selection=actual_selection)
  residue_torsions = collect_residue_torsion_angles(
                   pdb_hierarchy=pdb_hierarchy,
                   atom_selection=bool_pdbh_selection,
                   chi_angles_only=chi_angles_only)
  for residue_info in residue_torsions:
    for chi in residue_info.chis:
      i_seqs = chi.i_seqs
      sites = []
      for i_seq in i_seqs:
        sites.append(selection_to_sites_map[i_seq])
      di = geometry_restraints.dihedral(
             sites=sites, angle_ideal=0.0, weight=weight)
      angle_ideal = di.angle_model
      dp = geometry_restraints.dihedral_proxy(
        i_seqs=i_seqs,
        angle_ideal=angle_ideal,
        weight=weight,
        limit=limit,
        top_out=top_out_potential,
        origin_id=origin_id)
      torsion_proxies.append(dp)
  return torsion_proxies

def get_selection_to_sites_map(sites_cart, selection):
  sites_selection_map = {}
  assert isinstance(selection, flex.size_t)
  for i, i_seq in enumerate(selection):
    sites_selection_map[i_seq] = sites_cart[i]
  return sites_selection_map

def add_coordinate_restraints(
      sites_cart,
      selection=None,
      sigma=0.5,
      limit=1.0,
      top_out_potential=False):
  import boost_adaptbx.boost.python as bp
  ext_rcp = bp.import_ext("mmtbx_reference_coordinate_ext")
  result = ext_rcp.shared_reference_coordinate_proxy()
  if (selection is not None):
    if (isinstance(selection, flex.bool)):
      selection = selection.iselection()
  if selection is None:
    selection = flex.bool(
      len(sites_cart),
      True).iselection()
  # Not clear why this assertion should present. What if we want to restrain
  # only part of the molecule?
  # Very likely the mechanics is similar to the above procedure (generate_torsion_restraints),
  # selection is related to sites_cart.
  assert len(sites_cart) == len(selection)
  weight = 1.0 / (sigma**2)
  for k, i_seq in enumerate(selection):
    i_seqs = [i_seq]
    ref_sites = sites_cart[k]
    proxy = ext_rcp.reference_coordinate_proxy(
              i_seqs=i_seqs,
              ref_sites=ref_sites,
              weight=weight,
              limit=limit,
              top_out=top_out_potential)
    result.append(proxy)
  return result

def exclude_outliers_from_reference_restraints_selection(
    pdb_hierarchy,
    restraints_selection):
  from mmtbx.validation.ramalyze import ramalyze
  # the import below is SLOW!!!
  from mmtbx.rotamer.rotamer_eval import RotamerEval
  assert restraints_selection is not None
  # ramachandran plot outliers
  rama_outlier_selection = ramalyze(pdb_hierarchy=pdb_hierarchy,
    outliers_only=False).outlier_selection()
  rama_outlier_selection = flex.bool(restraints_selection.size(),
    rama_outlier_selection)
  # rotamer outliers
  rota_outlier_selection = flex.size_t()
  rotamer_manager = RotamerEval() # SLOW!!!
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for residue_group in chain.residue_groups():
        conformers = residue_group.conformers()
        if(len(conformers)>1): continue
        for conformer in residue_group.conformers():
          residue = conformer.only_residue()
          if(rotamer_manager.evaluate_residue(residue)=="OUTLIER"):
            rota_outlier_selection.extend(residue.atoms().extract_i_seq())
  rota_outlier_selection = flex.bool(restraints_selection.size(),
    rota_outlier_selection)
  outlier_selection = rama_outlier_selection | rota_outlier_selection
  return restraints_selection & (~outlier_selection)
