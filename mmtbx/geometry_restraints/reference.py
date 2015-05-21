
from __future__ import division
from cctbx.array_family import flex
from cctbx import geometry_restraints
from cctbx.geometry_restraints import shared_dihedral_proxy
import iotbx.pdb
import boost.python
ext = boost.python.import_ext("mmtbx_reference_coordinate_ext")


def add_coordinate_restraints(
      sites_cart,
      selection=None,
      sigma=0.5,
      limit=1.0,
      top_out_potential=False):
  import boost.python
  ext_rcp = boost.python.import_ext("mmtbx_reference_coordinate_ext")
  result = ext_rcp.shared_reference_coordinate_proxy()
  if (selection is not None):
    if (isinstance(selection, flex.bool)):
      selection = selection.iselection()
  if selection is None:
    selection = flex.bool(
      len(sites_cart),
      True).iselection()
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

class manager(object):
  def __init__(self):
    self.reference_torsion_proxies = None

  def add_torsion_restraints(
        self,
        pdb_hierarchy,
        sites_cart,
        selection=None,
        sigma=2.5,
        limit=15.0,
        chi_angles_only=False,
        top_out_potential=False):
    assert [atom.i_seq for atom in pdb_hierarchy.atoms()].count(0) <= 1
    from mmtbx.torsion_restraints.reference_model import build_torsion_proxies
    local_reference_torsion_proxies = \
        build_torsion_proxies(
          pdb_hierarchy=pdb_hierarchy,
          sites_cart=sites_cart,
          selection=selection,
          sigma=sigma,
          limit=limit,
          chi_angles_only=chi_angles_only,
          top_out_potential=top_out_potential)
    if self.reference_torsion_proxies is None:
      self.reference_torsion_proxies = local_reference_torsion_proxies
    else:
      self.add_existing_reference_torsion_proxies(local_reference_torsion_proxies)

  def add_existing_reference_torsion_proxies(self, proxies):
    if self.reference_torsion_proxies is None:
      self.reference_torsion_proxies = shared_dihedral_proxy()
    if proxies is not None:
      for dp in proxies:
        self.reference_torsion_proxies.append(dp)

  def remove_torsion_restraints(self, selection):
    self.reference_torsion_proxies = \
      self.reference_torsion_proxies.proxy_remove(selection=selection)

  def remove_chi_angle_restraints(self,
                                  pdb_hierarchy,
                                  selection=None):
    if self.reference_torsion_proxies is not None:
      awl = tuple(pdb_hierarchy.atoms_with_labels())
      updated_torsion_restraints = \
        geometry_restraints.shared_dihedral_proxy()
      for dp in self.reference_torsion_proxies:
        if not torsion_is_chi_angle(dp, awl):
          updated_torsion_restraints.append(dp)
      if len(updated_torsion_restraints) > 0:
        self.reference_torsion_proxies = updated_torsion_restraints
      else:
        self.reference_torsion_proxies = None

  def target_and_gradients(self,
                           sites_cart,
                           gradient_array,
                           unit_cell=None):
    target = 0.0
    if self.reference_torsion_proxies is not None:
      if unit_cell is None:
        target += geometry_restraints.dihedral_residual_sum(
                    sites_cart=sites_cart,
                    proxies=self.reference_torsion_proxies,
                    gradient_array=gradient_array)
      else:
        target += geometry_restraints.dihedral_residual_sum(
                    unit_cell=unit_cell,
                    sites_cart=sites_cart,
                    proxies=self.reference_torsion_proxies,
                    gradient_array=gradient_array)
    return target

  def select(self, n_seq, iselection):
    if self.reference_torsion_proxies is not None:
      self.reference_torsion_proxies = self.reference_torsion_proxies.\
          proxy_select(n_seq, iselection)

def torsion_is_chi_angle(dp, atoms_with_labels):
  get_class = iotbx.pdb.common_residue_names_get_class
  resname = None
  atoms = []
  for i_seq in dp.i_seqs:
    if resname is None:
      resname = atoms_with_labels[i_seq].resname
    else:
      if resname != atoms_with_labels[i_seq].resname:
        return False
    atoms.append(atoms_with_labels[i_seq])
  if (get_class(resname) != "common_amino_acid"):
    return False
  resseq = None
  for atom in atoms:
    if resseq is None:
      resseq = atom.resseq
    else:
      if resseq != atom.resseq:
        return False
  return True
