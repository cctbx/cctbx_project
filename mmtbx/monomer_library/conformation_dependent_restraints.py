from cctbx import geometry_restraints
import libtbx.load_env

if (not libtbx.env.has_module("conformation_dependent_geometry")):
  conformation_dependent_geometry = None
  is_available = False
else:
  import conformation_dependent_geometry.angles
  is_available = True

class pgd_server(object):
  def __init__(self):
    self.lib = conformation_dependent_geometry.angles.create_all_databases(
      conformation_dependent_geometry.angles.databases)
  def lookup(self, residue, next_residue, phi, psi):
    fields, geometry = conformation_dependent_geometry.angles.get_geometry(
                         dblist=self.lib,
                         residue=residue,
                         next_residue=next_residue,
                         phi=phi,
                         psi=psi)
    return geometry

if (not is_available):
  pgd_lib = None
else:
  pgd_lib = pgd_server()

class conformation_dependent_restraints(object):
  def __init__(self, residue_name, next_residue_name, conformation_proxies, i_phi_proxy, i_psi_proxy, i_dynamic_angles, i_dynamic_dihedrals):
    self.residue_name = residue_name
    self.next_residue_name = next_residue_name
    self.conformation_proxies = conformation_proxies
    self.i_phi_proxy = i_phi_proxy
    self.i_psi_proxy = i_psi_proxy
    self.i_dynamic_angles = i_dynamic_angles
    self.i_dynamic_dihedrals = i_dynamic_dihedrals

  def _get_dihedral(sites_cart, dihedral_proxies, i_proxy):
    if i_proxy is not None:
      # compute dihedral
        dihedral = geometry_restraints.dihedral(
                     sites_cart=sites_cart,
                     proxy=dihedral_proxies[i_proxy])
        return dihedral
    return None

  _get_dihedral = staticmethod(_get_dihedral)

  def _get_angle(sites_cart, angle_proxies, i_proxy):
    if i_proxy is not None:
      # compute angle
      angle = geometry_restraints.angle(
                sites_cart=sites_cart,
                proxy=angle_proxies[i_proxy])
      return angle
    return None

  _get_angle = staticmethod(_get_angle)

  def _get_average_and_weight(geometry, param_name):
    pgd_param_average_name = \
      conformation_dependent_geometry.angles.get_database_attribute_average_name(
      geometry_name=param_name)
    pgd_param_deviation_name = \
      conformation_dependent_geometry.angles.get_database_attribute_deviation_name(
      geometry_name=param_name)

    return getattr(geometry, pgd_param_average_name), \
           1/((getattr(geometry, pgd_param_deviation_name))**2)

  _get_average_and_weight = staticmethod(_get_average_and_weight)

  def update_restraints(self, sites_cart, dihedral_proxies, angle_proxies):
    if (1):
      from libtbx.utils import null_out
      log = null_out()
    else:
      import sys
      log = sys.stdout
    try:
      phi = self._get_dihedral(
                    sites_cart=sites_cart,
                    dihedral_proxies=self.conformation_proxies,
                    i_proxy=self.i_phi_proxy)
    except Exception: # XXX BAD
      phi = None
    try:
      psi = self._get_dihedral(
                    sites_cart=sites_cart,
                    dihedral_proxies=self.conformation_proxies,
                    i_proxy=self.i_psi_proxy)
    except Exception: # XXX BAD
      psi = None

    # Shows real dihedral value
    if phi is not None:
      print >> log, 'phi', phi.angle_model, phi.delta
    if psi is not None:
      print >> log, 'psi', psi.angle_model, psi.delta

    if phi is not None and psi is not None:

      # get restraint from our database here
      geometry = pgd_lib.lookup(
                   residue=self.residue_name,
                   next_residue=self.next_residue_name,
                   phi=phi.angle_model,
                   psi=psi.angle_model)

      # grab angles from our database
      # using zip() on i_dynamic_angles and our values
      # plug it into restraints in angle proxies
      angles = [self._get_angle(
                  sites_cart=sites_cart,
                  angle_proxies=angle_proxies,
                  i_proxy=i_proxy)
                for i_proxy in self.i_dynamic_angles]

      for angle, angle_name, i_proxy in zip(
                                          angles,
                                          conformation_dependent_geometry.angles.angle_names,
                                          self.i_dynamic_angles):
        # i_dynamic_angles contains None for angles/atoms that don't
        # exist so don't have restraints to update.
        if i_proxy is not None:

          new_angle_ideal, new_weight = \
                           self._get_average_and_weight(geometry, angle_name)

          # Create a new angle proxy here with our restraint
          new_angle_proxy = geometry_restraints.angle_proxy(
                              i_seqs=angle_proxies[i_proxy].i_seqs,
                              angle_ideal=new_angle_ideal,
                              weight=new_weight
                              )

          # Overwrite the old proxy
          angle_proxies[i_proxy] = new_angle_proxy

          # Show that we actually did update the proxy
          proxy = angle_proxies[i_proxy]
          print >> log, self.residue_name, self.next_residue_name,
          print >> log, angle_name, proxy.i_seqs, proxy.angle_ideal, proxy.weight

      # grab dihedrals from our database
      # using zip() on i_dynamic_dihedrals and our values
      # plug it into restraints in dihedral proxies
      dihedrals = [self._get_dihedral(
                  sites_cart=sites_cart,
                  dihedral_proxies=dihedral_proxies,
                  i_proxy=i_proxy)
                for i_proxy in self.i_dynamic_dihedrals]

      for dihedral, dihedral_name, i_proxy in zip(
                                          dihedrals,
                                          conformation_dependent_geometry.angles.dihedral_names,
                                          self.i_dynamic_dihedrals):
        # i_dynamic_dihedrals contains None for dihedrals/atoms that don't
        # exist so don't have restraints to update.
        if i_proxy is not None:

          new_angle_ideal, new_weight = \
                           self._get_average_and_weight(geometry, dihedral_name)

          new_angle_ideal = dihedral_proxies[i_proxy].angle_ideal
          new_weight = dihedral_proxies[i_proxy].weight

          # Create a new dihedral proxy here with our restraint
          new_dihedral_proxy = geometry_restraints.dihedral_proxy(
                              i_seqs=dihedral_proxies[i_proxy].i_seqs,
                              angle_ideal=new_angle_ideal,
                              weight=new_weight
                              )

          # Overwrite the old proxy
          dihedral_proxies[i_proxy] = new_dihedral_proxy

          # Show that we actually did update the proxy
          proxy = dihedral_proxies[i_proxy]
          print >> log, self.residue_name, self.next_residue_name,
          print >> log, dihedral_name, proxy.i_seqs, proxy.angle_ideal, proxy.weight

def build_conformation_dependent_angle_proxies(
      angle_proxy_registry,
      dihedral_proxy_registry,
      monomer_mappings,
      connectivity_i_j,
      connectivity_j_k,
      sites_cart):
  """
  Sets up conformation_dependent_restraints object.

  Finds atom indexes from registries

  Looks for C_prev in m_i, most atoms in m_j, N_next in m_k.
  """
  assert len(monomer_mappings) == 3

  if monomer_mappings[1] is None:
    residue_name = None
  else:
    residue_name = monomer_mappings[1].residue_name

  if monomer_mappings[2] is None:
    next_residue_name = None
  else:
    next_residue_name = monomer_mappings[2].residue_name

  dihedral_i_proxies = []
  conformation_proxies = geometry_restraints.shared_dihedral_proxy()
  for dihedral_number, dihedral_definition in enumerate(conformation_dependent_geometry.angles.dihedral_atoms):
    i_seqs = []
    for residue_index, atom_name in dihedral_definition:
      mm=monomer_mappings[residue_index]
      if mm is not None:
        i_seqs.append(getattr(mm.expected_atoms.get(atom_name), "i_seq", None))
    # account for missing atoms
    if len(i_seqs) == 4 and i_seqs.count(None) == 0:
      if dihedral_number == 0:
        # phi/psi: Restraints don't matter, we throw them away. This is just
        # so we can get eventually get the current value.
        phi = geometry_restraints.dihedral_proxy(
                i_seqs=i_seqs,
                angle_ideal=60.0,
                weight=1/20.0**2,
                periodicity=3
                )
        conformation_proxies.append(phi)
        dihedral_i_proxies.append(None)
        i_phi_proxy = 0
      elif dihedral_number == 1:
        psi = geometry_restraints.dihedral_proxy(
                i_seqs=i_seqs,
                angle_ideal=160.0,
                weight=1/30.0**2,
                periodicity=3
                )
        conformation_proxies.append(psi)
        dihedral_i_proxies.append(None)
        i_psi_proxy = 1
      elif dihedral_number == 2:
        # omega
        dihedral_i_proxy, dihedral_sign = dihedral_proxy_registry.lookup_i_proxy(i_seqs)
        dihedral_i_proxies.append(dihedral_i_proxy)
      else:
        pass
    # we're on a dihedral with missing length or atoms that are None
    else:
      if dihedral_number == 0:
        i_phi_proxy = None
      elif dihedral_number == 1:
        i_psi_proxy = None
      elif dihedral_number == 2:
        dihedral_i_proxies.append(None)
      else:
        pass

  angle_i_proxies = []
  for angle_definition in conformation_dependent_geometry.angles.angle_atoms:
    i_seqs = []
    for residue_index, atom_name in angle_definition:
      mm=monomer_mappings[residue_index]
      if mm is not None:
        i_seqs.append(getattr(mm.expected_atoms.get(atom_name), "i_seq", None))
    # account for missing atoms

    if len(i_seqs) == 3 and i_seqs.count(None) == 0:
      # go into angle_proxy_registry
      angle_i_proxy = angle_proxy_registry.lookup_i_proxy(i_seqs)
      angle_i_proxies.append(angle_i_proxy)
    else:
      # By filling in None for blanks, we can assume the lengths are equal
      # here and other places angle_atoms/angle_names are used
      angle_i_proxies.append(None)

  cfd = conformation_dependent_restraints(
          residue_name=residue_name,
          next_residue_name=next_residue_name,
          conformation_proxies=conformation_proxies,
          i_phi_proxy=i_phi_proxy,
          i_psi_proxy=i_psi_proxy,
          i_dynamic_angles=angle_i_proxies,
          i_dynamic_dihedrals=dihedral_i_proxies
          )

  return cfd
