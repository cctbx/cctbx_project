from __future__ import absolute_import, division, print_function
from cctbx import crystal
from cctbx import adp_restraints
from cctbx.array_family import flex


def covalent_pair_sym_table(xray_structure, buffer_thickness=3.5,
                            exclude_hydrogens=False):
  """ Which atoms are bonded to which.

  Hydrogen and deuterium are kept unless asked for otherwise. They have to be:
  which pairs of atoms get a restraint is read off this table, and dropping
  them would leave every hydrogen unrestrained -- which matters wherever
  hydrogen ADPs are refined rather than riding, as they are under an aspherical
  model.
  """
  asu_mappings = xray_structure.asu_mappings(buffer_thickness=buffer_thickness)
  pair_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
  scattering_types = xray_structure.scatterers().extract_scattering_types()
  if exclude_hydrogens:
    pair_asu_table.add_covalent_pairs(
      scattering_types,
      exclude_scattering_types=flex.std_string(("H", "D")))
  else:
    pair_asu_table.add_covalent_pairs(scattering_types)
  return pair_asu_table.extract_pair_sym_table()


def terminal_connectivity(xray_structure, buffer_thickness=3.5):
  """ Connectivity for deciding whether an atom is terminal, hydrogens aside.

  Being terminal governs the weight a restraint gets -- a terminal atom's
  displacement is the least well determined, so it is restrained more loosely.
  What decides that is the heavy-atom skeleton: a methyl carbon carries one
  carbon and three hydrogens, and it is the archetypal loosely-held group, yet
  counting its hydrogens makes it look four-coordinate and gives it the tight
  weight of a backbone atom. Likewise a hydroxyl oxygen is terminal whatever
  its hydrogen does.

  Separate from the table above on purpose: which atoms are restrained is one
  question and how tightly is another, and only the second one wants the
  hydrogens gone.

  Returns None when there is no structure to derive it from, and the caller
  then falls back on the connectivity it was given.
  """
  if xray_structure is None:
    return None
  return covalent_pair_sym_table(
    xray_structure, buffer_thickness,
    exclude_hydrogens=True).full_simple_connectivity()


class adp_similarity_restraints(object):
  def __init__(self, xray_structure=None, pair_sym_table=None, proxies=None,
               i_seqs=None, sigma=0.04, sigma_terminal=None,
               buffer_thickness=3.5, connectivity=None):
    assert [xray_structure, pair_sym_table].count(None) == 1
    scatterers = None
    if xray_structure is not None:
      scatterers = xray_structure.scatterers()

    def is_suitable(idx):
      """ Whether this scatterer's ADPs are being refined, so worth restraining.

      With no structure there are no flags to consult -- this class may be given
      a pair_sym_table and nothing else, which the assertion above allows for --
      and the answer is then that nothing is known against the scatterer, not
      that it is unsuitable. Answering False there rejects every scatterer and
      the caller gets no restraints at all.
      """
      if scatterers is None:
        return True
      if scatterers[idx].flags.use_u_iso() and scatterers[idx].flags.grad_u_iso():
        return True
      if scatterers[idx].flags.use_u_aniso() and scatterers[idx].flags.grad_u_aniso():
        return True
      return False

    if i_seqs is not None and len(i_seqs) == 0: i_seqs = None
    if sigma_terminal is None: sigma_terminal = 2 * sigma
    if proxies is None:
      proxies = adp_restraints.shared_adp_similarity_proxy()
    if pair_sym_table is None:
      pair_sym_table = covalent_pair_sym_table(xray_structure, buffer_thickness)
    if connectivity is None:
      connectivity = pair_sym_table.full_simple_connectivity()
    # hydrogens do not make an atom non-terminal; see terminal_connectivity
    terminal = terminal_connectivity(xray_structure, buffer_thickness)
    if terminal is None:
      terminal = connectivity

    for i_seq, j_seq_dict in enumerate(pair_sym_table):
      if i_seqs is not None and i_seq not in i_seqs: continue
      if not is_suitable(i_seq): continue
      for j_seq, sym_ops in j_seq_dict.items():
        if i_seqs is not None and j_seq not in i_seqs: continue
        if not is_suitable(j_seq): continue
        for sym_op in sym_ops:
          if sym_op.is_unit_mx():
            i_is_terminal = (terminal[i_seq].size() <= 1)
            j_is_terminal = (terminal[j_seq].size() <= 1)
            if i_is_terminal or j_is_terminal:
              weight = 1/(sigma_terminal*sigma_terminal)
            else:
              weight = 1/(sigma*sigma)
            proxies.append(adp_restraints.adp_similarity_proxy(
              i_seqs=(i_seq,j_seq),weight=weight))
          break
    self.proxies = proxies


# for use bu RIGU and DELU below
def build_proxies(proxies, proxy_type, sigma_12, sigma_13,
    xray_structure=None, pair_sym_table=None, i_seqs=None,
    buffer_thickness=3.5, connectivity=None):
  scatterers = None
  if xray_structure is not None:
    scatterers = xray_structure.scatterers()

  def is_suitable(idx):
    """ \\copydoc adp_similarity_restraints.is_suitable

    RIGU and DELU are restraints between anisotropic atoms, so an isotropic one
    is filtered out here whatever its flags say.
    """
    if scatterers is None:
      return True
    return scatterers[idx].flags.use_u_aniso() \
       and scatterers[idx].flags.grad_u_aniso()

  if pair_sym_table is None:
    pair_sym_table = covalent_pair_sym_table(xray_structure, buffer_thickness)
  if connectivity is None:
    connectivity = pair_sym_table.full_simple_connectivity()
  ij_seqs = set()
  for i_seq, j_seq_dict in enumerate(pair_sym_table):
    if i_seqs is not None and i_seq not in i_seqs: continue
    if not is_suitable(i_seq): continue
    for j_seq in connectivity[i_seq]:
      if i_seqs is not None and j_seq not in i_seqs: continue
      if not is_suitable(j_seq): continue
      if i_seq < j_seq:
        j_sym_ops = pair_sym_table[i_seq][j_seq]
      else:
        k_sym_ops = pair_sym_table[j_seq][i_seq]
      for sym_op in j_sym_ops:
        if (    sym_op.is_unit_mx()
            and i_seq < j_seq
            and (i_seq, j_seq) not in ij_seqs):
          ij_seqs.add((i_seq, j_seq))
          weight = 1/(sigma_12*sigma_12)
          proxies.append(proxy_type(
            i_seqs=(i_seq,j_seq),weight=weight))
          break
      if connectivity[j_seq].size() > 1:
        for k_seq in connectivity[j_seq]:
          if i_seqs is not None and k_seq not in i_seqs: continue
          if not is_suitable(k_seq): continue
          if k_seq != i_seq:
            for sym_op in j_sym_ops:
              if sym_op.is_unit_mx():
                if j_seq < k_seq:
                  k_sym_ops = pair_sym_table[j_seq][k_seq]
                else:
                  k_sym_ops = pair_sym_table[k_seq][j_seq]
                for sym_op in k_sym_ops:
                  if (    sym_op.is_unit_mx()
                      and i_seq < k_seq
                      and (i_seq, k_seq) not in ij_seqs):
                    ij_seqs.add((i_seq, k_seq))
                    weight = 1/(sigma_13*sigma_13)
                    proxies.append(proxy_type(
                      i_seqs=(i_seq,k_seq),weight=weight))
                    break
                break

class rigid_bond_restraints(object):
  def __init__(self, xray_structure=None, pair_sym_table=None, proxies=None,
               i_seqs=None, sigma_12=0.01, sigma_13=None,
               buffer_thickness=3.5, connectivity=None):
    """ sigma_12 and sigma_13 are the effective standard deviations used for
        1,2- and 1,3-distances respectively
    """
    assert [xray_structure, pair_sym_table].count(None) == 1
    if i_seqs is not None and len(i_seqs) == 0: i_seqs = None
    if sigma_13 is None: sigma_13 = sigma_12
    if proxies is None:
      proxies = adp_restraints.shared_rigid_bond_proxy()
    build_proxies(proxies, adp_restraints.rigid_bond_proxy, sigma_12, sigma_13,
      xray_structure=xray_structure, pair_sym_table=pair_sym_table,
      i_seqs=i_seqs, buffer_thickness=buffer_thickness, connectivity=connectivity)
    self.proxies = proxies

class rigu_restraints(object):
  def __init__(self, xray_structure=None, pair_sym_table=None, proxies=None,
               i_seqs=None, sigma_12=0.004, sigma_13=None,
               buffer_thickness=3.5, connectivity=None):
    """ sigma_12 and sigma_13 are the effective standard deviations used for
        1,2- and 1,3-distances respectively
    """
    assert [xray_structure, pair_sym_table].count(None) == 1
    if i_seqs is not None and len(i_seqs) == 0: i_seqs = None
    if sigma_13 is None: sigma_13 = sigma_12
    if proxies is None:
      proxies = adp_restraints.shared_rigu_proxy()

    build_proxies(proxies, adp_restraints.rigu_proxy, sigma_12, sigma_13,
      xray_structure=xray_structure, pair_sym_table=pair_sym_table,
      i_seqs=i_seqs, buffer_thickness=buffer_thickness, connectivity=connectivity)

    self.proxies = proxies

class isotropic_adp_restraints(object):
  def __init__(self, xray_structure, pair_sym_table=None, proxies=None,
               i_seqs=None, sigma=0.1, sigma_terminal=None,
                buffer_thickness=3.5, connectivity=None):
    if sigma_terminal is None: sigma_terminal = 2 * sigma
    if i_seqs is not None and len(i_seqs) == 0: i_seqs = None
    if proxies is None:
      proxies = adp_restraints.shared_isotropic_adp_proxy()
    scattering_types = xray_structure.scatterers().extract_scattering_types()
    use_u_aniso = xray_structure.scatterers().extract_use_u_aniso()
    if pair_sym_table is None:
      pair_sym_table = covalent_pair_sym_table(xray_structure, buffer_thickness)
    if connectivity is None:
      connectivity = pair_sym_table.full_simple_connectivity()
    # hydrogens do not make an atom non-terminal; see terminal_connectivity
    terminal = terminal_connectivity(xray_structure, buffer_thickness)
    if terminal is None:
      terminal = connectivity

    for i_seq, neighbours in enumerate(connectivity):
      if i_seqs is not None and i_seq not in i_seqs: continue
      elif not use_u_aniso[i_seq]: continue
      if terminal[i_seq].size() <= 1:
        weight = 1/(sigma_terminal*sigma_terminal)
      else:
        weight = 1/(sigma*sigma)
      proxies.append(adp_restraints.isotropic_adp_proxy(
        i_seqs=(i_seq,),weight=weight))
    self.proxies = proxies

class fixed_u_eq_adp_restraints(object):
  def __init__(self, xray_structure, u_eq_ideal, proxies=None,
               i_seqs=None, sigma=0.1):
    if proxies is None:
      proxies = adp_restraints.shared_fixed_u_eq_adp_proxy()
    weight = 1/(sigma*sigma)
    if i_seqs is None:
      i_seqs = [i for i, s in enumerate(xray_structure.scatterers())]
    for i_seq in i_seqs:
      proxies.append(adp_restraints.fixed_u_eq_adp_proxy(
        i_seqs=(i_seq,),weight=weight, u_eq_ideal=u_eq_ideal))
    self.proxies = proxies

class adp_u_eq_similarity_restraints(object):
  def __init__(self, xray_structure, proxies=None,
               i_seqs=None, sigma=0.1):
    if proxies is None:
      proxies = adp_restraints.shared_adp_u_eq_similarity_proxy()
    weight = 1/(sigma*sigma)
    if i_seqs is None:
      i_seqs = [i for i, s in enumerate(xray_structure.scatterers())]
    assert len(i_seqs) > 1
    proxies.append(adp_restraints.adp_u_eq_similarity_proxy(
      i_seqs=i_seqs, weight=weight))
    self.proxies = proxies

class adp_volume_similarity_restraints(object):
  def __init__(self, xray_structure, proxies=None,
               i_seqs=None, sigma=0.1):
    if proxies is None:
      proxies = adp_restraints.shared_adp_volume_similarity_proxy()
    weight = 1/(sigma*sigma)
    if i_seqs is None:
      i_seqs = [i for i, s in enumerate(xray_structure.scatterers())]
    assert len(i_seqs) > 1
    proxies.append(adp_restraints.adp_volume_similarity_proxy(
      i_seqs=i_seqs, weight=weight))
    self.proxies = proxies
