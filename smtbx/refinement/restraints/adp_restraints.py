from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from cctbx import crystal
from cctbx import adp_restraints

class adp_similarity_restraints(object):
  def __init__(self, xray_structure=None, pair_sym_table=None, proxies=None,
               i_seqs=None, sigma=0.04, sigma_terminal=None,
               buffer_thickness=3.5):
    assert [xray_structure, pair_sym_table].count(None) == 1
    if i_seqs is not None and len(i_seqs) == 0: i_seqs = None
    if sigma_terminal is None: sigma_terminal = 2 * sigma
    if proxies is None:
      proxies = adp_restraints.shared_adp_similarity_proxy()
    if pair_sym_table is None:
      asu_mappings = xray_structure.asu_mappings(buffer_thickness=buffer_thickness)
      pair_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
      scattering_types = xray_structure.scatterers().extract_scattering_types()
      pair_asu_table.add_covalent_pairs(
        scattering_types, exclude_scattering_types=flex.std_string(("H","D")))
      pair_sym_table = pair_asu_table.extract_pair_sym_table()
    connectivity = pair_sym_table.full_simple_connectivity()

    for i_seq, j_seq_dict in enumerate(pair_sym_table):
      if i_seqs is not None and i_seq not in i_seqs: continue
      for j_seq, sym_ops in j_seq_dict.items():
        if i_seqs is not None and j_seq not in i_seqs: continue
        for sym_op in sym_ops:
          if sym_op.is_unit_mx():
            i_is_terminal = (connectivity[i_seq].size() <= 1)
            j_is_terminal = (connectivity[j_seq].size() <= 1)
            if i_is_terminal or j_is_terminal:
              weight = 1/(sigma_terminal*sigma_terminal)
            else:
              weight = 1/(sigma*sigma)
            proxies.append(adp_restraints.adp_similarity_proxy(
              i_seqs=(i_seq,j_seq),weight=weight))
          break
    self.proxies = proxies

class rigid_bond_restraints(object):
  def __init__(self, xray_structure=None, pair_sym_table=None, proxies=None,
               i_seqs=None, sigma_12=0.01, sigma_13=None,
               buffer_thickness=3.5):
    """ sigma_12 and sigma_13 are the effective standard deviations used for
        1,2- and 1,3-distances respectively
    """
    assert [xray_structure, pair_sym_table].count(None) == 1
    if i_seqs is not None and len(i_seqs) == 0: i_seqs = None
    if sigma_13 is None: sigma_13 = sigma_12
    if proxies is None:
      proxies = adp_restraints.shared_rigid_bond_proxy()
    if pair_sym_table is None:
      asu_mappings = xray_structure.asu_mappings(buffer_thickness=buffer_thickness)
      pair_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
      scattering_types = xray_structure.scatterers().extract_scattering_types()
      pair_asu_table.add_covalent_pairs(
        scattering_types, exclude_scattering_types=flex.std_string(("H","D")))
      pair_sym_table = pair_asu_table.extract_pair_sym_table()
    connectivity = pair_sym_table.full_simple_connectivity()
    ij_seqs = set()

    for i_seq, j_seq_dict in enumerate(pair_sym_table):
      if i_seqs is not None and i_seq not in i_seqs: continue
      for j_seq in connectivity[i_seq]:
        if i_seqs is not None and j_seq not in i_seqs: continue
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
            proxies.append(adp_restraints.rigid_bond_proxy(
              i_seqs=(i_seq,j_seq),weight=weight))
            break
        if connectivity[j_seq].size() > 1:
          for k_seq in connectivity[j_seq]:
            if i_seqs is not None and k_seq not in i_seqs: continue
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
                      proxies.append(adp_restraints.rigid_bond_proxy(
                        i_seqs=(i_seq,k_seq),weight=weight))
                      break
                  break
    self.proxies = proxies

class rigu_restraints(object):
  def __init__(self, xray_structure=None, pair_sym_table=None, proxies=None,
               i_seqs=None, sigma_12=0.004, sigma_13=None,
               buffer_thickness=3.5):
    """ sigma_12 and sigma_13 are the effective standard deviations used for
        1,2- and 1,3-distances respectively
    """
    assert [xray_structure, pair_sym_table].count(None) == 1
    if i_seqs is not None and len(i_seqs) == 0: i_seqs = None
    if sigma_13 is None: sigma_13 = sigma_12
    if proxies is None:
      proxies = adp_restraints.shared_rigu_proxy()
    if pair_sym_table is None:
      asu_mappings = xray_structure.asu_mappings(buffer_thickness=buffer_thickness)
      pair_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
      scattering_types = xray_structure.scatterers().extract_scattering_types()
      pair_asu_table.add_covalent_pairs(
        scattering_types, exclude_scattering_types=flex.std_string(("H","D")))
      pair_sym_table = pair_asu_table.extract_pair_sym_table()
    connectivity = pair_sym_table.full_simple_connectivity()
    ij_seqs = set()

    for i_seq, j_seq_dict in enumerate(pair_sym_table):
      if i_seqs is not None and i_seq not in i_seqs: continue
      for j_seq in connectivity[i_seq]:
        if i_seqs is not None and j_seq not in i_seqs: continue
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
            proxies.append(adp_restraints.rigu_proxy(
              i_seqs=(i_seq,j_seq),weight=weight))
            break
        if connectivity[j_seq].size() > 1:
          for k_seq in connectivity[j_seq]:
            if i_seqs is not None and k_seq not in i_seqs: continue
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
                      proxies.append(adp_restraints.rigu_proxy(
                        i_seqs=(i_seq,k_seq),weight=weight))
                      break
                  break
    self.proxies = proxies

class isotropic_adp_restraints(object):
  def __init__(self, xray_structure, pair_sym_table=None, proxies=None,
               i_seqs=None, sigma=0.1, sigma_terminal=None, buffer_thickness=3.5):
    if sigma_terminal is None: sigma_terminal = 2 * sigma
    if i_seqs is not None and len(i_seqs) == 0: i_seqs = None
    if proxies is None:
      proxies = adp_restraints.shared_isotropic_adp_proxy()
    scattering_types = xray_structure.scatterers().extract_scattering_types()
    use_u_aniso = xray_structure.scatterers().extract_use_u_aniso()
    if pair_sym_table is None:
      asu_mappings = xray_structure.asu_mappings(buffer_thickness=buffer_thickness)
      pair_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
      pair_asu_table.add_covalent_pairs(
        scattering_types, exclude_scattering_types=flex.std_string(("H","D")))
      pair_sym_table = pair_asu_table.extract_pair_sym_table()
    connectivity = pair_sym_table.full_simple_connectivity()

    for i_seq, neighbours in enumerate(connectivity):
      if i_seqs is not None and i_seq not in i_seqs: continue
      elif scattering_types[i_seq] in ('H','D'): continue
      elif not use_u_aniso[i_seq]: continue
      if neighbours.size() <= 1:
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
      i_seqs = [i for i, s in enumerate(xray_structure.scatterers())
                if s.scattering_type not in ('H', 'D')]
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
      i_seqs = [i for i, s in enumerate(xray_structure.scatterers())
                if s.scattering_type not in ('H', 'D')]
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
      i_seqs = [i for i, s in enumerate(xray_structure.scatterers())
                if s.scattering_type not in ('H', 'D')]
    assert len(i_seqs) > 1
    proxies.append(adp_restraints.adp_volume_similarity_proxy(
      i_seqs=i_seqs, weight=weight))
    self.proxies = proxies
