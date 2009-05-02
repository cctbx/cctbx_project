from cctbx.array_family import flex
from cctbx import crystal
from cctbx import adp_restraints

class adp_similarity_restraints(object):
  def __init__(self, xray_structure, pair_sym_table=None, i_seqs=None,
               sigma=0.04, sigma_terminal=None, buffer_thickness=3.5):
    if sigma_terminal is None: sigma_terminal = 2 * sigma
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
    self.proxies = proxies

class rigid_bond_restraints(object):
  def __init__(self, xray_structure, pair_sym_table=None, i_seqs=None,
               sigma_12=0.01, sigma_13=None, buffer_thickness=3.5):
    """ sigma_12 and sigma_13 are the effective standard deviations used for
        1,2- and 1,3-distances respectively
    """
    if sigma_13 is None: sigma_13 = sigma_12
    proxies = adp_restraints.shared_rigid_bond_proxy()
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
      for j_seq in connectivity[i_seq]:
        if i_seqs is not None and j_seq not in i_seqs: continue
        if i_seq < j_seq:
          j_sym_ops = pair_sym_table[i_seq][j_seq]
        else:
          k_sym_ops = pair_sym_table[j_seq][i_seq]
        for sym_op in j_sym_ops:
          if sym_op.is_unit_mx() and i_seq < j_seq:
            weight = 1/(sigma_12*sigma_12)
            proxies.append(adp_restraints.rigid_bond_proxy(
              i_seqs=(i_seq,j_seq),weight=weight))
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
                    if sym_op.is_unit_mx() and i_seq < k_seq:
                      weight = 1/(sigma_13*sigma_13)
                      proxies.append(adp_restraints.rigid_bond_proxy(
                        i_seqs=(i_seq,k_seq),weight=weight))
    self.proxies = proxies

class isotropic_adp_restraints(object):
  def __init__(self, xray_structure, pair_sym_table=None, i_seqs=None,
               sigma=0.1, sigma_terminal=None, buffer_thickness=3.5):
    if sigma_terminal is None: sigma_terminal = 2 * sigma
    proxies = adp_restraints.shared_isotropic_adp_proxy()
    if pair_sym_table is None:
      asu_mappings = xray_structure.asu_mappings(buffer_thickness=buffer_thickness)
      pair_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
      scattering_types = xray_structure.scatterers().extract_scattering_types()
      pair_asu_table.add_covalent_pairs(
        scattering_types, exclude_scattering_types=flex.std_string(("H","D")))
      pair_sym_table = pair_asu_table.extract_pair_sym_table()
      connectivity = pair_sym_table.full_simple_connectivity()

    for i_seq , neighbours in enumerate(connectivity):
      if i_seqs is not None and i_seq not in i_seqs: continue
      elif scattering_types[i_seq] == 'H': continue
      if neighbours.size() <= 1:
        weight = 1/(sigma_terminal*sigma_terminal)
      else:
        weight = 1/(sigma*sigma)
      proxies.append(adp_restraints.isotropic_adp_proxy(
        i_seq=i_seq,weight=weight))
    self.proxies = proxies
