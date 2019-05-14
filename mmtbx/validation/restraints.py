
"""
Validation of models of any type against basic covalent geometry restraints.
By default this will flag all restrained atoms deviating by more than 4 sigma
from the target value.
"""

from __future__ import absolute_import, division, print_function
from mmtbx.validation import atoms, validation, get_atoms_info
from libtbx.str_utils import make_sub_header
from libtbx import slots_getstate_setstate
from math import sqrt
import sys

__restraint_attr__ = [
  "sigma",
  "target",
  "model",
  "delta",
  "residual",
] # XXX others?

class restraint(atoms):
  n_atoms = None
  """
  Base class for covalent sterochemistry restraint outliers (except for
  planarity, which is weird and different).  Unlike most of the other
  outlier implementations elsewhere in the validation module, the restraint
  outliers are printed on multiple lines to facilitate display of the atoms
  involved.
  """
  __slots__ = atoms.__slots__ + __restraint_attr__
  def __init__(self, **kwds):
    atoms.__init__(self, **kwds)
    if (self.n_atoms is not None):
      assert (len(self.atoms_info) == self.n_atoms)
    if (self.score is None):
      self.score = abs(self.delta / self.sigma)

  @staticmethod
  def header():
    return "%-20s  %7s  %7s  %7s  %6s  %6s  %10s" % ("atoms", "ideal", "model",
      "delta", "sigma", "residual", "deviation")

  def as_table_row_phenix(self):
    """
    Values for populating ListCtrl in Phenix GUI.
    """
    atoms_str = ", ".join([ a.id_str() for a in self.atoms_info ])
    return [ atoms_str, self.target, self.model, self.score ]

  def id_str(self, ignore_altloc=None):
    return ",".join([ a.id_str() for a in self.atoms_info ])

  def as_string(self, prefix=""):
    id_strs = [ a.id_str() for a in self.atoms_info ]
    id_len = max([ len(s) for s in id_strs ])
    lines = []
    for atom_str in id_strs :
      lines.append("%s%-20s" % (prefix, atom_str))
    lines[-1] += "  " + self.format_values()
    return "\n".join(lines)

  def format_values(self):
    return "%7.2f  %7.2f  %7.2f  %6.2e  %6.2e  %4.1f*sigma" % (self.target,
      self.model, self.delta, self.sigma, self.residual, self.score)

  def __cmp__(self, other):
    return cmp(other.score, self.score)

  def kinemage_key(self):
    atom0 = self.atoms_info[0]
    atom_names = [ a.name.strip().lower() for a in self.atoms_info ]
    kin_key = "%1s%3s%2s%4s%1s %s" % (self.get_altloc(),
      atom0.resname.lower(), atom0.chain_id, atom0.resseq, atom0.icode,
      "-".join(atom_names))
    return kin_key

class bond(restraint):
  n_atoms = 2
  __bond_attr__ = [
    "slack",
    "symop",
  ]
  __slots__ = restraint.__slots__ + __bond_attr__
  def as_table_row_phenix(self):
    return [ self.atoms_info[0].id_str(), self.atoms_info[1].id_str(),
             self.target, self.model, self.score ]

  @staticmethod
  def header():
    return "%-20s  %5s  %6s  %6s  %6s  %6s  %8s  %10s" % ("atoms", "ideal",
      "model", "delta", "sigma", "slack", "residual", "deviation")

  def formate_values(self):
    return "%5.3f  %6.2f  %6.3f  %6.3f  %6.2e  %8.2e  %4.1f*sigma" % \
      (self.target, self.model, self.delta, self.sigma, self.slack,
       self.residual, abs(self.score))

  def as_kinemage(self):
    from mmtbx.kinemage.validation import bond_outlier_as_kinemage
    return bond_outlier_as_kinemage(self)

class angle(restraint):
  n_atoms = 3
  def as_kinemage(self):
    from mmtbx.kinemage.validation import angle_outlier_as_kinemage
    return angle_outlier_as_kinemage(self)

class dihedral(restraint):
  n_atoms = 4
  def as_kinemage(self):
    return None

class chirality(restraint):
  def as_kinemage(self):
    return None

class planarity(restraint):
  __slots__ = atoms.__slots__ + [
    "rms_deltas",
    "delta_max",
    "residual",
  ]
  def as_table_row_phenix(self):
    atoms_str = ", ".join([ a.id_str() for a in self.atoms_info ])
    return [ atoms_str, self.delta_max, self.rms_deltas, self.score ]

  @staticmethod
  def header():
    return "%-20s  %10s  %10s  %10s  %10s" % ("atoms", "rms_deltas",
      "delta_max", "residual", "deviation")

  def format_values(self):
    return "%10.3f  %10.3f  %10.2f  %4.1f*sigma" % (self.rms_deltas,
      self.delta_max, self.residual, self.score)

  def as_kinemage(self):
    return None

class restraint_validation(validation):
  """
  Base class for collecting information about all restraints of a certain
  type, including overall statistics and individual outliers.
  """
  restraint_type = None
  kinemage_header = None
  gui_list_headers = ["Atoms","Ideal value","Model value","Deviation (sigmas)"]
  gui_formats = ["%s", "%.3f", "%.3f", "%.1f"]
  wx_column_widths = [500, 100, 100, 180]
  __restraints_attr__ = [
    "min",
    "max",
    "mean",
    "z_min",
    "z_max",
    "z_mean",
    "target",
  ]
  __slots__ = validation.__slots__ + __restraints_attr__
  def __init__(self,
      pdb_atoms,
      sites_cart,
      energies_sites,
      restraint_proxies,
      unit_cell,
      ignore_hd=True,
      sigma_cutoff=4.0,
      outliers_only=True,
      use_segids_in_place_of_chainids=False):
    validation.__init__(self)
    self.z_min = self.z_max = self.z_mean = None
    deviations_method = getattr(energies_sites, "%s_deviations" %
      self.restraint_type)
    self.min, self.max, self.mean = deviations_method()
    target = getattr(energies_sites, "%s_residual_sum" %
      self.restraint_type)
    self.n_total = getattr(energies_sites, "n_%s_proxies" %
      self.restraint_type)
    if (self.n_total > 0):
      self.target = target / self.n_total
    else :
      self.target = 0
    deviations_z_method = getattr(energies_sites, "%s_deviations_z" %
      self.restraint_type, None)
    if (deviations_z_method is not None):
      deviations_z = deviations_z_method()
      self.z_min, self.z_max, self.z_mean = deviations_z_method()
    self.results = sorted(self.get_outliers(
      proxies=restraint_proxies,
      unit_cell=unit_cell,
      sites_cart=sites_cart,
      pdb_atoms=pdb_atoms,
      sigma_cutoff=sigma_cutoff,
      outliers_only=outliers_only,
      use_segids_in_place_of_chainids=use_segids_in_place_of_chainids))
    self.n_outliers = len(self.results)

  def get_outliers(self, proxies, unit_cell, sites_cart, pdb_atoms,
      sigma_cutoff):
    raise NotImplementedError()

  def show_old_output(self, *args, **kwds):
    raise NotImplementedError()

  def show(self, out=sys.stdout, prefix="  ", verbose=True):
    if (len(self.results) > 0):
      print(prefix + self.get_result_class().header(), file=out)
      for result in self.results :
        print(result.as_string(prefix=prefix), file=out)
    self.show_summary(out=out, prefix=prefix)

  def show_summary(self, out=sys.stdout, prefix=""):
    if (self.n_total == 0):
      print(prefix + "No restraints of this type.", file=out)
      return
    elif (self.n_outliers == 0):
      print(prefix + \
        "All restrained atoms within 4.0 sigma of ideal values.", file=out)
    print("", file=out)
    if (self.z_mean is not None):
      print(prefix + "Min. delta:  %7.3f (Z=%7.3f)" % (self.min,
        self.z_min), file=out)
      print(prefix + "Max. delta:  %7.3f (Z=%7.3f)" % (self.max,
        self.z_max), file=out)
      print(prefix + "Mean delta:  %7.3f (Z=%7.3f)" % (self.mean,
        self.z_mean), file=out)
    else :
      print(prefix + "Min. delta:  %7.3f" % self.min, file=out)
      print(prefix + "Max. delta:  %7.3f" % self.max, file=out)
      print(prefix + "Mean delta:  %7.3f" % self.mean, file=out)

  def as_kinemage(self, chain_id=None):
    header = self.kinemage_header
    if (header is not None):
      kin_blocks = []
      for result in self.results :
        if (result.is_outlier()) and (result.is_in_chain(chain_id)):
          outlier_kin_txt = result.as_kinemage()
          if (outlier_kin_txt is not None):
            kin_blocks.append(outlier_kin_txt)
      return header + "\n".join(kin_blocks)
    return None

class bonds(restraint_validation):
  restraint_type = "bond"
  restraint_label = "Bond length"
  kinemage_header = "@subgroup {length devs} dominant\n"
  gui_list_headers = ["Atom 1","Atom 2","Ideal value","Model value",
                      "Deviation (sigmas)"]
  gui_formats = ["%s", "%s", "%.3f", "%.3f", "%.1f"]
  wx_column_widths = [150, 150, 100, 100, 180]

  def get_result_class(self) : return bond

  def get_outliers(self, proxies, unit_cell, sites_cart, pdb_atoms,
      sigma_cutoff, outliers_only=True,
      use_segids_in_place_of_chainids=False):
    from scitbx.array_family import flex
    from cctbx.geometry_restraints.linking_class import linking_class
    origin_ids = linking_class()
    site_labels = flex.bool(sites_cart.size(), True).iselection()
    sorted_table, not_shown = proxies.get_sorted(
      by_value="residual",
      sites_cart=sites_cart,
      site_labels=site_labels,
      origin_id=origin_ids.get_origin_id('covalent geometry'))
    # this can happen for C-alpha-only models, etc.
    if (sorted_table is None):
      return []
    outliers = []
    for restraint_info in sorted_table :
      (i_seq, j_seq, i_seqs, ideal, model, slack, delta, sigma, weight, residual, sym_op_j,
       rt_mx) = restraint_info
      bond_atoms = get_atoms_info(pdb_atoms, iselection=i_seqs,
        use_segids_in_place_of_chainids=use_segids_in_place_of_chainids)
      if sym_op_j:
        import scitbx
        m3 = rt_mx.r().as_double()
        m3 = scitbx.matrix.sqr(m3)
        t = rt_mx.t().as_double()
        t = scitbx.matrix.col((t[0],t[1],t[2]))
        xyz = unit_cell.fractionalize(flex.vec3_double([bond_atoms[1].xyz]))
        new_xyz = unit_cell.orthogonalize(m3.elems*xyz+t)
        bond_atoms[1].xyz = new_xyz[0]
      outlier = bond(
        atoms_info=bond_atoms,
        target=ideal,
        model=model,
        sigma=sigma,
        slack=slack,
        delta=delta,
        residual=residual,
        symop=sym_op_j,
        outlier=True,
        xyz=get_mean_xyz(bond_atoms))
      if (outlier.score > sigma_cutoff):
        outliers.append(outlier)
      elif (not outliers_only):
        outlier.outlier=False
        outliers.append(outlier)
    return outliers

class angles(restraint_validation):
  restraint_type = "angle"
  restraint_label = "Bond angle"
  kinemage_header = "@subgroup {geom devs} dominant\n"
  def get_result_class(self) : return angle

  def get_outliers(self, proxies, unit_cell, sites_cart, pdb_atoms,
      sigma_cutoff, outliers_only=True,
      use_segids_in_place_of_chainids=False):
    import cctbx.geometry_restraints
    sorted = _get_sorted(proxies,
      unit_cell=unit_cell,
      sites_cart=sites_cart,
      pdb_atoms=pdb_atoms,
      use_segids_in_place_of_chainids=use_segids_in_place_of_chainids)
    outliers = []
    for proxy, proxy_atoms in sorted :
      restraint = cctbx.geometry_restraints.angle(
        unit_cell=unit_cell,
        proxy=proxy,
        sites_cart=sites_cart)
      outlier = angle(
        atoms_info=proxy_atoms,
        target=restraint.angle_ideal,
        delta=restraint.delta,
        model=restraint.angle_model,
        sigma=cctbx.geometry_restraints.weight_as_sigma(restraint.weight),
        residual=restraint.residual(),
        outlier=True,
        xyz=proxy_atoms[1].xyz)
      if (outlier.score > sigma_cutoff):
        outliers.append(outlier)
      elif (not outliers_only):
        outlier.outlier=False
        outliers.append(outlier)
    return outliers

class dihedrals(restraint_validation):
  restraint_type = "dihedral"
  restraint_label = "Dihedral angle"
  def get_result_class(self) : return dihedral

  def get_outliers(self, proxies, unit_cell, sites_cart, pdb_atoms,
      sigma_cutoff, outliers_only=True,
      use_segids_in_place_of_chainids=False):
    import cctbx.geometry_restraints
    sorted = _get_sorted(proxies,
      unit_cell=unit_cell,
      sites_cart=sites_cart,
      pdb_atoms=pdb_atoms)
    outliers = []
    for proxy, proxy_atoms in sorted :
      restraint = cctbx.geometry_restraints.dihedral(
        unit_cell=unit_cell,
        proxy=proxy,
        sites_cart=sites_cart)
      outlier = dihedral(
        atoms_info=proxy_atoms,
        target=restraint.angle_ideal,
        delta=restraint.delta,
        model=restraint.angle_model,
        sigma=cctbx.geometry_restraints.weight_as_sigma(restraint.weight),
        residual=restraint.residual(),
        xyz=get_mean_xyz([proxy_atoms[1], proxy_atoms[2]]),
        outlier=True)
      if (outlier.score > sigma_cutoff):
        outliers.append(outlier)
      elif (not outliers_only):
        outlier.outlier=False
        outliers.append(outlier)
    return outliers

class chiralities(restraint_validation):
  restraint_type = "chirality"
  restraint_label = "Chiral volume"
  def get_result_class(self) : return chirality

  def get_outliers(self, proxies, unit_cell, sites_cart, pdb_atoms,
      sigma_cutoff, outliers_only=True,
      use_segids_in_place_of_chainids=False):
    import cctbx.geometry_restraints
    sorted = _get_sorted(proxies,
      unit_cell=None,
      sites_cart=sites_cart,
      pdb_atoms=pdb_atoms)
    outliers = []
    for proxy, proxy_atoms in sorted :
      restraint = cctbx.geometry_restraints.chirality(
        proxy=proxy,
        sites_cart=sites_cart)
      outlier = chirality(
        atoms_info=proxy_atoms,
        target=restraint.volume_ideal,
        delta=restraint.delta,
        model=restraint.volume_model,
        sigma=cctbx.geometry_restraints.weight_as_sigma(restraint.weight),
        residual=restraint.residual(),
        outlier=True,
        xyz=get_mean_xyz(proxy_atoms))
      if (outlier.score > sigma_cutoff):
        outliers.append(outlier)
      elif (not outliers_only):
        outlier.outlier=False
        outliers.append(outlier)
    return outliers

class planarities(restraint_validation):
  restraint_type = "planarity"
  restraint_label = "Planar group"
  gui_list_headers = ["Atoms", "Max. delta", "RMS(delta)", "Deviation (sigmas)"]
  gui_formats = ["%s", "%.3f", "%.3f", "%.1f"]
  wx_column_widths = [250, 100, 100, 130]

  def get_result_class(self) : return planarity

  def get_outliers(self, proxies, unit_cell, sites_cart, pdb_atoms,
      sigma_cutoff, outliers_only=True,
      use_segids_in_place_of_chainids=False):
    import cctbx.geometry_restraints
    from scitbx.array_family import flex
    site_labels = flex.bool(sites_cart.size(), True).iselection()
    sorted_table, n_not_shown = proxies.get_sorted(
      by_value="residual",
      sites_cart=sites_cart,
      site_labels=site_labels,
      unit_cell=unit_cell)
    if (sorted_table is None) : return []
    outliers = []
    for restraint_info in sorted_table :
      (plane_atoms, rms_delta, residual) = restraint_info
      i_seqs = [ a[0] for a in plane_atoms ]
      deviation = max([ a[1] / a[2] for a in plane_atoms ])
      plane_atoms_ = get_atoms_info(pdb_atoms, iselection=i_seqs)
      outlier = planarity(
        atoms_info=plane_atoms_,
        rms_deltas=rms_delta,
        residual=residual,
        delta_max=max([ a[1] for a in plane_atoms ]),
        score=deviation,
        outlier=True,
        xyz=get_mean_xyz(plane_atoms_))
      if (outlier.score > sigma_cutoff):
        outliers.append(outlier)
      elif (not outliers_only):
        outlier.outlier=False
        outliers.append(outlier)
    return outliers

def get_mean_xyz(atoms):
  from scitbx.matrix import col
  sum = col(atoms[0].xyz)
  for atom in atoms[1:] :
    sum += col(atom.xyz)
  return sum / len(atoms)

def _get_sorted(O,
        unit_cell,
        sites_cart,
        pdb_atoms,
        by_value="residual",
        use_segids_in_place_of_chainids=False):
  assert by_value in ["residual", "delta"]
  if (O.size() == 0): return []
  import cctbx.geometry_restraints
  from scitbx.array_family import flex
  from cctbx.geometry_restraints.linking_class import linking_class
  origin_ids = linking_class()

  deltas = flex.abs(O.deltas(sites_cart=sites_cart))
  residuals = O.residuals(sites_cart=sites_cart)
  if (by_value == "residual"):
    data_to_sort = residuals
  elif (by_value == "delta"):
    data_to_sort = deltas
  i_proxies_sorted = flex.sort_permutation(data=data_to_sort, reverse=True)
  sorted_table = []
  for i_proxy in i_proxies_sorted:
    proxy = O[i_proxy]
    if proxy.origin_id != origin_ids.get_origin_id('covalent geometry'):
      continue
    sigma = cctbx.geometry_restraints.weight_as_sigma(proxy.weight)
    score = sqrt(residuals[i_proxy]) / sigma
    proxy_atoms = get_atoms_info(pdb_atoms, iselection=proxy.i_seqs,
      use_segids_in_place_of_chainids=use_segids_in_place_of_chainids)
    sorted_table.append((proxy, proxy_atoms))
  return sorted_table

class combined(slots_getstate_setstate):
  """
  Container for individual validations of each of the five covalent restraint
  classes.
  """
  __geo_types__ = ["bonds", "angles", "dihedrals", "chiralities", "planarities"]
  __slots__ = __geo_types__ + ["_use_cdl"]
  def __init__(self,
      pdb_hierarchy,
      xray_structure,
      geometry_restraints_manager,
      ignore_hd=True,
      sigma_cutoff=4.0,
      outliers_only=True,
      use_segids_in_place_of_chainids=False,
      cdl=None):
    self._use_cdl = cdl
    from mmtbx import restraints
    restraints_manager = restraints.manager(
      geometry=geometry_restraints_manager)
    sites_cart = xray_structure.sites_cart()
    hd_selection = xray_structure.hd_selection()
    pdb_atoms = pdb_hierarchy.atoms()
    if (ignore_hd and hd_selection.count(True) > 0):
      restraints_manager = restraints_manager.select(selection = ~hd_selection)
      sites_cart = sites_cart.select(~hd_selection)
      pdb_atoms = pdb_atoms.select(~hd_selection)
    energies_sites = restraints_manager.energies_sites(
      sites_cart=sites_cart,
      compute_gradients=False).geometry
    for geo_type in self.__geo_types__ :
      restraint_validation_class = globals()[geo_type]
      if (geo_type == "bonds" ):
        restraint_proxies = restraints_manager.geometry.pair_proxies(
          sites_cart=sites_cart).bond_proxies
      else :
        restraint_proxies = getattr(restraints_manager.geometry,
          "%s_proxies" % restraint_validation_class.restraint_type)
      rv = restraint_validation_class(
        pdb_atoms=pdb_atoms,
        sites_cart=sites_cart,
        energies_sites=energies_sites,
        restraint_proxies=restraint_proxies,
        unit_cell=xray_structure.unit_cell(),
        ignore_hd=ignore_hd,
        sigma_cutoff=sigma_cutoff,
        outliers_only=outliers_only,
        use_segids_in_place_of_chainids=use_segids_in_place_of_chainids)
      setattr(self, geo_type, rv)

  def show(self, out=sys.stdout, prefix="", verbose=True):
    for geo_type in self.__geo_types__ :
      rv = getattr(self, geo_type)
      make_sub_header(rv.restraint_label + "s", out=out)
      if (geo_type == "angles") and getattr(self, "_use_cdl", False):
        print("  Using conformation-dependent library for mainchain "+\
                      "bond angle targets", file=out)
        print("", file=out)
      rv.show(out=out, prefix=prefix)

  def get_bonds_angles_rmsds(self):
    return (self.bonds.mean, self.angles.mean)

  def as_kinemage(self, chain_id=None):
    kin_txt = self.angles.as_kinemage(chain_id=chain_id)
    kin_txt += "\n"
    kin_txt += self.bonds.as_kinemage(chain_id=chain_id)
    return kin_txt
