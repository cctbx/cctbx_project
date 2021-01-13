
from __future__ import absolute_import, division, print_function
from cctbx import maptbx, miller
from cctbx.sgtbx import space_group_info
from iotbx.file_reader import any_file
from mmtbx.validation import atom, atom_info, validation
from mmtbx.validation import experimental
from libtbx.str_utils import format_value
import sys

class water(atom):
  """
  Container for information about a water atom, including electron density
  properties.
  """

  __slots__ = atom.__slots__ + experimental.__real_space_attr__ + [
    "anom",
    "nearest_contact",
    "nearest_atom",
    "n_hbonds",
    "fmodel"
  ]

  @property
  def cc(self):
    return self.score

  @staticmethod
  def header():
    return "%-20s  %6s  %4s  %6s  %6s  %6s  %5s" % ("atom", "b_iso", "occ",
      "2Fo-Fc", "Fo-Fc", "Anom", "CC")

  def id_str(self):
    return "%2s%4s%1s%1s %-4s" % (self.chain_id, self.resseq, self.icode,
      self.altloc, self.name)

  def as_string(self, prefix="", highlight_if_heavy=False):
    flag = ""
    if highlight_if_heavy and self.is_heavy_atom():
      flag = " ***"
    return "%-20s  %6.2f  %4.2f  %6.2f  %6.2f  %s  %5.3f%s" % (
      self.id_str(), self.b_iso, self.occupancy, self.two_fofc, self.fofc,
      format_value("%6.2f", self.anom, replace_none_with="---"), self.cc, flag)

  def is_bad_water(self):
    return ((self.cc < 0.8) or (self.fofc < -3.0) or (self.occupancy < 0.5)
            or (self.two_fofc < 1.0))

  def is_heavy_atom(self):
    return ((self.fofc is not None and self.fofc > 3.0)
            or (self.anom is not None and self.anom > 3.0)
            or (self.b_iso is not None and self.b_iso == 0))

  def as_table_row_phenix(self):
    return [ self.id_str(), self.b_iso, self.occupancy, self.two_fofc,
             self.fmodel, self.score]

class waters(validation):
  """
  Assess the properties of solvent atoms, including local environment and
  electron density.
  """

  __slots__ = validation.__slots__ + ["n_bad", "n_heavy"]

  def get_result_class(self) : return water

  def __init__(self, pdb_hierarchy, xray_structure, fmodel,
                distance_cutoff=4.0, collect_all=True,
                molprobity_map_params=None):
    validation.__init__(self)
    from mmtbx.real_space_correlation import extract_map_stats_for_single_atoms
    from cctbx import adptbx
    from scitbx.matrix import col
    self.n_bad = 0
    self.n_heavy = 0
    pdb_atoms = pdb_hierarchy.atoms()
    if(len(pdb_atoms)>1):
      assert (not pdb_atoms.extract_i_seq().all_eq(0))
    unit_cell = xray_structure.unit_cell()
    pair_asu_table = xray_structure.pair_asu_table(
      distance_cutoff = distance_cutoff)
    asu_mappings = pair_asu_table.asu_mappings()
    asu_table = pair_asu_table.table()
    u_isos = xray_structure.extract_u_iso_or_u_equiv()
    occupancies = xray_structure.scatterers().extract_occupancies()
    sites_frac = xray_structure.sites_frac()
    sel_cache = pdb_hierarchy.atom_selection_cache()
    water_sel = sel_cache.selection("water")

    if (molprobity_map_params is not None):
      # assume parameters have been validated (symmetry of pdb and map matches)
      two_fofc_map = None
      fc_map = None
      d_min = None
      crystal_gridding = None

      # read two_fofc_map
      if (molprobity_map_params.map_file_name is not None):
        f = any_file(molprobity_map_params.map_file_name)
        two_fofc_map = f.file_object.map_data()
        d_min = molprobity_map_params.d_min
        crystal_gridding = maptbx.crystal_gridding(
          f.file_object.unit_cell(),
          space_group_info=space_group_info(f.file_object.space_group_number),
          pre_determined_n_real=f.file_object.unit_cell_grid)

        pdb_atoms = pdb_hierarchy.atoms()
        xray_structure = pdb_hierarchy.extract_xray_structure(
          crystal_symmetry=f.crystal_symmetry())
        unit_cell = xray_structure.unit_cell()
        # check for origin shift
        # ---------------------------------------------------------------------
        soin = maptbx.shift_origin_if_needed(
          map_data         = two_fofc_map,
          sites_cart       = xray_structure.sites_cart(),
          crystal_symmetry = xray_structure.crystal_symmetry())
        two_fofc_map   = soin.map_data
        xray_structure.set_sites_cart(soin.sites_cart)
        # ---------------------------------------------------------------------
        pair_asu_table = xray_structure.pair_asu_table(
          distance_cutoff = distance_cutoff)
        asu_mappings = pair_asu_table.asu_mappings()
        asu_table = pair_asu_table.table()
        u_isos = xray_structure.extract_u_iso_or_u_equiv()
        occupancies = xray_structure.scatterers().extract_occupancies()
        sites_frac = xray_structure.sites_frac()
        sel_cache = pdb_hierarchy.atom_selection_cache()
        water_sel = sel_cache.selection("water")

      elif (molprobity_map_params.map_coefficients_file_name is not None):
        f = any_file(molprobity_map_params.map_coefficients_file_name)
        fourier_coefficients = f.file_server.get_miller_array(
          molprobity_map_params.map_coefficients_label)
        crystal_symmetry = fourier_coefficients.crystal_symmetry()
        d_min = fourier_coefficients.d_min()
        crystal_gridding = maptbx.crystal_gridding(
          crystal_symmetry.unit_cell(), d_min, resolution_factor=0.25,
          space_group_info=crystal_symmetry.space_group_info())
        two_fofc_map = miller.fft_map(
          crystal_gridding=crystal_gridding,
          fourier_coefficients=fourier_coefficients).apply_sigma_scaling().\
          real_map_unpadded()

      # calculate fc_map
      assert( (d_min is not None) and (crystal_gridding is not None) )
      f_calc = xray_structure.structure_factors(d_min=d_min).f_calc()
      fc_map = miller.fft_map(crystal_gridding=crystal_gridding,
                              fourier_coefficients=f_calc)
      fc_map = fc_map.apply_sigma_scaling().real_map_unpadded()

      map_stats = extract_map_stats_for_single_atoms(
        pdb_atoms=pdb_atoms,
        xray_structure=xray_structure,
        fmodel=None,
        selection=water_sel,
        fc_map=fc_map,
        two_fofc_map=two_fofc_map)
    else:
      map_stats = extract_map_stats_for_single_atoms(
        pdb_atoms=pdb_atoms,
        xray_structure=xray_structure,
        fmodel=fmodel,
        selection=water_sel)
    waters = []
    for i_seq, atom in enumerate(pdb_atoms):
      if (water_sel[i_seq]):
        rt_mx_i_inv = asu_mappings.get_rt_mx(i_seq, 0).inverse()
        self.n_total += 1
        asu_dict = asu_table[i_seq]
        nearest_atom = nearest_contact = None
        for j_seq, j_sym_groups in asu_dict.items():
          atom_j = pdb_atoms[j_seq]
          site_j = sites_frac[j_seq]
          # Filter out hydrogens
          if atom_j.element.upper().strip() in ["H", "D"]:
            continue
          for j_sym_group in j_sym_groups:
            rt_mx = rt_mx_i_inv.multiply(asu_mappings.get_rt_mx(j_seq,
              j_sym_group[0]))
            site_ji = rt_mx * site_j
            site_ji_cart = xray_structure.unit_cell().orthogonalize(site_ji)
            vec_i = col(atom.xyz)
            vec_ji = col(site_ji_cart)
            dxyz = abs(vec_i - vec_ji)
            if (nearest_contact is None) or (dxyz < nearest_contact):
              nearest_contact = dxyz
              nearest_atom = atom_info(pdb_atom=atom_j, symop=rt_mx)
        w = water(
          pdb_atom=atom,
          b_iso=adptbx.u_as_b(u_isos[i_seq]),
          occupancy=occupancies[i_seq],
          nearest_contact=nearest_contact,
          nearest_atom=nearest_atom,
          score=map_stats.two_fofc_ccs[i_seq],
          fmodel=map_stats.fmodel_values[i_seq],
          two_fofc=map_stats.two_fofc_values[i_seq],
          fofc=map_stats.fofc_values[i_seq],
          anom=map_stats.anom_values[i_seq],
          n_hbonds=None) # TODO
        if (w.is_bad_water()):
          w.outlier = True
          self.n_bad += 1
        elif (w.is_heavy_atom()):
          w.outlier = True
          self.n_heavy += 1
        if (w.outlier) or (collect_all):
          self.results.append(w)
    self.n_outliers = len(self.results)

  def show_summary(self, out=sys.stdout, prefix="  "):
    if (self.n_bad > 0):
      print("%sPoorly ordered waters:  %4d" % (prefix, self.n_bad), file=out)
    if (self.n_heavy > 0):
      print("%sMislabeled waters:      %4d" % (prefix, self.n_heavy), file=out)
    if (self.n_bad == 0) and (self.n_heavy == 0):
      print("%sAll waters okay." % prefix, file=out)

  def show(self, out=sys.stdout, prefix="  ", verbose=True):
    if (len(self.results) > 0):
      if (self.n_bad > 0):
        print(prefix + "Waters in poor density:", file=out)
        print(prefix + self.get_result_class().header(), file=out)
        for result in self.results :
          if (result.is_bad_water()):
            print(prefix + str(result), file=out)
        print("", file=out)
      if (self.n_heavy > 0):
        print(prefix + "Possibly mislabeled atoms:", file=out)
        print(prefix + self.get_result_class().header(), file=out)
        for result in self.results :
          if (result.is_heavy_atom()):
            print(prefix + str(result), file=out)
        print("", file=out)
    self.show_summary(out=out, prefix=prefix)
