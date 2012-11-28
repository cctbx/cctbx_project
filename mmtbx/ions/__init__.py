"""
Examines a structure for metal ions. Can iterate over all atoms, examining their density
and chemical environment to determine if they are correctly identified, or if there are
better candidate ions that they can be replaced with.

See build.build_ions to actually modify the structure, code in this module only prints out
messages to the log.
"""

from __future__ import division
from cctbx import crystal, adptbx
from cctbx.eltbx import sasaki, henke
from scitbx.matrix import col
from scitbx.array_family import flex
import libtbx.phil
import libtbx.load_env
from libtbx import group_args, adopt_init_args, Auto
from libtbx.str_utils import make_sub_header
from libtbx.utils import null_out, Sorry
from libtbx import easy_mp

from mmtbx.ions.parameters import get_charge, server, MetalParameters
from mmtbx.ions.geometry import find_coordination_geometry
from math import sqrt
import cStringIO
import sys

chloride_params_str = """
max_distance_to_amide_n = 3.5
  .type = float
max_distance_to_cation = 3.5
  .type = float
min_distance_to_anion = 3.5
  .type = float
min_distance_to_other_sites = 2.1
  .type = float
delta_amide_h_angle = 20
  .type = float
delta_planar_angle = 10
  .type = float
radius = 2.0
  .type = float
"""

ion_master_phil = libtbx.phil.parse("""
refine_occupancies = True
  .type = bool
starting_occupancy = 0.9
  .type = float
use_phaser = True
  .type = bool
require_valence = False
  .type = bool
ambiguous_valence_cutoff = 0.5
  .type = float
  .help = "Cutoff to uniquely select one of many metals by comparing its observed and expected valence"

water {
  min_2fofc_level = 1.8
    .type = float
  max_fofc_level = 3.0
    .type = float
  max_anom_level = 4.0
    .type = float
  max_occ = 1.0
    .type = float
  max_stddev_b_iso = 3
    .type = float
  min_frac_b_iso = 0.25
    .type = float
}
chloride {
  %s
}
phaser {
  distance_cutoff = 1.5
    .type = float
  distance_cutoff_same_site = 0.5
    .type = float
  fpp_ratio_min = 0.4
    .type = float
  use_llg_anom_map = False
    .type = bool
    .help = If True, the anomalous LLG map from Phaser will be used instead \
      of a conventional anomalous difference map.
}
""" % (chloride_params_str))

WATER_RES_NAMES = ["HOH", "WAT"]
DEFAULT_IONS = ["MG", "CA", "ZN", "CL"]
HALIDES = ["F", "CL", "BR", "I"]

# Signals a built water might be a ion:
# - Abnormal b-factors from nearby chain
# - Coordinated by other waters
# - Nearby carbon (< 3.0 A)
# - FoFc peak/hole
# - Phaser Anomalous signal
class Manager (object):
  """
  Wrapper object for extracting local environment and adding or modifying
  scatterers.
  """

  def __init__ (self,
                fmodel,
                pdb_hierarchy,
                xray_structure,
                params = None,
                wavelength = None,
                connectivity = None,
                nproc = 1,
                verbose = False,
                log = None) :
    self.fmodel = fmodel
    self.params = params
    self.wavelength = wavelength
    self.server = server()
    self.map_anom = None # optional
    self.nproc = nproc
    self.phaser_substructure = None
    self.fpp_from_phaser_ax_sites = None
    self.update_structure(
      pdb_hierarchy=pdb_hierarchy,
      xray_structure=xray_structure,
      connectivity=connectivity)
    self.update_maps()
    if ((params.use_phaser) and
        (libtbx.env.has_module("phaser")) and
        (fmodel.f_obs().anomalous_flag())) :
      self.phaser_substructure = find_anomalous_scatterers(
        fmodel, pdb_hierarchy, wavelength=wavelength).atoms()
      self.analyze_substructure(log=log, verbose=verbose)

  def update_structure (self, pdb_hierarchy, xray_structure,
      connectivity=None, log=None) :
    """
    Set the current atomic data: PDB hierarchy, Xray structure, and simple
    connectivity list.
    """
    self.pdb_hierarchy = pdb_hierarchy
    self.xray_structure = xray_structure
    self.fmodel.update_xray_structure(xray_structure, update_f_calc=True)
    self.sites_frac = self.xray_structure.sites_frac()
    self.connectivity = connectivity
    self.pdb_atoms = pdb_hierarchy.atoms()
    self.unit_cell = xray_structure.unit_cell()
    self.ax_chain = None

    # Extract some information about the structure, including B-factor
    # statistics for non-HD atoms and waters
    sctr_keys = xray_structure.scattering_type_registry().type_count_dict().keys()
    self.hd_present = ("H" in sctr_keys) or ("D" in sctr_keys)
    sel_cache = pdb_hierarchy.atom_selection_cache()
    water_sel = sel_cache.selection("(%s) and name O" %
                                    (" or ".join("resname " + i for i in WATER_RES_NAMES)))
    not_hd_sel = sel_cache.selection("not (element H or element D)")

    self.n_waters = len(water_sel)
    self.n_heavy = len(not_hd_sel)
    u_iso_all = xray_structure.extract_u_iso_or_u_equiv()
    u_iso_all_tmp = u_iso_all.select(not_hd_sel)
    self.b_mean_all = adptbx.u_as_b(flex.mean(u_iso_all_tmp))
    self.b_stddev_all = adptbx.u_as_b(
       u_iso_all_tmp.standard_deviation_of_the_sample())
    self.b_mean_hoh = self.b_stddev_hoh = None
    if (self.n_waters > 0) :
      u_iso_hoh = u_iso_all.select(water_sel)
      self.b_mean_hoh = adptbx.u_as_b(flex.mean(u_iso_hoh))
      self.b_stddev_hoh = adptbx.u_as_b(
        u_iso_hoh.standard_deviation_of_the_sample())
    if (self.phaser_substructure is not None) :
      self.analyze_substructure(log=log)

  def update_maps (self) :
    """
    Generate new maps for the current structure, including anomalous map if
    data are anomalous.

    If Phaser is installed and anomalous data are available, the anomalous
    log-likelihood gradient (LLG) map can be used instead of the conventional
    anomalous difference map. This is only really useful if the anomalous
    scattering of existing atoms is modeled (and ideally, refined).
    """
    def fft_map (map_coeffs, resolution_factor = 0.25):
      return map_coeffs.fft_map(resolution_factor = resolution_factor,
        ).apply_sigma_scaling().real_map_unpadded()
    f_map_coeffs = self.fmodel.map_coefficients(map_type="2mFo-DFc")
    df_map_coeffs = self.fmodel.map_coefficients(map_type="mFo-DFc")
    self.map_2fofc = fft_map(f_map_coeffs)
    self.map_fofc = fft_map(df_map_coeffs)
    self.anom_map = None
    if (self.fmodel.f_obs().anomalous_flag()) :
      if ((self.params.phaser.use_llg_anom_map) and
          (libtbx.env.has_module("phaser"))) :
        import mmtbx.map_tools
        anom_map_coeffs = mmtbx.map_tools.get_phaser_sad_llg_map_coefficients(
          fmodel=self.fmodel,
          pdb_hierarchy=self.pdb_hierarchy,
          log=None)
      else :
        anom_map_coeffs = self.fmodel.map_coefficients(map_type="anom")
      self.map_anom = fft_map(anom_map_coeffs)

  def find_nearby_atoms (self, i_seq, distance_cutoff = 3.0,
      filter_by_bonding = True):
    """
    Given site in the structure, return a list of nearby atoms with the
    supplied cutoff, and the vectors between them and the atom's site. Takes
    into account symmetry operations when finding nearby sites.
    """

    assert (i_seq < len(self.sites_frac))

    # Use pair_asu_table to find atoms within distance_cutoff of one another,
    # taking into account potential cell symmetry.
    # XXX should probably move this up, but it seems relatively fast
    unit_cell = self.xray_structure.unit_cell()
    pair_asu_table = self.xray_structure.pair_asu_table(
      distance_cutoff = distance_cutoff)
    asu_mappings = pair_asu_table.asu_mappings()
    asu_table = pair_asu_table.table()
    asu_dict = asu_table[i_seq]
    contacts = []
    site_i = self.sites_frac[i_seq]
    rt_mx_i_inv = asu_mappings.get_rt_mx(i_seq, 0).inverse()

    for j_seq, j_sym_groups in asu_dict.items():
      site_j = self.sites_frac[j_seq]

      for j_sym_group in j_sym_groups:
        rt_mx = rt_mx_i_inv.multiply(asu_mappings.get_rt_mx(j_seq,
          j_sym_group[0]))
        site_ji = rt_mx * site_j

        vec_i, vec_ji = col(unit_cell.orthogonalize(site_i)), \
                        col(unit_cell.orthogonalize(site_ji))
        assert abs(vec_i - vec_ji) < distance_cutoff + 0.5
        contacts.append((self.pdb_atoms[j_seq], vec_i - vec_ji))

    # Filter out atoms that are judged to be "contacts", but actually involved
    # in bidentate coordination of the target atom.  Currently this only
    # handles carboxyl groups.
    if filter_by_bonding and self.connectivity:
      filtered = []
      all_i_seqs = [atom.i_seq for atom, vector in contacts]

      for atom, vector in contacts:
        if _get_element(atom) not in ["C"]:
          filtered += (atom, vector),
          continue

        bonded_j_seqs = []
        for j_seq in self.connectivity[atom.i_seq] :
          if (j_seq in all_i_seqs) :
            bonded_j_seqs.append(j_seq)
        keep = True

        for j_seq in bonded_j_seqs:
          other_atom, other_vector = contacts[all_i_seqs.index(j_seq)]

          if _get_element(other_atom) in ["N", "O"]:
            if abs(other_vector) < abs(vector):
              keep = False
              break

        if keep:
          filtered += (atom, vector),

      contacts = filtered

    return contacts

  def principal_axes_of_inertia (self, site_cart, radius) :
    """
    Extracts the map grid points around a site, and calculates the axes of
    inertia (using the density values as weights).
    """
    # XXX not totally confident about how I'm weighting this...
    from cctbx import maptbx
    return maptbx.principal_axes_of_inertia(
      real_map=self.map_2fofc,
      site_cart=site_cart,
      unit_cell=self.xray_structure.unit_cell(),
      radius=radius)

  def get_map_sphere_variance (self, site_cart, radius,
      use_difference_map=False) :
    """
    Calculate the density levels for points on a sphere around a given point.
    This will give us some indication whether the density falls off in the
    expected manner around a single atom.
    """
    from cctbx import maptbx
    real_map = self.map_2fofc
    if (use_difference_map) :
      real_map = self.map_fofc
    return maptbx.spherical_variance_around_point(
      real_map=real_map,
      unit_cell=self.xray_structure.unit_cell(),
      site_cart=site_cart,
      radius=radius)

  def map_stats (self, i_seq) :
    """
    Given a site in the structure, find the signal of the 2FoFc, FoFc, and
    anomalous map (When available).
    """

    site = self.sites_frac[i_seq]
    return self.map_stats_at_point(site)

  def map_stats_at_point (self, site_frac) :
    """
    Given a fractional coordinate, find the signal of the 2FoFc, FoFc, and
    anomalous map (When available).
    """

    assert (self.map_2fofc is not None)
    assert (self.map_fofc is not None)

    value_2fofc = self.map_2fofc.tricubic_interpolation(site_frac)
    value_fofc = self.map_fofc.tricubic_interpolation(site_frac)
    value_anom = None

    if (self.map_anom is not None) :
      value_anom = self.map_anom.tricubic_interpolation(site_frac)

    return group_args(
      two_fofc=value_2fofc,
      fofc=value_fofc,
      anom=value_anom)

  def find_atoms_near_site (self,
      site_cart,
      distance_cutoff = 1.5,
      distance_cutoff_same_site = 0.5) :
    """
    Given an XYZ coordinate, finds atoms near that site, separating out those
    which are close enough to be essentially equivalent.  Used to analyze the
    anomalously-scattering substructure calculated by Phaser.

    Returns same_atoms, a list of atoms within distance_cutoff_same_site of
    site_cart, and and other_atoms, the rest of the atoms within
    distance_cutoff of site_cart.
    """

    site_frac = self.xray_structure.unit_cell().fractionalize(
      site_cart = site_cart)
    sites_frac = flex.vec3_double([site_frac])
    asu_mappings = self.xray_structure.asu_mappings(buffer_thickness =
      distance_cutoff + 0.1)
    asu_mappings.process_sites_frac(sites_frac,
      min_distance_sym_equiv = self.xray_structure.min_distance_sym_equiv())
    pair_generator = crystal.neighbors_fast_pair_generator(
      asu_mappings = asu_mappings,
      distance_cutoff = distance_cutoff)
    n_xray = self.xray_structure.scatterers().size()

    same_atoms = []
    other_atoms = []

    # Iterate through all interacting pairs of atoms in the structure
    # within distance_cutoff of each other
    for pair in pair_generator:
      # Find a pair where one's sequence ID is < n_xray and the other's
      # is >= n_xray and gather the site information about it
      if pair.i_seq < n_xray:
        if pair.j_seq < n_xray:
          continue

        rt_mx_i = asu_mappings.get_rt_mx_i(pair)
        rt_mx_j = asu_mappings.get_rt_mx_j(pair)
        rt_mx_ji = rt_mx_i.inverse().multiply(rt_mx_j)
        new_site_frac = rt_mx_ji * site_frac
        new_site_cart = self.unit_cell.orthogonalize(site_frac = new_site_frac)

        site_info = group_args(
          i_seq = pair.i_seq,
          rt_mx = rt_mx_ji,
          new_site_cart = new_site_cart,
          distance = sqrt(pair.dist_sq))
      else:
        if pair.j_seq >= n_xray:
          continue

        rt_mx_i = asu_mappings.get_rt_mx_i(pair)
        rt_mx_j = asu_mappings.get_rt_mx_j(pair)
        rt_mx_ij = rt_mx_j.inverse().multiply(rt_mx_i)
        new_site_frac = rt_mx_ij * site_frac
        new_site_cart = self.unit_cell.orthogonalize(site_frac = new_site_frac)

        site_info = group_args(
          i_seq = pair.j_seq,
          rt_mx = rt_mx_ij,
          new_site_cart = new_site_cart,
          distance = sqrt(pair.dist_sq))

      if site_info:
        if site_info.distance < distance_cutoff_same_site:
          same_atoms.append(site_info)
        else:
          other_atoms.append(site_info)

    same_atoms.sort(key = lambda x: x.distance)
    other_atoms.sort(key = lambda x: x.distance)

    return same_atoms, other_atoms

  def analyze_substructure (self, log=None, verbose=True) :
    """
    Given a list of AX pseudo-atoms placed by Phaser, finds the nearest real
    atoms, and if possible determines their equivalence.
    """
    assert (self.phaser_substructure is not None)
    self.fpp_from_phaser_ax_sites = flex.double(self.pdb_atoms.size(), -1)
    if (log is None) or (not verbose) : log = null_out()

    make_sub_header("Analyzing Phaser anomalous substructure", out=log)
    for atom in self.phaser_substructure :
      same_atoms, other_atoms = self.find_atoms_near_site(
        atom.xyz,
        distance_cutoff=self.params.phaser.distance_cutoff,
        distance_cutoff_same_site=self.params.phaser.distance_cutoff_same_site)

      def ax_atom_id (atom) :
        return "AX %s (fpp=%.3f)" % (atom.serial.strip(), atom.occ)

      if (len(same_atoms) == 0) :
        print >> log, "  No match for %s" % ax_atom_id(atom)
        for site_info in other_atoms :
          print >> log, "    %s (distance = %.3f)" % (
            self.pdb_atoms[site_info.i_seq].id_str(), site_info.distance)
      elif (len(same_atoms) == 1) :
        print >> log, "  %s maps to %s" % (ax_atom_id(atom),
          self.pdb_atoms[same_atoms[0].i_seq].id_str())
        self.fpp_from_phaser_ax_sites[same_atoms[0].i_seq] = atom.occ
      else :
        print >> log, "  ambiguous results for %s:" % ax_atom_id(atom)
        for a in same_atoms :
          print >> log, "    %s" % self.pdb_atoms[a.i_seq].id_str()
    print >> log, ""

  def get_fpp (self, i_seq) :
    """
    Retrieve the refined f'' for a site.
    """
    if (self.fpp_from_phaser_ax_sites is None) : return None
    fpp = self.fpp_from_phaser_ax_sites[i_seq]
    if (fpp < 0) : fpp = None
    return fpp

  # XXX approximate guess, needs more rational parameters
  def looks_like_halide_ion (self,
      i_seq,
      element = "CL",
      assume_hydrogens_all_missing = True):
    """
    Given a site, analyze the nearby atoms (taking geometry into account) to
    guess whether it might be a halide ion.  Among other things, halides will
    often be coordinated by amide hydrogens, and not in close contact with
    negatively charged atoms.  Note that this procedure has a very high
    false positive rate when used by itself, so additional information about
    the electron count (map, occ, b) is absolutely essential.
    """
    atom = self.pdb_atoms[i_seq]
    assert element.upper() in HALIDES
    # discard atoms with B-factors greater than mean-1sigma for waters
    if ((self.b_mean_hoh is not None) and
        (atom.b > self.b_mean_hoh - self.b_stddev_hoh)) :
      return False

    params = self.params.chloride
    nearby_atoms = self.find_nearby_atoms(
      i_seq,
      distance_cutoff = 4.0,
      filter_by_bonding = False)

    if assume_hydrogens_all_missing:
      xrs = self.xray_structure
      sctr_keys = xrs.scattering_type_registry().type_count_dict().keys()
      assume_hydrogens_all_missing = "H" not in sctr_keys and \
                                     "D" not in sctr_keys

    binds_backbone_amide = False
    near_cation = False
    near_arg_lys = False
    xyz = col(atom.xyz)

    fpp = self.get_fpp(i_seq)
    good_fpp = None
    if (fpp is not None) and (fpp > 0.1) :
      fpp_expected = sasaki.table(element).at_angstrom(
        self.wavelength).fdp()
      if (fpp_expected != 0) :
        fpp_ratio = fpp / fpp_expected
        if (fpp_ratio > 0.4) and (fpp_ratio < 1.05) :
          good_fpp = True
        else :
          return False

    for atom, vector in nearby_atoms:
      # to analyze local geometry, we use the target site mapped to be in the
      # same ASU as the interacting site
      resname = atom.fetch_labels().resname.strip().upper()
      atom_name = atom.name.strip()
      element = _get_element(atom)

      # XXX need to figure out exactly what this should be - CL has a
      # fairly large radius though (1.67A according to ener_lib.cif)
      if abs(vector) < params.min_distance_to_other_sites:
        return False

      if not element in ["C","N","H","O","S"]:
        charge = get_charge(element)

        if charge < 0 and abs(vector) <= params.min_distance_to_anion:
          # Nearby anion that is too close
          return False

        if charge > 0 and abs(vector) <= params.max_distance_to_cation:
          # Nearby cation
          near_cation = True
      elif ((abs(vector) < 3.2) and
            (atom_name in ["OD1","OD2","OE1","OE2"]) and
            (resname in ["GLU","ASP"])) :
        # Negatively charged sidechains
        # XXX this really needs to be more exhaustive (e.g. phosphate oxygens)
        return False

      elif (atom_name in ["NZ", "NE"] and
            resname in ["ARG","LYS"] and
            abs(vector) <= params.max_distance_to_cation):
        # Positively charged sidechains.
        # XXX actually, shouldn't it also be roughly coplanar with the
        # guanidinium group in Arg?
        near_arg_lys = True

      elif atom_name in ["H"]:
        # Backbone amide, explicit H
        # XXX can we make this more general for any amide H?
        xyz_h = col(atom.xyz)
        bonded_atoms = self.connectivity[i_seq]
        assert (len(bonded_atoms) == 1)
        xyz_n = col(self.pdb_atoms[bonded_atoms[0]].xyz)
        vec_hn = xyz_h - xyz_n
        vec_hx = xyz_h - xyz
        angle = abs(vec_hn.angle(vec_hx, deg = True))

        # If Cl, H, and N line up, Cl binds the amide group
        if abs(angle - 180) <= params.delta_amide_h_angle:
          binds_backbone_amide = True

      elif atom_name in ["N"] and assume_hydrogens_all_missing:
        # Backbone amide, implicit H
        xyz_n = col(atom.xyz)
        bonded_atoms = self.connectivity[atom.i_seq]
        ca_same = c_prev = None

        for j_seq in bonded_atoms :
          other = self.pdb_atoms[j_seq]
          if other.name.strip().upper() in ["CA"]:
            ca_same = col(other.xyz)
          elif other.name.strip().upper() in ["C"]:
            c_prev = col(other.xyz)

        if ca_same is not None and c_prev is not None:
          xyz_cca = (ca_same + c_prev) / 2
          vec_ncca = xyz_n - xyz_cca
          # 0.86 is the backbone N-H bond distance in geostd
          xyz_h = xyz_n + (vec_ncca.normalize() * 0.86)
          vec_nh = xyz_n - xyz_h
          vec_nx = xyz_n - xyz
          angle = abs(vec_nh.angle(vec_nx, deg = True))

          if abs(angle - 180) <= params.delta_amide_h_angle:
            binds_backbone_amide = True

      elif ((atom_name in ["HD1","HD2"] and resname in ["ASN"]) or
            (atom_name in ["HE1","HE2"] and resname in ["GLN"])):
        # ASN/GLN sidechain amide, explicit H
        bonded_atoms = self.connectivity[atom.i_seq]
        assert (len(bonded_atoms) == 1)
        xyz_n = col(self.pdb_atoms[bonded_atoms[0]].xyz)
        xyz_h = col(atom.xyz)
        vec_nh = xyz_n - xyz_h
        vec_xh = xyz - xyz_h
        angle = abs(vec_nh.angle(vec_xh, deg=True))

        if abs(angle - 180) <= params.delta_amide_h_angle:
          binds_sidechain_amide = True

      elif (assume_hydrogens_all_missing and
            ((atom.name in ["ND2"] and resname in ["ASN"]) or
             (atom.name in ["NE2"] and resname in ["GLN"]))):
        # ASN/GLN sidechain amide, implicit H
        bonded_atoms = self.connectivity[atom.i_seq]
        assert (len(bonded_atoms) == 1)
        c_i_seq = bonded_atoms[0]
        c_bonded_atoms = self.connectivity[c_i_seq]
        xyz_o = None

        # Grab the amide group's oxygen coordinate
        for j_seq in c_bonded_atoms:
          other = self.pdb_atoms[j_seq]
          if other.name in ["OD1", "OE1"]:
            xyz_o = col(other.xyz)
            break
        else:
          continue

        # Find the angle for X-N-C and N-C-O
        xyz_c = col(self.pdb_atoms[c_i_seq].xyz)
        xyz_n = col(atom.xyz)
        vec_nc = xyz_n - xyz_c
        vec_xn = xyz - xyz_n
        vec_oc = xyz_o - xyz_c
        angle_xnc = abs(vec_nc.angle(vec_xn, deg = True))
        angle_nco = abs(vec_nc.angle(vec_oc, deg = True))

        # 120 degrees is the angle between CG-ND2-HD* in ASN;
        # the ion should also be approximately co-planar
        if (abs(angle_xnc - 120) <= params.delta_amide_h_angle and
            abs(angle_nco % 180) <= params.delta_planar_angle):
          binds_sidechain_amide = True

    # TODO something smart...
    # the idea here is to determine whether the blob of density around the
    # site is approximately spherical and limited in extent.
    pai = self.principal_axes_of_inertia(atom.xyz,
      radius=self.params.chloride.radius)
    #print list(pai.eigensystem().values())
    map_variance = self.get_map_sphere_variance(atom.xyz,
      radius=self.params.chloride.radius)
    map_variance.show(prefix="  ")
    good_map_falloff = False
    if ((map_variance.mean < 1.0) and
        (map_variance.standard_deviation < 0.4)) : # XXX very arbitrary
      good_map_falloff = True

    # XXX probably need something more sophisticated here too.  a CL which
    # coordinates a metal (e.g. in 4aqi) may not have a clear blob, but
    # will still be detectable by other criteria.
    return binds_backbone_amide or near_cation or near_arg_lys or \
      (good_fpp and good_map_falloff)

  def check_water_properties (self, atom_props):
    """
    Checks the atom characteristics against what we would expect for a water.
    Updates atom_props with any inaccuracies noticed (Surpringly low b-factor,
    high occupancy, etc).
    """

    params = self.params.water
    inaccuracies = atom_props.inaccuracies["HOH"] = set()

    # Skip over water if the 2FoFc peak is low
    if atom_props.peak_2fofc < params.min_2fofc_level:
      return None

    # waters that have B-factors above the mean are presumed correct
    # XXX this may be problematic for some of our testing!
    if (atom_props.atom.b > self.b_mean_hoh) :
      return None

    if atom_props.peak_fofc > params.max_fofc_level:
      inaccuracies.add(atom_props.FOFC_PEAK)

    if atom_props.atom.occ > params.max_occ:
      inaccuracies.add(atom_props.HIGH_OCC)

    if (atom_props.peak_anom is not None and
        atom_props.peak_anom > params.max_anom_level):
      inaccuracies.add(atom_props.ANOM_PEAK)

    # Check the b-factor is within a specified number of standard
    # deviations from the mean and above a specified fraction of
    # the mean.
    z_value = (atom_props.atom.b - self.b_mean_hoh) / self.b_stddev_hoh

    if z_value < -params.max_stddev_b_iso:
      inaccuracies.add(atom_props.LOW_B)
    elif atom_props.atom.b < self.b_mean_hoh * params.min_frac_b_iso:
      inaccuracies.add(atom_props.LOW_B)

    return atom_props.is_correctly_identified()

  def analyze_water(self,
      i_seq,
      show_only_map_outliers = True,
      debug = True,
      candidates = Auto) :
    """
    Examines the environment around a single atom to determine if it is actually
    a misidentified metal ion.
    """
    atom = self.pdb_atoms[i_seq]
    if (candidates is Auto) :
      candidates = DEFAULT_IONS

    candidates = [i.strip().upper() for i in candidates]

    resname = atom.fetch_labels().resname.strip().upper()
    assert resname in WATER_RES_NAMES

    atom_props = AtomProperties(atom=atom, i_seq=i_seq, manager=self)
    final_choice = None
    # Gather some quick statistics on the atom and see if it looks like water
    looks_like_water = self.check_water_properties(atom_props)
    if (looks_like_water is None) : # not trustworthy, skip
      return None

    # Filter out metals based on whether they are more or less eletron-dense
    # in comparison with water
    filtered_candidates = []
    nuc_phosphate_site = atom_props.looks_like_nucleotide_phosphate_site()
    for symbol in candidates :
      elem = self.server.get_metal_parameters(symbol)
      if (elem is None) :
        if (symbol in HALIDES) : # halides are special!
          continue
        else :
          raise Sorry("Element '%s' not supported!" % symbol)
      if (nuc_phosphate_site) and (not symbol in ["MG", "CA"]) :
        continue
      atom_number = sasaki.table(symbol.upper()).atomic_number()
      if (((looks_like_water) and (atom_number < 12)) or
          ((not looks_like_water) and (atom_number >= 12))) :
        filtered_candidates.append(elem)
    if (len(filtered_candidates) == 0) :
      return None

    # Try each different ion and see what is reasonable
    def try_candidates (require_valence=True) :
      reasonable = []
      for elem_params in filtered_candidates:
        # Try the ion with check_ion()
        atom_props.check_ion(
          ion_params=elem_params,
          server=self.server,
          wavelength=self.wavelength,
          require_valence=require_valence)
        identity = str(elem_params)
        if atom_props.is_correctly_identified(identity) :
          reasonable.append((elem_params, atom_props.score[identity]))
      return reasonable

    # first try with bond valence included in criteria - then if that doesn't
    # yield any hits, try again without valences
    reasonable = try_candidates(require_valence=True)
    valence_used = True
    if (len(reasonable) == 0) and (not self.params.require_valence) :
      reasonable = try_candidates(require_valence=False)
      if (len(reasonable) > 0) :
        valence_used = False

    looks_like_halide = False
    if not reasonable and ((not looks_like_water) or (nuc_phosphate_site)):
      # try halides now
      candid_halides = set(candidates).intersection(HALIDES)
      filtered_halides = \
        [i for i in candid_halides
         if self.looks_like_halid_ion(i_seq = i_seq, element = i)]

      looks_like_halide = len(filtered_halides) > 0
      reasonable += [(MetalParameters(element = halide, charge = -1), 0)
                     for halide in filtered_halides]

    if len(reasonable) == 1:
      elem_params = reasonable[0][0]
      final_choice = elem_params

    return water_result(
      atom_props = atom_props,
      filtered_candidates = filtered_candidates,
      matching_candidates = reasonable,
      nuc_phosphate_site = nuc_phosphate_site,
      looks_like_water = looks_like_water,
      looks_like_halide = looks_like_halide,
      ambiguous_valence_cutoff = self.params.ambiguous_valence_cutoff,
      valence_used = valence_used,
      final_choice = final_choice)

  def analyze_waters (self, out = sys.stdout, debug = True,
      show_only_map_outliers = True, candidates = Auto):
    """
    Iterates through all of the waters in a model, examining the maps and local
    environment to check their identity and suggest likely ions where
    appropriate.
    """
    model = self.pdb_hierarchy.only_model()

    waters = []
    for chain in model.chains():
      for residue_group in chain.residue_groups():
        atom_groups  = residue_group.atom_groups()
        if (len(atom_groups) > 1) : # alt conf, skip
          continue
        for atom_group in atom_groups :
          # Check for non standard atoms in the residue
          # Or a label indicating the residue is a water
          elements = set(i.strip().upper()
                         for i in atom_group.atoms().extract_element())
          elements.difference_update(["H", "N", "C", "O", "S", "P", "SE"])
          resname = atom_group.resname.strip().upper()

          if (resname in WATER_RES_NAMES) :
            atoms = atom_group.atoms()
            if (len(atoms) == 1) : # otherwise it probably has hydrogens, skip
              waters.append(atoms[0].i_seq)
    print >> out, ""
    print >> out, "%d waters to analyze" % len(waters)
    if (len(waters) == 0) : return
    nproc = easy_mp.get_processes(self.nproc)
    ions = []
    if (nproc == 1) :
      print >> out, ""
      for water_i_seq in waters :
        water_props = self.analyze_water(
          i_seq=water_i_seq,
          debug = debug,
          candidates = candidates,
          show_only_map_outliers = show_only_map_outliers)
        if (water_props is not None) :
          water_props.show_summary(out=out, debug=debug)
          ions.append((water_i_seq, [water_props.final_choice]))
    else :
      print >> out, "Parallelizing across %d processes" % nproc
      print >> out, ""
      analyze_water = _analyze_water_wrapper(manager=self,
        debug=debug,
        candidates=candidates,
        show_only_map_outliers=show_only_map_outliers)
      results = easy_mp.pool_map(
        fixed_func=analyze_water,
        args=waters,
        processes=nproc)
      for water_i_seq, final_choice, result_str in results :
        if (result_str is not None) :
          print >> out, result_str
        if final_choice is not None:
          ions.append((water_i_seq, [final_choice]))
    return ions

class _analyze_water_wrapper (object) :
  """
  Simple wrapper for calling manager.analyze_water with keyword arguments
  in a parallelized loop.  Because the water_result object is not pickle-able,
  we only return the final ion choice and summary string.
  """
  def __init__ (self, manager, **kwds) :
    self.manager = manager
    self.kwds = dict(kwds)

  def __call__ (self, i_seq) :
    result = self.manager.analyze_water(i_seq, **(self.kwds))
    out = cStringIO.StringIO()
    if (result is not None) :
      result.show_summary(out=out,
        debug=self.kwds.get("debug", False))

    result_str = out.getvalue()
    if (result_str == "") : result_str = None
    return i_seq, getattr(result, "final_choice", None), result_str

class water_result (object) :
  """
  Container for storing the results of manager.analyze_water for later display
  and retrieval.
  """
  def __init__ (self,
      atom_props,
      filtered_candidates,
      matching_candidates,
      nuc_phosphate_site,
      looks_like_water,
      looks_like_halide,
      ambiguous_valence_cutoff,
      valence_used,
      final_choice) :
    adopt_init_args(self, locals())

  def show_summary (self, out=None, debug=False) :
    if (out is None) : out = sys.stdout
    results = self.matching_candidates
    if (len(results) > 0) :
      self.atom_props.show_properties(identity="HOH", out = out)
      if (self.nuc_phosphate_site) :
        print >> out, "  appears to be nucleotide coordination site"
      if (self.final_choice is not None) :
        # We have one result that we are reasonably certain of
        elem_params, score = results[0]
        self.atom_props.show_ion_results(
          identity=str(self.final_choice),
          out=out,
          valence_used=self.valence_used,
          confirmed=True)
        print >> out, ""
      elif (len(results) > 1) :
        # We have a couple possible identities for the atom
        below_cutoff = [ elem_params for elem_params, score in results
                        if score < self.ambiguous_valence_cutoff]
        if len(below_cutoff) == 1:
          elem_params = below_cutoff[0]
          print >> out, "  ambigous results, best valence from %s" % \
            str(elem_params)
          self.atom_props.show_ion_results(
            identity = str(elem_params),
            out=out,
            valence_used = True)
          print >> out, ""
        else:
          ions = [str(i[0]) for i in sorted(results, key=lambda x: x[1])]
          print >> out, "  ambiguous results, could be %s" % \
            (", ".join(ions))
          for elem_params, score in results :
            self.atom_props.show_ion_results(identity=str(elem_params), out=out)
          print >> out, ""
    else:
      if (not self.looks_like_water) or (self.nuc_phosphate_site) :
        self.atom_props.show_properties(identity="HOH", out = out)
        if (self.nuc_phosphate_site) :
          print >> out, "  appears to be nucleotide coordination site"
        # try anions now
        if (self.looks_like_halide) :
          print >> out, "  Probable element: %s" % str(self.final_choice)
          print >> out, ""
        else:
          # atom is definitely not water, but no reasonable candidates found
          # print out why all the metals we tried failed
          if (debug) :
            print >> out, "  insufficient data to identify atom"
            for params in self.filtered_candidates:
              if (self.atom_props.has_compatible_ligands(str(params))) :
                self.atom_props.show_ion_results(identity = str(params),
                  out = out)
            print >> out, ""

class AtomProperties (object):
  """
  Collect physical attributes of an atom, including B, occupancy, and map
  statistics, and track those which are at odds with its chemical identity.
  """

  LOW_B, HIGH_B, LOW_OCC, HIGH_OCC, NO_2FOFC_PEAK, FOFC_PEAK, FOFC_HOLE, \
    ANOM_PEAK, NO_ANOM_PEAK, BAD_GEOMETRY, NO_GEOMETRY, BAD_VECTORS, \
    BAD_VALENCES, TOO_FEW_NON_WATERS, TOO_FEW_COORD, TOO_MANY_COORD, \
    LIKE_COORD, BAD_COORD_ATOM, BAD_FPP, BAD_COORD_RESIDUE, VERY_BAD_VALENCES \
    = range(21)

  def __init__(self, atom, i_seq, manager):
    self.atom = atom
    self.i_seq = i_seq
    self.resname = atom.parent().resname

    # Grab all the atoms within 3.5 Angstroms
    nearby_atoms = manager.find_nearby_atoms(
      i_seq = i_seq,
      distance_cutoff = 3.5)

    self.nearby_atoms = [i for i in nearby_atoms
                         if not _same_atom_different_altloc(i[0], atom) and
                         _get_element(i[0]) not in ["H", "D"]]

    self.residue_counts = count_coordinating_residues(self.nearby_atoms)
    self.geometries = find_coordination_geometry(nearby_atoms)

    map_stats = manager.map_stats(i_seq)
    self.peak_2fofc = map_stats.two_fofc
    self.peak_fofc = map_stats.fofc
    self.peak_anom = map_stats.anom

    self.inaccuracies = {}
    self.ignored = {}
    self.vectors = {}
    self.valence_sum = {}
    self.vector_sum = {}
    self.score = {}
    self.fpp_expected = {}
    self.expected_params = {}
    self.bad_coords = {}

    # Determine the f'' value if possible using phaser
    self.fpp = manager.get_fpp(i_seq)
    self.fpp_ratios = {}

  def is_correctly_identified(self, identity = "HOH"):
    """
    Returns whether factors indicate that the atom was correctly identified.
    """
    if identity is None:
      identity = _identity(self.atom)

    return len(self.inaccuracies[identity]) == 0

  def has_compatible_ligands (self, identity) :
    """
    Indicates whether the coordinating atoms are of the allowed type (e.g.
    no N or S atoms coordinating CA, etc.)
    """
    return len(self.bad_coords[identity]) == 0

  # XXX obsolete, delete?
  def atom_weight (self):
    """
    Evaluates whether factors indicate that the atom is lighter, heavier, or
    isoelectric to what it is currently identified as.

    Returns -1 if lighter, 0 if isoelectronic, and 1 if heavier.
    """
    identity = _identity(self.atom)
    if (self.resname == "HOH") :
      # waters that don't have B-factors at least 1 stddev below the mean are
      # presumed to be correct
      if (self.atom.b > manager.b_mean_hoh - manager.b_stddev_hoh) :
        return 0
      identity = "HOH"

    if self.is_correctly_identified(identity):
      return 0

    # B-factors/occupancies?
    if self.FOFC_PEAK in self.inaccuracies[identity] or self.atom.b < 1:
      return 1

    if self.FOFC_HOLE in self.inaccuracies[identity]:
      return -1

    return 0

  def check_ion(self,
      ion_params,
      server,
      wavelength = None,
      require_valence=True):
    """
    Checks whether or not a given ion, satisfies the metal-coordination
    parameters, specified by ion_params, such as valence sum, geometry, etc.
    """

    identity = _identity(ion_params)
    inaccuracies = self.inaccuracies[identity] = set()
    self.expected_params[identity] = ion_params
    ignored = self.ignored[identity] = set()
    import iotbx.pdb
    get_class = iotbx.pdb.common_residue_names_get_class

    # if the atom is clearly not a water, optionally relax some rules.  this
    # will be more sensitive for transition metals, without finding a lot of
    # spurious Mg/Na sites.
    strict_rules = require_valence or self.is_correctly_identified("HOH")

    # Check for all non-overlapping atoms within 3 A of the metal
    coord_atoms = [pair for index, pair in enumerate(self.nearby_atoms)
                   if abs(pair[1]) < 3 and
                   all(abs(pair[1] - next_pair[1]) > 0.3
                       for next_pair in self.nearby_atoms[index + 1:])]

    if len(coord_atoms) < ion_params.coord_num_lower:
      inaccuracies.add(self.TOO_FEW_COORD)

    if len(coord_atoms) > ion_params.coord_num_upper:
      inaccuracies.add(self.TOO_MANY_COORD)

    # Coordinating atoms closer than 3.0 A are not positively charged
    n_non_water = 0
    self.bad_coords[identity] = []

    for other_atom, vector in self.nearby_atoms:
      other_name = other_atom.name.strip().upper()
      other_resname = other_atom.fetch_labels().resname.strip().upper()
      other_element = _get_element(other_atom)

      if (not other_resname in WATER_RES_NAMES) :
        n_non_water += 1
      else:
        # Everything can potentially be coordinated by water
        continue

      if abs(vector) < 3.0:
        # Anion coordinated with anion, or cation coordinating with cation
        # (Make sure atom charges are on the opposite sides of 0)
        other_charge = other_atom.charge_as_int()

        if other_charge is None:
          other_charge = get_charge(other_element)

        if (ion_params.allowed_coordinating_atoms and
            other_element not in ion_params.allowed_coordinating_atoms):
          # Check if atom is of an allowed element, if restricted
          self.bad_coords[identity].append((other_atom, vector))
          inaccuracies.add(self.BAD_COORD_ATOM)
        elif (ion_params.allowed_coordinating_residues and
              other_resname not in ion_params.allowed_coordinating_residues and
              get_class(other_resname) == "common_amino_acid" and
              other_name not in ["C", "O", "N", "CA"]) : # XXX probably just O
          # Check if atom is of an allowed residue type, if part of a sidechain
          self.bad_coords[identity].append((other_atom, vector))
          inaccuracies.add(self.BAD_COORD_RESIDUE)
        elif (cmp(0, other_atom.charge_as_int()) == cmp(0, ion_params.charge)) :
          # Check if coordinating atom is of opposite charge
          self.bad_coords[identity].append((other_atom, vector))
          inaccuracies.add(self.LIKE_COORD)
        elif (ion_params.charge > 0 and
            other_element in ["N"] and
            other_resname in ["LYS", "ARG", "ASN", "GLN"]):
          # Coordinating nitrogen most likely positive.
          #
          # Ignore nitrogens without a charge label that are on positively
          # charged amino acids.
          self.bad_coords[identity].append((other_atom, vector))
          inaccuracies.add(self.LIKE_COORD)

    # Check the number of coordinating waters
    if (n_non_water < ion_params.min_coordinating_non_waters) :
      inaccuracies.add(self.TOO_FEW_NON_WATERS)

    # Check the geometry of the coordinating atoms
    if (ion_params.allowed_geometries) and (strict_rules) :
      allowed = [i[0] in ion_params.allowed_geometries
                 for i in self.geometries]

      if not self.geometries:
        if strict_rules:
          inaccuracies.add(self.NO_GEOMETRY)
      elif not any(allowed):
        inaccuracies.add(self.BAD_GEOMETRY)

    # Check for reasonable vector/valence values
    vectors = server.calculate_valences(ion_params, self.nearby_atoms)
    self.vectors[identity] = vectors

    self.valence_sum[identity] = sum([abs(i) for i in vectors])
    self.vector_sum[identity] = abs(sum(vectors, col((0, 0, 0))))

    if self.vector_sum[identity] > ion_params.vec_sum_cutoff:
      if (strict_rules) :
        inaccuracies.add(self.BAD_VECTORS)
      else :
        ignored.add(self.BAD_VECTORS)

    if (self.valence_sum[identity] < ion_params.cvbs_expected * 0.75 or
        self.valence_sum[identity] > ion_params.cvbs_expected * 1.25):
      inaccuracies.add(self.VERY_BAD_VALENCES)
    else:
      if (self.valence_sum[identity] < ion_params.cvbs_lower or
          self.valence_sum[identity] > ion_params.cvbs_upper):
        if (strict_rules) :
          inaccuracies.add(self.BAD_VALENCES)
        else :
          ignored.add(self.BAD_VALENCES)

    self.score[identity] = abs(self.valence_sum[identity] -
                               ion_params.cvbs_expected)

    # Check the f'' values if available
    # XXX in theory the fpp_ratio should be no more than 1.0 unless we are
    # right on the peak wavelength.  in practice Phaser can overshoot a little
    # bit, so we need to be more tolerant.  picking the maximum f'' from the
    # Sasaki and Henke tables will also limit the ratio.
    fpp_ratio_max = 1.05
    fpp_ratio_min = 0.4
    if self.fpp is not None and wavelength is not None:
      fpp_expected_sasaki = sasaki.table(ion_params.element).at_angstrom(
        wavelength).fdp()
      fpp_expected_henke = henke.table(ion_params.element).at_angstrom(
        wavelength).fdp()
      self.fpp_expected[identity] = max(fpp_expected_sasaki,fpp_expected_henke)
      self.fpp_ratios[identity] = self.fpp / self.fpp_expected[identity]
      if ((self.fpp_ratios[identity] > fpp_ratio_max) or
          ((self.fpp >= 0.2) and (self.fpp_ratios[identity] < fpp_ratio_min))) :
        inaccuracies.add(self.BAD_FPP)

  def show_properties (self, identity, out=sys.stdout) :
    """
    Show atomic properties that are independent of the suspected identity.
    """
    print >> out, "%s:" % self.atom.id_str()
    b_flag = ""
    if (self.LOW_B in self.inaccuracies[identity]) :
      b_flag = " <<<"
    elif (self.HIGH_B in self.inaccuracies[identity]) :
      b_flag = " !!!"
    print >> out, "  B-factor:      %6.2f%s" % (self.atom.b, b_flag)
    occ_flag = ""
    if (self.LOW_OCC in self.inaccuracies[identity]) :
      occ_flag = " !!!"
    elif (self.HIGH_OCC in self.inaccuracies[identity]) :
      occ_flag = " <<<"
    print >> out, "  Occupancy:     %6.2f%s" % (self.atom.occ, occ_flag)
    twofofc_flag = ""
    if (self.NO_2FOFC_PEAK in self.inaccuracies[identity]) :
      twofofc_flag = " !!!"
    print >> out, "  2mFo-DFc map:  %6.2f%s" % (self.peak_2fofc, twofofc_flag)
    fofc_flag = ""
    if (self.FOFC_PEAK in self.inaccuracies[identity]) :
      fofc_flag = " <<<"
    elif (self.FOFC_HOLE in self.inaccuracies[identity]) :
      fofc_flag = " !!!"
    print >> out, "  mFo-DFc map:   %6.2f%s" % (self.peak_fofc, fofc_flag)
    if (self.peak_anom is not None) :
      anom_flag = ""
      if (self.ANOM_PEAK in self.inaccuracies[identity]) :
        anom_flag = " <<<"
      elif (self.NO_ANOM_PEAK in self.inaccuracies[identity]) :
        anom_flag = " !!!"
      print >> out, "  Anomalous map: %6.2f%s" % (self.peak_anom, anom_flag)
    if self.fpp is not None:
      fpp_flag = ""
      if (self.fpp >= 0.2) :
        fpp_flag = " <<<"
      print >> out, "  f'':           %6.2f%s" % (self.fpp, fpp_flag)
    if self.nearby_atoms is not None:
      print >> out, "  Nearby atoms:"
      angstrom = u"\u00C5".encode("utf-8", "strict").strip()
      for atom, vector in self.nearby_atoms :
        print >> out, "    %s (%5.3f %s)" % (atom.id_str(), abs(vector),
          angstrom)
      if self.geometries:
        print >> out, "  Coordinating geometry:"
        degree = u"\N{DEGREE SIGN}".encode("utf-8", "strict")
        for geometry, deviation in self.geometries:
          print >> out, "    %-15s (average deviation: %.3f%s)" % \
            (geometry, deviation, degree)

  def show_ion_results (self, identity = None, out = sys.stdout,
      confirmed = False, valence_used = True):
    """
    Show statistics for a proposed element identity.
    """

    if not identity:
      identity = _identity(self.atom)

    inaccuracies = self.inaccuracies[identity]
    ignored = self.ignored[identity]

    if identity != _identity(self.atom):
      if (confirmed) :
        print >> out, "  Probable element: %s" % identity
      else :
        print >> out, "  Atom as %s:" % identity
    else:
      print >> out, "    %s:" % self.atom.id_str()

    if identity in self.vector_sum and self.vector_sum[identity] is not None:
      problem = ((self.BAD_VECTORS in inaccuracies) or
                 (self.BAD_VECTORS in ignored))

      print >> out, "    Vector sum:  %6.3f" % self.vector_sum[identity],
      print >> out, " !!!" if problem else ""

    if identity in self.valence_sum and self.valence_sum[identity] is not None:
      problem = inaccuracies.union(ignored).intersection(
        [self.BAD_VALENCES, self.VERY_BAD_VALENCES])

      print >> out, "    Valence sum: %6.3f" % self.valence_sum[identity],
      if valence_used:
        print >> out, "(expected: %6.3f)" % \
          self.expected_params[identity].cvbs_expected,
      print >> out, " !!!" if problem else ""

    if self.NO_GEOMETRY in inaccuracies:
      print >> out, "    No distinct geometry !!!"

    if self.BAD_GEOMETRY in inaccuracies:
      print >> out, "    Unexpected geometry  !!!"

    bad_coord = [self.LIKE_COORD, self.BAD_COORD_ATOM, self.BAD_COORD_RESIDUE]
    if inaccuracies.intersection(bad_coord):
      print >> out, "    Bad coordinating atom%s:" % \
        ("s" if len(self.bad_coords[identity]) != 1 else "")
      angstrom = u"\u00C5".encode("utf-8", "strict").strip()
      for atom, vector in self.bad_coords[identity]:
        print >> out, "    %s (%5.3f %s) !!!" % (atom.id_str(), abs(vector),
          angstrom)

    if self.TOO_FEW_NON_WATERS in inaccuracies:
      print >> out, "    Too few coordinating waters !!!"
    if self.TOO_FEW_COORD in inaccuracies:
      print >> out, "    Too few coordinating atoms !!!"
    if self.TOO_MANY_COORD in inaccuracies:
      print >> out, "    Too many coordinating atoms !!!"

    if (self.fpp is not None) :
      print >> out, "    f'' ratio:   %6.3f" % (self.fpp_ratios[identity]),
      if self.BAD_FPP in inaccuracies:
        print >> out, " !!!"
      else :
        print >> out, ""

  # XXX can we get away with just one oxygen?
  def looks_like_nucleotide_phosphate_site (self,
      min_phosphate_oxygen_atoms=2,
      distance_cutoff=2.5) : # XXX wild guess
    """
    Decide whether the atom is coordinating phosphate oxygens from a
    nucleotide, based on common atom names.
    """
    n_phosphate_oxygens = 0
    for atom, vector in self.nearby_atoms :
      if (len(atom.name) < 3) or (_get_element(atom) not in ["O"]) :
        continue
      if ((atom.name[1:3] in ["O1","O2","O3"]) and
          (atom.name[3] in ["A","B","G"])) :
        if (abs(vector) <= distance_cutoff) :
          n_phosphate_oxygens += 1
    return (n_phosphate_oxygens == min_phosphate_oxygen_atoms)

def _get_element(atom):
  """
  Grabs an atom's element, stripping off whitespace and making sure the
  letters are capitalized.
  """
  return atom.element.strip().upper()

def _identity(atom):
  """
  Covers an atom into a string representing its element and charge.
  """
  charge = atom.charge

  if not isinstance(charge, int):
    charge = atom.charge_as_int()

  return "%s(%+d)" % (_get_element(atom), charge)

def find_anomalous_scatterers (
    fmodel,
    pdb_hierarchy,
    wavelength,
    prefix="phaser_AX",
    verbose=True,
    log=None) :
  """
  Runs phaser to collect a list of anomalous scatterers around the unit cell.
  """
  if (log is None) : log = sys.stdout
  if (wavelength is None) :
    raise Sorry("Please specify the X-ray wavelength.")
  if (not libtbx.env.has_module("phaser")) :
    print "Phaser not available"
    return None
  from phaser.phenix_adaptors import sad_target
  import phaser
  assert (fmodel.f_obs().anomalous_flag())
  adaptor = sad_target.data_adaptor(
    f_obs=fmodel.f_obs(),
    r_free_flags=fmodel.r_free_flags(),
    verbose=True)
  phaser_input = adaptor.set_ep_auto_data(
    xray_structure=fmodel.xray_structure,
    pdb_hierarchy=pdb_hierarchy)
  phaser_input.setLLGC_COMP(True)
  phaser_input.addLLGC_SCAT("AX")
  phaser_input.setWAVE(wavelength)
  phaser_input.setROOT(prefix)
  phaser_input.setMUTE(not verbose)
  phaser_input.setHKLO(False)
  phaser_input.setXYZO(False)
  out = phaser.Output()
  if (not verbose) :
    out.setPackagePhenix(null_out())
  else :
    out.setPackagePhenix(log)
  result = phaser.runEP_AUTO(phaser_input, out)
  return result.getIotbx() # pdb hierarchy object

def create_manager (
    pdb_hierarchy,
    geometry_restraints_manager,
    fmodel,
    wavelength,
    params,
    resolution_factor=0.25,
    nproc=Auto,
    verbose=False,
    log=None) :
  connectivity = \
    geometry_restraints_manager.shell_sym_tables[0].full_simple_connectivity()
  manager = Manager(
    fmodel = fmodel,
    pdb_hierarchy = pdb_hierarchy,
    xray_structure = fmodel.xray_structure,
    connectivity = connectivity,
    wavelength = wavelength,
    params = params,
    nproc = nproc,
    verbose = verbose,
    log = log)
  return manager

def _same_atom_different_altloc(atom1, atom2):
  """
  Determines whether atom1 and atom2 differ only by their alternate location.
  """

  label1, label2 = [i.fetch_labels() for i in [atom1, atom2]]
  name1, name2 = atom1.name.strip(), atom2.name.strip()
  chain1, chain2 = label1.chain_id, label2.chain_id
  res1, res2 = label1.resid(), label2.resid()

  return name1 == name2 and chain1 == chain2 and res1 == res2

def count_coordinating_residues (nearby_atoms, distance_cutoff=3.0) :
  """
  Count the number of residues of each type involved in the coordination
  sphere.  This may yield additional clues to the identity of ions, e.g. only
  Zn will have 4 Cys residues.
  """
  unique_residues = []
  residue_counts = {}
  for atom, vector in nearby_atoms :
    if (abs(vector) <= distance_cutoff) :
      parent = atom.parent()
      for residue in unique_residues :
        if (residue == parent) :
          break
      else :
        resname = parent.resname
        if (not resname in residue_counts) :
          residue_counts[resname] = 0
        residue_counts[resname] += 1
        unique_residues.append(parent)
  return residue_counts
