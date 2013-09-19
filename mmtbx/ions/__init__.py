# -*- coding: utf-8; py-indent-offset: 2 -*-
"""
Examines a structure for metal ions. Can iterate over all atoms, examining
their density and chemical environment to determine if they are correctly
identified, or if there are better candidate ions that they can be replaced
with.

See mmtbx.ions.build to actually modify the structure, code in this module
only prints out messages to the log.
"""

from __future__ import division
from mmtbx.ions.parameters import server, MetalParameters
from mmtbx.ions.geometry import find_coordination_geometry
from cctbx import crystal, adptbx
from cctbx.eltbx import sasaki, henke
from scitbx.matrix import col, distance_from_plane
from scitbx.array_family import flex
from libtbx import group_args, adopt_init_args, Auto, slots_getstate_setstate
from libtbx.str_utils import make_sub_header, format_value, framed_output
from libtbx.utils import null_out, Sorry
from libtbx import easy_mp
import libtbx.phil
import libtbx.load_env
from math import sqrt
import cStringIO
import time
import sys

chloride_params_str = """
max_distance_to_amide_n = 3.5
  .type = float
  .input_size = 80
  .help = Max distance to amide N atom
max_distance_to_cation = 3.5
  .type = float
  .input_size = 80
min_distance_to_anion = 3.5
  .type = float
  .input_size = 80
min_distance_to_other_sites = 2.1
  .type = float
  .input_size = 80
max_distance_to_hydroxyl = 3.5
  .type = float
delta_amide_h_angle = 20
  .type = float
  .input_size = 80
  .help = Allowed deviation (in degrees) from ideal N-H-X angle
delta_planar_angle = 10
  .type = float
  .input_size = 80
max_deviation_from_plane = 0.8
  .type = float
radius = 2.0
  .type = float
  .input_size = 80
  .short_caption = Sampling radius
"""

ion_master_phil = libtbx.phil.parse("""
require_valence = False
  .type = bool
  .help = Toggles the use of valence calculations
  .short_caption = Require good bond valence
ambiguous_valence_cutoff = 0.5
  .type = float
  .help = "Cutoff to uniquely select one of many metals by comparing its observed and expected valence"
d_min_strict_valence = 1.5
  .type = float
  .short_caption = Resolution limit for strict valence rules
anom_map_type = *residual simple llg
  .type = choice
  .help = Type of anomalous difference map to use.  Default is a residual \
    map, showing only unmodeled scattering.
find_anomalous_substructure = Auto
  .type = bool
use_phaser = True
  .type = bool
  .help = Toggles the use of Phaser for calculation of f-prime values.
  .short_caption = Use Phaser to calculate f-double-prime
aggressive = False
  .type = bool
  .help = Toggles more permissive settings for flagging waters as heavier \
    elements.  Not recommended for automated use.
water
  .short_caption = Water filtering
  .style = box auto_align
{
  min_2fofc_level = 1.8
    .type = float
    .input_size = 80
    .help = Minimum water 2mFo-DFc map value.  Waters below this cutoff will \
      not be analyzed further.
    .short_caption = Min. allowed 2mFo-DFc map value
  max_fofc_level = 3.0
    .type = float
    .input_size = 80
    .help = Maximum water mFo-DFc map value
    .short_caption = Max. expected mFo-DFc map value
  max_anom_level = 3.0
    .type = float
    .input_size = 80
    .help = Maximum water anomalous map value
    .short_caption = Max. expected anomalous map value
  max_occ = 1.0
    .type = float
    .input_size = 80
    .help = Maximum water occupancy
    .short_caption = Max. expected occupancy
  min_b_iso = 1.0
    .type = float
    .input_size = 80
    .help = Minimum water isotropic B-factor
    .short_caption = Min. expected B-factor
  fp_max = 1.0
    .type = float
  fpp_max = 0
    .type = float
  max_stddev_b_iso = 5
    .type = float
    .input_size = 80
    .short_caption = Max. standard deviations below mean B-iso
  min_frac_b_iso = 0.2
    .type = float
    .input_size = 80
    .short_caption = Min. fraction of mean B-iso
  min_frac_calpha_b_iso = 0.75
    .type = float
  max_frac_calpha_2fofc = 1.2
    .type = float
  min_2fofc_coordinating = 0.9
    .type = float
    .help = Minimum 2mFo-DFc for waters involved in coordination shells.
}
chloride
  .short_caption = Halide ions
  .style = box auto_align
{
  %s
}
phaser
  .short_caption = Phaser options
  .caption = These parameters control the identification of anomalous \
    scatterers using Phaser substructure completion.
  .style = box auto_align
{
  llgc_ncycles = None
    .type = int
    .short_caption = Number of cycles
  distance_cutoff = 1.5
    .type = float
    .input_size = 80
    .short_caption = Min. required separation from other atoms
  distance_cutoff_same_site = 0.7
    .type = float
    .input_size = 80
    .short_caption = Max. separation from mapped atom
  fpp_ratio_min = 0.2
    .type = float
    .input_size = 80
    .help = Minimum ratio of refined/theoretical f-double-prime.
    .short_caption = Min. f-double-prime ratio
  fpp_ratio_max = 1.1
    .type = float
    .input_size = 80
    .help = Maximum ratio of refined/theoretical f-double-prime.
    .short_caption = Max. f-double-prime ratio
}
""" % (chloride_params_str))

from iotbx.pdb import common_residue_names_water as WATER_RES_NAMES
DEFAULT_IONS = ["MG", "CA", "ZN", "CL"]
HALIDES = ["F", "CL", "BR", "IOD"]
TRANSITION_METALS = ["MN", "FE", "CO", "CU", "NI", "ZN", "PT"]
NUC_PHOSPHATE_BINDING = ["MG", "CA", "MN"]
SUPPORTED = TRANSITION_METALS + HALIDES + ["NA", "MG", "K", "CA", "CD", "HG"]
WATER, WATER_POOR, LIGHT_ION, HEAVY_ION = range(4)

def check_supported (elements) :
  if (elements is None) :
    raise Sorry("No elements specified for ion picking - must be either "+
      "'Auto' or a comma-separated list of element symbols.")
  elif (elements is not Auto) :
    # XXX somehow comma-separation of phil strings fields doesn't work
    if (isinstance(elements, list)) and (len(elements) == 1) :
      elements = elements[0].split(",")
    if (elements == ['X']) : # XXX hack for testing - X is "dummy" element
      return True
    for elem in elements :
      if (not elem.strip().upper() in SUPPORTED) :
        raise Sorry("Element '%s' not supported!  Choices are: %s" %
                    (elem, " ".join(SUPPORTED)))
  return True

class atom_contact (slots_getstate_setstate) :
  """
  Container for information about an interacting atom.  Most of the methods
  are simply wrappers for frequently called operations on the atom object, but
  symmetry-aware.
  """
  __slots__ = ["atom", "vector", "site_cart", "rt_mx", "server"]
  def __init__ (self, atom, vector, site_cart, rt_mx, server) :
    self.atom = atom
    self.vector = vector
    self.site_cart = site_cart
    self.rt_mx = rt_mx
    self.server = server

  def distance (self) :
    """Actual distance from target atom"""
    return abs(self.vector)

  def distance_from (self, other) :
    """Distance from another contact"""
    return abs(self.vector - other.vector)

  def id_str (self, suppress_rt_mx = False) :
    if (not self.rt_mx.is_unit_mx()) and (not suppress_rt_mx) :
      return self.atom.id_str() + " " + str(self.rt_mx)
    else :
      return self.atom.id_str()

  def atom_name (self) :
    return self.atom.name.strip()

  def resname (self) :
    return self.atom.fetch_labels().resname.strip().upper()

  @property
  def element (self) :
    return self.server.get_element(self.atom)

  @property
  def charge (self) :
    return self.server.get_charge(self.atom)

  @property
  def occ (self) :
    return self.atom.occ

  def atom_i_seq (self) :
    return self.atom.i_seq

  def altloc (self) :
    return self.atom.fetch_labels().altloc.strip()

  def atom_id_no_altloc (self, suppress_rt_mx = False) :
    """Unique identifier for an atom, ignoring the altloc but taking the
    symmetry operator (if any) into account."""
    labels = self.atom.fetch_labels()
    base_id = labels.chain_id + labels.resid() + self.atom.name
    if (self.rt_mx.is_unit_mx()) or (suppress_rt_mx) :
      return base_id
    else :
      return base_id + " " + str(self.rt_mx)

  def __eq__ (self, other) :
    """Equality operator, taking symmetry into account but ignoring the
    altloc identifier."""
    return (other.atom_id_no_altloc() == self.atom_id_no_altloc())

  def __abs__ (self) :
    return self.distance()

  def __str__ (self) :
    return self.id_str()

# - Phaser Anomalous signal
class Manager (object) :
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
    if (log is None) : log = null_out()
    self.fmodel = fmodel
    self.atoms_to_props = {}
    self.params = params
    self.wavelength = wavelength
    self.server = server()
    self.nproc = nproc
    self.phaser_substructure = None
    self.flag_refine_substructure = False
    self.fpp_from_phaser_ax_sites = None
    self.site_fp = None
    self.site_fdp = None
    self._map_values = {}
    self.update_structure(
      pdb_hierarchy = pdb_hierarchy,
      xray_structure = xray_structure,
      connectivity = connectivity,
      log = log)
    self.update_maps()
    # The default behavior is to refine the anomalous structure when possible;
    # this can use either the CCTBX anomalous group refinement or Phaser's
    # substructure completion (but not both).  Some additional machinery is
    # required to handle this distinction.
    refine_substructure = \
      params.find_anomalous_substructure if params is not None else False
    if (refine_substructure is Auto) :
      refine_substructure = fmodel.f_obs().anomalous_flag() and \
                            (wavelength is not None)
    if (refine_substructure) :
      if (wavelength is None) :
        raise Sorry("Wavelength required when "+
                    "find_anomalous_substructure=True.")
      elif (not fmodel.f_obs().anomalous_flag()) :
        raise Sorry("Anomalous data required when "+
                    "find_anomalous_substructure=True.")
      if ((params.use_phaser) and (libtbx.env.has_module("phaser"))) :
        print >> log, "  Running Phaser substructure completion..."
        t1 = time.time()
        #open("tmp1.pdb", "w").write(pdb_hierarchy.as_pdb_string(xray_structure))
        phaser_result = find_anomalous_scatterers(
          fmodel = fmodel,
          pdb_hierarchy = pdb_hierarchy,
          wavelength = wavelength,
          verbose = verbose,
          log = log,
          n_cycles = params.phaser.llgc_ncycles)
        t2 = time.time()
        print >> log, "    time: %.1fs" % (t2-t1)
        if (phaser_result is None) :
          print >> log, "  ERROR: Phaser substructure completion failed!"
        else :
          self.phaser_substructure = phaser_result.atoms()
          if (len(self.phaser_substructure) == 0) :
            print >> log, "  No anomalous scatterers found!"
          else :
            print >> log, "  %d anomalous scatterers found" % \
              len(self.phaser_substructure)
            self.analyze_substructure(log = log, verbose = True)
      else :
        self.flag_refine_substructure = True
        self.refine_anomalous_substructure(log = log)

  def refine_anomalous_substructure (self, log) :
    """
    Run simple "substructure completion" implemented using CCTBX tools (the
    command-line equivalent is mmtbx.refine_anomalous_substructure).  A faster
    alternative to Phaser, but not clear how effective this is at the moment.
    """
    from mmtbx.refinement import anomalous_scatterer_groups
    fmodel_tmp = self.fmodel.deep_copy()
    # XXX should we only be refining f''?
    anom_groups = anomalous_scatterer_groups.refine_anomalous_substructure(
      fmodel = fmodel_tmp,
      pdb_hierarchy = self.pdb_hierarchy,
      wavelength = self.wavelength,
      reset_water_u_iso = True,
      verbose = True,
      use_all_anomalous = True,
      out = log)
    scatterers = fmodel_tmp.xray_structure.scatterers()
    self.use_fdp = flex.bool(scatterers.size(), False)
    self.site_fp = flex.double(scatterers.size(), 0)
    self.site_fdp = flex.double(scatterers.size(), 0)
    for group in anom_groups :
      for i_seq in group.iselection :
        sc = scatterers[i_seq]
        self.use_fdp[i_seq] = True
        self.site_fp[i_seq] = sc.fp
        self.site_fdp[i_seq] = sc.fdp

  def update_structure (self, pdb_hierarchy, xray_structure,
      connectivity = None, log = None, refine_if_necessary = True) :
    """
    Set the current atomic data: PDB hierarchy, Xray structure, and simple
    connectivity list.
    """
    self.pdb_hierarchy = pdb_hierarchy
    self.xray_structure = xray_structure
    if self.fmodel:
      self.fmodel.update_xray_structure(xray_structure, update_f_calc = True)
    self.sites_frac = self.xray_structure.sites_frac()
    self.sites_cart = self.xray_structure.sites_cart()
    self.connectivity = connectivity
    self.pdb_atoms = pdb_hierarchy.atoms()
    self.unit_cell = xray_structure.unit_cell()
    self.use_fdp = flex.bool(xray_structure.scatterers().size(), False)
    self.ax_chain = None
    self._pair_asu_cache = {}
    self.u_iso_all = xray_structure.extract_u_iso_or_u_equiv()
    # Extract some information about the structure, including B-factor
    # statistics for non-HD atoms and waters
    sctr_keys = xray_structure.scattering_type_registry().type_count_dict()\
      .keys()
    self.hd_present = ("H" in sctr_keys) or ("D" in sctr_keys)
    sel_cache = pdb_hierarchy.atom_selection_cache()
    self.carbon_sel = sel_cache.selection("element C").iselection()
    self.calpha_sel = sel_cache.selection("name CA and element C").iselection()
    assert (len(self.carbon_sel) > 0)
    water_sel = sel_cache.selection("({}) and element O".format(
      " or ".join("resname " + i for i in WATER_RES_NAMES))).iselection()
    self.n_waters = len(water_sel)

    not_hd_sel = sel_cache.selection("not (element H or element D)").iselection()
    self.n_heavy = len(not_hd_sel)
    u_iso_all_tmp = self.u_iso_all.select(not_hd_sel)
    self.b_mean_all = adptbx.u_as_b(flex.mean(u_iso_all_tmp))
    self.b_stddev_all = adptbx.u_as_b(
       u_iso_all_tmp.standard_deviation_of_the_sample())
    self.b_mean_calpha = 0
    if (len(self.calpha_sel) > 0) :
      self.b_mean_calpha = adptbx.u_as_b(flex.mean(self.u_iso_all.select(
        self.calpha_sel)))
    self.b_mean_hoh = self.b_stddev_hoh = None
    if (self.n_waters > 0) :
      u_iso_hoh = self.u_iso_all.select(water_sel)
      self.b_mean_hoh = adptbx.u_as_b(flex.mean(u_iso_hoh))
      self.b_stddev_hoh = adptbx.u_as_b(
        u_iso_hoh.standard_deviation_of_the_sample())
    if (self.phaser_substructure is not None) :
      self.analyze_substructure(log = log)
    elif (self.flag_refine_substructure) and (refine_if_necessary) :
      self.refine_anomalous_substructure(log = log)

  def get_initial_b_iso (self) :
    if (getattr(self, "b_mean_hoh", None) is not None) :
      return self.b_mean_hoh
    else :
      return self.b_mean_all

  def get_map (self, map_type) :
    map_coeffs = self.fmodel.map_coefficients(
      map_type = map_type,
      exclude_free_r_reflections = True,
      fill_missing = True,
      pdb_hierarchy = self.pdb_hierarchy)
    if (map_coeffs is None) :
      return None
    return map_coeffs.fft_map(resolution_factor = 0.25,
      ).apply_sigma_scaling().real_map_unpadded()

  def update_maps (self) :
    """
    Generate new maps for the current structure, including anomalous map if
    data are anomalous.

    If Phaser is installed and anomalous data are available, the anomalous
    log-likelihood gradient (LLG) map can be used instead of the conventional
    anomalous difference map. This is only really useful if the anomalous
    scattering of existing atoms is modeled (and ideally, refined).
    """
    if self.fmodel is None:
      return
    def fft_map (map_coeffs, resolution_factor = 0.25):
      return map_coeffs.fft_map(resolution_factor = resolution_factor,
        ).apply_sigma_scaling().real_map_unpadded()
    map_types = ["2mFo-DFc", "mFo-DFc"]
    map_keys = ["2mFo-DFc", "mFo-DFc"]
    if (self.fmodel.f_obs().anomalous_flag()) :
      if (self.params.anom_map_type == "phaser") :
        map_types.append("llg")
      elif (self.params.anom_map_type == "residual") :
        map_types.append("anom_residual")
      else :
        map_types.append("anom")
      map_keys.append("anom")
    # To save memory, we sample atomic positions immediately and throw out
    # the actual maps (instead of keeping up to 3 in memory)
    sites_frac = self.xray_structure.sites_frac()
    sites_cart = self.xray_structure.sites_cart()
    self._principal_axes_of_inertia = [ None ] * len(sites_frac)
    self._map_variances = [ None ] * len(sites_frac)
    self.calpha_mean_two_fofc = 0
    for map_type, map_key in zip(map_types, map_keys) :
      real_map = self.get_map(map_type)
      if (real_map is not None) :
        # Gather values for map peaks at each site
        self._map_values[map_key] = flex.double(sites_frac.size(), 0)
        for i_seq, site_frac in enumerate(sites_frac) :
          atom = self.pdb_atoms[i_seq]
          resname = atom.fetch_labels().resname.strip().upper()
          if (resname in WATER_RES_NAMES + SUPPORTED or
              atom.segid.strip().upper() in ["ION"]):
            value = real_map.eight_point_interpolation(site_frac)
            self._map_values[map_key][i_seq] = value

        if map_type in ["2mFo-DFc"]:
          # Gather values on map variance and principal axes of interia
          from cctbx import maptbx
          for i_seq, site_cart in enumerate(sites_cart) :
            resname = self.pdb_atoms[i_seq].fetch_labels().resname.strip()
            if resname in WATER_RES_NAMES + SUPPORTED:
              # XXX not totally confident about how I'm weighting this...
              p_a_i = maptbx.principal_axes_of_inertia(
                real_map = real_map,
                site_cart = site_cart,
                unit_cell = self.unit_cell,
                radius = self.params.chloride.radius)
              self._principal_axes_of_inertia[i_seq] = p_a_i
              variance = maptbx.spherical_variance_around_point(
                real_map = real_map,
                unit_cell = self.unit_cell,
                site_cart = site_cart,
                radius = self.params.chloride.radius)
              self._map_variances[i_seq] = variance
            elif (i_seq in self.calpha_sel) :
              # Also collect some info in average C_alpha 2FoFc peak heights
              self.calpha_mean_two_fofc += real_map.eight_point_interpolation(
                sites_frac[i_seq])
        del real_map

    if (self.calpha_mean_two_fofc > 0) :
      n_calpha = len(self.calpha_sel)
      assert (n_calpha > 0)
      self.calpha_mean_two_fofc /= n_calpha

    # Gather info on carbons' average Fo peak height for use in estimating other
    # sites' atomic weight
    self.carbon_fo_values = None
    if (len(self.carbon_sel) > 0) :
      self.carbon_fo_values = flex.double()
      self._map_values["mFo"] = flex.double(sites_frac.size(), 0)
      fo_map = fft_map(self.fmodel.map_coefficients(
        map_type = "mFo",
        exclude_free_r_reflections = True,
        fill_missing = True))

      for i_seq, site_frac in enumerate(sites_frac) :
        resname = self.pdb_atoms[i_seq].fetch_labels().resname.strip()
        element = self.pdb_atoms[i_seq].element.strip()
        if (element == "C") or ((element == "O") and (resname in WATER_RES_NAMES)):
          map_value = fo_map.eight_point_interpolation(site_frac)
          self._map_values["mFo"][i_seq] = map_value
          if (element == "C") :
            self.carbon_fo_values.append(map_value)
      del fo_map

  def show_current_scattering_statistics (self, out = sys.stdout) :
    print >> out, ""
    print >> out, "Model and map statistics:"
    print >> out, "  mean mFo map height @ carbon: %s" % format_value("%.2f",
      flex.max(self.carbon_fo_values))
    if (self.calpha_mean_two_fofc > 0) :
      print >> out, "  mean 2mFo-DFc map height @ C-alpha: %s" % format_value(
        "%.2f", self.calpha_mean_two_fofc)
    print >> out, "  mean B-factor: %s" % format_value("%.2f", self.b_mean_all)
    if (self.b_mean_calpha > 0) :
      print >> out, "  mean C-alpha B-factor: %s" % format_value("%.2f",
        self.b_mean_calpha)
    print >> out, "  mean water B-factor: %s" % format_value("%.2f",
      self.b_mean_hoh)
    print >> out, ""

  def get_strict_valence_flag (self) :
    d_min = self.fmodel.f_obs().d_min()
    return (d_min < self.params.d_min_strict_valence)

  def find_nearby_atoms (self, i_seq, far_distance_cutoff = 3.0,
      near_distance_cutoff = 1.5, filter_by_bonding = True,
      filter_by_two_fofc = True):
    """
    Given site in the structure, return a list of nearby atoms with the
    supplied cutoff, and the vectors between them and the atom's site. Takes
    into account symmetry operations when finding nearby sites.
    """

    assert (i_seq < len(self.sites_frac))

    # Use pair_asu_table to find atoms within distance_cutoff of one another,
    # taking into account potential cell symmetry.
    # XXX should probably move this up, but it seems relatively fast
    pair_asu_table, asu_mappings, asu_table = self._pair_asu_cache.get(
      far_distance_cutoff, (None,None,None))
    if (pair_asu_table is None) :
      pair_asu_table = self.xray_structure.pair_asu_table(
        distance_cutoff = far_distance_cutoff)
      asu_mappings = pair_asu_table.asu_mappings()
      asu_table = pair_asu_table.table()
      self._pair_asu_cache[far_distance_cutoff] = (pair_asu_table, asu_mappings,
        asu_table)
    asu_dict = asu_table[i_seq]
    contacts = []
    site_i = self.sites_frac[i_seq]
    site_i_cart = self.sites_cart[i_seq]
    rt_mx_i_inv = asu_mappings.get_rt_mx(i_seq, 0).inverse()
    atom_i = self.pdb_atoms[i_seq]

    # Create the primary list of contacts
    for j_seq, j_sym_groups in asu_dict.items():
      site_j = self.sites_frac[j_seq]
      atom_j = self.pdb_atoms[j_seq]

      # Filter out hydrogens
      if atom_j.element.upper().strip() in ["H", "D"]:
        continue

      # Filter out alternate conformations of this atom
      if same_atom_different_altloc(atom_i, atom_j):
        continue

      # Gather up contacts with all symmetric copies
      for j_sym_group in j_sym_groups:
        rt_mx = rt_mx_i_inv.multiply(asu_mappings.get_rt_mx(j_seq,
          j_sym_group[0]))
        site_ji = rt_mx * site_j
        site_ji_cart = self.unit_cell.orthogonalize(site_ji)
        vec_i = col(site_i_cart)
        vec_ji = col(site_ji_cart)
        assert abs(vec_i - vec_ji) < far_distance_cutoff + 0.5
        contact = atom_contact(
          atom = atom_j,
          vector = vec_i - vec_ji,
          site_cart = site_ji_cart,
          rt_mx = rt_mx,
          server = self.server)
        # XXX I have no idea why the built-in handling of special positions
        # doesn't catch this for us
        if (j_seq == i_seq) and (not rt_mx.is_unit_mx()) :
          continue
        if contact.distance() < near_distance_cutoff:
          continue
        contacts.append(contact)

    # Filter out carbons that are judged to be "contacts", but are actually
    # just bonded to genuine coordinating atoms.  This is basically just a
    # way to handle sidechains such as His, Asp/Glu, or Cys where the carbons
    # may be relatively close to the metal site.
    if filter_by_bonding and self.connectivity:
      filtered = []
      all_i_seqs = [contact.atom_i_seq() for contact in contacts]

      for contact in contacts:
        # oxygen is always allowed, other elements may not be
        if contact.element not in ["C", "P", "S", "N"]:
          filtered.append(contact)
          continue

        # Remove atoms within 1.9 A contact distance
        if any(other_contact != contact and
               contact.distance_from(other_contact) < 1.9 and
               contact.distance() > other_contact.distance()
               for other_contact in contacts):
          continue

        # Examine the mode connectivity to catch more closely bonded atoms
        bonded_j_seqs = []
        for j_seq in self.connectivity[contact.atom_i_seq()]:
          if (j_seq in all_i_seqs) :
            bonded_j_seqs.append(j_seq)

        for j_seq in bonded_j_seqs:
          other_contact = contacts[all_i_seqs.index(j_seq)]

          if other_contact.element in ["N", "O", "S"] and \
            abs(other_contact) < abs(contact):
            break
        else:
          filtered.append(contact)

      contacts = filtered

    # Discard waters in poor density
    if (filter_by_two_fofc) :
      filtered = []
      for contact in contacts :
        if (contact.resname() in WATER_RES_NAMES) :
          two_fofc = self._map_values["2mFo-DFc"][contact.atom_i_seq()]
          if two_fofc < self.params.water.min_2fofc_coordinating:
            continue
        filtered.append(contact)
      contacts = filtered

    return contacts

  def guess_b_iso_real (self, i_seq) :
    """
    Guess an approximate B_iso for an atom by averaging the values for
    coordinating atoms.  This should partially compensate for the deflation of
    the B-factor due to incorrect scattering type, and consequently make the
    occupancy refinement more accurate.
    """
    contacts = self.find_nearby_atoms(i_seq, far_distance_cutoff = 3.0)
    if (len(contacts) == 0) :
      return adptbx.u_as_b(self.u_iso_all[i_seq])
    u_iso_sum = 0
    for contact in contacts :
      u_iso_sum += self.u_iso_all[contact.atom_i_seq()]
    return adptbx.u_as_b(u_iso_sum / len(contacts))

  def principal_axes_of_inertia (self, i_seq) :
    """
    Extracts the map grid points around a site, and calculates the axes of
    inertia (using the density values as weights). This is used to calculate
    the sphericity of the blob of density around a site.
    """
    return self._principal_axes_of_inertia[i_seq]

  def get_map_sphere_variance (self, i_seq) :
    """
    Calculate the density levels for points on a sphere around a given point.
    This will give us some indication whether the density falls off in the
    expected manner around a single atom.
    """
    return self._map_variances[i_seq]

  def get_b_iso (self, i_seq) :
    sc = self.xray_structure.scatterers()[i_seq]
    return adptbx.u_as_b(sc.u_iso_or_equiv(self.unit_cell))

  def map_stats (self, i_seq) :
    """
    Given a site in the structure, find the signal of the 2FoFc, FoFc, and
    anomalous map (When available).
    """
    value_2fofc = self._map_values["2mFo-DFc"][i_seq]
    value_fofc = self._map_values["mFo-DFc"][i_seq]
    value_anom = None
    if ("anom" in self._map_values) :
      value_anom = self._map_values["anom"][i_seq]
    return group_args(
      two_fofc = value_2fofc,
      fofc = value_fofc,
      anom = value_anom)

  def guess_molecular_weight (self, i_seq) :
    """
    Guesses the molecular weight of a site by scaling the signal from
    the mFo map by the average signal from carbon atoms.
    """
    map_values = self._map_values.get("mFo", None)
    if (map_values is None) : return None
    height = map_values[i_seq]
    mean_carbon = flex.mean(self.carbon_fo_values)
    assert (mean_carbon > 0)
    return 6 * height / mean_carbon

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

    site_frac = self.unit_cell.fractionalize(site_cart = site_cart)
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
      if (pair.i_seq < n_xray and pair.j_seq < n_xray) or \
        (pair.i_seq >= n_xray and pair.j_seq >= n_xray):
        continue

      rt_mx_i = asu_mappings.get_rt_mx_i(pair)
      rt_mx_j = asu_mappings.get_rt_mx_j(pair)
      rt_mx_ji = rt_mx_i.inverse().multiply(rt_mx_j)
      new_site_frac = rt_mx_ji * site_frac
      new_site_cart = self.unit_cell.orthogonalize(site_frac = new_site_frac)
      atom_seq = pair.i_seq if pair.i_seq < n_xray else pair.j_seq

      site_info = group_args(
        i_seq = atom_seq,
        rt_mx = rt_mx_ji,
        new_site_cart = new_site_cart,
        distance = sqrt(pair.dist_sq))

      if site_info.distance < distance_cutoff_same_site:
        same_atoms.append(site_info)
      else:
        other_atoms.append(site_info)

    same_atoms.sort(key = lambda x: x.distance)
    other_atoms.sort(key = lambda x: x.distance)

    return same_atoms, other_atoms

  def analyze_substructure (self, log = None, verbose = True) :
    """
    Given a list of AX pseudo-atoms placed by Phaser, finds the nearest real
    atoms, and if possible determines their equivalence.
    """
    def ax_atom_id (atom) :
      return "AX %s (fpp=%.3f)" % (atom.serial.strip(), atom.occ)

    assert (self.phaser_substructure is not None)
    self.fpp_from_phaser_ax_sites = flex.double(self.pdb_atoms.size(), -1)
    if (log is None) or (not verbose) : log = null_out()

    def _filter_site_infos(site_infos):
      # Organize a dictionary, keyed with each site's atom i_seq
      from collections import OrderedDict
      i_seqs = OrderedDict()
      for site_info in site_infos:
        if site_info.i_seq in i_seqs:
          i_seqs[site_info.i_seq].append(site_info)
        else:
          i_seqs[site_info.i_seq] = [site_info]

      # If we are picking from multiple of the same i_seq, select the closest
      for val in i_seqs.values():
        # Prefer unit operators
        if any(i.rt_mx.is_unit_mx() for i in val):
          for site_info in val:
            if not site_info.rt_mx.is_unit_mx():
              val.remove(site_info)
        val.sort(key = lambda x: x.distance)

      return [i[0] for i in i_seqs.values()]

    make_sub_header("Analyzing Phaser anomalous substructure", out = log)
    for atom in self.phaser_substructure :
      print >> log, ax_atom_id(atom)
      same_atoms, other_atoms = self.find_atoms_near_site(
        atom.xyz,
        distance_cutoff = self.params.phaser.distance_cutoff,
        distance_cutoff_same_site = self.params.phaser.distance_cutoff_same_site)

      same_atoms, other_atoms = \
        _filter_site_infos(same_atoms), _filter_site_infos(other_atoms)

      if len(same_atoms) == 0:
        print >> log, "  No match for %s" % ax_atom_id(atom)
        for other_atom in other_atoms:
          print >> log, "    %s (distance = %.3f)" % \
                (self.pdb_atoms[other_atom.i_seq].id_str(), other_atom.distance)
      elif len(same_atoms) == 1:
        print >> log, "  %s maps to %s" % \
          (ax_atom_id(atom), self.pdb_atoms[same_atoms[0].i_seq].id_str())
        self.use_fdp[same_atoms[0].i_seq] = True
        self.fpp_from_phaser_ax_sites[same_atoms[0].i_seq] = atom.occ
      else :
        print >> log, "  ambiguous results for %s:" % ax_atom_id(atom)
        for same_atom in same_atoms:
          print >> log, "    %s" % self.pdb_atoms[same_atom.i_seq].id_str()
    print >> log, ""

  def get_fpp (self, i_seq) :
    """
    Retrieve the refined f'' for a site.  Because this can come from either the
    built-in anomalous refinement or Phaser, it is handled differently than
    the f' retrieval.
    """
    if (not self.use_fdp[i_seq]) :
      return None
    if (self.fpp_from_phaser_ax_sites is None) :
      return self.site_fdp[i_seq]
    fpp = self.fpp_from_phaser_ax_sites[i_seq]
    if (fpp < 0) : fpp = None
    return fpp

  def get_fp (self, i_seq) :
    """
    Retrieve the refined f' for a site.
    """
    if (self.site_fp is None) or (not self.use_fdp[i_seq]) :
      return None
    return self.site_fp[i_seq]

  # XXX approximate guess, needs more rational parameters
  def looks_like_halide_ion (self,
      i_seq,
      element = "CL",
      assume_hydrogens_all_missing = Auto):
    """
    Given a site, analyze the nearby atoms (taking geometry into account) to
    guess whether it might be a halide ion.  Among other things, halides will
    often be coordinated by amide hydrogens, and not in close contact with
    negatively charged atoms.  Note that this procedure has a very high
    false positive rate when used by itself, so additional information about
    the electron count (map, occ, b) is absolutely essential.
    """
    atom = self.pdb_atoms[i_seq]
    sites_frac = self.xray_structure.sites_frac()
    #print "checking for halide:", atom.id_str()
    assert element.upper() in HALIDES
    # discard atoms with B-factors greater than mean-1sigma for waters
    if ((self.b_mean_hoh is not None) and
        (atom.b > self.b_mean_hoh - self.b_stddev_hoh)) :
      #print "Bad b-factor"
      return False

    params = self.params.chloride
    nearby_atoms = self.find_nearby_atoms(
      i_seq,
      far_distance_cutoff = 4.0,
      filter_by_bonding = False)

    if (assume_hydrogens_all_missing in [Auto,None]) :
      xrs = self.xray_structure
      sctr_keys = xrs.scattering_type_registry().type_count_dict().keys()
      assume_hydrogens_all_missing = "H" not in sctr_keys and \
                                     "D" not in sctr_keys

    binds_amide_hydrogen = False
    near_cation = False
    near_lys = False
    near_hydroxyl = False
    xyz = col(atom.xyz)
    min_distance_to_cation = max(self.unit_cell.parameters()[0:3])
    min_distance_to_hydroxyl = min_distance_to_cation

    for contact in nearby_atoms:
      # to analyze local geometry, we use the target site mapped to be in the
      # same ASU as the interacting site
      def get_site (k_seq) :
        return self.unit_cell.orthogonalize(
          site_frac = (contact.rt_mx * sites_frac[k_seq]))
      other = contact.atom
      resname = contact.resname()
      atom_name = contact.atom_name()
      element = contact.element
      distance = abs(contact)
      j_seq = other.i_seq

      # XXX need to figure out exactly what this should be - CL has a
      # fairly large radius though (1.67A according to ener_lib.cif)
      if (distance < params.min_distance_to_other_sites) :
        return False

      if not element in ["C", "N", "H", "O", "S"]:
        charge = self.server.get_charge(element)

        if charge < 0 and distance <= params.min_distance_to_anion:
          # Nearby anion that is too close
          return False

        if charge > 0 and distance <= params.max_distance_to_cation:
          # Nearby cation
          near_cation = True
          if (distance < min_distance_to_cation) :
            min_distance_to_cation = distance
      # Lysine sidechains (can't determine planarity)
      elif (atom_name in ["NZ"] and #, "NE", "NH1", "NH2"] and
            resname in ["LYS"] and
            distance <= params.max_distance_to_cation):
        near_lys = True
        if (distance < min_distance_to_cation) :
          min_distance_to_cation = distance
      # sidechain amide groups, no hydrogens (except Arg)
      # XXX this would be more reliable if we also calculate the expected
      # hydrogen positions and use the vector method below
      elif (atom_name in ["NZ","NH1","NH2","ND2","NE2"] and
            resname in ["ARG","ASN","GLN"] and
            (assume_hydrogens_all_missing or resname == "ARG") and
            distance <= params.max_distance_to_cation) :
        if (is_coplanar_with_sidechain(atom, other.parent(),
              distance_cutoff = params.max_deviation_from_plane)) :
          binds_amide_hydrogen = True
          if (resname == "ARG") and (distance < min_distance_to_cation) :
            min_distance_to_cation = distance
      # hydroxyl groups - note that the orientation of the hydrogen is usually
      # arbitrary and we can't determine precise bonding
      elif ((atom_name in ["OG1", "OG2", "OH1"]) and
            (resname in ["SER", "THR", "TYR"]) and
            (distance <= params.max_distance_to_hydroxyl)) :
        near_hydroxyl = True
        if (distance < min_distance_to_hydroxyl) :
          min_distance_to_hydroxyl = distance

      # Backbone amide, explicit H
      elif atom_name in ["H"]:
        # TODO make this more general for any amide H?
        xyz_h = col(contact.site_cart)
        bonded_atoms = self.connectivity[j_seq]
        if (len(bonded_atoms) != 1) :
          continue
        xyz_n = col(get_site(bonded_atoms[0]))
        vec_hn = xyz_h - xyz_n
        vec_hx = xyz_h - xyz
        angle = abs(vec_hn.angle(vec_hx, deg = True))

        # If Cl, H, and N line up, Cl binds the amide group
        if abs(angle - 180) <= params.delta_amide_h_angle:
          binds_amide_hydrogen = True
        else :
          pass #print "%s N-H-X angle: %s" % (atom.id_str(), angle)
      # Backbone amide, implicit H
      elif atom_name in ["N"] and assume_hydrogens_all_missing:
        xyz_n = col(contact.site_cart)
        bonded_atoms = self.connectivity[j_seq]
        ca_same = c_prev = None

        for k_seq in bonded_atoms :
          other2 = self.pdb_atoms[k_seq]
          if other2.name.strip().upper() in ["CA"]:
            ca_same = col(get_site(k_seq))
          elif other2.name.strip().upper() in ["C"]:
            c_prev = col(get_site(k_seq))

        if ca_same is not None and c_prev is not None:
          xyz_cca = (ca_same + c_prev) / 2
          vec_ncca = xyz_n - xyz_cca
          # 0.86 is the backbone N-H bond distance in geostd
          xyz_h = xyz_n + (vec_ncca.normalize() * 0.86)
          vec_nh = xyz_n - xyz_h
          vec_nx = xyz_n - xyz
          angle = abs(vec_nh.angle(vec_nx, deg = True))

          if abs(angle - 180) <= params.delta_amide_h_angle:
            binds_amide_hydrogen = True

      # sidechain NH2 groups, explicit H
      elif ((atom_name in ["HD1","HD2"] and resname in ["ASN"]) or
            (atom_name in ["HE1","HE2"] and resname in ["GLN"])) :
            # XXX not doing this for Arg because it can't handle the bidentate
            # coordination
            #(atom_name in ["HH11","HH12","HH21","HH22"] and resname == "ARG")):
        bonded_atoms = self.connectivity[j_seq]
        assert (len(bonded_atoms) == 1)
        xyz_n = col(get_site(bonded_atoms[0]))
        xyz_h = col(contact.site_cart)
        vec_nh = xyz_n - xyz_h
        vec_xh = xyz - xyz_h
        angle = abs(vec_nh.angle(vec_xh, deg = True))

        if abs(angle - 180) <= params.delta_amide_h_angle:
          binds_amide_hydrogen = True
        else :
          pass #print "%s amide angle: %s" % (atom.id_str(), angle)

    # now check again for negatively charged sidechain (etc.) atoms (e.g.
    # carboxyl groups), but with some leeway if a cation is also nearby.
    # backbone carbonyl atoms are also excluded.
    for contact in nearby_atoms:
      if (contact.altloc() not in ["", "A"]) :
        continue
      resname = contact.resname()
      atom_name = contact.atom_name()
      distance = abs(contact)
      if ((distance < 3.2) and
          (distance < (min_distance_to_cation + 0.2)) and
          is_negatively_charged_oxygen(atom_name, resname)) :
        #print contact.id_str(), distance
        return False

    # TODO something smart...
    # the idea here is to determine whether the blob of density around the
    # site is approximately spherical and limited in extent.
    pai = self.principal_axes_of_inertia(i_seq)
    #print list(pai.eigensystem().values())
    map_variance = self.get_map_sphere_variance(i_seq)
    #map_variance.show(prefix = "  ")
    good_map_falloff = False
    #print map_variance.mean, map_variance.standard_deviation
    if ((map_variance is not None) and
        (map_variance.mean < 1.5) and
        (map_variance.standard_deviation < 0.4)) : # XXX very arbitrary
      good_map_falloff = True

    # XXX probably need something more sophisticated here too.  a CL which
    # coordinates a metal (e.g. in 4aqi) may not have a clear blob, but
    # will still be detectable by other criteria.
    return (binds_amide_hydrogen or near_cation or near_lys)
           # good_map_falloff)

  def analyze_water(self,
      i_seq,
      debug = True,
      candidates = Auto,
      no_final = False,
      out = sys.stdout):
    """
    Examines the environment around a single atom to determine if it is actually
    a misidentified metal ion.

    no_final is used internally, when candidates is not Auto, to try all of the
    Auto without selecting any as a final choice when none of the elements in
    candidates appear to probable ion identities.
    """

    atom = self.pdb_atoms[i_seq]
    # Keep track of this in case we find nothing from the user-specific candidates
    auto_candidates = candidates is Auto
    if auto_candidates:
      candidates = DEFAULT_IONS
    candidates = [i.strip().upper() for i in candidates]
    if (candidates == ['X']) : # XXX hack for testing - X is "dummy" element
      candidates = []

    resname = atom.fetch_labels().resname.strip().upper()
    assert resname in WATER_RES_NAMES

    atom_props = self.atoms_to_props[i_seq]
    final_choice = None

    # Gather some quick statistics on the atom and see if it looks like water
    atom_type = atom_props.get_atom_type(params=self.params.water,
      aggressive=self.params.aggressive)
    if (atom_type == WATER_POOR) : # not trustworthy, skip
      #PRINT_DEBUG("SKIPPING %s" % atom.id_str())
      return None

    # Filter out metals based on whether they are more or less eletron-dense
    # in comparison with water
    filtered_candidates = []
    halide_candidates = False
    nuc_phosphate_site = atom_props.looks_like_nucleotide_phosphate_site()
    max_carbon_fo_map_value = sys.maxint
    if (self.carbon_fo_values is not None) :
      max_carbon_fo_map_value = flex.max(self.carbon_fo_values)

    for symbol in candidates :
      elem = self.server.get_metal_parameters(symbol)
      if (elem is None) :
        if (symbol in HALIDES) : # halides are special!
          halide_candidates = True
          continue
        else :
          raise Sorry("Element '%s' not supported!" % symbol)

      # If we definitely look like water or a light ion, only look at the
      # metals isoelectronic with water.
      if ((atom_type in [WATER, LIGHT_ION]) and
          (symbol not in ["NA", "MG", "F", "NE"])) :
        continue

      if (nuc_phosphate_site) and (not symbol in NUC_PHOSPHATE_BINDING) :
        continue

      n_elec = sasaki.table(symbol.upper()).atomic_number() - elem.charge
      # lighter elements are not expected to have any anomalous signal
      if (n_elec <= 12) and (atom_props.fpp > 0.1) :
        continue
      mass_ratio = atom_props.estimated_weight / max(n_elec, 1)
      # note that anomalous peaks are more important than the 2mFo-DFc level
      if (mass_ratio < 0.4) and (atom_type == WATER) :
        continue

      filtered_candidates.append(elem)

    # if len(filtered_candidates) == 0 and not halide_candidates:
    #   return None

    # Try each different ion and see what is reasonable
    def try_candidates (require_valence = True) :
      reasonable = []
      unreasonable = []
      for elem_params in filtered_candidates:
        # Try the ion with check_ion_environment()
        atom_props.check_ion_environment(
          ion_params = elem_params,
          server = self.server,
          wavelength = self.wavelength,
          require_valence = require_valence)
        atom_props.check_fpp_ratio(
          ion_params = elem_params,
          wavelength = self.wavelength,
          fpp_ratio_min = self.params.phaser.fpp_ratio_min,
          fpp_ratio_max = self.params.phaser.fpp_ratio_max)
        identity = str(elem_params)
        if atom_props.is_correctly_identified(identity = identity):
          reasonable.append((elem_params, atom_props.score[identity]))
        else :
          unreasonable.append((elem_params, atom_props.score[identity]))
      return reasonable, unreasonable

    # first try with bond valence included in criteria - then if that doesn't
    # yield any hits, try again without valences
    reasonable, unreasonable = try_candidates(require_valence = True)
    valence_used = True
    if (len(reasonable) == 0) and (not self.params.require_valence) :
      reasonable, unreasonable = try_candidates(require_valence = False)
      if (len(reasonable) > 0) :
        valence_used = False

    looks_like_halide = False
    if (not (reasonable or atom_type < HEAVY_ION or nuc_phosphate_site)) :
      # try halides now
      candidate_halides = set(candidates).intersection(HALIDES)
      filtered_halides = []
      for element in candidate_halides :
        fpp_ratio = atom_props.check_fpp_ratio(
          ion_params = MetalParameters(element = element, charge = -1),
          wavelength = self.wavelength,
          fpp_ratio_min = self.params.phaser.fpp_ratio_min,
          fpp_ratio_max = self.params.phaser.fpp_ratio_max)
        # XXX chlorides are tricky, because they tend to be partial occupancy
        # anyway and the f'' is already small, so prone to error here
        if (fpp_ratio is not None) and (element != "CL") :
          if ((fpp_ratio < self.params.phaser.fpp_ratio_min) or
              (fpp_ratio > self.params.phaser.fpp_ratio_max)) :
            #print "fpp_ratio:", fpp_ratio
            continue
        if (self.looks_like_halide_ion(i_seq = i_seq, element = element)) :
          filtered_halides.append(element)

      looks_like_halide = (len(filtered_halides) > 0)
      reasonable += [(MetalParameters(element = halide, charge = -1), 0)
                     for halide in filtered_halides]

    # If we can't find anything reasonable, relax our constraints...
    # Look for something that is a compatible ligand, an apparent fpp,
    # and compatible geometries
    if (len(reasonable) == 0) :
      compatible = [params for params in filtered_candidates
                    if atom_props.has_compatible_ligands(str(params))]
      for ion_params in compatible :
        atomic_number = sasaki.table(ion_params.element).atomic_number()
        weight_ratio = 0
        if (atom_props.estimated_weight is not None) :
          weight_ratio = atom_props.estimated_weight / (atomic_number - \
            ion_params.charge)
        # special handling for transition metals, but only if the user has
        # explicitly requested one (and there is no ambiguity)
        if ((ion_params.element in TRANSITION_METALS) and
            (atom_type == HEAVY_ION) and
            #(atom_props.peak_2fofc > max_carbon_fo_map_value) and
            (not auto_candidates) and
            (atom_props.is_compatible_site(ion_params))) :
          n_good_res = atom_props.number_of_favored_ligand_residues(ion_params,
            distance = 2.7)
          n_total_coord_atoms = atom_props.number_of_atoms_within_radius(2.8)
          # if we see at least one favorable residue coordinating the atom
          # and no more than six atoms total, accept the current guess
          if ((n_good_res >= 1) and (2 <= n_total_coord_atoms <= 6)) :
            reasonable.append((ion_params, 0))
          else : pass
            #print "n_good_res = %d, n_total_coord_atoms = %d" % (n_good_res,
            #  n_total_coord_atoms)
        elif ((ion_params.element in ["K","CA"]) and
              (atom_type == HEAVY_ION) and
              (not auto_candidates) and
              (atom_props.is_compatible_site(ion_params))) :
          n_good_res = atom_props.number_of_favored_ligand_residues(ion_params,
            distance = 2.9, exclude_atoms = ["O"])
          n_bb_oxygen = atom_props.number_of_backbone_oxygens(
            distance_cutoff = 2.9)
          n_total_coord_atoms = atom_props.number_of_atoms_within_radius(
            distance_cutoff = 2.9)
          if ((n_good_res + n_bb_oxygen) >= 2) and (n_total_coord_atoms >= 4) :
            reasonable.append((ion_params, 0))
        # another special case: very heavy ions, which are probably not binding
        # physiologically
        elif ((atomic_number > 30) and (atom_type == HEAVY_ION) and
              (((weight_ratio > 0.5) and (weight_ratio < 1.05)) or
               (atom_props.fp > 10)) and
              (not auto_candidates) and
              (atom_props.is_compatible_site(ion_params,
                ignore_valence=not self.get_strict_valence_flag()))) :
          n_total_coord_atoms = atom_props.number_of_atoms_within_radius(2.8)
          if (n_total_coord_atoms >= 3) :
            reasonable.append((ion_params, 0))
        else : pass
          #print atom.id_str(), atomic_number, looks_like_water, weight_ratio, \
          #  atom_props.is_compatible_site(ion_params)
      if (len(compatible) == 1) and (not self.get_strict_valence_flag()) :
        inaccuracies = atom_props.inaccuracies[str(ion_params)]
        if (compatible[0] in atom_props.fpp_ratios and
            atom_props.BAD_FPP not in inaccuracies and
            not inaccuracies.intersection([atom_props.BAD_GEOMETRY,
                                           atom_props.NO_GEOMETRY])):
          PRINT_DEBUG(atom.id_str())
          PRINT_DEBUG(filtered_candidates)
          PRINT_DEBUG(compatible)

    if (len(reasonable) == 1) :
      final_choice = reasonable[0][0]

    if (not reasonable) and (not auto_candidates):
      # Couldn't find anything from what the user suggested, try the other
      # default candidates and just let the user know about them
      result = self.analyze_water(
        i_seq = i_seq,
        debug = debug,
        no_final = True,
        out = out)
      if (result.final_choice is not None) :
        return result
    return water_result(
      atom_props = atom_props,
      filtered_candidates = filtered_candidates,
      matching_candidates = reasonable,
      rejected_candidates = unreasonable,
      nuc_phosphate_site = nuc_phosphate_site,
      atom_type = atom_type,
      looks_like_halide = looks_like_halide,
      ambiguous_valence_cutoff = self.params.ambiguous_valence_cutoff,
      valence_used = valence_used,
      final_choice = final_choice,
      wavelength = self.wavelength,
      no_final = no_final)

  def analyze_waters (self, out = sys.stdout, debug = True, candidates = Auto):
    """
    Iterates through all of the waters in a model, examining the maps and local
    environment to check their identity and suggest likely ions where
    appropriate.

    Returns a list of (i_seq, MetalParameter), which indicate that waters of
    i_seq would be better changed to the new metal identity.
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
    print >> out, "  %d waters to analyze" % len(waters)
    if (len(waters) == 0) : return
    nproc = easy_mp.get_processes(self.nproc)
    ions = []
    self.atoms_to_props = dict((i_seq, AtomProperties(i_seq, self))
                               for i_seq in waters)
    if (nproc == 1) :
      print >> out, ""
      for water_i_seq in waters :
        t1 = time.time()
        water_props = self.analyze_water(
          i_seq = water_i_seq,
          debug = debug,
          candidates = candidates,
          out = out)
        if (water_props is not None) :
          water_props.show_summary(out = out, debug = debug)
          map_stats = self.map_stats(water_i_seq)
          if ((water_props.final_choice is not None) and
              (not water_props.no_final)) :
            ions.append((water_i_seq, [water_props.final_choice],
              map_stats.two_fofc))
        t2 = time.time()
        #print "%s: %.3fs" % (self.pdb_atoms[water_i_seq].id_str(), t2-t1)
    else :
      print >> out, "  Parallelizing across %d processes" % nproc
      print >> out, ""
      analyze_water = _analyze_water_wrapper(manager = self,
        debug = debug,
        candidates = candidates,
        out = out)
      results = easy_mp.pool_map(
        fixed_func = analyze_water,
        args = waters,
        processes = nproc)
      for result in results :
        if (result is None) :
          continue
        water_i_seq, final_choice, result_str = result
        if (result_str is not None) :
          print >> out, result_str
        if final_choice is not None :
          map_stats = self.map_stats(water_i_seq)
          ions.append((water_i_seq, [final_choice], map_stats.two_fofc))

    return sorted(ions, lambda a, b: cmp(b[2], a[2]))

  def validate_ion (self, i_seq, out = sys.stdout, debug = True):
    """
    Examines one site in the model and determines if it was correctly modelled,
    returning a boolean indicating correctness.
    """

    atom_props = self.atoms_to_props[i_seq]
    element = self.server.get_element(atom_props.atom)
    elem_params = self.server.get_metal_parameters(element)

    if elem_params is not None:
      atom_type = atom_props.get_atom_type(params=self.params.water)
      atom_props.check_ion_environment(
        ion_params = elem_params,
        server = self.server,
        wavelength = self.wavelength,
        require_valence = self.params.require_valence)
      atom_props.check_fpp_ratio(
        ion_params = elem_params,
        wavelength = self.wavelength,
        fpp_ratio_min = self.params.phaser.fpp_ratio_min,
        fpp_ratio_max = self.params.phaser.fpp_ratio_max)
    elif element in HALIDES:
      identity = atom_props.identity()
      atom_props.inaccuracies[identity] = set()

      if not self.looks_like_halide_ion(i_seq = i_seq, element = element):
        atom_props.inaccuracies[identity].add(atom_props.BAD_HALIDE)
    else:
      raise Sorry("Element '%s' not supported:\n%s" %
                  (element, atom_props.atom.format_atom_record()))

    return atom_props

  def validate_ions (self, out = sys.stdout, debug = True, segid = None):
    """
    Iterate over all the ions built into the model by this module and double
    check their correctness. Looks for tell-tale signs such as negative peaks
    in the mFo-DFc and anomalous difference maps, abnormal b-factor, and
    disallowed chemical environments.

    Prints out a table of bad ions and their information and returns a list of
    their sites' i_seqs.
    """

    if segid is None:
      ions = []
      for model in self.pdb_hierarchy.models():
        for chain in model.chains():
          for residue_group in chain.residue_groups():
            for atom_group in residue_group.atom_groups():
              if atom_group.resname.strip() in SUPPORTED:
                atoms = atom_group.atoms()
                assert (len(atoms) == 1)
                for atom in atoms:
                  ions += atom.i_seq,
    else:
      sel_cache = self.pdb_hierarchy.atom_selection_cache()
      ions = sel_cache.selection("segid {0}".format(segid)).iselection()

    if len(ions) == 0:
      return
    ion_status = []

    self.atoms_to_props = dict((i_seq, AtomProperties(i_seq, self))
                               for i_seq in ions)
    for i_seq in ions:
      atom_props = self.validate_ion(
        i_seq = i_seq,
        debug = debug)
      ion_status.append((atom_props, atom_props.is_correctly_identified()))
    scatterers = self.xray_structure.scatterers()
    headers = ("atom", "occ", "b_iso", "2mFo-DFc", "mFo-DFc", "fp", "fdp",
               "ratio", "BVS", "VECSUM")
    fmt = "%-15s %-5s  %-5s  %-8s  %-7s  %-5s  %-5s  %-5s %-5s  %-6s"
    box = framed_output(out, title = "Validating new ions", width = 80)
    print >> box, fmt % headers
    print >> box, " " + ("-" * 75)
    for props, okay_flag in ion_status :
      i_seq = props.atom.i_seq
      sc = scatterers[i_seq]
      fp = fdp = None
      if (sc.flags.use_fp_fdp) :
        fp = sc.fp
        fdp = sc.fdp
      # XXX props.atom.b does not work here!
      b_iso = adptbx.u_as_b(sc.u_iso_or_equiv(unit_cell = self.unit_cell))
      identity = props.identity()
      def ff (fs, val) : return format_value(fs, val, replace_none_with = "---")
      print >> box, fmt % (props.atom.id_str(suppress_segid = True)[5:-1],
        ff("%.2f", props.atom.occ), ff("%.2f", b_iso),
        ff("%.2f", props.peak_2fofc), ff("%.2f", props.peak_fofc),
        ff("%.2f", fp), ff("%.2f", fdp),
        ff("%.2f", None if self.wavelength is None else fdp / sasaki.table(
            self.server.get_element(props.atom)).at_angstrom(self.wavelength).fdp()),
        ff("%.2f", props.valence_sum.get(identity)),
        ff("%.2f", props.vector_sum.get(identity)))
      # print >> box, props.geometries
      if not okay_flag:
        print >> box, "\n".join("!!! " + props.error_strs[i]
                                for i in props.inaccuracies[identity])
      print >> box, ""
    box.close()

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
    try:
      result = self.manager.analyze_water(i_seq, **(self.kwds))
      out = cStringIO.StringIO()
      if (result is not None) :
        result.show_summary(out = out,
          debug = self.kwds.get("debug", False))
      result_str = out.getvalue()
      if (result_str == "") : result_str = None
      final_choice = None
      # only accept final choice if it's in the original list of elements
      if (not getattr(result, "no_final", False)) :
        final_choice = getattr(result, "final_choice", None)
        #PRINT_DEBUG("final_choice = %s" % final_choice)
      return i_seq, final_choice, result_str
    except KeyboardInterrupt :
      return (None, None, None)

class _validate_ion_wrapper (object) :
  """
  Simple wrapper for calling manager.validate_ion with keyword arguments
  in a parallelized loop.  Because the water_result object is not pickle-able,
  we only return the final ion choice and summary string.
  """
  def __init__ (self, manager, **kwds) :
    self.manager = manager
    self.kwds = dict(kwds)

  def __call__ (self, i_seq) :
    try:
      correct = self.manager.validate_ion(i_seq, **(self.kwds))
      return i_seq, correct
    except KeyboardInterrupt :
      return (None, None)

class water_result (object):
  """
  Container for storing the results of manager.analyze_water for later display
  and retrieval.
  """
  def __init__ (self,
      atom_props,
      filtered_candidates,
      matching_candidates,
      rejected_candidates,
      nuc_phosphate_site,
      atom_type,
      looks_like_halide,
      ambiguous_valence_cutoff,
      valence_used,
      final_choice,
      wavelength=None,
      no_final = False) :
    adopt_init_args(self, locals())

  def show_summary (self, out = None, debug = False) :
    if (out is None) : out = sys.stdout
    results = self.matching_candidates
    if (len(results) > 0) :
      self.atom_props.show_properties(identity = "HOH", out = out)
      if (self.nuc_phosphate_site) :
        print >> out, "  appears to be nucleotide coordination site"
      if (self.no_final) :
        print >> out, "  Found potential ion%s outside of specified set:" % \
          ("s" if len(results) > 1 else "")
      if (self.final_choice is not None) :
        # We have one result that we are reasonably certain of
        elem_params, score = results[0]
        if elem_params.element not in HALIDES:
          self.atom_props.show_ion_results(
            identity = str(self.final_choice),
            out = out,
            valence_used = self.valence_used,
            confirmed = True)
        else:
          print >> out, "  Probable anion:", str(elem_params)
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
            out = out,
            valence_used = True)
          print >> out, ""
        else:
          ions = [str(i[0]) for i in sorted(results, key = lambda x: x[1])]
          print >> out, "  ambiguous results, could be %s" % ", ".join(ions)
          for elem_params, score in results :
            self.atom_props.show_ion_results(identity = str(elem_params),
              out = out)
          print >> out, ""
    else:
      if (self.atom_type != WATER) or (self.nuc_phosphate_site) :
        self.atom_props.show_properties(identity = "HOH", out = out)
        if (self.nuc_phosphate_site) :
          print >> out, "  appears to be nucleotide coordination site"
        # try anions now
        if (self.looks_like_halide) :
          print >> out, "  Probable cation: %s" % str(self.final_choice)
          print >> out, ""
        else:
          # atom is definitely not water, but no reasonable candidates found
          # print out why all the metals we tried failed
          if (debug) and (len(self.filtered_candidates) > 0) :
            print >> out, "  insufficient data to identify atom"
            possible = True
            for params in self.filtered_candidates:
              if (self.atom_props.has_compatible_ligands(str(params))) :
                if possible:
                  print >> out, "  possible candidates:"
                  possible = False
                self.atom_props.show_ion_results(identity = str(params),
                  out = out)
              else :
                print >> out, "  incompatible ligands for %s" % str(params)
            #print >> out, "  rejected as unsuitable:"
            #for params in self.rejected_candidates:
            #  if (self.atom_props.has_compatible_ligands(str(params))) :
            #    self.atom_props.show_ion_results(identity = str(params),
            #      out = out)
            #  else :
            #    print >> out, "  incompatible ligands for %s" % str(params)
          print >> out, ""

class AtomProperties (object) :
  """
  Collect physical attributes of an atom, including B, occupancy, and map
  statistics, and track those which are at odds with its chemical identity.
  """

  LOW_B, HIGH_B, LOW_OCC, HIGH_OCC, NO_2FOFC_PEAK, FOFC_PEAK, FOFC_HOLE, \
    ANOM_PEAK, NO_ANOM_PEAK, BAD_GEOMETRY, NO_GEOMETRY, BAD_VECTORS, \
    BAD_VALENCES, TOO_FEW_NON_WATERS, TOO_FEW_COORD, TOO_MANY_COORD, \
    LIKE_COORD, BAD_COORD_ATOM, BAD_FPP, BAD_COORD_RESIDUE, VERY_BAD_VALENCES, \
    BAD_HALIDE, HIGH_2FOFC, COORDING_GEOMETRY, CLOSE_CONTACT \
    = range(25)

  error_strs = {
    LOW_B: "Abnormally low b-factor",
    HIGH_B: "Abnormally high b-factor",
    LOW_OCC: "Abnormally low occupancy",
    HIGH_OCC: "Abnormally high occupancy",
    NO_2FOFC_PEAK: "No 2mFo-DFc map peak",
    FOFC_PEAK: "Peak in mFo-DFc map",
    FOFC_HOLE: "Negative peak in dFo-mFc map",
    ANOM_PEAK: "Peak in the anomalous map",
    NO_ANOM_PEAK: "No peak in the anomalous map",
    BAD_GEOMETRY: "Unexpected geometry for coordinating atoms",
    NO_GEOMETRY: "No distinct geometry for coordinating atoms",
    BAD_VECTORS: "VECSUM above cutoff",
    BAD_VALENCES: "BVS above or below cutoff",
    TOO_FEW_NON_WATERS: "Too few non-water coordinating atoms",
    TOO_FEW_COORD: "Too few coordinating atoms",
    TOO_MANY_COORD: "Too many coordinating atoms",
    LIKE_COORD: "Like charge coordinating like",
    BAD_COORD_ATOM: "Disallowed coordinating atom",
    BAD_FPP: "Bad refined f'' value",
    BAD_COORD_RESIDUE: "Disallowed coordinating residue",
    VERY_BAD_VALENCES: "BVS far above or below cutoff",
    BAD_HALIDE: "Bad halide site",
    HIGH_2FOFC: "Unexpectedly high 2mFo-DFc value",
    COORDING_GEOMETRY: "No distinct geometry and coordinating another atom with distinct geometry",
    CLOSE_CONTACT : "Close contact to oxygen atom",
    }

  def __init__(self, i_seq, manager):
    self.i_seq = i_seq
    self.atom = manager.pdb_atoms[i_seq]
    self.resname = self.atom.parent().resname.strip().upper()
    self.d_min = manager.fmodel.f_obs().d_min()
    self.anomalous_flag = manager.fmodel.f_obs().anomalous_flag()
    self.strict_valence = manager.get_strict_valence_flag()

    # Grab all the atoms within 3.5 Angstroms
    nearby_atoms_unfiltered = manager.find_nearby_atoms(
      i_seq = i_seq,
      far_distance_cutoff = 3.5)
    self.nearby_atoms = []
    nearby_atoms_no_alts = []
    for contact in nearby_atoms_unfiltered :
      if (contact.element not in ["H", "D"]) :
        self.nearby_atoms.append(contact)
        for other in nearby_atoms_no_alts :
          if (other == contact) :
            break
        else :
          nearby_atoms_no_alts.append(contact)

    self.residue_counts = count_coordinating_residues(self.nearby_atoms)
    self.geometries = find_coordination_geometry(nearby_atoms_no_alts)

    map_stats = manager.map_stats(i_seq)
    self.peak_2fofc = map_stats.two_fofc
    self.peak_fofc = map_stats.fofc
    self.peak_anom = map_stats.anom
    self.estimated_weight = manager.guess_molecular_weight(i_seq)
    self.b_iso = manager.get_b_iso(i_seq)
    # cache global properties
    self.b_mean_hoh = manager.b_mean_hoh
    self.b_stddev_hoh = manager.b_stddev_hoh
    self.b_mean_calpha = manager.b_mean_calpha
    self.calpha_mean_two_fofc = manager.calpha_mean_two_fofc

    self.inaccuracies = {}
    self.ignored = {}
    self.vectors = {}
    self.valence_sum = {}
    self.vector_sum = {}
    self.score = {}
    self.fpp_expected = {}
    self.expected_params = {}
    self.bad_coords = {}

    # Determine the f'' value if possible using phaser or anomalous refinement
    self.fpp = manager.get_fpp(i_seq)
    self.fpp_ratios = {}
    self.fp = manager.get_fp(i_seq)
    self.manager = manager

  def is_correctly_identified(self, identity = None):
    """
    Returns whether factors indicate that the atom was correctly identified.
    """
    if identity is None:
      identity = self.identity()

    return len(self.inaccuracies[identity]) == 0

  def has_compatible_ligands (self, identity) :
    """
    Indicates whether the coordinating atoms are of the allowed type (e.g.
    no N or S atoms coordinating CA, etc.) and residue (e.g. Ser is not an
    appropriate ligand for ZN).
    """
    return ((len(self.bad_coords[identity]) == 0) and
            (not self.BAD_COORD_RESIDUE in self.inaccuracies[identity]))

  def is_compatible_site (self, ion_params, require_anom = True,
      ignore_valence=False) :
    """
    More minimal criteria for determining whether a site is chemically
    compatible, allowing for incomplete coordination shells.
    """
    inaccuracies = self.inaccuracies[str(ion_params)]
    anom_allowed = self.is_compatible_anomalous_scattering(ion_params) or \
                   (not require_anom)
    return (self.has_compatible_ligands(str(ion_params)) and
            (not self.TOO_MANY_COORD in inaccuracies) and
            (ignore_valence or (not self.VERY_BAD_VALENCES in inaccuracies))
            and (not self.BAD_COORD_RESIDUE in inaccuracies) and
            (anom_allowed))

  def number_of_favored_ligand_residues (self, ion_params, distance = 3.0,
      exclude_atoms = ()) :
    """
    Counts the number of preferred residues coordinating the atom.  Used for
    approximate detection of transition-metal binding sites.
    """
    n_res = 0
    resids = []
    for contact in self.nearby_atoms:
      if (contact.atom_name() in exclude_atoms) :
        continue
      if (contact.distance() < distance) :
        labels = contact.atom.fetch_labels()
        other_resname = contact.resname()
        other_resid = labels.chain_id + labels.resid()
        if ((ion_params.allowed_coordinating_residues is not None) and
            (other_resname in ion_params.allowed_coordinating_residues) and
            (not other_resid in resids)) :
          n_res += 1
          resids.append(other_resid)
    return n_res

  def number_of_atoms_within_radius (self, distance_cutoff) :
    n_atoms = 0
    atom_ids = []
    for contact in self.nearby_atoms:
      other_id = contact.atom_id_no_altloc()
      if (not other_id in atom_ids) :
        if (contact.distance() < distance_cutoff) :
          n_atoms += 1
        atom_ids.append(other_id) # check for alt confs.
    return n_atoms

  def number_of_backbone_oxygens (self, distance_cutoff = 3.0) :
    n_bb_ox = 0
    for contact in self.nearby_atoms :
      if (contact.atom_name() == "O") :
        if (contact.distance() <= distance_cutoff) :
          if (not contact.resname() in WATER_RES_NAMES) :
            n_bb_ox += 1
    return n_bb_ox

  # FIXME needs to be refactored and combined with check_fpp_ratio
  def is_compatible_anomalous_scattering (self, ion_params) :
    # lighter elements should have effectively no anomalous scattering
    if (ion_params.element.upper() in ["MG", "NA"]) :
      return ((self.fpp is None) and (self.peak_anom is not None) and
              (self.peak_anom < 1.0))
    else :
      # XXX somewhat dangerous - we really need f'' for this to work reliably
      if (self.fpp is None) :
        return (self.peak_anom is not None) and (self.peak_anom > 3.0)
      identity = self.identity(ion = ion_params)
      if (identity in self.fpp_ratios) :
        return (not self.BAD_FPP in self.inaccuracies[identity])
    return False

  # XXX obsolete, delete?
  def atom_weight (self, manager):
    """
    Evaluates whether factors indicate that the atom is lighter, heavier, or
    isoelectric to what it is currently identified as.

    Returns -1 if lighter, 0 if isoelectronic, and 1 if heavier.
    """
    identity = "HOH" if self.resname in WATER_RES_NAMES else self.identity()

    # Waters that don't have B-factors at least 1 stddev below the mean are
    # presumed to be correct
    if (identity == "HOH" and
        (self.atom.b > manager.b_mean_hoh - manager.b_stddev_hoh)):
      return 0

    if self.is_correctly_identified(identity = identity):
      return 0

    # B-factors/occupancies?
    if self.FOFC_PEAK in self.inaccuracies[identity] or self.atom.b < 1:
      return 1

    if self.FOFC_HOLE in self.inaccuracies[identity]:
      return -1

    return 0

  def check_ion_environment (self,
      ion_params,
      server,
      wavelength = None,
      require_valence = True):
    """
    Checks whether or not the specified ion satisfies the metal-coordination
    parameters, specified by ion_params, such as valence sum, geometry, etc.
    The criteria used here are quite strict, but many of the analyses are
    saved for later if we want to use looser critera.
    """
    from iotbx.pdb import common_residue_names_get_class as get_class

    identity = self.identity(ion_params)
    inaccuracies = self.inaccuracies[identity] = set()
    self.expected_params[identity] = ion_params
    ignored = self.ignored[identity] = set()

    # if the atom is clearly not a water, optionally relax some rules.  this
    # will be more sensitive for transition metals, without finding a lot of
    # spurious Mg/Na sites.
    strict_rules = require_valence or \
      self.is_correctly_identified(identity = "HOH") or \
      self.strict_valence or \
      ion_params.element in ["NA","MG"]

    # Check for all non-overlapping atoms within 3 A of the metal
    n_closest = 0
    coord_atoms = []
    for i_pair, contact1 in enumerate(self.nearby_atoms) :
      distance = contact1.distance()
      if (distance < 3.0) :
        for contact2 in self.nearby_atoms[(i_pair+1):] :
          if ((contact1 == contact2) or
              (contact1.distance_from(contact2) <= 0.3)) :
            break
        else :
          coord_atoms.append(contact1)
          if (distance < 2.7) :
            n_closest += 1

    if len(coord_atoms) < ion_params.coord_num_lower:
      inaccuracies.add(self.TOO_FEW_COORD)

    if n_closest > ion_params.coord_num_upper:
      inaccuracies.add(self.TOO_MANY_COORD)

    # Coordinating atoms closer than 3.0 A are not positively charged
    n_non_water = 0
    self.bad_coords[identity] = []

    for contact in self.nearby_atoms:
      other_name = contact.atom_name()
      other_resname = contact.resname()
      other_element = contact.element

      if (not other_resname in WATER_RES_NAMES) :
        n_non_water += 1
      else:
        # Everything can potentially be coordinated by water
        continue

      if (contact.distance() < 3.0) :
        # XXX: So, we have a a fair number of rules restricting nitrogens and
        # nitrogen-containing residues from coordinating a number of cations.
        #
        # However, this rule is dependent on the protonation of the nitrogen,
        # if the pKa is low at the site, it is possible for a metal to
        # coordinate the residue fine.
        #
        # We want a complex rule that takes into account coordinating geometry,
        # density signal, and the presence of other coordinating atoms that
        # might drop the site's pKa enough to lose the hydrogen.
        if ((ion_params.allowed_coordinating_atoms is not None) and
            (other_element not in ion_params.allowed_coordinating_atoms)) :
          self.bad_coords[identity].append(contact)
          inaccuracies.add(self.BAD_COORD_ATOM)
        if (get_class(other_resname) == "common_amino_acid") :
          # limit elements allowed to bind to backbone atoms (mainly carbonyl
          # oxygen)
          if ((other_name in ["C","N","O","CA","H","HA"]) and
              ((ion_params.allowed_backbone_atoms is None) or
               (not other_name in ion_params.allowed_backbone_atoms))) :
            if (other_name == "O") and (is_carboxy_terminus(contact.atom)) :
              pass # C-terminal carboxyl group is allowed
            else :
              self.bad_coords[identity].append(contact)
              inaccuracies.add(self.BAD_COORD_ATOM)
          # Check if atom is of an allowed residue type, if part of a sidechain
          if (ion_params.allowed_coordinating_residues is not None) :
            allowed = ion_params.allowed_coordinating_residues
            if ((not other_resname in allowed) and
                (other_name not in ["C", "O", "N", "CA", "OXT"])) :
                # XXX probably just O
              self.bad_coords[identity].append(contact)
              inaccuracies.add(self.BAD_COORD_RESIDUE)
        elif (cmp(0, server.get_charge(contact.atom)) ==
              cmp(0, ion_params.charge)) :
          # Check if coordinating atom is of opposite charge
          self.bad_coords[identity].append(contact)
          inaccuracies.add(self.LIKE_COORD)
        elif (ion_params.charge > 0 and
            other_element in ["N"] and
            other_resname in ["LYS", "ARG", "ASN", "GLN"]):
          # Coordinating nitrogen most likely positive.
          #
          # Ignore nitrogens without a charge label that are on positively
          # charged amino acids.
          self.bad_coords[identity].append(contact)
          inaccuracies.add(self.LIKE_COORD)

    # Check the number of coordinating waters
    if (n_non_water < ion_params.min_coordinating_non_waters) :
      inaccuracies.add(self.TOO_FEW_NON_WATERS)

    # Check the geometry of the coordinating atoms
    if ion_params.allowed_geometries and strict_rules:
      allowed = [i[0] in ion_params.allowed_geometries
                 for i in self.geometries]
      if "any" in ion_params.allowed_geometries:
        pass
      elif not self.geometries:
        if strict_rules:
          inaccuracies.add(self.NO_GEOMETRY)
      elif not any(allowed):
        inaccuracies.add(self.BAD_GEOMETRY)
      else:
        strict_rules = False

    # If no distinct geometry, check that none of the coordinating have distinct
    # geometry, either
    if self.geometries == []:
      for contact in self.nearby_atoms:
        o_atom = contact.atom
        if o_atom.i_seq in self.manager.atoms_to_props:
          o_geometry = self.manager.atoms_to_props[o_atom.i_seq].geometries
          if o_geometry != []:
            inaccuracies.add(self.COORDING_GEOMETRY)

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

    # XXX I am not sure how low a valence sum we want to allow, but many
    # structures with non-physiological cation binding have partial and/or
    # irregular coordination shells
    if (self.valence_sum[identity] < ion_params.cvbs_expected * 0.25 or
        self.valence_sum[identity] > ion_params.cvbs_expected * 1.25):
      inaccuracies.add(self.VERY_BAD_VALENCES)
    else:
      if (self.valence_sum[identity] < ion_params.cvbs_lower or
          self.valence_sum[identity] > ion_params.cvbs_upper):
        if strict_rules:
          inaccuracies.add(self.BAD_VALENCES)
        else :
          ignored.add(self.BAD_VALENCES)

    self.score[identity] = abs(self.valence_sum[identity] -
                               ion_params.cvbs_expected)

  # FIXME this really needs to be refactored and combined with the method
  # is_compatible_anomalous_scattering
  def check_fpp_ratio (self,
      ion_params,
      wavelength,
      fpp_ratio_min = 0.3,
      fpp_ratio_max = 1.05) :
    """Compare the refined and theoretical f'' values if available"""
    identity = str(ion_params)
    inaccuracies = self.inaccuracies.get(identity, None)
    if (inaccuracies is None) :
      inaccuracies = self.inaccuracies[identity] = set()
    if (ion_params.element.upper() in ["MG", "NA"]) :
      if (self.fpp is not None) or (self.peak_anom > 1) :
        inaccuracies.add(self.BAD_FPP)
    else :
      # XXX in theory the fpp_ratio should be no more than 1.0 unless we are
      # right on the peak wavelength.  in practice Phaser can overshoot a little
      # bit, so we need to be more tolerant.  picking the maximum f'' from the
      # Sasaki and Henke tables will also limit the ratio.
      if (wavelength is not None) and (self.anomalous_flag) :
        fpp_expected_sasaki = sasaki.table(ion_params.element).at_angstrom(
          wavelength).fdp()
        fpp_expected_henke = henke.table(ion_params.element).at_angstrom(
          wavelength).fdp()
        self.fpp_expected[identity] = max(fpp_expected_sasaki,
          fpp_expected_henke)
        if (self.fpp is not None) and (self.fpp_expected[identity] != 0) :
          self.fpp_ratios[identity] = self.fpp / self.fpp_expected[identity]
          if ((self.fpp_ratios[identity] > fpp_ratio_max) or
              ((self.fpp >= 0.2) and
               (self.fpp_ratios[identity] < fpp_ratio_min))) :
            inaccuracies.add(self.BAD_FPP)
        elif (self.fpp_expected[identity] > 0.75) and (self.peak_anom < 2) :
          inaccuracies.add(self.BAD_FPP)
    return self.fpp_ratios.get(identity)

  def show_properties (self, identity, out = sys.stdout) :
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
    elif (self.HIGH_2FOFC in self.inaccuracies[identity]) :
      twofofc_flag = " <<<"
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
    if (self.estimated_weight is not None) :
      print >> out, "  Approx. mass:  %6d" % self.estimated_weight
    if self.fpp is not None:
      fpp_flag = ""
      if (self.fpp >= 0.2) :
        fpp_flag = " <<<"
      print >> out, "  f'':           %6.2f%s" % (self.fpp, fpp_flag)
      print >> out, "  f'' ratio:     %s" % format_value("%6.2f",
        self.fpp_ratios.get(identity))
    if self.nearby_atoms is not None:
      angstrom = u"\N{ANGSTROM SIGN}".encode("utf-8", "strict")
      degree = u"\N{DEGREE SIGN}".encode("utf-8", "strict")

      print >> out, "  Nearby atoms: (%d within 3.0 %s)" % \
        (len([i for i in self.nearby_atoms if i.distance() < 3]), angstrom)

      for contact in self.nearby_atoms :
        print >> out, "    %s (%5.3f %s)" % \
        (contact.id_str(), contact.distance(), angstrom)

      if self.geometries:
        print >> out, "  Coordinating geometry:"
        for geometry, deviation in self.geometries:
          print >> out, "    %-15s (average deviation: %.3f%s)" % \
            (geometry, deviation, degree)

  def show_ion_results (self, identity = None, out = sys.stdout,
      confirmed = False, valence_used = True):
    """
    Show statistics for a proposed element identity.
    """

    if not identity:
      identity = self.identity(self.atom)

    inaccuracies = self.inaccuracies.get(identity, set([]))
    ignored = self.ignored.get(identity, set([]))

    if identity != self.identity():
      if (confirmed) :
        print >> out, "  Probable cation: %s" % identity
      else :
        print >> out, "  Atom as %s:" % identity
    else:
      print >> out, "    %s:" % self.atom.id_str()

    if identity in self.vector_sum and self.vector_sum[identity] is not None:
      problem = ((self.BAD_VECTORS in inaccuracies) or
                 (self.BAD_VECTORS in ignored))

      print >> out, "    Vector sum:  %6.3f %s" % \
        (self.vector_sum[identity], " !!!" if problem else "")

    if identity in self.valence_sum and self.valence_sum[identity] is not None:
      problem = inaccuracies.union(ignored).intersection(
        [self.BAD_VALENCES, self.VERY_BAD_VALENCES])

      print >> out, "    Valence sum: %6.3f" % self.valence_sum[identity]
      if valence_used:
        print >> out, "(expected: %6.3f) %s" % \
        (self.expected_params[identity].cvbs_expected, " !!!" if problem else ""),

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
        print >> out, "    %s (%5.3f %s) !!!" % \
          (atom.id_str(), abs(vector), angstrom)

    if self.TOO_FEW_NON_WATERS in inaccuracies:
      print >> out, "    Too few coordinating non-waters !!!"
    if self.TOO_FEW_COORD in inaccuracies:
      print >> out, "    Too few coordinating atoms !!!"
    if self.TOO_MANY_COORD in inaccuracies:
      print >> out, "    Too many coordinating atoms !!!"

    if (self.fpp is not None) and (identity in self.fpp_ratios) :
      print >> out, "    f'' ratio:   %6.3f%s" % \
         (self.fpp_ratios[identity], " !!!" if self.BAD_FPP in inaccuracies else "")

  # XXX can we get away with just one oxygen?
  def looks_like_nucleotide_phosphate_site (self,
      min_phosphate_oxygen_atoms = 2,
      distance_cutoff = 2.5) : # XXX wild guess
    """
    Decide whether the atom is coordinating phosphate oxygens from a
    nucleotide, based on common atom names.
    """
    n_phosphate_oxygens = 0
    for contact in self.nearby_atoms :
      atom_name = contact.atom_name()
      if (len(atom_name) < 3) or (contact.element not in ["O"]) :
        continue
      if ((atom_name[0:2] in ["O1","O2","O3"]) and
          (atom_name[2] in ["A","B","G"])) :
        if (contact.distance() <= distance_cutoff) :
          n_phosphate_oxygens += 1
    return (n_phosphate_oxygens == min_phosphate_oxygen_atoms)

  def identity(self, ion = None):
    """
    Covers an atom into a string representing its element and charge.
    """
    if ion is None:
      ion = self.atom
    element = self.manager.server.get_element(ion)
    charge = self.manager.server.get_charge(ion)
    return "{}{:+}".format(element, charge)

  def get_atom_type (self, params, aggressive=False) :
    """
    Checks the atom characteristics against what we would expect for a water.
    Updates self with any inaccuracies noticed (Surpringly low b-factor,
    high occupancy, etc).  Note that HEAVY_ION does not necessarily rule out
    NA/MG, as these often have mFo-DFc peaks at high resolution.
    """
    inaccuracies = self.inaccuracies["HOH"] = set()
    atom_type = WATER
    # Skip over water if the 2mFo-DFc or mFo-DFc value is too low
    if ((self.peak_2fofc < params.min_2fofc_level) or
        (self.peak_fofc < -2.0)) :
      return WATER_POOR
    if (self.fpp is not None) and (self.fpp > params.fpp_max) :
      return HEAVY_ION
    if (self.fp is not None) and (self.fp > params.fp_max) :
      return HEAVY_ION
    if (self.peak_anom > params.max_anom_level) :
      inaccuracies.add(self.ANOM_PEAK)
      atom_type = HEAVY_ION
    if self.peak_fofc > params.max_fofc_level:
      inaccuracies.add(self.FOFC_PEAK)
      atom_type = HEAVY_ION
    if self.atom.occ > params.max_occ: # this will probably never happen...
      inaccuracies.add(self.HIGH_OCC)
      atom_type = HEAVY_ION
    # very low B-factors automatically trigger a check
    if (self.atom.b < params.min_b_iso) :
      inaccuracies.add(self.LOW_B)
      atom_type = HEAVY_ION
    elif (self.b_stddev_hoh is not None) and (self.b_stddev_hoh > 0) :
      # high B-factor relative to other waters
      z_value = (self.b_iso - self.b_mean_hoh) / self.b_stddev_hoh
      if z_value < -params.max_stddev_b_iso:
        inaccuracies.add(self.LOW_B)
      elif self.atom.b < self.b_mean_hoh * params.min_frac_b_iso:
        inaccuracies.add(self.LOW_B)
      if (atom_type == WATER) and (self.LOW_B in inaccuracies) :
        atom_type = LIGHT_ION
    if (aggressive) and (atom_type == WATER) :
      # high B-factor relative to C-alpha
      if (self.b_mean_calpha > 0) :
        relative_b = self.b_iso / self.b_mean_calpha
        if (relative_b < params.min_frac_calpha_b_iso) :
          inaccuracies.add(self.LOW_B)
          atom_type = LIGHT_ION
      # high 2Fo-Fc relative to C-alpha
      if (self.calpha_mean_two_fofc > 0) :
        relative_2fofc = self.peak_2fofc / self.calpha_mean_two_fofc
        if (relative_2fofc > params.max_frac_calpha_2fofc) :
          inaccuracies.add(self.HIGH_2FOFC)
          atom_type = LIGHT_ION
    # check for close contacts
    for i_pair, contact1 in enumerate(self.nearby_atoms) :
      if (contact1.element.strip() == "O") :
        distance = contact1.distance()
        if (distance < 2.4) :
          inaccuracies.add(self.CLOSE_CONTACT)
          if (atom_type < LIGHT_ION) :
            atom_type = LIGHT_ION
    return atom_type


def find_anomalous_scatterers (*args, **kwds) :
  """
  Wrapper for corresponding method in phaser.substructure, if phaser is
  available and configured.
  """
  if (not libtbx.env.has_module("phaser")) :
    if "log" in kwds:
      print >> kwds["log"], "Phaser not available"
    return None
  from phaser import substructure
  return substructure.find_anomalous_scatterers(*args, **kwds)

def create_manager (
    pdb_hierarchy,
    geometry_restraints_manager,
    fmodel,
    wavelength,
    params,
    resolution_factor = 0.25,
    nproc = Auto,
    verbose = False,
    log = None) :
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

def same_atom_different_altloc(atom1, atom2):
  """
  Determines whether atom1 and atom2 differ only by their alternate location.
  """

  label1, label2 = [i.fetch_labels() for i in [atom1, atom2]]
  name1, name2 = atom1.name.strip(), atom2.name.strip()
  chain1, chain2 = label1.chain_id, label2.chain_id
  res1, res2 = label1.resid(), label2.resid()

  return name1 == name2 and chain1 == chain2 and res1 == res2

def count_coordinating_residues (nearby_atoms, distance_cutoff = 3.0) :
  """
  Count the number of residues of each type involved in the coordination
  sphere.  This may yield additional clues to the identity of ions, e.g. only
  Zn will have 4 Cys residues.
  """
  unique_residues = []
  residue_counts = {}
  for contact in nearby_atoms :
    if (contact.distance() <= distance_cutoff) :
      parent = contact.atom.parent()
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

def PRINT_DEBUG (*args) :
  print >> sys.stderr, sys.stderr, args

def is_carboxy_terminus (pdb_object) :
  atoms = None
  if (type(pdb_object).__name__ == "atom") :
    atoms = pdb_object.parent().atoms()
  else :
    assert (type(pdb_object).__name__ in ['residue_group','atom_group'])
    atoms = pdb_object.atoms()
  for atom in atoms :
    if (atom.name.strip() == "OXT") :
      return True
  return False

def is_negatively_charged_oxygen (atom_name, resname) :
  """
  Determine whether the oxygen atom of interest is either negatively charged
  (usually a carboxyl group or sulfate/phosphate), or has a lone pair (and
  no hydrogen atom) that would similarly repel anions.
  """
  if ((atom_name in ["OD1","OD2","OE1","OE2"]) and
      (resname in ["GLU","ASP","GLN","ASN"])) :
    return True
  elif ((atom_name == "O") and (not resname in WATER_RES_NAMES)) :
    return True # sort of - the lone pair acts this way
  elif ((len(atom_name) == 3) and (atom_name[0:2] in ["O1","O2","O3"]) and
        (atom_name[2] in ["A","B","G"])) :
    return True
  elif (resname in ["SO4","PO4"]) :
    return True
  return False

# XXX distance cutoff may be too generous, but 0.5 is too strict
def is_coplanar_with_sidechain (atom, residue, distance_cutoff = 0.75) :
  """
  Given an isolated atom and an interacting residue with one or more amine
  groups, determine whether the atom is approximately coplanar with the terminus
  of the sidechain (and thus interacting with the amine hydrogen(s) along
  approximately the same axis as the N-H bond).
  """
  sidechain_sites = []
  resname = residue.resname
  for other in residue.atoms() :
    name = other.name.strip()
    if (resname == "ARG") and (name in ["NH1","NH2","NE"]) :
      sidechain_sites.append(other.xyz)
    elif (resname == "GLN") and (name in ["OE1","NE2","CD"]) :
      sidechain_sites.append(other.xyz)
    elif (resname == "ASN") and (name in ["OD1","ND2","CG"]) :
      sidechain_sites.append(other.xyz)
  if (len(sidechain_sites) != 3) : # XXX probably shouldn't happen
    return False
  D = distance_from_plane(atom.xyz, sidechain_sites)
  #print atom.id_str(), D
  return (D <= distance_cutoff)
