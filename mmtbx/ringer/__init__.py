
# TODO unit tests wouldn't hurt...

"""
Implementation of the Ringer method for torsion-angle sampling of sidechain
electron density to screen for alternate conformations.

Reference:
  Lang PT, Ng HL, Fraser JS, Corn JE, Echols N, Sales M, Holton JM, Alber T.
  Automated electron-density sampling reveals widespread conformational
  polymorphism in proteins. Protein Sci. 2010 Jul;19(7):1420-31. PubMed PMID:
  20499387
"""

from __future__ import absolute_import, division, print_function
from libtbx import adopt_init_args, Auto
from libtbx.utils import Sorry
from libtbx import easy_mp
import sys
from six.moves import range

ringer_phil_str = """
sampling_angle = 5
  .type = int
  .input_size = 80
scaling = *sigma volume
  .type = choice(multi=False)
skip_alt_confs = True
  .type = bool
  .short_caption = Skip existing alternate conformations
nproc = 1
  .type = int
  .short_caption = Processors
  .style = renderer:draw_nproc_widget
"""

class ringer_chi(object):
  """
  Sampling results for a single sidechain Chi angle.
  """
  def __init__(self, id, angle_current, densities, sampling,
      fofc_densities=None):
    adopt_init_args(self, locals())
    assert len(densities) > 0
    if (angle_current < 0):
      self.angle_current = 360 + angle_current
    self.peak_chi, self.peak_rho = self.find_peaks(densities)
    self.deviation = self.deviate(self.peak_chi)
    self.rho_mean = sum(densities) / len(densities)
    # Add a tiny number to avoid dividing by 0 (which shouldn't happen anyway)
    self.rho_rel = self.peak_rho/(self.rho_mean+.000000000000000001)

  def format_csv(self, fofc=False):
    densities = [ "%.3f" % x for x in self.densities ]
    if (fofc) and (self.fofc_densities is not None):
      densities = [ "%.3f" % x for x in self.fofc_densities ]
    return "chi%d,%.1f,%s" % (self.id, self.angle_current, ",".join(densities))

  def find_peaks(self, densities):
    rho_max = max(densities)
    for i, rho in enumerate(densities):
      if rho == max(densities):
        i = i * (360 / len(densities))
        return i, rho
    # This should never happen, but just in case, dump this in
    # place of throwing an error.
    return 0,0

  def deviate(self, chi):
    return min(abs(chi-i) for i in [60, 180, 300])

class ringer_residue(object):
  """
  Container of sampling results for a single residue with at least one Chi
  angle.
  """
  def __init__(self, resname, chain_id, resid, altloc, n_chi, xyz=None):
    adopt_init_args(self, locals())
    self._angles = {}

  def format(self):
    if (self.altloc == ""):
      return "%s%2s%s" % (self.resname, self.chain_id, self.resid)
    else :
      return "%s%2s%s (conformer %s)" % (self.resname, self.chain_id,
        self.resid, self.altloc)

  def format_csv(self, include_map_label=True):
    if (self.altloc == ""):
      prefix = "%s%2s%s," % (self.resname, self.chain_id, self.resid)
    else :
      prefix = "%s%2s%s %s," % (self.resname, self.chain_id, self.resid,
        self.altloc)
    lines = []
    for i in range(1, self.n_chi+1):
      chi = self.get_angle(i)
      if (chi is not None):
        if include_map_label :
          lines.append(prefix + "2mFo-DFc," + chi.format_csv())
        else :
          lines.append(prefix + chi.format_csv())
        if (chi.fofc_densities is not None):
          if include_map_label :
            lines.append(prefix + "mFo-DFc," + chi.format_csv(fofc=True))
          else :
            lines.append(prefix + chi.format_csv(fofc=True))
    return "\n".join(lines)

  def add_angle(self, **kwds):
    chi = ringer_chi(**kwds)
    self._angles[chi.id] = chi

  def get_angle(self, id):
    return self._angles.get(id, None)

def sample_angle(
    i_seqs,
    sites_cart,
    map_coeffs,
    real_map,
    difference_map,
    sigma,
    angle_start,
    params,
    sampling_method="linear",
    unit_cell=None):
  """
  Given a set of four sites defining a rotatable dihedral angle, sample the
  density at the fourth site in small angular increments.

  returns: a tuple of lists containing the sampled density values (floats) for
           the primary map and optional difference map.
  """
  frac_matrix = None
  if (unit_cell is None):
    assert (map_coeffs is not None)
    unit_cell = map_coeffs.unit_cell()
  frac_matrix = unit_cell.fractionalization_matrix()
  assert (sampling_method != "direct") or (map_coeffs is not None)
  from cctbx import maptbx
  from scitbx.matrix import rotate_point_around_axis
  point = rotate_point_around_axis(
    axis_point_1=sites_cart[1],
    axis_point_2=sites_cart[2],
    point=sites_cart[3],
    angle=-angle_start,
    deg=True)
  # TODO: present option to have point (sites_cart[3]) be generated based on
  # idealized geometry.
  n_degrees = 0
  densities = []
  difference_densities = []
  while (n_degrees < 360):
    point = rotate_point_around_axis(
      axis_point_1=sites_cart[1],
      axis_point_2=sites_cart[2],
      point=point,
      angle=params.sampling_angle,
      deg=True)
    point_frac = unit_cell.fractionalize(site_cart=point)
    rho = rho_fofc = None
    if (sampling_method == "spline") and (map_coeffs is not None):
      rho = real_map.tricubic_interpolation(point_frac)
      if (difference_map is not None):
        rho_fofc = difference_map.tricubic_interpolation(point_frac)
    elif (sampling_method == "linear") or (map_coeffs is None):
      if (map_coeffs is None):
        rho = maptbx.non_crystallographic_eight_point_interpolation(
          map=real_map,
          gridding_matrix=frac_matrix,
          site_cart=point)
          #allow_out_of_bounds=True)
      else :
        rho = real_map.eight_point_interpolation(point_frac)
        if (difference_map is not None):
          rho_fofc = difference_map.eight_point_interpolation(point_frac)
    else :
      rho = map_coeffs.direct_summation_at_point(
        site_frac=point_frac,
        sigma=sigma).real
    densities.append(rho)
    if (rho_fofc is not None):
      difference_densities.append(rho_fofc)
    n_degrees += params.sampling_angle
  #print densities
  return densities, difference_densities

class iterate_over_residues(object):
  """
  Given a PDB hierarchy and electron density, run Ringer analysis for all
  applicable amino acid residues in the model.  Defaults to examining all chi
  angles but this can be overriden.  Implemented as a class to facilitate
  parallelization, but the instantiated object can be discarded after the
  'results' attribute is retrieved.
  """
  def __init__(self,
                pdb_hierarchy,
                params,
                map_coeffs=None,
                difference_map_coeffs=None,
                map_data=None,
                unit_cell=None,
                grid_spacing=0.2,
                sampling_method="linear",
                n_chi_max=4,
                log=None):
    if (log is None) : log = sys.stdout
    adopt_init_args(self, locals())
    models = pdb_hierarchy.models()
    if (len(models) > 1):
      raise Sorry("Multi-model PDB files not supported.")
    self.sigma = self.real_map = self.difference_map = None
    if (map_coeffs is not None):
      self.unit_cell = map_coeffs.unit_cell()
      if (params.sampling_method == "direct"):
        self.map_coeffs = self.map_coeffs.expand_to_p1()
        if (not map_coeffs.anomalous_flag()):
          self.map_coeffs = self.map_coeffs.generate_bijvoet_mates()
      if (sampling_method != "direct") or (params.scaling == "sigma"):
        fft_map = self.map_coeffs.fft_map(resolution_factor=grid_spacing)
        if (params.scaling == "sigma"):
          self.sigma = fft_map.statistics().sigma()
          fft_map.apply_sigma_scaling()
        else :
          fft_map.apply_volume_scaling()
        self.real_map = fft_map.real_map_unpadded()
    else :
      assert (map_data is not None)
      self.unit_cell = unit_cell
      self.real_map = map_data
#      space_group_number = ccp4_map.space_group_number
#      from cctbx import crystal
#      crystal_symmetry_map = crystal.symmetry(ccp4_map.unit_cell().parameters(), space_group_number)
#      if not crystal_symmetry_model:
#        print >> self.log, """Warning: the model does not contain symmetry information. Using map information."""
#      elif not crystal_symmetry_map.is_similar_symmetry(crystal_symmetry_model):
#        print >> self.log, """Warning: The map and model appear to have different crystal symmetry information.
#          EMRinger will assume the map symmetry data is correct and process."""
#      # If map space group is P1, then check that model space group is also either not present or is P1.
#      # If both are p1 or model symmetry is not present, then do the shift. Otherwise, no shift.
#      if space_group_number == 1 and not (crystal_symmetry_model and crystal_symmetry_model.space_group() and crystal_symmetry_model.space_group_number() != 1):
#        import mmtbx.utils
#        shift_manager = mmtbx.utils.shift_origin(
#        map_data = ccp4_map.data.as_double(),
#        pdb_hierarchy = pdb_hierarchy,
#        crystal_symmetry = crystal_symmetry_map)
#        if not shift_manager.shift_cart == None:
#          print >> self.log, "Warning: Model and Map use different origin. Applying origin shift to compensate."
#        pdb_hierarchy = shift_manager.pdb_hierarchy # gives you shifted model
#
#        self.real_map = shift_manager.map_data # gives you shifted map
#      else:
#        print >> self.log, """Warning: Structure is not P1, so automatic origin shifts cannot currently be applied"""
#        self.real_map = ccp4_map.data.as_double()
      # XXX assume that the map is already scaled properly (in the original
      # unit cell)
      #models = pdb_hierarchy.models()
      self.sigma = 1 #ccp4_map.statistics().sigma()
      # XXX the unit cell that we need for the non-crystallographic
      # interpolation is not what comes out of the map - it's the
      #self.unit_cell = ccp4_map.grid_unit_cell()
    if (difference_map_coeffs is not None):
      if (sampling_method == "direct"):
        self.difference_map_coeffs = self.difference_map_coeffs.expand_to_p1()
        if (not difference_map_coeffs.anomalous_flag()):
          self.difference_map_coeffs = \
            self.difference_map_coeffs.generate_bijvoet_mates()
      if (sampling_method != "direct") or (params.scaling == "sigma"):
        fft_map = self.difference_map_coeffs.fft_map(
          resolution_factor=params.grid_spacing)
        if (params.scaling == "sigma"):
          fft_map.apply_sigma_scaling()
        else :
          fft_map.apply_volume_scaling()
        self.difference_map = fft_map.real_map_unpadded()
    results = []
    from mmtbx.rotamer import sidechain_angles
    self.angle_lookup = sidechain_angles.SidechainAngles(False)
    self.sites_cart = pdb_hierarchy.atoms().extract_xyz()
    self.residue_groups = []
    for chain in models[0].chains():
      self.residue_groups.extend(chain.residue_groups())
    if (params.nproc in [None,Auto]) or (params.nproc > 1):
      # this will be a list of lists
      results_ = easy_mp.pool_map(
        processes=params.nproc,
        fixed_func=self.__sample_density,
        args=list(range(len(self.residue_groups))))
      # now flatten it out
      self.results = []
      for result_list in results_ : self.results.extend(result_list)
    else :
      self.results = []
      for i_res in range(len(self.residue_groups)):
        self.results.extend(self.__sample_density(i_res, verbose=True))
    if len(self.results) == 0:
      raise Sorry("""No residues could be scanned by EMRinger, so scores cannot be generated.
      There are a few problems that can lead to this, including not having
      modeled side chains (poly-A or poly-G models), mismatches between the map
      and model grid, or corrupted map density values. These problems can often
      be assessed with molecular graphics tools such as pymol or coot.""")

  def __sample_density(self, i_res, verbose=False):
    import iotbx.pdb
    get_class = iotbx.pdb.common_residue_names_get_class
    residue_group = self.residue_groups[i_res]
    conformers = residue_group.conformers()
    results = []
    for i_conf, conformer in enumerate(residue_group.conformers()):
      if (i_conf > 0) and (self.params.skip_alt_confs):
        continue
      residue = conformer.only_residue()
      if (get_class(residue.resname) == "common_amino_acid"):
        n_chi = int(self.angle_lookup.chisPerAA.get(residue.resname.lower(),0))
        if (n_chi == 0) : continue
        xyz = None
        for atom in residue.atoms():
          if (atom.name.strip() == "CA"):
            xyz = atom.xyz
            break
        res_out = ringer_residue(
          resname=residue.resname,
          chain_id=residue_group.parent().id,
          # resid=residue.resid(),
          resid=residue.resseq_as_int(),
          altloc=conformer.altloc,
          n_chi=n_chi,
          xyz=xyz)
        if (verbose):
          print("  %s:" % residue.id_str(), file=self.log)
        for i in range(1, min(self.n_chi_max+1, n_chi+1)):
          try :
            atoms = self.angle_lookup.extract_chi_atoms("chi%d" % i, residue)
          except AttributeError as e :
            print("Warning: Could not load chi {} atoms".format(i), file=self.log)
            pass
          else :
            try :
              if (atoms is None):
                print("Warning: No side chain atoms detected in model", file=self.log)
                break
              i_seqs = [ atom.i_seq for atom in atoms ]
              sites_chi = [ self.sites_cart[i_seq] for i_seq in i_seqs ]
              from cctbx.geometry_restraints import dihedral
              chi = dihedral(
                sites=sites_chi,
                angle_ideal=0,
                weight=0)
              if (verbose):
                print("    chi%d = %.1f" % (i, chi.angle_model), file=self.log)
              densities, fofc_densities = sample_angle(
                i_seqs=i_seqs,
                sites_cart=sites_chi,
                map_coeffs=self.map_coeffs,
                real_map=self.real_map,
                difference_map=self.difference_map,
                unit_cell=self.unit_cell,
                angle_start=chi.angle_model,
                sigma=self.sigma,
                params=self.params,
                sampling_method=self.sampling_method)
              if (len(fofc_densities) == 0):
                fofc_densities = None
              else :
                assert (len(fofc_densities) == len(densities))
              if (verbose) : pass
              res_out.add_angle(
                id=i,
                angle_current=chi.angle_model,
                densities=densities,
                fofc_densities=fofc_densities,
                sampling=self.params.sampling_angle)
            except Exception as e :
              pass
        if not len(res_out._angles) == 0:
          results.append(res_out)
    return results

class Peak(object):
  """
  Container for information about the sampling angle where density is at the
  global maximum for the given Chi angle.  Used for EMRinger.
  """
  # The peak object, should eventually get moved into ringer I suspect.
  def __init__(self, resname, resid, chain_id, n_chi, chi_value, rho_value):
    adopt_init_args(self, locals())
    self.chi_value=chi_value%360

  def __repr__(self):
    return "\n%s\t%s\t%s\t%s\t%d\t%f" % (self.resname,self.resid,self.chain_id,self.n_chi,self.chi_value*5,self.rho_value)

class Peaklist(object):
  # Right now this is just a slightly specialized list. I may add functionality
  # later, however.
  def __init__(self):
    self.peaks=[]

  def sorted(self, *key):
    return sorted(self.peaks,*key)

  def append_lists(self,other_peaklist):
    self.peaks = self.peaks+ other_peaklist.peaks

  def add_new(self,resname, resid, chain_id, n_chi, chi_value, rho_value):
    self.peaks.append(Peak(resname, resid, chain_id, n_chi, chi_value, rho_value))

  def get_peaks(self):
    return self.peaks

  def __len__(self):
    return len(self.peaks)

  def __repr__(self):
    return str(sorted(self.peaks,key=lambda peak: peak.chi_value))
