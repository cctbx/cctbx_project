import cctbx.eltbx.caasf

from scitbx.python_utils import misc
ext = misc.import_ext("cctbx_boost.xray_ext")
misc.import_regular_symbols(globals(), ext.__dict__)
del misc

from cctbx import crystal
from cctbx import miller
from cctbx import adptbx
from cctbx import maptbx
from cctbx.eltbx.caasf import wk1995
from cctbx.array_family import flex
from cctbx import sgtbx
from cctbx import matrix
from scitbx import fftpack
from scitbx.python_utils import dicts
from scitbx.python_utils.misc import adopt_init_args
import sys

def scatterer(label="",
              site=(0,0,0),
              u=0,
              occupancy=1,
              caasf="",
              fp_fdp=0j):
  """\
Python wrapper for C++ constructor.
"""
  if (type(caasf) == type("")):
    if (caasf == ""):
      caasf = wk1995(label, 0)
    else:
      caasf = wk1995(caasf, 1)
  return ext.scatterer(label, site, u, occupancy, caasf, fp_fdp)

class structure(crystal.special_position_settings):

  def __init__(self, special_position_settings, scatterers=None):
    crystal.special_position_settings._copy_constructor(
      self, special_position_settings)
    self._scatterers = flex.xray_scatterer()
    self._special_position_indices = flex.size_t()
    if (scatterers != None):
      self.add_scatterers(scatterers)

  def _copy_constructor(self, other):
    crystal.special_position_settings._copy_constructor(
      self, special_position_settings)
    self._scatterers = other._scatterers
    self._special_position_indices = other._special_position_indices

  def deep_copy_scatterers(self):
    cp = structure(self)
    cp._scatterers = self._scatterers.deep_copy()
    cp._special_position_indices = self._special_position_indices.deep_copy()
    return cp

  def scatterers(self):
    return self._scatterers

  def special_position_indices(self):
    return self._special_position_indices

  def apply_symmetry(self, i):
    site_symmetry = self._scatterers[i].apply_symmetry(
      self.unit_cell(),
      self.space_group(),
      self.min_distance_sym_equiv(),
      self.u_star_tolerance(),
      self.assert_is_positive_definite(),
      self.assert_min_distance_sym_equiv())
    if (not site_symmetry.is_point_group_1()):
      self.special_position_indices().append(i)

  def all_apply_symmetry(self):
    self._special_position_indices = apply_symmetry(
      self.unit_cell(),
      self.space_group(),
      self.scatterers(),
      self.min_distance_sym_equiv(),
      self.u_star_tolerance(),
      self.assert_is_positive_definite(),
      self.assert_min_distance_sym_equiv())

  def add_scatterer(self, scatterer):
    i = self.scatterers().size()
    self._scatterers.append(scatterer)
    self.apply_symmetry(i)

  def add_scatterers(self, scatterers):
    #XXX crash with RedHat 7.3/gcc 2.96 when running phenix translation search
    #for scatterer in scatterers:
    #  self.add_scatterer(scatterer)
    for i in xrange(scatterers.size()):
      self.add_scatterer(scatterers[i])

  def structure_factors_direct(self, anomalous_flag=None, d_min=None):
    miller_set = miller.build_set(self, anomalous_flag, d_min)
    return structure_factors_direct(self, miller_set)

  def show_summary(self, f=sys.stdout):
    print >> f, "Number of scatterers:", self.scatterers().size()
    print >> f, "At special positions:", self.special_position_indices().size()
    crystal.symmetry.show_summary(self, f)
    return self

  def show_scatterers(self, f=sys.stdout):
    print >> f, "Label  M  Coordinates            Occ  Uiso or Ustar"
    for scatterer in self.scatterers():
      print >> f, "%-4s" % (scatterer.label,),
      print >> f, "%3d" % (scatterer.multiplicity(),),
      print >> f, "%7.4f %7.4f %7.4f" % scatterer.site,
      print >> f, "%4.2f" % (scatterer.occupancy,),
      if (not scatterer.anisotropic_flag):
        print >> f, "%6.4f" % (scatterer.u_iso,),
      else:
        print >> f, ("%6.3f " * 5 + "%6.3f") % adptbx.u_star_as_u_cart(
          self.unit_cell(), scatterer.u_star),
      print >> f
      if (abs(scatterer.fp_fdp) != 0):
        print >> f, "     fp,fdp = %6.4f,%6.4f" % (
          scatterer.fp_fdp.real,
          scatterer.fp_fdp.imag)
    return self

  def apply_special_position_ops_d_target_d_site(self, d_target_d_site):
    for i in self.special_position_indices():
      site_symmetry = self.site_symmetry(self.scatterers()[i].site)
      r = matrix.sqr(site_symmetry.special_op().r().as_double())
      d_target_d_site[i] = (matrix.row(d_target_d_site[i]) * r).elems

  def expand_to_p1(self):
    new_structure = structure(
      crystal.special_position_settings(
        crystal.symmetry.cell_equivalent_p1(self)))
    for scatterer in self.scatterers():
      assert not scatterer.anisotropic_flag, "Not implemented." # XXX
      site_symmetry = self.site_symmetry(scatterer.site)
      equiv_sites = sgtbx.sym_equiv_sites(site_symmetry)
      new_scatterer = scatterer.copy()
      for site in equiv_sites.coordinates():
        new_scatterer.site = site
        new_structure.add_scatterer(new_scatterer)
    return new_structure

  def change_basis(self, cb_op):
    new_structure = structure(
      crystal.special_position_settings.change_basis(self, cb_op))
    for scatterer in self.scatterers():
      assert not scatterer.anisotropic_flag, "Not implemented." # XXX
      new_scatterer = scatterer.copy()
      new_scatterer.site = cb_op(new_scatterer.site)
      new_structure.add_scatterer(new_scatterer)
    return new_structure

  def as_emma_model(self):
    from cctbx import euclidean_model_matching as emma
    positions = []
    for scatterer in self.scatterers():
      positions.append(emma.position(scatterer.label, scatterer.site))
    return emma.model(self, positions)

class structure_factors_direct:

  def __init__(self, xray_structure,
                     miller_set,
                     d_target_d_f_calc=None,
                     d_site_flag=00000,
                     d_u_iso_flag=00000,
                     d_u_star_flag=00000,
                     d_occupancy_flag=00000,
                     d_fp_flag=00000,
                     d_fdp_flag=00000):
    assert xray_structure.unit_cell().is_similar_to(miller_set.unit_cell())
    assert xray_structure.space_group() == miller_set.space_group()
    self._xray_structure = xray_structure
    self._miller_set = miller_set
    self._d_target_d_f_calc = d_target_d_f_calc
    if (d_target_d_f_calc == None):
      d_target_d_f_calc = flex.complex_double()
    self._results = structure_factors_direct_with_first_derivatives(
      self._miller_set.unit_cell(),
      self._miller_set.space_group(),
      self._miller_set.indices(),
      self._xray_structure.scatterers(),
      d_target_d_f_calc,
      d_site_flag,
      d_u_iso_flag,
      d_u_star_flag,
      d_occupancy_flag,
      d_fp_flag,
      d_fdp_flag)

  def xray_structure(self):
    return self._xray_structure

  def miller_set(self):
    return self._miller_set

  def d_target_d_f_calc(self):
    return self._d_target_d_f_calc

  def f_calc_array(self):
    return miller.array(self._miller_set, self._results.f_calc())

  def d_target_d_site(self):
    d_target_d_site = self._results.d_target_d_site()
    xray_structure = self.xray_structure()
    assert d_target_d_site.size() == xray_structure.scatterers().size()
    return d_target_d_site

  def d_target_d_u_iso(self):
    d_target_d_u_iso = self._results.d_target_d_u_iso()
    xray_structure = self.xray_structure()
    assert d_target_d_u_iso.size() == xray_structure.scatterers().size()
    return d_target_d_u_iso

  def d_target_d_u_star(self):
    d_target_d_u_star = self._results.d_target_d_u_star()
    xray_structure = self.xray_structure()
    assert d_target_d_u_star.size() == xray_structure.scatterers().size()
    return d_target_d_u_star

  def d_target_d_occupancy(self):
    d_target_d_occupancy = self._results.d_target_d_occupancy()
    xray_structure = self.xray_structure()
    assert d_target_d_occupancy.size() == xray_structure.scatterers().size()
    return d_target_d_occupancy

  def d_target_d_fp(self):
    d_target_d_fp = self._results.d_target_d_fp()
    xray_structure = self.xray_structure()
    assert d_target_d_fp.size() == xray_structure.scatterers().size()
    return d_target_d_fp

  def d_target_d_fdp(self):
    d_target_d_fdp = self._results.d_target_d_fdp()
    xray_structure = self.xray_structure()
    assert d_target_d_fdp.size() == xray_structure.scatterers().size()
    return d_target_d_fdp

  def d_target_d_site_inplace_frac_as_cart(self, d_target_d_site):
    structure_factors_d_target_d_site_in_place_frac_as_cart(
      self.miller_set().unit_cell(), d_target_d_site)

class structure_factors_fft:

  def __init__(self, xray_structure,
                     miller_set,
                     grid_resolution_factor=1./3,
                     symmetry_flags=maptbx.use_space_group_symmetry,
                     quality_factor=100,
                     wing_cutoff=1.e-3,
                     exp_table_one_over_step_size=-100,
                     max_prime=5):
    assert xray_structure.unit_cell().is_similar_to(miller_set.unit_cell())
    assert xray_structure.space_group() == miller_set.space_group()
    assert symmetry_flags.use_space_group_symmetry()
    self._xray_structure = xray_structure
    self._miller_set = miller_set
    d_min = miller_set.d_min()
    n_real = miller_set.determine_gridding(
      resolution_factor=grid_resolution_factor,
      d_min=d_min,
      symmetry_flags=symmetry_flags,
      max_prime=max_prime)
    rfft = fftpack.real_to_complex_3d(n_real)
    u_extra = calc_u_extra(d_min, grid_resolution_factor, quality_factor)
    force_complex = 00000
    electron_density_must_be_positive = 1
    sampled_density = sampled_model_density(
      xray_structure.unit_cell(),
      xray_structure.scatterers(),
      rfft.n_real(),
      rfft.m_real(),
      u_extra,
      wing_cutoff,
      exp_table_one_over_step_size,
      force_complex,
      electron_density_must_be_positive)
    tags = maptbx.grid_tags(rfft.n_real())
    symmetry_flags = maptbx.symmetry_flags(use_space_group_symmetry=0001)
    tags.build(xray_structure.space_group_info().type(), symmetry_flags)
    sampled_density.apply_symmetry(tags)
    if (not sampled_density.anomalous_flag()):
      map = sampled_density.real_map()
      sf_map = rfft.forward(map)
      collect_conj = 1
    else:
      cfft = fftpack.complex_to_complex_3d(rfft.n_real())
      map = sampled_density.complex_map()
      sf_map = cfft.backward(map)
      collect_conj = 0
    self._f_calc_data = maptbx.structure_factors.from_map(
      sampled_density.anomalous_flag(),
      miller_set.indices(),
      sf_map,
      collect_conj).data()
    sampled_density.eliminate_u_extra_and_normalize(
      miller_set.indices(),
      self._f_calc_data)

  def xray_structure(self):
    return self._xray_structure

  def miller_set(self):
    return self._miller_set

  def f_calc_array(self):
    return miller.array(self.miller_set(), self._f_calc_data)

def structure_factors(xray_structure, miller_set, method=None):
  assert method in (None, "fft", "direct")
  if (method == None):
    approx_number_of_atoms = xray_structure.scatterers().size() \
                           * xray_structure.space_group().order_z()
    if (approx_number_of_atoms > 30):
      method = "fft"
    else:
      method = "direct"
  if (method == "fft"):
    result = structure_factors_fft(xray_structure, miller_set)
  else:
    result = structure_factors_direct(xray_structure, miller_set)
  result.f_calc_method = method
  return result

class _target_functor_base:

  def __call__(self, f_calc_array, compute_derivatives):
    assert f_calc_array.unit_cell().is_similar_to(
           self.f_obs_array().unit_cell())
    assert f_calc_array.space_group() == self.f_obs_array().space_group()
    if (self.weights() != None):
      return self._target_calculator(self.f_obs_array().data(),
                                     self.weights(),
                                     f_calc_array.data(),
                                     compute_derivatives)
    else:
      return self._target_calculator(self.f_obs_array().data(),
                                     f_calc_array.data(),
                                     compute_derivatives)

class _least_squares_residual(_target_functor_base):

  def __init__(self, f_obs_array, weights=None,
               use_sigmas_as_weights=00000):
    adopt_init_args(self, locals(), hide=0001)
    assert self._weights == None or self._use_sigmas_as_weights == 00000
    self._target_calculator = targets_least_squares_residual
    if (self._use_sigmas_as_weights):
      self._weights = self._f_obs_array.sigmas().data()

  def f_obs_array(self):
    return self._f_obs_array

  def weights(self):
    return self._weights

  def use_sigmas_as_weights(self):
    return self._use_sigmas_as_weights

class _intensity_correlation(_target_functor_base):

  def __init__(self, f_obs_array, weights=None,
               use_multiplicities_as_weights=00000):
    adopt_init_args(self, locals(), hide=0001)
    assert self._weights==None or self._use_multiplicities_as_weights==00000
    self._target_calculator = targets_intensity_correlation
    if (self._use_multiplicities_as_weights):
      self._weights = self._f_obs_array.multiplicities().data()

  def f_obs_array(self):
    return self._f_obs_array

  def weights(self):
    return self._weights

  def use_multiplicities_as_weights(self):
    return self._use_multiplicities_as_weights

target_functors = dicts.easy(
  least_squares_residual=_least_squares_residual,
  intensity_correlation=_intensity_correlation)
