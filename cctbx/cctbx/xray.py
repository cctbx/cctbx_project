from scitbx.python_utils.misc import import_regular_symbols
from cctbx_boost import xray_ext as ext
import_regular_symbols(globals(), ext.__dict__)
del import_regular_symbols

from cctbx import crystal
from cctbx import miller
from cctbx import adptbx
from cctbx.eltbx.caasf import wk1995
from cctbx.array_family import flex
from cctbx import matrix
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
    for scatterer in scatterers:
      self.add_scatterer(scatterer)

  def structure_factors_direct(self, anomalous_flag=None, d_min=None):
    miller_set = miller.build_set(self, anomalous_flag, d_min)
    return structure_factors_direct(self, miller_set)

  def show_summary(self, f=sys.stdout):
    print "Number of scatterers:", self.scatterers().size()
    print "At special positions:", self.special_position_indices().size()
    crystal.symmetry.show_summary(self, f)
    return self

  def show_scatterers(self, f=sys.stdout):
    print "Label  M  Coordinates            Occ  Uiso or Ustar"
    for scatterer in self.scatterers():
      print "%-4s" % (scatterer.label,),
      print "%3d" % (scatterer.multiplicity(),),
      print "%7.4f %7.4f %7.4f" % scatterer.site,
      print "%4.2f" % (scatterer.occupancy,),
      if (not scatterer.anisotropic_flag):
        print "%6.4f" % (scatterer.u_iso,),
      else:
        print ("%6.3f " * 5 + "%6.3f") % adptbx.u_star_as_u_cart(
          self.unit_cell(), scatterer.u_star),
      print
      if (abs(scatterer.fp_fdp) != 0):
        print "     fp,fdp = %6.4f,%6.4f" % (
          scatterer.fp_fdp.real,
          scatterer.fp_fdp.imag)
    return self

  def apply_special_position_ops_d_target_d_site(self, d_target_d_site):
    for i in self.special_position_indices():
      site_symmetry = self.site_symmetry(self.scatterers()[i].site)
      r = matrix.sqr(site_symmetry.special_op().r().as_double())
      d_target_d_site[i] = (matrix.row(d_target_d_site[i]) * r).elems

class structure_factors_direct:

  def __init__(self, xray_structure,
                     miller_set,
                     d_target_d_f_calc=None,
                     d_site_flag=False,
                     d_u_iso_flag=False,
                     d_u_star_flag=False,
                     d_occupancy_flag=False,
                     d_fp_flag=False,
                     d_fdp_flag=False):
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

class target_functor_base:

  def __call__(self, f_calc_array, compute_derivatives):
    assert f_calc_array.unit_cell().is_similar_to(self.f_obs_array.unit_cell())
    assert f_calc_array.space_group() == self.f_obs_array.space_group()
    if (self.weights):
      return self.target_calculator(self.f_obs_array.data(),
                                    self.weights,
                                    f_calc_array.data(),
                                    compute_derivatives)
    else:
      return self.target_calculator(self.f_obs_array.data(),
                                    f_calc_array.data(),
                                    compute_derivatives)

class least_squares_residual(target_functor_base):

  def __init__(self, f_obs_array, weights=None,
               use_sigmas_as_weights=False):
    adopt_init_args(self, locals())
    assert self.weights == None or self.use_sigmas_as_weights == False
    self.target_calculator = targets_least_squares_residual
    if (self.use_sigmas_as_weights):
      self.weights = self.sigmas()

class intensity_correlation(target_functor_base):

  def __init__(self, f_obs_array, weights=None,
               use_multiplicities_as_weights=False):
    adopt_init_args(self, locals())
    assert self.weights == None or self.use_multiplicities_as_weights == False
    self.target_calculator = targets_intensity_correlation
    if (self.use_multiplicities_as_weights):
      self.weights = self.multiplicities()

target_functors = dicts.easy()
target_functors.least_squares_residual = least_squares_residual
target_functors.intensity_correlation = intensity_correlation
