import cctbx.eltbx.caasf

from scitbx.python_utils import misc
ext = misc.import_ext("cctbx_boost.xray_ext")
misc.import_regular_symbols(globals(), ext.__dict__)
del misc

from cctbx.xray import structure_factors
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
import types
import sys

def scatterer(label="",
              site=(0,0,0),
              u=None,
              occupancy=1,
              caasf="",
              fp_fdp=0j,
              b=None):
  """\
Python wrapper for C++ constructor.
"""
  assert u is None or b is None
  if   (b is not None): u = adptbx.b_as_u(b)
  elif (u is None): u = 0
  if (type(caasf) == type("")):
    if (caasf == ""):
      caasf = wk1995(label, 0)
    else:
      caasf = wk1995(caasf, 1)
  return ext.scatterer(label, site, u, occupancy, caasf, fp_fdp)

def _scatterer_copy(self,
                    label=None,
                    site=None,
                    u=None,
                    b=None,
                    occupancy=None,
                    caasf=None,
                    fp_fdp=None):
  assert u is None or b is None
  if (b is not None): u = adptbx.b_as_u(b)
  if (label is None): label = self.label
  if (site is None): site = self.site
  if (u is None):
    if (self.anisotropic_flag): u = self.u_star
    else: u = self.u_iso
  if (occupancy is None): occupancy = self.occupancy
  if (caasf is None): caasf = self.caasf
  if (fp_fdp is None): fp_fdp = self.fp_fdp
  return scatterer(
    label=label,
    site=site,
    u=u,
    occupancy=occupancy,
    caasf=caasf,
    fp_fdp=fp_fdp)

ext.scatterer.copy = _scatterer_copy

class structure(crystal.special_position_settings):

  def __init__(self, special_position_settings, scatterers=None):
    crystal.special_position_settings._copy_constructor(
      self, special_position_settings)
    if (scatterers is None):
      self._scatterers = flex.xray_scatterer()
      self._special_position_indices = flex.size_t()
    else:
      self._scatterers = scatterers.deep_copy()
      self.all_apply_symmetry()

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

  def __getitem__(self, slice_object):
    assert type(slice_object) == types.SliceType
    assert self.scatterers() is not None
    return structure(
      special_position_settings=self,
      scatterers=self.scatterers().__getitem__(slice_object))

  def add_scatterer(self, scatterer):
    i = self.scatterers().size()
    self._scatterers.append(scatterer)
    self.apply_symmetry(i)

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

  def add_scatterers(self, scatterers):
    n = self.scatterers().size()
    special_position_indices = self._all_apply_symmetry(scatterers) + n
    self._scatterers.append(scatterers)
    self._special_position_indices.append(special_position_indices)

  def all_apply_symmetry(self):
    self._special_position_indices =self._all_apply_symmetry(self.scatterers())

  def _all_apply_symmetry(self, scatterers):
    return apply_symmetry(
      self.unit_cell(),
      self.space_group(),
      scatterers,
      self.min_distance_sym_equiv(),
      self.u_star_tolerance(),
      self.assert_is_positive_definite(),
      self.assert_min_distance_sym_equiv())

  def structure_factors(self, anomalous_flag=None, d_min=None,
                              direct=00000, fft=00000):
    miller_set = miller.build_set(self, anomalous_flag, d_min)
    return structure_factors.from_scatterers(
      crystal_symmetry=self,
      d_min=d_min)(
        xray_structure=self,
        miller_set=miller_set,
        direct=direct,
        fft=fft)

  def show_summary(self, f=sys.stdout):
    print >> f, "Number of scatterers:", self.scatterers().size()
    print >> f, "At special positions:", self.special_position_indices().size()
    crystal.symmetry.show_summary(self, f)
    return self

  def show_scatterers(self, f=sys.stdout):
    print >> f, "Label  M  Coordinates            Occ  Uiso or Ucart"
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
      new_structure.add_scatterer(scatterer.copy(site=cb_op(scatterer.site)))
    return new_structure

  def change_hand(self):
    ch_op = self.space_group_info().type().change_of_hand_op()
    return self.change_basis(ch_op)

  def apply_shift(self, shift):
    shifted_scatterers = self.scatterers().deep_copy()
    shifted_scatterers.set_sites(
      shifted_scatterers.extract_sites() + shift)
    return structure(
      special_position_settings=self,
      scatterers=shifted_scatterers)

  def sort(self, by_value="occupancy", reverse=00000):
    assert by_value in ("occupancy",)
    assert reverse in (00000, 0001)
    p = flex.sort_permutation(
      self.scatterers().extract_occupancies(),
      reverse)
    return structure(
      special_position_settings=self,
      scatterers=self.scatterers().shuffle(p))

  def as_emma_model(self):
    from cctbx import euclidean_model_matching as emma
    positions = []
    for scatterer in self.scatterers():
      positions.append(emma.position(scatterer.label, scatterer.site))
    return emma.model(self, positions)

class _target_functor_base:

  def __call__(self, f_calc, compute_derivatives):
    assert f_calc.unit_cell().is_similar_to(
           self.f_obs().unit_cell())
    assert f_calc.space_group() == self.f_obs().space_group()
    if (self.weights() is not None):
      return self._target_calculator(self.f_obs().data(),
                                     self.weights(),
                                     f_calc.data(),
                                     compute_derivatives)
    else:
      return self._target_calculator(self.f_obs().data(),
                                     f_calc.data(),
                                     compute_derivatives)

class _least_squares_residual(_target_functor_base):

  def __init__(self, f_obs, weights=None,
               use_sigmas_as_weights=00000):
    adopt_init_args(self, locals(), hide=0001)
    assert self._weights is None or self._use_sigmas_as_weights == 00000
    self._target_calculator = targets_least_squares_residual
    if (self._use_sigmas_as_weights):
      self._weights = self._f_obs.sigmas().data()

  def f_obs(self):
    return self._f_obs

  def weights(self):
    return self._weights

  def use_sigmas_as_weights(self):
    return self._use_sigmas_as_weights

class _intensity_correlation(_target_functor_base):

  def __init__(self, f_obs, weights=None,
               use_multiplicities_as_weights=00000):
    adopt_init_args(self, locals(), hide=0001)
    assert self._weights is None or self._use_multiplicities_as_weights==00000
    self._target_calculator = targets_intensity_correlation
    if (self._use_multiplicities_as_weights):
      self._weights = self._f_obs.multiplicities().data()

  def f_obs(self):
    return self._f_obs

  def weights(self):
    return self._weights

  def use_multiplicities_as_weights(self):
    return self._use_multiplicities_as_weights

target_functors = dicts.easy(
  least_squares_residual=_least_squares_residual,
  intensity_correlation=_intensity_correlation)
