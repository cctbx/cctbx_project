from __future__ import absolute_import, division, print_function
from cctbx.eltbx.neutron import neutron_news_1992_table
from cctbx import miller
from cctbx import crystal
from cctbx import sgtbx
from cctbx import adptbx
from scitbx.array_family import flex
import math
import sys
from six.moves import range

class scatterer(object):

  def __init__(self, label="",
                     site=(0,0,0),
                     u=None,
                     occupancy=1,
                     scattering_info="",
                     b=None):
    assert u is None or b is None
    if   (b is not None): u = adptbx.b_as_u(b)
    elif (u is None): u = 0
    if (type(scattering_info) == type("")):
      if (scattering_info == ""):
        scattering_info = neutron_news_1992_table(label, 0)
      else:
        scattering_info = neutron_news_1992_table(scattering_info, 1)
    self.label = label
    self.site = site
    try:
      self.u_iso = float(u)
      self.anisotropic_flag = False
    except Exception:
      assert len(u) == 6
      self.anisotropic_flag = True
      self.u_star = tuple([float(uij) for uij in u])
    self.occupancy = occupancy
    self.scattering_info = scattering_info
    self._multiplicity = 0
    self._weight = 0

  def multiplicity(self):
    return self._multiplicity

  def weight(self):
    return self._weight

  def update_weight(self, space_group_order_z):
    self._weight = self.occupancy * self._multiplicity \
                 / float(space_group_order_z)

  def apply_symmetry(self, unit_cell,
                           space_group,
                           min_distance_sym_equiv=0.5,
                           u_star_tolerance=0,
                           assert_min_distance_sym_equiv=True):
    site_symmetry = sgtbx.site_symmetry(
      unit_cell,
      space_group,
      self.site,
      min_distance_sym_equiv,
      assert_min_distance_sym_equiv)
    self.site = site_symmetry.exact_site()
    self._multiplicity = site_symmetry.multiplicity()
    self.update_weight(space_group.order_z())
    if (self.anisotropic_flag):
      if (u_star_tolerance > 0.):
        assert site_symmetry.is_compatible_u_star(self.u_star,u_star_tolerance)
      self.u_star = site_symmetry.average_u_star(self.u_star)
    return site_symmetry

class structure(crystal.special_position_settings):

  def __init__(self, special_position_settings, scatterers=None):
    crystal.special_position_settings._copy_constructor(
      self, special_position_settings)
    self._scatterers = []
    self._special_position_indices = []
    if (scatterers is not None):
      self.add_scatterers(scatterers)

  def _copy_constructor(self, other):
    crystal.special_position_settings._copy_constructor(
      self, special_position_settings)
    self._scatterers = other._scatterers
    self._special_position_indices = other._special_position_indices

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
      self.assert_min_distance_sym_equiv())
    if (not site_symmetry.is_point_group_1()):
      self.special_position_indices().append(i)

  def add_scatterer(self, scatterer):
    i = len(self.scatterers())
    self._scatterers.append(scatterer)
    self.apply_symmetry(i)

  def add_scatterers(self, scatterers):
    for i in range(len(scatterers)):
      self.add_scatterer(scatterers[i])

  def structure_factors(self, d_min):
    print("WARNING: RESULTS NOT VERIFIED") # XXX
    miller_set = miller.build_set(
      crystal_symmetry=self,
      anomalous_flag=True, # XXX always True?
      d_min=d_min)
    f_calc = flex.complex_double()
    for h in miller_set.indices():
      fc = 0j
      for scatterer in self.scatterers():
        site_symmetry = self.site_symmetry(scatterer.site)
        equiv_sites = sgtbx.sym_equiv_sites(site_symmetry)
        sum_exp_j_two_pi_hx = 0j
        for i_symop,x in enumerate(equiv_sites.coordinates()):
          sum_hx = 0
          for i in range(3):
            sum_hx += h[i] * x[i]
          phase = 2 * math.pi * sum_hx
          exp_j_two_pi_hx = complex(math.cos(phase), math.sin(phase))
          if (scatterer.anisotropic_flag):
            r = self.space_group()(i_symop).r()
            hr = h
            dw = adptbx.debye_waller_factor_u_star(hr, scatterer.u_star)
            exp_j_two_pi_hx *= dw
          sum_exp_j_two_pi_hx += exp_j_two_pi_hx
        b_j = scatterer.scattering_info.bound_coh_scatt_length()
        fc_site = scatterer.weight() * b_j * sum_exp_j_two_pi_hx
        if (not scatterer.anisotropic_flag):
          d_star_sq = self.unit_cell().d_star_sq(h)
          dw = adptbx.debye_waller_factor_u_iso(d_star_sq/4, scatterer.u_iso)
          fc_site *= dw
        fc += fc_site
      f_calc.append(fc)
    return miller.array(
      miller_set=miller_set,
      data=f_calc)

  def show_summary(self, f=None):
    if (f is None): f = sys.stdout
    print("Number of scatterers:", len(self.scatterers()), file=f)
    print("At special positions:", len(self.special_position_indices()), file=f)
    crystal.symmetry.show_summary(self, f)
    return self

  def show_scatterers(self, f=None):
    if (f is None): f = sys.stdout
    print("Label  M  Coordinates            Occ  Uiso or Ucart", file=f)
    for scatterer in self.scatterers():
      print("%-4s" % (scatterer.label,), end=' ', file=f)
      print("%3d" % (scatterer.multiplicity(),), end=' ', file=f)
      print("%7.4f %7.4f %7.4f" % scatterer.site, end=' ', file=f)
      print("%4.2f" % (scatterer.occupancy,), end=' ', file=f)
      if (not scatterer.anisotropic_flag):
        print("%6.4f" % (scatterer.u_iso,), end=' ', file=f)
      else:
        print(("%6.3f " * 5 + "%6.3f") % adptbx.u_star_as_u_cart(
          self.unit_cell(), scatterer.u_star), end=' ', file=f)
      print(file=f)
    return self
