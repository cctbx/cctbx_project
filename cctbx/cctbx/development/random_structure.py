from cctbx import sgtbx
from cctbx import crystal
from cctbx import adptbx
from cctbx import xray
from cctbx import eltbx
from scitbx.python_utils.misc import adopt_init_args
import random

def have_suitable_hetero_distance(existing_sites,
                                  sym_equiv_sites_of_other_site,
                                  min_hetero_distance):
  for existing_site in existing_sites:
    if (sgtbx.min_sym_equiv_distance_info(
          sym_equiv_sites_of_other_site, existing_site).dist()
        < min_hetero_distance):
      return 00000
  return 0001

def random_site(special_position_settings,
                existing_sites,
                min_hetero_distance=1.5,
                general_position_only=00000,
                grid=None,
                t_centre_of_inversion=None,
                max_trials=100):
  for trial in xrange(max_trials):
    if (grid is None):
      site = (random.random(), random.random(), random.random())
    else:
      site = [random.randrange(g) / float(g) for g in grid]
    site_symmetry = special_position_settings.site_symmetry(site)
    if (general_position_only and not site_symmetry.is_point_group_1()):
      continue
    sym_equiv_sites = sgtbx.sym_equiv_sites(site_symmetry)
    if (not have_suitable_hetero_distance(
              existing_sites, sym_equiv_sites, min_hetero_distance)):
      continue
    site = site_symmetry.exact_site()
    if (t_centre_of_inversion is None):
      return site
    site_inv = [-x+t for x,t in zip(site, t_centre_of_inversion)]
    site_symmetry_inv = special_position_settings.site_symmetry(site_inv)
    if (general_position_only and not site_symmetry_inv.is_point_group_1()):
      continue
    sym_equiv_sites_inv = sgtbx.sym_equiv_sites(site_symmetry_inv)
    if (not have_suitable_hetero_distance(
              existing_sites + [site],
              sym_equiv_sites_inv, min_hetero_distance)):
      continue
    return site, site_symmetry_inv.exact_site()
  return None

def random_sites(special_position_settings,
                 existing_sites,
                 n_new,
                 min_hetero_distance=1.5,
                 general_positions_only=00000,
                 grid=None,
                 t_centre_of_inversion=None,
                 max_trials=100,
                 max_back_track=100):
  n_loop = n_new
  if (t_centre_of_inversion is not None):
    assert n_new % 2 == 0
    n_loop /= 2
  for i_back_track in xrange(max_back_track):
    all_sites = existing_sites[:]
    for i_new in xrange(n_loop):
      site = random_site(special_position_settings,
                         all_sites,
                         min_hetero_distance,
                         general_positions_only,
                         grid,
                         t_centre_of_inversion,
                         max_trials)
      if (site is None):
        break
      if (t_centre_of_inversion is None):
        all_sites.append(site)
      else:
        all_sites.extend(site)
    if (len(all_sites) == len(existing_sites) + n_new):
      return all_sites
  raise RuntimeError, "Cannot find sites matching all constraints."

def random_modify_site(special_position_settings, site, gauss_sigma,
                       max_distance=0,
                       vary_z_only=00000,
                       max_trials=100):
  site_symmetry = special_position_settings.site_symmetry(site)
  assert site_symmetry.distance_moved() < 1.e-5
  unit_cell = special_position_settings.unit_cell()
  site_cart = list(unit_cell.orthogonalize(site))
  for trial in xrange(max_trials):
    if (vary_z_only):
      modified_site_cart = site_cart[:2] \
                         + [random.gauss(site_cart[2], gauss_sigma)]
    else:
      modified_site_cart = [random.gauss(x, gauss_sigma) for x in site_cart]
    modified_site = site_symmetry.special_op() \
                  * unit_cell.fractionalize(modified_site_cart)
    if (max_distance > 0):
      distance = unit_cell.distance(site, modified_site)
      if (distance > max_distance): continue
    modified_site_symmetry = special_position_settings.site_symmetry(
      modified_site)
    if (modified_site_symmetry.special_op() != site_symmetry.special_op()):
      continue
    return modified_site
  raise RuntimeError, "Cannot find suitable site."

class xray_structure(xray.structure):

  def __init__(self,
               space_group_info,
               elements=None,
               n_scatterers=None,
               volume_per_atom=50.,
               min_distance=1.5,
               general_positions_only=00000,
               random_f_prime_d_min=0,
               random_f_prime_scale=0.6,
               random_f_double_prime=0,
               random_f_double_prime_scale=0.6,
               random_u_iso=00000,
               random_u_iso_scale=0.3,
               u_iso=0,
               anisotropic_flag=00000,
               random_u_cart_scale=0.3,
               random_occupancy=00000):
    assert elements is None or n_scatterers is None
    assert not (elements is None and n_scatterers is None)
    adopt_init_args(self, locals(), exclude=("space_group_info",))
    if (elements is not None):
      self.n_scatterers = len(elements)
    crystal_symmetry = crystal.symmetry(
      space_group_info.any_compatible_unit_cell(
        self.n_scatterers
        * volume_per_atom
        * space_group_info.group().order_z()),
      space_group_info=space_group_info)
    special_position_settings = crystal.special_position_settings(
      crystal_symmetry,
      min_distance_sym_equiv=min_distance,
      u_star_tolerance=0,
      assert_is_positive_definite=00000,
      assert_min_distance_sym_equiv=0001)
    xray.structure.__init__(self, special_position_settings)
    if (elements is not None):
      self.build_scatterers(elements)

  def build_scatterers(self, elements, grid=None, t_centre_of_inversion=None):
    all_sites = random_sites(
      special_position_settings=self,
      existing_sites=[scatterer.site for scatterer in self.scatterers()],
      n_new=len(elements),
      min_hetero_distance=self.min_distance,
      general_positions_only=self.general_positions_only,
      grid=grid,
      t_centre_of_inversion=t_centre_of_inversion)
    assert len(all_sites) <= self.n_scatterers
    sf_dict = {}
    for element in elements:
      if (not sf_dict.has_key(element)):
        sf_dict[element] = eltbx.caasf.best_approximation(element)
    fp = 0
    fdp = 0
    n_existing = self.scatterers().size()
    i_label = n_existing
    for element,site in zip(elements, all_sites[n_existing:]):
      i_label += 1
      scatterer = xray.scatterer(
        label=element + str(i_label),
        scattering_type=element,
        site=site)
      site_symmetry = scatterer.apply_symmetry(
        self.unit_cell(),
        self.space_group(),
        self.min_distance_sym_equiv())
      if (self.random_f_prime_d_min):
        f0 = sf_dict[element].at_d_star_sq(1./self.random_f_prime_d_min**2)
        assert f0 > 0
        fp = -f0 * random.random() * self.random_f_prime_scale
      if (self.random_f_double_prime):
        f0 = sf_dict[element].at_d_star_sq(0)
        fdp = f0 * random.random() * self.random_f_double_prime_scale
      scatterer.fp = fp
      scatterer.fdp = fdp
      if (not self.anisotropic_flag):
        scatterer.anisotropic_flag = 00000
        u_iso = self.u_iso
        if (not u_iso and self.random_u_iso):
          u_iso = random.random() * self.random_u_iso_scale
        scatterer.u_iso = u_iso
      else:
        scatterer.anisotropic_flag = 0001
        run_away_counter = 0
        while 1:
          run_away_counter += 1
          assert run_away_counter < 100
          u_cart = [random.random()*self.random_u_cart_scale
                    for i in xrange(3)] + [0.,0.,0.]
          u_cart = adptbx.random_rotate_ellipsoid(u_cart)
          scatterer.u_star = site_symmetry.average_u_star(
                               adptbx.u_cart_as_u_star(
                                 self.unit_cell(), u_cart))
          u_cart = adptbx.u_star_as_u_cart(self.unit_cell(), scatterer.u_star)
          eigenvalues = adptbx.eigenvalues(u_cart)
          if (min(eigenvalues) > 0.001):
            break
      if (self.random_occupancy):
        scatterer.occupancy = max(0.1, min(1.0, random.gauss(0.5, 0.2)))
      self.add_scatterer(scatterer)

  def random_modify_site(self, site, gauss_sigma,
                         max_distance=0,
                         vary_z_only=00000,
                         max_trials=100):
    return random_modify_site(
      self, site, gauss_sigma, max_distance, vary_z_only, max_trials)

  def random_modify_u_iso(self, u_iso, gauss_sigma):
    return max(0.1, random.gauss(u_iso, gauss_sigma))

  def random_modify_u_star(self, u_star, gauss_sigma,
                                 max_relative_difference=1./3,
                                 max_trials=100):
    for trial in xrange(max_trials):
      modified_u_star = []
      for i in xrange(len(u_star)):
        u = u_star[i]
        max_diff = u * max_relative_difference
        modified_u = random.gauss(u, gauss_sigma)
        if (modified_u - u > u + max_diff):
          modified_u = u + max_diff
        elif (u - modified_u > u + max_diff):
          modified_u = u - max_diff
        modified_u_star.append(modified_u)
      u_cart = adptbx.u_star_as_u_cart(self.unit_cell(), modified_u_star)
      eigenvalues = adptbx.eigenvalues(u_cart)
      if (min(eigenvalues) > 0.001):
        return modified_u_star
    raise RuntimeError, "Cannot find suitable u_star."

  def random_modify_occupancy(self, occupancy, gauss_sigma):
    return max(0.1, occupancy - abs(random.gauss(0, gauss_sigma)))

  def random_modify_fp(self, fp, gauss_sigma):
    assert fp < 0
    return min(-0.1, random.gauss(fp, gauss_sigma))

  def random_modify_fdp(self, fdp, gauss_sigma):
    assert fdp > 0
    return max(0.1, random.gauss(fdp, gauss_sigma))

  def random_modify_parmeters(self, parameter_name, gauss_sigma=0.1,
                                    vary_z_only=00000):
    modified_structure = self.deep_copy_scatterers()
    for scatterer in modified_structure.scatterers():
      if (parameter_name == "site"):
        scatterer.site = \
          self.random_modify_site(scatterer.site, gauss_sigma,
                                  vary_z_only=vary_z_only)
      elif (parameter_name == "u_iso"):
        scatterer.u_iso = \
          self.random_modify_u_iso(scatterer.u_iso, gauss_sigma)
      elif (parameter_name == "u_star"):
        scatterer.u_star = \
          self.random_modify_u_star(scatterer.u_star, gauss_sigma)
      elif (parameter_name == "occupancy"):
        scatterer.occupancy = \
          self.random_modify_occupancy(scatterer.occupancy, gauss_sigma)
        scatterer.update_weight(self.space_group().order_z())
      elif (parameter_name == "fp"):
        scatterer.fp = self.random_modify_fp(scatterer.fp, gauss_sigma)
      elif (parameter_name == "fdp"):
        scatterer.fdp = self.random_modify_fdp(scatterer.fdp, gauss_sigma)
      else:
        raise RuntimeError
    return modified_structure
