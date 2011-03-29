from __future__ import division

import boost.python
ext = boost.python.import_ext("cctbx_symmetry_search_ext")
from cctbx_symmetry_search_ext import ls_with_scale_and_bias

from cctbx import miller
from cctbx import sgtbx
from cctbx.sgtbx import lattice_symmetry
from cctbx import maptbx
from libtbx import adopt_optional_init_args
from scitbx import matrix as mat
import scitbx.math
from scitbx.math import clustering
from cctbx.array_family import flex
import itertools

class symmetrised_shifted_structure_factors(object):

  def __init__(self, miller_set, f_c_in_p1, x, compute_gradient=False):
    assert miller_set.unit_cell().is_similar_to(f_c_in_p1.unit_cell())
    sssft = ext.symmetrised_shifted_structure_factors(
      miller_set.space_group(),
      miller_set.indices(),
      miller.f_calc_map(f_c_in_p1.indices(), f_c_in_p1.data(),
                        miller_set.anomalous_flag()),
      x,
      compute_gradient)
    self.f_x = miller.array(miller_set, sssft.f_x)
    self.grad_f_x = sssft.grad_f_x

  def misfit(self, f_o_or_f_o_sq):
    assert self.f_x.is_similar_symmetry(f_o_or_f_o_sq)
    assert self.f_x.indices().all_eq(f_o_or_f_o_sq.indices())
    assert (f_o_or_f_o_sq.is_xray_amplitude_array()
            or f_o_or_f_o_sq.is_xray_intensity_array())
    if f_o_or_f_o_sq.is_xray_amplitude_array():
      f_o_sq = f_o_or_f_o_sq.f_as_f_sq()
    else:
      f_o_sq = f_o_or_f_o_sq
    return ls_with_scale_and_bias(self.f_x.data(),
                                  self.grad_f_x,
                                  f_o_sq.data(),
                                  self.f_x.multiplicities().data().as_double())


class shift_refinement(object):

  def __init__(self, f_obs, fc_in_p1, initial_shift):
    self.f_obs = f_obs
    self.fc_in_p1 = fc_in_p1
    self.initial_shift = initial_shift
    self.x = flex.double(self.initial_shift)
    self.initial_goos = None
    scitbx.lbfgs.run(self,
                     scitbx.lbfgs.termination_parameters(
                       traditional_convergence_test_eps=0.01))
    self.shift = mat.col(self.x)
    del self.x

  def compute_functional_and_gradients(self):
    sssf = symmetrised_shifted_structure_factors(self.f_obs,
                                                 self.fc_in_p1,
                                                 tuple(self.x),
                                                 compute_gradient=True)
    self.symmetrised_shifted_sf = sssf
    goos = sssf.misfit(self.f_obs)
    if self.initial_goos is None:
      self.initial_goos = goos
    self.goos = goos
    return goos.value, flex.double(goos.gradient)


# the denominator used throughout for space-group symmetries
sg_t_den = 12

class structure_factor_symmetry(object):
  """ Crystallographic symmetry of complex structure factors.

  All attributes featuring crystallographic elements (origin, spacegroup,
  etc) are in a primitive unit cell. """

  grid_resolution_factor = 0.4
  cross_correlation_cutoff_for_centring = 0.75
  phi_sym_acceptance_cutoff = 0.25
  phi_sym_rejection_cutoff = 0.5

  def __init__(self, f_in_p1, **kwds):
    isinstance(f_in_p1, miller.array)
    assert f_in_p1.space_group().type().hall_symbol() == ' P 1'
    self.f_in_p1 = f_in_p1
    adopt_optional_init_args(self, kwds)

    self.search_parameters = maptbx.peak_search_parameters(
      peak_search_level=3,
      interpolate=True,
      min_distance_sym_equiv=0.25,
    )

    self.space_group = sgtbx.space_group('P 1', t_den=sg_t_den)
    self.origin = None
    self.symmetry_pool = []

    self.find_centring_translations()

    if self.space_group.order_z() > 1:
      f_in_centered = self.f_in_p1.customized_copy(
        space_group_info=sgtbx.space_group_info(group=self.space_group)
        ).eliminate_sys_absent().merge_equivalents().array()
      self.cb_op_to_primitive = \
          f_in_centered.change_of_basis_op_to_primitive_setting()
      self.f_in_p1 = f_in_centered.change_basis(self.cb_op_to_primitive)
      self.space_group = self.f_in_p1.space_group()
    else:
      self.cb_op_to_primitive = sgtbx.change_of_basis_op()

    self.f_in_p1 = self.f_in_p1.generate_bijvoet_mates()
    self.find_space_group()
    self.space_group = sgtbx.space_group(self.space_group.type().hall_symbol(),
                                         t_den=sgtbx.sg_t_den)
    self.space_group_info = sgtbx.space_group_info(group=self.space_group)

  def centring_translation_peak_sites(self):
    f = self.f_in_p1
    cc_sf = f * f.conjugate().data() / f.sum_sq()
    cc_map = cc_sf.fft_map(
      symmetry_flags=maptbx.use_space_group_symmetry,
      resolution_factor=self.grid_resolution_factor)

    heights = flex.double()
    sites = []
    cc_map_peaks = cc_map.peak_search(self.search_parameters)
    zero = mat.zeros(3)
    for cc_peak_info in itertools.islice(cc_map_peaks, 5):
      if abs(mat.col(cc_peak_info.site) - zero) < 0.01: continue
      heights.append(cc_peak_info.height)
      sites.append(cc_peak_info.site)

    if not heights: return ()
    clusters = clustering.two_means(heights)
    if clusters.highest_stat < self.cross_correlation_cutoff_for_centring:
      return ()
    return sites[:clusters.cut]

  def find_centring_translations(self):
    for d in self.centring_translation_peak_sites():
      t = [ scitbx.math.continued_fraction.from_real(x, 1e-2).as_rational()
            for x in d ]
      if t.count(0) > 1: continue
      unique_denominators = dict(
        [ (r.denominator(), 1) for r in t if r.numerator() != 0 ]).keys()
      assert len(unique_denominators) in (0, 1)
      if len(unique_denominators) == 1:
        den = unique_denominators[0]
        num = [ r.numerator() for r in t ]
        self.space_group.expand_ltr(
          sgtbx.tr_vec(num, den).new_denominator(sg_t_den))

  def possible_point_group_generators(self):
    lattice_group = lattice_symmetry.group(self.f_in_p1.unit_cell(),
                                           max_delta=1)
    lattice_group.expand_inv(sgtbx.tr_vec((0,0,0)))
    rot_parts = set()
    decorated_rot_parts = []
    for op in lattice_group:
      r = op.r()
      if r.is_unit_mx(): continue
      if op.inverse() in rot_parts: continue
      r_info = sgtbx.rot_mx_info(r)
      if r_info.type() < -2: continue
      rot_parts.add(op)
      decorated_rot_parts.append(
        (r_info.type() == -1, # inversion shall come first,
         list(r_info.ev()).count(0), # axes // to unit cell shall come first
                                     # note Python 2.5- compatibility
         r_info.type() == -2, # mirrors preferred.
         r.order(), # higher order preferred
         op))
    decorated_rot_parts.sort()
    decorated_rot_parts.reverse()
    for item in decorated_rot_parts: yield item[-1]

  def cross_correlation_peaks(self):
    f0 = self.f_in_p1
    for op in self.possible_point_group_generators():
      f, op_times_f = f0.common_sets(
        f0.change_basis(sgtbx.change_of_basis_op(op)),
        assert_is_similar_symmetry=False)
      #assert f.size() == f0.size()
      #XXX a better sanity check is needed here to check the amount of overlap
      #XXX between transformed indices
      cc_sf = f * op_times_f.conjugate().data() / f.sum_sq()
      cc_map = cc_sf.fft_map(
        symmetry_flags=maptbx.use_space_group_symmetry,
        resolution_factor=self.grid_resolution_factor)
      if 0: # display 3D map
        from crys3d.qttbx import map_viewer
        map_viewer.display(window_title=op.as_xyz(),
                           fft_map = cc_map,
                           iso_level_positive_range_fraction=0.8,
                           wires=True,
                           orthographic=True)
      cc_map_peaks = cc_map.peak_search(self.search_parameters)
      peak = cc_map_peaks.next()
      yield (op.r(), mat.col(peak.site))

  def find_space_group(self):
    decorated_symmetry_pool = []
    denominator = 12**3
    for i, (r, d) in enumerate(self.cross_correlation_peaks()):
      t = sgtbx.tr_vec((d*denominator).as_int(), tr_den=denominator)
      cb_op = sgtbx.change_of_basis_op(sgtbx.rt_mx(r, t))
      phi_sym = self.f_in_p1.symmetry_agreement_factor(
        cb_op, assert_is_similar_symmetry=False)
      if phi_sym < self.phi_sym_acceptance_cutoff:
        status = possible_symmetry.accepted
      elif phi_sym < self.phi_sym_rejection_cutoff:
        status = possible_symmetry.unsure
      else:
        status = possible_symmetry.rejected
      decorated_symmetry_pool.append(
        (-status, i, possible_symmetry(r, d, phi_sym, status)))
    decorated_symmetry_pool.sort()
    self.symmetry_pool = [ item[-1] for item in decorated_symmetry_pool ]

    self.origin = mat.mutable_zeros(3)
    for symm in self.symmetry_pool:
      if symm.status != symm.accepted: continue
      symm.set_components_of_global_origin(self.origin)
      if self.origin.elems.count(0) == 0: break
    for symm in self.symmetry_pool:
      if symm.status != symm.accepted: continue
      symm.change_origin(self.origin)
      self.space_group.expand_smx(symm.rt)

  def __str__(self):
    return '\n'.join([ str(e) for e in self.symmetry_pool ])

  def symmetrised_structure_factors(self, x=None, delta=None):
    assert x is None or delta is None
    if delta is not None: x = -self.origin - delta
    elif x is None: x = -self.origin
    f_o = self.f_in_p1\
        .as_amplitude_array()\
        .customized_copy(space_group_info=self.space_group_info)\
        .merge_equivalents().array()
    ssf = symmetrised_shifted_structure_factors(f_o, self.f_in_p1, x)
    gos = ssf.misfit(f_o)
    return gos, ssf.f_x


class possible_symmetry(object):

  accepted, unsure, rejected = xrange(3) # possible values of self.status

  def __init__(self, r, d, symmetry_agreement, status):
    assert r.den() == 1
    self.r = r
    order = r.order()
    self.r_info = sgtbx.rot_mx_info(r)
    type = self.r_info.type()
    axis = self.r_info.ev()

    self.symmetry_agreement = symmetry_agreement
    self.status = status

    # compute intrinsic and location part of d, using p, which is
    # the order times the projector onto r's invariant space
    p = mat.sqr(r.accumulate().as_double())
    t_i_num = (p*d).as_int()
    t_l = d - t_i_num/order
    t_i = sgtbx.tr_vec(sg_t_den//order*t_i_num, tr_den=sg_t_den)

    # compute the origin corresponding to t_l by solving
    # (1 - r) o = t_l
    one_minus_r = -mat.sqr(self.r.minus_unit_mx().num())
    one_minus_r_row_echelon = one_minus_r.as_flex_int_matrix()
    q = mat.identity(3).as_flex_int_matrix()
    rank = scitbx.math.row_echelon_form_t(one_minus_r_row_echelon, q)
    qd = flex.double(mat.sqr(q)*t_l)[:rank]
    o = flex.double((0,0,0))
    scitbx.math.row_echelon_back_substitution_float(
      one_minus_r_row_echelon, qd, o)

    # construct object state
    self.t_i, self.raw_origin = t_i, mat.col(o)
    self.origin = None
    self.one_minus_r = one_minus_r

  def set_components_of_global_origin(self, origin):
    for i in xrange(3):
      if origin[i] == 0 and self.raw_origin[i] != 0:
        origin[i] = self.raw_origin[i]

  def change_origin(self, origin):
    o = self.raw_origin - origin
    t_l = self.one_minus_r*o
    t_l = sgtbx.tr_vec((sg_t_den*t_l).as_int(), tr_den=sg_t_den)
    self.rt = sgtbx.rt_mx(self.r, self.t_i.plus(t_l).new_denominator(sg_t_den))
    tr_info = sgtbx.translation_part_info(self.rt)
    self.origin = tr_info.origin_shift()
    self.t_i = tr_info.intrinsic_part()

  fmt3vec = ','.join(["%.3f"]*3)

  def __str__(self):
    return "[%s] % i |(%s) +(%s) @(%s)<=(%s): %.4f" % (
      {self.accepted: '*', self.unsure: '~', self.rejected: ' '}[self.status],
      self.r.type(), "% i, % i, % i" % self.r_info.ev(), self.t_i, self.origin,
      self.fmt3vec % tuple(self.raw_origin), self.symmetry_agreement)
