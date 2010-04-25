from __future__ import division

import boost.python
ext = boost.python.import_ext("cctbx_symmetry_search_ext")


class goodness_of_symmetry(ext.goodness_of_symmetry):

  def __init__(self, f_o_or_f_o_sq, f_c_in_p1, x, compute_gradient=False):
    assert (f_o_or_f_o_sq.is_xray_amplitude_array()
            or f_o_or_f_o_sq.is_xray_intensity_array())
    if f_o_or_f_o_sq.is_xray_amplitude_array():
      f_o_sq = f_o_or_f_o_sq.f_as_f_sq()
    else:
      f_o_sq = f_o_or_f_o_sq
    super(goodness_of_symmetry, self).__init__(
      f_o_sq.space_group(),
      f_o_sq.indices(),
      f_o_sq.data(),
      miller.f_calc_map(f_c_in_p1.indices(), f_c_in_p1.data(),
                        f_o_sq.anomalous_flag()),
      x,
      compute_gradient)


from cctbx import miller
from cctbx import sgtbx
import sgtbx.cosets
from cctbx.sgtbx import lattice_symmetry
from cctbx import maptbx
from libtbx import adopt_optional_init_args
from libtbx import itertbx
from scitbx import matrix as mat
import scitbx.math
from scitbx.math import clustering
from cctbx.array_family import flex
from libtbx import containers
from libtbx import wingide


class structure_factor_symmetry(object):

  grid_resolution_factor = 1/3
  cross_correlation_cutoff_for_centring = 0.75
  min_completeness_under_symmetry = 0.9
  cross_correlation_symmetry_rejection_cutoff = 0.5
  cross_correlation_symmetry_acceptance_cutoff = 0.75

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
    self.log = []
    self.space_group = sgtbx.space_group()
    self.n_ltr = None

  def find_centring_translations(self):
    f = self.f_in_p1
    cc_sf = f * f.conjugate().data() / f.sum_sq()
    cc_map = cc_sf.fft_map(
      symmetry_flags=maptbx.use_space_group_symmetry,
      resolution_factor=self.grid_resolution_factor)

    heights = flex.double()
    sites = []
    cc_map_peaks = cc_map.peak_search(self.search_parameters)
    zero = mat.zeros(3)
    for cc_peak_info in itertbx.islice(cc_map_peaks, 5):
      if abs(mat.col(cc_peak_info.site) - zero) < 0.01: continue
      heights.append(cc_peak_info.height)
      sites.append(cc_peak_info.site)

    if not heights: return
    clusters = clustering.two_means(heights)
    if clusters.highest_mean < self.cross_correlation_cutoff_for_centring:
      return

    for d in sites[:clusters.cut]:
      t = [ scitbx.math.continued_fraction.from_real(x, 1e-2).as_rational()
            for x in d ]
      unique_denominators = dict(
        [ (r.denominator(), 1) for r in t if r.numerator() != 0 ]).keys()
      assert len(unique_denominators) in (0, 1)
      if len(unique_denominators) == 1:
        den = unique_denominators[0]
        num = [ r.numerator() for r in t ]
        self.space_group.expand_ltr(
          sgtbx.tr_vec(num, den).new_denominator(sgtbx.sg_t_den))

  def tentative_rotation_parts(self):
    self.lattice_group = lattice_symmetry.group(self.f_in_p1.unit_cell(),
                                                max_delta=1)
    self.lattice_group.expand_inv(sgtbx.tr_vec((0,0,0)))
    queue = containers.hashed_queue().push(
      sgtbx.space_group('P 1').make_tidy())
    rot_parts = set()
    while queue:
      low_pg = queue.pull()
      r_den = low_pg.r_den()
      t_den = low_pg.t_den()
      cosets = sgtbx.cosets.left_decomposition_point_groups_only(
        g=self.lattice_group,
        h=low_pg)
      if len(cosets.partitions) == 1: continue
      higher_symmetries = [ coset[0] for coset in cosets.partitions[1:] ]
      while higher_symmetries:
        op = higher_symmetries.pop()
        higher_pg = sgtbx.space_group(low_pg)\
                  .expand_smx(op.new_denominators(r_den, t_den))\
                  .make_tidy()
        if higher_pg not in queue:
          queue.push(higher_pg)
          if op.inverse() in rot_parts: continue
          if op.r().type() < -2: continue
          rot_parts.add(op)

    rot_part_info = [ rot_mx_info(op.r()) for op in rot_parts ]
    rot_part_info.sort(reverse=True)
    return rot_part_info

  def cross_correlation_peaks(self):
    f0 = self.f_in_p1
    for r_info in self.tentative_rotation_parts():
      r = r_info.r
      if r.is_unit_mx(): continue
      op = sgtbx.rt_mx(r)
      f, op_times_f = f0.common_sets(
        f0.change_basis(sgtbx.change_of_basis_op(op)))
      if f.size() / f0.size() < self.min_completeness_under_symmetry:
        self.log.append(
          "%s: too small completeness under symmetry" % op.as_xyz())
        continue
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
      # For each rotation part, the max number of peaks to search is the
      # number of lattice translations. We search a few more for
      # clustering to work
      heights = flex.double()
      sites = []
      for peak in itertbx.islice(cc_map_peaks, 5):
        heights.append(peak.height)
        sites.append(peak.site)
      yield (r_info, sites, heights)

  def find_space_group(self):
    self.symmetry_pool = []
    self.origin = mat.mutable_zeros(3)
    for r_info, sites, heights in self.cross_correlation_peaks():
      symm = possible_symmetry(r_info.r, sites[0], heights[0])
      self.symmetry_pool.append(symm)
      if symm.cc < self.cross_correlation_symmetry_rejection_cutoff:
        symm.quality = possible_symmetry.black
        continue
      if symm.cc > self.cross_correlation_symmetry_acceptance_cutoff:
        symm.quality = possible_symmetry.white
        symm.set_components_of_global_origin(self.origin)
      else:
        symm.quality = possible_symmetry.grey

    self.point_group = sgtbx.space_group()
    for symm in self.symmetry_pool:
      symm.change_origin(self.origin)
      if symm.quality != possible_symmetry.white: continue
      r0 = sgtbx.rt_mx(symm.r)
      if r0 in self.point_group:
        symm.already_in_spacegroup = True
        continue
      self.space_group.expand_smx(
        symm.rt.new_denominators(r_den=1, t_den=sgtbx.sg_t_den))
      self.point_group.expand_smx(r0)
    self.space_group_info = sgtbx.space_group_info(group=self.space_group)

  def __str__(self):
    import cStringIO
    result = cStringIO.StringIO()
    print >>result, "Rejected symmetries"
    print >>result, "-------------------"
    for symm in self.symmetry_pool:
      if symm.quality == possible_symmetry.black:
        print >>result, symm
        print >>result
    print >>result
    print >>result, "Grey symmetries"
    print >>result, "---------------"
    for symm in self.symmetry_pool:
      if symm.quality == possible_symmetry.grey:
        print >>result, symm
        print >>result
    print >>result
    print >>result, "Accepted symmetries generated by previous ones"
    print >>result, "----------------------------------------------"
    for symm in self.symmetry_pool:
      if (symm.quality == possible_symmetry.white
          and symm.already_in_spacegroup):
        print >>result, symm
        print >>result
    print >>result
    print >>result, "New origin: (%.3f, %.3f, %.3f)" % tuple(self.origin)
    print >>result, ("Space group: %s"
           % self.space_group.type().hall_symbol())
    print >>result
    print >>result, "Accepted Generators"
    print >>result, "-------------------"
    for symm in self.symmetry_pool:
      if (symm.quality == possible_symmetry.white
          and not symm.already_in_spacegroup):
        print >>result, symm
        print >>result
    return result.getvalue()


class rot_mx_info(object):

  def __init__(self, r):
    self.r = r
    self.r_info = sgtbx.rot_mx_info(r)
    self.type = self.r_info.type()
    self.axis = self.r_info.ev()
    self.order = r.order()
    self.ordering = (self.type == -1,
                     self.axis.count(0),
                     self.order,
                     self.type == -2 or self.type > 0,
                     self.axis == (0,0,1),)
  def __str__(self):
    return "%s [%s]" % (self.r_info, self.r.as_xyz())

  def __gt__(self, other):
    return self.ordering > other.ordering

  def __eq__(self, other):
    return self.ordering == other.ordering


class possible_symmetry(rot_mx_info):

  mult = int(1e6)

  white, grey, black = xrange(3)

  def __init__(self, r, d, cc):
    assert r.den() == 1
    super(possible_symmetry, self).__init__(r)
    d = mat.col(d)
    self.proj_r_invariant = mat.sqr(r.accumulate().as_double())
    t_i = (self.proj_r_invariant*d).as_int()
    self.t_i = sgtbx.tr_vec(12//r.order()*t_i, tr_den=12)
    self.d = d - t_i/r.order()
    r_minus_one = flex.int(r.minus_unit_mx().num())
    r_minus_one.reshape(flex.grid(3,3))
    q = mat.identity(3).as_flex_int_matrix()
    rank = scitbx.math.row_echelon_form_t(r_minus_one, q)
    qd = flex.double(mat.sqr(q)*self.d)[:rank]
    o = flex.double((0,0,0))
    scitbx.math.row_echelon_back_substitution_float(
      r_minus_one, qd, o)
    o = tuple(o)
    self.shifted_origin = mat.col(o)
    self.origin = None
    self.cc = cc
    self.quality = possible_symmetry.black
    self.already_in_spacegroup = False

  def set_components_of_global_origin(self, origin):
    for i in xrange(3):
      if origin[i] == 0 and self.shifted_origin[i] != 0:
        origin[i] = self.shifted_origin[i]

  def change_origin(self, origin):
    self.origin = self.shifted_origin - origin

  def t_l(self):
    return sgtbx.tr_vec(
      (mat.sqr(self.r.minus_unit_mx().num())*mat.col(self.origin)*12).as_int(),
      tr_den=12)
  t_l = property(t_l)

  def rt(self):
    return sgtbx.rt_mx(self.r, self.t_i.plus(self.t_l))
  rt = property(rt)

  def __str__(self):
    fmt_3vec = ",".join(["% .3f"]*3)
    result = [ super(possible_symmetry, self).__str__(),
               ("\td=(" + ",".join(["% .3f"]*3) + ")") % self.d.elems,
               ("\tt_i=(%s)") % self.t_i, ]
    if self.origin is not None:
      result.append(("\to=(%s)" % fmt_3vec) % tuple(self.origin))
      result.append("\tt_l=(%s)" % self.t_l)
    result.extend([
      ("\to_shifted=(%s)" % fmt_3vec) % self.shifted_origin.elems,
      "\tcc=%.2f" % self.cc ])
    return '\n'.join(result)
