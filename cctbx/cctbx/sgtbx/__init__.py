from cctbx import uctbx

import boost.python
ext = boost.python.import_ext("cctbx_sgtbx_ext")
from cctbx_sgtbx_ext import *

class empty: pass

from cctbx.array_family import flex
from scitbx import matrix
from boost import rational
import random
import sys

class _space_group(boost.python.injector, ext.space_group):

  def expand_smx(self, smx):
    if (isinstance(smx, str)):
      smx = rt_mx(smx)
    self.raw_expand_smx(smx)

class space_group_info:

  __safe_for_unpickling__ = 0001

  def __init__(self, symbol=None, table_id=None, group=None, number=None):
    assert [symbol, group, number].count(None) >= 2
    if (number is not None):
      symbol = str(number)
    if (symbol is None):
      assert table_id is None
      self._group = group
    else:
      assert group is None
      if (table_id is None):
        self._group = space_group(space_group_symbols(symbol))
      else:
        if (isinstance(symbol, int)): symbol = str(symbol)
        self._group = space_group(space_group_symbols(symbol, table_id))
    if (self._group is not None):
      self._group.make_tidy()
    self._space_group_info_cache = empty()

  def _copy_constructor(self, other):
    self._group = other._group
    self._space_group_info_cache = other._space_group_info_cache

  def __getinitargs__(self):
    return (str(self),)

  def __getstate__(self):
    return None

  def __setstate__(self, state):
    pass

  def group(self):
    return self._group

  def type(self, tidy_cb_op=0001, r_den=cb_r_den, t_den=cb_t_den):
    cache = self._space_group_info_cache
    if (not hasattr(cache, "_type")
        or cache._type_parameters != (tidy_cb_op, r_den, t_den)):
      cache._type_parameters = (tidy_cb_op, r_den, t_den)
      cache._type = space_group_type(self._group, tidy_cb_op, r_den, t_den)
    return cache._type

  def reciprocal_space_asu(self):
    cache = self._space_group_info_cache
    if (not hasattr(cache, "_reciprocal_space_asu")):
      cache._reciprocal_space_asu = reciprocal_space_asu(self.type())
    return cache._reciprocal_space_asu

  def direct_space_asu(self):
    from cctbx.sgtbx.direct_space_asu import reference_table
    cache = self._space_group_info_cache
    if (not hasattr(cache, "_direct_space_asu")):
      reference_asu = reference_table.get_asu(self.type().number())
      cache._direct_space_asu = reference_asu.change_basis(
        self.type().cb_op().inverse())
    return cache._direct_space_asu

  def brick(self):
    cache = self._space_group_info_cache
    if (not hasattr(cache, "_brick")):
      cache._brick = brick(self.type())
    return cache._brick

  def wyckoff_table(self):
    cache = self._space_group_info_cache
    if (not hasattr(cache, "_wyckoff_table")):
      cache._wyckoff_table = wyckoff_table(self.type())
    return cache._wyckoff_table

  def structure_seminvariants(self):
    cache = self._space_group_info_cache
    if (not hasattr(cache, "_structure_seminvariants")):
      cache._structure_seminvariants = structure_seminvariants(self._group)
    return cache._structure_seminvariants

  def reference_setting(self):
    return space_group_info(symbol=self.type().number())

  def is_reference_setting(self):
    return self.type().cb_op().is_identity_op()

  def as_reference_setting(self):
    return self.change_basis(self.type().cb_op())

  def change_basis(self, cb_op):
    if (isinstance(cb_op, str)):
      cb_op = change_of_basis_op(cb_op)
    return space_group_info(group=self.group().change_basis(cb_op))

  def change_hand(self):
    return self.change_basis(self.type().change_of_hand_op())

  def primitive_setting(self):
    return self.change_basis(self.group().z2p_op())

  def __str__(self):
    cache = self._space_group_info_cache
    if (not hasattr(cache, "_lookup_symbol")):
      cache._lookup_symbol = self.type().lookup_symbol()
    return cache._lookup_symbol

  def show_summary(self, f=None, prefix="Space group: "):
    if (f is None): f = sys.stdout
    print >> f, "%s%s (No. %d)" % (
      prefix, str(self), self.type().number())

  def any_compatible_unit_cell(self, volume):
    sg_number = self.type().number()
    if   (sg_number <   3):
      params = (1., 1.3, 1.7, 83, 109, 129)
    elif (sg_number <  16):
      params = (1., 1.3, 1.7, 90, 109, 90)
    elif (sg_number <  75):
      params = (1., 1.3, 1.7, 90, 90, 90)
    elif (sg_number < 143):
      params = (1., 1., 1.7, 90, 90, 90)
    elif (sg_number < 195):
      params = (1., 1., 1.7, 90, 90, 120)
    else:
      params = (1., 1., 1., 90, 90, 90)
    unit_cell = self.type().cb_op().inverse().apply(uctbx.unit_cell(params))
    f = (volume / unit_cell.volume())**(1/3.)
    params = list(unit_cell.parameters())
    for i in xrange(3): params[i] *= f
    return uctbx.unit_cell(params)

def row_echelon_back_substitution(rt_mx, v=None, sol=None, indep=None):
  return ext.row_echelon_back_substitution(rt_mx, v, sol, indep)

class _tr_vec(boost.python.injector, tr_vec):

  def as_rational(self):
    return matrix.col(rational.vector(self.num(), self.den()))

class _rot_mx(boost.python.injector, rot_mx):

  def as_rational(self):
    return matrix.sqr(rational.vector(self.num(), self.den()))

class _rt_mx(boost.python.injector, rt_mx):

  def as_rational(self):
    return matrix.rt((self.r().as_rational(), self.t().as_rational()))

class _search_symmetry_flags(boost.python.injector, ext.search_symmetry_flags):

  def show_summary(self, f=None):
    if (f is None): f = sys.stdout
    print >> f, "use_space_group_symmetry:", self.use_space_group_symmetry()
    print >> f, "use_space_group_ltr:", self.use_space_group_ltr()
    print >> f, "use_normalizer_k2l:", self.use_normalizer_k2l()
    print >> f, "use_normalizer_l2n:", self.use_normalizer_l2n()
    print >> f, "use_seminvariants:", self.use_seminvariants()

class _wyckoff_table(boost.python.injector, wyckoff_table):

  def random_site_symmetry(self,
        special_position_settings,
        i_position,
        unit_shift_range=(-5,6),
        tolerance=1.e-8):
    position = self.position(i_position)
    run_away_counter = 0
    while 1:
      run_away_counter += 1
      assert run_away_counter < 1000
      site = position.special_op() * [random.random() for i in xrange(3)]
      if (unit_shift_range is not None):
        site = [x + random.randrange(*unit_shift_range) for x in site]
      site_symmetry = special_position_settings.site_symmetry(site)
      if (site_symmetry.distance_moved() < tolerance):
        assert site_symmetry.multiplicity() == position.multiplicity()
        return site_symmetry
