from cctbx import uctbx

from scitbx.python_utils import misc
ext = misc.import_ext("cctbx_boost.sgtbx_ext")
misc.import_regular_symbols(globals(), ext.__dict__)
del misc

class empty: pass

from cctbx.array_family import flex
from cctbx import matrix
from scitbx import rational
from scitbx.boost_python_utils import injector
import sys

class space_group_info:

  __safe_for_unpickling__ = 0001

  def __init__(self, symbol=None, table_id=None, group=None):
    if (symbol is None):
      assert table_id is None
      self._group = group
    else:
      assert group is None
      if (table_id is None):
        self._group = space_group(space_group_symbols(symbol))
      else:
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

  def structure_seminvariant(self):
    cache = self._space_group_info_cache
    if (not hasattr(cache, "_structure_seminvariant")):
      cache._structure_seminvariant = structure_seminvariant(self._group)
    return cache._structure_seminvariant

  def reference_setting(self):
    return space_group_info(symbol=self.type().number())

  def is_reference_setting(self):
    return self.type().cb_op().is_identity_op()

  def as_reference_setting(self):
    return self.change_basis(self.type().cb_op())

  def change_basis(self, cb_op):
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

  def show_summary(self, f=sys.stdout, prefix="Space group: "):
    print >> f, "%s%s (No. %d)" % (
      prefix, str(self), self.type().number())

  def any_compatible_unit_cell(self, volume):
    sg_number = self.type().number()
    if   (sg_number <   3):
      params = (1., 1.3, 1.7, 83, 109, 129)
    elif (sg_number <  15):
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

class _tr_vec(injector, tr_vec):

  def as_rational(self):
    return matrix.col(rational.vector(self.num(), self.den()))

class _rot_mx(injector, rot_mx):

  def as_rational(self):
    return matrix.sqr(rational.vector(self.num(), self.den()))

class _rt_mx(injector, rt_mx):

  def as_rational(self):
    return matrix.rt((self.r().as_rational(), self.t().as_rational()))
