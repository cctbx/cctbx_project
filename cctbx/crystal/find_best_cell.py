from __future__ import absolute_import, division, print_function
import functools
from cctbx import sgtbx
from six.moves import range
from six.moves import zip
#from cctbx import crystal

# --------------------------------------------------------------------------------------
# Per-space-group "plan" cache
#
# The set of candidate change-of-basis operations, and *which* of them preserve the input
# space group, depend ONLY on the input space group (its number and its setting) -- never
# on the unit cell. find_best_cell.__init__ is called millions of times (it is the
# dominant cost of dials.stills_process indexing) but across only a few hundred distinct
# space-group settings, so this space-group-level work is constant per setting and is
# wastefully recomputed on every call.
#
# We therefore compute, once per distinct input space group, a "plan": which branch applies
# (the cheap early-return cases, monoclinic, or orthorhombic) plus the list of pre-VALIDATED
# base change-of-basis ops (those that map the space group to itself). The per-call path then
# does only cell-level work -- apply each validated op to the actual unit cell, compare, keep
# the best -- with no candidate generation, no op construction, no space-group filtering, and
# no space_group_type build (the last per-call bottleneck once the rest was hoisted).
#
# The plan is keyed on the input space group OBJECT. crystal.symmetry always tidies its
# group, so it is hashable, and hashing it (~0.5 us) is ~460x cheaper than building
# space_group_info().type() (~230 us) just to obtain a Hall-symbol key. Equal groups hash and
# compare equal regardless of how they were built, so the keying is exact.
#
# Validity of an op is `input_sgi.change_basis(op).group() == input_sg`, which is purely a
# space-group test (cell-independent) -- identical to the per-iteration filter the old code
# ran inline -- so caching it is behavior-preserving. new_denominators() returns a fresh op
# and does not mutate its receiver, so the cached base ops are safe to reuse across calls.
#
# Both caches are functools.lru_cache(maxsize=None): unbounded (the key space is only a few
# hundred space groups) and process-lifetime. Use _best_cell_plan.cache_clear() /
# .cache_info() to reset or inspect.
# --------------------------------------------------------------------------------------

@functools.lru_cache(maxsize=None)
def _standard_info_for(space_group_number):
  standard_info = sgtbx.space_group_info(
    symbol=space_group_number,
    table_id="A1983")
  cb_op_std_ref = standard_info.type().cb_op()
  ends_with_colon_2 = str(standard_info).endswith(":2")
  return (standard_info, cb_op_std_ref, ends_with_colon_2)

@functools.lru_cache(maxsize=None)
def _best_cell_plan(input_sg):
  # Keyed on the input space group object itself: it is hashable (crystal.symmetry always
  # tidies its group) and hashing it costs ~0.5 us, versus ~230 us to build
  # space_group_info().type(). Keying here means that space_group_type build -- the per-call
  # bottleneck once everything else was hoisted -- happens only on a cache miss, i.e. once
  # per distinct space group rather than on every one of the millions of calls.
  input_sgi = sgtbx.space_group_info(group=input_sg)
  space_group_number = input_sgi.type().number()
  if (space_group_number == 1):
    return ("niggli", None, None)
  if (space_group_number < 3 or space_group_number >= 75):
    return ("identity", None, None)
  standard_info, cb_op_std_ref, std_ends_with_colon_2 = _standard_info_for(space_group_number)
  cb_op_inp_ref = input_sgi.type().cb_op()
  cb_op_std_inp = cb_op_inp_ref.inverse() * cb_op_std_ref
  assert standard_info.group().change_basis(cb_op_std_inp) == input_sg
  if (space_group_number <= 15):
    branch = "monoclinic"
    two_fold_info = sgtbx.rot_mx_info(input_sg(1).r())
    assert abs(two_fold_info.type()) == 2
    ev = list(two_fold_info.ev())
    assert ev.count(0) == 2
    unique_axis = ev.index(1)
    candidate_mxs = sgtbx.find_affine(input_sg).cb_mx()
  else:
    branch = "orthorhombic"
    assert not std_ends_with_colon_2
    unique_axis = None
    candidate_mxs = sgtbx.space_group("P 4 3*").change_basis(cb_op_std_inp)
  # Build each base op once (from the xyz string -- the rt_mx matrix overload makes
  # boost.python probe numpy's mat3/vec3 converters) and keep only the ops that preserve
  # the space group. The filter is purely a space-group test, so it is the same selection
  # the old per-iteration `alt_symmetry.space_group() == input_sg` check made.
  valid_base_ops = []
  for mx in candidate_mxs:
    base_op = sgtbx.change_of_basis_op(mx.as_xyz())
    if (input_sgi.change_basis(base_op).group() == input_sg):
      valid_base_ops.append(base_op)
  return (branch, valid_base_ops, unique_axis)

class find_best_cell(object):

  def __init__(self,
        input_symmetry,
        angular_tolerance=None,
        best_monoclinic_beta=True):
    if (angular_tolerance is None):
      angular_tolerance = 3
    # symmetry() and all_cells() are built lazily (see those accessors). The dominant
    # production caller -- crystal.symmetry.change_of_basis_op_to_best_cell -- uses only
    # cb_op(), so the per-call work below avoids constructing any crystal.symmetry object:
    # the best cell is selected by comparing UNIT CELLS, computed as
    # input_sg.average_unit_cell(input_uc.change_basis(cb_op)). That is bit-identical to the
    # unit cell crystal.symmetry.change_basis would produce (its force_compatible_unit_cell
    # path calls the same average_unit_cell), so the selection is unchanged; the full
    # crystal.symmetry objects are only materialized if symmetry()/all_cells() are called.
    self._input_symmetry = input_symmetry
    self._cb_ops = None        # per-iteration cb_ops, in order -> recipe for all_cells()
    self._best_index = None    # index into _cb_ops of the best cell; None => input is best
    self._symmetry = None      # built lazily, or set directly for the early-return cases
    self._all_cells = None     # built lazily, or set directly for the early-return cases
    input_sg = input_symmetry.space_group()
    branch, valid_base_ops, unique_axis = _best_cell_plan(input_sg)
    if (branch == "niggli"):
      self._cb_op = input_symmetry.change_of_basis_op_to_niggli_cell()
      self._symmetry = input_symmetry.change_basis(self._cb_op)
      self._all_cells = [self._symmetry]
      return
    if (branch == "identity"):
      self._cb_op = sgtbx.change_of_basis_op()
      self._symmetry = input_symmetry
      self._all_cells = [input_symmetry]
      return
    input_uc = input_symmetry.unit_cell()
    best_cb_op = sgtbx.change_of_basis_op()
    best_uc = input_uc
    cb_ops = []
    if (branch == "monoclinic"):
      do_best_beta = (best_monoclinic_beta and unique_axis == 1)
      for base_op in valid_base_ops:
        cb_op = base_op.new_denominators(best_cb_op)
        alt_uc = input_sg.average_unit_cell(input_uc.change_basis(cb_op))
        if (do_best_beta):
          cb_op_best_beta = alt_uc.change_of_basis_op_for_best_monoclinic_beta()
          if (not cb_op_best_beta.is_identity_op()):
            cb_op.update(cb_op_best_beta)
            alt_uc = input_sg.average_unit_cell(input_uc.change_basis(cb_op))
        cb_ops.append(cb_op)
        cmp_result = best_uc.compare_monoclinic(
          other=alt_uc,
          unique_axis=unique_axis,
          angular_tolerance=angular_tolerance)
        if (cmp_result > 0):
          best_cb_op = cb_op
          best_uc = alt_uc
          self._best_index = len(cb_ops) - 1
    else:
      for base_op in valid_base_ops:
        cb_op = base_op.new_denominators(best_cb_op)
        alt_uc = input_sg.average_unit_cell(input_uc.change_basis(cb_op))
        cb_ops.append(cb_op)
        cmp_result = best_uc.compare_orthorhombic(alt_uc)
        if (cmp_result > 0):
          best_cb_op = cb_op
          best_uc = alt_uc
          self._best_index = len(cb_ops) - 1
    self._cb_op = best_cb_op
    self._cb_ops = cb_ops

  def cb_op(self):
    return self._cb_op

  def symmetry(self):
    # The best symmetry: input_symmetry itself when no alternative beat it (_best_index is
    # None, best_cb_op is identity), else the best alternative -- shared with the
    # corresponding all_cells() entry so the two return the same object, as before.
    if (self._symmetry is None):
      if (self._best_index is None):
        self._symmetry = self._input_symmetry
      elif (self._all_cells is not None):
        self._symmetry = self._all_cells[self._best_index]
      else:
        self._symmetry = self._input_symmetry.change_basis(self._cb_op)
    return self._symmetry

  def all_cells(self):
    if (self._all_cells is None):
      cells = [self._input_symmetry.change_basis(cb) for cb in self._cb_ops]
      if (self._best_index is not None):
        if (self._symmetry is not None):
          cells[self._best_index] = self._symmetry   # keep identity with symmetry()
        else:
          self._symmetry = cells[self._best_index]
      self._all_cells = cells
    return self._all_cells


# this class is an implementation of an alternative method
# of finding the best cell. The above algorithm only works well for
# point groups. When translations are available in the SG,
# space group changes might occur.
# This class does the following: it determines axis permutiatons
# that do not change the spacegroup (apart from an origin shift maybe)


class alternative_find_best_cell(object):
  def __init__(self,
               unit_cell,
               space_group):
    from cctbx import crystal

    self.unit_cell = unit_cell
    self.xs = crystal.symmetry( self.unit_cell, space_group=space_group)
    self.sg_info = sgtbx.space_group_info(group=space_group)
    self.hall_symbol = self.sg_info.type().hall_symbol()
    self.best_cell = None
    self.best_cb_op = sgtbx.change_of_basis_op( 'x,y,z' )

    # now we have to go to the reference setting
    tmp = self.sg_info.change_of_basis_op_to_reference_setting()
    self.xs = self.xs.change_basis( tmp )
    self.unit_cell = self.xs.unit_cell()
    self.sg_info=self.sg_info.change_basis( tmp )
    self.best_cb_op = tmp*self.best_cb_op


    # note that the order in which the operators are checked is important!
    self.axes_permut = [ sgtbx.change_of_basis_op( 'x,y,z' ),# this will work

                         sgtbx.change_of_basis_op( '-x,z,y' ),# this might for pxxy
                         sgtbx.change_of_basis_op( 'y,x,-z' ),# this might for pxxy
                         sgtbx.change_of_basis_op( 'z,-y,x' ),# this might for pxxy

                         sgtbx.change_of_basis_op( 'z,x,y' ),# this not for pxxy
                         sgtbx.change_of_basis_op( 'y,z,x' ) # this not for pxxy
                         ]

    fix_flags = ['All',0,2,1,None,None]

    self.allowed_cb_ops = [sgtbx.change_of_basis_op( 'x,y,z' )]
    self.allowed_cb_ops_to_ref = [sgtbx.change_of_basis_op( 'x,y,z' )]

    fixed_element = None

    identity_op = sgtbx.change_of_basis_op( 'x,y,z' ).c().r()

    for cb_op, fixed in zip(self.axes_permut[1:],
                            fix_flags[1:]): # leave alone the first op
      sg_new = self.sg_info.change_basis( cb_op )
      sg_new_hall_symbol = sg_new.type().hall_symbol()
      #print self.sg_info, "--(", cb_op.as_xyz(),")-->>", sg_new, ";",
      cp_op_to_ref_rotational_part = \
        sg_new.change_of_basis_op_to_reference_setting().c().r()
      #print cp_op_to_ref_rotational_part.as_xyz()

      if  cp_op_to_ref_rotational_part == identity_op :
        # cb_op leaves sg invariantapart from may an origin shift
        self.allowed_cb_ops.append( cb_op )
        self.allowed_cb_ops_to_ref.append( sg_new.change_of_basis_op_to_reference_setting() )
        fixed_element = fixed

    # we now have allowed cb ops, as well as an indication which axis is fixed
    # we can now very easely generate all p;ossible cells
    # and check for order of cell axes among non fixed cel constants
    # the one that has all constants in proper order, is the things we are interested in
    # cell axis order is only an issue in
    # triclinic space groups
    # monoclinic spacegroups
    # orthorhombic spacegroups
    # for other sg's the order is fixed by the symmetry

    self.unit_cell_array = [ unit_cell ]
    self.order_check_array =  [self.order_check( unit_cell, fixed_element )]

    for cb_op in  self.allowed_cb_ops[1:]:
      tmp_uc = self.xs.change_basis( cb_op ).unit_cell()
      self.unit_cell_array.append( tmp_uc )
      self.order_check_array.append( self.order_check( tmp_uc, fixed_element ) )
    self.find_it()


  def order_check(self, unit_cell, fixed=None):
    prefered_order = False
    abc = [ unit_cell.parameters()[0],
            unit_cell.parameters()[1],
            unit_cell.parameters()[2] ]
    if fixed == None:
      if abc[0] <= abc[1]:
        if abc[1] <= abc[2]:
          prefered_order=True

    else:
      assert fixed >= 0
      assert fixed <= 2
      tmp_abc = []
      for ii in range(3):
        if ii != fixed:
          tmp_abc.append( abc[ ii ] )
      if tmp_abc[0] <= tmp_abc[1]:
        prefered_order = True

    return prefered_order

  def find_it(self):
    # check how many trues we have
    n_true = ( self.order_check_array ).count(True)
    best_index = None
    if n_true == 1: # there is only one solution
      best_index = self.order_check_array.index( True )
    else: # there is more then one possible solution, use the first solution one encounters
      for order, ii in zip( self.order_check_array,
                            range(len(self.order_check_array)) ):
        if order:
          best_index = ii
          break

    # Needed for cases like C 1 2 1
    if best_index == None:
      best_index = 0

    self.best_cb_op = (self.allowed_cb_ops_to_ref[ best_index ]*
                       self.allowed_cb_ops[ best_index ] *
                       self.best_cb_op)
    self.best_xs = self.xs.change_basis( self.allowed_cb_ops[ best_index ] * self.allowed_cb_ops_to_ref[ best_index ] )

    self.best_cell = self.best_xs.unit_cell()

  def return_best_cell(self):
    return self.best_cell

  def return_change_of_basis_op_to_best_cell(self):
    return self.best_cb_op

  def return_best_xs(self):
    return self.best_xs


def exercise_alternative():
  from cctbx import crystal
  cb_op = sgtbx.change_of_basis_op("x,y,z")
  for space_group_number in range(1,231):
    sgi = sgtbx.space_group_info(symbol=space_group_number)
    uc = sgi.any_compatible_unit_cell(volume=1000)
    xs = crystal.symmetry(unit_cell=uc, space_group_info=sgi).change_basis(cb_op)

    abest = alternative_find_best_cell( xs.unit_cell(),
                                        sgi.group() )
    print("Space group : ", sgi)
    print("Unit cell   : ", abest.return_best_cell().parameters())
    print("cb op       : ", abest.return_change_of_basis_op_to_best_cell().as_xyz())
    print("N cb op     : ", len( abest.allowed_cb_ops ))
    print()
    print()



def exercise():
  from cctbx import crystal
  cb_op = sgtbx.change_of_basis_op("y,z,x")
  for space_group_number in [1] + list(range(3,73)):
    sgi = sgtbx.space_group_info(symbol=space_group_number)
    uc = sgi.any_compatible_unit_cell(volume=1000)

    xs = crystal.symmetry(unit_cell=uc, space_group_info=sgi).change_basis(cb_op)

    best = find_best_cell(xs,angular_tolerance=3)
    best.symmetry().show_summary()
    print(best.cb_op().as_xyz())
    print(best.cb_op().c())
    print(len(best.all_cells()))
    print()

if (__name__ == "__main__"):
  exercise_alternative()
  exercise()
