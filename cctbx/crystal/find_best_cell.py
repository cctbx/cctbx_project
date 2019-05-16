from __future__ import absolute_import, division, print_function
from cctbx import sgtbx
from six.moves import range
from six.moves import zip
#from cctbx import crystal

class find_best_cell(object):

  def __init__(self,
        input_symmetry,
        angular_tolerance=None,
        best_monoclinic_beta=True):
    if (angular_tolerance is None):
      angular_tolerance = 3
    self._all_cells = []
    space_group_number = input_symmetry.space_group_info().type().number()
    if (space_group_number == 1):
      self._cb_op = input_symmetry.change_of_basis_op_to_niggli_cell()
      self._symmetry = input_symmetry.change_basis(self._cb_op)
      self._all_cells.append(self._symmetry)
      return
    if (space_group_number < 3 or space_group_number >= 75):
      self._cb_op = sgtbx.change_of_basis_op()
      self._symmetry = input_symmetry
      self._all_cells.append(self._symmetry)
      return
    standard_info = sgtbx.space_group_info(
      symbol=space_group_number,
      table_id="A1983")
    cb_op_inp_ref = input_symmetry.space_group_info().type().cb_op()
    cb_op_std_ref = standard_info.type().cb_op()
    cb_op_std_inp = cb_op_inp_ref.inverse() * cb_op_std_ref
    assert standard_info.group().change_basis(cb_op_std_inp) == input_symmetry.space_group()
    best_cb_op = sgtbx.change_of_basis_op()
    best_symmetry = input_symmetry
    if (space_group_number <= 15):
      two_fold_info = sgtbx.rot_mx_info(input_symmetry.space_group()(1).r())
      assert abs(two_fold_info.type()) == 2
      ev = list(two_fold_info.ev())
      assert ev.count(0) == 2
      unique_axis = ev.index(1)
      affine = sgtbx.find_affine(input_symmetry.space_group())
      for cb_mx in affine.cb_mx():
        cb_op = sgtbx.change_of_basis_op(cb_mx).new_denominators(best_cb_op)
        alt_symmetry = input_symmetry.change_basis(cb_op)
        if (alt_symmetry.space_group() == input_symmetry.space_group()):
          if (best_monoclinic_beta and unique_axis == 1):
            cb_op_best_beta = alt_symmetry.unit_cell() \
              .change_of_basis_op_for_best_monoclinic_beta()
            if (not cb_op_best_beta.is_identity_op()):
              cb_op.update(cb_op_best_beta)
              alt_symmetry = input_symmetry.change_basis(cb_op)
          self._all_cells.append(alt_symmetry)
          cmp_result = best_symmetry.unit_cell().compare_monoclinic(
            other=alt_symmetry.unit_cell(),
            unique_axis=unique_axis,
            angular_tolerance=angular_tolerance)
          if (cmp_result > 0):
            best_cb_op = cb_op
            best_symmetry = alt_symmetry
    else:
      assert not str(standard_info).endswith(":2")
      affine_group = sgtbx.space_group("P 4 3*").change_basis(
        cb_op_std_inp)
      for affine_s in affine_group:
        cb_op = sgtbx.change_of_basis_op(affine_s) \
          .new_denominators(best_cb_op)
        alt_symmetry = input_symmetry.change_basis(cb_op)
        if (alt_symmetry.space_group() == input_symmetry.space_group()):
          self._all_cells.append(alt_symmetry)
          cmp_result = best_symmetry.unit_cell().compare_orthorhombic(
            alt_symmetry.unit_cell())
          if (cmp_result > 0):
            best_cb_op = cb_op
            best_symmetry = alt_symmetry
    self._cb_op = best_cb_op
    self._symmetry = best_symmetry

  def cb_op(self):
    return self._cb_op

  def symmetry(self):
    return self._symmetry

  def all_cells(self):
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
