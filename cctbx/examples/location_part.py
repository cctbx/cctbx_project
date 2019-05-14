from __future__ import absolute_import, division, print_function
from cctbx import sgtbx

class symmetries_with_nonzero_location_parts(object):

  def browse(self):
    for self.symbol in sgtbx.space_group_symbol_iterator():
      self.sg = sgtbx.space_group(self.symbol.hall()).make_tidy()
      self.z2p_op = self.sg.z2p_op()
      self.sg_p = self.sg.change_basis(self.z2p_op)
      self.on_new_space_group()
      for self.op in self.sg_p:
        self.rot_info = self.op.r().info()
        self.tr_info = sgtbx.translation_part_info(self.op)
        if self.tr_info.origin_shift().is_zero(): continue
        self.on_new_symmetry()

  def on_new_space_group(self):
    self.space_group_printout = (
      "%s (%i) [ %s ]"
      % (self.symbol.hall(), self.symbol.number(), self.z2p_op.as_xyz()))


  def on_new_symmetry(self):
    if self.space_group_printout is not None:
      print()
      print(self.space_group_printout)
      self.space_group_printout = None
    print("\t% i |%s +(%s) @(%s)" % (self.rot_info.type(), self.rot_info.ev(),
                                    self.tr_info.intrinsic_part(),
                                    self.tr_info.origin_shift()))


class a_theorem_in_primitive_settings(symmetries_with_nonzero_location_parts):
  """ For any space group in a primitive setting, if it were to contain
      two elements (R|t) and (R|t+d) where t is the intrinsic part,
      then d is a lattice translation.
  """

  def on_new_symmetry(self):
    r = self.op.r()
    self.n_translations.setdefault(r, 0)
    self.n_translations[r] += 1
    if self.n_translations[r] > 1:
      yield super(a_theorem_in_primitive_settings, self).on_new_symmetry()

  def on_new_space_group(self):
    self.n_translations.clear()
    super(a_theorem_in_primitive_settings, self).on_new_space_group()

  def verify(self):
    print(self.__class__.__doc__)
    print()
    print("Let's try to find a counter-example ...")
    self.n_translations = {}
    self.browse()
    print("The search is over!")


class symmetries_with_both_nonzero_location_and_intrinsic_parts(
  symmetries_with_nonzero_location_parts):

  def on_new_symmetry(self):
    if (not self.tr_info.intrinsic_part().is_zero()
        and not self.tr_info.location_part().is_zero()):
      super(symmetries_with_both_nonzero_location_and_intrinsic_parts,
            self).on_new_symmetry()

def run():
  import sys
  if sys.argv[1] == "browse":
    symmetries_with_nonzero_location_parts().browse()
  elif sys.argv[1] == "verify":
    a_theorem_in_primitive_settings().verify()
  elif sys.argv[1] == "nonzero-location-and-intrinsic":
    symmetries_with_both_nonzero_location_and_intrinsic_parts().browse()

if __name__ == '__main__':
  run()
