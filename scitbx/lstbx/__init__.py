import boost.python
ext = boost.python.import_ext("scitbx_lstbx_ext")
from scitbx_lstbx_ext import *

class normal_equations_extension(boost.python.injector, normal_equations):

  def __iter__(self):
    yield self.normal_matrix_packed_u
    yield self.right_hand_side
