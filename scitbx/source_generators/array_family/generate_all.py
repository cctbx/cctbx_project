from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from . import generate_reductions
from . import generate_std_imports
from . import generate_operator_functors
from . import generate_operator_traits_builtin
from . import generate_algebras
from . import generate_apply
import os

def refresh(array_family):
  assert os.path.isdir(array_family)
  array_family_detail = os.path.join(array_family, "detail")
  assert os.path.isdir(array_family_detail)
  print('  Generating C++ header files in:\n    "%s"' % array_family)
  generate_reductions.run(array_family)
  generate_std_imports.run(array_family)
  generate_operator_functors.run(array_family)
  generate_operator_traits_builtin.run(array_family)
  generate_algebras.run(array_family)
  generate_apply.run(array_family)

assert __name__ != "__main__"
