from __future__ import absolute_import, division, print_function

import os

import scitbx.source_generators.array_family.generate_reductions
import scitbx.source_generators.array_family.generate_std_imports
import scitbx.source_generators.array_family.generate_operator_functors
import scitbx.source_generators.array_family.generate_operator_traits_builtin
import scitbx.source_generators.array_family.generate_algebras
import scitbx.source_generators.array_family.generate_apply

def refresh(array_family):
  assert os.path.isdir(array_family)
  array_family_detail = os.path.join(array_family, "detail")
  assert os.path.isdir(array_family_detail)
  print('  Generating C++ header files in:\n    "%s"' % array_family)
  scitbx.source_generators.array_family.generate_reductions.run(array_family)
  scitbx.source_generators.array_family.generate_std_imports.run(array_family)
  scitbx.source_generators.array_family.generate_operator_functors.run(array_family)
  scitbx.source_generators.array_family.generate_operator_traits_builtin.run(array_family)
  scitbx.source_generators.array_family.generate_algebras.run(array_family)
  scitbx.source_generators.array_family.generate_apply.run(array_family)

assert __name__ != "__main__"
