from __future__ import absolute_import, division, print_function
from scitbx.source_generators.array_family import generate_algebras
from scitbx.source_generators import utils

this = "scitbx.source_generators.array_family.generate_apply"

def one_type(target_dir, array_type_name):
  f = utils.join_open(target_dir, "%s_apply.h" % array_type_name, "w")
  utils.write_this_is_auto_generated(f, this)
  include_array_type_name = array_type_name
  if (array_type_name == "ref"):
    include_array_type_name = "versa"
  generic_include = "functors"
  if (generate_algebras.base_array_type_name(array_type_name) == "tiny"):
    generic_include = "operators"
  print("""\
#ifndef SCITBX_ARRAY_FAMILY_%s_APPLY_H
#define SCITBX_ARRAY_FAMILY_%s_APPLY_H

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#include <scitbx/type_holder.h>
#include <scitbx/array_family/%s.h>
#include <scitbx/array_family/detail/generic_array_%s.h>

namespace scitbx { namespace af {
""" % ((array_type_name.upper(),) * 2 + (
    include_array_type_name, generic_include)), file=f)

  generate_algebras.generate_unary_apply(f, array_type_name)

  print("""}} // namespace scitbx::af

#endif // DOXYGEN_SHOULD_SKIP_THIS

#endif // SCITBX_ARRAY_FAMILY_%s_APPLY_H""" % (array_type_name.upper(),), file=f)
  f.close()

def run(target_dir):
  for array_type_name in (
    "tiny_plain", "tiny",
    "small_plain", "small",
    "shared_plain", "shared",
    "versa_plain", "versa",
    "ref"):
    one_type(target_dir, array_type_name)

if (__name__ == "__main__"):
  run(".")
