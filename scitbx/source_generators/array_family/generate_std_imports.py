from scitbx.source_generators.array_family import generate_operator_functors
from scitbx.source_generators import utils

this = "scitbx.source_generators.array_family.generate_std_imports"

cmath_1arg = (
  'acos', 'cos', 'tan',
  'asin', 'cosh', 'tanh',
  'atan', 'exp', 'sin',
  'fabs', 'log', 'sinh',
  'ceil', 'floor', 'log10', 'sqrt',
)

cmath_2arg = (
  'fmod', 'pow', 'atan2',
)

cstdlib_1arg = (
  'abs',
)

algorithm_2arg = (
  'each_min', 'each_max',
)

complex_1arg = (
# "cos",
# "cosh",
# "exp",
# "log",
# "log10",
# "sin",
# "sinh",
# "sqrt",
# "tan",
# "tanh",
  "conj",
)

complex_special = (
("ElementType", "real", "std::complex<ElementType>"),
("ElementType", "imag", "std::complex<ElementType>"),
("ElementType", "abs", "std::complex<ElementType>"),
("ElementType", "arg", "std::complex<ElementType>"),
("ElementType", "norm", "std::complex<ElementType>"),
("std::complex<ElementType>", "pow", "std::complex<ElementType>",
                                     "int"),
("std::complex<ElementType>", "pow", "std::complex<ElementType>",
                                     "ElementType"),
("std::complex<ElementType>", "pow", "std::complex<ElementType>",
                                     "std::complex<ElementType>"),
("std::complex<ElementType>", "pow", "ElementType",
                                     "std::complex<ElementType>"),
("std::complex<ElementType>", "polar", "ElementType",
                                       "ElementType"),
)

complex_special_addl_1arg = ("real", "imag", "arg", "norm")
complex_special_addl_2arg = ("polar",)

def filter_function_name(name):
  if name in ('each_min', 'each_max'):
    return name[-3:]
  return name

def generate_1arg(f):
  for function_name in (
    cmath_1arg + cstdlib_1arg + complex_1arg + complex_special_addl_1arg):
    generate_operator_functors.generate_unary(f,
      function_name, function_name + "(x)")

def generate_2arg(f):
  for function_name in cmath_2arg + algorithm_2arg + complex_special_addl_2arg:
    generate_operator_functors.generate_binary(f,
      function_name, filter_function_name(function_name) + "(x, y)")

def run(target_dir):
  f = utils.join_open(target_dir, "detail/std_imports.h", "w")
  utils.write_this_is_auto_generated(f, this)
  print >> f, """\
#ifndef SCITBX_ARRAY_FAMILY_STD_IMPORTS_H
#define SCITBX_ARRAY_FAMILY_STD_IMPORTS_H

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#include <cmath>
#include <cstdlib>
#include <complex>

namespace scitbx { namespace fn {
"""

  all_function_names = []
  for function_name in (cmath_1arg + cmath_2arg + cstdlib_1arg
                        + algorithm_2arg + complex_1arg):
    if (not function_name in all_function_names):
      all_function_names.append(function_name)
  for entry in complex_special:
    function_name = entry[1]
    if (not function_name in all_function_names):
      all_function_names.append(function_name)

  for function_name in all_function_names:
    print >> f, "  using std::" + filter_function_name(function_name) + ";"

  generate_1arg(f)
  generate_2arg(f)

  print >> f, """
}} // namespace scitbx::af

#endif // DOXYGEN_SHOULD_SKIP_THIS

#endif // SCITBX_ARRAY_FAMILY_STD_IMPORTS_H"""
  f.close()

if (__name__ == "__main__"):
  run(".")
