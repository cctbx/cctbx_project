from __future__ import absolute_import, division, print_function
from scitbx.source_generators import utils
import sys

this = "scitbx.source_generators.array_family.generate_operator_traits_builtin"

# Signed types only, to avoid the pitfalls of signed/unsigned conversions.
types_ordered = (
  "signed char",
  "short",
  "int",
  "long",
  "long long",
  "float",
  "double",
  "std::complex<float>",
  "std::complex<double>",
)

unsigned = (
  "unsigned char",
  "unsigned short",
  "unsigned int",
  "unsigned long",
  "unsigned long long"
)

floating = ("float", "double")

special_pairs = {
  ("double", "std::complex<float>"): "std::complex<double>",
  ("std::complex<float>", "double"): "std::complex<double>",
}

def build_pairs():
  op_types = []
  result_type = []
  for i, itype in enumerate(types_ordered):
    for j, jtype in enumerate(types_ordered):
      op_types.append((itype, jtype))
      if (i >= j):
        result_type.append(0)
      else:
        result_type.append(jtype)
  for unsigned_t in unsigned:
    for floating_t in floating:
      op_types.append((unsigned_t, floating_t))
      result_type.append(floating_t)
      op_types.append((floating_t, unsigned_t))
      result_type.append(floating_t)
  for op_t in special_pairs:
    result_type[op_types.index(op_t)] = special_pairs[op_t]
  return op_types, result_type

def run(target_dir):
  op_types, result_type = build_pairs()
  assert len(op_types) == len(result_type)
  if ("--Raw" in sys.argv):
    for i, optype in enumerate(op_types):
      print("%s + %s = %s" % (
        optype[0], optype[1], result_type[i]), file=f)
  else:
    f = utils.join_open(target_dir, "operator_traits_builtin.h", "w")
    utils.write_this_is_auto_generated(f, this)
    print("""\
#ifndef SCITBX_ARRAY_FAMILY_OPERATOR_TRAITS_BUILTIN_H
#define SCITBX_ARRAY_FAMILY_OPERATOR_TRAITS_BUILTIN_H

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#include <complex>

namespace scitbx { namespace af {

  // The default traits: the result type is the type of the lhs argument.
  template<typename TypeLHS, typename TypeRHS>
  struct binary_operator_traits {
    typedef TypeLHS arithmetic;
  };

  // The remainder of this file defines the traits where the
  // result type is the type of the rhs argument.
""", file=f)

    for i, optype in enumerate(op_types):
      if (result_type[i]):
        print("""  template<>
  struct binary_operator_traits<%s, %s > {
    typedef %s arithmetic;
  };
""" % (optype[0], optype[1], result_type[i]), file=f)

    print("}} // namespace scitbx::af", file=f)
    print(file=f)
    print("#endif // DOXYGEN_SHOULD_SKIP_THIS", file=f)
    print(file=f)
    print("#endif // SCITBX_ARRAY_FAMILY_OPERATOR_TRAITS_BUILTIN_H", file=f)
    f.close()

if (__name__ == "__main__"):
  run(".")
