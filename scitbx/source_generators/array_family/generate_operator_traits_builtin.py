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

class pair(object):

  def __init__(self, lhs, rhs):
    self.lhs = lhs
    self.rhs = rhs
  def __cmp__(self, other):
    if (self.lhs == other.lhs and self.rhs == other.rhs): return 0
    return 1

special_pairs = (
  (pair("double", "std::complex<float>"), "std::complex<double>"),
  (pair("std::complex<float>", "double"), "std::complex<double>"),
)

def build_pairs():
  op_types = []
  result_type = []
  for i in xrange(len(types_ordered)):
    for j in xrange(len(types_ordered)):
      op_types.append(pair(types_ordered[i], types_ordered[j]))
      if (i >= j):
        result_type.append(0)
      else:
        result_type.append(types_ordered[j])
  for unsigned_t in unsigned:
    for floating_t in floating:
      op_types.append(pair(unsigned_t, floating_t))
      result_type.append(floating_t)
      op_types.append(pair(floating_t, unsigned_t))
      result_type.append(floating_t)
  for op_t, result_t in special_pairs:
    result_type[op_types.index(op_t)] = result_t
  return op_types, result_type

def run(target_dir):
  op_types, result_type = build_pairs()
  assert len(op_types) == len(result_type)
  if ("--Raw" in sys.argv):
    for i in xrange(len(op_types)):
      print >> f, "%s + %s = %s" % (
        op_types[i].lhs, op_types[i].rhs, result_type[i])
  else:
    f = utils.join_open(target_dir, "operator_traits_builtin.h", "w")
    utils.write_this_is_auto_generated(f, this)
    print >> f, """\
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
"""

    for i in xrange(len(op_types)):
      if (result_type[i]):
        print >> f, """  template<>
  struct binary_operator_traits<%s, %s > {
    typedef %s arithmetic;
  };
""" % (op_types[i].lhs, op_types[i].rhs, result_type[i])

    print >> f, "}} // namespace scitbx::af"
    print >> f
    print >> f, "#endif // DOXYGEN_SHOULD_SKIP_THIS"
    print >> f
    print >> f, "#endif // SCITBX_ARRAY_FAMILY_OPERATOR_TRAITS_BUILTIN_H"
    f.close()

if (__name__ == "__main__"):
  run(".")
