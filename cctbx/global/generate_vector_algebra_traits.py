def write_copyright():
  print """// This is an automatically generated file. Do not edit.
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Dec 2001: Created (R.W. Grosse-Kunstleve)
 */"""

# Signed types only, to avoid the pitfalls of signed/unsigned conversions.
types_ordered = (
  "signed char",
  "short",
  "int",
  "long",
  "float",
  "double",
  "std::complex<float>",
  "std::complex<double>",
)

integer_types = (
  "signed char",
  "short",
  "int",
  "long",
)

class pair:

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

import sys

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
  for op_t, result_t in special_pairs:
    result_type[op_types.index(op_t)] = result_t
  return op_types, result_type

def run():
  op_types, result_type = build_pairs()
  assert len(op_types) == len(result_type)
  if ("--Raw" in sys.argv):
    for i in xrange(len(op_types)):
      print "%s + %s = %s" % (op_types[i].lhs, op_types[i].rhs, result_type[i])
  else:
    f = open("algebra_traits.h", "w")
    sys.stdout = f
    write_copyright()
    print """
#ifndef CCTBX_VECTOR_ALGEBRA_TRAITS_H
#define CCTBX_VECTOR_ALGEBRA_TRAITS_H

#include <complex>

namespace cctbx { namespace vector {

  template <typename T>
  struct algebra_constructor : T {
    algebra_constructor(std::size_t sz) : T(sz) {}
  };

  struct qualifier_v_v {}; // lhs = vector type, rhs = vector type
  struct qualifier_v_s {}; // lhs = vector type, rhs = scalar type
  struct qualifier_s_v {}; // etc.
  struct qualifier_s_s {};

  // Traits for scalar operators.

  // The default traits: the result type is the type of the lhs argument.
  template<typename TypeLHS, typename TypeRHS>
  struct algebra_traits {
    typedef TypeLHS value_type;
    typedef TypeLHS promotion_type;
    typedef bool bool_type;
    typedef qualifier_s_s type_qualifier;
  };

  // The remainder of this file defines the traits where the
  // result type is the type of the rhs argument.
"""

    for i in xrange(len(op_types)):
      if (result_type[i]):
        print """  template<>
  struct algebra_traits<%s, %s > {
    typedef %s value_type;
    typedef %s promotion_type;
    typedef bool bool_type;
    typedef qualifier_s_s type_qualifier;
  };
""" % (op_types[i].lhs, op_types[i].rhs, result_type[i], result_type[i])

    print "}} // namespace cctbx::vector"
    print ""
    print "#endif // CCTBX_VECTOR_ALGEBRA_TRAITS_H"
    sys.stdout = sys.__stdout__
    f.close()

if (__name__ == "__main__"):
  run()
