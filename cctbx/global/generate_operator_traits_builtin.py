def write_copyright():
  print """// This is an automatically generated file. Do not edit.
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jan 2002: Created, based on generate_vector_algebra_traits.py (rwgk)
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
    f = open("operator_traits_builtin.h", "w")
    sys.stdout = f
    write_copyright()
    print """
#ifndef CCTBX_ARRAY_FAMILY_OPERATOR_TRAITS_BUILTIN_H
#define CCTBX_ARRAY_FAMILY_OPERATOR_TRAITS_BUILTIN_H

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#include <complex>

namespace cctbx { namespace af {

  template<typename T>
  struct unary_operator_traits{
    typedef T arithmetic;
    typedef bool logical;
  };

  // The default traits: the result type is the type of the lhs argument.
  template<typename TypeLHS, typename TypeRHS>
  struct binary_operator_traits {
    typedef TypeLHS arithmetic;
    typedef bool logical;
    typedef bool boolean;
  };

  // The remainder of this file defines the traits where the
  // result type is the type of the rhs argument.
"""

    for i in xrange(len(op_types)):
      if (result_type[i]):
        print """  template<>
  struct binary_operator_traits<%s, %s > {
    typedef %s arithmetic;
    typedef bool logical;
    typedef bool boolean;
  };
""" % (op_types[i].lhs, op_types[i].rhs, result_type[i])

    print "}} // namespace cctbx::af"
    print ""
    print "#endif // DOXYGEN_SHOULD_SKIP_THIS"
    print ""
    print "#endif // CCTBX_ARRAY_FAMILY_OPERATOR_TRAITS_BUILTIN_H"
    sys.stdout = sys.__stdout__
    f.close()

if (__name__ == "__main__"):
  run()
