# $Id$

# Enumeration of type promotions for binary vector operations.
# This approach does not rely on partial specialization.

def write_copyright():
  print """// This is an automatically generated file. Do not edit.
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Dec 2001: Created (R.W. Grosse-Kunstleve)
 */"""

types_ordered = (
# "bool",
# "signed char",
# "short",
  "int",
# "long",
  "float",
  "double",
  "std::complex<float>",
  "std::complex<double>",
)

integer_types = (
  "bool",
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
    write_copyright()
    print """
#ifndef CCTBX_VECTOR_PROMOTIONS_H
#define CCTBX_VECTOR_PROMOTIONS_H

#include <complex>

namespace cctbx { namespace vector {

  // The default promotion: the result type is the type of the lhs argument.
  template<typename TypeLHS, typename TypeRHS>
  struct promotion_trait {
    typedef TypeLHS promotion_type;
  };

}} // namespace cctbx::vector"

// The remainder of this file defines the promotions where the
// result type is the type of the rhs argument.

# define CCTBX_VECTOR_PROMOTIONS_DECLARE(VectorTemplate) \\
\\"""
    for i in xrange(len(op_types)):
      if (result_type[i]):
        print """  template<>\\
  struct promotion_trait<VectorTemplate<%s >, VectorTemplate<%s > > {\\
    typedef VectorTemplate<%s > promotion_type;\\
  };\\
  template<>\\
  struct promotion_trait<VectorTemplate<%s >, %s > {\\
    typedef VectorTemplate<%s > promotion_type;\\
  };\\
  template<>\\
  struct promotion_trait<%s, VectorTemplate<%s > > {\\
    typedef VectorTemplate<%s > promotion_type;\\
  };\\
\\""" % (op_types[i].lhs, op_types[i].rhs, result_type[i],
         op_types[i].lhs, op_types[i].rhs, result_type[i],
         op_types[i].lhs, op_types[i].rhs, result_type[i])

    print ""
    print "#endif // CCTBX_VECTOR_PROMOTIONS_H"

if (__name__ == "__main__"):
  run()
