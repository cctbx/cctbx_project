def write_copyright():
  print """// This is an automatically generated file. Do not edit.
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jan 2002: Created, based on generate_af_operators.py (rwgk)
 */"""

arithmetic_unary_ops = ("-")
arithmetic_binary_ops = ("+", "-", "*", "/", "%")
arithmetic_in_place_binary_ops = ("+=", "-=", "*=", "/=", "%=")
logical_unary_ops = ("!")
logical_binary_ops = ("&&", "||")
boolean_ops = ("==", "!=", ">", "<", ">=", "<=")

class empty: pass

def op_vars(type_flags):
  v = empty()
  v.param = ["ValueTypeLhs", "ValueTypeRhs"]
  v.dotv = ["", ""]
  if (type_flags[0]):
    v.param[0] = "flagged_value<ValueTypeLhs>"
    v.dotv[0] = ".v"
    v.have_both_test = "lhs.f"
  if (type_flags[1]):
    v.param[1] = "flagged_value<ValueTypeRhs>"
    v.dotv[1] = ".v"
    v.have_both_test = "rhs.f"
  if (type_flags[0] and type_flags[1]):
    v.have_both_test = "lhs.f && rhs.f"
  return v

def elementwise_binary_op(op_class, op_symbol, type_flags):
  v = op_vars(type_flags)
  print """  template <typename ValueTypeLhs, typename ValueTypeRhs>
  inline
  flagged_value<
    typename binary_operator_traits<
      ValueTypeLhs, ValueTypeRhs>::%s>
  operator%s(
    const %s& lhs,
    const %s& rhs) {
    flagged_value<
      typename binary_operator_traits<
        ValueTypeLhs, ValueTypeRhs>::%s> result;
    if (%s) {
      result.v = lhs%s %s rhs%s;
      result.f = true;
    }
    return result;
  }
""" % (op_class, op_symbol, v.param[0], v.param[1], op_class,
       v.have_both_test, v.dotv[0], op_symbol, v.dotv[1])

def elementwise_inplace_binary_op(op_symbol, type_flags):
  v = op_vars(type_flags)
  if (type_flags == (1,1)):
    action = """    if (lhs.f) {
      if (rhs.f) lhs.v %s rhs.v;
      else lhs.f = false;
    }""" % (op_symbol,)
  else:
    action = "    if (lhs.f) lhs.v %s rhs;" % (op_symbol,)
  print """  template <typename ValueTypeLhs, typename ValueTypeRhs>
  inline
  flagged_value<ValueTypeLhs>&
  operator%s(
    flagged_value<ValueTypeLhs>& lhs,
    const %s& rhs) {
%s
    return lhs;
  }
""" % (op_symbol, v.param[1], action)

def generate_elementwise_binary_op(op_class, op_symbol):
  for type_flags in ((1,1), (1,0), (0,1)):
    elementwise_binary_op(op_class, op_symbol, type_flags)
  if (op_class == "arithmetic"):
    for type_flags in ((1,1), (1,0)):
      elementwise_inplace_binary_op(op_symbol + "=", type_flags)

def run():
  import sys
  f = open("flagged_value_operators.h", "w")
  sys.stdout = f
  write_copyright()
  print """
#ifndef CCTBX_ARRAY_FAMILY_FLAGGED_VALUE_OPERATORS_H
#define CCTBX_ARRAY_FAMILY_FLAGGED_VALUE_OPERATORS_H

#include <cctbx/array_family/operator_traits.h>

namespace cctbx { namespace af {
"""

  #generate_unary_ops(array_type_name)
  for op_symbol in arithmetic_binary_ops:
    generate_elementwise_binary_op("arithmetic", op_symbol)
  for op_symbol in logical_binary_ops:
    generate_elementwise_binary_op("logical", op_symbol)
  for op_symbol in boolean_ops:
    generate_elementwise_binary_op("boolean", op_symbol)

  print """}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_FLAGGED_VALUE_OPERATORS_H"""
  sys.stdout = sys.__stdout__
  f.close()

if (__name__ == "__main__"):
  run()
