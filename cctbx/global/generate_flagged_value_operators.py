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

def op_vars(type_flags, op_class):
  v = empty()
  v.param = ["ValueTypeLhs", "ValueTypeRhs"]
  v.dotv = ["", ""]
  if (type_flags[0]):
    v.template_types = "typename ValueTypeLhs"
    v.result_type = "  flagged_value<ValueTypeLhs>"
    v.param[0] = "flagged_value<ValueTypeLhs>"
    v.param[1] = "ValueTypeLhs"
    v.dotv[0] = ".v"
    v.have_both_test = "lhs.f"
  if (type_flags[1]):
    v.template_types = "typename ValueTypeRhs"
    v.result_type = "  flagged_value<ValueTypeRhs>"
    v.param[1] = "flagged_value<ValueTypeRhs>"
    v.param[0] = "ValueTypeRhs"
    v.dotv[1] = ".v"
    v.have_both_test = "rhs.f"
  if (type_flags[0] and type_flags[1]):
    v.template_types = "typename ValueTypeLhs, typename ValueTypeRhs"
    v.param[0] = "flagged_value<ValueTypeLhs>"
    v.param[1] = "flagged_value<ValueTypeRhs>"
    v.have_both_test = "lhs.f && rhs.f"
    v.result_type = """  flagged_value<
    typename binary_operator_traits<
      ValueTypeLhs, ValueTypeRhs>::%s>""" % (op_class,)
  return v

def elementwise_binary_op(op_class, op_symbol, type_flags):
  v = op_vars(type_flags, op_class)
  print """  template <%s>
  inline
%s
  operator%s(
    const %s& lhs,
    const %s& rhs) {
%s result;
    if (%s) {
      result.v = lhs%s %s rhs%s;
      result.f = true;
    }
    return result;
  }
""" % (v.template_types, v.result_type,
       op_symbol, v.param[0], v.param[1], v.result_type,
       v.have_both_test, v.dotv[0], op_symbol, v.dotv[1])

def elementwise_inplace_binary_op(op_symbol, type_flags):
  v = op_vars(type_flags, "n/a")
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

def generate_unary_ops():
  for op_class, op_symbol in (("arithmetic", "-"),
                              ("logical", "!")):
    print """  template <typename ValueType>
  inline
  flagged_value<
    typename unary_operator_traits<
      ValueType>::%s>
  operator%s(
    const flagged_value<ValueType>& fv) {
    flagged_value<
      typename unary_operator_traits<
        ValueType>::%s> result;
    if (fv.f) {
      result.v = %sfv.v;
      result.f = true;
    }
    return result;
  }
""" % (op_class, op_symbol, op_class, op_symbol)

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

  generate_unary_ops()
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
