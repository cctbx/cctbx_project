def write_copyright():
  print """// This is an automatically generated file. Do not edit.
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     Jan 2002: Created, based on generate_tiny_operators.py (rwgk)
 */"""

arithmetic_unary_ops = ("-")
arithmetic_binary_ops = ("+", "-", "*", "/", "%")
arithmetic_in_place_binary_ops = ("+=", "-=", "*=", "/=", "%=")
logical_unary_ops = ("!")
logical_binary_ops = ("&&", "||")
boolean_ops = ("==", "!=", ">", "<", ">=", "<=")

class empty: pass

def form_param(array_type_name, access, Xhs, distinct_N):
  if (access == ""): return "ElementType" + Xhs
  NXhs = ""
  if (distinct_N):
    NXhs = ", N" + Xhs
  else:
    NXhs = ", N"
  return "%s<ElementType%s%s>" % (array_type_name, Xhs, NXhs)

def op_vars(array_type_name, access_lhs, access_rhs):
  v = empty()
  v.result_constructor_args = ""
  v.loop_n = "N"
  v.size_assert = ""
  if (array_type_name != "tiny"):
    if (access_lhs != ""):
      v.result_constructor_args = "(lhs.size())"
      v.loop_n = "lhs.size()"
    else:
      v.result_constructor_args = "(rhs.size())"
      v.loop_n = "rhs.size()"
    if (access_lhs != "" and access_rhs != ""):
      v.size_assert = """if (lhs.size() != rhs.size()) throw_range_error();
    """
  if (array_type_name != "tiny" and access_lhs != "" and access_rhs != ""):
    v.template_head = \
"""template <typename ElementTypeLhs, std::size_t NLhs,
            typename ElementTypeRhs, std::size_t NRhs>"""
    v.Nresult = "((NLhs<NRhs)?NLhs:NRhs)"
    distinct_N = 1
  else:
    v.template_head = \
"""template <typename ElementTypeLhs,
            typename ElementTypeRhs, std::size_t N>"""
    v.Nresult = "N"
    distinct_N = 0
  v.param_lhs = form_param(array_type_name, access_lhs, "Lhs", distinct_N)
  v.param_rhs = form_param(array_type_name, access_rhs, "Rhs", distinct_N)
  return v

def elementwise_binary_op(
      op_class, op_symbol, function_name,
      array_type_name, access_lhs, access_rhs):
  v = op_vars(array_type_name, access_lhs, access_rhs)
  print """  %s
  inline
  %s<
    typename binary_operator_traits<
      ElementTypeLhs, ElementTypeRhs>::%s, %s>
  %s(
    const %s& lhs,
    const %s& rhs) {
    %s<
      typename binary_operator_traits<
        ElementTypeLhs, ElementTypeRhs>::%s, %s>
    result%s;
    %sfor(std::size_t i=0;i<%s;i++) result[i] = lhs%s %s rhs%s;
    return result;
  }
""" % (v.template_head, array_type_name, op_class, v.Nresult,
       function_name, v.param_lhs, v.param_rhs,
       array_type_name, op_class, v.Nresult,
       v.result_constructor_args, v.size_assert, v.loop_n,
       access_lhs, op_symbol, access_rhs)

def elementwise_inplace_binary_op(
      op_symbol, function_name,
      array_type_name, access_rhs):
  v = op_vars(array_type_name, "[i]", access_rhs)
  print """  %s
  inline
  %s<ElementTypeLhs, %s>
  %s(
    const %s<ElementTypeLhs, %s>& lhs,
    const %s& rhs) {
    %sfor(std::size_t i=0;i<%s;i++) lhs[i] %s= rhs%s;
    return lhs;
  }
""" % (v.template_head, array_type_name, v.Nresult,
       function_name, array_type_name, v.Nresult,
       v.param_rhs, v.size_assert, v.loop_n, op_symbol, access_rhs)

def generate_elementwise_binary_op(
      array_type_name,
      inplace, op_class, op_symbol, function_name = None):
  if (function_name == None):
    function_name = "operator" + op_symbol
  elementwise_binary_op(op_class, op_symbol, function_name,
    array_type_name, "[i]", "[i]")
  elementwise_binary_op(op_class, op_symbol, function_name,
    array_type_name, "[i]", "")
  elementwise_binary_op(op_class, op_symbol, function_name,
    array_type_name, "", "[i]")
  if (inplace):
    elementwise_inplace_binary_op(op_symbol, "operator" + op_symbol + "=",
      array_type_name, "[i]")
    elementwise_inplace_binary_op(op_symbol, "operator" + op_symbol + "=",
      array_type_name, "")

def reducing_boolean_op(
      op_symbol,
      array_type_name, access_lhs, access_rhs,
      truth_test_type):
  v = op_vars(array_type_name, access_lhs, access_rhs)
  if (op_symbol == "=="):
    tests = (
"""      if (lhs%s != rhs%s) return %s() != %s();"""
    % (access_lhs, access_rhs, truth_test_type, truth_test_type))
    final_op = "=="
  elif (op_symbol == "!="):
    tests = (
"""      if (lhs%s != rhs%s) return %s() == %s();"""
    % (access_lhs, access_rhs, truth_test_type, truth_test_type))
    final_op = "!="
  elif (op_symbol in ("<", ">")):
    tests = (
"""      if (lhs%s %s rhs%s) return %s() == %s();
      if (rhs%s %s lhs%s) return %s() != %s();"""
    % (access_lhs, op_symbol, access_rhs, truth_test_type, truth_test_type,
       access_rhs, op_symbol, access_lhs, truth_test_type, truth_test_type))
    final_op = "!="
  elif (op_symbol in ("<=", ">=")):
    tests = (
"""      if (!(lhs%s %s rhs%s)) return %s() != %s();"""
    % (access_lhs, op_symbol, access_rhs, truth_test_type, truth_test_type))
    final_op = "=="
  print """  %s
  inline
  typename binary_operator_traits<
    ElementTypeLhs, ElementTypeRhs>::boolean
  operator%s(
    const %s& lhs,
    const %s& rhs) {
    %sfor(std::size_t i=0;i<%s;i++) {
%s
    }
    return %s() %s %s();
  }
""" % (v.template_head, op_symbol, v.param_lhs, v.param_rhs,
       v.size_assert, v.loop_n, tests,
       truth_test_type, final_op, truth_test_type)

def generate_reducing_boolean_op(array_type_name, op_symbol):
  reducing_boolean_op(op_symbol,
    array_type_name, "[i]", "[i]", "ElementTypeLhs")
  reducing_boolean_op(op_symbol,
    array_type_name, "[i]", "", "ElementTypeLhs")
  reducing_boolean_op(op_symbol,
    array_type_name, "", "[i]", "ElementTypeRhs")

def generate_unary_ops(array_type_name):
  result_constructor_args = ""
  if (array_type_name != "tiny"): result_constructor_args = "(a.size())"
  for op_class, op_symbol in (("arithmetic", "-"),
                              ("logical", "!")):
    print """  template <typename ElementType, std::size_t N>
  inline
  %s<
    typename unary_operator_traits<
      ElementType>::%s, N>
  operator%s(const %s<ElementType, N>& a) {
    %s<
      typename unary_operator_traits<
        ElementType>::%s, N> result%s;
    for(std::size_t i=0;i<a.size();i++) result[i] = %sa[i];
    return result;
  }
""" % (array_type_name, op_class, op_symbol, array_type_name,
       array_type_name, op_class, result_constructor_args,
       op_symbol)

def one_type(array_type_name):
  import sys
  f = open("%s_operators.h" % (array_type_name,), "w")
  sys.stdout = f
  write_copyright()
  print """
#ifndef CCTBX_ARRAY_FAMILY_%s_OPERATORS_H
#define CCTBX_ARRAY_FAMILY_%s_OPERATORS_H

#include <cctbx/array_family/operator_traits.h>

namespace cctbx { namespace af {
""" % ((array_type_name.upper(),) * 2)

  generate_unary_ops(array_type_name)
  for op_symbol in arithmetic_binary_ops:
    generate_elementwise_binary_op(
      array_type_name, 1, "arithmetic", op_symbol)
  for op_symbol in logical_binary_ops:
    generate_elementwise_binary_op(
      array_type_name, 0, "logical", op_symbol)
  for op_symbol, function_name in (
    ("==", "equal_to"),
    ("!=", "not_equal_to"),
    (">", "greater"),
    ("<", "less"),
    (">", "greater_equal"),
    ("<", "less_equal")):
    generate_elementwise_binary_op(
      array_type_name, 0, "boolean", op_symbol, function_name)
  for op_symbol in boolean_ops:
    generate_reducing_boolean_op(array_type_name, op_symbol)

  print """}} // namespace cctbx::af

#endif // CCTBX_ARRAY_FAMILY_%s_OPERATORS_H""" % (array_type_name.upper(),)
  sys.stdout = sys.__stdout__
  f.close()

def run():
  one_type("tiny")
  one_type("small")

if (__name__ == "__main__"):
  run()
