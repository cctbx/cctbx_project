import sys
import generate_vector_promotions

unary_operators = ("-")
arithmetic_operators = ("+", "-", "*", "/", "%")
boolean_operators = ("==", "!=", "<", "<=", ">", ">=")

def one_binop(op_symbol,
              TypeLHS, LHS_is_vector,
              TypeRHS, RHS_is_vector,
              TypeResult, ValueTypeResult):
  r_subscript = ""
  l_subscript = ""
  if (RHS_is_vector):
    v_size = "r"
    r_subscript = "[i]"
  if (LHS_is_vector):
    v_size = "l"
    l_subscript = "[i]"
  print """  inline %s\\
  operator%s(const %s& lhs, const %s& rhs) {\\
    %s result(%shs.size());\\
    for(std::size_t i=0;i<%shs.size();i++)\\
      result[i] = %s(lhs%s) %s %s(rhs%s);\\
    return result;\\
  }\\""" % (TypeResult, op_symbol, TypeLHS, TypeRHS,
            TypeResult, v_size,
            v_size,
            ValueTypeResult, l_subscript, op_symbol,
            ValueTypeResult, r_subscript)

def generate_arithmetic_binop(op_symbol,
                              ValueTypeLHS,
                              ValueTypeRHS,
                              ValueTypeResult):
  for LHS_is_vector, RHS_is_vector in ((1,1), (1,0), (0,1)):
    TypeLHS = ValueTypeLHS
    TypeRHS = ValueTypeRHS
    if (LHS_is_vector): TypeLHS = "VectorTemplate<" + TypeLHS + " >"
    if (RHS_is_vector): TypeRHS = "VectorTemplate<" + TypeRHS + " >"
    one_binop(op_symbol,
              TypeLHS, LHS_is_vector,
              TypeRHS, RHS_is_vector,
              "VectorTemplate<%s >" % (ValueTypeResult,), ValueTypeResult)

def run():
  f = open("algebra.h", "w")
  sys.stdout = f
  generate_vector_promotions.write_copyright()
  print """
#ifndef CCTBX_VECTOR_ALGEBRA_H
#define CCTBX_VECTOR_ALGEBRA_H

#include <complex>

#define CCTBX_VECTOR_ALGEBRA_DECLARE(VectorTemplate) \\
namespace cctbx { namespace vector {\\"""
  op_types, result_type = generate_vector_promotions.build_pairs()
  for op_symbol in arithmetic_operators:
    for i in xrange(len(op_types)):
      if (    op_symbol == "%"
          and (not result_type[i]
                   in generate_vector_promotions.integer_types)): continue
      if (result_type[i] == 0):
        generate_arithmetic_binop(
          op_symbol, op_types[i].lhs, op_types[i].rhs, op_types[i].lhs)
      else:
        generate_arithmetic_binop(
          op_symbol, op_types[i].lhs, op_types[i].rhs, result_type[i])
  print """}} // namespace cctbx::vector

#endif // CCTBX_VECTOR_ALGEBRA_H"""
  sys.stdout = sys.__stdout__
  f.close()

if (__name__ == "__main__"):
  run()
