import sys
from string import upper
import generate_vector_algebra_traits

arithmetic_unary_ops = ("-")
arithmetic_binary_ops = ("+", "-", "*", "/", "%")
arithmetic_in_place_binary_ops = ("+=", "-=", "*=", "/=", "%=")
boolean_ops = ("==", "!=", ">", "<", ">=", "<=")
logical_unary_ops = ("!")
logical_binary_ops = ("&&", "||")
logical_in_place_binary_ops = ("&&=", "||=")

std_abinop_function_objects = {
  "+": "std::plus",
  "-": "std::minus",
  "*": "std::multiplies",
  "/": "std::divides",
  "%": "std::modulus",
}

aipbinop_function_objects = {
  "+=": "ip_plus",
  "-=": "ip_minus",
  "*=": "ip_multiplies",
  "/=": "ip_divides",
  "%=": "ip_modulus",
}

std_boolop_function_objects = {
  "==": "std::equal_to",
  "!=": "std::not_equal_to",
  ">":  "std::greater",
  "<":  "std::less",
  ">=": "std::greater_equal",
  "<=": "std::less_equal",
}

boolop_function_name = {
  "==": "vector_equal_to",
  "!=": "vector_not_equal_to",
  ">":  "vector_greater",
  "<":  "vector_less",
  ">=": "vector_greater_equal",
  "<=": "vector_less_equal",
}

std_lbinop_function_objects = {
  "&&": "logical_and",
  "||": "logical_or",
}

std_lipbinop_function_objects = {
  "&&=": "ip_logical_and",
  "||=": "ip_logical_or",
}

def get_fillin(LHS_is_vector, RHS_is_vector):
  v_size = None
  r_qf = "s"
  l_qf = "s"
  r_subscript = ""
  l_subscript = ""
  if (RHS_is_vector):
    v_size = "r"
    r_subscript = "[i]"
    r_qf = "v"
  if (LHS_is_vector):
    v_size = "l"
    l_subscript = "[i]"
    l_qf = "v"
  return l_qf, r_qf, v_size, l_subscript, r_subscript

def one_abinop(LHS_is_vector, RHS_is_vector):
  l_qf, r_qf, v_size, l_subscript, r_subscript = get_fillin(
    LHS_is_vector, RHS_is_vector)
  print """
  template <typename FunctionObjectType, typename TypeLHS, typename TypeRHS>
  typename algebra_traits<TypeLHS, TypeRHS>::promotion_type
  abinop(const FunctionObjectType& op,
         qualifier_%s_%s, const TypeLHS& lhs, const TypeRHS& rhs)
  {""" % (l_qf, r_qf)
  if (LHS_is_vector and RHS_is_vector):
    print "    assert(lhs.size() == rhs.size());"
  if (LHS_is_vector or RHS_is_vector):
    print \
"""    typedef typename algebra_traits<TypeLHS, TypeRHS>::promotion_type
    promotion_type;
    algebra_constructor<promotion_type> result(%shs.size());
    for(typename Type%sHS::size_type i=0;i<%shs.size();i++)
      result[i] = op(lhs%s, rhs%s);
    return result;
  }""" % (v_size, upper(v_size), v_size, l_subscript, r_subscript)
  else:
    print \
"""    return op(lhs, rhs);
  }"""

def generate_abinops():
  for LHS_is_vector, RHS_is_vector in ((1,1), (1,0), (0,1), (0,0)):
    one_abinop(LHS_is_vector, RHS_is_vector)
  for op in arithmetic_binary_ops:
    print """
  template <typename TypeLHS, typename TypeRHS>
  typename algebra_traits<TypeLHS, TypeRHS>::promotion_type
  operator%s(const TypeLHS& lhs, const TypeRHS& rhs) {
    typedef typename algebra_traits<TypeLHS, TypeRHS>::value_type value_type;
    typedef typename algebra_traits<TypeLHS, TypeRHS>::type_qualifier qf;
    return abinop(%s<value_type>(), qf(), lhs, rhs);
  }""" % (op, std_abinop_function_objects[op])

def one_aipbinop(LHS_is_vector, RHS_is_vector):
  l_qf, r_qf, v_size, l_subscript, r_subscript = get_fillin(
    LHS_is_vector, RHS_is_vector)
  print """
  template <typename FunctionObjectType, typename TypeLHS, typename TypeRHS>
  TypeLHS&
  aipbinop(const FunctionObjectType& op,
           qualifier_%s_%s, TypeLHS& lhs, const TypeRHS& rhs)
  {""" % (l_qf, r_qf)
  if (LHS_is_vector and RHS_is_vector):
    print "    assert(lhs.size() == rhs.size());"
  print \
"""    typedef typename algebra_traits<TypeLHS, TypeLHS>::value_type
    value_type;"""
  if (LHS_is_vector or RHS_is_vector):
    print \
"""    for(typename Type%sHS::size_type i=0;i<%shs.size();i++)
      op(lhs%s, value_type(rhs%s));
    return lhs;
  }""" % (upper(v_size), v_size, l_subscript, r_subscript)
  else:
    print \
"""    op(lhs, value_type(rhs));
    return lhs;
  }"""

def generate_aipbinop_function_objects():
  for op in arithmetic_in_place_binary_ops:
    print """
  template <typename T>
  struct %s : std::binary_function<T, T, T>
  {
    T& operator()(T& x, const T& y) const { x %s y; return x; }
  };""" % (aipbinop_function_objects[op], op)

def generate_aipbinops():
  generate_aipbinop_function_objects()
  for LHS_is_vector, RHS_is_vector in ((1,1), (1,0), (0,1), (0,0)):
    one_aipbinop(LHS_is_vector, RHS_is_vector)
  for op in arithmetic_in_place_binary_ops:
    print """
  template <typename TypeLHS, typename TypeRHS>
  TypeLHS&
  operator%s(TypeLHS& lhs, const TypeRHS& rhs) {
    typedef typename algebra_traits<TypeLHS, TypeLHS>::value_type value_type;
    typedef typename algebra_traits<TypeLHS, TypeRHS>::type_qualifier qf;
    return aipbinop(%s<value_type>(), qf(), lhs, rhs);
  }""" % (op, aipbinop_function_objects[op])

def one_boolop(LHS_is_vector, RHS_is_vector):
  l_qf, r_qf, v_size, l_subscript, r_subscript = get_fillin(
    LHS_is_vector, RHS_is_vector)
  print """
  template <typename FunctionObjectType, typename TypeLHS, typename TypeRHS>
  typename algebra_traits<TypeLHS, TypeRHS>::bool_type
  boolop(const FunctionObjectType& op,
         qualifier_%s_%s, const TypeLHS& lhs, const TypeRHS& rhs)
  {""" % (l_qf, r_qf)
  if (LHS_is_vector and RHS_is_vector):
    print "    assert(lhs.size() == rhs.size());"
  if (LHS_is_vector or RHS_is_vector):
    print \
"""    typedef typename algebra_traits<TypeLHS, TypeRHS>::bool_type
    bool_type;
    algebra_constructor<bool_type> result(%shs.size());
    for(typename Type%sHS::size_type i=0;i<%shs.size();i++)
      result[i] = op(lhs%s, rhs%s);
    return result;
  }""" % (v_size, upper(v_size), v_size, l_subscript, r_subscript)
  else:
    print \
"""    return op(lhs, rhs);
  }"""

def generate_boolops():
  for LHS_is_vector, RHS_is_vector in ((1,1), (1,0), (0,1), (0,0)):
    one_boolop(LHS_is_vector, RHS_is_vector)
  for op in boolean_ops:
    print """
  template <typename TypeLHS, typename TypeRHS>
  typename algebra_traits<TypeLHS, TypeRHS>::bool_type
  %s(const TypeLHS& lhs, const TypeRHS& rhs) {
    typedef typename algebra_traits<TypeLHS, TypeRHS>::value_type value_type;
    typedef typename algebra_traits<TypeLHS, TypeRHS>::type_qualifier qf;
    return boolop(%s<value_type>(), qf(), lhs, rhs);
  }""" % (boolop_function_name[op], std_boolop_function_objects[op])

def run():
  f = open("algebra_operators.h", "w")
  sys.stdout = f
  generate_vector_algebra_traits.write_copyright()
  print """
#ifndef CCTBX_VECTOR_ALGEBRA_OPERATORS_H
#define CCTBX_VECTOR_ALGEBRA_OPERATORS_H

#include <complex>
#include <cassert>

namespace cctbx { namespace vector {"""
  generate_abinops()
  generate_aipbinops()
  generate_boolops()
  print """
}} // namespace cctbx::vector

#endif // CCTBX_VECTOR_ALGEBRA_OPERATORS_H"""
  sys.stdout = sys.__stdout__
  f.close()

if (__name__ == "__main__"):
  run()
