from scitbx.source_generators.array_family import operator_functor_info
from scitbx.source_generators.array_family import generate_std_imports
from scitbx.source_generators import utils

this = "scitbx.source_generators.array_family.generate_algebras"

misc_functions_a = (
  "absolute", "pow2",
)

misc_functions_a_a = (
  "fmod_positive",
)

misc_functions_x_x_s = (
  (["approx_equal"], 1,
   ["ElementType const& tolerance", "tolerance"]),
)

class empty(object):
  def join(self, other):
    self.__dict__.update(other.__dict__)

def format_header(indent, str, max_line_length = 79):
  maxlen = max_line_length - len(indent)
  extra_indent = ""
  lei = len(extra_indent)
  result = ""
  rest = str.strip()
  while (lei + len(rest) > maxlen):
    if (lei == 0):
      i = rest.index("<")
    else:
      i = rest.index(",")
      try: i += rest[i+1:].index(",") + 1
      except Exception: pass
    result += indent + extra_indent + rest[:i+1] + '\n'
    extra_indent = "  "
    lei = 2
    rest = rest[i+1:].strip()
  result += indent + extra_indent + rest
  return result

def format_list(indent, list):
  r = ""
  for line in list[:-1]:
    r += indent + line + "\n"
  return r + indent + list[-1]

def base_array_type_name(array_type_name):
  return array_type_name.replace("_plain", "")

def get_template_args(array_type_name):
  result = [["typename", "ElementType"]]
  if (base_array_type_name(array_type_name) in ("tiny", "small")):
    result.append(["std::size_t", "N"])
  elif (base_array_type_name(array_type_name) in ("ref","const_ref", "versa")):
    result.append(["typename", "AccessorType"])
  return result

def get_numbered_template_args(array_type_name, n_params, equal_element_type):
  single = get_template_args(array_type_name)
  if (n_params == 1): return [single]
  result = []
  for i in xrange(1, n_params+1):
    if (equal_element_type):
      single_numbered = [single[0]]
    else:
      single_numbered = [[single[0][0], single[0][1]+str(i)]]
    if (base_array_type_name(array_type_name) == "tiny"):
      single_numbered.append(single[1])
    else:
      for s in single[1:]:
        single_numbered.append([s[0], s[1]+str(i)])
    result.append(single_numbered)
  return result

def get_template_header_args(numbered_template_args):
  result = []
  for p in numbered_template_args:
    for d in p:
      if (not d in result): result.append(d)
  return result

def get_template_parameter_args(numbered_template_args):
  from string import join
  result = []
  for p in numbered_template_args:
    result.append(join([d[1] for d in p], ", "))
  return result

def get_template_header(numbered_template_args):
  from string import join
  ha = get_template_header_args(numbered_template_args)
  result = "template<" + join([join(x) for x in ha], ", ") + ">"
  return result

def get_template_parameters(array_type_name, template_parameter_args):
  result = []
  for p in template_parameter_args:
    result.append(array_type_name + "<" + p + ">")
  return result

def get_template_header_and_parameters(
  array_type_name, n_params, equal_element_type = 0
):
  result = empty()
  result.nta = get_numbered_template_args(
    array_type_name, n_params, equal_element_type)
  result.tpa = get_template_parameter_args(result.nta)
  result.header = get_template_header(result.nta)
  result.params = get_template_parameters(array_type_name, result.tpa)
  return result

def derive_return_array_type_simple(param):
  if (not param.startswith("small")): return param
  return param.replace("N1", "SCITBX_ARRAY_FAMILY_SMALL_ALGEBRA_MIN_N1_N2")

def wrap_element_type(array_type_name, element_type, addl):
  from string import join
  r = array_type_name + "<" + join([element_type] + addl, ", ") + ">"
  if (r.endswith(">>")): return r[:-1] + " >"
  return r

def adjust_array_type_names(array_type_name):
  r = empty()
  r.return_array_type_name = array_type_name
  r.param_array_type_name = array_type_name
  if (array_type_name == "ref"):
    r.return_array_type_name = "versa"
    r.param_array_type_name = "const_ref"
  return r

def special_decl_params(array_type_name, special_def):
  r = adjust_array_type_names(array_type_name)
  r.return_elment_type = special_def[0]
  r.function_name = special_def[1]
  r.arg_element_types = special_def[2:]
  r.nta = get_numbered_template_args(r.param_array_type_name, 1, 0)
  addl = []
  if (len(r.nta[0]) == 2): addl = [r.nta[0][1][1]]
  r.header = get_template_header(r.nta)
  r.return_array_type = wrap_element_type(
    r.return_array_type_name, r.return_elment_type, addl)
  r.arg_array_types = []
  for aet in r.arg_element_types:
    r.arg_array_types.append(wrap_element_type(
      r.param_array_type_name, aet, addl))
  return r

def operator_decl_params(array_type_name, op_type, op_class, type_flags,
  equal_element_type = 0
):
  r = adjust_array_type_names(array_type_name)
  if (type_flags != (1,1)):
    r.join(get_template_header_and_parameters(r.param_array_type_name, 1,
      equal_element_type))
    r.params.insert(type_flags[0], "ElementType")
    r.return_element_type = ["ElementType"]
    if (op_class in ("boolean", "bool_result", "logical")):
      r.return_element_type = ["bool"]
  else:
    r.join(get_template_header_and_parameters(
      r.param_array_type_name, 2, equal_element_type))
    r.return_element_type = [r.nta[0][0][1]]
    if (op_class in ("boolean", "bool_result", "logical")):
      r.return_element_type = ["bool"]
    elif (op_class != "n/a"):
      assert op_class == "arithmetic"
      r.return_element_type = [
        "typename binary_operator_traits<",
        "  ElementType1, ElementType2>::" + op_class]
  if (len(r.return_element_type) == 1):
    r.return_array_type = [
      r.return_array_type_name + "<" + r.return_element_type[0]]
  else:
    r.return_array_type = [
      r.return_array_type_name + "<",
      "  " + r.return_element_type[0],
      "    " + r.return_element_type[1]]
  if (len(r.nta[0]) == 2):
    if (r.nta[0][1][1] == "N1"):
      r.return_array_type[-1] +=", SCITBX_ARRAY_FAMILY_SMALL_ALGEBRA_MIN_N1_N2"
    else:
      r.return_array_type[-1] += ", " + r.nta[0][1][1]
  r.return_array_type[-1] += ">"
  r.typedef_return_array_type = (["typedef " +  r.return_array_type[0]]
    + r.return_array_type[1:])
  r.element_types = ["ElementType", "ElementType"]
  if (type_flags == (1,1)):
    r.element_types = ["ElementType1", "ElementType2"]
  return r

def get_result_constructor_args(array_type_name, type_flags = None):
  arg_name = "a"
  if (type_flags is not None):
    arg_name = "a%d" % ((type_flags[0] + 1) % 2 + 1,)
  if (base_array_type_name(array_type_name) == "versa"):
    return "%s.accessor()" % (arg_name,)
  return "%s.size()" % (arg_name,)

def binary_operator_algo_params(array_type_name, type_flags):
  r = empty()
  r.loop_n = "N"
  r.size_assert = ""
  r.result_constructor_args = ""
  if (base_array_type_name(array_type_name) != "tiny"):
    r.loop_n = "a%d.size()" % ((type_flags[0] + 1) % 2 + 1,)
    if (type_flags == (1,1)):
      r.size_assert = """if (a1.size() != a2.size()) throw_range_error();
    """
    r.result_constructor_args = get_result_constructor_args(
      array_type_name, type_flags)
  r.begin = ["", ""]
  for i in xrange(2):
    if (type_flags[i]): r.begin[i] = ".begin()"
  r.type_flags_code = "sa"[type_flags[0]] + "_" + "sa"[type_flags[1]]
  return r

def generate_unary_ops(f, array_type_name):
  result_constructor_args = get_result_constructor_args(array_type_name)
  for op_class, op_symbol in (("arithmetic", "-"),
                              ("logical", "!")):
    d = operator_decl_params(array_type_name, "unary", op_class, (1,0))
    print >> f, """%s
  inline
%s
  operator%s(%s const& a) {
%s
    return_array_type;
    typedef typename return_array_type::value_type return_element_type;""" % (
      format_header("  ", d.header),
      format_list("  ", d.return_array_type),
      op_symbol, d.params[0],
      format_list("    ", d.typedef_return_array_type))
    if (base_array_type_name(array_type_name) == "tiny"):
      print >> f, """    return_array_type result;
    array_operation_a(fn::functor_%s<
        return_element_type,
        ElementType>(),
      a.begin(), result.begin(), a.size(), true_type());
    return result;
  }
""" % (operator_functor_info.unary_functors[op_symbol],)
    else:
      print >> f, """    return return_array_type(%s,
      make_init_functor(make_array_functor_a(
        fn::functor_%s<
          return_element_type,
          ElementType>(), a.begin())));
  }
""" % (result_constructor_args,
       operator_functor_info.unary_functors[op_symbol])

def generate_unary_apply(f, array_type_name):
  result_constructor_args = get_result_constructor_args(array_type_name)
  d = operator_decl_params(array_type_name, "unary", "n/a", (1,0))
  nta = get_numbered_template_args(array_type_name, 1, 0)
  addl = []
  if (len(nta[0]) == 2): addl = [nta[0][1][1]]
  return_array_type = wrap_element_type(
    d.return_array_type_name, "ReturnElementType", addl)
  nta.insert(0, [["typename", "UnaryOperation"]])
  header2 = get_template_header(nta)
  nta.append([["typename", "ReturnElementType"]])
  header1 = get_template_header(nta)
  print >> f, """%s
  inline
  %s
  apply(UnaryOperation const& op,
        %s const& a,
        type_holder<ReturnElementType> /*result_type_holder*/) {
    typedef %s return_array_type;""" % (
       format_header("  ", header1),
       return_array_type,
       d.params[0],
       return_array_type)
  if (base_array_type_name(array_type_name) == "tiny"):
    print >> f, """    return_array_type result;
    array_operation_a(op,
      a.begin(), result.begin(), a.size(), true_type());
    return result;
  }"""
  else:
    print >> f, """    return return_array_type(%s,
      make_init_functor(
        array_functor_a<UnaryOperation, ElementType, ReturnElementType>(
          op, a.begin())));
  }""" % (result_constructor_args,)
  print >> f, """
%s
  inline
  %s
  apply(UnaryOperation const& op,
        %s const& a) {
    return apply(op, a, type_holder<typename UnaryOperation::result_type>());
  }
""" % (format_header("  ", header2),
       return_array_type.replace(
         "ReturnElementType", "typename UnaryOperation::result_type"),
       d.params[0])

def elementwise_binary_op(f,
      array_type_name, op_class, op_symbol, type_flags, function_name):
  d = operator_decl_params(array_type_name, "binary", op_class, type_flags)
  a = binary_operator_algo_params(array_type_name, type_flags)
  print >> f, """%s
  inline
%s
  %s(
    %s const& a1,
    %s const& a2) {
%s
    return_array_type;
    typedef typename return_array_type::value_type return_element_type;""" % (
      format_header("  ", d.header),
      format_list("  ", d.return_array_type),
      function_name, d.params[0], d.params[1],
      format_list("    ", d.typedef_return_array_type))
  if (base_array_type_name(array_type_name) == "tiny"):
    print >> f, """    return_array_type result;
    array_operation_%s(fn::functor_%s<
        return_element_type,
        %s,
        %s>(),
      a1%s, a2%s, result.begin(), %s, true_type());
    return result;
  }
""" % (a.type_flags_code,
       operator_functor_info.binary_functors[op_symbol],
       d.element_types[0], d.element_types[1],
       a.begin[0], a.begin[1], a.loop_n)
  else:
    print >> f, """    %sreturn return_array_type(%s,
      make_init_functor(make_array_functor_%s(
        fn::functor_%s<
          return_element_type,
          %s,
          %s>(), a1%s, a2%s)));
  }
""" % (a.size_assert,
       a.result_constructor_args,
       a.type_flags_code,
       operator_functor_info.binary_functors[op_symbol],
       d.element_types[0], d.element_types[1],
       a.begin[0], a.begin[1])

def elementwise_inplace_binary_op(f,
      array_type_name, op_class, op_symbol, type_flags):
  d = operator_decl_params(array_type_name, "binary", "n/a", type_flags)
  d.params[0] = d.params[0].replace("const_ref", "ref")
  a = binary_operator_algo_params(array_type_name, type_flags)
  print >> f, """%s
  inline
  %s&
  operator%s(
    %s& a1,
    %s const& a2) {
    %sarray_operation_in_place_%s(fn::functor_%s<
        %s,
        %s>(),
      a1.begin(), a2%s, %s);
    return a1;
  }
""" % (format_header("  ", d.header),
       d.params[0],
       op_symbol, d.params[0], d.params[1],
       a.size_assert,
       a.type_flags_code,
       operator_functor_info.in_place_binary_functors[op_symbol],
       d.return_element_type[0],
       d.element_types[1],
       a.begin[1], a.loop_n);

def generate_elementwise_binary_op(f, array_type_name, op_class, op_symbol):
  function_name = "operator" + op_symbol
  for type_flags in ((1,1), (1,0), (0,1)):
    elementwise_binary_op(f,
      array_type_name, op_class, op_symbol, type_flags, function_name)

def generate_elementwise_inplace_binary_op(f,
      array_type_name, op_class, op_symbol):
  for type_flags in ((1,1), (1,0)):
    elementwise_inplace_binary_op(f,
      array_type_name, op_class, op_symbol, type_flags)

def generate_1arg_element_wise(f, array_type_name, function_names):
  result_constructor_args = get_result_constructor_args(array_type_name)
  d = operator_decl_params(array_type_name, "unary", "arithmetic", (1,0))
  for function_name in function_names:
    print >> f, """%s
  inline
%s
  %s(%s const& a) {
%s
    return_array_type;
    typedef typename return_array_type::value_type return_element_type;""" % (
      format_header("  ", d.header),
      format_list("  ", d.return_array_type),
      function_name, d.params[0],
      format_list("    ", d.typedef_return_array_type))
    if (base_array_type_name(array_type_name) == "tiny"):
      print >> f, """    return_array_type result;
    array_operation_a(fn::functor_%s<return_element_type, ElementType>(),
      a.begin(), result.begin(), a.size(), true_type());
    return result;
  }
""" % (function_name,)
    else:
      print >> f, """    return return_array_type(%s,
      make_init_functor(make_array_functor_a(
        fn::functor_%s<return_element_type, ElementType>(), a.begin())));
  }
""" % (result_constructor_args, function_name)

def generate_2arg_element_wise(f,
  array_type_name, function_names,
  equal_element_type = 0
):
  for function_name in function_names:
    for type_flags in ((1,1), (1,0), (0,1)):
      d = operator_decl_params(
        array_type_name, "binary", "n/a", type_flags, equal_element_type)
      a = binary_operator_algo_params(array_type_name, type_flags)
      print >> f, """%s
  inline
%s
  %s(
    %s const& a1,
    %s const& a2) {
%s
    return_array_type;
    typedef typename return_array_type::value_type return_element_type;""" % (
       format_header("  ", d.header),
       format_list("  ", d.return_array_type),
       function_name, d.params[0], d.params[1],
       format_list("    ", d.typedef_return_array_type))
      if (base_array_type_name(array_type_name) == "tiny"):
        print >> f, """
    return_array_type result;
    array_operation_%s(fn::functor_%s<
        return_element_type, %s, %s>(),
      a1%s, a2%s, result.begin(), %s, true_type());
    return result;
  }
""" % (a.type_flags_code,
       function_name, d.element_types[0], d.element_types[1],
       a.begin[0], a.begin[1], a.loop_n)
      else:
        print >> f, """    %sreturn return_array_type(%s,
      make_init_functor(make_array_functor_%s(
        fn::functor_%s<return_element_type,
          %s, %s>(),
        a1%s, a2%s)));
  }
""" % (a.size_assert,
       a.result_constructor_args,
       a.type_flags_code,
       function_name, d.element_types[0], d.element_types[1],
       a.begin[0], a.begin[1])

def generate_x_x_s_element_wise(f,
  array_type_name, function_names,
  equal_element_type,
  addl_args
):
  for function_name in function_names:
    for type_flags in ((1,1), (1,0), (0,1)):
      d = operator_decl_params(
        array_type_name, "binary", "bool_result",
        type_flags, equal_element_type)
      a = binary_operator_algo_params(array_type_name, type_flags)
      print >> f, """%s
  inline
%s
  %s(
    %s const& a1,
    %s const& a2,
    %s) {
%s
    return_array_type;
    typedef typename return_array_type::value_type return_element_type;""" % (
      (format_header("  ", d.header),
       format_list("  ", d.return_array_type),
       function_name, d.params[0], d.params[1], addl_args[0],
       format_list("    ", d.typedef_return_array_type)))
      if (base_array_type_name(array_type_name) == "tiny"):
        print >> f, """
    return_array_type result;
    array_operation_%s_s(fn::functor_%s<
        return_element_type, %s, %s, %s>(),
      a1%s, a2%s, %s,
      result.begin(), %s, true_type());
    return result;
  }
""" % (a.type_flags_code,
       function_name, "ElementType", "ElementType", "ElementType",
       a.begin[0], a.begin[1], addl_args[1], a.loop_n)
      else:
        print >> f, """    %sreturn return_array_type(%s,
      make_init_functor(make_array_functor_%s_s(
        fn::functor_%s<return_element_type,
          %s, %s, %s>(),
        a1%s, a2%s, %s)));
  }
""" % (a.size_assert,
       a.result_constructor_args,
       a.type_flags_code,
       function_name, "ElementType", "ElementType", "ElementType",
       a.begin[0], a.begin[1], addl_args[1])

def generate_element_wise_special(f,
  array_type_name, special_def
):
  p = special_decl_params(array_type_name, special_def)
  if (len(p.arg_array_types) == 1):
    print >> f, """%s
  inline
  %s
  %s(%s const& a) {
    typedef %s return_array_type;""" % (
      format_header("  ", p.header),
      p.return_array_type,
      p.function_name, p.arg_array_types[0],
      p.return_array_type)
    if (base_array_type_name(array_type_name) == "tiny"):
      print >> f, """    return_array_type result;
    array_operation_a(fn::functor_%s<
      %s,
      %s >(),
      a.begin(), result.begin(), a.size(), true_type());
    return result;
  }
""" % (p.function_name,
       special_def[0], special_def[2])
    else:
      result_constructor_args = get_result_constructor_args(
        array_type_name)
      print >> f, """    return return_array_type(%s,
      make_init_functor(make_array_functor_a(
        fn::functor_%s<
          %s,
          %s >(), a.begin())));
  }
""" % (result_constructor_args,
       p.function_name,
       special_def[0], special_def[2])
  else:
    for type_flags in ((1,1), (1,0), (0,1)):
      a = binary_operator_algo_params(array_type_name, type_flags)
      params = []
      for i in xrange(2):
        if (type_flags[i]):
          params.append(p.arg_array_types[i])
        else:
          params.append(p.arg_element_types[i])
      print >> f, """%s
  inline
  %s
  %s(
    %s const& a1,
    %s const& a2) {
    typedef %s return_array_type;""" % (
      format_header("  ", p.header), p.return_array_type,
      p.function_name, params[0], params[1],
      p.return_array_type)
      if (base_array_type_name(array_type_name) == "tiny"):
        print >> f, """    return_array_type result;
    array_operation_%s(fn::functor_%s<
        %s,
        %s,
        %s >(),
      a1%s, a2%s, result.begin(), %s, true_type());
    return result;
  }
""" % (a.type_flags_code,
       p.function_name,
       special_def[0], special_def[2], special_def[3],
       a.begin[0], a.begin[1], a.loop_n)
      else:
        print >> f, """    %sreturn return_array_type(%s,
      make_init_functor(make_array_functor_%s(
        fn::functor_%s<
          %s,
          %s,
          %s >(), a1%s, a2%s)));
  }
""" % (a.size_assert,
       a.result_constructor_args,
       a.type_flags_code,
       p.function_name,
       special_def[0], special_def[2], special_def[3],
       a.begin[0], a.begin[1])

def one_type(target_dir, array_type_name):
  f = utils.join_open(target_dir, "%s_algebra.h" % array_type_name, "w")
  utils.write_this_is_auto_generated(f, this)
  include_array_type_name = array_type_name
  if (array_type_name == "ref"):
    include_array_type_name = "versa"
  generic_include = "functors"
  if (base_array_type_name(array_type_name) == "tiny"):
    generic_include = "operators"
  print >> f, """\
#ifndef SCITBX_ARRAY_FAMILY_%s_ALGEBRA_H
#define SCITBX_ARRAY_FAMILY_%s_ALGEBRA_H

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#include <scitbx/array_family/%s.h>
""" % ((array_type_name.upper(),) * 2 + (include_array_type_name,))
  if (array_type_name == "small"):
    print >> f, """#if (defined(BOOST_MSVC) && BOOST_MSVC <= 1300) // VC++ 7.0
#define SCITBX_ARRAY_FAMILY_SMALL_ALGEBRA_MIN_N1_N2 N1
#else
#define SCITBX_ARRAY_FAMILY_SMALL_ALGEBRA_MIN_N1_N2 (N1<N2?N1:N2)
#endif

"""
  print >> f, """#include <scitbx/array_family/operator_traits_builtin.h>
#include <scitbx/array_family/detail/operator_functors.h>
#include <scitbx/array_family/detail/generic_array_%s.h>
#include <scitbx/array_family/detail/std_imports.h>
#include <scitbx/array_family/misc_functions.h>

namespace scitbx { namespace af {
""" % (generic_include,)

  generate_unary_ops(f, array_type_name)
  for op_symbol in operator_functor_info.arithmetic_binary_ops:
    generate_elementwise_binary_op(f,
      array_type_name, "arithmetic", op_symbol)
    generate_elementwise_inplace_binary_op(f,
      array_type_name, "arithmetic", op_symbol + "=")
  for op_symbol in operator_functor_info.logical_binary_ops:
    generate_elementwise_binary_op(f,
      array_type_name, "logical", op_symbol)
  for op_symbol in operator_functor_info.boolean_binary_ops:
    generate_elementwise_binary_op(f,
      array_type_name, "boolean", op_symbol)
  generate_1arg_element_wise(f,
    array_type_name,
    misc_functions_a
    + generate_std_imports.cmath_1arg
    + generate_std_imports.cstdlib_1arg
    + generate_std_imports.complex_1arg)
  generate_2arg_element_wise(f,
    array_type_name,
    misc_functions_a_a
    + generate_std_imports.cmath_2arg
    + generate_std_imports.algorithm_2arg)
  for special_def in generate_std_imports.complex_special:
    generate_element_wise_special(f, array_type_name, special_def)
  for args in misc_functions_x_x_s:
    apply(generate_x_x_s_element_wise, (f, array_type_name) + args)

  print >> f, "}} // namespace scitbx::af"
  print >> f
  print >> f, "#endif // DOXYGEN_SHOULD_SKIP_THIS"
  print >> f
  print >> f, "#endif // SCITBX_ARRAY_FAMILY_%s_ALGEBRA_H" % (
    array_type_name.upper(),)

  f.close()

def run(target_dir):
  for array_type_name in ("ref", "tiny", "small", "shared", "versa"):
    one_type(target_dir, array_type_name)

if (__name__ == "__main__"):
  run(".")
