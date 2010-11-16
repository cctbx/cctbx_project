from scitbx.source_generators.array_family import operator_functor_info
from scitbx.source_generators import utils

this = "scitbx.source_generators.array_family.generate_operator_functors"

def generate_unary(f, name, op):
  print >> f, """
  template <typename ResultType,
            typename ArgumentType>
  struct functor_%s {
    typedef ResultType result_type;
    ResultType operator()(ArgumentType const& x) const {
      return ResultType(%s);
    }
  };""" % (name, op)

def generate_binary(f, name, op):
  print >> f, """
  template <typename ResultType,
            typename ArgumentType1,
            typename ArgumentType2>
  struct functor_%s {
    typedef ResultType result_type;
    ResultType operator()(ArgumentType1 const& x,
                          ArgumentType2 const& y) const {
      return ResultType(%s);
    }
  };""" % (name, op)

def generate_in_place_binary(f, name, op):
  print >> f, """
  template <typename ArgumentType1,
            typename ArgumentType2>
  struct functor_%s {
    ArgumentType1& operator()(ArgumentType1& x,
                              ArgumentType2 const& y) const {
      %s;
      return x;
    }
  };""" % (name, op)

def run(target_dir):
  f = utils.join_open(target_dir, "detail/operator_functors.h", "w")
  utils.write_this_is_auto_generated(f, this)
  print >> f, """\
#ifndef SCITBX_ARRAY_FAMILY_OPERATOR_FUNCTORS_H
#define SCITBX_ARRAY_FAMILY_OPERATOR_FUNCTORS_H

namespace scitbx { namespace fn {"""

  for op, ftor_name in operator_functor_info.unary_functors.items():
    generate_unary(f, ftor_name, op + "x")
  for op, ftor_name in operator_functor_info.binary_functors.items():
    generate_binary(f, ftor_name, "x " + op + " y")
  for op, ftor_name in operator_functor_info.in_place_binary_functors.items():
    generate_in_place_binary(f, ftor_name, "x " + op + " y")

  print >> f, """
}} // namespace scitbx::fn

#endif // SCITBX_ARRAY_FAMILY_OPERATOR_FUNCTORS_H"""
  f.close()

if (__name__ == "__main__"):
  run(".")
