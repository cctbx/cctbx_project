from scitbx.source_generators import utils

this = "scitbx.source_generators.array_family.generate_reductions"

def make_dict(**kwds):
  return kwds

def substitute(subs, template):
  for key, value in subs.items():
    template = template.replace("${"+key+"}", value)
  assert template.find("${") < 0, "Incomplete substitutions."
  return template[1:]

def generate_order(f, subs):
  print >> f, substitute(subs, """
  template <typename ElementType1${templ_decl_2_1},
            typename ElementType2${templ_decl_2_2}>
  int
  inline
  order(
    ${array_type_plain}<ElementType1${templ_inst_2_1}> const& a1,
    ${array_type_plain}<ElementType2${templ_inst_2_2}> const& a2)
  {
    return order(a1.const_ref(), a2.const_ref());
  }
""")

def generate_max_index_etc(f, subs):
  for func_name in ("max_index", "min_index"):
    subs["func_name"] = func_name
    print >> f, substitute(subs, """
  template <typename ElementType${templ_decl_2}>
  inline
  std::size_t
  ${func_name}(${array_type_plain}<ElementType${templ_inst_2}> const& a)
  {
    return ${func_name}(a.const_ref());
  }
""")

def generate_first_index_etc(f, subs):
  for func_name in ("first_index", "last_index"):
    subs["func_name"] = func_name
    print >> f, substitute(subs, """
  template <typename ElementType${templ_decl_2}, class PredicateType>
  inline
  boost::optional<std::size_t>
  ${func_name}(${array_type_plain}<ElementType${templ_inst_2}> const& a,
               PredicateType p)
  {
    return ${func_name}(a.const_ref(), p);
  }
""")

def generate_max_etc(f, subs):
  for func_name in ("max", "min", "max_absolute",
                    "sum", "sum_sq",
                    "product", "mean", "mean_sq"):
    subs["func_name"] = func_name
    print >> f, substitute(subs, """
  template <typename ElementType${templ_decl_2}>
  inline
  ElementType
  ${func_name}(${array_type_plain}<ElementType${templ_inst_2}> const& a)
  {
    return ${func_name}(a.const_ref());
  }
""")

def generate_mean_weighted_etc(f, subs):
  for func_name in ("mean_weighted", "mean_sq_weighted"):
    subs["func_name"] = func_name
    print >> f, substitute(subs, """
  template <typename ElementTypeValues${templ_decl_2_1eq},
            typename ElementTypeWeights${templ_decl_2_2eq}>
  inline
  ElementTypeValues
  ${func_name}(
    ${array_type_plain}<ElementTypeValues${templ_inst_2_1eq}> const& values,
    ${array_type_plain}<ElementTypeWeights${templ_inst_2_2eq}> const& weights)
  {
    return ${func_name}(values.const_ref(), weights.const_ref());
  }
""")

def one_type(target_dir, subs):
  array_type = subs["array_type"]
  subs["array_type_plain"] = array_type + "_plain"
  subs["ARRAY_TYPE"] = array_type.upper()
  f = utils.join_open(target_dir, "%s_reductions.h" % array_type, "w")
  utils.write_this_is_auto_generated(f, this)
  print >> f, substitute(subs, """
#ifndef SCITBX_ARRAY_FAMILY_${ARRAY_TYPE}_REDUCTIONS_H
#define SCITBX_ARRAY_FAMILY_${ARRAY_TYPE}_REDUCTIONS_H

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#include <scitbx/array_family/ref_reductions.h>
#include <scitbx/array_family/${array_type_plain}.h>

namespace scitbx { namespace af {
""")

  generate_order(f, subs)
  generate_max_index_etc(f, subs)
  generate_max_etc(f, subs)
  generate_mean_weighted_etc(f, subs)
  generate_first_index_etc(f, subs)

  print >> f, substitute(subs, """
}} // namespace scitbx::af

#endif // DOXYGEN_SHOULD_SKIP_THIS

#endif // SCITBX_ARRAY_FAMILY_${ARRAY_TYPE}_REDUCTIONS_H
""")

  f.close()

def run(target_dir):
  tiny_subs = make_dict(
    array_type="tiny",
    templ_decl_2=", std::size_t N",
    templ_decl_2_1=", std::size_t N1",
    templ_decl_2_2=", std::size_t N2",
    templ_decl_2_1eq=", std::size_t N",
    templ_decl_2_2eq="",
    templ_inst_2=", N",
    templ_inst_2_1=", N1",
    templ_inst_2_2=", N2",
    templ_inst_2_1eq=", N",
    templ_inst_2_2eq=", N")
  small_subs = make_dict(
    array_type="small",
    templ_decl_2=", std::size_t N",
    templ_decl_2_1=", std::size_t N1",
    templ_decl_2_2=", std::size_t N2",
    templ_decl_2_1eq=", std::size_t N1",
    templ_decl_2_2eq=", std::size_t N2",
    templ_inst_2=", N",
    templ_inst_2_1=", N1",
    templ_inst_2_2=", N2",
    templ_inst_2_1eq=", N1",
    templ_inst_2_2eq=", N2")
  shared_subs = make_dict(
    array_type="shared",
    templ_decl_2="",
    templ_decl_2_1="",
    templ_decl_2_2="",
    templ_decl_2_1eq="",
    templ_decl_2_2eq="",
    templ_inst_2="",
    templ_inst_2_1="",
    templ_inst_2_2="",
    templ_inst_2_1eq="",
    templ_inst_2_2eq="")
  versa_subs = make_dict(
    array_type="versa",
    templ_decl_2=", typename AccessorType",
    templ_decl_2_1=", typename AccessorType1",
    templ_decl_2_2=", typename AccessorType2",
    templ_decl_2_1eq=", typename AccessorType1",
    templ_decl_2_2eq=", typename AccessorType2",
    templ_inst_2=", AccessorType",
    templ_inst_2_1=", AccessorType1",
    templ_inst_2_2=", AccessorType2",
    templ_inst_2_1eq=", AccessorType1",
    templ_inst_2_2eq=", AccessorType2")
  for subs in (tiny_subs, small_subs, shared_subs, versa_subs):
    one_type(target_dir, subs)

if (__name__ == "__main__"):
  run(".")
