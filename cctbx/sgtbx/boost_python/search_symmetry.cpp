#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <cctbx/sgtbx/search_symmetry.h>

namespace cctbx { namespace sgtbx { namespace boost_python {

namespace {

  struct search_symmetry_flags_wrappers
  {
    typedef search_symmetry_flags w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      class_<w_t>("search_symmetry_flags", no_init)
        .def(init<bool, optional<int, bool, bool, bool> >(
          (arg_("use_space_group_symmetry"),
           arg_("use_space_group_ltr")=0,
           arg_("use_seminvariant")=false,
           arg_("use_normalizer_k2l")=false,
           arg_("use_normalizer_l2n")=false)))
        .def("use_space_group_symmetry", &w_t::use_space_group_symmetry)
        .def("use_space_group_ltr", &w_t::use_space_group_ltr)
        .def("use_seminvariant", &w_t::use_seminvariant)
        .def("use_normalizer_k2l", &w_t::use_normalizer_k2l)
        .def("use_normalizer_l2n", &w_t::use_normalizer_l2n)
        .def("__eq__", &w_t::operator==)
        .def("__ne__", &w_t::operator!=)
      ;
    }
  };

  struct search_symmetry_wrappers
  {
    typedef search_symmetry w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      typedef return_value_policy<copy_const_reference> ccr;
      typedef return_internal_reference<> rir;
      class_<w_t>("search_symmetry", no_init)
        .def(init<search_symmetry_flags const&,
                  space_group_type const&>(
          (arg_("flags"),
           arg_("space_group_type"))))
        .def(init<search_symmetry_flags const&,
                  space_group_type const&,
                  structure_seminvariant const&>(
          (arg_("flags"),
           arg_("space_group_type"),
           arg_("seminvariant"))))
        .def("flags", &w_t::flags, rir())
        .def("group", &w_t::group, rir())
        .def("continuous_shifts", &w_t::continuous_shifts, ccr())
      ;
    }
  };

} // namespace <anoymous>

  void wrap_search_symmetry()
  {
    search_symmetry_flags_wrappers::wrap();
    search_symmetry_wrappers::wrap();
  }

}}} // namespace cctbx::sgtbx::boost_python
