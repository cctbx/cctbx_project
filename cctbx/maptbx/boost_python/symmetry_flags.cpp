#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/maptbx/symmetry_flags.h>
#include <boost/python/class.hpp>

namespace cctbx { namespace maptbx { namespace boost_python {

namespace {

  struct symmetry_flags_wrappers
  {
    typedef symmetry_flags w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      class_<w_t>("symmetry_flags", no_init)
        .def(init<bool, optional<bool, bool> >(
          (arg_("use_space_group_symmetry"),
           arg_("use_normalizer_k2l")=false,
           arg_("use_structure_seminvariants")=false)))
        .def("use_space_group_symmetry", &w_t::use_space_group_symmetry)
        .def("use_normalizer_k2l", &w_t::use_normalizer_k2l)
        .def("use_structure_seminvariants", &w_t::use_structure_seminvariants)
        .def("select_sub_space_group", &w_t::select_sub_space_group,
          (arg_("space_group_type")))
        .def("__eq__", &w_t::operator==)
        .def("__ne__", &w_t::operator!=)
      ;
    }
  };

} // namespace <anoymous>

  void wrap_symmetry_flags()
  {
    symmetry_flags_wrappers::wrap();
  }

}}} // namespace cctbx::maptbx::boost_python
