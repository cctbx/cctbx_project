#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <cctbx/translation_search/symmetry_flags.h>

namespace cctbx { namespace translation_search { namespace boost_python {

namespace {

  struct symmetry_flags_wrappers
  {
    typedef symmetry_flags w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      class_<w_t, bases<sgtbx::search_symmetry_flags> >(
        "symmetry_flags", no_init)
        .def(init<bool, bool>(
          (arg_("is_isotropic_search_model"), arg_("have_f_part"))))
        .def("is_isotropic_search_model", &w_t::is_isotropic_search_model)
        .def("have_f_part", &w_t::have_f_part)
      ;
    }
  };

} // namespace <anoymous>

  void wrap_symmetry_flags()
  {
    symmetry_flags_wrappers::wrap();
  }

}}} // namespace cctbx::translation_search::boost_python
