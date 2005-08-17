#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <cctbx/sgtbx/brick.h>

namespace cctbx { namespace sgtbx { namespace boost_python {

namespace {

  struct brick_wrappers
  {
    typedef brick w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("brick", no_init)
        .def(init<space_group_type const&>((arg_("space_group_type"))))
        .def("as_string", &w_t::as_string)
        .def("__str__", &w_t::as_string)
        .def("is_inside",
          (bool(w_t::*)(tr_vec const&) const) &w_t::is_inside, (
            arg_("point")))
      ;
    }
  };

} // namespace <anoymous>

  void wrap_brick()
  {
    brick_wrappers::wrap();
  }

}}} // namespace cctbx::sgtbx::boost_python
