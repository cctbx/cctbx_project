#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <cctbx/xray/refinement_flags.h>

namespace cctbx { namespace xray { namespace boost_python {

namespace {

  struct refinement_flags_wrappers
  {
    typedef refinement_flags w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("refinement_flags", no_init)
        .def(init<optional<
          bool, bool, bool, bool, bool, bool, bool, float> >((
            arg_("site")=false, arg_("u_iso")=false, arg_("u_aniso")=false,
            arg_("occupancy")=false, arg_("fp")=false, arg_("fdp")=false,
            arg_("tan_u_iso")=false, arg_("param")=0)))
        .def(init<w_t const&>()) // copy constructor
        .def("site", &w_t::site)
        .def("set_site", &w_t::set_site, (arg_("state")))
        .def("u_iso", &w_t::u_iso)
        .def("set_u_iso", &w_t::set_u_iso, (arg_("state")))
        .def("u_aniso", &w_t::u_aniso)
        .def("set_u_aniso", &w_t::set_u_aniso, (arg_("state")))
        .def("occupancy", &w_t::occupancy)
        .def("set_occupancy", &w_t::set_occupancy, (arg_("state")))
        .def("fp", &w_t::fp)
        .def("set_fp", &w_t::set_fp, (arg_("state")))
        .def("fdp", &w_t::fdp)
        .def("set_fdp", &w_t::set_fdp, (arg_("state")))
        .def("tan_u_iso", &w_t::tan_u_iso)
        .def("set_tan_u_iso", &w_t::set_tan_u_iso, (arg_("state")))
        .def_readonly("bits", &w_t::bits)
        .def_readwrite("param", &w_t::param)
      ;
    }
  };

} // namespace <anoymous>

  void wrap_refinement_flags()
  {
    refinement_flags_wrappers::wrap();
  }

}}} // namespace cctbx::xray::boost_python
