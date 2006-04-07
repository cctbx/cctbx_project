#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_arg.hpp>
#include <cctbx/xray/scatterer.h>

namespace cctbx { namespace xray { namespace boost_python {

namespace {

  struct scatterer_flags_wrappers
  {
    typedef scatterer_flags w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("scatterer_flags", no_init)
        .def(init<>())
        .def(init<w_t const&>()) // copy constructor
        .def("use", &w_t::use)
        .def("use_u_iso", &w_t::use_u_iso)
        .def("use_u_aniso", &w_t::use_u_aniso)
        .def("grad_site", &w_t::grad_site)
        .def("grad_u_iso", &w_t::grad_u_iso)
        .def("grad_u_aniso", &w_t::grad_u_aniso)
        .def("grad_occupancy", &w_t::grad_occupancy)
        .def("grad_fp", &w_t::grad_fp)
        .def("grad_fdp", &w_t::grad_fdp)
        .def("curv_site_site", &w_t::curv_site_site)
        .def("curv_site_u_iso", &w_t::curv_site_u_iso)
        .def("curv_site_u_aniso", &w_t::curv_site_u_aniso)
        .def("curv_site_occupancy", &w_t::curv_site_occupancy)
        .def("curv_site_fp", &w_t::curv_site_fp)
        .def("curv_site_fdp", &w_t::curv_site_fdp)
        .def("curv_u_iso_u_iso", &w_t::curv_u_iso_u_iso)
        .def("curv_u_iso_u_aniso", &w_t::curv_u_iso_u_aniso)
        .def("curv_u_iso_occupancy", &w_t::curv_u_iso_occupancy)
        .def("curv_u_iso_fp", &w_t::curv_u_iso_fp)
        .def("curv_u_iso_fdp", &w_t::curv_u_iso_fdp)
        .def("curv_u_aniso_u_aniso", &w_t::curv_u_aniso_u_aniso)
        .def("curv_u_aniso_occupancy", &w_t::curv_u_aniso_occupancy)
        .def("curv_u_aniso_fp", &w_t::curv_u_aniso_fp)
        .def("curv_u_aniso_fdp", &w_t::curv_u_aniso_fdp)
        .def("curv_occupancy_occupancy", &w_t::curv_occupancy_occupancy)
        .def("curv_occupancy_fp", &w_t::curv_occupancy_fp)
        .def("curv_occupancy_fdp", &w_t::curv_occupancy_fdp)
        .def("curv_fp_fp", &w_t::curv_fp_fp)
        .def("curv_fp_fdp", &w_t::curv_fp_fdp)
        .def("curv_fdp_fdp", &w_t::curv_fdp_fdp)
        .def("tan_u_iso", &w_t::tan_u_iso)
        .def("set_use", &w_t::set_use, (arg_("state")), return_self<>())
        .def("set_use_u_iso", &w_t::set_use_u_iso, (arg_("state")), return_self<>())
        .def("set_use_u_aniso", &w_t::set_use_u_aniso, (arg_("state")), return_self<>())
        .def("set_grad_site", &w_t::set_grad_site, (arg_("state")), return_self<>())
        .def("set_grad_u_iso", &w_t::set_grad_u_iso, (arg_("state")), return_self<>())
        .def("set_grad_u_aniso", &w_t::set_grad_u_aniso, (arg_("state")), return_self<>())
        .def("set_grad_occupancy", &w_t::set_grad_occupancy, (arg_("state")), return_self<>())
        .def("set_grad_fp", &w_t::set_grad_fp, (arg_("state")), return_self<>())
        .def("set_grad_fdp", &w_t::set_grad_fdp, (arg_("state")), return_self<>())
        .def("set_curv_site_site", &w_t::set_curv_site_site, (arg_("state")), return_self<>())
        .def("set_curv_site_u_iso", &w_t::set_curv_site_u_iso, (arg_("state")), return_self<>())
        .def("set_curv_site_u_aniso", &w_t::set_curv_site_u_aniso, (arg_("state")), return_self<>())
        .def("set_curv_site_occupancy", &w_t::set_curv_site_occupancy, (arg_("state")), return_self<>())
        .def("set_curv_site_fp", &w_t::set_curv_site_fp, (arg_("state")), return_self<>())
        .def("set_curv_site_fdp", &w_t::set_curv_site_fdp, (arg_("state")), return_self<>())
        .def("set_curv_u_iso_u_iso", &w_t::set_curv_u_iso_u_iso, (arg_("state")), return_self<>())
        .def("set_curv_u_iso_u_aniso", &w_t::set_curv_u_iso_u_aniso, (arg_("state")), return_self<>())
        .def("set_curv_u_iso_occupancy", &w_t::set_curv_u_iso_occupancy, (arg_("state")), return_self<>())
        .def("set_curv_u_iso_fp", &w_t::set_curv_u_iso_fp, (arg_("state")), return_self<>())
        .def("set_curv_u_iso_fdp", &w_t::set_curv_u_iso_fdp, (arg_("state")), return_self<>())
        .def("set_curv_u_aniso_u_aniso", &w_t::set_curv_u_aniso_u_aniso, (arg_("state")), return_self<>())
        .def("set_curv_u_aniso_occupancy", &w_t::set_curv_u_aniso_occupancy, (arg_("state")), return_self<>())
        .def("set_curv_u_aniso_fp", &w_t::set_curv_u_aniso_fp, (arg_("state")), return_self<>())
        .def("set_curv_u_aniso_fdp", &w_t::set_curv_u_aniso_fdp, (arg_("state")), return_self<>())
        .def("set_curv_occupancy_occupancy", &w_t::set_curv_occupancy_occupancy, (arg_("state")), return_self<>())
        .def("set_curv_occupancy_fp", &w_t::set_curv_occupancy_fp, (arg_("state")), return_self<>())
        .def("set_curv_occupancy_fdp", &w_t::set_curv_occupancy_fdp, (arg_("state")), return_self<>())
        .def("set_curv_fp_fp", &w_t::set_curv_fp_fp, (arg_("state")), return_self<>())
        .def("set_curv_fp_fdp", &w_t::set_curv_fp_fdp, (arg_("state")), return_self<>())
        .def("set_curv_fdp_fdp", &w_t::set_curv_fdp_fdp, (arg_("state")), return_self<>())
        .def("set_tan_u_iso", &w_t::set_tan_u_iso, (arg_("state")), return_self<>())
        .def_readonly("bits", &w_t::bits)
        .def_readwrite("param", &w_t::param)
      ;
    }
  };

  struct scatterer_grad_flags_counts_wrappers
  {
    typedef scatterer_grad_flags_counts w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("scatterer_grad_flags_counts", no_init)
        .def(init<scitbx::af::const_ref<scatterer<> > const&>())
        .def_readonly("site", &w_t::site)
        .def_readonly("u_iso", &w_t::u_iso)
        .def_readonly("u_aniso", &w_t::u_aniso)
        .def_readonly("occupancy", &w_t::occupancy)
        .def_readonly("fp", &w_t::fp)
        .def_readonly("fdp", &w_t::fdp)
        .def_readonly("tan_u_iso", &w_t::tan_u_iso)
        .def("all_zero", &w_t::all_zero)
      ;
    }
  };

  BOOST_PYTHON_FUNCTION_OVERLOADS(set_scatterer_grad_flags_overloads,
    set_scatterer_grad_flags,1,9)

} // namespace <anoymous>

  void wrap_scatterer_flags()
  {
    using namespace boost::python;
    scatterer_flags_wrappers::wrap();
    scatterer_grad_flags_counts_wrappers::wrap();
    def("set_scatterer_grad_flags",(void
      (*)(scitbx::af::ref<scatterer<> > const&,
                           bool,
                           bool,
                           bool,
                           bool,
                           bool,
                           bool,
                           bool,
                           int)) set_scatterer_grad_flags<scatterer<> >,
         (arg_("scatterers"),
          arg_("site")=false,
          arg_("u_iso")=false,
          arg_("u_aniso")=false,
          arg_("occupancy")=false,
          arg_("fp")=false,
          arg_("fdp")=false,
          arg_("tan_u_iso")=false,
          arg_("param")=0));
  }

}}} // namespace cctbx::xray::boost_python
