#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_arg.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/make_constructor.hpp>
#include <cctbx/xray/scatterer.h>
#include <scitbx/array_family/boost_python/shared_wrapper.h>

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
        .def("implies", &w_t::implies)
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
        .def("use_fp_fdp", &w_t::use_fp_fdp)
        .def("use_u_iso_only", &w_t::use_u_iso_only)
        .def("use_u_aniso_only", &w_t::use_u_aniso_only)
        .def("set_use", &w_t::set_use, (arg("state")), return_self<>())
        .def("set_use_u_iso", &w_t::set_use_u_iso, (
          arg("state")), return_self<>())
        .def("set_use_u_aniso", &w_t::set_use_u_aniso, (
          arg("state")), return_self<>())
        .def("set_grad_site", &w_t::set_grad_site, (
          arg("state")), return_self<>())
        .def("set_grad_u_iso", &w_t::set_grad_u_iso, (
          arg("state")), return_self<>())
        .def("set_grad_u_aniso", &w_t::set_grad_u_aniso, (
          arg("state")), return_self<>())
        .def("set_grad_occupancy", &w_t::set_grad_occupancy, (
          arg("state")), return_self<>())
        .def("set_grad_fp", &w_t::set_grad_fp, (
          arg("state")), return_self<>())
        .def("set_grad_fdp", &w_t::set_grad_fdp, (
          arg("state")), return_self<>())
        .def("set_curv_site_site", &w_t::set_curv_site_site, (
          arg("state")), return_self<>())
        .def("set_curv_site_u_iso", &w_t::set_curv_site_u_iso, (
          arg("state")), return_self<>())
        .def("set_curv_site_u_aniso", &w_t::set_curv_site_u_aniso, (
          arg("state")), return_self<>())
        .def("set_curv_site_occupancy", &w_t::set_curv_site_occupancy, (
          arg("state")), return_self<>())
        .def("set_curv_site_fp", &w_t::set_curv_site_fp, (
          arg("state")), return_self<>())
        .def("set_curv_site_fdp", &w_t::set_curv_site_fdp, (
          arg("state")), return_self<>())
        .def("set_curv_u_iso_u_iso", &w_t::set_curv_u_iso_u_iso, (
          arg("state")), return_self<>())
        .def("set_curv_u_iso_u_aniso", &w_t::set_curv_u_iso_u_aniso, (
          arg("state")), return_self<>())
        .def("set_curv_u_iso_occupancy", &w_t::set_curv_u_iso_occupancy, (
          arg("state")), return_self<>())
        .def("set_curv_u_iso_fp", &w_t::set_curv_u_iso_fp, (
          arg("state")), return_self<>())
        .def("set_curv_u_iso_fdp", &w_t::set_curv_u_iso_fdp, (
          arg("state")), return_self<>())
        .def("set_curv_u_aniso_u_aniso", &w_t::set_curv_u_aniso_u_aniso, (
          arg("state")), return_self<>())
        .def("set_curv_u_aniso_occupancy", &w_t::set_curv_u_aniso_occupancy, (
          arg("state")), return_self<>())
        .def("set_curv_u_aniso_fp", &w_t::set_curv_u_aniso_fp, (
          arg("state")), return_self<>())
        .def("set_curv_u_aniso_fdp", &w_t::set_curv_u_aniso_fdp, (
          arg("state")), return_self<>())
        .def("set_curv_occupancy_occupancy",
          &w_t::set_curv_occupancy_occupancy, (
            arg("state")), return_self<>())
        .def("set_curv_occupancy_fp", &w_t::set_curv_occupancy_fp, (
          arg("state")), return_self<>())
        .def("set_curv_occupancy_fdp", &w_t::set_curv_occupancy_fdp, (
          arg("state")), return_self<>())
        .def("set_curv_fp_fp", &w_t::set_curv_fp_fp, (
          arg("state")), return_self<>())
        .def("set_curv_fp_fdp", &w_t::set_curv_fp_fdp, (
          arg("state")), return_self<>())
        .def("set_curv_fdp_fdp", &w_t::set_curv_fdp_fdp, (
          arg("state")), return_self<>())
        .def("set_tan_u_iso", &w_t::set_tan_u_iso, (
          arg("state")), return_self<>())
        .def("set_use_fp_fdp", &w_t::set_use_fp_fdp, arg("state"),
             return_self<>())
        .def("set_use_u_iso_only", &w_t::set_use_u_iso_only)
        .def("set_use_u_aniso_only", &w_t::set_use_u_aniso_only)
        .def("set_use_u", &w_t::set_use_u, (arg("iso"), arg("aniso")))
        .def("set_grads", &w_t::set_grads, (arg("state")))
        .def_readonly("bits", &w_t::bits)
        .def_readwrite("param", &w_t::param)
      ;
    }
  };

  struct scatterer_flags_array_wrappers
  {
    typedef af::shared<scatterer_flags> wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt> wrapper = scitbx::af::boost_python::shared_wrapper<
        scatterer_flags,
        return_internal_reference<> >::wrap("shared_scatterer_flags");
      wrapper
        .def("__init__",
             make_constructor(from_scatterers,
                              default_call_policies(),
                              arg("scatterers")))
        .def("n_parameters", n_parameters)
        .def("assign_to", assign_to, return_self<>())
        ;
    }

    static wt*
    from_scatterers(af::const_ref<scatterer<> > const & scatterers)
    {
      wt *result = new wt();
      result->reserve(scatterers.size());
      for(std::size_t i=0; i < scatterers.size(); ++i) {
        result->push_back(scatterers[i].flags);
      }
      return result;
    }

    static std::size_t n_parameters(wt const &self) {
      return grad_flags_counts(self.const_ref()).n_parameters();
    }

    static
    wt const&
    assign_to(
      wt const &self, af::ref<scatterer<> > const &scatterers)
    {
      CCTBX_ASSERT(self.size() == scatterers.size());
      for (std::size_t i=0; i < self.size(); ++i) {
        scatterers[i].flags = self[i];
      }
      return self;
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
        .def_readonly("use_u_iso", &w_t::use_u_iso)
        .def_readonly("use_u_aniso", &w_t::use_u_aniso)
        .def("n_parameters", &w_t::n_parameters)
      ;
    }
  };

} // namespace <anoymous>

  void wrap_scatterer_flags()
  {
    using namespace boost::python;
    scatterer_flags_wrappers::wrap();
    scatterer_flags_array_wrappers::wrap();
    scatterer_grad_flags_counts_wrappers::wrap();
    def("set_scatterer_grad_flags",
        set_scatterer_grad_flags<scatterer<> >,
         (arg("scatterers"),
          arg("site")=false,
          arg("u_iso")=false,
          arg("u_aniso")=false,
          arg("occupancy")=false,
          arg("fp")=false,
          arg("fdp")=false,
          arg("tan_u_iso")=false,
          arg("param")=0));
  }

}}} // namespace cctbx::xray::boost_python
