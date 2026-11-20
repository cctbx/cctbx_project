#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <cctbx/adp_restraints/npd_adp.h>
#include <scitbx/boost_python/container_conversions.h>


namespace cctbx { namespace adp_restraints {

namespace {

  struct npd_adp_proxy_wrappers
  {
    typedef npd_adp_proxy w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t, bases<adp_restraint_proxy<1> > >
        ("npd_adp_proxy", no_init)
        .def(init<
           af::tiny<unsigned, 1> const &,
           double,
           double>(
          (arg("i_seqs"),
           arg("weight"),
           arg("s_target"))))
        .def_readwrite("s_target", &w_t::s_target)
      ;
      {
        scitbx::af::boost_python::shared_wrapper<w_t>::wrap(
          "shared_npd_adp_proxy")
        ;
      }
    }
  };

  struct npd_adp_wrappers
  {
    typedef npd_adp w_t;

    static void
    wrap() {
      using namespace boost::python;
      class_<w_t>("npd_adp", no_init)
        .def(init<
            scitbx::sym_mat3<double> const &,
            double,
            double>(
          (arg("u_cart"),
           arg("weight"),
           arg("s_target"))))
        .def(init<
            adp_restraint_params<double> const &,
            npd_adp_proxy const &>(
          (arg("params"),
           arg("proxy"))))
        .def("active", &w_t::active)
        .def("s", &w_t::s)
        .def("delta", &w_t::delta)
        .def("residual", &w_t::residual)
        .def("rms_deltas", &w_t::rms_deltas)
        .def("gradients", &w_t::gradients)
      ;
    }
  };

  void wrap_all() {
    npd_adp_wrappers::wrap();
    npd_adp_proxy_wrappers::wrap();
  }

}

namespace boost_python {

  void
  wrap_npd_adp() { wrap_all(); }

}

}} // cctbx::adp_restraints
