#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <cctbx/adp_restraints/adp_similarity.h>
#include <scitbx/boost_python/container_conversions.h>


namespace cctbx { namespace adp_restraints {

namespace {

  struct adp_similarity_proxy_wrappers {
    typedef adp_similarity_proxy w_t;

    static void wrap() {
      using namespace boost::python;
      class_<w_t, bases<adp_restraint_proxy<2> > >
        ("adp_similarity_proxy", no_init)
        .def(init<af::tiny<unsigned, 2> const&, double>((
           arg("i_seqs"),
           arg("weight"))))
      ;
      {
        scitbx::af::boost_python::shared_wrapper<w_t>::wrap(
          "shared_adp_similarity_proxy")
        ;
      }
    }
  };

  struct adp_similarity_wrappers {
    typedef adp_similarity w_t;

    static void wrap() {
      using namespace boost::python;
      class_<w_t, bases<adp_restraint_base<2> > >
        ("adp_similarity", no_init)
        .def(init<
           af::tiny<scitbx::sym_mat3<double>, 2> const&,
           double>(
          (arg("u_cart"),
           arg("weight"))))
        .def(init<
           af::tiny<double, 2> const&,
           double>(
          (arg("u_iso"),
           arg("weight"))))
        .def(init<
           scitbx::sym_mat3<double> const&,
           double,
           double>(
          (arg("u_cart"),
           arg("u_iso"),
           arg("weight"))))
        .def(init<
            double,
            scitbx::sym_mat3<double> const&,
           double>(
          (arg("u_iso"),
           arg("u_cart"),
           arg("weight"))))
        .def(init<
            adp_restraint_params<double> const&,
            adp_similarity_proxy const&>(
          (arg("params"),
           arg("proxy"))))
        .def("gradients2", &w_t::gradients2)
      ;
    }
  };

  struct adp_u_eq_similarity_proxy_wrappers {
    typedef adp_u_eq_similarity_proxy w_t;

    static void wrap() {
      using namespace boost::python;
      class_<w_t, bases<adp_restraint_proxy<2> > >
        ("adp_u_eq_similarity_proxy", no_init)
        .def(init<af::tiny<unsigned, 2> const&, double>((
           arg("i_seqs"),
           arg("weight"))))
      ;
      {
        scitbx::af::boost_python::shared_wrapper<w_t>::wrap(
          "shared_adp_u_eq_similarity_proxy")
        ;
      }
    }
  };

  struct adp_u_eq_similarity_wrappers {
    typedef adp_u_eq_similarity w_t;

    static void wrap() {
      using namespace boost::python;
      class_<w_t, bases<adp_restraint_base<2> > >
        ("adp_u_eq_similarity", no_init)
        .def(init<
           af::tiny<scitbx::sym_mat3<double>, 2> const&,
           double>(
          (arg("u_cart"),
           arg("weight"))))
        .def(init<
           af::tiny<double, 2> const&,
           double>(
          (arg("u_iso"),
           arg("weight"))))
        .def(init<
           scitbx::sym_mat3<double> const&,
           double,
           double>(
          (arg("u_cart"),
           arg("u_iso"),
           arg("weight"))))
        .def(init<
            double,
            scitbx::sym_mat3<double> const&,
           double>(
          (arg("u_iso"),
           arg("u_cart"),
           arg("weight"))))
        .def(init<
            adp_restraint_params<double> const&,
            adp_u_eq_similarity_proxy const&>(
          (arg("params"),
           arg("proxy"))))
        .def("gradients2", &w_t::gradients2)
      ;
    }
  };

  struct adp_volume_similarity_proxy_wrappers {
    typedef adp_volume_similarity_proxy w_t;

    static void wrap() {
      using namespace boost::python;
      class_<w_t, bases<adp_restraint_proxy<2> > >
        ("adp_volume_similarity_proxy", no_init)
        .def(init<af::tiny<unsigned, 2> const&, double>((
           arg("i_seqs"),
           arg("weight"))))
      ;
      {
        scitbx::af::boost_python::shared_wrapper<w_t>::wrap(
          "shared_adp_volume_similarity_proxy")
        ;
      }
    }
  };

  struct adp_volume_similarity_wrappers {
    typedef adp_volume_similarity w_t;

    static void wrap() {
      using namespace boost::python;
      class_<w_t, bases<adp_restraint_base<2> > >
        ("adp_volume_similarity", no_init)
        .def(init<
           af::tiny<scitbx::sym_mat3<double>, 2> const&,
           double>(
          (arg("u_cart"),
           arg("weight"))))
        .def(init<
           af::tiny<double, 2> const&,
           double>(
          (arg("u_iso"),
           arg("weight"))))
        .def(init<
           scitbx::sym_mat3<double> const&,
           double,
           double>(
          (arg("u_cart"),
           arg("u_iso"),
           arg("weight"))))
        .def(init<
            double,
            scitbx::sym_mat3<double> const&,
           double>(
          (arg("u_iso"),
           arg("u_cart"),
           arg("weight"))))
        .def(init<
            adp_restraint_params<double> const&,
            adp_volume_similarity_proxy const&>(
          (arg("params"),
           arg("proxy"))))
        .def("gradients2", &w_t::gradients2)
      ;
    }
  };

  void wrap_all() {
    using namespace boost::python;
    adp_similarity_wrappers::wrap();
    adp_similarity_proxy_wrappers::wrap();
    adp_u_eq_similarity_wrappers::wrap();
    adp_u_eq_similarity_proxy_wrappers::wrap();
    adp_volume_similarity_wrappers::wrap();
    adp_volume_similarity_proxy_wrappers::wrap();
  }

}

namespace boost_python {

  void wrap_adp_similarity() { wrap_all(); }

}}}
