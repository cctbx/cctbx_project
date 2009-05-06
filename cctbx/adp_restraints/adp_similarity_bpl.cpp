#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/return_by_value.hpp>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <cctbx/adp_restraints/adp_similarity.h>
#include <scitbx/boost_python/container_conversions.h>


namespace cctbx { namespace adp_restraints {

namespace {


  struct adp_similarity_proxy_wrappers
  {
    typedef adp_similarity_proxy w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("adp_similarity_proxy", no_init)
        .def(init<af::tiny<unsigned, 2> const&, double>((
           arg_("i_seqs"),
           arg_("weight"))))
        .add_property("i_seqs", make_getter(&w_t::i_seqs, rbv()))
        .add_property("weight", make_getter(&w_t::weight, rbv()))
      ;
      {
        scitbx::af::boost_python::shared_wrapper<w_t>::wrap(
          "shared_adp_similarity_proxy")
        ;
      }
    }
  };

  struct adp_similarity_wrappers
  {
    typedef adp_similarity w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("adp_similarity", no_init)
        .def(init<
           af::tiny<scitbx::sym_mat3<double>, 2> const&,
           af::tiny<double, 2> const&,
           af::tiny<bool, 2> const&,
           double>(
          (arg_("u_cart"),
           arg_("u_iso"),
           arg_("use_u_aniso"),
           arg_("weight"))))
        .def(init<
           af::const_ref<scitbx::sym_mat3<double> > const&,
           af::const_ref<double> const&,
           af::const_ref<bool> const&,
           adp_similarity_proxy const&>(
          (arg_("u_cart"),
           arg_("u_iso"),
           arg_("use_u_aniso"),
           arg_("proxy"))))
        .add_property("u_cart", make_getter(&w_t::u_cart, rbv()))
        .add_property("u_iso", make_getter(&w_t::u_iso, rbv()))
        .add_property("use_u_aniso", make_getter(&w_t::use_u_aniso, rbv()))
        .add_property("weight", make_getter(&w_t::weight, rbv()))
        .def("deltas", &w_t::deltas)
        .def("rms_deltas", &w_t::rms_deltas)
        .def("residual", &w_t::residual)
        .def("gradients", &w_t::gradients)
      ;
    }
  };

  void
  wrap_all()
  {
    using namespace scitbx::boost_python::container_conversions;
    tuple_mapping_fixed_size<scitbx::af::tiny<scitbx::sym_mat3<double>, 2> >();
    tuple_mapping_fixed_size<scitbx::af::tiny<bool, 2> >();

    using namespace boost::python;
    adp_similarity_wrappers::wrap();
    adp_similarity_proxy_wrappers::wrap();
    def("adp_similarity_residual_sum", adp_similarity_residual_sum,
      (arg_("u_cart"),
       arg_("u_iso"),
       arg_("use_u_aniso"),
       arg_("proxies"),
       arg_("gradients_aniso_cart"),
       arg_("gradients_iso")));
    def("adp_similarity_residuals", adp_similarity_residuals,
      (arg_("u_cart"),
       arg_("u_iso"),
       arg_("use_u_aniso"),
       arg_("proxies")));
    def("adp_similarity_deltas_rms", adp_similarity_deltas_rms,
      (arg_("u_cart"),
       arg_("u_iso"),
       arg_("use_u_aniso"),
       arg_("proxies")));
  }

}

namespace boost_python {

  void
  wrap_adp_similarity() { wrap_all(); }

}}}
