#include <cctbx/boost_python/flex_fwd.h>
#include <cctbx/geometry_restraints/shared_wrapper_pickle.hpp>

#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/return_by_value.hpp>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <scitbx/boost_python/is_polymorphic_workaround.h>
#include <cctbx/geometry_restraints/chirality.h>
#include <cctbx/geometry_restraints/proxy_select.h>

SCITBX_BOOST_IS_POLYMORPHIC_WORKAROUND(cctbx::geometry_restraints::chirality)

namespace cctbx { namespace geometry_restraints {
namespace {

  struct chirality_proxy_wrappers : boost::python::pickle_suite
  {
    typedef chirality_proxy w_t;

    static boost::python::tuple
      getinitargs(w_t const& self)
    {
      return boost::python::make_tuple(self.i_seqs,
        self.sym_ops,
        self.volume_ideal,
        self.both_signs,
        self.weight,
        self.origin_id
        );
    }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("chirality_proxy", no_init)
        .def(init<af::tiny<unsigned, 4> const&,
          optional_container<af::shared<sgtbx::rt_mx> > const&,
          double, bool, double,
          unsigned char>((
          arg("i_seqs"),
          arg("sym_ops"),
          arg("volume_ideal"),
          arg("both_signs"),
          arg("weight"),
          arg("origin_id")=0)))
        .def(init<af::tiny<unsigned, 4> const&, double, bool, double,
          unsigned char>((
          arg("i_seqs"),
          arg("volume_ideal"),
          arg("both_signs"),
          arg("weight"),
          arg("origin_id")=0)))
        .def(init<af::tiny<unsigned, 4> const&, w_t const&>((
          arg("i_seqs"),
          arg("proxy"))))
        .def("scale_weight", &w_t::scale_weight, (arg("factor")))
        .def("sort_i_seqs", &w_t::sort_i_seqs)
        .add_property("i_seqs", make_getter(&w_t::i_seqs, rbv()))
        .add_property("sym_ops", make_getter(&w_t::sym_ops, rbv()))
        .def_readonly("volume_ideal", &w_t::volume_ideal)
        .def_readonly("both_signs", &w_t::both_signs)
        .def_readwrite("weight", &w_t::weight)
        .def_readwrite("origin_id", &w_t::origin_id)
        .def_pickle(chirality_proxy_wrappers())
        ;
      {
        typedef scitbx::af::boost_python::shared_wrapper<w_t> shared_w_t;
        scitbx::af::boost_python::shared_wrapper<w_t>::wrap(
          "shared_chirality_proxy")
          .def("proxy_select",
            (af::shared<w_t>(*)(
              af::const_ref<w_t> const&,
              std::size_t,
              af::const_ref<std::size_t> const&))
                shared_proxy_select, (
            arg("n_seq"), arg("iselection")))
          .def("proxy_select",
            (af::shared<w_t>(*)(
              af::const_ref<w_t> const&,
              unsigned char))
                shared_proxy_select_origin, (
            arg("origin_id")))
          .def("proxy_remove",
            (af::shared<w_t>(*)(
              af::const_ref<w_t> const&,
              af::const_ref<bool> const&))
                shared_proxy_remove, (
            arg("selection")))
          .def("proxy_remove",
            (af::shared<w_t>(*)(
              af::const_ref<w_t> const&,
              unsigned char))
                shared_proxy_remove, (
            arg("origin_id")))
          .def_pickle(shared_wrapper_pickle_suite< shared_w_t::w_t >())
        ;
      }
    }
  };

  struct chirality_wrappers : boost::python::pickle_suite
  {
    typedef chirality w_t;

    static boost::python::tuple
      getinitargs(w_t const& self)
    {
      return boost::python::make_tuple(self.sites,
        self.volume_ideal,
        self.both_signs,
        self.weight
        );
    }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("chirality", no_init)
        .def(init<af::tiny<scitbx::vec3<double>, 4> const&,
                  double, bool, double>(
          (arg("sites"),
           arg("volume_ideal"),
           arg("both_signs"),
           arg("weight"))))
        .def(init<uctbx::unit_cell const&,
            af::const_ref<scitbx::vec3<double> > const&,
            chirality_proxy const&>((
            arg("unit_cell"), arg("sites_cart"), arg("proxy"))))
        .def(init<af::const_ref<scitbx::vec3<double> > const&,
                  chirality_proxy const&>(
          (arg("sites_cart"), arg("proxy"))))
        .add_property("sites", make_getter(&w_t::sites, rbv()))
        .add_property("sym_ops", make_getter(&w_t::sym_ops, rbv()))
        .def_readonly("volume_ideal", &w_t::volume_ideal)
        .def_readonly("both_signs", &w_t::both_signs)
        .def_readonly("weight", &w_t::weight)
        .def_readonly("volume_model", &w_t::volume_model)
        .def_readonly("delta_sign", &w_t::delta_sign)
        .def_readonly("delta", &w_t::delta)
        .def("residual", &w_t::residual)
        .def("gradients", &w_t::gradients)
        .def_pickle(chirality_wrappers())
        ;
    }
  };

  void
  wrap_all()
  {
    using namespace boost::python;
    chirality_proxy_wrappers::wrap();
    chirality_wrappers::wrap();
    def("chirality_deltas",
      (af::shared<double>(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<chirality_proxy> const&))
        chirality_deltas,
      (arg("sites_cart"), arg("proxies")));
    def("chirality_residuals",
      (af::shared<double>(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<chirality_proxy> const&))
        chirality_residuals,
      (arg("sites_cart"), arg("proxies")));
    def("chirality_residual_sum",
      (double(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<chirality_proxy> const&,
        af::ref<scitbx::vec3<double> > const&))
        chirality_residual_sum,
      (arg("sites_cart"), arg("proxies"), arg("gradient_array")));
    def("chirality_deltas",
      (af::shared<double>(*)(
        uctbx::unit_cell const&,
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<chirality_proxy> const&))
        chirality_deltas,
      (arg("unit_cell"), arg("sites_cart"), arg("proxies")));
    def("chirality_residuals",
      (af::shared<double>(*)(
        uctbx::unit_cell const&,
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<chirality_proxy> const&))
        chirality_residuals,
      (arg("unit_cell"), arg("sites_cart"), arg("proxies")));
    def("chirality_residual_sum",
      (double(*)(
        uctbx::unit_cell const&,
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<chirality_proxy> const&,
        af::ref<scitbx::vec3<double> > const&))
        chirality_residual_sum,
      (arg("unit_cell"), arg("sites_cart"), arg("proxies"), arg("gradient_array")));
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_chirality() { wrap_all(); }

}}} // namespace cctbx::geometry_restraints::boost_python
