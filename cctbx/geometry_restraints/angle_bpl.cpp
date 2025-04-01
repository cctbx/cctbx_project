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
#include <cctbx/geometry_restraints/angle.h>
#include <cctbx/geometry_restraints/proxy_select.h>

SCITBX_BOOST_IS_POLYMORPHIC_WORKAROUND(cctbx::geometry_restraints::angle)

namespace cctbx { namespace geometry_restraints {
namespace {

  struct angle_proxy_wrappers : boost::python::pickle_suite
  {
    typedef angle_proxy w_t;

    static boost::python::tuple
      getinitargs(w_t const& self)
    {
      return boost::python::make_tuple(self.i_seqs,
        self.angle_ideal,
        self.weight,
        self.slack,
        self.origin_id
        );
    }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("angle_proxy", no_init)
        .def(init<
          af::tiny<unsigned, 3> const&, double, double,
          double, unsigned char>((
          arg("i_seqs"), arg("angle_ideal"), arg("weight"),
          arg("slack")=0, arg("origin_id")=0)))
        .def(init<
          af::tiny<unsigned, 3> const&,
          optional_container<af::shared<sgtbx::rt_mx> > const&,
          double,
          double,
          double,
          unsigned char>((
            arg("i_seqs"),
            arg("sym_ops"),
            arg("angle_ideal"),
            arg("weight"),
            arg("slack")=0,
            arg("origin_id")=0)))
        .def(init<af::tiny<unsigned, 3> const&, w_t const&>((
          arg("i_seqs"), arg("proxy"))))
        .def("scale_weight", &w_t::scale_weight, (arg("factor")))
        .def("sort_i_seqs", &w_t::sort_i_seqs)
        .add_property("i_seqs", make_getter(&w_t::i_seqs, rbv()))
        .add_property("sym_ops", make_getter(&w_t::sym_ops, rbv()))
        .def_readwrite("angle_ideal", &w_t::angle_ideal)
        .def_readwrite("weight", &w_t::weight)
        .def_readwrite("slack", &w_t::slack)
        .def_readwrite("origin_id", &w_t::origin_id)
        .def_pickle(angle_proxy_wrappers())
        ;
      {
        typedef return_internal_reference<> rir;
        typedef scitbx::af::boost_python::shared_wrapper<w_t, rir> shared_w_t;
        scitbx::af::boost_python::shared_wrapper<w_t, rir>::wrap(
          "shared_angle_proxy")
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

  struct angle_wrappers : boost::python::pickle_suite
  {
    typedef angle w_t;

    static boost::python::tuple
      getinitargs(w_t const& self)
    {
      return boost::python::make_tuple(self.sites,
        self.angle_ideal,
        self.weight,
        self.slack,
        self.origin_id
        );
    }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t, bases<cctbx::geometry::angle<double> > >("angle", no_init)
        .def(init<af::tiny<scitbx::vec3<double>, 3> const&,
          double,
          double,
          double,
          unsigned char>(
          (arg("sites"), arg("angle_ideal"), arg("weight"),
          arg("slack")=0, arg("origin_id")=0)))
        .def(init<af::const_ref<scitbx::vec3<double> > const&,
                  angle_proxy const&>(
          (arg("sites_cart"), arg("proxy"))))
        .def(init<uctbx::unit_cell const&,
                  af::const_ref<scitbx::vec3<double> > const&,
                  angle_proxy const&>(
          (arg("unit_cell"), arg("sites_cart"), arg("proxy"))))
        .add_property("sites", make_getter(&w_t::sites, rbv()))
        .def_readonly("angle_ideal", &w_t::angle_ideal)
        .def_readonly("weight", &w_t::weight)
        .def_readonly("slack", &w_t::slack)
        .def_readonly("have_angle_model", &w_t::have_angle_model)
        .def_readonly("angle_model", &w_t::angle_model)
        .def_readonly("delta", &w_t::delta)
        .def_readonly("delta_slack", &w_t::delta_slack)
        .def_readonly("origin_id", &w_t::origin_id)
        .def("residual", &w_t::residual)
        .def("grads_and_curvs", &w_t::grads_and_curvs, (
          arg("epsilon")=1e-100))
        .def("gradients", &w_t::gradients, (
          arg("epsilon")=1e-100))
        .def_pickle(angle_wrappers())
        ;
    }
  };

  void
  wrap_all()
  {
    using namespace boost::python;
    angle_proxy_wrappers::wrap();
    angle_wrappers::wrap();
    def("angle_deltas",
      (af::shared<double>(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<angle_proxy> const&))
      angle_deltas,
      (arg("sites_cart"), arg("proxies")));
    def("angle_deltas",
      (af::shared<double>(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<angle_proxy> const&,
        unsigned char))
      angle_deltas,
      (arg("sites_cart"), arg("proxies"), arg("origin_id")));
    def("angle_residuals",
      (af::shared<double>(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<angle_proxy> const&))
      angle_residuals,
      (arg("sites_cart"), arg("proxies")));
    def("angle_residual_sum",
      (double(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<angle_proxy> const&,
        af::ref<scitbx::vec3<double> > const&))
      angle_residual_sum,
      (arg("sites_cart"), arg("proxies"), arg("gradient_array")));
    def("angle_deltas",
      (af::shared<double>(*)(
        uctbx::unit_cell const&,
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<angle_proxy> const&))
      angle_deltas,
      (arg("unit_cell"), arg("sites_cart"), arg("proxies")));
    def("angle_deltas",
      (af::shared<double>(*)(
        uctbx::unit_cell const&,
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<angle_proxy> const&,
        unsigned char))
      angle_deltas,
      (arg("unit_cell"), arg("sites_cart"), arg("proxies"), arg("origin_id")));
    def("angle_residuals",
      (af::shared<double>(*)(
        uctbx::unit_cell const&,
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<angle_proxy> const&))
      angle_residuals,
      (arg("unit_cell"), arg("sites_cart"), arg("proxies")));
    def("angle_residual_sum",
      (double(*)(
        uctbx::unit_cell const&,
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<angle_proxy> const&,
        af::ref<scitbx::vec3<double> > const&))
      angle_residual_sum,
      (arg("unit_cell"), arg("sites_cart"), arg("proxies"), arg("gradient_array")));
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_angle() { wrap_all(); }

}}} // namespace cctbx::geometry_restraints::boost_python
