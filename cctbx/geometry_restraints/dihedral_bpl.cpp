#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/return_by_value.hpp>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <cctbx/geometry_restraints/dihedral.h>
#include <cctbx/geometry_restraints/proxy_select.h>

namespace cctbx { namespace geometry_restraints {
namespace {

  struct dihedral_proxy_wrappers
  {
    typedef dihedral_proxy w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      typedef default_call_policies dcp;
      object none;
      class_<w_t>("dihedral_proxy", no_init)
        .def(init<
          af::tiny<unsigned, 4> const&, double, double,
          int, alt_angle_ideals_type, double, bool, double>((
            arg("i_seqs"), arg("angle_ideal"), arg("weight"),
            arg("periodicity")=0, arg("alt_angle_ideals")=none,
            arg("limit")=-1.0, arg("top_out")=false,
            arg("slack")=0.0)))
        .def(init<
          af::tiny<unsigned, 4> const&,
          optional_container<af::shared<sgtbx::rt_mx> > const&,
          double,
          double,
          int,
          alt_angle_ideals_type,double,bool,double>((
            arg("i_seqs"), arg("sym_ops"), arg("angle_ideal"), arg("weight"),
            arg("periodicity")=0,
            arg("alt_angle_ideals")=none,
            arg("limit")=-1.0, arg("top_out")=false, arg("slack")=-1.0)))
        .def(init<af::tiny<unsigned, 4> const&, w_t const&>((
          arg("i_seqs"), arg("proxy"))))
        .def("scale_weight", &w_t::scale_weight, (arg("factor")))
        .def("sort_i_seqs", &w_t::sort_i_seqs)
        .add_property("i_seqs", make_getter(&w_t::i_seqs, rbv()))
        .def_readwrite("angle_ideal", &w_t::angle_ideal)
        .def_readwrite("weight", &w_t::weight)
        .def_readwrite("periodicity", &w_t::periodicity)
        .def_readwrite("limit", &w_t::limit)
        .def_readwrite("top_out", &w_t::top_out)
        .def_readwrite("slack", &w_t::slack)
        .add_property("alt_angle_ideals",
          make_getter(&w_t::alt_angle_ideals, rbv()),
          make_setter(&w_t::alt_angle_ideals, dcp()))
        .add_property("sym_ops", make_getter(&w_t::sym_ops, rbv()))
      ;
      {
        typedef return_internal_reference<> rir;
        scitbx::af::boost_python::shared_wrapper<w_t, rir>::wrap(
          "shared_dihedral_proxy")
          .def("count_harmonic", dihedral_count_harmonic)
          .def("proxy_select",
            (af::shared<w_t>(*)(
              af::const_ref<w_t> const&,
              std::size_t,
              af::const_ref<std::size_t> const&))
              shared_proxy_select, (
            arg("n_seq"), arg("iselection")))
          .def("proxy_remove",
            (af::shared<w_t>(*)(
              af::const_ref<w_t> const&,
              af::const_ref<bool> const&))
                shared_proxy_remove, (
            arg("selection")))
        ;
      }
    }
  };

  struct dihedral_wrappers
  {
    typedef dihedral w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      typedef default_call_policies dcp;
      object none;
      class_<w_t>("dihedral", no_init)
        .def(init<
          af::tiny<scitbx::vec3<double>, 4> const&, double, double,
          int, alt_angle_ideals_type const&, double, bool, double>((
            arg("sites"), arg("angle_ideal"), arg("weight"),
            arg("periodicity")=0, arg("alt_angle_ideals")=none,
            arg("limit")=-1.0, arg("top_out")=false, arg("slack")=0.0)))
        .def(init<
          af::const_ref<scitbx::vec3<double> > const&,
          dihedral_proxy const&>((
            arg("sites_cart"), arg("proxy"))))
        .def(init<uctbx::unit_cell const&,
          af::const_ref<scitbx::vec3<double> > const&,
          dihedral_proxy const&>((
            arg("unit_cell"), arg("sites_cart"), arg("proxy"))))
        .add_property("sites", make_getter(&w_t::sites, rbv()))
        .def_readonly("angle_ideal", &w_t::angle_ideal)
        .def_readonly("weight", &w_t::weight)
        .def_readonly("periodicity", &w_t::periodicity)
        .def_readonly("limit", &w_t::limit)
        .def_readonly("top_out", &w_t::top_out)
        .def_readonly("slack", &w_t::slack)
        .add_property("alt_angle_ideals",
          make_getter(&w_t::alt_angle_ideals, rbv()),
          make_setter(&w_t::alt_angle_ideals, dcp()))
        .def_readonly("have_angle_model", &w_t::have_angle_model)
        .def_readonly("angle_model", &w_t::angle_model)
        .def_readonly("delta", &w_t::delta)
        .def("residual", &w_t::residual)
        .def("gradients", &w_t::gradients, (arg("epsilon")=1e-100))
      ;
    }
  };

  void
  wrap_all()
  {
    using namespace boost::python;
    dihedral_proxy_wrappers::wrap();
    dihedral_wrappers::wrap();
    def("dihedral_deltas",
      (af::shared<double>(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<dihedral_proxy> const&))
        dihedral_deltas,
      (arg("sites_cart"), arg("proxies")));
    def("dihedral_residuals",
      (af::shared<double>(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<dihedral_proxy> const&))
        dihedral_residuals,
      (arg("sites_cart"), arg("proxies")));
    def("dihedral_residual_sum",
      (double(*)(
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<dihedral_proxy> const&,
        af::ref<scitbx::vec3<double> > const&))
        dihedral_residual_sum,
      (arg("sites_cart"), arg("proxies"), arg("gradient_array")));
    def("dihedral_deltas",
      (af::shared<double>(*)(
        uctbx::unit_cell const&,
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<dihedral_proxy> const&))
        dihedral_deltas,
      (arg("unit_cell"), arg("sites_cart"), arg("proxies")));
    def("dihedral_residuals",
      (af::shared<double>(*)(
        uctbx::unit_cell const&,
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<dihedral_proxy> const&))
        dihedral_residuals,
      (arg("unit_cell"), arg("sites_cart"), arg("proxies")));
    def("dihedral_residual_sum",
      (double(*)(
        uctbx::unit_cell const&,
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<dihedral_proxy> const&,
        af::ref<scitbx::vec3<double> > const&))
        dihedral_residual_sum,
      (arg("unit_cell"), arg("sites_cart"), arg("proxies"), arg("gradient_array")));
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_dihedral() { wrap_all(); }

}}} // namespace cctbx::geometry_restraints::boost_python
