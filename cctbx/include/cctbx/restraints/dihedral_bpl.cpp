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
#include <cctbx/restraints/dihedral.h>

namespace cctbx { namespace restraints {
namespace {

  struct dihedral_proxy_wrappers
  {
    typedef dihedral_proxy w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("dihedral_proxy", no_init)
        .def(init<af::tiny<unsigned, 4> const&, double, double,
                  optional<int> >(
          (arg_("i_seqs"), arg_("angle_ideal"), arg_("weight"),
           arg_("periodicity")=0)))
        .add_property("i_seqs", make_getter(&w_t::i_seqs, rbv()))
        .def_readonly("angle_ideal", &w_t::angle_ideal)
        .def_readonly("weight", &w_t::weight)
        .def_readonly("periodicity", &w_t::periodicity)
      ;
      {
        scitbx::af::boost_python::shared_wrapper<w_t>::wrap(
          "shared_dihedral_proxy");
      }
    }
  };

  struct dihedral_wrappers
  {
    typedef dihedral w_t;

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      gradients_overloads, gradients, 0, 1)

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("dihedral", no_init)
        .def(init<af::tiny<scitbx::vec3<double>, 4> const&, double, double,
                  optional<int> >(
          (arg_("sites"), arg_("angle_ideal"), arg_("weight"),
           arg_("periodicity")=0)))
        .def(init<af::const_ref<scitbx::vec3<double> > const&,
                  dihedral_proxy const&>(
          (arg_("sites_cart"), arg_("proxy"))))
        .add_property("sites", make_getter(&w_t::sites, rbv()))
        .def_readonly("angle_ideal", &w_t::angle_ideal)
        .def_readonly("weight", &w_t::weight)
        .def_readonly("periodicity", &w_t::periodicity)
        .def_readonly("have_angle_model", &w_t::have_angle_model)
        .def_readonly("angle_model", &w_t::angle_model)
        .def_readonly("delta", &w_t::delta)
        .def("residual", &w_t::residual)
        .def("gradients", &w_t::gradients, gradients_overloads(
          (arg_("epsilon"))))
      ;
    }
  };

  void
  wrap_all()
  {
    using namespace boost::python;
    typedef boost::python::arg arg_; // gcc 2.96 workaround
    dihedral_proxy_wrappers::wrap();
    dihedral_wrappers::wrap();
    def("dihedral_deltas", dihedral_deltas,
      (arg_("sites_cart"), arg_("proxies")));
    def("dihedral_residuals", dihedral_residuals,
      (arg_("sites_cart"), arg_("proxies")));
    def("dihedral_residual_sum", dihedral_residual_sum,
      (arg_("sites_cart"), arg_("proxies"), arg_("gradient_array")));
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_dihedral() { wrap_all(); }

}}} // namespace cctbx::restraints::boost_python
