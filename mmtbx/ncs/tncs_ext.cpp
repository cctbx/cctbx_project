#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <mmtbx/ncs/tncs.h>
#include <boost/python/return_value_policy.hpp>
#include <boost/python.hpp>

namespace mmtbx { namespace ncs {
  namespace bp = boost::python;

namespace {
  boost::python::tuple
  getinitargs(pair<> const& self)
  {
    return boost::python::make_tuple(self.r, self.t);
  }

  using namespace boost::python;
  using boost::python::arg;

  template <typename FloatType>
  struct tncs_eps_factor_refinery_wrapper
  {
    typedef tncs_eps_factor_refinery<FloatType> w_t;
    static void wrap() {
      using namespace boost::python;
      class_<w_t>("tncs_eps_factor_refinery", no_init)
        .def(init<
             bp::list const&,
             scitbx::mat3<double>,
             af::shared<mat3<double> >
                  >((arg("xxx"),
                     arg("ccc"),
                     arg("ddd"))))
        .def("tncs_epsfac", &w_t::tncs_epsfac_result)
      ;
    }
  };

  void init_module()
  {
    tncs_eps_factor_refinery_wrapper<double>::wrap();

    typedef return_value_policy<return_by_value> rbv;
    class_<pair<> >("pair")
      .def(init<scitbx::mat3<double> const&,
                scitbx::vec3<double> const&,
                double,
                double,
                double,
                af::shared<double> >((
                  arg("r"),
                  arg("t"),
                  arg("radius"),
                  arg("weight"),
                  arg("fracscat"),
                  arg("rho_mn"))))
      .add_property("r",      make_getter(&pair<>::r,      rbv()))
      .add_property("t",      make_getter(&pair<>::t,      rbv()))
      .add_property("radius", make_getter(&pair<>::radius, rbv()))
      .add_property("weight", make_getter(&pair<>::weight, rbv()))
      .add_property("fracscat", make_getter(&pair<>::fracscat, rbv()))
      .add_property("rho_mn", make_getter(&pair<>::rho_mn, rbv()))
      .enable_pickling()
      .def("__getinitargs__", getinitargs)
    ;

  }

} // namespace <anonymous>
}} // namespace mmtbx::ncs

BOOST_PYTHON_MODULE(mmtbx_ncs_ext)
{
  mmtbx::ncs::init_module();
}
