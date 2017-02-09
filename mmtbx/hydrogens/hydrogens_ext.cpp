#include <boost/python/class.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python.hpp>
#include <mmtbx/hydrogens/hydrogens.h>

//SCITBX_BOOST_IS_POLYMORPHIC_WORKAROUND(mmtbx::hydrogens::common)

namespace mmtbx { namespace hydrogens {
  namespace bp = boost::python;

namespace {
  boost::python::tuple
  getinitargs(riding_coefficients const& self)
  {
    return boost::python::make_tuple(self.a0);
  }

  void init_module()
  {
    using namespace boost::python;
    using boost::python::arg;
    typedef return_value_policy<return_by_value> rbv;

    class_<riding_coefficients>("riding_coefficients")
      .def(init<string const&, int const&, int const&, int const&, int const&,
                double const&, double const&,
                double const&, int const&, double const& >(
                    (arg("htype"),arg("a0"),arg("a1"),arg("a2"), arg("a3"),
                     arg("a"), arg("b"), arg("h"), arg("n"), arg("disth"))))
      .add_property("htype", make_getter(&riding_coefficients::htype, rbv()))
      .add_property("a0", make_getter(&riding_coefficients::a0, rbv()))
      .add_property("a1", make_getter(&riding_coefficients::a1, rbv()))
      .add_property("a2", make_getter(&riding_coefficients::a2, rbv()))
      .add_property("a3", make_getter(&riding_coefficients::a3, rbv()))
      .add_property("a", make_getter(&riding_coefficients::a, rbv()))
      .add_property("b", make_getter(&riding_coefficients::b, rbv()))
      .add_property("h", make_getter(&riding_coefficients::h, rbv()))
      .add_property("n", make_getter(&riding_coefficients::n, rbv()))
      .add_property("disth", make_getter(&riding_coefficients::disth, rbv()))
      .def("__getinitargs__", getinitargs)
    ;

//    def("compute_H_position",
//         (vec3<double>(*)
//               (riding_coefficients,
//                af::shared<vec3<double> > const&,
//                int const&)) compute_H_position,
//                  (arg("riding_coefficients"),
//                   arg("sites_cart"),
//                   arg("ih")))
//    ;

  }

} // namespace <anonymous>
}} // namespace mmtbx::hydrogens

BOOST_PYTHON_MODULE(mmtbx_hydrogens_ext)
{
  mmtbx::hydrogens::init_module();
}
