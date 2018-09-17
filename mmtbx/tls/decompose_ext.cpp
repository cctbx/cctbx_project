#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <scitbx/boost_python/is_polymorphic_workaround.h>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python.hpp>

#include <scitbx/array_family/shared.h>
#include <mmtbx/tls/decompose.h>

SCITBX_BOOST_IS_POLYMORPHIC_WORKAROUND(mmtbx::tls::common)

namespace mmtbx { namespace tls { namespace decompose {
  namespace bp = boost::python;
  namespace af = scitbx::af;

namespace {

  boost::python::tuple get_v_amplitudes(const decompose_tls_matrices& a)
  { return boost::python::tuple(a.v_amplitudes); }

  boost::python::tuple get_l_amplitudes(const decompose_tls_matrices& a)
  { return boost::python::tuple(a.l_amplitudes); }

  boost::python::tuple get_s_amplitudes(const decompose_tls_matrices& a)
  { return boost::python::tuple(a.s_amplitudes); }

  af::shared<scitbx::vec3<double>> get_v_axis_directions(const decompose_tls_matrices& a)
  {
      af::shared<scitbx::vec3<double>> s(3);
      s[0] = a.v1_M;
      s[1] = a.v2_M;
      s[2] = a.v3_M;
      return s;
  }

  af::shared<scitbx::vec3<double>> get_l_axis_directions(const decompose_tls_matrices& a)
  {
      af::shared<scitbx::vec3<double>> s(3);
      s[0] = a.l1_M;
      s[1] = a.l2_M;
      s[2] = a.l3_M;
      return s;
  }

  af::shared<scitbx::vec3<double>> get_l_axis_intersections(const decompose_tls_matrices& a)
  {
      af::shared<scitbx::vec3<double>> s(3);
      s[0] = a.w1_M;
      s[1] = a.w2_M;
      s[2] = a.w3_M;
      return s;
  }

  boost::python::tuple get_R_ML(const decompose_tls_matrices& a)
  { return boost::python::tuple(a.R_ML); }

  boost::python::tuple get_R_MV(const decompose_tls_matrices& a)
  { return boost::python::tuple(a.R_MV); }

  void init_module()
  {
    using namespace boost::python;
    using boost::python::arg;

    class_<decompose_tls_matrices>("decompose_tls_matrices",
            init< scitbx::sym_mat3<double> const&,
                  scitbx::sym_mat3<double> const&,
                  scitbx::mat3<double> const&,
                  bool, bool,
                  double, double,
                  std::string, double>(
                      (arg("T"),
                       arg("L"),
                       arg("S"),
                       arg("l_and_s_in_degrees")=true,
                       arg("verbose")=false,
                       arg("tol")=1.e-6,
                       arg("eps")=1.e-8,
                       arg("t_S_formula")="11",
                       arg("t_S_value")=0.0)))

        .def("is_valid", &decompose_tls_matrices::is_valid)
        .def("error", &decompose_tls_matrices::error)

        .add_property("v_amplitudes", get_v_amplitudes)
        .add_property("l_amplitudes", get_l_amplitudes)
        .add_property("s_amplitudes", get_s_amplitudes)

        .add_property("v_axis_directions",    get_v_axis_directions)
        .add_property("l_axis_directions",    get_l_axis_directions)
        .add_property("l_axis_intersections", get_l_axis_intersections)

        .add_property("R_ML", get_R_ML)
        .add_property("R_MV", get_R_MV)

    ;

  }

} // namespace <anonymous>
}}} // namespace mmtbx::tls::decompose

BOOST_PYTHON_MODULE(mmtbx_tls_decompose_ext)
{
  mmtbx::tls::decompose::init_module();
}
