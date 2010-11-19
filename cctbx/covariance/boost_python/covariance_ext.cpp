#include <boost/python/module.hpp>
#include <boost/python/def.hpp>

#include <cctbx/covariance/covariance.h>

namespace cctbx { namespace covariance {
  namespace boost_python {

    void wrap_covariance_matrix() {
      using namespace boost::python;

      def("extract_covariance_matrix_for_sites", (
        af::versa<double, af::packed_u_accessor>(*)(
        af::const_ref<std::size_t> const &,
        af::const_ref<double, af::packed_u_accessor> const &,
        cctbx::xray::parameter_map<cctbx::xray::scatterer<double> > const &))
        extract_covariance_matrix_for_sites,
        (arg("i_seqs"), arg("matrix"), arg("parameter_map")));
      def("extract_covariance_matrix_for_u_aniso", (
        af::versa<double, af::packed_u_accessor>(*)(
        std::size_t,
        af::const_ref<double, af::packed_u_accessor> const &,
        cctbx::xray::parameter_map<cctbx::xray::scatterer<double> > const &))
        extract_covariance_matrix_for_u_aniso,
        (arg("i_seq"), arg("matrix"), arg("parameter_map")));
      def("variance_for_u_iso", (
        double(*)(
        std::size_t,
        af::const_ref<double, af::packed_u_accessor> const &,
        cctbx::xray::parameter_map<cctbx::xray::scatterer<double> > const &))
        variance_for_u_iso,
        (arg("i_seq"), arg("matrix"), arg("parameter_map")));
      def("orthogonalize_covariance_matrix", (
        af::versa<double, af::packed_u_accessor>(*)(
        af::const_ref<double, af::packed_u_accessor> const &,
        cctbx::uctbx::unit_cell const &,
        cctbx::xray::parameter_map<cctbx::xray::scatterer<double> > const &))
        orthogonalize_covariance_matrix,
        (arg("matrix"), arg("unit_cell"), arg("parameter_map")));
      def("covariance_orthogonalization_matrix", (
       scitbx::sparse::matrix<double>(*)(
        cctbx::uctbx::unit_cell const &,
        cctbx::xray::parameter_map<cctbx::xray::scatterer<double> > const &))
        covariance_orthogonalization_matrix,
        (arg("unit_cell"), arg("parameter_map")));
    }

    void init_module()
    {
      wrap_covariance_matrix();
    }



}}} // namespace cctbx::covariance::boost_python

BOOST_PYTHON_MODULE(cctbx_covariance_ext)
{
  cctbx::covariance::boost_python::init_module();
}
