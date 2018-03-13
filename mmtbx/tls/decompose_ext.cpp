#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <mmtbx/tls/tls.h>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <scitbx/boost_python/is_polymorphic_workaround.h>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python.hpp>

SCITBX_BOOST_IS_POLYMORPHIC_WORKAROUND(mmtbx::tls::common)

namespace mmtbx { namespace tls { namespace decompose {
  namespace bp = boost::python;

namespace {

  void init_module()
  {
    using namespace boost::python;
    using boost::python::arg;

   class_<decompose_tls_matrices>("decompose_tls_matrices",
           init< sym_mat3<double> const&,
                 sym_mat3<double> const&,
                 mat3<double> const&,
                 vec3<double> const& >())
       .def("is_valid", &decompose_tls_matrices::is_valid)
   ;

  }

} // namespace <anonymous>
}}} // namespace mmtbx::tls::decompose

BOOST_PYTHON_MODULE(mmtbx_tls_decompose_ext)
{
  mmtbx::tls::decompose::init_module();
}
