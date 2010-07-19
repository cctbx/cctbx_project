#include <boost/python/module.hpp>

#if defined(SCITBX_LAPACK_FEM)
namespace lapack_fem { namespace boost_python {
  void wrap();
}}
#endif

namespace scitbx { namespace matrix { namespace boost_python {

  void wrap_matrix();
  void wrap_householder();
  void wrap_svd();
  void wrap_eigensystem();
  void wrap_cholesky();

  namespace {
    void init_module() {
      wrap_matrix();
      wrap_householder();
      wrap_svd();
      wrap_eigensystem();
      wrap_cholesky();
#if defined(SCITBX_LAPACK_FEM)
      lapack_fem::boost_python::wrap();
#endif
    }
  }
}}}

BOOST_PYTHON_MODULE(scitbx_linalg_ext)
{
  scitbx::matrix::boost_python::init_module();
}
