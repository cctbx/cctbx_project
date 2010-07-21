#include <boost/python/module.hpp>

namespace scitbx { namespace lapack { namespace boost_python {
  void wrap();
}}}

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
      lapack::boost_python::wrap();
    }
  }
}}}

BOOST_PYTHON_MODULE(scitbx_linalg_ext)
{
  scitbx::matrix::boost_python::init_module();
}
