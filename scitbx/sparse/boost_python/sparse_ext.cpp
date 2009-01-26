#include <boost/python/module.hpp>

namespace scitbx { namespace sparse { namespace boost_python {

void wrap_vector();
void wrap_matrix();
void wrap_lu_factorization();
void wrap_approx_equal();
void wrap_random();

void init_module() {
  wrap_vector();
  wrap_matrix();
  wrap_random();
  wrap_approx_equal();
  wrap_lu_factorization();
}

}}}


BOOST_PYTHON_MODULE(scitbx_sparse_ext)
{
  scitbx::sparse::boost_python::init_module();
}
