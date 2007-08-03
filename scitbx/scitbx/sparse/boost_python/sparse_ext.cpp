#include <boost/python/module.hpp>

namespace scitbx { namespace sparse { namespace boost_python {

void wrap_vector();
void wrap_matrix();

void init_module() {
  wrap_vector();
  wrap_matrix();
}

}}}


BOOST_PYTHON_MODULE(scitbx_sparse_ext)
{
  scitbx::sparse::boost_python::init_module();
}
