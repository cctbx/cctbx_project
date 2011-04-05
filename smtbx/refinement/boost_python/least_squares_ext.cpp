#include <boost/python/module.hpp>

namespace smtbx { namespace refinement { namespace least_squares {
  namespace boost_python {

  void wrap_weighting_schemes();
  void wrap_least_squares();
  namespace {
    void init_module() {
      wrap_weighting_schemes();
      wrap_least_squares();
    }
  }
}}}}

BOOST_PYTHON_MODULE(smtbx_refinement_least_squares_ext)
{
  smtbx::refinement::least_squares::boost_python::init_module();
}
