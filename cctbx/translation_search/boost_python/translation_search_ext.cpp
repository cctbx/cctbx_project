#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>

namespace cctbx { namespace translation_search { namespace boost_python {

  void wrap_fast_nv1995();
  void wrap_fast_terms();
  void wrap_symmetry_flags();

namespace {

  void init_module()
  {
    wrap_fast_nv1995();
    wrap_fast_terms();
    wrap_symmetry_flags();
  }

} // namespace <anonymous>
}}} // namespace cctbx::translation_search::boost_python

BOOST_PYTHON_MODULE(translation_search_ext)
{
  cctbx::translation_search::boost_python::init_module();
}
