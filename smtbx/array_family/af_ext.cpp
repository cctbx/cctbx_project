#include <boost_adaptbx/tuple_conversion.h>
#include <cctbx/coordinates.h>

namespace smtbx { namespace af { namespace boost_python {
namespace {

void init_module() {
  {
    using boost_adaptbx::tuple_conversion::to_python;
    typedef cctbx::cartesian<double> cart_t;
    // used in smtbx/refinement/constraints/geometric_hydrogen.h
    to_python<boost::tuple<cart_t, cart_t, cart_t> >();
  }
}

}}}}

BOOST_PYTHON_MODULE(smtbx_array_family_ext)
{
  smtbx::af::boost_python::init_module();
}
