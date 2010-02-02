#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>

#include <scitbx/wigner/boost_python/wigner3j.h>

namespace scitbx { namespace wigner {
namespace boost_python{

  void wrap_wigner3j_fast();
  void wrap_wigner3j();
  void wrap_wigner3j_zero();
  namespace {

    void init_module()
    {
      using namespace boost::python;
      wrap_wigner3j_fast();
      wrap_wigner3j();
      wrap_wigner3j_zero();
    }

}

}}} // scitbx::wigner::<anonymous>

BOOST_PYTHON_MODULE(scitbx_wigner_ext)
{
  scitbx::wigner::boost_python::init_module();
}
