#include <boost/python/module.hpp>
#include <scitbx/stl/map_fwd.h>
#include <scitbx/stl/map_wrapper.h>
#include <cctbx/xray/scatterer_flags.h>
#include <map>

namespace smtbx { namespace stl { namespace boost_python {
namespace {

void init_module() {
  typedef boost::python::return_internal_reference<> rir;
  // used by smtbx.refinement.constraints
  scitbx::stl::boost_python::map_wrapper<
    std::map<std::size_t,
             cctbx::xray::scatterer_flags> >::wrap(
               "std_size_t_xray_scatterer_flags");
}

}}}}

BOOST_PYTHON_MODULE(smtbx_stl_map_ext)
{
  smtbx::stl::boost_python::init_module();
}
