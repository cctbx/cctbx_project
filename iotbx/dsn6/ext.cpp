
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <iotbx/dsn6/dsn6.hpp>

namespace iotbx { namespace dsn6 { namespace boost_python {

  void
  init_module()
  {
    using namespace boost::python;
    def("write_dsn6_map", write_dsn6_map, (
        arg("file_name"),
        arg("unit_cell"),
        arg("gridding_first"),
        arg("gridding_last"),
        arg("map_data")));
  }

}}} // namespace iotbx::ccp4_map

BOOST_PYTHON_MODULE(iotbx_dsn6_map_ext)
{
  iotbx::dsn6::boost_python::init_module();
}
