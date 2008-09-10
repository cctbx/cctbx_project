#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

#include <cmaplib.h>

namespace iotbx { namespace ccp4_map {

  void
  test_read(
    const char* file_name)
  {
    CMap_io::CMMFile* mfile = static_cast<CMap_io::CMMFile*>(
      CMap_io::ccp4_cmap_open(file_name, O_RDONLY));
    CMap_io::ccp4_cmap_close(mfile);
  }

  void
  init_module()
  {
    using namespace boost::python;
    def("test_read", test_read, (arg_("file_name")));
  }

}} // namespace iotbx::ccp4_map

BOOST_PYTHON_MODULE(iotbx_ccp4_map_ext)
{
  iotbx::ccp4_map::init_module();
}
