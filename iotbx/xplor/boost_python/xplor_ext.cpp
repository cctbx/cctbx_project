#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/def.hpp>

#include "iotbx/xplor/map_reader.h"
#include "iotbx/xplor/map_writer.h"

#include <cctbx/uctbx.h>

namespace af = scitbx::af;

namespace iotbx { namespace boost_python { namespace xplor_ext {

  void init_module()
  {
    using namespace boost::python;
    using namespace iotbx::xplor;

    class_<map_reader>("map_reader", no_init)
      .def(init<std::string const&, std::size_t, af::flex_grid<> const&>(
        (arg("file_name"), arg("n_header_lines"), arg("grid"))))
      .def_readonly("data", &map_reader::data)
      .def_readonly("average", &map_reader::average)
      .def_readonly("standard_deviation", &map_reader::standard_deviation)
    ;

    def("map_writer", map_writer_box,
      (arg("file_name"), arg("unit_cell"),
       arg("data"),
       arg("average"), arg("standard_deviation")));
    def("map_writer", map_writer_p1_cell,
      (arg("file_name"), arg("unit_cell"),
       arg("gridding_first"), arg("gridding_last"), arg("data"),
       arg("average"), arg("standard_deviation")));
  }

}}} // namespace iotbx::boost_python::xplor_ext

BOOST_PYTHON_MODULE(iotbx_xplor_ext)
{
  iotbx::boost_python::xplor_ext::init_module();
}
