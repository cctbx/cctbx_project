#include <cctbx/boost_python/flex_fwd.h>
#include <iotbx/mtzwriter.h>
#include <boost/python/class.hpp>

namespace iotbx { namespace mtz { namespace boost_python {

  void wrap_MtzWriter()
  {
    using namespace boost::python;
    class_<MtzWriter>("MtzWriter", init<>())
      .def("setTitle",      &MtzWriter::setTitle)
      .def("raw_setSpaceGroup", &MtzWriter::setSpaceGroup)
      .def("oneCrystal",    &MtzWriter::oneCrystal)
      .def("oneDataset",    &MtzWriter::oneDataset)
      .def("addColumn",     &MtzWriter::addColumn)
      .def("write",         &MtzWriter::write)
    ;
  }

}}} // namespace iotbx::mtz::boost_python
