#include <boost/python.hpp>
#include <cctbx/boost_python/flex_fwd.h>
#include <scitbx/boost_python/utils.h>
#include "iotbx/cppmtz.h"
using namespace boost::python;

BOOST_PYTHON_MODULE(mtz)
{
    class_<iotbx::mtz::Mtz>("Mtz", init<std::string>())
      .def("title",      &iotbx::mtz::Mtz::title)
      .def("SpaceGroup", &iotbx::mtz::Mtz::SpaceGroup)
      .def("size",       &iotbx::mtz::Mtz::size,
          return_value_policy<copy_non_const_reference>())
      .def("ncrystals",  &iotbx::mtz::Mtz::ncrystals,
          return_value_policy<copy_non_const_reference>())
      .def("columns",    &iotbx::mtz::Mtz::columns)
      .def("history",    &iotbx::mtz::Mtz::history)
      .def("getColumn",  &iotbx::mtz::Mtz::getColumn)
      .def("getShared",  &iotbx::mtz::Mtz::getShared)
      .def("getCrystal", &iotbx::mtz::Mtz::getCrystal)
      .def("printHeader",&iotbx::mtz::Mtz::printHeader)
      .def("UnitCell",   &iotbx::mtz::Mtz::UnitCell)
      .def("ndatasets",  &iotbx::mtz::Mtz::ndatasets)
      .def("ncolumns",   &iotbx::mtz::Mtz::ncolumns,
           return_value_policy<copy_non_const_reference>())
      .def("MIx",        &iotbx::mtz::Mtz::MIx)
      .def("printHeaderAdv",&iotbx::mtz::Mtz::printHeaderAdv)
    ;

    class_<iotbx::mtz::Crystal>("Crystal", no_init)
      .def("crystal_name",&iotbx::mtz::Crystal::crystal_name)
      .def("project_name",&iotbx::mtz::Crystal::project_name)
      .def("ndatasets",   &iotbx::mtz::Crystal::ndatasets)
      .def("UnitCell",    &iotbx::mtz::Crystal::UnitCell)
      .def("getDataset",  &iotbx::mtz::Crystal::getDataset)
    ;

    class_<iotbx::mtz::Dataset>("Dataset", no_init)
      .def("dataset_name",&iotbx::mtz::Dataset::dataset_name)
      .def("wavelength",  &iotbx::mtz::Dataset::wavelength)
      .def("ncolumns",    &iotbx::mtz::Dataset::ncolumns)
      .def("getColumn",   &iotbx::mtz::Dataset::getColumn)
    ;

    class_<iotbx::mtz::Column>("Column", no_init)
      .def("label",       &iotbx::mtz::Column::label)
      .def("type",        &iotbx::mtz::Column::type)
      .def("__getitem__", &iotbx::mtz::Column::lookup)
      .def("__call__",    &iotbx::mtz::Column::lookup)
    ;

    class_<iotbx::mtz::Foo>("Foo")
      .def("value",&iotbx::mtz::Foo::value)
    ;
}
