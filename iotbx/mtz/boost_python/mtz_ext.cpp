#include <cctbx/boost_python/flex_fwd.h>
#include <iotbx/cppmtz.h>
#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_non_const_reference.hpp>

namespace iotbx { namespace mtz { namespace boost_python {

  void wrap_Mtz()
  {
    using namespace boost::python;
    class_<Mtz>("Mtz", init<std::string>())
      .def("title",      &Mtz::title)
      .def("SpaceGroup", &Mtz::SpaceGroup)
      .def("getSgtbxSpaceGroup", &Mtz::getSgtbxSpaceGroup)
      .def("nsym",       &Mtz::nsym)
      .def("size",       &Mtz::size,
          return_value_policy<copy_non_const_reference>())
      .def("ncrystals",  &Mtz::ncrystals,
          return_value_policy<copy_non_const_reference>())
      .def("columns",    &Mtz::columns)
      .def("history",    &Mtz::history)
      .def("getColumn",  &Mtz::getColumn)
      .def("getCrystal", &Mtz::getCrystal)
      .def("columnToCrystal", &Mtz::columnToCrystal)
      .def("printHeader",&Mtz::printHeader)
      .def("UnitCell",   &Mtz::UnitCell)
      .def("ndatasets",  &Mtz::ndatasets)
      .def("ncolumns",   &Mtz::ncolumns,
           return_value_policy<copy_non_const_reference>())
      .def("MIx",        &Mtz::MIx)
      .def("printHeaderAdv",&Mtz::printHeaderAdv)
      .def("valid_indices", &Mtz::valid_indices)
      .def("valid_indices", &Mtz::valid_indices_anomalous)
      .def("valid_values", &Mtz::valid_values)
      .def("valid_values", &Mtz::valid_values_anomalous)
      .def("valid_complex", &Mtz::valid_complex)
      .def("valid_complex", &Mtz::valid_complex_anomalous)
      .def("valid_hl", &Mtz::valid_hl)
      .def("valid_integers", &Mtz::valid_integers)
    ;
  }

  void wrap_Crystal()
  {
    using namespace boost::python;
    class_<Crystal>("Crystal", no_init)
      .def("crystal_name",&Crystal::crystal_name)
      .def("project_name",&Crystal::project_name)
      .def("ndatasets",   &Crystal::ndatasets)
      .def("UnitCell",    &Crystal::UnitCell)
      .def("getDataset",  &Crystal::getDataset)
    ;
  }

  void wrap_Dataset()
  {
    using namespace boost::python;
    class_<Dataset>("Dataset", no_init)
      .def("dataset_name",&Dataset::dataset_name)
      .def("wavelength",  &Dataset::wavelength)
      .def("ncolumns",    &Dataset::ncolumns)
      .def("getColumn",   &Dataset::getColumn)
    ;
  }

  void wrap_Column()
  {
    using namespace boost::python;
    class_<Column>("Column", no_init)
      .def("label",       &Column::label)
      .def("type",        &Column::type)
      .def("__getitem__", &Column::lookup) // XXX potential crash!
      .def("__call__",    &Column::lookup) // XXX potential crash!
      .def("isnan",       &Column::ccp4_isnan) // XXX potential crash!
    ;
  }

  void wrap_MtzWriter();

}}} // namespace iotbx::mtz::boost_python

BOOST_PYTHON_MODULE(iotbx_mtz_ext)
{
  using namespace iotbx::mtz::boost_python;
  wrap_Mtz();
  wrap_Crystal();
  wrap_Dataset();
  wrap_Column();
  wrap_MtzWriter();
}
