#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <iotbx/mtz/column.h>
#include <scitbx/array_family/boost_python/shared_wrapper.h>

namespace iotbx { namespace mtz {
namespace {

  struct dataset_wrappers
  {
    typedef dataset w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("dataset", no_init)
        .def(init<crystal const&, int>((
          arg_("mtz_crystal"), arg_("i_dataset"))))
        .def("mtz_crystal", &w_t::mtz_crystal)
        .def("i_dataset", &w_t::i_dataset)
        .def("mtz_object", &w_t::mtz_object)
        .def("id", &w_t::id)
        .def("name", &w_t::name)
        .def("wavelength", &w_t::wavelength)
        .def("n_columns", &w_t::n_columns)
        .def("columns", &w_t::columns)
        .def("add_column", &w_t::add_column, (
          arg_("label"), arg_("type")))
      ;
      {
        scitbx::af::boost_python::shared_wrapper<w_t>::wrap(
          "shared_dataset");
      }
    }
  };

  void
  wrap_all()
  {
    dataset_wrappers::wrap();
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_dataset() { wrap_all(); }

}}} // namespace iotbx::mtz::boost_python
