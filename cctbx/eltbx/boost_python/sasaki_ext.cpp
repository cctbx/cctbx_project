#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <scitbx/boost_python/iterator_wrappers.h>
#include <cctbx/eltbx/sasaki.h>

namespace cctbx { namespace eltbx { namespace sasaki { namespace boost_python {

namespace {

  struct table_wrappers
  {
    typedef table w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("table", no_init)
        .def(init<std::string const&, optional<bool> >())
        .def("label", &w_t::label)
        .def("atomic_number", &w_t::atomic_number)
        .def("at_ev", &w_t::at_ev)
        .def("at_kev", &w_t::at_kev)
        .def("at_angstrom", &w_t::at_angstrom)
      ;
    }
  };

  void init_module()
  {
    table_wrappers::wrap();
    scitbx::boost_python::iterator_wrappers<
      table, table_iterator>::wrap("table_iterator");
  }

} // namespace <anonymous>
}}}} // namespace cctbx::eltbx::sasaki::boost_python

BOOST_PYTHON_MODULE(cctbx_eltbx_sasaki_ext)
{
  cctbx::eltbx::sasaki::boost_python::init_module();
}
