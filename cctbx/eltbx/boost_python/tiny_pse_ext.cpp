#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <scitbx/boost_python/iterator_wrappers.h>
#include <cctbx/eltbx/tiny_pse.h>

namespace cctbx { namespace eltbx { namespace tiny_pse {
namespace boost_python {

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
        .def(init<int>())
        .def("atomic_number", &w_t::atomic_number)
        .def("symbol", &w_t::symbol)
        .def("name", &w_t::name)
        .def("weight", &w_t::weight)
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
}}}} // namespace cctbx::eltbx::tiny_pse::boost_python

BOOST_PYTHON_MODULE(cctbx_eltbx_tiny_pse_ext)
{
  cctbx::eltbx::tiny_pse::boost_python::init_module();
}
