#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <scitbx/boost_python/iterator_wrappers.h>
#include <cctbx/eltbx/neutron.h>

namespace cctbx { namespace eltbx { namespace neutron {
namespace boost_python {

namespace {

  struct neutron_news_1992_table_wrappers
  {
    typedef neutron_news_1992_table w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("neutron_news_1992_table", no_init)
        .def(init<std::string const&, optional<bool> >())
        .def("label", &w_t::label)
        .def("bound_coh_scatt_length", &w_t::bound_coh_scatt_length)
        .def("abs_cross_sect", &w_t::abs_cross_sect)
      ;
    }
  };

  void init_module()
  {
    neutron_news_1992_table_wrappers::wrap();
    scitbx::boost_python::iterator_wrappers<
      neutron_news_1992_table,
      neutron_news_1992_table_iterator>::wrap(
        "neutron_news_1992_table_iterator");
  }

} // namespace <anonymous>
}}}} // namespace cctbx::eltbx::neutron::boost_python

BOOST_PYTHON_MODULE(cctbx_eltbx_neutron_ext)
{
  cctbx::eltbx::neutron::boost_python::init_module();
}
