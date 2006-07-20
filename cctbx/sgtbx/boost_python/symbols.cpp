#include <scitbx/boost_python/iterator_wrappers.h>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <cctbx/sgtbx/symbols.h>

namespace cctbx { namespace sgtbx { namespace boost_python {

namespace {

  struct space_group_symbols_wrappers
  {
    typedef space_group_symbols w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("space_group_symbols", no_init)
        .def(init<std::string const&, optional<std::string const&> >((
          arg_("symbol"),
          arg_("table_id")="")))
        .def(init<int, optional<std::string const&, std::string const&> >((
          arg_("space_group_number"),
          arg_("extension")="",
          arg_("table_id")="")))
        .def("number", &w_t::number)
        .def("schoenflies", &w_t::schoenflies, ccr())
        .def("qualifier", &w_t::qualifier, ccr())
        .def("hermann_mauguin", &w_t::hermann_mauguin, ccr())
        .def("extension", &w_t::extension)
        .def("change_of_basis_symbol", &w_t::change_of_basis_symbol, ccr())
        .def("universal_hermann_mauguin",
          &w_t::universal_hermann_mauguin, ccr())
        .def("hall", &w_t::hall, ccr())
        .def("point_group_type", &w_t::point_group_type)
        .def("laue_group_type", &w_t::laue_group_type)
        .def("crystal_system", &w_t::crystal_system)
      ;
    }
  };

} // namespace <anoymous>

  void wrap_symbols()
  {
    space_group_symbols_wrappers::wrap();
    scitbx::boost_python::iterator_wrappers<
      space_group_symbols,
      space_group_symbol_iterator>::wrap(
        "space_group_symbol_iterator");
  }

}}} // namespace cctbx::sgtbx::boost_python
