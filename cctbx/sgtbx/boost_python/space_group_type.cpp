#include <boost/python/class.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <cctbx/sgtbx/space_group_type.h>

namespace cctbx { namespace sgtbx { namespace boost_python {

namespace {

  struct space_group_type_wrappers : boost::python::pickle_suite
  {
    typedef space_group_type w_t;

    static boost::python::tuple
    getinitargs(w_t const& self)
    {
      return boost::python::make_tuple(
        "Hall: " + self.hall_symbol(false),
        "",
        self.cb_op_is_tidy());
    }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      typedef return_internal_reference<> rir;
      class_<w_t>("space_group_type")
        .def(init<std::string const&, optional<std::string const&, bool> >((
          arg("symbol"),
          arg("table_id")="",
          arg("tidy_cb_op")=true)))
        .def(init<space_group const&, optional<bool, int, int> >((
          arg("group"),
          arg("tidy_cb_op")=true,
          arg("r_den")=cb_r_den,
          arg("t_den")=cb_t_den)))
        .def("group", &w_t::group, rir())
        .def("number", &w_t::number)
        .def("cb_op", &w_t::cb_op, ccr())
        .def("cb_op_is_tidy", &w_t::cb_op_is_tidy)
        .def("addl_generators_of_euclidean_normalizer",
          &w_t::addl_generators_of_euclidean_normalizer, (
            arg("flag_k2l"), arg("flag_l2n")))
        .def("expand_addl_generators_of_euclidean_normalizer",
          &w_t::expand_addl_generators_of_euclidean_normalizer, (
            arg("flag_k2l"), arg("flag_l2n")))
        .def("is_enantiomorphic", &w_t::is_enantiomorphic)
        .def("is_symmorphic", &w_t::is_symmorphic)
        .def("change_of_hand_op", &w_t::change_of_hand_op)
        .def("hall_symbol", &w_t::hall_symbol, (arg("tidy_cb_op")=true))
        .def("universal_hermann_mauguin_symbol",
          &w_t::universal_hermann_mauguin_symbol, (
            arg("tidy_cb_op")=true))
        .def("lookup_symbol", &w_t::lookup_symbol, (arg("ad_hoc_1992")=false))
        .def_pickle(space_group_type_wrappers())
      ;
    }
  };

} // namespace <anoymous>

  void wrap_space_group_type()
  {
    space_group_type_wrappers::wrap();
  }

}}} // namespace cctbx::sgtbx::boost_python
