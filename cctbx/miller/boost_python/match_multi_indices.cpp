// -*- Mode: C++; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-

#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/miller/match_multi_indices.h>
#include <boost/python/class.hpp>

namespace cctbx { namespace miller { namespace boost_python {

namespace {

  struct match_multi_indices_wrappers
  {
    typedef match_multi_indices w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("match_multi_indices", no_init)
        .def(init<af::shared<index<> > const&,
                  af::shared<index<> > const&>(
                    (arg_("miller_indices_unique"),
                     arg_("miller_indices"))))
        .def("number_of_matches", &w_t::number_of_matches,
             arg_("i_array"))
        .def("have_singles", &w_t::have_singles)
        .def("pairs", &w_t::pairs)
        .def("singles", &w_t::singles,
             arg_("i_array"))
        .def("pair_selection", &w_t::pair_selection,
             arg_("i_array"))
        .def("single_selection", &w_t::single_selection,
             arg_("i_array"))
        .def("paired_miller_indices", &w_t::paired_miller_indices,
             arg_("i_array"));
    }
  };

} // namespace <anoymous>

  void wrap_match_multi_indices()
  {
    match_multi_indices_wrappers::wrap();
  }

}}} // namespace cctbx::miller::boost_python
