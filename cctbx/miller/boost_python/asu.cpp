/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Created (rwgk)
 */

#include <cctbx/miller/asu.h>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>

namespace cctbx { namespace miller { namespace boost_python {

namespace {

  struct index_table_layout_adaptor_wrappers
  {
    typedef index_table_layout_adaptor w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t,bases<sym_equiv_index> >("index_table_layout_adaptor",no_init)
        .def("h", &w_t::h)
        .def("i_column", &w_t::i_column)
      ;
    }
  };

  struct asym_index_wrappers
  {
    typedef asym_index w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t, bases<sym_equiv_index> >("asym_index", no_init)
        .def(init<sgtbx::space_group const&,
                  sgtbx::reciprocal_space::asu const&,
                  index<> const&>())
        .def(init<sgtbx::space_group const&,
                  index<> const&>())
        .def(init<sym_equiv_indices const&>())
        .def("one_column", &w_t::one_column)
        .def("two_column", &w_t::two_column)
      ;
    }
  };

  template <typename ValueType>
  struct map_to_asu_wrappers
  {
    static
    void
    no_bool(
      sgtbx::space_group_type const& sg_type,
      bool anomalous_flag,
      af::ref<index<> > const& miller_indices,
      af::ref<ValueType> const& data_array)
    {
      map_to_asu(sg_type, anomalous_flag, miller_indices, data_array);
    }

    static
    void
    with_bool(
      sgtbx::space_group_type const& sg_type,
      bool anomalous_flag,
      af::ref<index<> > const& miller_indices,
      af::ref<ValueType> const& data_array,
      bool deg)
    {
      map_to_asu(sg_type, anomalous_flag, miller_indices, data_array, deg);
    }
  };

} // namespace <anoymous>

  void wrap_asu()
  {
    using namespace boost::python;
    index_table_layout_adaptor_wrappers::wrap();
    asym_index_wrappers::wrap();
    def("map_to_asu", map_to_asu_wrappers<double>::no_bool);
    def("map_to_asu", map_to_asu_wrappers<double>::with_bool);
    def("map_to_asu", map_to_asu_wrappers<std::complex<double> >::no_bool);
    def("map_to_asu", map_to_asu_wrappers<hendrickson_lattman<> >::no_bool);
  }

}}} // namespace cctbx::miller::boost_python
