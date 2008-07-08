#ifndef SMTBX_CONSTRAINTS_ARRAY_BP
#define SMTBX_CONSTRAINTS_ARRAY_BP

#include <boost/python/class.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <smtbx/refinement/constraints/constraint_array.h>

namespace smtbx { namespace refinement { namespace constraints {
namespace boost_python {

template<class ConstraintType>
struct constraint_array_wrapper
{
  typedef constraint_array<ConstraintType, af::boost_python::flex_1d> wt;

  static void delitem(wt &self, int i) {
    self.erase(&self[i]);
  }

  static void wrap(char const *name) {
    using namespace boost::python;
    return_internal_reference<> rir;
    class_<wt>(name, no_init)
      .def(init<uctbx::unit_cell const &, // 2th
                sgtbx::site_symmetry_table const &, // 3nd
                af::shared<typename wt::xray_scatterer_type>, // 4rd
                typename wt::parameter_map_type const &, // 5th
                af::shared<xray::scatterer_flags> >( // 6th
           (arg_("unit_cell"),
            arg_("site_symmetry_table"),
            arg_("scatterers"),
            arg_("crystallographic_parameter_map"),
            arg_("constraint_flags")))
           [with_custodian_and_ward<1,2,
            with_custodian_and_ward<1,3,
            with_custodian_and_ward<1,5> > >()])
      .def("__getitem__", &wt::operator[], rir)
      .def("__len__", &wt::size)
      .def("append", &wt::push_back)
      .def("__delitem__", delitem)
      .def("compute_gradients", &wt::compute_gradients)
      .def("apply_shifts", &wt::apply_shifts)
      .def("place_constrained_scatterers", &wt::place_constrained_scatterers)
      .add_property("already_constrained", &wt::already_constrained)
      ;
  }
};

}}}}

#endif // GUARD
