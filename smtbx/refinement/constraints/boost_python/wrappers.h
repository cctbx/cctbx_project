#ifndef SMTBX_REFINEMENT_CONSTRAINTS_BOOST_PYTHON_WRAPPERS_H
#define SMTBX_REFINEMENT_CONSTRAINTS_BOOST_PYTHON_WRAPPERS_H

#include <smtbx/boost_python/flex_fwd.h>
#include <boost/python/class.hpp>
#include <boost/python/with_custodian_and_ward.hpp>
#include <boost/python/return_internal_reference.hpp>

#include <scitbx/array_family/boost_python/passing_flex_by_reference.h>

#include <smtbx/refinement/constraints/constraint_array.h>

namespace smtbx { namespace refinement { namespace constraints {
namespace boost_python {

template <class WrappedType>
struct constrained_scatterers_wrapper
{
  typedef WrappedType wt;

  static boost::python::class_<wt> wrap(char const *name) {
    using namespace boost::python;

    return
    class_<wt>(name, no_init)
      .def(init<uctbx::unit_cell const &, // 2nd
                sgtbx::site_symmetry_table const &, // 3rd
                af::shared<typename wt::xray_scatterer_type>, // 4th
                typename wt::parameter_map_type const &, // 5th
                af::shared<xray::scatterer_flags> >( // 6th
           (arg("unit_cell"),
            arg("site_symmetry_table"),
            arg("scatterers"),
            arg("crystallographic_parameter_map"),
            arg("constraint_flags")))
           [with_custodian_and_ward<1,2,
            with_custodian_and_ward<1,3,
            with_custodian_and_ward<1,5> > >()])
      .add_property("already_constrained", &wt::already_constrained)
      .def("n_reparametrization_variables",
           &wt::n_reparametrization_variables)
      .def("compute_gradients", &wt::compute_gradients)
      .def("apply_shifts", &wt::apply_shifts)
      .def("place_constrained_scatterers", &wt::place_constrained_scatterers)
      ;
  }
};

template<class ConstraintType>
struct constraint_array_wrapper
{
  typedef constraint_array<ConstraintType, af::boost_python::flex_1d> wt;

  static void delitem(wt &self, std::size_t i) {
    self.erase(&self[i]);
  }

  static void wrap(char const *name) {
    using namespace boost::python;
    return_internal_reference<> rir;
    constrained_scatterers_wrapper<wt>::wrap(name)
      .def("__getitem__", &wt::operator[], rir)
      .def("__len__", &wt::size)
      .def("append", &wt::push_back)
      .def("__delitem__", delitem)
      ;
  }

};


}}}} // smtbx::refinement::constraints::boost_python

#endif // GUARD
