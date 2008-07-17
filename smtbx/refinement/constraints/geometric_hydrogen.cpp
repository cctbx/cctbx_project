#include <smtbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/slice.hpp>
#include <boost/python/return_value_policy.hpp>

#include <scitbx/array_family/boost_python/passing_flex_by_reference.h>

#include <cctbx/xray/scatterer.h>

#include <smtbx/refinement/constraints/constraints_array_bp.h>
#include <smtbx/refinement/constraints/geometric_hydrogen.h>

namespace smtbx { namespace refinement { namespace constraints {
namespace boost_python {



template<template<typename FloatType_bis,
                  class XrayScattererType_bis,
                  template<typename FloatType_ter> class SharedArray1D_>
         class ConstraintTemplate,
         typename FloatType=double,
         class XrayScattererType=xray::scatterer<> >
struct few_scatterer_constraint_wrapper
{
  typedef ConstraintTemplate<FloatType,
                             XrayScattererType,
                             af::boost_python::flex_1d> wt;

  boost::python::class_<wt> wrapped;

  few_scatterer_constraint_wrapper(char const *name)
    : wrapped(name, boost::python::no_init)
  {
    using namespace boost::python;
    std::string array_name = std::string(name) + std::string("_array");
    constraint_array_wrapper<wt>::wrap(array_name.c_str());

    wrapped
      .def("compute_gradients", &wt::compute_gradients)
      .def("apply_shifts", &wt::apply_shifts)
      .def("place_constrained_scatterers", &wt::place_constrained_scatterers)
      ;
  }
};

template<template<typename FloatType_bis,
                  class XrayScattererType_bis,
                  template<typename FloatType_ter> class SharedArray1D>
         class ConstraintTemplate,
         typename FloatType=double,
         class XrayScattererType=xray::scatterer<> >
struct hydrogen_constraint_wrapper
  : few_scatterer_constraint_wrapper<ConstraintTemplate,
                                     FloatType,
                                     XrayScattererType>
{
  typedef few_scatterer_constraint_wrapper<ConstraintTemplate,
                                           FloatType,
                                           XrayScattererType> base_t;
  typedef typename base_t::wt wt;
  using base_t::wrapped;

  hydrogen_constraint_wrapper(char const *name)
    : base_t(name)
  {
    using namespace boost::python;
    wrapped
      .add_property("pivot", &wt::pivot)
      .add_property("hydrogens", &wt::hydrogens)
      .add_property("bond_length", &wt::bond_length, &wt::set_bond_length)
      .add_property("stretching", &wt::stretching, &wt::set_stretching)
      ;
  }
};

template<template<typename FloatType_bis,
                  class XrayScattererType_bis,
                  template<typename FloatType_ter> class SharedArray1D>
         class ConstraintTemplate,
         int n_neighbours, int n_hydrogens,
         typename FloatType=double,
         class XrayScattererType=xray::scatterer<> >
struct common_hydrogen_constraint_wrapper
  : hydrogen_constraint_wrapper<ConstraintTemplate,
                                FloatType,
                                XrayScattererType>
{
  typedef hydrogen_constraint_wrapper<ConstraintTemplate,
                                      FloatType,
                                      XrayScattererType> base_t;
  typedef typename base_t::wt wt;
  typedef typename wt::float_type float_type;
  using base_t::wrapped;

  common_hydrogen_constraint_wrapper(char const *name)
    : base_t(name)
  {
    using namespace boost::python;
    wrapped
      .def(init<int,
                af::tiny<int, n_neighbours>,
                af::tiny<int, n_hydrogens>,
                float_type,
                bool>((
            arg("pivot"), arg("pivot_neighbours"),
            arg("hydrogens"),
            arg("bond_length"),
            arg("stretching")=false)))
      ;
  }
};

template<typename FloatType=double,
         class XrayScattererType=xray::scatterer<> >
struct terminal_X_Hn_wrapper
  : hydrogen_constraint_wrapper<terminal_X_Hn,
                                FloatType,
                                XrayScattererType>
{
  typedef hydrogen_constraint_wrapper<terminal_X_Hn,
                                      FloatType,
                                      XrayScattererType> base_t;
  typedef typename base_t::wt wt;
  typedef typename wt::float_type float_type;
  using base_t::wrapped;

  terminal_X_Hn_wrapper(char const *name)
    : base_t(name)
  {
    using namespace boost::python;
    wrapped
      .def(init<int, int, af::small<int, 3>, float_type, float_type, bool, bool>
          ((arg("pivot"), arg("pivot_neighbour"),
            arg("hydrogens"),
            arg("azimuth"), arg("bond_length"),
            arg("rotating")=true, arg("stretching")=false)))
      .add_property("rotating", &wt::rotating, &wt::set_rotating)
      .add_property("azimuth", &wt::azimuth, &wt::set_azimuth)
      .add_property("local_cartesian_frame", &wt::local_cartesian_frame)
      ;
  }
};

template<typename FloatType=double,
         class XrayScattererType=xray::scatterer<> >
struct secondary_CH2_wrapper
  : common_hydrogen_constraint_wrapper<secondary_CH2, 2, 2,
                                       FloatType,
                                       XrayScattererType>
{
  typedef common_hydrogen_constraint_wrapper<secondary_CH2, 2, 2,
                                             FloatType,
                                             XrayScattererType>
          base_t;
  typedef typename base_t::wt wt;
  typedef typename wt::float_type float_type;
  using base_t::wrapped;

  secondary_CH2_wrapper(char const *name)
    : base_t(name)
  {
    using namespace boost::python;
    wrapped
      .def_readwrite("theta0", wt::theta0)
      .def_readwrite("dtheta_over_dXY_sq", wt::dtheta_over_dXY_sq)
      ;
  }
};

void wrap_geometric_hydrogen() {
  terminal_X_Hn_wrapper<>("terminal_X_Hn");
  secondary_CH2_wrapper<>("secondary_CH2");
  common_hydrogen_constraint_wrapper<tertiary_CH, 3, 1>("tertiary_CH");
}

}}}} // smtbx::refinement::constraints::boost_python
