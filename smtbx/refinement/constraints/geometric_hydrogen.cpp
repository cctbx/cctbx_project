#include <boost/python/class.hpp>
#include <boost/python/slice.hpp>
#include <boost/python/return_value_policy.hpp>

#include <scitbx/array_family/boost_python/passing_flex_by_reference.h>

#include <cctbx/xray/scatterer.h>

#include <smtbx/refinement/constraints/constraints_array_bp.h>
#include <smtbx/refinement/constraints/geometric_hydrogen.h>

namespace smtbx { namespace refinement { namespace constraints {
namespace boost_python {

template<typename FloatType, class XrayScattererType>
struct stretchable_rotatable_riding_ch3_wrapper
{
  typedef stretchable_rotatable_riding_terminal_X_Hn<FloatType,
                                                     XrayScattererType,
                                                     af::boost_python::flex_1d>
           wt;

  typedef typename wt::float_type float_type;

  static void wrap() {
    using namespace boost::python;
    std::string name("stretchable_rotatable_riding_terminal_X_Hn");
    std::string array_name = name + std::string("_array");

    constraint_array_wrapper<wt>::wrap(array_name.c_str());

    class_<wt>(name.c_str(), no_init)
      .def(init<int, int, af::small<int, 3>,
                float_type, float_type,
                bool, bool>((
            arg("pivot"), arg("pivot_neighbour"),
            arg("hydrogens"),
            arg("azimuth"), arg("bond_length"),
            arg("rotating")=true, arg("stretching")=false)))
      .add_property("pivot", &wt::pivot)
      .add_property("hydrogens", &wt::hydrogens)
      .add_property("rotating", &wt::rotating, &wt::set_rotating)
      .add_property("stretching", &wt::stretching, &wt::set_stretching)
      .add_property("azimuth", &wt::azimuth, &wt::set_azimuth)
      .add_property("bond_length", &wt::bond_length, &wt::set_bond_length)
      .add_property("local_cartesian_frame", &wt::local_cartesian_frame)
      .def("compute_gradients", &wt::compute_gradients)
      .def("apply_shifts", &wt::apply_shifts)
      .def("place_constrained_scatterers", &wt::place_constrained_scatterers)
      ;
  }
};

void wrap_geometric_hydrogen() {
  stretchable_rotatable_riding_ch3_wrapper<double, xray::scatterer<> >::wrap();
}

}}}} // smtbx::refinement::constraints::boost_python
