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

struct common_wrapper
{
  template<class wt>
  static void wrap(boost::python::class_<wt> &wrapped, char const *name) {
    std::string array_name = std::string(name) + std::string("_array");
    constraint_array_wrapper<wt>::wrap(array_name.c_str());

    wrapped
      .add_property("pivot", &wt::pivot)
      .add_property("hydrogens", &wt::hydrogens)
      .add_property("stretching", &wt::stretching, &wt::set_stretching)
      .add_property("bond_length", &wt::bond_length, &wt::set_bond_length)
      .add_property("initialise_in_context", &wt::initialise_in_context)
      .add_property("place_constrained_scatterers",
                    &wt::place_constrained_scatterers)
      .add_property("compute_gradients", &wt::compute_gradients)
      .add_property("apply_shifts", &wt::apply_shifts)
      ;
  }
};

template<typename FloatType=double,
         class XrayScattererType=xray::scatterer<> >
struct terminal_X_Hn_wrapper
{
  typedef terminal_X_Hn<FloatType, XrayScattererType, af::boost_python::flex_1d>
          wt;
  typedef typename wt::float_type float_type;

  static void wrap(char const *name) {
    using namespace boost::python;
    class_<wt> wrapped(name, no_init);
    common_wrapper::wrap(wrapped, name);
    wrapped
      .def(init<int, int, af::small<int,3>, float_type, float_type, bool, bool>
          ((arg("pivot"), arg("pivot_neighbour"), arg("hydrogens"),
            arg("azimuth"), arg("bond_length"), arg("rotating")=true,
            arg("stretching")=false)))
      .add_property("rotating", &wt::rotating, &wt::set_rotating)
      .add_property("azimuth", &wt::azimuth, &wt::set_azimuth)
      .add_property("local_cartesian_frame", &wt::local_cartesian_frame)
      ;
    }
};

template<typename FloatType=double,
         class XrayScattererType=xray::scatterer<> >
struct secondary_CH2_wrapper
{
  typedef secondary_CH2<FloatType, XrayScattererType, af::boost_python::flex_1d>
          wt;
  typedef typename wt::float_type float_type;

  static void wrap(char const *name) {
    using namespace boost::python;
    class_<wt> wrapped(name, no_init);
    common_wrapper::wrap(wrapped, name);
    wrapped
      .def(init<int, af::tiny<int,2>, af::tiny<int,2>, float_type, bool>
          ((arg("pivot"), arg("pivot_neighbours"), arg("hydrogens"),
            arg("bond_length"), arg("stretching")=false)))
      .def_readonly("theta0", wt::theta0)
      .def_readonly("dtheta_over_dXY_sq", &wt::dtheta_over_dXY_sq)
      ;
  }
};

template<typename FloatType=double,
         class XrayScattererType=xray::scatterer<> >
struct tertiary_CH_wrapper
{
  typedef tertiary_CH<FloatType, XrayScattererType, af::boost_python::flex_1d>
          wt;
  typedef typename wt::float_type float_type;

  static void wrap(char const *name) {
    using namespace boost::python;
    class_<wt> wrapped(name, no_init);
    common_wrapper::wrap(wrapped, name);
    wrapped
      .def(init<int, af::tiny<int,3>, int, float_type, bool>
          ((arg("pivot"), arg("pivot_neighbours"), arg("hydrogen"),
            arg("bond_length"), arg("stretching")=false)))
      ;
  }
};

template<typename FloatType=double,
         class XrayScattererType=xray::scatterer<> >
struct aromatic_CH_or_amide_NH_wrapper
{
  typedef aromatic_CH_or_amide_NH<FloatType, XrayScattererType,
                                  af::boost_python::flex_1d>
          wt;
  typedef typename wt::float_type float_type;

  static void wrap(char const *name) {
    using namespace boost::python;
    class_<wt> wrapped(name, no_init);
    common_wrapper::wrap(wrapped, name);
    wrapped
      .def(init<int, af::tiny<int,2>, int, float_type, bool>
          ((arg("pivot"), arg("pivot_neighbours"), arg("hydrogen"),
            arg("bond_length"), arg("stretching")=false)))
      ;
  }
};


void wrap_geometric_hydrogen() {
  terminal_X_Hn_wrapper<>::wrap("terminal_X_Hn");
  secondary_CH2_wrapper<>::wrap("secondary_CH2");
  tertiary_CH_wrapper<>::wrap("tertiary_CH");
  aromatic_CH_or_amide_NH_wrapper<>::wrap("aromatic_CH_or_amide_NH");
}

}}}} // smtbx::refinement::constraints::boost_python
