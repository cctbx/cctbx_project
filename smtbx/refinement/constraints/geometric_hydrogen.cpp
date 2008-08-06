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
      .def("n_reparametrization_variables",
                    &wt::n_reparametrization_variables)
      .def("initialise_in_context", &wt::initialise_in_context)
      .def("place_constrained_scatterers",
           &wt::place_constrained_scatterers)
      .def("compute_gradients", &wt::compute_gradients)
      .def("apply_shifts", &wt::apply_shifts)
      ;
  }
};

template<typename FloatType=double,
         class XrayScattererType=xray::scatterer<> >
struct terminal_tetrahedral_XHn_wrapper
{
  typedef terminal_tetrahedral_XHn<FloatType, XrayScattererType,
                                   af::boost_python::flex_1d>
          wt;
  typedef typename wt::float_type float_type;

  static void wrap(char const *name) {
    using namespace boost::python;
    class_<wt> wrapped(name, no_init);
    common_wrapper::wrap(wrapped, name);
    wrapped
      .def(init<std::size_t, std::size_t, af::small<std::size_t,3>, float_type,
                float_type, bool, bool>
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
      .def(init<std::size_t, af::tiny<std::size_t,2>, af::tiny<std::size_t,2>,
                float_type, bool>
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
      .def(init<std::size_t, af::tiny<std::size_t,3>, std::size_t,
                float_type, bool>
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
      .def(init<std::size_t, af::tiny<std::size_t,2>, std::size_t,
                float_type, bool>
          ((arg("pivot"), arg("pivot_neighbours"), arg("hydrogen"),
            arg("bond_length"), arg("stretching")=false)))
      ;
  }
};

template<typename FloatType=double,
         class XrayScattererType=xray::scatterer<> >
struct terminal_trihedral_XH2_wrapper
{
  typedef terminal_trihedral_XH2<FloatType, XrayScattererType,
                                 af::boost_python::flex_1d>
          wt;
  typedef typename wt::float_type float_type;

  static void wrap(char const *name) {
    using namespace boost::python;
    class_<wt> wrapped(name, no_init);
    common_wrapper::wrap(wrapped, name);
    wrapped
      .def(init<std::size_t, std::size_t, std::size_t, af::tiny<std::size_t,2>,
                float_type, bool>
          ((arg("pivot"),
            arg("pivot_neighbour"),
            arg("pivot_neighbour_substituent"),
            arg("hydrogens"),
            arg("bond_length"), arg("stretching")=false)))
      ;
  }
};

template<typename FloatType=double,
         class XrayScattererType=xray::scatterer<> >
struct acetylenic_CH_wrapper
{
  typedef acetylenic_CH<FloatType, XrayScattererType, af::boost_python::flex_1d>
          wt;
  typedef typename wt::float_type float_type;

  static void wrap(char const *name) {
    using namespace boost::python;
    class_<wt> wrapped(name, no_init);
    common_wrapper::wrap(wrapped, name);
    wrapped
      .def(init<std::size_t, std::size_t, std::size_t, float_type, bool>
          ((arg("pivot"), arg("pivot_neighbour"), arg("hydrogen"),
            arg("bond_length"), arg("stretching")=false)))
      ;
  }
};

template<typename FloatType=double,
         class XrayScattererType=xray::scatterer<> >
struct polyhedral_BH_wrapper
{
  typedef polyhedral_BH<FloatType, XrayScattererType, af::boost_python::flex_1d>
          wt;
  typedef typename wt::float_type float_type;

  static void wrap(char const *name) {
    using namespace boost::python;
    class_<wt> wrapped(name, no_init);
    common_wrapper::wrap(wrapped, name);
    wrapped
      .def(init<std::size_t, af::small<std::size_t,5> const &, std::size_t, bool,
                float_type, bool>
          ((arg("pivot"), arg("pivot_neighbours"), arg("hydrogen"),
            arg("missing_fifth"),
            arg("bond_length"), arg("stretching")=false)))
      .def(init<std::size_t, af::small<std::size_t,5> const &, std::size_t,
                float_type, bool>
          ((arg("pivot"), arg("pivot_neighbours"), arg("hydrogen"),
            arg("bond_length"), arg("stretching")=false)))
      ;
  }
};

void wrap_geometric_hydrogen() {
  terminal_tetrahedral_XHn_wrapper<>::wrap("terminal_tetrahedral_XHn");
  secondary_CH2_wrapper<>::wrap("secondary_CH2");
  tertiary_CH_wrapper<>::wrap("tertiary_CH");
  aromatic_CH_or_amide_NH_wrapper<>::wrap("aromatic_CH_or_amide_NH");
  terminal_trihedral_XH2_wrapper<>::wrap("terminal_trihedral_XH2");
  acetylenic_CH_wrapper<>::wrap("acetylenic_CH");
  polyhedral_BH_wrapper<>::wrap("polyhedral_BH");
}

}}}} // smtbx::refinement::constraints::boost_python
