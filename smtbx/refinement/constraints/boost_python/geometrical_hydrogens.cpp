#include <boost/python/class.hpp>
#include <boost/python/implicit.hpp>

#include <scitbx/boost_python/container_conversions.h>

#include <smtbx/refinement/constraints/geometrical_hydrogens.h>

#include <sstream>

namespace smtbx { namespace refinement { namespace constraints {
namespace boost_python {

  template <int n_hydrogens, bool staggered>
  struct terminal_tetrahedral_xhn_sites_wrapper
  {
    typedef terminal_tetrahedral_xhn_sites<n_hydrogens, staggered> wt;

    static void wrap() {
      using namespace boost::python;
      std::ostringstream sname;
      if (staggered) sname << "staggered_";
      sname << "terminal_tetrahedral_xh";
      if (n_hydrogens > 1) sname << n_hydrogens;
      sname << "_site";
      if (n_hydrogens > 1) sname << "s";
      std::string name = sname.str();
      if (staggered) {
        class_<wt,
               bases<asu_parameter>,
               std::auto_ptr<wt> >(name.c_str(), no_init)
          .def(init<site_parameter *,
                    site_parameter *,
                    site_parameter *,
                    independent_scalar_parameter *,
                    af::tiny<typename wt::scatterer_type *,
                             n_hydrogens> const &>
               ((arg("pivot"), arg("pivot_neighbour"), arg("stagger_on"),
                 arg("length"),
                 arg("hydrogen"))))
          ;
      }
      else {
        class_<wt,
               bases<asu_parameter>,
               std::auto_ptr<wt> >(name.c_str(), no_init)
          .def(init<site_parameter *,
                    site_parameter *,
                    independent_scalar_parameter *,
                    independent_scalar_parameter *,
                    cart_t const &,
                    af::tiny<typename wt::scatterer_type *,
                             n_hydrogens> const &>
               ((arg("pivot"), arg("pivot_neighbour"),
                 arg("azimuth"), arg("length"),
                 arg("e_zero_azimuth"),
                 arg("hydrogen"))))
          ;
      }
      implicitly_convertible<std::auto_ptr<wt>, std::auto_ptr<parameter> >();
    }
  };

  struct angle_starting_tetrahedral_wrapper
  {
    typedef angle_starting_tetrahedral wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt,
             bases<scalar_parameter>,
             std::auto_ptr<wt> >("angle_starting_tetrahedral", no_init)
        .def(init<bool>(arg("variable")))
        ;
      implicitly_convertible<std::auto_ptr<wt>, std::auto_ptr<parameter> >();
    }
  };

  struct secondary_ch2_sites_wrapper
  {
    typedef secondary_ch2_sites wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt,
             bases<asu_parameter>,
             std::auto_ptr<wt> >("secondary_ch2_sites", no_init)
        .def(init<site_parameter *,
                site_parameter *,
                site_parameter *,
                independent_scalar_parameter *,
                angle_starting_tetrahedral *,
                wt::scatterer_type *,
                wt::scatterer_type *>
           ((arg("pivot"), arg("pivot_neighbour_0"), arg("pivot_neighbour_1"),
             arg("length"), arg("h_c_h_angle"),
             arg("hydrogen_0"), arg("hydrogen_1"))));
      implicitly_convertible<std::auto_ptr<wt>, std::auto_ptr<parameter> >();
    }
  };

  struct tertiary_ch_site_wrapper
  {
    typedef tertiary_ch_site wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt,
             bases<asu_parameter>,
             std::auto_ptr<wt> >("tertiary_ch_site", no_init)
        .def(init<site_parameter *,
                  site_parameter *,
                  site_parameter *,
                  site_parameter *,
                  independent_scalar_parameter *,
                  wt::scatterer_type *>
             ((arg("pivot"), arg("pivot_neighbour_0"), arg("pivot_neighbour_1"),
               arg("pivot_neighbour_2"), arg("length"),
               arg("hydrogen"))));
      implicitly_convertible<std::auto_ptr<wt>, std::auto_ptr<parameter> >();
    }
  };

  struct secondary_planar_xh_site_wrapper
  {
    typedef secondary_planar_xh_site wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt,
             bases<asu_parameter>,
             std::auto_ptr<wt> >("secondary_planar_xh_site", no_init)
        .def(init<site_parameter *,
                  site_parameter *,
                  site_parameter *,
                  independent_scalar_parameter *,
                  wt::scatterer_type *>
             ((arg("pivot"), arg("pivot_neighbour_0"), arg("pivot_neighbour_1"),
               arg("length"),
               arg("hydrogen"))));
      implicitly_convertible<std::auto_ptr<wt>, std::auto_ptr<parameter> >();
    }
  };

  struct terminal_planar_xh2_sites_wrapper
  {
    typedef terminal_planar_xh2_sites wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt,
             bases<asu_parameter>,
             std::auto_ptr<wt> >("terminal_planar_xh2_sites", no_init)
        .def(init<site_parameter *,
                  site_parameter *,
                  site_parameter *,
                  independent_scalar_parameter *,
                  wt::scatterer_type *, wt::scatterer_type *>
             ((arg("pivot"), arg("pivot_neighbour"),
               arg("pivot_neighbour_substituent"), arg("length"),
               arg("hydrogen_0"), arg("hydrogen_1"))));
      implicitly_convertible<std::auto_ptr<wt>, std::auto_ptr<parameter> >();
    }
  };


  struct terminal_linear_ch_site_wrapper
  {
    typedef terminal_linear_ch_site wt;

    static void wrap() {
      using namespace boost::python;
      class_<wt,
             bases<asu_parameter>,
             std::auto_ptr<wt> >("terminal_linear_ch_site", no_init)
        .def(init<site_parameter *,
                  site_parameter *,
                  independent_scalar_parameter *,
                  wt::scatterer_type *>
             ((arg("pivot"), arg("pivot_neighbour"), arg("length"),
               arg("hydrogen"))));
      implicitly_convertible<std::auto_ptr<wt>, std::auto_ptr<parameter> >();
    }
  };

  void wrap_geometrical_hydrogens() {
    {
      using namespace scitbx::boost_python::container_conversions;
      tuple_mapping_fixed_size<af::tiny<asu_parameter::scatterer_type *, 1> >();
      tuple_mapping_fixed_size<af::tiny<asu_parameter::scatterer_type *, 2> >();
      tuple_mapping_fixed_size<af::tiny<asu_parameter::scatterer_type *, 3> >();
    }
    //                                    #H  #staggered?
    terminal_tetrahedral_xhn_sites_wrapper<1, false>::wrap();
    terminal_tetrahedral_xhn_sites_wrapper<2, false>::wrap();
    terminal_tetrahedral_xhn_sites_wrapper<3, false>::wrap();
    terminal_tetrahedral_xhn_sites_wrapper<1, true>::wrap();
    terminal_tetrahedral_xhn_sites_wrapper<2, true>::wrap();
    terminal_tetrahedral_xhn_sites_wrapper<3, true>::wrap();

    angle_starting_tetrahedral_wrapper::wrap();
    secondary_ch2_sites_wrapper::wrap();
    tertiary_ch_site_wrapper::wrap();
    secondary_planar_xh_site_wrapper::wrap();
    terminal_planar_xh2_sites_wrapper::wrap();
    terminal_linear_ch_site_wrapper::wrap();
  }


}}}}
