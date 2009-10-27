#include <smtbx/refinement/least_squares.h>
#include <smtbx/refinement/weighting_schemes.h>

#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/with_custodian_and_ward.hpp>


namespace smtbx { namespace refinement { namespace least_squares {
  namespace boost_python {

  template <typename FloatType>
  struct floating_origin_restraints_wrapper
  {
    typedef floating_origin_restraints<FloatType> wt;
    typedef typename wt::scatterer_t scatterer_t;
    typedef typename wt::scalar_t scalar_t;

    static boost::python::tuple singular_directions(wt const &self) {
      using namespace boost::python;
      switch (self.singular_directions.size()) {
        case 3:
          return make_tuple(self.singular_directions[0],
                            self.singular_directions[1],
                            self.singular_directions[2]);
        case 2:
          return make_tuple(self.singular_directions[0],
                            self.singular_directions[1]);
        case 1:
          return make_tuple(self.singular_directions[0]);
        case 0:
          return boost::python::tuple();
      }
      // It can't get there: just to prevent compiler warnings.
      return boost::python::tuple();
    }

    static void wrap(char const *name) {
      using namespace boost::python;
      class_<wt>(name, no_init)
        .def(init<sgtbx::space_group const &,
                  sgtbx::site_symmetry_table const &,
                  af::shared<scatterer_t> const &,
                  scalar_t>
             ((arg("space_group"),
               arg("site_symmetry_table"),
               arg("scatterers"),
               arg("floating_origin_restraint_relative_weight"))))
        .def("add_to", &wt::add_to, arg("normal_equations"))
        .add_property("singular_directions", singular_directions)
        ;
    }
  };

  template <typename FloatType,
            template<typename> class NormalEquations,
            template<typename> class WeightingScheme,
            class OneMillerIndexLinearisation,
            class ConstraintsType>
  struct normal_equation_building
  {
    static void
    delegate(NormalEquations<FloatType> const &normal_equations,
             af::const_ref<miller::index<> > const &miller_indices,
             af::const_ref<FloatType> const &data,
             af::const_ref<FloatType> const &sigmas,
             WeightingScheme<FloatType> const &weighting_scheme,
             FloatType scale_factor,
             OneMillerIndexLinearisation const &one_h_linearisation,
             ConstraintsType const &constraints)
    {
      smtbx::refinement::least_squares::build_normal_equations(
        const_cast<NormalEquations<FloatType> &>(normal_equations),
        miller_indices,
        data,
        sigmas,
        weighting_scheme,
        scale_factor,
        const_cast<OneMillerIndexLinearisation &>(one_h_linearisation),
        constraints);
    }

    static void wrap() {
      using namespace boost::python;
      def("build_normal_equations", delegate,
          (arg("normal_equations"),
           arg("miller_indices"),
           arg("data"),
           arg("sigmas"),
           arg("weighting_scheme"),
           arg("scale_factor"),
           arg("one_h_linearisation"),
           arg("constraints")));
    }
  };

  template <typename FloatType>
  struct special_position_constraints_wrapper
  {
    typedef special_position_constraints<FloatType> wt;
    typedef typename wt::scatterer_t scatterer_t;

    static void wrap(char const *name) {
      using namespace boost::python;
      class_<wt>(name, no_init)
        .def(init<uctbx::unit_cell const &,
                  sgtbx::site_symmetry_table const &,
                  af::shared<scatterer_t> const &>
             ((arg("unit_cell"),
               arg("site_symmetry_table"),
               arg("scatterers")))
              [with_custodian_and_ward<1, 2,
               with_custodian_and_ward<1, 3,
               with_custodian_and_ward<1, 4> > >()])
        .def_readonly("n_independent_params",
                      &wt::n_independent_params)
        .def_readonly("n_crystallographic_params",
                      &wt::n_crystallographic_params)
        .def("apply_shifts", &wt::apply_shifts,
             (arg("shifts"),
              arg("scatterers")))
        ;
  }
  };



  void wrap_least_squares() {
    using namespace boost::python;

    floating_origin_restraints_wrapper<
      double>::wrap("floating_origin_restraints");

    normal_equation_building<
      double,
      scitbx::lstbx::normal_equations_separating_scale_factor,
      unit_weighting,
      structure_factors::direct::one_h_linearisation::std_trigonometry<
        double,
        true,
        structure_factors::direct::one_h_linearisation::modulus_squared
      >,
      special_position_constraints<double>
    >::wrap();

    special_position_constraints_wrapper<
      double>::wrap("special_position_constraints");
  }


}}}}
