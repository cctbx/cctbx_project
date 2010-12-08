#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/with_custodian_and_ward.hpp>

#include <smtbx/refinement/least_squares.h>
#include <smtbx/refinement/weighting_schemes.h>


namespace smtbx { namespace refinement { namespace least_squares {
  namespace boost_python {

  template <typename FloatType>
  struct floating_origin_restraints_wrapper
  {
    typedef floating_origin_restraints<FloatType> wt;
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
                  af::const_ref<constraints::scatterer_parameters> const &,
                  scitbx::sparse::matrix<FloatType> const &,
                  scalar_t>
             ((arg("space_group"),
               arg("all_scatterer_parameters"),
               arg("jacobian_transpose_matching_grad_fc"),
               arg("floating_origin_restraint_relative_weight"))))
        .def("add_to", &wt::add_to, arg("normal_equations"))
        .add_property("singular_directions", singular_directions)
        ;
    }
  };

  template <typename FloatType,
            template<typename> class NormalEquations,
            class OneMillerIndexLinearisation>
  struct normal_equation_building
  {
    typedef build_normal_equations<FloatType> wt;

    static void wrap(char const *name) {
      using namespace boost::python;
      class_<wt>(name, no_init)
        .def(init<NormalEquations<FloatType> &,
                  af::const_ref<miller::index<> > const &,
                  af::const_ref<FloatType> const &,
                  af::const_ref<FloatType> const &,
                  af::const_ref<std::complex<FloatType> > const &,
                  mainstream_shelx_weighting<FloatType> const &,
                  FloatType,
                  OneMillerIndexLinearisation &,
                  scitbx::sparse::matrix<FloatType> const &>
                  ())
        .def(init<NormalEquations<FloatType> &,
                  af::const_ref<miller::index<> > const &,
                  af::const_ref<FloatType> const &,
                  af::const_ref<FloatType> const &,
                  af::const_ref<std::complex<FloatType> > const &,
                  unit_weighting<FloatType> const &,
                  FloatType,
                  OneMillerIndexLinearisation &,
                  scitbx::sparse::matrix<FloatType> const &>
                  ())
        .def(init<NormalEquations<FloatType> &,
                  af::const_ref<miller::index<> > const &,
                  af::const_ref<FloatType> const &,
                  af::const_ref<FloatType> const &,
                  af::const_ref<std::complex<FloatType> > const &,
                  sigma_weighting<FloatType> const &,
                  FloatType,
                  OneMillerIndexLinearisation &,
                  scitbx::sparse::matrix<FloatType> const &>
                  ())
        .def("f_calc", &wt::f_calc)
        .def("weights", &wt::weights)
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
      structure_factors::direct::one_h_linearisation::std_trigonometry<
        double,
        true,
        structure_factors::direct::one_h_linearisation::modulus_squared
      >
    >::wrap("build_normal_equations");
  }


}}}}
