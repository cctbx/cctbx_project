#include <scitbx/lstbx/normal_equations.h>

#include <boost_adaptbx/optional_conversions.h>

#include <boost/python/class.hpp>
#include <boost/python/tuple.hpp>


namespace scitbx { namespace lstbx { namespace boost_python {

  template <typename FloatType>
  struct normal_equations_wrapper
  {
    typedef normal_equations<FloatType> wt;
    typedef typename wt::scalar_t scalar_t;
    typedef typename wt::symmetric_matrix_t symmetric_matrix_t;
    typedef typename wt::vector_t vector_t;

    static void wrap(char const *name) {
      using namespace boost::python;

      boost_adaptbx::optional_conversions::to_and_from_python<
        boost::optional<wt> >();

      class_<wt>(name, no_init)
        .def(init<int>(arg("n_parameters")))
        .def(init<symmetric_matrix_t const &, vector_t const &>(
             (arg("normal_matrix"), arg("right_hand_side"))))
        .def("add_equation",
             (void (wt::*)(scalar_t b, af::const_ref<scalar_t> const &, scalar_t))
             &wt::add_equation,
             (arg("right_hand_side"), arg("design_matrix_row"), arg("weight")))
        .def("add_equations", &wt::add_equations,
             (arg("right_hand_side"), arg("design_matrix"), arg("weights"),
              arg("negate_right_hand_side")=false))
        .def("solve", &wt::solve)
        .add_property("solved", &wt::solved)
        /* We use 'def' instead of add_property for those because they may
           throw if called on an instanced which is not finalised.
           On the Python side, an attribute lookup which may throw is a
           source of confusion (e.g. hasattr does not work correctly for those).
         */
        .def("normal_matrix_packed_u", &wt::normal_matrix)
        .def("right_hand_side", &wt::right_hand_side)
        .def("cholesky_factor_packed_u", &wt::cholesky_factor)
        .def("solution", &wt::solution)
        ;
    }
  };


  template <typename FloatType>
  struct normal_equations_separating_scale_factor_wrapper
  {
    typedef normal_equations_separating_scale_factor<FloatType> wt;
    typedef typename wt::scalar_t scalar_t;

    static void add_equation(wt &self,
                             scalar_t yc, af::const_ref<scalar_t> const &grad_yc,
                             scalar_t yo, scalar_t w)
    {
      self.add_equation(yc, grad_yc, yo, w);
    }

    static void wrap(char const *name) {
      using namespace boost::python;
      class_<wt>(name, no_init)
        .def(init<int, bool>((arg("n_parameters"), arg("normalised")=true)))
        .def("add_equation", add_equation,
             (arg("y_calc"), arg("grad_y_calc"), arg("y_obs"), arg("weight")))
        .def("finalise", &wt::finalise)
        .add_property("finalised", &wt::finalised)
        .def("reset", &wt::reset)
        /* We use 'def' instead of add_property for those because they may
           throw if called on an instanced which is not finalised.
           On the Python side, an attribute lookup which may throw is a
           source of confusion (e.g. hasattr does not work correctly for those).
         */
        .def("optimal_scale_factor", &wt::optimal_scale_factor)
        .def("sum_w_yo_sq", &wt::sum_w_yo_sq)
        .def("objective", &wt::objective)
        .def("step_equations", &wt::step_equations)
        ;
    }
  };

  void wrap_normal_equations() {
    normal_equations_wrapper<double>::wrap("normal_equations");
    normal_equations_separating_scale_factor_wrapper<double>
      ::wrap("normal_equations_separating_scale_factor");
  }

}}}
