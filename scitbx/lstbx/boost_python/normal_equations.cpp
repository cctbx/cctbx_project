#include <scitbx/lstbx/normal_equations.h>

#include <boost_adaptbx/optional_conversions.h>

#include <boost/python/class.hpp>
#include <boost/python/tuple.hpp>


namespace scitbx { namespace lstbx { namespace boost_python {

  template <typename FloatType>
  struct normal_equations_wrapper
  {
    typedef normal_equations<FloatType> wt;
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
        .add_property("normal_matrix_packed_u", &wt::normal_matrix)
        .add_property("right_hand_side", &wt::right_hand_side)
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
        .def(init<int>(arg("n_parameters")))
        .def("add_equation", add_equation,
             (arg("y_calc"), arg("grad_y_calc"), arg("y_obs"), arg("weight")))
        .def("optimised_scale_factor", &wt::optimised_scale_factor)
        .def("objective", &wt::objective)
        .def("gradient", &wt::gradient)
        .def("equations", &wt::equations);
    }
  };

  void wrap_normal_equations() {
    normal_equations_wrapper<double>::wrap("normal_equations");
    normal_equations_separating_scale_factor_wrapper<double>
      ::wrap("normal_equations_separating_scale_factor");
  }

}}}
