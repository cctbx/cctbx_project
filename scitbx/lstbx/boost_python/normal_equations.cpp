#include <scitbx/lstbx/normal_equations.h>

#include <boost/python/class.hpp>
#include <boost/python/tuple.hpp>


namespace scitbx { namespace lstbx { namespace boost_python {

  template <typename FloatType>
  struct normal_equations_separating_scale_factor_wrapper
  {
    typedef normal_equations_separating_scale_factor<FloatType> wt;
    typedef typename wt::scalar_t scalar_t;

    static void add_datum(wt &self,
                          scalar_t yc, af::const_ref<scalar_t> const &grad_yc,
                          scalar_t yo, scalar_t w)
    {
      self.add_datum(yc, grad_yc, yo, w);
    }

    static boost::python::tuple equations(wt &self) {
      typename wt::symmetric_matrix_t a;
      typename wt::vector_t b;
      boost::tie(a, b) = self.equations();
      af::flex_grid<> acc(a.accessor().size_1d());
      af::versa<scalar_t, af::flex_grid<> > a_(a, acc);
      return boost::python::make_tuple(a_, b);
    }

    static void wrap(char const *name) {
      using namespace boost::python;
      class_<wt>(name, no_init)
        .def(init<int>(arg("n_parameters")))
        .def("add_datum", add_datum,
             (arg("y_calc"), arg("grad_y_calc"), arg("y_obs"), arg("weight")))
        .def("optimised_scale_factor", &wt::optimised_scale_factor)
        .def("normal_matrix_packed_u_and_rhs", equations);
    }
  };

  void wrap_normal_equations() {
    normal_equations_separating_scale_factor_wrapper<double>
      ::wrap("normal_equations_separating_scale_factor");
  }

}}}
