#include <boost/python/class.hpp>
#include <boost/python/def.hpp>

#include <scitbx/matrix/eigensystem.h>

namespace scitbx { namespace matrix { namespace boost_python {

  struct eigensystem_real_symmetric_wrappers
  {
    typedef eigensystem::real_symmetric<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("eigensystem_real_symmetric", no_init)
        .def(init<
          af::const_ref<double, af::c_grid<2> > const&, double, double>((
            arg("m"),
            arg("relative_epsilon")=1.e-10,
            arg("absolute_epsilon")=0)))
        .def(init<
          scitbx::sym_mat3<double> const&, double, double>((
            arg("m"),
            arg("relative_epsilon")=1.e-10,
            arg("absolute_epsilon")=0)))
        .def("min_abs_pivot", &w_t::min_abs_pivot)
        .def("vectors", &w_t::vectors)
        .def("values", &w_t::values)
        .def("generalized_inverse_as_packed_u",
          &w_t::generalized_inverse_as_packed_u)
      ;
    }
  };

  // simlar to time_dsyev_*()
  vec3<double>
  time_eigensystem_real_symmetric(
    sym_mat3<double> const& m, std::size_t n_repetitions)
  {
    SCITBX_ASSERT(n_repetitions % 2 == 0);
    vec3<double> result(0,0,0);
    for(std::size_t i=0;i<n_repetitions/2;i++) {
      result += vec3<double>(
        eigensystem::real_symmetric<>(m).values().begin());
      result -= vec3<double>(
        eigensystem::real_symmetric<>(m).values().begin());
    }
    return result / static_cast<double>(n_repetitions);
  }

  void wrap_eigensystem() {
    using namespace boost::python;
    eigensystem_real_symmetric_wrappers::wrap();
    def("time_eigensystem_real_symmetric", time_eigensystem_real_symmetric);
  }

}}}
