#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/return_value_policy.hpp>

#include <scitbx/random.h>
#include <scitbx/matrix/householder.h>

namespace scitbx { namespace matrix { namespace boost_python {

  template <class TriangularDecompositionType>
  struct householder_triangular_decomposition_wrapper
  {
    typedef TriangularDecompositionType wt;
    typedef typename wt::scalar_t scalar_t;

    static void wrap(char const *name) {
      using namespace boost::python;
      class_<wt>(name, no_init)
      .def(init<af::ref<scalar_t, af::mat_grid> const &, optional<bool> >(
           (arg("matrix"), arg("may_accumulate_q"))))
      .def("q", &wt::q, arg("thin")=true)
      .def("accumulate_q_in_place", &wt::accumulate_q_in_place)
      ;
    }
  };

  template <typename FloatType>
  struct householder_bidiagonalisation_wrapper
  {
    typedef householder::bidiagonalisation<FloatType> wt;
    typedef typename wt::scalar_t scalar_t;

    static void wrap(char const *name) {
      using namespace boost::python;
      class_<wt>(name, no_init)
      .def(init<af::ref<scalar_t, af::mat_grid> const &>())
      .def("u", &wt::u, arg("thin")=true)
      .def("v", &wt::v, arg("thin")=true)
      ;
    }
  };

  template <typename FloatType, class UniformRandomNumberGenerator>
  struct random_normal_matrix_generator_wrapper
  {
    typedef householder::random_normal_matrix_generator<
      FloatType, UniformRandomNumberGenerator> wt;

    static af::shared<std::size_t> get_state(wt const &self) {
      return self.normal_gen.engine().getstate();
    }

    static void set_state(wt &self,
                          af::const_ref<std::size_t> const &state)
    {
      self.normal_gen.engine().setstate(state);
    }

    static void wrap(char const *name) {
      using namespace boost::python;
      class_<wt>(name, no_init)
      .def(init<int, int>(args("rows", "columns")))
      .def("normal_matrix", &wt::normal_matrix)
      .def("matrix_with_singular_values", &wt::matrix_with_singular_values)
      .def("symmetric_matrix_with_eigenvalues",
           &wt::symmetric_matrix_with_eigenvalues)
      .add_property("state", get_state, set_state)
      ;
    }
  };

  void wrap_householder() {
    using namespace matrix::boost_python;
    householder_triangular_decomposition_wrapper<
    scitbx::matrix::householder::qr_decomposition<double> >::wrap(
      "householder_qr_decomposition");
    householder_triangular_decomposition_wrapper<
    scitbx::matrix::householder::lq_decomposition<double> >::wrap(
      "householder_lq_decomposition");
    householder_bidiagonalisation_wrapper<double>::wrap(
      "householder_bidiagonalisation");
    random_normal_matrix_generator_wrapper<
      double,
      boost_random::mt19937
    >::wrap("random_normal_matrix_generator");
  }

}}}
