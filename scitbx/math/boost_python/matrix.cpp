#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/def.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost_adaptbx/easy_overloads.h>
#include <scitbx/random.h>
#include <scitbx/matrix/householder.h>
#include <scitbx/matrix/svd.h>
#include <scitbx/matrix/tests.h>

namespace scitbx {

  namespace matrix { namespace boost_python {

    template <class TriangularDecompositionType>
    struct householder_triangular_decomposition_wrapper
    {
      typedef TriangularDecompositionType wt;
      typedef typename wt::scalar_t scalar_t;

      static void wrap(char const *name) {
        using namespace boost::python;
        return_value_policy<return_by_value> rbv;
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
        return_value_policy<return_by_value> rbv;
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
          .def(init<int, int>(
               args("rows", "columns")))
          .def("normal_matrix", &wt::normal_matrix)
          .def("matrix_with_singular_values", &wt::matrix_with_singular_values)
          .add_property("state", get_state, set_state)
          ;
        }
    };



    template <typename FloatType>
    struct bidiagonal_matrix_svd_decomposition_wrapper
    {
      typedef svd::bidiagonal_decomposition<FloatType> wt;
      typedef typename wt::scalar_t scalar_t;

      static void wrap(char const *name) {
        using namespace boost::python;
        enum_<svd::bidiagonal_kind>("bidiagonal_matrix_kind")
          .value("upper_diagonal", svd::upper_bidiagonal_kind)
          .value("lower_diagonal", svd::lower_bidiagonal_kind)
          ;
        class_<wt>(name, no_init)
          .def(init<af::ref<scalar_t> const&, af::ref<scalar_t> const&, int,
                    af::ref<scalar_t, af::mat_grid> const &, bool,
                    af::ref<scalar_t, af::mat_grid> const &, bool,
                    optional<scalar_t,
                             int> >(
                   (arg("diagonal"), arg("off_diagonal"), arg("kind"),
                    arg("u"), arg("accumulate_u"),
                    arg("v"), arg("accumulate_v"),
                    arg("epsilon"),
                    arg("max_iteration_multiplier"))))
          .def("compute", &wt::compute)
          .def("sort", &wt::sort, arg("reverse")=true)
          .def_readonly("has_converged", &wt::has_converged)
          ;
      }
    };

    template <typename FloatType>
    BOOST_ADAPTBX_FUNCTION_OVERLOADS(matrix_normality_ratio_overloads,
                                     normality_ratio<FloatType>, 1, 2,
                                     boost::python::args("a", "epsilon"));

    template <typename FloatType>
    BOOST_ADAPTBX_FUNCTION_OVERLOADS(matrix_equality_ratio_overloads,
                                     equality_ratio<FloatType>, 2, 3,
                                     boost::python::args("a", "b", "epsilon"));

  }} // matrix::boost_python


  namespace math { namespace boost_python {

    void wrap_matrix() {
      using namespace matrix::boost_python;
      householder_triangular_decomposition_wrapper<
        scitbx::matrix::householder::qr_decomposition<double> >::wrap(
                "householder_qr_decomposition");
      householder_triangular_decomposition_wrapper<
        scitbx::matrix::householder::lq_decomposition<double> >::wrap(
                "householder_lq_decomposition");
      householder_bidiagonalisation_wrapper<double>::wrap(
        "householder_bidiagonalisation");
      bidiagonal_matrix_svd_decomposition_wrapper<double>::wrap(
        "svd_decomposition_of_bidiagonal_matrix");
      matrix_normality_ratio_overloads<double>::wrap("matrix_normality_ratio");
      matrix_equality_ratio_overloads<double>::wrap("matrix_equality_ratio");

      random_normal_matrix_generator_wrapper<
        double,
        boost_random::mt19937
      >::wrap("random_normal_matrix_generator");

      using namespace boost::python;
      def("reconstruct_svd", matrix::svd::reconstruct<double>);
    }

  }} // math::boost_python

} // scitbx
