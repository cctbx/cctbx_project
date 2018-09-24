#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/enum.hpp>

#include <scitbx/matrix/svd.h>

namespace scitbx { namespace matrix { namespace boost_python {

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
                optional<scalar_t, int> >(
                  (arg("diagonal"), arg("off_diagonal"), arg("kind"),
                   arg("u"), arg("accumulate_u"),
                   arg("v"), arg("accumulate_v"),
                   arg("epsilon"),
                   arg("max_iteration_multiplier"))))
      .def("compute", &wt::compute)
      .def("sort", &wt::sort)
      .def("numerical_rank", &wt::numerical_rank)
      .def_readonly("has_converged", &wt::has_converged)
      ;
    }
  };

  template <typename FloatType>
  struct matrix_svd_decomposition_wrapper {
    typedef svd::decompose<FloatType> wt;

    static void wrap(char const *name) {
      using namespace boost::python;

      class_<wt>(name, no_init)
        .def(init<af::ref<FloatType, af::mat_grid> const &,
          optional<FloatType, bool, bool> >(
          (arg("matrix"),
            arg("crossover")=5./3,
            arg("accumulate_u") = false,
            arg("accumulate_v") = false)))
        .add_property("u", &wt::getU)
        .add_property("v", &wt::getV)
        .add_property("sigma", &wt::getSigma)
        .def("numerical_rank", &wt::numerical_rank)
        .def("reconstruct", &wt::reconstruct)
        ;
    }
  };

  void wrap_svd() {
    using namespace matrix::boost_python;
    bidiagonal_matrix_svd_decomposition_wrapper<double>::wrap(
      "svd_decomposition_of_bidiagonal_matrix");
    matrix_svd_decomposition_wrapper<double>::wrap(
      "svd_decompose");
    using namespace boost::python;
    def("reconstruct_svd", matrix::svd::reconstruct<double>);
  }

}}}
