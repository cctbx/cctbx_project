#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python.hpp>

#include <lapack_fem/dsyev.hpp>
#include <scitbx/sym_mat3.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <boost/scoped_array.hpp>

namespace lapack_fem { namespace boost_python {

  namespace af = scitbx::af;

  lapack_fem::common cmn;

  int
  dsyev_wrapper(
    std::string const& jobz,
    std::string const& uplo,
    af::ref<double, af::c_grid<2> > const& a,
    af::ref<double> const& w)
  {
    ASSERTBX(a.accessor().is_square());
    int n = a.accessor()[0];
    ASSERTBX(w.size() == n);
    int lwork = std::max(1, 3*n-1);
    boost::scoped_array<double> work(new double[lwork]);
    int info;
    dsyev(
      cmn,
      fem::str_cref(jobz.data(), jobz.size()),
      fem::str_cref(uplo.data(), uplo.size()),
      n,
      a[0],
      /*lda*/ n,
      w[0],
      work[0],
      lwork,
      info);
    return info;
  }

  scitbx::vec3<double>
  time_dsyev(
    scitbx::sym_mat3<double> const& m, std::size_t n_repetitions)
  {
    SCITBX_ASSERT(n_repetitions % 2 == 0);
    scitbx::vec3<double> result(0,0,0);
    for(std::size_t i=0;i<n_repetitions/2;i++) {
      for(std::size_t j=0;j<2;j++) {
        scitbx::vec3<double> w;
        scitbx::mat3<double> a(m);
        dsyev_wrapper("V", "U",
          af::ref<double, af::c_grid<2> >(a.begin(), af::c_grid<2>(3,3)),
          w.ref());
        if (j == 0) result += w;
        else        result -= w;
      }
    }
    return result / static_cast<double>(n_repetitions);
  }

  void
  wrap()
  {
    using namespace boost::python;
    def("lapack_dsyev", dsyev_wrapper, (
      arg("jobz"), arg("uplo"), arg("a"), arg("w")));
    def("time_lapack_dsyev", time_dsyev);
  }

}} // namespace lapack_fem::boost_python
