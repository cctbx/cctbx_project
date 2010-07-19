#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python.hpp>

#include <lapack_fem/dsyev.hpp>
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

  void
  wrap()
  {
    using namespace boost::python;
    def("lapack_dsyev", dsyev_wrapper, (
      arg("jobz"), arg("uplo"), arg("a"), arg("w")));
  }

}} // namespace lapack_fem::boost_python
