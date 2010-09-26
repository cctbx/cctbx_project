#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python.hpp>

#include <scitbx/sym_mat3.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <boost/scoped_array.hpp>

namespace boost_python_meta_ext { struct holder {}; }

extern "C" {

  void
  dgesdd_(
    char const* jobz,
    int const* m,
    int const* n,
    double* a,
    int const* lda,
    double* s,
    double* u,
    int const* ldu,
    double* vt,
    int const* ldvt,
    double* work,
    int const* lwork,
    int* iwork,
    int* info,
    int jobz_len);

  void
  dgesvd_(
    char const* jobu,
    char const* jobvt,
    int const* m,
    int const* n,
    double* a,
    int const* lda,
    double* s,
    double* u,
    int const* ldu,
    double* vt,
    int const* ldvt,
    double* work,
    int const* lwork,
    int* info,
    int jobu_len,
    int jobvt_len);

  void
  dsyev_(
    char const* jobz,
    char const* uplo,
    int const* n,
    double* a,
    int const* lda,
    double* w,
    double* work,
    int const* lwork,
    int* info,
    int jobz_len,
    int uplo_len);

} // extern "C"

#if defined(SCITBX_LAPACK_FEM)
#  include <lapack_fem/selected.hpp>
#endif

namespace scitbx { namespace lapack { namespace boost_python {

#if defined(SCITBX_LAPACK_FEM)
  lapack_fem::common cmn;
#endif

  boost::python::object
  dgesdd_wrapper(
    af::ref<double, af::c_grid<2> > const& a,
    bool use_fortran=false)
  {
    int m = a.accessor()[1];
    int n = a.accessor()[0];
    SCITBX_ASSERT(m > 0);
    SCITBX_ASSERT(n > 0);
    boost::python::object result;
#if defined(SCITBX_LAPACK_FEM) || defined(SCITBX_LAPACK_FOR)
    int p = std::min(m,n);
    af::shared<double> s(p, 0.);
    af::versa<double, af::c_grid<2> > u(af::c_grid<2>(p, m), 0.);
    af::versa<double, af::c_grid<2> > vt(af::c_grid<2>(n, p), 0.);
    boost::scoped_array<int> iwork(new int[8*p]);
    int lwork = -1;
    int info;
    for(unsigned i_pass=0;i_pass<2;i_pass++) {
      boost::scoped_array<double> work(new double[std::max(1,lwork)]);
#endif
      if (!use_fortran) {
#if defined(SCITBX_LAPACK_FEM)
        lapack_fem::dgesdd(
          cmn,
          /*jobz*/ "S",
          m,
          n,
          a[0],
          /*lda*/ m,
          s[0],
          u[0],
          /*ldu*/ m,
          vt[0],
          /*ldvt*/ p,
          work[0],
          lwork,
          iwork[0],
          info);
#else
        return result;
#endif
      }
      else {
#if defined(SCITBX_LAPACK_FOR)
        dgesdd_(
          /*jobz*/ "S",
          &m,
          &n,
          &a[0],
          /*lda*/ &m,
          &s[0],
          &u[0],
          /*ldu*/ &m,
          &vt[0],
          /*ldvt*/ &p,
          &work[0],
          &lwork,
          &iwork[0],
          &info,
          /*jobz_len*/ 1);
#else
        return result;
#endif
      }
#if defined(SCITBX_LAPACK_FEM) || defined(SCITBX_LAPACK_FOR)
      if (i_pass == 0) {
        lwork = static_cast<int>(work[0]);
      }
    }
    result = boost::python::object(boost_python_meta_ext::holder());
    result.attr("s") = s;
    result.attr("u") = u;
    result.attr("vt") = vt;
    result.attr("info") = info;
#endif
    return result;
  }

  boost::python::object
  dgesvd_wrapper(
    af::ref<double, af::c_grid<2> > const& a,
    bool use_fortran=false)
  {
    int m = a.accessor()[1];
    int n = a.accessor()[0];
    SCITBX_ASSERT(m > 0);
    SCITBX_ASSERT(n > 0);
    boost::python::object result;
#if defined(SCITBX_LAPACK_FEM) || defined(SCITBX_LAPACK_FOR)
    int p = std::min(m,n);
    af::shared<double> s(p, 0.);
    af::versa<double, af::c_grid<2> > u(af::c_grid<2>(p, m), 0.);
    af::versa<double, af::c_grid<2> > vt(af::c_grid<2>(n, p), 0.);
    int lwork = -1;
    int info;
    for(unsigned i_pass=0;i_pass<2;i_pass++) {
      boost::scoped_array<double> work(new double[std::max(1,lwork)]);
#endif
      if (!use_fortran) {
#if defined(SCITBX_LAPACK_FEM)
        lapack_fem::dgesvd(
          cmn,
          /*jobu*/ "S",
          /*jobvt*/ "S",
          m,
          n,
          a[0],
          /*lda*/ m,
          s[0],
          u[0],
          /*ldu*/ m,
          vt[0],
          /*ldvt*/ p,
          work[0],
          lwork,
          info);
#else
        return result;
#endif
      }
      else {
#if defined(SCITBX_LAPACK_FOR)
        dgesvd_(
          /*jobu*/ "S",
          /*jobvt*/ "S",
          &m,
          &n,
          &a[0],
          /*lda*/ &m,
          &s[0],
          &u[0],
          /*ldu*/ &m,
          &vt[0],
          /*ldvt*/ &p,
          &work[0],
          &lwork,
          &info,
          /*jobu_len*/ 1,
          /*jobvt_len*/ 1);
#else
        return result;
#endif
      }
#if defined(SCITBX_LAPACK_FEM) || defined(SCITBX_LAPACK_FOR)
      if (i_pass == 0) {
        lwork = static_cast<int>(work[0]);
      }
    }
    result = boost::python::object(boost_python_meta_ext::holder());
    result.attr("s") = s;
    result.attr("u") = u;
    result.attr("vt") = vt;
    result.attr("info") = info;
#endif
    return result;
  }

  int
  dsyev_wrapper(
    std::string const& jobz,
    std::string const& uplo,
    af::ref<double, af::c_grid<2> > const& a,
    af::ref<double> const& w,
    bool use_fortran=false)
  {
    SCITBX_ASSERT(a.accessor().is_square());
    int n = a.accessor()[0];
    SCITBX_ASSERT(w.size() == n);
    int info = 99;
#if defined(SCITBX_LAPACK_FEM) || defined(SCITBX_LAPACK_FOR)
    int lwork = -1;
    for(unsigned i_pass=0;i_pass<2;i_pass++) {
      boost::scoped_array<double> work(new double[std::max(1,lwork)]);
#endif
      if (!use_fortran) {
#if defined(SCITBX_LAPACK_FEM)
        lapack_fem::dsyev(
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
#endif
      }
      else {
#if defined(SCITBX_LAPACK_FOR)
        dsyev_(
          jobz.data(),
          uplo.data(),
          &n,
          &a[0],
          /*lda*/ &n,
          &w[0],
          &work[0],
          &lwork,
          &info,
          jobz.size(),
          uplo.size());
#endif
      }
#if defined(SCITBX_LAPACK_FEM) || defined(SCITBX_LAPACK_FOR)
      if (i_pass == 0) {
        lwork = static_cast<int>(work[0]);
      }
    }
#endif
    return info;
  }

  // simlar to time_eigensystem_real_symmetric()
  scitbx::vec3<double>
  time_dsyev(
    scitbx::sym_mat3<double> const& m,
    std::size_t n_repetitions,
    bool use_fortran=false)
  {
    SCITBX_ASSERT(n_repetitions % 2 == 0);
    scitbx::vec3<double> result(0,0,0);
    int info = 99;
    for(std::size_t i=0;i<n_repetitions/2;i++) {
      for(std::size_t j=0;j<2;j++) {
        scitbx::vec3<double> w;
        scitbx::mat3<double> a(m);
        info = dsyev_wrapper("V", "U",
          af::ref<double, af::c_grid<2> >(a.begin(), af::c_grid<2>(3,3)),
          w.ref(), use_fortran);
        if (j == 0) result += w;
        else        result -= w;
      }
    }
    SCITBX_ASSERT(info == 0);
    return result / static_cast<double>(n_repetitions);
  }

  bool
  fem_is_available()
  {
#if defined(SCITBX_LAPACK_FEM)
    return true;
#else
    return false;
#endif
  }

  bool
  for_is_available()
  {
#if defined(SCITBX_LAPACK_FOR)
    return true;
#else
    return false;
#endif
  }

  void
  wrap()
  {
    using namespace boost::python;

    def("fem_is_available", fem_is_available);
    def("for_is_available", for_is_available);

    def("lapack_dgesdd", dgesdd_wrapper, (
      arg("a"), arg("use_fortran")=false));
    def("lapack_dgesvd", dgesvd_wrapper, (
      arg("a"), arg("use_fortran")=false));
    def("lapack_dsyev", dsyev_wrapper, (
      arg("jobz"), arg("uplo"), arg("a"), arg("w"), arg("use_fortran")=false));

    def("time_lapack_dsyev", time_dsyev, (
      arg("m"), arg("n_repetitions"), arg("use_fortran")=false));
  }

}}} // namespace scitbx::lapack::boost_python
