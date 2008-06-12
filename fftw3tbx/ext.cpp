#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <scitbx/array_family/boost_python/utils.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>
#include <string>

#include <fftw3.h>

namespace {

  namespace af = scitbx::af;

  void
  complex_to_complex_in_place(
    af::ref<std::complex<double> > const& data,
    int exp_sign)
  {
    SCITBX_ASSERT(exp_sign == FFTW_FORWARD || exp_sign == FFTW_BACKWARD);
    int n = static_cast<int>(data.size());
    fftw_complex *in = reinterpret_cast<fftw_complex*>(data.begin());
    fftw_plan p = fftw_plan_dft_1d(
      n, in, in, exp_sign, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
  }

  void
  complex_to_complex_3d_in_place(
    af::ref<std::complex<double>, af::c_grid<3> > const& data,
    int exp_sign)
  {
    SCITBX_ASSERT(exp_sign == FFTW_FORWARD || exp_sign == FFTW_BACKWARD);
    int nx = static_cast<int>(data.accessor()[0]);
    int ny = static_cast<int>(data.accessor()[1]);
    int nz = static_cast<int>(data.accessor()[2]);
    fftw_complex *in = reinterpret_cast<fftw_complex*>(data.begin());
    fftw_plan p = fftw_plan_dft_3d(
      nx, ny, nz, in, in, exp_sign, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
  }

  af::versa<std::complex<double>, af::flex_grid<> >
  real_to_complex_in_place(
    af::versa<double, af::flex_grid<> >& data)
  {
    af::boost_python::assert_0_based_1d(data.accessor());
    int m = static_cast<int>(data.accessor().all()[0]);
    int n = static_cast<int>(data.accessor().focus()[0]);
    SCITBX_ASSERT(m == 2*(n/2+1));
    double* in = data.begin();
    fftw_complex *out = reinterpret_cast<fftw_complex*>(in);
    fftw_plan p = fftw_plan_dft_r2c_1d(
      n, in, out, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
    return af::versa<std::complex<double>, af::flex_grid<> >(
      data.handle(), af::flex_grid<>(m/2));
  }

  af::versa<double, af::flex_grid<> >
  complex_to_real_in_place(
    af::versa<std::complex<double>, af::flex_grid<> >& data,
    int n)
  {
    af::boost_python::assert_0_based_1d(data.accessor());
    SCITBX_ASSERT(!data.accessor().is_padded());
    int m = static_cast<int>(data.accessor().all()[0]) * 2;
    SCITBX_ASSERT(m == 2*(n/2+1));
    fftw_complex *in = reinterpret_cast<fftw_complex*>(data.begin());
    double* out = reinterpret_cast<double*>(in);
    fftw_plan p = fftw_plan_dft_c2r_1d(
      n, in, out, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
    return af::versa<double, af::flex_grid<> >(
      data.handle(), af::flex_grid<>(m).set_focus(n));
  }

  af::versa<std::complex<double>, af::flex_grid<> >
  real_to_complex_3d_in_place(
    af::versa<double, af::flex_grid<> >& data)
  {
    af::boost_python::assert_0_based_3d(data.accessor());
    int mx = static_cast<int>(data.accessor().all()[0]);
    int my = static_cast<int>(data.accessor().all()[1]);
    int mz = static_cast<int>(data.accessor().all()[2]);
    int nx = static_cast<int>(data.accessor().focus()[0]);
    int ny = static_cast<int>(data.accessor().focus()[1]);
    int nz = static_cast<int>(data.accessor().focus()[2]);
    SCITBX_ASSERT(mx == nx);
    SCITBX_ASSERT(my == ny);
    SCITBX_ASSERT(mz == 2*(nz/2+1));
    double* in = data.begin();
    fftw_complex *out = reinterpret_cast<fftw_complex*>(in);
    fftw_plan p = fftw_plan_dft_r2c_3d(
      nx, ny, nz, in, out, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
    return af::versa<std::complex<double>, af::flex_grid<> >(
      data.handle(),
      af::flex_grid<>((af::adapt(af::tiny<int, 3>(mx,my,mz/2)))));
  }

  af::versa<double, af::flex_grid<> >
  complex_to_real_3d_in_place(
    af::versa<std::complex<double>, af::flex_grid<> >& data,
    af::tiny<int, 3> const& n)
  {
    af::boost_python::assert_0_based_3d(data.accessor());
    SCITBX_ASSERT(!data.accessor().is_padded());
    int mx = static_cast<int>(data.accessor().all()[0]);
    int my = static_cast<int>(data.accessor().all()[1]);
    int mz = static_cast<int>(data.accessor().all()[2]) * 2;
    int nx = n[0];
    int ny = n[1];
    int nz = n[2];
    SCITBX_ASSERT(mx == nx);
    SCITBX_ASSERT(my == ny);
    SCITBX_ASSERT(mz == 2*(nz/2+1));
    fftw_complex *in = reinterpret_cast<fftw_complex*>(data.begin());
    double* out = reinterpret_cast<double*>(in);
    fftw_plan p = fftw_plan_dft_c2r_3d(
      nx, ny, nz, in, out, FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
    return af::versa<double, af::flex_grid<> >(
      data.handle(),
      af::flex_grid<>((af::adapt(af::tiny<int, 3>(nx,ny,mz))))
        .set_focus(af::adapt(n)));
  }

  void
  wrap_fftw3()
  {
    using namespace boost::python;
    scope().attr("fftw_version") = std::string(fftw_version);
    def("complex_to_complex_in_place", complex_to_complex_in_place, (
      arg_("data"), arg_("exp_sign")));
    def("complex_to_complex_3d_in_place", complex_to_complex_3d_in_place, (
      arg_("data"), arg_("exp_sign")));
    def("real_to_complex_in_place", real_to_complex_in_place, (
      arg_("data")));
    def("complex_to_real_in_place", complex_to_real_in_place, (
      arg_("data"), arg_("n")));
    def("real_to_complex_3d_in_place", real_to_complex_3d_in_place, (
      arg_("data")));
    def("complex_to_real_3d_in_place", complex_to_real_3d_in_place, (
      arg_("data"), arg_("n")));
  }

} // namespace <anonymous>

BOOST_PYTHON_MODULE(fftw3tbx_ext)
{
  wrap_fftw3();
}
