#include <iostream>
#include <vector>
#include <cctbx/fftbx/complex_to_complex_3d.h>
#include <cctbx/fftbx/real_to_complex_3d.h>

int main(void)
{
  std::size_t i;

  cctbx::fftbx::complex_to_complex<double> cfft(10);
  std::vector<std::complex<double> > vc(cfft.N());
  for(i=0;i<cfft.N();i++) {
    vc[i] = std::complex<double>(2.*i, 2.*i+1.);
  }
  cfft.forward(vc.begin());
  for(i=0;i<cfft.N();i++) {
    std::cout << vc[i].real() << " " << vc[i].imag() << std::endl;
  }
  cfft.backward(vc.begin());
  for(i=0;i<cfft.N();i++) {
    std::cout << vc[i].real() << " " << vc[i].imag() << std::endl;
  }

  cctbx::fftbx::real_to_complex<double> rfft(10);
  std::vector<double> vr(2 * rfft.Ncomplex());
  for(i=0;i<rfft.Nreal();i++) {
    vr[i] = 1.*i;
  }
  rfft.forward(vr.begin());
  for(i=0;i<2*rfft.Ncomplex();i++) {
    std::cout << vr[i] << std::endl;
  }
  rfft.backward(vr.begin());
  for(i=0;i<rfft.Nreal();i++) {
    std::cout << vr[i] << std::endl;
  }

  cctbx::fftbx::complex_to_complex_3d<double> cfft3d(2, 3, 5);
  cctbx::dimension_end<3> dim_c3d(cfft3d.N());
  std::vector<std::complex<double> > vc3d(dim_c3d.N1d());
  cctbx::ndim_accessor<
    cctbx::dimension_end<3>,
    std::vector<std::complex<double> >::iterator,
    std::vector<std::complex<double> >::value_type >
  c3dmap(dim_c3d, vc3d.begin());
  cfft3d.forward(c3dmap);
  cfft3d.backward(c3dmap);

  cctbx::fftbx::real_to_complex_3d<double> rfft3d(3, 4, 5);
  cctbx::dimension_end<3> dim_r3d(rfft3d.Mreal());
  std::vector<double> vr3d(dim_r3d.N1d());
  cctbx::ndim_accessor<
    cctbx::dimension_end<3>,
    std::vector<double>::iterator,
    std::vector<double>::value_type>
  r3dmap(dim_r3d, vr3d.begin());
  rfft3d.forward(r3dmap);
  rfft3d.backward(r3dmap);
#ifdef NEVER_DEFINED
#endif

  return 0;
}
