#include <iostream>
#include <vector>
#include <cctbx/fftbx/complex_to_complex_3d.h>
#include <cctbx/fftbx/real_to_complex_3d.h>
#include <cctbx/array_family/versa.h>

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
  cctbx::af::versa<std::complex<double>, cctbx::af::grid<3> >
  c3dmap(cctbx::af::grid<3>(cfft3d.N()));
  cfft3d.forward(c3dmap.ref());
  cfft3d.backward(c3dmap.ref());

  cctbx::fftbx::real_to_complex_3d<double> rfft3d(3, 4, 5);
  cctbx::af::versa<double, cctbx::af::grid<3> >
  r3dmap(cctbx::af::grid<3>(rfft3d.Mreal()));
  rfft3d.forward(r3dmap.ref());
  rfft3d.backward(r3dmap.ref());
#ifdef NEVER_DEFINED
#endif

  return 0;
}
