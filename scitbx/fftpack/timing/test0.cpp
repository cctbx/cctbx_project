/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Ported from cctbx (rwgk)
     2001 Nov: Created (R.W. Grosse-Kunstleve)
 */

#include <iostream>
#include <vector>
#include <scitbx/fftpack/complex_to_complex_3d.h>
#include <scitbx/fftpack/real_to_complex_3d.h>
#include <scitbx/array_family/versa.h>

int main(void)
{
  std::size_t i;

  scitbx::fftpack::complex_to_complex<double> cfft(10);
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

  scitbx::fftpack::real_to_complex<double> rfft(10);
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

  scitbx::fftpack::complex_to_complex_3d<double> cfft3d(2, 3, 5);
  scitbx::af::versa<std::complex<double>, scitbx::af::grid<3> >
  c3dmap(scitbx::af::grid<3>(cfft3d.N()));
  cfft3d.forward(c3dmap.ref());
  cfft3d.backward(c3dmap.ref());

  scitbx::fftpack::real_to_complex_3d<double> rfft3d(3, 4, 5);
  scitbx::af::versa<double, scitbx::af::grid<3> >
  r3dmap(scitbx::af::grid<3>(rfft3d.Mreal()));
  rfft3d.forward(r3dmap.ref());
  rfft3d.backward(r3dmap.ref());
#ifdef NEVER_DEFINED
#endif

  return 0;
}
