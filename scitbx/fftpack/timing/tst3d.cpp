/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Aug: Ported from cctbx (rwgk)
     2001 Nov: Created (R.W. Grosse-Kunstleve)
 */

#include <stdlib.h>
#include <iostream>
#include <string>
#include <cassert>

#include <scitbx/fftpack/complex_to_complex_3d.h>
#include <scitbx/fftpack/real_to_complex_3d.h>

#include <fftw.h>
#include <rfftw.h>

using namespace scitbx;

namespace {

  typedef af::shared<double> shared_real_array;
  typedef af::shared<std::complex<double> > shared_complex_array;

  shared_complex_array init_cseq(const af::int3& n)
  {
    shared_complex_array cseq(af::product(n), af::reserve_flag());
    for(int i=0;i<cseq.capacity(); i++) {
      cseq.push_back(shared_complex_array::value_type(
        double(37-i)/(af::max(n)+11),
        double(i-73)/(af::max(n)+13)));
    }
    return cseq;
  }

  shared_real_array init_rseq(const af::int3& n)
  {
    shared_real_array rseq(af::product(n), af::reserve_flag());
    for(int i=0;i<rseq.capacity(); i++) {
      rseq.push_back(double(37-i)/(af::max(n)+11));
    }
    return rseq;
  }

  shared_complex_array tst_complex_to_complex(char dir, const af::int3& n)
  {
    fftpack::complex_to_complex_3d<double> fft(n);
    shared_complex_array cseq = init_cseq(n);
    af::ref<shared_complex_array::value_type, af::grid<3> >
    cmap(cseq.begin(), af::grid<3>(n));
    if (dir == 'f') {
      fft.forward(cmap);
    }
    else {
      fft.backward(cmap);
    }
    return cseq;
  }

  shared_real_array tst_real_to_complex(char dir, const af::int3& n)
  {
    fftpack::real_to_complex_3d<double> fft(n);
    shared_real_array rseq = init_rseq(fft.m_real());
    af::ref<shared_real_array::value_type, af::grid<3> >
    rmap(rseq.begin(), af::grid<3>(fft.m_real()));
    if (dir == 'f') {
      fft.forward(rmap);
    }
    else {
      fft.forward(rmap); // complex values have some symmetry
      fft.backward(rmap);
    }
    return rseq;
  }

  void show_fftw_complex(const fftw_complex *cseq, int n)
  {
    for(int i=0;i<n;i++) {
      std::cout << "re " << cseq[i].re << std::endl;
      std::cout << "im " << cseq[i].im << std::endl;
    }
  }

  shared_complex_array tst_fftw(char dir, const af::int3& n)
  {
    shared_complex_array cseq = init_cseq(n);
    fftwnd_plan Plan;
    if (dir == 'f') {
      Plan = fftw3d_create_plan(
        n[0], n[1], n[2], FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
    }
    else {
      Plan = fftw3d_create_plan(
        n[0], n[1], n[2], FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
    }
    fftwnd_one(Plan, (fftw_complex *) cseq.begin(), 0);
    fftwnd_destroy_plan(Plan);
    return cseq;
  }

  shared_real_array tst_rfftw(char dir, const af::int3& n)
  {
    fftpack::real_to_complex_3d<double> fft(n);
    shared_real_array rseq = init_rseq(fft.m_real());
    rfftwnd_plan Plan;
    if (dir == 'f') {
      Plan = rfftw3d_create_plan(
        n[0], n[1], n[2], FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
      rfftwnd_one_real_to_complex(Plan, (fftw_real *) rseq.begin(), 0);
    }
    else {
      af::ref<shared_real_array::value_type, af::grid<3> >
      rmap(rseq.begin(), af::grid<3>(fft.m_real()));
      fft.forward(rmap); // complex values have some symmetry
      Plan = rfftw3d_create_plan(
        n[0], n[1], n[2], FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
      rfftwnd_one_complex_to_real(Plan, (fftw_complex *) rseq.begin(), 0);
    }
    rfftwnd_destroy_plan(Plan);
    return rseq;
  }

  void show_cseq(const shared_complex_array& cseq)
  {
    for(std::size_t i=0;i<cseq.size();i++) {
      std::cout << cseq[i].real() << std::endl;
      std::cout << cseq[i].imag() << std::endl;
    }
  }

  void show_rseq(const shared_real_array& rseq, const af::int3& n)
  {
    af::int3 m = fftpack::m_real_from_n_real(n);
    assert(rseq.size() == af::product(m));
    af::c_index_1d_calculator<3> i1d;
    af::int3 I;
    for(I[0]=0;I[0]<n[0];I[0]++)
    for(I[1]=0;I[1]<n[1];I[1]++)
    for(I[2]=0;I[2]<n[2];I[2]++) {
      std::cout << rseq[i1d.get(m, I)] << std::endl;
    }
  }

  std::ostream& operator<<(std::ostream& os, const af::int3& t) {
    std::cout << t[0] << " " << t[1] << " " << t[2];
    return os;
  }

  void timing_complex_to_complex_3d(char dir,
                                    const af::int3& n,
                                    std::size_t loop_iterations)
  {
    shared_complex_array cseq(af::product(n));
    af::ref<shared_complex_array::value_type, af::grid<3> >
    cmap(cseq.begin(), af::grid<3>(n));
    fftpack::complex_to_complex_3d<double> fft(n);
    if (dir == 'f') {
      std::cout << "timing_complex_to_complex_3d forward " << n << std::endl;
      for (std::size_t i=0;i<loop_iterations;i++) {
        fft.forward(cmap);
      }
    }
    else {
      std::cout << "timing_complex_to_complex_3d backward " << n << std::endl;
      for (std::size_t i=0;i<loop_iterations;i++) {
        fft.backward(cmap);
      }
    }
  }

  void timing_real_to_complex_3d(char dir,
                                 const af::int3& n,
                                 std::size_t loop_iterations)
  {
    fftpack::real_to_complex_3d<double> fft(n);
    shared_real_array rseq(af::product(fft.m_real()));
    af::ref<shared_real_array::value_type, af::grid<3> >
    rmap(rseq.begin(), af::grid<3>(fft.m_real()));
    if (dir == 'f') {
      std::cout << "timing_real_to_complex_3d forward " << n << std::endl;
      for (std::size_t i=0;i<loop_iterations;i++) {
        fft.forward(rmap);
      }
    }
    else {
      std::cout << "timing_real_to_complex_3d backward " << n << std::endl;
      for (std::size_t i=0;i<loop_iterations;i++) {
        fft.backward(rmap);
      }
    }
  }

  void timing_fftw_3d(char dir,
                      const af::int3& n,
                      std::size_t loop_iterations)
  {
    shared_complex_array cseq(af::product(n));
    if (dir == 'f') {
      std::cout << "timing_fftw_3d forward " << n << std::endl;
      fftwnd_plan Plan = fftw3d_create_plan(
        n[0], n[1], n[2], FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
      for (std::size_t i=0;i<loop_iterations;i++) {
        fftwnd_one(Plan, (fftw_complex *) cseq.begin(), 0);
      }
      fftwnd_destroy_plan(Plan);
    }
    else {
      std::cout << "timing_fftw_3d backward " << n << std::endl;
      fftwnd_plan Plan = fftw3d_create_plan(
        n[0], n[1], n[2], FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
      for (std::size_t i=0;i<loop_iterations;i++) {
        fftwnd_one(Plan, (fftw_complex *) cseq.begin(), 0);
      }
      fftwnd_destroy_plan(Plan);
    }
  }

  void timing_rfftw_3d(char dir,
                       const af::int3& n,
                       std::size_t loop_iterations)
  {
    af::int3 m = fftpack::m_real_from_n_real(n);
    shared_real_array rseq(af::product(m));
    if (dir == 'f') {
      std::cout << "timing_rfftw_3d forward " << n << std::endl;
      rfftwnd_plan Plan = rfftw3d_create_plan(
        n[0], n[1], n[2], FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
      for (std::size_t i=0;i<loop_iterations;i++) {
        rfftwnd_one_real_to_complex(Plan, (fftw_real *) rseq.begin(), 0);
      }
      rfftwnd_destroy_plan(Plan);
    }
    else {
      std::cout << "timing_rfftw_3d backward " << n << std::endl;
      rfftwnd_plan Plan = rfftw3d_create_plan(
        n[0], n[1], n[2], FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
      for (std::size_t i=0;i<loop_iterations;i++) {
        rfftwnd_one_complex_to_real(Plan, (fftw_complex *) rseq.begin(), 0);
      }
      rfftwnd_destroy_plan(Plan);
    }
  }

}

void usage() {
  std::cerr
    << "usage: tst3d fftpack|fftw cf|cb|rf|rb nx ny nz iter"
    << std::endl;
  exit(1);
}

int main(int argc, const char* argv[])
{
  if (argc != 7) usage();
  std::string package = std::string(argv[1]);
  if (package != "fftpack" && package != "fftw") usage();
  std::string type_and_dir = std::string(argv[2]);
  if (type_and_dir.size() != 2) usage();
  if (type_and_dir[0] != 'c' && type_and_dir[0] != 'r') usage();
  if (type_and_dir[1] != 'f' && type_and_dir[1] != 'b') usage();
  af::int3 n(atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));
  int iter = atoi(argv[6]);
  if (iter < 0) {
    if (package == "fftpack") {
      if (type_and_dir[0] == 'c') {
        shared_complex_array cseq = tst_complex_to_complex(type_and_dir[1],n);
        show_cseq(cseq);
      }
      else {
        shared_real_array rseq = tst_real_to_complex(type_and_dir[1], n);
        show_rseq(rseq, n);
      }
    }
    else {
      if (type_and_dir[0] == 'c') {
        shared_complex_array cseq = tst_fftw(type_and_dir[1], n);
        show_cseq(cseq);
      }
      else {
        shared_real_array rseq = tst_rfftw(type_and_dir[1], n);
        show_rseq(rseq, n);
      }
    }
  }
  else {
    if (package == "fftpack") {
      if (type_and_dir[0] == 'c') {
        timing_complex_to_complex_3d(type_and_dir[1], n, iter);
      }
      else {
        timing_real_to_complex_3d(type_and_dir[1], n, iter);
      }
    }
    else {
      if (type_and_dir[0] == 'c') {
        timing_fftw_3d(type_and_dir[1], n, iter);
      }
      else {
        timing_rfftw_3d(type_and_dir[1], n, iter);
      }
    }
  }
  return 0;
}
