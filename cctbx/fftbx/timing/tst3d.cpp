// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Nov 12: Created (R.W. Grosse-Kunstleve)
 */

#include <stdlib.h>
#include <iostream>
#include <string>
#include <cassert>

#include <boost/smart_ptr.hpp>

#include <cctbx/fftbx/complex_to_complex_3d.h>
#include <cctbx/fftbx/real_to_complex_3d.h>

#include <fftw.h>
#include <rfftw.h>

using namespace cctbx;

namespace {

  typedef boost::shared_ptr<std::vector<double> >
  shared_real_vector;
  typedef boost::shared_ptr<std::vector<std::complex<double> > >
  shared_complex_vector;

  typedef cctbx::array<std::size_t, 3> triple;

  shared_complex_vector init_cseq(const triple& N)
  {
    shared_complex_vector cseq(new shared_complex_vector::element_type);
    for(int i=0;i<cctbx::vector::product(N); i++) {
      cseq->push_back(shared_complex_vector::element_type::value_type(
        double(37-i)/(cctbx::vector::max(N)+11),
        double(i-73)/(cctbx::vector::max(N)+13)));
    }
    return cseq;
  }

  shared_real_vector init_rseq(const triple& N)
  {
    shared_real_vector rseq(new shared_real_vector::element_type);
    for(int i=0;i<cctbx::vector::product(N); i++) {
      rseq->push_back(double(37-i)/(cctbx::vector::max(N)+11));
    }
    return rseq;
  }

  shared_complex_vector tst_complex_to_complex(char dir, const triple& N)
  {
    fftbx::complex_to_complex_3d<double> fft(N);
    shared_complex_vector cseq = init_cseq(N);
    vecrefnd<
      shared_complex_vector::element_type::value_type,
      dimension<3> >
    cmap(cseq->begin(), N);
    if (dir == 'f') {
      fft.forward(cmap);
    }
    else {
      fft.backward(cmap);
    }
    return cseq;
  }

  shared_real_vector tst_real_to_complex(char dir, const triple& N)
  {
    fftbx::real_to_complex_3d<double> fft(N);
    shared_real_vector rseq = init_rseq(fft.Mreal());
    vecrefnd<
      shared_real_vector::element_type::value_type,
      dimension<3> >
    rmap(rseq->begin(), fft.Mreal());
    if (dir == 'f') {
      fft.forward(rmap);
    }
    else {
      fft.forward(rmap); // complex values have some symmetry
      fft.backward(rmap);
    }
    return rseq;
  }

  void show_fftw_complex(const fftw_complex *cseq, int N)
  {
    for(int i=0;i<N;i++) {
      std::cout << "re " << cseq[i].re << std::endl;
      std::cout << "im " << cseq[i].im << std::endl;
    }
  }

  shared_complex_vector tst_fftw(char dir, const triple& N)
  {
    shared_complex_vector cseq = init_cseq(N);
    fftwnd_plan Plan;
    if (dir == 'f') {
      Plan = fftw3d_create_plan(
        N[0], N[1], N[2], FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
    }
    else {
      Plan = fftw3d_create_plan(
        N[0], N[1], N[2], FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
    }
    fftwnd_one(Plan, (fftw_complex *) &(*cseq)[0], 0);
    fftwnd_destroy_plan(Plan);
    return cseq;
  }

  shared_real_vector tst_rfftw(char dir, const triple& N)
  {
    fftbx::real_to_complex_3d<double> fft(N);
    shared_real_vector rseq = init_rseq(fft.Mreal());
    rfftwnd_plan Plan;
    if (dir == 'f') {
      Plan = rfftw3d_create_plan(
        N[0], N[1], N[2], FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
      rfftwnd_one_real_to_complex(Plan, (fftw_real *) &(*rseq)[0], 0);
    }
    else {
      vecrefnd<
        shared_real_vector::element_type::value_type,
        dimension<3> >
      rmap(rseq->begin(), fft.Mreal());
      fft.forward(rmap); // complex values have some symmetry
      Plan = rfftw3d_create_plan(
        N[0], N[1], N[2], FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
      rfftwnd_one_complex_to_real(Plan, (fftw_complex *) &(*rseq)[0], 0);
    }
    rfftwnd_destroy_plan(Plan);
    return rseq;
  }

  void show_cseq(const shared_complex_vector& cseq)
  {
    for(std::size_t i=0;i<cseq->size();i++) {
      std::cout << (*cseq)[i].real() << std::endl;
      std::cout << (*cseq)[i].imag() << std::endl;
    }
  }

  void show_rseq(const shared_real_vector& rseq, const triple& N)
  {
    triple M = fftbx::Mreal_from_Nreal(N);
    assert(rseq->size() == cctbx::vector::product(M));
    c_index_1d<3> i1d;
    triple I;
    for(I[0]=0;I[0]<N[0];I[0]++)
    for(I[1]=0;I[1]<N[1];I[1]++)
    for(I[2]=0;I[2]<N[2];I[2]++) {
      std::cout << (*rseq)[i1d(M, I)] << std::endl;
    }
  }

  std::ostream& operator<<(std::ostream& os, const triple& t) {
    std::cout << t[0] << " " << t[1] << " " << t[2];
    return os;
  }

  void timing_complex_to_complex_3d(char dir,
                                    const triple& N,
                                    std::size_t loop_iterations)
  {
    shared_complex_vector cseq(new shared_complex_vector::element_type);
    cseq->resize(cctbx::vector::product(N));
    vecrefnd<
      shared_complex_vector::element_type::value_type,
      dimension<3> >
    cmap(cseq->begin(), N);
    fftbx::complex_to_complex_3d<double> fft(N);
    if (dir == 'f') {
      std::cout << "timing_complex_to_complex_3d forward " << N << std::endl;
      for (std::size_t i=0;i<loop_iterations;i++) {
        fft.forward(cmap);
      }
    }
    else {
      std::cout << "timing_complex_to_complex_3d backward " << N << std::endl;
      for (std::size_t i=0;i<loop_iterations;i++) {
        fft.backward(cmap);
      }
    }
  }

  void timing_real_to_complex_3d(char dir,
                                 const triple& N,
                                 std::size_t loop_iterations)
  {
    fftbx::real_to_complex_3d<double> fft(N);
    shared_real_vector rseq(new shared_real_vector::element_type);
    rseq->resize(cctbx::vector::product(fft.Mreal()));
    vecrefnd<
      shared_real_vector::element_type::value_type,
      dimension<3> >
    rmap(rseq->begin(), fft.Mreal());
    if (dir == 'f') {
      std::cout << "timing_real_to_complex_3d forward " << N << std::endl;
      for (std::size_t i=0;i<loop_iterations;i++) {
        fft.forward(rmap);
      }
    }
    else {
      std::cout << "timing_real_to_complex_3d backward " << N << std::endl;
      for (std::size_t i=0;i<loop_iterations;i++) {
        fft.backward(rmap);
      }
    }
  }

  void timing_fftw_3d(char dir,
                      const triple& N,
                      std::size_t loop_iterations)
  {
    std::vector<std::complex<double> > cseq(cctbx::vector::product(N));
    if (dir == 'f') {
      std::cout << "timing_fftw_3d forward " << N << std::endl;
      fftwnd_plan Plan = fftw3d_create_plan(
        N[0], N[1], N[2], FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
      for (std::size_t i=0;i<loop_iterations;i++) {
        fftwnd_one(Plan, (fftw_complex *) &cseq[0], 0);
      }
      fftwnd_destroy_plan(Plan);
    }
    else {
      std::cout << "timing_fftw_3d backward " << N << std::endl;
      fftwnd_plan Plan = fftw3d_create_plan(
        N[0], N[1], N[2], FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
      for (std::size_t i=0;i<loop_iterations;i++) {
        fftwnd_one(Plan, (fftw_complex *) &cseq[0], 0);
      }
      fftwnd_destroy_plan(Plan);
    }
  }

  void timing_rfftw_3d(char dir,
                       const triple& N,
                       std::size_t loop_iterations)
  {
    triple M = fftbx::Mreal_from_Nreal(N);
    std::vector<double> rseq(cctbx::vector::product(M));
    if (dir == 'f') {
      std::cout << "timing_rfftw_3d forward " << N << std::endl;
      rfftwnd_plan Plan = rfftw3d_create_plan(
        N[0], N[1], N[2], FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
      for (std::size_t i=0;i<loop_iterations;i++) {
        rfftwnd_one_real_to_complex(Plan, (fftw_real *) &rseq[0], 0);
      }
      rfftwnd_destroy_plan(Plan);
    }
    else {
      std::cout << "timing_rfftw_3d backward " << N << std::endl;
      rfftwnd_plan Plan = rfftw3d_create_plan(
        N[0], N[1], N[2], FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
      for (std::size_t i=0;i<loop_iterations;i++) {
        rfftwnd_one_complex_to_real(Plan, (fftw_complex *) &rseq[0], 0);
      }
      rfftwnd_destroy_plan(Plan);
    }
  }

}

void usage() {
  std::cerr
    << "usage: tst3d fftbx|fftw cf|cb|rf|rb Nx Ny Nz iter"
    << std::endl;
  exit(1);
}

int main(int argc, const char* argv[])
{
  if (argc != 7) usage();
  std::string package = std::string(argv[1]);
  if (package != "fftbx" && package != "fftw") usage();
  std::string type_and_dir = std::string(argv[2]);
  if (type_and_dir.size() != 2) usage();
  if (type_and_dir[0] != 'c' && type_and_dir[0] != 'r') usage();
  if (type_and_dir[1] != 'f' && type_and_dir[1] != 'b') usage();
  triple N(atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));
  int iter = atoi(argv[6]);
  if (iter < 0) {
    if (package == "fftbx") {
      if (type_and_dir[0] == 'c') {
        shared_complex_vector cseq = tst_complex_to_complex(type_and_dir[1], N);
        show_cseq(cseq);
      }
      else {
        shared_real_vector rseq = tst_real_to_complex(type_and_dir[1], N);
        show_rseq(rseq, N);
      }
    }
    else {
      if (type_and_dir[0] == 'c') {
        shared_complex_vector cseq = tst_fftw(type_and_dir[1], N);
        show_cseq(cseq);
      }
      else {
        shared_real_vector rseq = tst_rfftw(type_and_dir[1], N);
        show_rseq(rseq, N);
      }
    }
  }
  else {
    if (package == "fftbx") {
      if (type_and_dir[0] == 'c') {
        timing_complex_to_complex_3d(type_and_dir[1], N, iter);
      }
      else {
        timing_real_to_complex_3d(type_and_dir[1], N, iter);
      }
    }
    else {
      if (type_and_dir[0] == 'c') {
        timing_fftw_3d(type_and_dir[1], N, iter);
      }
      else {
        timing_rfftw_3d(type_and_dir[1], N, iter);
      }
    }
  }
  return 0;
}
