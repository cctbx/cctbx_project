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

#include <boost/smart_ptr.hpp>

#include <cctbx/fftbx/complex_to_complex_3d.h>
#include <cctbx/fftbx/real_to_complex_3d.h>

#include <fftw.h>
#include <rfftw.h>

using namespace cctbx;

namespace {

  typedef boost::shared_ptr<std::vector<double> > shared_vector;

  shared_vector init_cseq(const fftbx::triple& N)
  {
    shared_vector cseq(new std::vector<double>);
    for(int i=0;i<N.product(); i++) {
      cseq->push_back(double(37-i)/(N.max()+11));
      cseq->push_back(double(i-73)/(N.max()+13));
    }
    return cseq;
  }

  shared_vector tst_complex_to_complex(char dir, const fftbx::triple& N)
  {
    shared_vector cseq = init_cseq(N);
    fftbx::complex_to_complex_3d<std::vector<double> > fft(N);
    if (dir == 'f') {
      fft.forward(*cseq);
    }
    else {
      fft.backward(*cseq);
    }
    return cseq;
  }

  shared_vector tst_real_to_complex(char dir, const fftbx::triple& N)
  {
    fftbx::triple M = fftbx::Ncomplex_from_Nreal(N);
    shared_vector cseq = init_cseq(M);
    fftbx::real_to_complex_3d<std::vector<double> > fft(N);
    if (dir == 'f') {
      fft.forward(*cseq);
    }
    else {
      fft.backward(*cseq);
    }
    return cseq;
  }

  void show_fftw_complex(const fftw_complex *cseq, int N)
  {
    for(int i=0;i<N;i++) {
      std::cout << "re " << cseq[i].re << std::endl;
      std::cout << "im " << cseq[i].im << std::endl;
    }
  }

  shared_vector tst_fftw(char dir, const fftbx::triple& N)
  {
    shared_vector cseq = init_cseq(N);
    fftwnd_plan Plan;
    if (dir == 'f') {
      Plan = fftw3d_create_plan(
        N[0], N[1], N[2], FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
    }
    else {
      Plan = fftw3d_create_plan(
        N[0], N[1], N[2], FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
    }
    //show_fftw_complex((fftw_complex *) &(*cseq)[0], N.product());
    fftwnd_one(Plan, (fftw_complex *) &(*cseq)[0], 0);
    fftwnd_destroy_plan(Plan);
    return cseq;
  }

  shared_vector tst_rfftw(char dir, const fftbx::triple& N)
  {
    //std::cout << N[0] << " " << N[1] << " " << N[2] << std::endl;
    fftbx::triple M = fftbx::Ncomplex_from_Nreal(N);
    //std::cout << M[0] << " " << M[1] << " " << M[2] << std::endl;
    shared_vector cseq = init_cseq(M);
    //std::cout << cseq->size() << std::endl;
    rfftwnd_plan Plan;
    if (dir == 'f') {
      Plan = rfftw3d_create_plan(
        N[0], N[1], N[2], FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
      rfftwnd_one_real_to_complex(Plan, (fftw_real *) &(*cseq)[0], 0);
    }
    else {
      Plan = rfftw3d_create_plan(
        N[0], N[1], N[2], FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
      rfftwnd_one_complex_to_real(Plan, (fftw_complex *) &(*cseq)[0], 0);
    }
    rfftwnd_destroy_plan(Plan);
    return cseq;
  }

  void show_cseq(const shared_vector& cseq)
  {
    for(std::size_t i=0;i<cseq->size();i++) {
      std::cout << (*cseq)[i] << std::endl;
    }
  }

  std::ostream& operator<<(std::ostream& os, const fftbx::triple& t) {
    std::cout << t[0] << " " << t[1] << " " << t[2];
    return os;
  }

  void timing_complex_to_complex_3d(char dir,
                                    const fftbx::triple& N,
                                    std::size_t loop_iterations)
  {
    std::vector<double> cseq(2 * N.product());
    fftbx::complex_to_complex_3d<std::vector<double> > fft(N);
    if (dir == 'f') {
      std::cout << "timing_complex_to_complex_3d forward " << N << std::endl;
      for (std::size_t i=0;i<loop_iterations;i++) {
        fft.forward(cseq);
      }
    }
    else {
      std::cout << "timing_complex_to_complex_3d backward " << N << std::endl;
      for (std::size_t i=0;i<loop_iterations;i++) {
        fft.backward(cseq);
      }
    }
  }

  void timing_real_to_complex_3d(char dir,
                                 const fftbx::triple& N,
                                 std::size_t loop_iterations)
  {
    fftbx::triple M = fftbx::Ncomplex_from_Nreal(N);
    std::vector<double> cseq(2 * M.product());
    fftbx::real_to_complex_3d<std::vector<double> > fft(N);
    if (dir == 'f') {
      std::cout << "timing_real_to_complex_3d forward " << N << std::endl;
      for (std::size_t i=0;i<loop_iterations;i++) {
        fft.forward(cseq);
      }
    }
    else {
      std::cout << "timing_real_to_complex_3d backward " << N << std::endl;
      for (std::size_t i=0;i<loop_iterations;i++) {
        fft.backward(cseq);
      }
    }
  }

  void timing_fftw_3d(char dir,
                      const fftbx::triple& N,
                      std::size_t loop_iterations)
  {
    std::vector<double> cseq(2 * N.product());
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
                       const fftbx::triple& N,
                       std::size_t loop_iterations)
  {
    fftbx::triple M = fftbx::Ncomplex_from_Nreal(N);
    std::vector<double> cseq(2 * M.product());
    if (dir == 'f') {
      std::cout << "timing_rfftw_3d forward " << N << std::endl;
      rfftwnd_plan Plan = rfftw3d_create_plan(
        N[0], N[1], N[2], FFTW_FORWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
      for (std::size_t i=0;i<loop_iterations;i++) {
        rfftwnd_one_real_to_complex(Plan, (fftw_real *) &cseq[0], 0);
      }
      rfftwnd_destroy_plan(Plan);
    }
    else {
      std::cout << "timing_rfftw_3d backward " << N << std::endl;
      rfftwnd_plan Plan = rfftw3d_create_plan(
        N[0], N[1], N[2], FFTW_BACKWARD, FFTW_ESTIMATE | FFTW_IN_PLACE);
      for (std::size_t i=0;i<loop_iterations;i++) {
        rfftwnd_one_complex_to_real(Plan, (fftw_complex *) &cseq[0], 0);
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
  fftbx::triple N(atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));
  int iter = atoi(argv[6]);
  shared_vector cseq;
  if (iter < 0) {
    if (package == "fftbx") {
      if (type_and_dir[0] == 'c') {
        cseq = tst_complex_to_complex(type_and_dir[1], N);
      }
      else {
        cseq = tst_real_to_complex(type_and_dir[1], N);
      }
    }
    else {
      if (type_and_dir[0] == 'c') {
        cseq = tst_fftw(type_and_dir[1], N);
      }
      else {
        cseq = tst_rfftw(type_and_dir[1], N);
      }
    }
    show_cseq(cseq);
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
