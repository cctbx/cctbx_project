// $Id$
/* Copyright (c) 2001 The Regents of the University of California through
   E.O. Lawrence Berkeley National Laboratory, subject to approval by the
   U.S. Department of Energy. See files COPYRIGHT.txt and
   cctbx/LICENSE.txt for further details.

   Revision history:
     2001 Nov 08: Created (R.W. Grosse-Kunstleve)
 */

#include <stdlib.h>
#include <iostream>
#include <string>

#include <cctbx/fftbx/complex_to_complex.h>
#include <cctbx/fftbx/real_to_complex.h>

using namespace cctbx;

namespace {

  void
  timing_complex_to_complex(std::size_t transform_size,
                            std::size_t loop_iterations)
  {
    using namespace fftbx;
    std::vector<double> cseq(2 * transform_size);
    complex_to_complex<double> fft(transform_size);
    for (std::size_t i=0;i<loop_iterations;i++) {
      fft.forward(cseq.begin());
      fft.backward(cseq.begin());
    }
  }

  void
  timing_real_to_complex(std::size_t transform_size,
                         std::size_t loop_iterations)
  {
    using namespace fftbx;
    std::vector<double> rseq(transform_size);
    real_to_complex<double> fft(transform_size);
    for (std::size_t i=0;i<loop_iterations;i++) {
      fft.forward(rseq.begin());
      fft.backward(rseq.begin());
    }
  }

}

int main(int argc, const char* argv[])
{
  const char* usage = "need four arguments: cc|rc N iter factor";

  if (argc != 5) {
    std::cerr << usage << std::endl;
    return 1;
  }
  std::string fft_type;
  if (std::string(argv[1]) == "cc") {
    fft_type = "complex-to-complex";
  }
  else if (std::string(argv[1]) == "rc") {
    fft_type = "real-to-complex";
  }
  else {
    std::cerr << usage << std::endl;
    return 1;
  }
  std::size_t transform_size = atoi(argv[2]);
  std::size_t loop_iterations = atoi(argv[3]);
  std::size_t factor = atoi(argv[4]);
  loop_iterations *= factor;
  std::cout << "fftbx " << fft_type << std::endl;
  std::cout << "Transform size: " << transform_size << std::endl;
  std::cout << "Loop iterations: " << loop_iterations << std::endl;
  if (fft_type == "complex-to-complex") {
    timing_complex_to_complex(transform_size, loop_iterations);
  }
  else {
    timing_real_to_complex(transform_size, loop_iterations);
  }
  return 0;
}
