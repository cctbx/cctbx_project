#define ACTIVE_BLITZ
#define ACTIVE_SHARED
//#define MUL_ADD

#include <cstdlib>
#include <iostream>
#include <boost/timer.hpp>
#include <cctbx/array_family/shared_algebra.h>
#include <cctbx/array_family/tiny.h> // XXX
#include <cctbx/array_family/ref_algebra.h>
#include <cctbx/array_family/simple_io.h>

#if defined(ACTIVE_BLITZ)
#include <blitz/array.h>
#endif

int main(int argc, char* argv[])
{
  std::size_t power = 5;
  if (argc > 1) {
    power = atoi(argv[1]);
  }
  std::size_t max_n = std::pow(10, power);
  std::cout << "max_n = " << max_n << std::endl;
  boost::timer t;
  for (std::size_t n = 1; n <= max_n; n *= 10) {
    std::size_t n_iter = 3 * max_n / n;
    if (n_iter > 1000000) continue;
    std::cout << "n_iter = " << n_iter << std::endl;
    cctbx::af::shared<double> sa(n, 1.);
    cctbx::af::shared<double> sb(n, 2.);
    cctbx::af::shared<double> sc(n, 3.);
    cctbx::af::shared<double> sr;
#if defined(ACTIVE_BLITZ)
    blitz::Array<double, 1> ba(sa.begin(), n, blitz::duplicateData);
    blitz::Array<double, 1> bb(sb.begin(), n, blitz::duplicateData);
    blitz::Array<double, 1> bc(sc.begin(), n, blitz::duplicateData);
    blitz::Array<double, 1> br(n);
    t.restart();
    for (std::size_t i = 0; i < n_iter; i++) {
#if defined(MUL_ADD)
      br = ba * bb + bc;
#else
      br = ba * bb;
#endif
    }
    std::cout << "blitz::Array(" << n << ") " << t.elapsed() << std::endl;
#endif
#if defined(ACTIVE_SHARED)
    t.restart();
    for (std::size_t i = 0; i < n_iter; i++) {
#if defined(MUL_ADD)
      sr = sa * sb + sc;
#else
      sr = sa * sb;
#endif
    }
    std::cout << "  af::shared(" << n << ") " << t.elapsed() << std::endl;
#endif
#if defined(ACTIVE_BLITZ) && defined(ACTIVE_SHARED)
    cctbx::af::ref<double> br_ref(&*br.begin(), br.size());
    if (br_ref.size() != sr.size()) {
      std::cout << "ERROR size" << std::endl;
    }
    if (br_ref != sr.ref()) {
      std::cout << "ERROR data" << std::endl;
    }
    //std::cout << bc_ref << std::endl;
    //std::cout << sc.ref() << std::endl;
#endif
  }

  return 0;
}
