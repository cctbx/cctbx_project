#include <vector>
#include <boost/numeric/ublas/config.h>
#include <boost/numeric/ublas/vector.h>

int main()
{
  numerics::vector<float, std::vector<float> > ua(1000);
  numerics::vector<float, std::vector<float> > ub(1000);
  numerics::vector<float, std::vector<float> > uc(1000);
  uc = ua + ub;
  return 0;
}

/*
#define ACTIVE_UBLAS
//#define ACTIVE_SHARED
//#define MUL_ADD

#if defined(ACTIVE_SHARED)
#include <cctbx/array_family/shared_algebra.h>
#include <cctbx/array_family/tiny.h> // XXX
#include <cctbx/array_family/ref_algebra.h>
#include <cctbx/array_family/simple_io.h>
#endif

#if defined(ACTIVE_UBLAS)
#include <boost/numeric/ublas/config.h>
#include <boost/numeric/ublas/vector.h>
//#include <cctbx/array_family/shared_plain.h>
#include <vector>
#endif

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <boost/timer.hpp>

int main(int argc, char* argv[])
{
  int power = 5;
  if (argc > 1) {
    power = std::atoi(argv[1]);
  }
  std::size_t max_n = std::pow(10., power) + 0.5;
  std::cout << "max_n = " << max_n << std::endl;
  boost::timer t;
  for (std::size_t n = 1; n <= max_n; n *= 10) {
    std::size_t n_iter = 3 * max_n / n;
    if (n_iter > 1000000) continue;
    std::cout << "n_iter = " << n_iter << std::endl;
#if defined(ACTIVE_SHARED)
    cctbx::af::shared<double> sa(n, 1.);
    cctbx::af::shared<double> sb(n, 2.);
    cctbx::af::shared<double> sc(n, 3.);
    cctbx::af::shared<double> sr;
#endif
#if defined(ACTIVE_UBLAS)
    typedef std::vector<double> base_array;
    numerics::vector<double, base_array> ua(n);
    numerics::vector<double, base_array> ub(n);
    numerics::vector<double, base_array> uc(n);
    numerics::vector<double, base_array> ur(n);
    t.restart();
    for (std::size_t i = 0; i < n_iter; i++) {
#if defined(MUL_ADD)
      ur = ua * ub + uc;
#else
      ur = ua * ub;
#endif
    }
    std::cout << "ublas_shared(" << n << ") " << t.elapsed() << std::endl;
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
#if defined(XXX) && defined(ACTIVE_UBLAS) && defined(ACTIVE_SHARED)
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
*/
