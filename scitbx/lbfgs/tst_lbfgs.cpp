#include <iostream>
#include <stdlib.h>
#include <scitbx/lbfgs.h>

int main(int argc, const char* argv[])
{
  int passes = 1;
  if (argc > 1) {
    passes = atoi(argv[1]);
  }
  int n = 100;
  std::vector<double> x(n); double* xb = &(*(x.begin()));
  std::vector<double> g(n); double* gb = &(*(g.begin()));
  for (int ipass = 0; ipass < passes; ipass++) {
    for(int j=0;j<n;j+=2) {
      x[j] = -1.2;
      x[j+1] = 1.;
    }
    scitbx::lbfgs::minimizer<double, int> minimizer(n);
    scitbx::lbfgs::traditional_convergence_test<double, int> is_converged(n);
    if (passes == 1) {
      std::cout << "n: " << minimizer.n() << std::endl;
      std::cout << "m: " << minimizer.m() << std::endl;
      std::cout << "xtol: " << minimizer.xtol() << std::endl;
    }
    for(;;) {
      double f = 0.;
      for(int j=0;j<n;j+=2) {
        double t1 = 1.e0 - x[j];
        double t2 = 1.e1 * (x[j+1] - x[j] * x[j]);
        g[j+1] = 2.e1 * t2;
        g[j] = -2.e0 * (x[j] * g[j+1] + t1);
        f = f + t1 * t1 + t2 * t2;
      }
      if (minimizer.run(xb, f, gb)) continue;
      if (passes == 1) {
        std::cout << "f: " << f;
        std::cout << " gnorm: " << minimizer.euclidean_norm(&(*(g.begin())));
        std::cout << std::endl;
        std::cout << " " << minimizer.iter();
        std::cout << " " << minimizer.nfun();
        std::cout << " " << minimizer.stp();
        std::cout << std::endl;
      }
      if (is_converged(xb, gb)) break;
      if (minimizer.nfun() > 2000) break;
      minimizer.run(xb, f, gb);
    }
  }
  return 0;
}
