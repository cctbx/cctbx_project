#include <iostream>
#include <stdlib.h>
#include <cctbx/lbfgs.h>

int main(int argc, const char* argv[])
{
  int passes = 1;
  if (argc > 1) {
    passes = atoi(argv[1]);
  }
  int n = 100;
  std::vector<double> x(n);
  std::vector<double> g(n);
  for (int ipass = 0; ipass < passes; ipass++) {
    for(int j=0;j<n;j+=2) {
      x[j] = -1.2;
      x[j+1] = 1.;
    }
    cctbx::lbfgs::minimizer<double> minimizer(n);
    // We allow at most 2000 evaluations of f and g
    while (!minimizer.is_converged() && minimizer.nfun() < 2000) {
      double f = 0.;
      for(int j=0;j<n;j+=2) {
        double t1 = 1.e0 - x[j];
        double t2 = 1.e1 * (x[j+1] - x[j] * x[j]);
        g[j+1] = 2.e1 * t2;
        g[j] = -2.e0 * (x[j] * g[j+1] + t1);
        f = f + t1 * t1 + t2 * t2;
      }
      if (passes == 1) {
        if (minimizer.nfun() == 0) {
          std::cout << "n: " << minimizer.n() << std::endl;
          std::cout << "m: " << minimizer.m() << std::endl;
          std::cout << "eps: " << minimizer.eps() << std::endl;
          std::cout << "xtol: " << minimizer.xtol() << std::endl;
        }
        std::cout << "f: " << f;
        std::cout << " gnorm: " << minimizer.gnorm(&(*(g.begin())));
        std::cout << std::endl;
      }
      minimizer.run(&(*(x.begin())), f, &(*(g.begin())));
      if (passes == 1) {
        std::cout << " " << minimizer.iter();
        std::cout << " " << minimizer.nfun();
        std::cout << " " << minimizer.stp();
        std::cout << std::endl;
      }
    }
  }
  return 0;
}
