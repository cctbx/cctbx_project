#include <boost/python/module.hpp>
#include <boost/python/def.hpp>

namespace {

  long
  factorial(int n)
  {
    long fact = 1;
    for(int i=2;i<=n;i++) {
      fact *= i;
    }
    return fact;
  }

  double
  power(double x, int n)
  {
    double pow = x;
    for(int i=1;i<n;i++) {
      pow *= x;
    }
    return pow;
  }

  double
  sin(double x, int n_terms)
  {
    double result = x;
    int sign = -1;
    int pow = 3;
    for(int i=1;i<n_terms;i++) {
      result += power(x,pow)/(sign*factorial(pow));
      sign *= -1;
      pow += 2;
    }
    return result;
  }

  double
  run_c_plus_plus(int n, int n_terms)
  {
    double result = 0;
    while (n--) {
      for(int i=0;i<180;i++) {
        result += sin(i * 3.14159265359/180, n_terms);
      }
    }
    return result;
  }

} // namespace anonymous

BOOST_PYTHON_MODULE(boost_python_hybrid_times_ext)
{
  using namespace boost::python;
  def("run_c_plus_plus", run_c_plus_plus);
}
