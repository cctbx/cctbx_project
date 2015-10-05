#include <boost/python/module.hpp>
#include <boost/python/def.hpp>

#include <boost_adaptbx/tuple_conversion.h>
#include <boost/tuple/tuple.hpp>
#include <boost/thread.hpp>
#include <boost/ref.hpp>

namespace boost_adaptbx { namespace python { namespace {

const int N=100000;

struct pi_computation {
  double pi;
  void operator()() {
      pi = 0;
      double a=1, b=3;
      for(int i=1; i<2*N; ++i) {
        pi += 1/a - 1/b;
        a += 4;
        b += 4;
      }
      pi *= 4;
  }
};

struct e_computation {
  double e;
  void operator()() {
      e = 1;
      double a = 1;
      for(int i=1; i<N; ++i) {
        a /= i;
        e += a;
      }
  }
};

boost::tuple<double, double> test_boost_thread() {
  double e, pi;
  pi_computation compute_pi;
  e_computation compute_e;
  tuple_conversion::to_python<boost::tuple<double, double> >();
  boost::thread pi_worker(boost::ref(compute_pi));
  boost::thread e_worker(boost::ref(compute_e));
  pi_worker.join(); e_worker.join();
  return boost::make_tuple(compute_pi.pi, compute_e.e);
}

void wrap_all() {
  using namespace boost::python;
  def("test_boost_thread", test_boost_thread);
}

}}}

BOOST_PYTHON_MODULE(boost_adaptbx_boost_thread_test_ext)
{
  boost_adaptbx::python::wrap_all();
}
