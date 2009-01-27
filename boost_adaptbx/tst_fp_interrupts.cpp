#include <boost/python/module.hpp>
#include <boost/python/def.hpp>

double division(double x, double y) {
  return x/y;
}

BOOST_PYTHON_MODULE(tst_fp_interrupts_ext)
{
  using namespace boost::python;
  def("division", division);
}
