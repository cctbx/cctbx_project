#include <boost/python/class.hpp>

#include <scitbx/math/numeric_limits.h>

namespace scitbx { namespace math { namespace boost_python {

template<typename FloatType>
struct numeric_limits_wrapper
{
  typedef numeric_limits<FloatType> wt;

  // Workaround for a mysterious runtime failure on Windows with VC 8
  // when wrap() is executed on import of this Boost.Python extension
  static int radix()                 { return wt::radix;             }
  static FloatType min()             { return wt::min();             }
  static FloatType max()             { return wt::max();             }
  static int min_exponent()          { return wt::min_exponent;      }
  static int max_exponent()          { return wt::max_exponent;      }
  static int min_exponent10()        { return wt::min_exponent10;    }
  static int max_exponent10()        { return wt::max_exponent10;    }
  static int digits()                { return wt::digits;            }
  static int digits10()              { return wt::digits10;          }
  static FloatType epsilon()         { return wt::epsilon();         }
  static FloatType epsilon_x_radix() { return wt::epsilon_x_radix(); }
  static FloatType safe_min()        { return wt::safe_min();        }

  static void wrap(char const *name) {
    using namespace boost::python;
    typedef numeric_limits_wrapper<FloatType> wt;
    class_<wt>(name, no_init)
      .add_static_property("radix", &wt::radix)
      .add_static_property("min", &wt::min)
      .add_static_property("max", &wt::max)
      .add_static_property("min_exponent", &wt::min_exponent)
      .add_static_property("max_exponent", &wt::max_exponent)
      .add_static_property("min_exponent10", &wt::min_exponent10)
      .add_static_property("max_exponent10", &wt::max_exponent10)
      .add_static_property("digits", &wt::digits)
      .add_static_property("digits10", &wt::digits10)
      .add_static_property("epsilon", &wt::epsilon)
      .add_static_property("epsilon_x_radix", &wt::epsilon_x_radix)
      .add_static_property("safe_min", &wt::safe_min)
      ;
  }
};

void wrap_numeric_limits() {
  numeric_limits_wrapper<double>::wrap("double_numeric_limits");
}

}}} // scitbx::af::boost_python
