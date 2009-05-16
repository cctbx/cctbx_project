#ifndef BOOST_ADAPTBX_EASY_OVERLOADS_H
#define BOOST_ADAPTBX_EASY_OVERLOADS_H

#include <boost/python/overloads.hpp>

/// Easier definition of overloaded functions for Boost.Python.
/** Synopsis:
    \code
    double foo(int x, long y=1, double z=2);

    ...

    BOOST_ADAPTBX_FUNCTION_OVERLOADS(foo_overloads,
                                     foo, 1, 3,
                                     boost::python::args("x", "y", "z"));

    ...

    foo_overloads::wrap("py_foo");
    \endcode

    Templates are easily handled as
    \code
    template <typename Tx, typename Ty, typename Tz>
    double foo(Tx x, Ty y=1, Tz z=2);

    ...

    template <typename Tx, typename Ty, typename Tz>
    BOOST_ADAPTBX_FUNCTION_OVERLOADS(foo_overloads,
                                     foo<Tx, Ty, Tz>, 1, 3,
                                     boost::python::args("x", "y", "z"));

    ...

    foo_overloads<int, long, double>::wrap("py_foo");
    \endcode
*/
#define BOOST_ADAPTBX_FUNCTION_OVERLOADS(name, func, minargs, maxargs, keywords)\
  struct name                                                                   \
  {                                                                             \
    BOOST_PYTHON_FUNCTION_OVERLOADS(overloads, func, minargs, maxargs);         \
    static void wrap(char const *func_name)                                     \
    {                                                                           \
      boost::python::def(func_name, func, overloads(keywords));                 \
    }                                                                           \
  }


#endif // GUARD
