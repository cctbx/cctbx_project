#ifndef BOOST_ADAPTBX_STD_PAIR_CONVERSION_H
#define BOOST_ADAPTBX_STD_PAIR_CONVERSION_H

#include <boost/python/tuple.hpp>
#include <boost/python/to_python_converter.hpp>

namespace boost_adaptbx { namespace std_pair_conversions {

  namespace detail {
    template <typename T, typename U>
    struct to_tuple
    {
      static PyObject* convert(std::pair<T,U> const& p) {
        using namespace boost::python;
        return incref(boost::python::make_tuple(p.first, p.second).ptr());
      }

      static PyTypeObject const *get_pytype() { return &PyTuple_Type; }
    };
  }

  template <typename T, typename U>
  struct to_tuple
  {
    to_tuple() {
      using namespace boost::python;
      to_python_converter<std::pair<T,U>, detail::to_tuple<T,U>
#ifdef BOOST_PYTHON_SUPPORTS_PY_SIGNATURES
                                    , true
#endif
      >();
    }
  };

}}

#endif
