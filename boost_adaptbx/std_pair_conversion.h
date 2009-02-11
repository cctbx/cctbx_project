#ifndef BOOST_ADAPTBX_STD_PAIR_CONVERSION_H
#define BOOST_ADAPTBX_STD_PAIR_CONVERSION_H

#include <boost/python/tuple.hpp>
#include <boost/python/extract.hpp>
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

  template <typename T, typename U>
  struct from_tuple
  {
    from_tuple() {
      using namespace boost::python::converter;
      registry::push_back(&convertible,
                          &construct,
                          boost::python::type_id<std::pair<T,U> >()
                          #ifdef BOOST_PYTHON_SUPPORTS_PY_SIGNATURES
                          , get_pytype
                          #endif
      );
    }

    static const PyTypeObject *get_pytype() { return &PyTuple_Type; }

    static void *convertible(PyObject *o) {
      using namespace boost::python;
      if (!PyTuple_Check(o) || PyTuple_GET_SIZE(o) != 2) return 0;
      return o;
    }

    static void construct(
      PyObject *o,
      boost::python::converter::rvalue_from_python_stage1_data *data)
    {
      using boost::python::extract;
      using namespace boost::python::converter;
      PyObject *first  = PyTuple_GET_ITEM(o, 0),
               *second = PyTuple_GET_ITEM(o, 1);
      void *storage =
        ((rvalue_from_python_storage<std::pair<T,U> >*) data)->storage.bytes;
      new (storage) std::pair<T,U>(extract<T>(first), extract<U>(second));
      data->convertible = storage;
    }
  };

  template <typename T, typename U>
  struct to_and_from_tuple
  {
    to_and_from_tuple() {
      to_tuple<T,U>();
      from_tuple<T,U>();
    }
  };


}}

#endif
