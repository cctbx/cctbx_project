#pragma once
#include <boost/python/tuple.hpp>

namespace scitbx { namespace boost_python {
  namespace bp = boost::python;

  // used to be pair_as_tuple.hpp
  template<class Pair> struct PairToTupleConverter {
    static PyObject* convert(Pair const& pair) {
      return boost::python::incref(
        boost::python::make_tuple(pair.first, pair.second).ptr()
      );
    }
  };

  //https://stackoverflow.com/questions/16497889/how-to-expose-stdpair-to-python-using-boostpython
  template<typename T1, typename T2>
  struct PairToPythonConverter {
    static PyObject* convert(const std::pair<T1, T2>& pair) {
      return bp::incref(bp::make_tuple(pair.first, pair.second).ptr());
    }
  };

  template<typename T1, typename T2>
  struct PythonToPairConverter {
    PythonToPairConverter() {
      bp::converter::registry::push_back(&convertible, &construct,
        bp::type_id<std::pair<T1, T2> >());
    }
    static void* convertible(PyObject* obj) {
      if (!PyTuple_CheckExact(obj)) {
        return 0;
      }
      if (PyTuple_Size(obj) != 2) {
        return 0;
      }
      return obj;
    }
    static void construct(PyObject* obj,
      bp::converter::rvalue_from_python_stage1_data* data)
    {
      bp::tuple tuple(bp::borrowed(obj));
      void* storage = ((bp::converter::rvalue_from_python_storage<std::pair<T1, T2> >*) data)->storage.bytes;
      new (storage) std::pair<T1, T2>(bp::extract<T1>(tuple[0]), bp::extract<T2>(tuple[1]));
      data->convertible = storage;
    }
  };

  template <typename T1, typename T2>
  struct RegisterPyPair {
    static void to_py() {
      bp::to_python_converter<std::pair<T1, T2>, PairToPythonConverter<T1, T2> >();
    }
    static void from_py() {
      PythonToPairConverter<T1, T2>();
    }
    RegisterPyPair() {
      to_py();
      from_py();
    }
  };

}} // namespace scitbx::boost_python
