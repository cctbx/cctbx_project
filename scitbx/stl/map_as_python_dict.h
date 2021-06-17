#ifndef SCITBX_BOOST_PYTHON_STL_MAP_AS_PYTHON_DICT_H
#define SCITBX_BOOST_PYTHON_STL_MAP_AS_PYTHON_DICT_H

#include <boost/python/dict.hpp>
#include <boost/python/list.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/detail/api_placeholder.hpp>

#if PY_MAJOR_VERSION >= 3
#define IS_PY3K
#endif

namespace scitbx { namespace stl { namespace boost_python {

  template <typename MapType>
  boost::python::dict
  map_as_dict(MapType const& m)
  {
    boost::python::dict result;
    for(typename MapType::const_iterator
          pair=m.begin(); pair!=m.end(); pair++) {
      result[pair->first] = pair->second;
    }
    return result;
  }

  template <typename MapType>
  void
  update_map_from_dict(MapType& m, boost::python::dict d)
  {
    typedef typename MapType::key_type key_type;
    typedef typename MapType::mapped_type mapped_type;
    using namespace boost::python;
    #ifdef IS_PY3K
      PyObject *key_py3, *value_py3;
      Py_ssize_t pos = 0;
      PyObject* d_py3 = boost::python::object(d).ptr();
      while (PyDict_Next(d_py3, &pos, &key_py3, &value_py3)) {
        key_type key = extract<key_type>(key_py3)();
        mapped_type value = extract<mapped_type>(value_py3)();
        m[key]=value;
      }
    #else
      list keys = extract<list>(d.attr("keys")())();
      list values = extract<list>(d.attr("values")())();
      std::size_t len_keys = len(keys);
      for(std::size_t i=0;i<len_keys;i++) {
        key_type key = extract<key_type>(keys[i])();
        mapped_type value = extract<mapped_type>(values[i])();
        m[key] = value;
      }
    #endif
  }

}}} // namespace scitbx::stl::boost_python

#endif // SCITBX_BOOST_PYTHON_STL_MAP_AS_PYTHON_DICT_H
