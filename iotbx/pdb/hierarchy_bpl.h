#ifndef IOTBX_PDB_HIERARCHY_BPL_H
#define IOTBX_PDB_HIERARCHY_BPL_H

#include <boost/python/str.hpp>
#include <boost/optional.hpp>
#include <iotbx/pdb/hybrid_36_c.h>

#define IOTBX_PDB_HIERARCHY_DATA_WRAPPERS_SMALL_STR_GET(attr) \
    static \
    boost::python::str \
    get_##attr(w_t const& self) \
    { \
      return boost::python::str(self.data->attr.elems); \
    }

#define IOTBX_PDB_HIERARCHY_DATA_WRAPPERS_SMALL_STR_SET(attr) \
    static \
    void \
    set_##attr(w_t& self, const char* value) \
    { \
      self.data->attr.replace_with(value); \
    }

#define IOTBX_PDB_HIERARCHY_DATA_WRAPPERS_SMALL_STR_GET_SET(attr) \
  IOTBX_PDB_HIERARCHY_DATA_WRAPPERS_SMALL_STR_GET(attr) \
  IOTBX_PDB_HIERARCHY_DATA_WRAPPERS_SMALL_STR_SET(attr)

#define IOTBX_PDB_HIERARCHY_WRAPPERS_SET_HY36( \
          attr, data_attr, width, value_min, value_max) \
    static void \
    set_##attr(w_t& self, boost::python::object const& value) \
    { \
      PyObject* ptr = value.ptr(); \
      if (ptr == Py_None) { \
        self.data_attr.replace_with(0); \
        return; \
      } \
      if (PyString_Check(ptr)) { \
        self.data_attr.replace_with(PyString_AS_STRING(ptr)); \
        return; \
      } \
      if (PyInt_Check(ptr)) { \
        long v = PyInt_AS_LONG(ptr); \
        if (v < value_min) { \
          PyErr_SetString( \
            PyExc_ValueError, "value is less than " #value_min); \
          boost::python::throw_error_already_set(); \
        } \
        if (v > value_max) { \
          PyErr_SetString( \
            PyExc_ValueError, "value is greater than " #value_max); \
          boost::python::throw_error_already_set(); \
        } \
        const char* errmsg = hy36encode( \
          width, static_cast<int>(v), self.data_attr.elems); \
        if (errmsg != 0) { \
          PyErr_SetString(PyExc_ValueError, errmsg); \
          boost::python::throw_error_already_set(); \
        } \
        return; \
      } \
      PyErr_SetString(PyExc_TypeError, "value must be a Python str or int."); \
      boost::python::throw_error_already_set(); \
    }

  template <typename ChildType, typename ParentType>
  struct get_parent
  {
    static
    boost::python::object
    wrapper(ChildType const& child, bool optional)
    {
      boost::optional<ParentType> parent = child.parent(optional);
      if (!parent) return boost::python::object();
      return boost::python::object(*parent);
    }
  };

#endif // IOTBX_PDB_HIERARCHY_BPL_H
