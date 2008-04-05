#ifndef IOTBX_PDB_HIERARCHY_BPL_H
#define IOTBX_PDB_HIERARCHY_BPL_H

#include <boost/python/str.hpp>
#include <boost/optional.hpp>

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

  template <typename ChildType, typename ParentType>
  struct get_parent
  {
    static
    boost::python::object
    wrapper(ChildType const& child)
    {
      boost::optional<ParentType> parent = child.parent();
      if (!parent) return boost::python::object();
      return boost::python::object(*parent);
    }
  };

#endif // IOTBX_PDB_HIERARCHY_BPL_H
