#ifndef SCITBX_BOOST_PYTHON_STL_MAP_AS_DICT_H
#define SCITBX_BOOST_PYTHON_STL_MAP_AS_DICT_H

#include <boost/python/dict.hpp>

namespace scitbx { namespace boost_python {

  template <typename MapType>
  boost::python::dict
  stl_map_as_dict(MapType const& m)
  {
    boost::python::dict result;
    typedef typename MapType::const_iterator mci;
    for(mci i=m.begin(); i!= m.end(); i++) {
      result[i->first] = i->second;
    }
    return result;
  }

}} // scitbx::boost_python

#endif // SCITBX_BOOST_PYTHON_STL_MAP_AS_DICT_H
