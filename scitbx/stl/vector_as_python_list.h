#ifndef SCITBX_BOOST_PYTHON_STL_VECTOR_AS_PYTHON_LIST_H
#define SCITBX_BOOST_PYTHON_STL_VECTOR_AS_PYTHON_LIST_H

#include <boost/python/list.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/detail/api_placeholder.hpp>

namespace scitbx { namespace stl { namespace boost_python {

  template <typename VectorType>
  boost::python::list
  vector_as_list(VectorType const& v)
  {
    boost::python::list result;
    for(std::size_t i=0;i<v.size();i++) {
      result.append(v[i]);
    }
    return result;
  }

  template <typename VectorType>
  void
  update_vector_from_list(VectorType& v, boost::python::list values)
  {
    typedef typename VectorType::value_type value_type;
    using namespace boost::python;
    std::size_t len_values = len(values);
    v.reserve(len_values);
    for(std::size_t i=0;i<len_values;i++) {
      value_type value = extract<value_type>(values[i])();
      v.push_back(value);
    }
  }

}}} // namespace scitbx::stl::boost_python

#endif // SCITBX_BOOST_PYTHON_STL_VECTOR_AS_PYTHON_LIST_H
