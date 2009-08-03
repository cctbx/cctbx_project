#ifndef SCITBX_BOOST_PYTHON_SEQUENCE_AS_ARRAY_H
#define SCITBX_BOOST_PYTHON_SEQUENCE_AS_ARRAY_H

#include <boost/python/extract.hpp>
#include <boost/python/object.hpp>
#include <boost/numeric/conversion/cast.hpp>

namespace scitbx { namespace boost_python {

  template <typename ArrayType>
  ArrayType
  sequence_as(
    boost::python::object const& seq)
  {
    namespace bp = boost::python;
    bp::ssize_t n = bp::len(seq);
    ArrayType result;
    result.reserve(boost::numeric_cast<std::size_t>(n));
    typedef typename ArrayType::value_type value_type;
    for(bp::ssize_t i=0;i<n;i++) {
      result.push_back(bp::extract<value_type>(seq[i])());
    }
    return result;
  }

}} // scitbx::boost_python

#endif // GUARD
