#ifndef SCITBX_BOOSTPYTHON_PAIRASTUPLE_HPP_
#define SCITBX_BOOSTPYTHON_PAIRASTUPLE_HPP_

#include <boost/python/tuple.hpp>

namespace scitbx
{

namespace boost_python
{

template<class Pair>
struct PairToTupleConverter
{
  static PyObject* convert(Pair const& pair)
  {
    return boost::python::incref(
      boost::python::make_tuple( pair.first, pair.second ).ptr()
      );
  }
};

} // namespace boost_python
} // namespace scitbx

#endif // SCITBX_BOOSTPYTHON_PAIRASTUPLE_HPP_

