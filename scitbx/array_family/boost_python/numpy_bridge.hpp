#ifndef NUMPY_BRIDGE_H
#define NUMPY_BRIDGE_H

#include <boost/python/numpy.hpp>
#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/array_family/versa.h>

namespace scitbx {
namespace af {
namespace boost_python {

/// Do any numpy API imports that might be required (only call once)
void import_numpy_api_if_available();

/// Convert a flex versa-grid object to a numpy array
///
/// @tparam ElementType   The value-type for the flex array. A numpy datatype
///                       will be searched for to match this type.
/// @param  from          The versa-grid array to convert from
/// @param  optional      Ignored if numpy is available. If numpy is missing,
///                       and optional is True, the returned python object is
///                       None. Otherwise, if optional is False then an
///                       std::invalid_argument exception is thrown.
/// @returns A new numpy array object
///
/// @throws std::runtime_error      If no matching numpy type could be found,
///                                 or numpy unavailable and optional=false
template <typename ElementType>
boost::python::object
flex_as_numpy_array(ref<ElementType, flex_grid<> > const &from, bool optional);

/// Convert a numpy array to a flex object.
///
/// The contents of the numpy array will be converted to a numpy datatype
/// matching the target flex array before copying.
///
/// @tparam   ElementType         The value-type for the output flex array
/// @param    array               A numpy array wrapped by boost::python
/// @returns  A flex versa-grid array
///
/// @throws std::runtime_error    If numpy is not available, or if no
///                               matching numpy data type could be found.
/// @throws std::invalid_argument If the numpy array is non-contiguous
template <typename ElementType>
versa<ElementType, flex_grid<> > *
flex_from_numpy_array(boost::python::numpy::ndarray const &array);

} // namespace boost_python
} // namespace af
} // namespace scitbx

#endif
