#include <stdint.h>
#include <boost/python/numpy.hpp>

namespace scitbx { namespace af { namespace boost_python {

#define SCITBX_LOC(pyname, ElementType) \
  boost::python::object \
  flex_##pyname##_as_numpy_array( \
    ref<ElementType, flex_grid<> > const& O, \
    bool optional=false); \
 \
  versa<ElementType, flex_grid<> >* \
  flex_##pyname##_from_numpy_array( \
    boost::python::numpy::ndarray const& arr_obj);

  SCITBX_LOC(bool, bool)
  SCITBX_LOC(int, int)
  SCITBX_LOC(long, long)
  SCITBX_LOC(float, float)
  SCITBX_LOC(double, double)
#define SCITBX_LOC2 std::complex<double>
  SCITBX_LOC(complex_double, SCITBX_LOC2)
#undef SCITBX_LOC2
  SCITBX_LOC(size_t, std::size_t)
  SCITBX_LOC(int8, int8_t)
  SCITBX_LOC(int16, int16_t)
  // SCITBX_LOC(int32, int32_t)
  #if defined(_MSC_VER)
  SCITBX_LOC(int64, int64_t)
  #endif
  SCITBX_LOC(uint8, uint8_t)
  SCITBX_LOC(uint16, uint16_t)
  SCITBX_LOC(uint32, uint32_t)
  // SCITBX_LOC(uint64, uint64_t)
#undef SCITBX_LOC

}}} // namespace scitbx::af::boost_python
