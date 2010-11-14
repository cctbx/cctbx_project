#include <boost/python/numeric.hpp>

namespace scitbx { namespace af { namespace boost_python {

#define SCITBX_LOC(pyname, ElementType) \
  boost::python::object \
  flex_##pyname##_as_numpy_array( \
    ref<ElementType, flex_grid<> > const& O); \
 \
  versa<ElementType, flex_grid<> >* \
  flex_##pyname##_from_numpy_array( \
    boost::python::numeric::array const& arr_obj);

  SCITBX_LOC(bool, bool)
  SCITBX_LOC(int, int)
  SCITBX_LOC(long, long)
  SCITBX_LOC(float, float)
  SCITBX_LOC(double, double)
#define SCITBX_LOC2 std::complex<double>
  SCITBX_LOC(complex_double, SCITBX_LOC2)
#undef SCITBX_LOC2
  SCITBX_LOC(size_t, std::size_t)

#undef SCITBX_LOC

}}} // namespace scitbx::af::boost_python
