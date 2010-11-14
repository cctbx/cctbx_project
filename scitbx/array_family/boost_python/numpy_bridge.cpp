#include <boost/python.hpp>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/flex_grid.h>
#include <boost_adaptbx/type_id_eq.h>

#if defined(SCITBX_HAVE_NUMPY_INCLUDE)
#include <numpy/arrayobject.h>
#endif

namespace scitbx { namespace af { namespace boost_python {

  void
  import_numpy_api_if_available()
  {
#if defined(SCITBX_HAVE_NUMPY_INCLUDE)
    import_array();
#endif
  }

  template <typename ElementType>
  boost::python::object
  ref_flex_as_numpy_array(
    ref<ElementType, flex_grid<> > const& O,
    int typenum)
  {
    namespace bp = boost::python;
    bp::object result;
#if defined(SCITBX_HAVE_NUMPY_INCLUDE)
    npy_intp dims[flex_grid<>::index_type::capacity()];
    flex_grid<> const& grid = O.accessor();
    int nd = static_cast<int>(grid.nd());
    for(unsigned i=0;i<nd;i++) {
      dims[i] = static_cast<npy_intp>(grid.all()[i]);
    }
    result = bp::object(bp::handle<>(
      PyArray_SimpleNewFromData(nd, dims, typenum, O.begin())));
#endif
    return result;
  }

#define SCITBX_LOC(ElementType, typenum) \
  boost::python::object \
  ref_flex_as_numpy_array( \
    ref<ElementType, flex_grid<> > const& O) \
  { \
    return ref_flex_as_numpy_array(O, typenum); \
  }

#if !defined(SCITBX_HAVE_NUMPY_INCLUDE)
#define NPY_BOOL 0
#define NPY_INT 0
#define NPY_LONG 0
#define NPY_FLOAT 0
#define NPY_DOUBLE 0
#define NPY_CDOUBLE 0
#define NPY_ULONG 0
#define NPY_ULONGLONG 0
#endif

  SCITBX_LOC(bool, NPY_BOOL);
  SCITBX_LOC(int, NPY_INT);
  SCITBX_LOC(long, NPY_LONG);
  SCITBX_LOC(float, NPY_FLOAT);
  SCITBX_LOC(double, NPY_DOUBLE);
  SCITBX_LOC(std::complex<double>, NPY_CDOUBLE);

#if defined(BOOST_ADAPTBX_TYPE_ID_SIZE_T_EQ_UNSIGNED_LONG)
  SCITBX_LOC(std::size_t, NPY_ULONG);
#elif defined(BOOST_ADAPTBX_TYPE_ID_SIZE_T_EQ_UNSIGNED_LONG_LONG)
  SCITBX_LOC(std::size_t, NPY_ULONGLONG);
#else
# error Unknown size_t: review boost_adaptbx/type_id_eq.h
#endif

#undef SCITBX_LOC

}}} // namespace scitbx::af::boost_python
