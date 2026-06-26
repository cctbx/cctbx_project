#include <stdint.h>
#include <boost/python/numpy.hpp>
#include <boost/python.hpp>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/flex_grid.h>
#include <boost_adaptbx/type_id_eq.h>
#include <boost_adaptbx/floating_point_exceptions.h>

#if defined(SCITBX_HAVE_NUMPY_INCLUDE)
// https://docs.scipy.org/doc/numpy/reference/c-api.array.html
#define PY_ARRAY_UNIQUE_SYMOL FLEX_ARRAY_API
// #define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>
#else
#define NPY_BOOL 0
#define NPY_INT 0
#define NPY_LONG 0
#define NPY_FLOAT 0
#define NPY_DOUBLE 0
#define NPY_CDOUBLE 0
#define NPY_USHORT 0
#define NPY_UINT 0
#define NPY_ULONG 0
#define NPY_ULONGLONG 0
#define NPY_INT8 0
#define NPY_INT16 0
#define NPY_INT32 0
#define NPY_INT64 0
#define NPY_UINT8 0
#define NPY_UINT16 0
#define NPY_UINT32 0
#define NPY_UINT64 0
#endif

#if PY_MAJOR_VERSION >= 3
#define IS_PY3K
#endif

namespace scitbx { namespace af { namespace boost_python {

// import_array returns NULL in python3, and somehow that makes this function
// need a return type as well.
#ifdef IS_PY3K
  void*
#else
  void
#endif
  import_numpy_api_if_available()
  {
#if defined(SCITBX_HAVE_NUMPY_INCLUDE)
    {
    // Fix segmentation faults
    // http://boostorg.github.io/python/doc/html/numpy/tutorial/simple.html
    using namespace boost_adaptbx::floating_point;
    exception_trapping guard(exception_trapping::dont_trap);
    Py_Initialize();
    boost::python::numpy::initialize();
    import_array();
    }
#endif
#ifdef IS_PY3K
    return NULL;
#endif
  }

  char const*
  numpy_api_not_available() { return "numpy API not available"; }

  template <typename ElementType>
  boost::python::object
  ref_flex_as_numpy_array(
    ref<ElementType, flex_grid<> > const& O,
    bool optional,
    int type_num)
  {
    namespace bp = boost::python;
    bp::object result;
#if !defined(SCITBX_HAVE_NUMPY_INCLUDE)
    if (!optional) {
      throw std::runtime_error(numpy_api_not_available());
    }
#else
    npy_intp dims[flex_grid<>::index_type::capacity_value];
    flex_grid<> const& grid = O.accessor();
    int nd = static_cast<int>(grid.nd());
    for(unsigned i=0;i<nd;i++) {
      dims[i] = static_cast<npy_intp>(grid.all()[i]);
    }
    result = bp::object(bp::handle<>(PyArray_SimpleNew(nd, dims, type_num)));
    ElementType* result_data = reinterpret_cast<ElementType*>(
      PyArray_DATA(reinterpret_cast<PyArrayObject*>(result.ptr())));
    std::copy(O.begin(), O.end(), result_data);
#endif
    return result;
  }

  template <typename SourceType, typename TargetType>
  void
  copy_data_with_cast(
    std::size_t n,
    SourceType const* s,
    TargetType* t)
  {
    for(std::size_t i=0;i<n;i++) {
      t[i] = static_cast<TargetType>(static_cast<SourceType>(s[i]));
    }
  }

  template <typename ElementType>
  versa<ElementType, flex_grid<> >
  versa_flex_from_numpy_array(
    boost::python::numpy::ndarray const& arr_obj)
  {
    namespace bp = boost::python;
#if !defined(SCITBX_HAVE_NUMPY_INCLUDE)
    throw std::runtime_error(numpy_api_not_available());
#else
    PyObject* obj_ptr = arr_obj.ptr();
    PyArrayObject* arr_ptr = reinterpret_cast<PyArrayObject*>(obj_ptr);
    if (!PyArray_Check(arr_ptr)) {
      throw std::invalid_argument(
        "Expected a numpy.ndarray instance");
    }
    if (!PyArray_ISCONTIGUOUS(arr_ptr)) {
      throw std::invalid_argument(
        "numpy.ndarray instance is not contiguous");
    }
    typedef typename flex_grid<>::index_type fgit;
    typedef typename fgit::value_type fgivt;
    fgit all;
    npy_intp ndim = PyArray_NDIM(arr_ptr);
    SCITBX_ASSERT(ndim <= all.capacity());
    npy_intp* dims = PyArray_DIMS(arr_ptr);
    for(unsigned i=0;i<ndim;i++) {
      all.push_back(static_cast<fgivt>(dims[i]));
    }
    flex_grid<> grid(all);
    SCITBX_ASSERT(grid.size_1d() == PyArray_Size(obj_ptr));
    versa<ElementType, flex_grid<> > result(
      grid, init_functor_null<ElementType>());
    void* data = PyArray_DATA(arr_ptr);
    int type_num = PyArray_TYPE(arr_ptr);
#define SCITBX_LOC(tn, t) \
    (type_num == tn) { \
      copy_data_with_cast<t, ElementType>( \
        grid.size_1d(), \
        reinterpret_cast<t*>(data), \
        result.begin()); \
    }
#define SCITBX_LOC_COMPLEX(tn) \
    (type_num == tn) { \
      copy_data_with_cast<ElementType, ElementType>( \
        grid.size_1d(), \
        reinterpret_cast<ElementType*>(data), \
        result.begin()); \
    }
    if      SCITBX_LOC(NPY_BOOL, npy_bool)
    else if SCITBX_LOC(NPY_INT, npy_int)
    else if SCITBX_LOC(NPY_LONG, npy_long)
    else if SCITBX_LOC(NPY_FLOAT, npy_float)
    else if SCITBX_LOC(NPY_DOUBLE, npy_double)
    else if SCITBX_LOC_COMPLEX(NPY_CDOUBLE)
    else if SCITBX_LOC(NPY_USHORT, npy_ushort)
    else if SCITBX_LOC(NPY_UINT, npy_uint)
    else if SCITBX_LOC(NPY_ULONG, npy_ulong)
    else if SCITBX_LOC(NPY_ULONGLONG, npy_ulonglong)
    else if SCITBX_LOC(NPY_INT8, npy_int8)
    else if SCITBX_LOC(NPY_INT16, npy_int16)
    else if SCITBX_LOC(NPY_INT32, npy_int32)
    else if SCITBX_LOC(NPY_INT64, npy_int64)
    else if SCITBX_LOC(NPY_UINT8, npy_uint8)
    else if SCITBX_LOC(NPY_UINT16, npy_uint16)
    else if SCITBX_LOC(NPY_UINT32, npy_uint32)
    else if SCITBX_LOC(NPY_UINT64, npy_uint64)
    else {
      throw std::runtime_error(
        "Unsupported numpy.ndarray element type");
    }
#undef SCITBX_LOC_COMPLEX
#undef SCITBX_LOC
    return result;
#endif
  }

#define SCITBX_LOC(pyname, ElementType, type_num) \
  boost::python::object \
  flex_##pyname##_as_numpy_array( \
    ref<ElementType, flex_grid<> > const& O, \
    bool optional) \
  { \
    return ref_flex_as_numpy_array(O, optional, type_num); \
  } \
 \
  versa<ElementType, flex_grid<> >* \
  flex_##pyname##_from_numpy_array( \
    boost::python::numpy::ndarray const& arr_obj) \
  { \
    return new versa<ElementType, flex_grid<> >( \
      versa_flex_from_numpy_array<ElementType >(arr_obj)); \
  }

  SCITBX_LOC(bool, bool, NPY_BOOL);
  SCITBX_LOC(int, int, NPY_INT);
  SCITBX_LOC(long, long, NPY_LONG);
  SCITBX_LOC(float, float, NPY_FLOAT);
  SCITBX_LOC(double, double, NPY_DOUBLE);
  SCITBX_LOC(complex_double, std::complex<double>, NPY_CDOUBLE);
  SCITBX_LOC(int8, int8_t, NPY_INT8);
  SCITBX_LOC(int16, int16_t, NPY_INT16);
  // SCITBX_LOC(int32, int32_t, NPY_INT32);
  #if defined(_MSC_VER)
  SCITBX_LOC(int64, int64_t, NPY_INT64);
  #endif
  SCITBX_LOC(uint8, uint8_t, NPY_UINT8);
  SCITBX_LOC(uint16, uint16_t, NPY_UINT16);
  SCITBX_LOC(uint32, uint32_t, NPY_UINT32);
  SCITBX_LOC(uint64, uint64_t, NPY_UINT64);

#if defined(BOOST_ADAPTBX_TYPE_ID_SIZE_T_EQ_UNSIGNED_SHORT)
  SCITBX_LOC(size_t, std::size_t, NPY_USHORT);
#elif defined(BOOST_ADAPTBX_TYPE_ID_SIZE_T_EQ_UNSIGNED)
  SCITBX_LOC(size_t, std::size_t, NPY_UINT);
#elif defined(BOOST_ADAPTBX_TYPE_ID_SIZE_T_EQ_UNSIGNED_LONG)
  SCITBX_LOC(size_t, std::size_t, NPY_ULONG);
#elif defined(BOOST_ADAPTBX_TYPE_ID_SIZE_T_EQ_UNSIGNED_LONG_LONG)
  SCITBX_LOC(size_t, std::size_t, NPY_ULONGLONG);
#else
# error Unknown size_t: review boost_adaptbx/type_id_eq.h
#endif

#undef SCITBX_LOC

}}} // namespace scitbx::af::boost_python
