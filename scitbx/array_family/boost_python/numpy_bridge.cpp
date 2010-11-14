#include <boost/python.hpp>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/flex_grid.h>
#include <boost_adaptbx/type_id_eq.h>

#if defined(SCITBX_HAVE_NUMPY_INCLUDE)
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
#endif

namespace scitbx { namespace af { namespace boost_python {

  void
  import_numpy_api_if_available()
  {
#if defined(SCITBX_HAVE_NUMPY_INCLUDE)
    import_array();
    boost::python::numeric::array::set_module_and_type("numpy", "ndarray");
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
    npy_intp dims[flex_grid<>::index_type::capacity()];
    flex_grid<> const& grid = O.accessor();
    int nd = static_cast<int>(grid.nd());
    for(unsigned i=0;i<nd;i++) {
      dims[i] = static_cast<npy_intp>(grid.all()[i]);
    }
    result = bp::object(bp::handle<>(
      PyArray_SimpleNewFromData(nd, dims, type_num, O.begin())));
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
      t[i] = static_cast<TargetType>(s[i]);
    }
  }

  template <typename ElementType>
  versa<ElementType, flex_grid<> >
  versa_flex_from_numpy_array(
    boost::python::numeric::array const& arr_obj)
  {
    namespace bp = boost::python;
#if !defined(SCITBX_HAVE_NUMPY_INCLUDE)
    throw std::runtime_error(numpy_api_not_available());
#else
    PyObject* arr_ptr = arr_obj.ptr();
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
    SCITBX_ASSERT(grid.size_1d() == PyArray_Size(arr_ptr));
    versa<ElementType, flex_grid<> > result(
      grid, init_functor_null<ElementType>());
    void* data = PyArray_DATA(arr_ptr);
    int type_num = PyArray_TYPE(arr_ptr);
#define SCITBX_LOC(tn) \
    (type_num == tn) { \
      copy_data_with_cast( \
        grid.size_1d(), \
        reinterpret_cast<ElementType*>(data), \
        result.begin()); \
    }
    if      SCITBX_LOC(NPY_BOOL)
    else if SCITBX_LOC(NPY_INT)
    else if SCITBX_LOC(NPY_LONG)
    else if SCITBX_LOC(NPY_FLOAT)
    else if SCITBX_LOC(NPY_DOUBLE)
    else if SCITBX_LOC(NPY_CDOUBLE)
    else if SCITBX_LOC(NPY_USHORT)
    else if SCITBX_LOC(NPY_UINT)
    else if SCITBX_LOC(NPY_ULONG)
    else if SCITBX_LOC(NPY_ULONGLONG)
    else {
      throw std::runtime_error(
        "Unsupported numpy.ndarray element type");
    }
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
    boost::python::numeric::array const& arr_obj) \
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
