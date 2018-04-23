#include <Python.h>
#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include <boost_adaptbx/floating_point_exceptions.h>
#include <boost_adaptbx/type_id_eq.h>
#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/array_family/versa.h>

using namespace scitbx::af;

// Only include the actual conversion functions if we have numpy
#if defined(SCITBX_HAVE_NUMPY_INCLUDE)

// Declare how we're using the numpy API.
// This is the unique symbol that boost numpy uses - since we defer
// activation of the numpy API to boost, we need to use the same symbol.
#define PY_ARRAY_UNIQUE_SYMBOL BOOST_NUMPY_ARRAY_API
// Numpy initialization done by boost python, so we don't want import_array
#define NO_IMPORT_ARRAY
#include <numpy/arrayobject.h>

/// Returns a numpy type description object describing a built-in C++ type.
///
/// This uses the boost::python::numpy code to properly abstract over platform
/// differences, and lets us tie into the boost::python::numpy API, but as
/// minimally as possible.
///
/// @tparam   T   The C++ type to get the numpy equivalent for
/// @returns  A new reference to the equivalent numpy type object
template <typename T> PyArray_Descr *get_builtin_dtype() {
  boost::python::object dtype = boost::python::numpy::dtype::get_builtin<T>();
  if (PyObject_IsInstance(
          dtype.ptr(), reinterpret_cast<PyObject *>(&PyArrayDescr_Type)) != 1) {
    throw std::runtime_error("Didn't get expected type object");
  }
  PyArray_Descr *dtp = reinterpret_cast<PyArray_Descr *>(dtype.ptr());
  Py_INCREF(dtp);
  return dtp;
}

/// Do the actual conversion from flex to numpy
template <typename ElementType>
boost::python::object
ref_flex_as_numpy_array(scitbx::af::ref<ElementType, flex_grid<> > const &O) {
  // Extract the number and size of dimensions
  npy_intp dims[flex_grid<>::index_type::capacity_value];
  flex_grid<> const &grid = O.accessor();
  int nd = static_cast<int>(grid.nd());
  for (unsigned i = 0; i < nd; i++) {
    dims[i] = static_cast<npy_intp>(grid.all()[i]);
  }

  // Get the equivalent dtype for this C++ type
  PyArray_Descr *dtype = get_builtin_dtype<ElementType>();
  PyArrayObject *result = reinterpret_cast<PyArrayObject *>(
      PyArray_SimpleNewFromDescr(nd, dims, dtype));

  // Copy over the data from our flex array
  ElementType *result_data =
      reinterpret_cast<ElementType *>(PyArray_DATA(result));
  std::copy(O.begin(), O.end(), result_data);
  // Return the python object
  boost::python::handle<> handle(reinterpret_cast<PyObject *>(result));
  return boost::python::object(handle);
}

/// Do actual conversion from numpy to a flex versa-grid
template <typename ElementType>
versa<ElementType, flex_grid<> >
versa_flex_from_numpy_array(boost::python::numpy::ndarray const &arr_obj) {
  // Convert this to a plain API pointer to work with the API
  PyArrayObject *arr_ptr = reinterpret_cast<PyArrayObject *>(arr_obj.ptr());
  // Make sanity checks on the input array
  if (!PyArray_ISCONTIGUOUS(arr_ptr)) {
    throw std::invalid_argument("numpy.ndarray instance is not contiguous");
  }

  // Build a dimensions array in the correct type for flex_grid
  flex_grid<>::index_type all; // e.g. a small<long, 10>
  npy_intp ndim = PyArray_NDIM(arr_ptr);
  SCITBX_ASSERT(ndim <= all.capacity());
  npy_intp *dims = PyArray_DIMS(arr_ptr);
  for (unsigned i = 0; i < ndim; i++) {
    all.push_back(static_cast<flex_grid<>::index_type::value_type>(dims[i]));
  }

  // Create the new flex_grid to hold all the data
  flex_grid<> grid(all);
  // Double-check that we have the sizes correct
  const npy_intp n_items = PyArray_SIZE(arr_ptr);
  SCITBX_ASSERT(grid.size_1d() == n_items);

  // Now, use this flex_grid to create the final array object
  versa<ElementType, flex_grid<> > result(grid,
                                          init_functor_null<ElementType>());

  // Now, convert the original object to the same data type as the flex
  // array so that we can do a direct copy
  PyArray_Descr *dtype = get_builtin_dtype<ElementType>();
  PyArrayObject *converted_arr =
      reinterpret_cast<PyArrayObject *>(PyArray_CastToType(arr_ptr, dtype, 0));

  // Now we want to copy over the data from the type-correct array
  ElementType *data =
      reinterpret_cast<ElementType *>(PyArray_DATA(converted_arr));
  std::copy(data, data + (size_t)n_items, result.begin());
  Py_DECREF(converted_arr);

  return result;
}

#endif

char const *numpy_api_not_available() { return "numpy API not available"; }

///////////////////////////////////////////////////////////////////////////////
// Define all the public functions with numpy/no-numpy specific behaviour

namespace scitbx {
namespace af {
namespace boost_python {

void import_numpy_api_if_available() {
// If no numpy, then this function does nothing
#if defined(SCITBX_HAVE_NUMPY_INCLUDE)
  /* On Snow Leopard, numpy is shipped with the system
   and the call import_array() triggers a floating-point exception,
   which then crashes the program thanks to the default unforgiving
   FP trapping policy of the cctbx.
   Temporarily changing to a more lenient policy is the simplest
   possible fix.
   */
  using namespace boost_adaptbx::floating_point;
  exception_trapping guard(exception_trapping::dont_trap);
  boost::python::numpy::initialize();
#endif
}

template <typename ElementType>
boost::python::object
flex_as_numpy_array(ref<ElementType, flex_grid<> > const &from, bool optional) {
#if defined(SCITBX_HAVE_NUMPY_INCLUDE)
  return ref_flex_as_numpy_array(from);
#else
  // If no numpy, we have the option to return None instead of throwing
  if (!optional) {
    throw std::runtime_error(numpy_api_not_available());
  }
  // Return None
  return boost::python::object();
#endif
}

template <typename ElementType>
versa<ElementType, flex_grid<> > *
flex_from_numpy_array(boost::python::numpy::ndarray const &array) {
#if defined(SCITBX_HAVE_NUMPY_INCLUDE)
  return new versa<ElementType, flex_grid<> >(
      versa_flex_from_numpy_array<ElementType>(array));
#else
  throw std::runtime_error(numpy_api_not_available());
#endif
}

// Create a convenience macro that explicitly instantiates these functions
// for a specific type.
#define EXPLICIT_INSTANTIATE_FLEX_CONVERTERS(fortype)                          \
  template boost::python::object flex_as_numpy_array<fortype>(                 \
      ref<fortype, flex_grid<> > const &from, bool optional);                  \
  template versa<fortype, flex_grid<> > *flex_from_numpy_array<fortype>(       \
      boost::python::numpy::ndarray const &array);

// Explicitly instantiate converters for all flex types
EXPLICIT_INSTANTIATE_FLEX_CONVERTERS(bool)
EXPLICIT_INSTANTIATE_FLEX_CONVERTERS(int)
EXPLICIT_INSTANTIATE_FLEX_CONVERTERS(long)
EXPLICIT_INSTANTIATE_FLEX_CONVERTERS(float)
EXPLICIT_INSTANTIATE_FLEX_CONVERTERS(double)
EXPLICIT_INSTANTIATE_FLEX_CONVERTERS(std::complex<double>)
EXPLICIT_INSTANTIATE_FLEX_CONVERTERS(std::size_t)

} // namespace boost_python
} // namespace af
} // namespace scitbx
