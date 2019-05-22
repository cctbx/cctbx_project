#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/maptbx/accessors/c_grid_padded_p1.h>
#include <cctbx/maptbx/statistics.h>
#include <boost/python/def.hpp>

#if PY_MAJOR_VERSION >= 3
#define IS_PY3K

void destroy_pyobject(PyObject* obj_ptr) {
  free(PyCapsule_GetPointer(obj_ptr, NULL));
}

#endif

namespace cctbx { namespace maptbx { namespace boost_python {

namespace {

  // Based on code by N.K. Sauter
  template <typename InpFloatType,
            typename OutFloatType>
  struct as_CObjectZYX
  {
    typedef c_grid_padded_p1<3>::index_type index_type;

    static boost::python::object
    convert(
      af::const_ref<InpFloatType, c_grid_padded_p1<3> > const& map_unit_cell,
      index_type const& first,
      index_type const& last,
      bool apply_sigma_scaling)
    {
      InpFloatType mean = 0;
      InpFloatType sigma = 0;
      if (apply_sigma_scaling) {
        statistics<InpFloatType>
          map_statistics(af::const_ref<InpFloatType, af::flex_grid<> >(
            map_unit_cell.begin(),
            map_unit_cell.accessor().as_flex_grid()));
        mean = map_statistics.mean();
        sigma = map_statistics.sigma();
        if (sigma == 0) sigma = 1;
      }
      OutFloatType* out_mem = reinterpret_cast<OutFloatType*>(
        malloc(out_size(first, last) * sizeof(OutFloatType)));
      OutFloatType* out_ptr = out_mem;
      index_type out_pt;
      for (out_pt[2] = first[2]; out_pt[2] <= last[2]; out_pt[2]++) {
      for (out_pt[1] = first[1]; out_pt[1] <= last[1]; out_pt[1]++) {
      for (out_pt[0] = first[0]; out_pt[0] <= last[0]; out_pt[0]++) {
        InpFloatType val = map_unit_cell(out_pt);
        if (apply_sigma_scaling) val = (val - mean) / sigma;
        *out_ptr++ = static_cast<OutFloatType>(val);
      }}}
      return boost::python::object(boost::python::handle<>(
#ifdef IS_PY3K
        PyCapsule_New(out_mem, NULL, destroy_pyobject)));
#else
        PyCObject_FromVoidPtr(out_mem, free)));
#endif
    }

    static std::size_t
    out_size(index_type const& first, index_type const& last)
    {
      std::size_t result = 1;
      for(std::size_t i=0;i<3;i++) {
        CCTBX_ASSERT(last[i] >= first[i]);
        result *= (last[i] - first[i] + 1);
      }
      return result;
    }
  };

} // namespace <anoymous>

  void wrap_pymol_interface()
  {
    using namespace boost::python;
    def("as_CObjectZYX", as_CObjectZYX<float, float>::convert);
    def("as_CObjectZYX", as_CObjectZYX<double, float>::convert);
  }

}}} // namespace cctbx::maptbx::boost_python
