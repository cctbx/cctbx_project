#include <boost/python/def.hpp>
#include <boost/python/list.hpp>
#include <boost/python/tuple.hpp>

#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/auto_array.h>
#include <boost/shared_array.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <vector>

/* silence warnings triggered by bad coding of the standard library
   on MacOS X, at least with XCode 3.2.2
 */
#pragma clang diagnostic ignored "-Wmismatched-tags"
#include <valarray>

namespace scitbx { namespace af { namespace boost_python {

namespace {

#define SCITBX_LOC_ST(S, T) \
    result.append(make_tuple(S, sizeof(T)));

#define SCITBX_LOC_T(T) SCITBX_LOC_ST(#T, T)

#define SCITBX_LOC_empty_container_sizes(T) \
  boost::python::list \
  empty_container_sizes_##T() \
  { \
    using namespace boost::python; \
    list result; \
    SCITBX_LOC_T(T) \
    SCITBX_LOC_T(af::ref<T>) \
    SCITBX_LOC_T(scitbx::auto_array<T>) \
    SCITBX_LOC_T(boost::shared_array<T>) \
    SCITBX_LOC_T(boost::ptr_vector<T>) \
    SCITBX_LOC_T(std::vector<T>) \
    SCITBX_LOC_T(af::shared<T>) \
    result.append( \
      make_tuple("af::shared [cumulative]", \
                  sizeof(af::shared<T>) + sizeof(af::sharing_handle))); \
    SCITBX_LOC_T(std::valarray<T>) \
    SCITBX_LOC_T(af::versa<T>) \
    typedef af::versa<T, af::flex_grid<> > versa_flex_grid; \
    SCITBX_LOC_ST("af::versa<" #T ", af::flex_grid<> >", versa_flex_grid) \
    result.append( \
      make_tuple("af::versa<" #T ", af::flex_grid<> > [cumulative]", \
                 sizeof(versa_flex_grid) + sizeof(af::sharing_handle))); \
    return result; \
  }

  SCITBX_LOC_empty_container_sizes(int)
  SCITBX_LOC_empty_container_sizes(double)

#undef SCITBX_LOC_empty_container_sizes
#undef SCITBX_LOC_T
#undef SCITBX_LOC_ST

} // namespace anonymous

  void
  wrap_empty_container_sizes()
  {
    using namespace boost::python;
    def("empty_container_sizes_int", empty_container_sizes_int);
    def("empty_container_sizes_double", empty_container_sizes_double);
  }

}}} // namespace scitbx::af::boost_python
