#include <boost/python/def.hpp>
#include <boost/python/list.hpp>

#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/flex_grid.h>
#include <scitbx/array_family/shared.h>
#include <boost/shared_array.hpp>
#include <boost/format.hpp>
#include <vector>

namespace scitbx { namespace af { namespace boost_python {

namespace {

#define SCITBX_LOC_ST(S, T) \
    result.append((boost::format("sizeof(" S ") = %lu") % sizeof(T)).str());

#define SCITBX_LOC_T(T) SCITBX_LOC_ST(#T, T)

#define SCITBX_LOC_show_sizes(T) \
  boost::python::list \
  show_sizes_##T() \
  { \
    boost::python::list result; \
    SCITBX_LOC_T(T) \
    SCITBX_LOC_T(boost::shared_array<T>) \
    SCITBX_LOC_T(std::vector<T>) \
    SCITBX_LOC_T(af::shared<T>) \
    SCITBX_LOC_T(af::versa<T>) \
    typedef af::versa<T, af::flex_grid<> > versa_flex_grid; \
    SCITBX_LOC_ST("af::versa<" #T ", af::flex_grid<> >", versa_flex_grid) \
    return result; \
  }

  SCITBX_LOC_show_sizes(int)
  SCITBX_LOC_show_sizes(double)

#undef SCITBX_LOC_show_sizes
#undef SCITBX_LOC_T
#undef SCITBX_LOC_ST

} // namespace anonymous

  void
  wrap_show_sizes()
  {
    using namespace boost::python;
    def("show_sizes_int", show_sizes_int);
    def("show_sizes_double", show_sizes_double);
  }

}}} // namespace scitbx::af::boost_python
