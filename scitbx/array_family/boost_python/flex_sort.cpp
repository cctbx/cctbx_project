#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/array_family/sort.h>
#include <boost/python/def.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/args.hpp>

namespace scitbx { namespace af { namespace boost_python {

namespace {

  BOOST_PYTHON_FUNCTION_OVERLOADS(
    sort_permutation_overloads, sort_permutation, 1, 2)

  template <typename ElementType>
  struct sort_permutation_wrapper
  {
    static void
    wrap()
    {
      using namespace boost::python;
      boost::python::def("sort_permutation",
        (shared<std::size_t>(*)(const_ref<ElementType> const&, bool)) 0,
        sort_permutation_overloads((
          arg_("data"), arg_("reverse")=false)));
      }
  };

} // namespace <anonymous>

  void wrap_flex_sort()
  {
    sort_permutation_wrapper<int>::wrap();
    sort_permutation_wrapper<std::size_t>::wrap();
    sort_permutation_wrapper<double>::wrap();
  }

}}} // namespace scitbx::af::boost_python
