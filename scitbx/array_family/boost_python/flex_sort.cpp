#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/array_family/sort.h>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

namespace scitbx { namespace af { namespace boost_python {

namespace {

  template <typename ElementType>
  void
  sort_permutation_wrapper()
  {
    using namespace boost::python;
    def("sort_permutation",
      // cast needed for Intel C++ 9.1
      (shared<std::size_t>(*)(const_ref<ElementType> const&, bool, bool))
        af::sort_permutation, (
          arg("data"), arg("reverse")=false, arg("stable")=true));
  }

} // namespace <anonymous>

  void wrap_flex_sort()
  {
    sort_permutation_wrapper<int>();
    sort_permutation_wrapper<std::size_t>();
    sort_permutation_wrapper<double>();
  }

}}} // namespace scitbx::af::boost_python
