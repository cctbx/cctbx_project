#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/boost_python/flex_pickle_single_buffered.h>

namespace scitbx { namespace af { namespace boost_python {
namespace {

  bool
  next_permutation(af::ref<std::size_t> const& a)
  {
    return std::next_permutation(a.begin(), a.end());
  }

} // namespace <anonymous>

  void wrap_flex_size_t()
  {
    flex_wrapper<std::size_t>::integer("size_t", boost::python::scope())
      .def_pickle(flex_pickle_single_buffered<std::size_t>())
      .def("next_permutation", next_permutation);
    ;
  }

}}} // namespace scitbx::af::boost_python
