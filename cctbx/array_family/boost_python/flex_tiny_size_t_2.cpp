#include <cctbx/boost_python/flex_fwd.h>

#include <scitbx/array_family/boost_python/flex_wrapper.h>

namespace scitbx { namespace af { namespace boost_python {
namespace {

  shared<size_t>
  column(
    const_ref<tiny<std::size_t, 2> > const& a,
    std::size_t i_column)
  {
    SCITBX_ASSERT(i_column < 2);
    shared<size_t> result((reserve(a.size())));
    for(std::size_t i=0;i<a.size();i++) {
      result.push_back(a[i][i_column]);
    }
    return result;
  }

} // namespace <anonymous>

  template <>
  struct flex_default_element<tiny<std::size_t, 2> >
  {
    static tiny<std::size_t, 2>
    get() { return tiny<std::size_t, 2>(0,0); }
  };

  void wrap_flex_tiny_size_t_2()
  {
    flex_wrapper<tiny<std::size_t, 2> >::plain("tiny_size_t_2")
      .def("column", column);
    ;
  }

}}} // namespace scitbx::af::boost_python
