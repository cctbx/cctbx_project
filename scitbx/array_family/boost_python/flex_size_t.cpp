#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/boost_python/flex_pickle_single_buffered.h>
#include <boost/python/make_constructor.hpp>

namespace scitbx { namespace af { namespace boost_python {
namespace {

  flex<std::size_t>::type*
  from_stl_vector_unsigned(std::vector<unsigned> const& v)
  {
    af::shared<std::size_t> result(af::reserve(v.size()));
    for(std::size_t i=0;i<v.size();i++) {
      result.push_back(v[i]);
    }
    return new flex<std::size_t>::type(result, result.size());
  }

  bool
  next_permutation(af::ref<std::size_t> const& a)
  {
    return std::next_permutation(a.begin(), a.end());
  }

  af::shared<std::size_t>
  inverse_permutation(af::const_ref<std::size_t> const& self)
  {
    af::shared<std::size_t> result(self.size());
    af::ref<std::size_t> r = result.ref();
    for(std::size_t i=0;i<self.size();i++) {
      SCITBX_ASSERT(self[i] < self.size());
      r[self[i]] = i;
    }
    return result;
  }

} // namespace <anonymous>

  void wrap_flex_size_t()
  {
    using namespace boost::python;
    flex_wrapper<std::size_t>::integer("size_t", boost::python::scope())
      .def_pickle(flex_pickle_single_buffered<std::size_t>())
      .def("__init__", make_constructor(
        from_stl_vector_unsigned, default_call_policies()))
      .def("next_permutation", next_permutation)
      .def("inverse_permutation", inverse_permutation)
    ;
  }

}}} // namespace scitbx::af::boost_python
