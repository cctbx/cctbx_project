#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/boost_python/flex_pickle_single_buffered.h>
#include <scitbx/array_family/boost_python/byte_str.h>
#include <scitbx/array_family/boost_python/range_wrappers.h>
#include <scitbx/array_family/boost_python/numpy_bridge.hpp>
#include <scitbx/array_family/counts.h>
#include <scitbx/array_family/selections.h>
#include <scitbx/stl/map_fwd.h>
#include <boost/python/make_constructor.hpp>
#include <boost/python/args.hpp>
#include <vector>
#include <map>

namespace scitbx { namespace af { namespace boost_python {
namespace {

  flex<std::size_t>::type*
  from_stl_vector_unsigned(std::vector<unsigned> const& v)
  {
    shared<std::size_t> result(reserve(v.size()));
    for(std::size_t i=0;i<v.size();i++) {
      result.push_back(v[i]);
    }
    return new flex<std::size_t>::type(result, result.size());
  }

  shared<int>
  as_int(
    const_ref<std::size_t> const& O)
  {
    shared<int> result(O.size(), init_functor_null<int>());
    for(std::size_t i=0;i<O.size();i++) {
      result[i] = static_cast<int>(O[i]);
    }
    return result;
  }

  bool
  next_permutation(ref<std::size_t> const& a)
  {
    return std::next_permutation(a.begin(), a.end());
  }

  shared<std::size_t>
  inverse_permutation(const_ref<std::size_t> const& self)
  {
    shared<std::size_t> result(self.size());
    ref<std::size_t> r = result.ref();
    for(std::size_t i=0;i<self.size();i++) {
      SCITBX_ASSERT(self[i] < self.size());
      r[self[i]] = i;
    }
    return result;
  }

  std::size_t
  increment_and_track_up_from_zero(
    ref<std::size_t> const& O,
    const_ref<std::size_t> const& iselection)
  {
    std::size_t result = 0;
    for(std::size_t i=0;i<iselection.size();i++) {
      std::size_t ii = iselection[i];
      SCITBX_ASSERT(ii < O.size());
      if (O[ii]++ == 0) result++;
    }
    return result;
  }

  boost::python::tuple
  intersection_i_seqs(
    const_ref<std::size_t> const& self,
    const_ref<std::size_t> const& other)
  {
    intersection_with_tracking<std::size_t, std::size_t> proxy(
      self, other, /*track_matching_elements*/ false, /*track_i_seqs*/ true);
    return boost::python::make_tuple(proxy.self_i_seqs, proxy.other_i_seqs);
  }

} // namespace <anonymous>

  void wrap_flex_size_t()
  {
    using namespace boost::python;
    using boost::python::arg;
    flex_wrapper<std::size_t>::integer("size_t", scope())
      .def_pickle(flex_pickle_single_buffered<std::size_t>())
      .def("__init__", make_constructor(
        from_stl_vector_unsigned, default_call_policies()))
      .def("__init__", make_constructor(
        flex_size_t_from_numpy_array, default_call_policies()))
      .def("copy_to_byte_str",
        copy_to_byte_str<versa<std::size_t, flex_grid<> > >)
      .def("as_int", as_int)
      .def("intersection",
        (shared<std::size_t>(*)(
          const_ref<std::size_t> const&,
          const_ref<std::size_t> const&))
        intersection, (arg("self"), arg("other")))
      .def("intersection_i_seqs", intersection_i_seqs, (
        arg("self"), arg("other")))
      .def("counts", counts<std::size_t, std::map<long, long> >::unlimited)
      .def("counts", counts<std::size_t, std::map<long, long> >::limited, (
        arg("max_keys")))
      .def("next_permutation", next_permutation)
      .def("inverse_permutation", inverse_permutation)
      .def("increment_and_track_up_from_zero",
        increment_and_track_up_from_zero, (arg("iselection")))
      .def("as_numpy_array", flex_size_t_as_numpy_array, (
        arg("optional")=false))
    ;
    def(
      "size_t_from_byte_str",
      shared_from_byte_str<std::size_t>,
      (arg("byte_str")));
    range_wrappers<std::size_t, long, range_args::unsigned_check>::wrap(
      "size_t_range");
  }

}}} // namespace scitbx::af::boost_python
