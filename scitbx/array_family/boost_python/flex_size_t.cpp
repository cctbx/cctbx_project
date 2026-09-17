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

  inline //identify as inline to avoid an unused-function compiler warning
  flex<std::size_t>::type*
  from_stl_vector_unsigned(std::vector<unsigned> const& v)
  {
    shared<std::size_t> result(reserve(v.size()));
    for(std::size_t i=0;i<v.size();i++) {
      result.push_back(v[i]);
    }
    return new flex<std::size_t>::type(result, result.size());
  }

  template<typename uIntType>
  shared<int>
  as_int(
    const_ref<uIntType> const& O)
  {
    shared<int> result(O.size(), init_functor_null<int>());
    for(std::size_t i=0;i<O.size();i++) {
      result[i] = static_cast<int>(O[i]);
    }
    return result;
  }

  template<typename uIntType>
  bool
  next_permutation(ref<uIntType> const& a)
  {
    return std::next_permutation(a.begin(), a.end());
  }

  template<typename uIntType>
  shared<uIntType>
  inverse_permutation(const_ref<uIntType> const& self)
  {
    shared<uIntType> result(self.size());
    ref<uIntType> r = result.ref();
    for(std::size_t i=0;i<self.size();i++) {
      SCITBX_ASSERT(self[i] < self.size());
      r[self[i]] = i;
    }
    return result;
  }

  template<typename uIntType>
  uIntType
  increment_and_track_up_from_zero(
    ref<uIntType> const& O,
    const_ref<uIntType> const& iselection)
  {
    uIntType result = 0;
    for(std::size_t i=0;i<iselection.size();i++) {
      uIntType ii = iselection[i];
      SCITBX_ASSERT(ii < O.size());
      if (O[ii]++ == 0) result++;
    }
    return result;
  }

  template<typename uIntType>
  boost::python::tuple
  intersection_i_seqs(
    const_ref<uIntType> const& self,
    const_ref<uIntType> const& other)
  {
    intersection_with_tracking<uIntType, uIntType> proxy(
      self, other, /*track_matching_elements*/ false, /*track_i_seqs*/ true);
    return boost::python::make_tuple(proxy.self_i_seqs, proxy.other_i_seqs);
  }

} // namespace <anonymous>

  // wrap signed integer types
  // ==========================================================================
  #define WRAP_FLEX(pyname, uIntType) \
  namespace {\
  flex<uIntType>::type* \
  from_stl_vector_unsigned_to_##uIntType(std::vector<unsigned> const& v) \
  { \
    shared<uIntType> result(reserve(v.size())); \
    for(std::size_t i=0;i<v.size();i++) { \
      result.push_back(v[i]); \
    } \
    return new flex<uIntType>::type(result, result.size()); \
  }} \
  void wrap_flex_##uIntType() \
  { \
    using namespace boost::python; \
    using boost::python::arg; \
    flex_wrapper<uIntType>::integer(#pyname, scope()) \
      .def_pickle(flex_pickle_single_buffered<uIntType>()) \
      .def("__init__", make_constructor( \
        from_stl_vector_unsigned_to_##uIntType, default_call_policies())) \
      .def("__init__", make_constructor( \
        flex_##pyname##_from_numpy_array, default_call_policies())) \
      .def("copy_to_byte_str", \
        copy_to_byte_str<versa<uIntType, flex_grid<> > >) \
      .def("as_int", as_int<uIntType>) \
      .def("intersection", \
        (shared<uIntType>(*)(  \
          const_ref<uIntType> const&, \
          const_ref<uIntType> const&)) \
        intersection<uIntType>, (arg("self"), arg("other"))) \
      .def("intersection_i_seqs", intersection_i_seqs<uIntType>, ( \
        arg("self"), arg("other"))) \
      .def("counts", counts<uIntType, std::map<long, long> >::unlimited) \
      .def("counts", counts<uIntType, std::map<long, long> >::limited, ( \
        arg("max_keys"))) \
      .def("next_permutation", next_permutation<uIntType>) \
      .def("inverse_permutation", inverse_permutation<uIntType>) \
      .def("increment_and_track_up_from_zero", \
        increment_and_track_up_from_zero<uIntType>, (arg("iselection"))) \
      .def("as_numpy_array", flex_##pyname##_as_numpy_array, ( \
        arg("optional")=false)) \
    ; \
    def( \
      #pyname"_from_byte_str", \
      shared_from_byte_str<uIntType>, \
      (arg("byte_str"))); \
    range_wrappers<uIntType, int64_t, range_args::unsigned_check>::wrap( \
      #pyname"_range"); \
  }

  // --------------------------------------------------------------------------
  WRAP_FLEX(size_t, size_t);
  WRAP_FLEX(uint8, uint8_t);
  WRAP_FLEX(uint16, uint16_t);
  WRAP_FLEX(uint32, uint32_t);
  WRAP_FLEX(uint64, uint64_t);  // wrapped as size_t

}}} // namespace scitbx::af::boost_python
