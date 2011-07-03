#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/boost_python/flex_pickle_single_buffered.h>
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
    af::shared<std::size_t> result(af::reserve(v.size()));
    for(std::size_t i=0;i<v.size();i++) {
      result.push_back(v[i]);
    }
    return new flex<std::size_t>::type(result, result.size());
  }

  af::shared<int>
  as_int(
    af::const_ref<std::size_t> const& O)
  {
    af::shared<int> result(O.size(), af::init_functor_null<int>());
    for(std::size_t i=0;i<O.size();i++) {
      result[i] = static_cast<int>(O[i]);
    }
    return result;
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

  std::size_t
  increment_and_track_up_from_zero(
    af::ref<std::size_t> const& O,
    af::const_ref<std::size_t> const& iselection)
  {
    std::size_t result = 0;
    for(std::size_t i=0;i<iselection.size();i++) {
      std::size_t ii = iselection[i];
      SCITBX_ASSERT(ii < O.size());
      if (O[ii]++ == 0) result++;
    }
    return result;
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
      .def("as_int", as_int)
      .def("intersection",
        (af::shared<std::size_t>(*)(
          af::const_ref<std::size_t> const&,
          af::const_ref<std::size_t> const&))
        af::intersection, (arg("self"), arg("other")))
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
    range_wrappers<std::size_t, long, range_args::unsigned_check>::wrap(
      "size_t_range");
  }

}}} // namespace scitbx::af::boost_python
