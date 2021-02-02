#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <boost/python/module.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost_adaptbx/type_id_eq.h>
#include <vector>
#include <set>
#include <scitbx/mat3.h>
#include <scitbx/array_family/accessors/c_grid_padded.h>

namespace scitbx { namespace af { namespace boost_python {
namespace {

  void
  append_union_of_selected_arrays(
    af::shared<std::set<unsigned> >& self,
    af::const_ref<std::vector<unsigned> > const& arrays,
    af::const_ref<std::size_t> const& selection)
  {
    self.push_back(std::set<unsigned>());
    std::set<unsigned>& u = self.back();
    for(std::size_t i=0;i<selection.size();i++) {
      unsigned selection_i =  selection[i];
      SCITBX_ASSERT(selection_i < arrays.size());
      std::vector<unsigned> const& v = arrays[selection_i];
      u.insert(v.begin(), v.end());
    }
  }

  void init_module()
  {
    using namespace boost::python;
    using boost::python::arg;
    typedef return_internal_reference<> rir;
// #if !defined(BOOST_ADAPTBX_TYPE_ID_SIZE_T_EQ_UNSIGNED)
//     shared_wrapper<unsigned>::wrap("unsigned");
// #endif
    shared_wrapper<std::vector<unsigned>, rir>::wrap("stl_vector_unsigned");
    shared_wrapper<std::vector<double>, rir>::wrap("stl_vector_double");
    shared_wrapper<std::set<unsigned>, rir>::wrap("stl_set_unsigned")
      .def("append_union_of_selected_arrays",
        append_union_of_selected_arrays, (
          arg("arrays"), arg("selection")))
      .enable_pickling()
    ;
    shared_wrapper<mat3<int> >::wrap("mat3_int");
    shared_wrapper<tiny<int, 3> >::wrap("tiny_int_3");
    shared_wrapper<tiny<int, 4> >::wrap("tiny_int_4");
    shared_wrapper<tiny<int, 2> >::wrap("tiny_int_2");
      // used by scitbx.iso_surface
  }

}}}} // namespace scitbx::af::boost_python::<anonymous>

BOOST_PYTHON_MODULE(scitbx_array_family_shared_ext)
{
  scitbx::af::boost_python::init_module();
}
