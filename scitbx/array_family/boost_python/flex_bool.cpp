#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/boost_python/flex_pickle_single_buffered.h>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>

namespace scitbx { namespace af { namespace boost_python {

namespace {

  af::shared<std::size_t>
  iselection(
    af::const_ref<bool, flex_grid<> > const& a,
    bool test_value=true)
  {
    af::shared<std::size_t> result;
    for(std::size_t i=0;i<a.size();i++) {
      if (a[i] == test_value) result.push_back(i);
    }
    return result;
  }

  BOOST_PYTHON_FUNCTION_OVERLOADS(iselection_overloads, iselection, 1, 2)

} // namespace <anonymous>

  void wrap_flex_bool()
  {
    typedef boost::python::arg arg_; // gcc 2.96 workaround
    flex_wrapper<bool>::logical("bool", boost::python::scope())
      .def_pickle(flex_pickle_single_buffered<bool>())
      .def("iselection", iselection,
        iselection_overloads((arg_("self"), arg_("test_value")=true)))
    ;
  }

}}} // namespace scitbx::af::boost_python
