#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/xray/scattering_dictionary.h>
#include <boost/python/class.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

namespace cctbx { namespace xray { namespace boost_python {

namespace {

  struct scatterer_group_wrappers
  {
    typedef scatterer_group w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("scatterer_group", no_init)
        .def_readonly("coefficients", &w_t::coefficients)
        .add_property("member_indices",
          make_getter(&w_t::member_indices, rbv()))
      ;
    }
  };

  struct scattering_dictionary_wrappers
  {
    typedef scattering_dictionary w_t;

    static boost::python::dict
    dict(w_t const& o)
    {
      boost::python::dict result;
      w_t::dict_type::const_iterator e;
      for(e=o.dict().begin();e!=o.dict().end();e++) {
        result[e->first] = e->second;
      }
      return result;
    }

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("scattering_dictionary", no_init)
        .def(init<af::const_ref<scatterer<> > const&>())
        .def("size", &w_t::size)
        .def("dict", dict)
        .def("assign", &w_t::assign)
        .def("assign_from_table", &w_t::assign_from_table)
      ;
    }
  };

} // namespace <anoymous>

  void wrap_scattering_dictionary()
  {
    scatterer_group_wrappers::wrap();
    scattering_dictionary_wrappers::wrap();
  }

}}} // namespace cctbx::xray::boost_python
