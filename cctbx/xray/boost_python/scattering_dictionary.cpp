#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/xray/scattering_dictionary.h>
#include <boost/python/class.hpp>
#include <boost/python/list.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/detail/api_placeholder.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_by_value.hpp>

namespace cctbx { namespace xray { namespace boost_python {

namespace {

  struct scatterer_group_wrappers : boost::python::pickle_suite
  {
    typedef scatterer_group w_t;

    static
    boost::python::tuple
    getstate(boost::python::object w_obj)
    {
      using namespace boost::python;
      w_t const& w = extract<w_t const&>(w_obj)();
      return make_tuple(
        w_obj.attr("__dict__"),
        w.coefficients,
        w.member_indices);
    }

    static
    void
    setstate(boost::python::object w_obj, boost::python::tuple state)
    {
      using namespace boost::python;
      w_t& w = extract<w_t&>(w_obj)();
      if (len(state) != 3) {
        PyErr_SetObject(PyExc_ValueError,
          ("expected 3-item tuple in call to __setstate__; got %s"
           % state).ptr());
        throw_error_already_set();
      }
      // restore the object's __dict__
      dict d = extract<dict>(w_obj.attr("__dict__"))();
      d.update(state[0]);
      // restore the internal state of the C++ object
      w.coefficients = extract<eltbx::caasf::custom const&>(state[1])();
      w.member_indices = extract<af::shared<std::size_t> >(state[2])();
    }

    static bool getstate_manages_dict() { return true; }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("scatterer_group")
        .def_readonly("coefficients", &w_t::coefficients)
        .add_property("member_indices",
          make_getter(&w_t::member_indices, rbv()))
        .def_pickle(scatterer_group_wrappers())
      ;
    }
  };

  struct scattering_dictionary_wrappers : boost::python::pickle_suite
  {
    typedef scattering_dictionary w_t;

    static boost::python::dict
    as_python_dict(w_t const& o)
    {
      boost::python::dict result;
      w_t::dict_type::const_iterator e;
      for(e=o.dict().begin();e!=o.dict().end();e++) {
        result[e->first] = e->second;
      }
      return result;
    }

    static std::size_t
    dict_size(w_t const& o) { return o.dict().size(); }

    static
    boost::python::tuple
    getstate(boost::python::object w_obj)
    {
      using namespace boost::python;
      w_t const& w = extract<w_t const&>(w_obj)();
      return make_tuple(
        w_obj.attr("__dict__"),
        as_python_dict(w));
    }

    static
    void
    setstate(boost::python::object w_obj, boost::python::tuple state)
    {
      using namespace boost::python;
      w_t& w = extract<w_t&>(w_obj)();
      if (len(state) != 2) {
        PyErr_SetObject(PyExc_ValueError,
          ("expected 2-item tuple in call to __setstate__; got %s"
           % state).ptr());
        throw_error_already_set();
      }
      // restore the object's __dict__
      dict d = extract<dict>(w_obj.attr("__dict__"))();
      d.update(state[0]);
      // restore the internal state of the C++ object
      list keys = extract<list>(state[1].attr("keys")())();
      list values = extract<list>(state[1].attr("values")())();
      CCTBX_ASSERT(len(keys) == len(values));
      for(std::size_t i=0;i<len(keys);i++) {
        std::string label = extract<std::string>(keys[i])();
        scatterer_group group = extract<scatterer_group>(values[i])();
        w.dict()[label] = group;
      }
    }

    static bool getstate_manages_dict() { return true; }

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("scattering_dictionary")
        .def(init<af::const_ref<scatterer<> > const&>())
        .def("n_scatterers", &w_t::n_scatterers)
        .def("dict", as_python_dict)
        .def("dict_size", dict_size)
        .def("find_undefined", &w_t::find_undefined)
        .def("lookup", &w_t::lookup, ccr())
        .def("scatterer_permutation", &w_t::scatterer_permutation)
        .def("assign", &w_t::assign)
        .def("assign_from_table", &w_t::assign_from_table)
        .def_pickle(scattering_dictionary_wrappers())
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
