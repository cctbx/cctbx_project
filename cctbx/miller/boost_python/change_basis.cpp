#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <cctbx/miller/change_basis.h>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

namespace cctbx { namespace miller { namespace boost_python {

namespace {

  template <typename DataType, typename ChangeBasisPolicy>
  struct change_basis_wrappers
  {
    typedef change_basis<DataType, ChangeBasisPolicy> w_t;

    static void
    wrap(const char* python_name)
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>(python_name, no_init)
        .def(init<sgtbx::change_of_basis_op const&,
                  af::const_ref<index<> > const&,
                  af::const_ref<DataType> const&,
                  optional<bool> >((
          arg_("cb_op"),
          arg_("indices_in"),
          arg_("data_in"),
          arg_("deg")=false)))
        .add_property("indices", make_getter(&w_t::indices, rbv()))
        .add_property("data", make_getter(&w_t::data, rbv()))
      ;
    }
  };

}

  void wrap_change_basis()
  {
    change_basis_wrappers<
      double,
      change_basis_phase_policy<double> >::wrap(
        "change_basis_phases_double");
    change_basis_wrappers<
      std::complex<double>,
      change_basis_complex_policy<double> >::wrap(
        "change_basis_complex_double");
    change_basis_wrappers<
      hendrickson_lattman<>,
      change_basis_hendrickson_lattman_policy<double> >::wrap(
        "change_basis_hendrickson_lattman");
  }

}}} // namespace cctbx::miller::boost_python
