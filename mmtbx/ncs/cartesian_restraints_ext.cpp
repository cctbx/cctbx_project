#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/register_ptr_to_python.hpp>
#include <mmtbx/ncs/cartesian_restraints.h>
#include <scitbx/boost_python/container_conversions.h>

namespace mmtbx { namespace ncs { namespace cartesian_restraints {
namespace {

  void
  register_to_tuple()
  {
    using namespace scitbx::boost_python::container_conversions;
    to_tuple_mapping<af::tiny<af::shared<std::size_t>, 2> >();
    to_tuple_mapping<std::vector<af::tiny<af::shared<std::size_t>, 2> > >();
  }

  struct pair_registry_wrappers : boost::python::pickle_suite
  {
    typedef pair_registry w_t;

    static boost::python::tuple
      getinitargs(w_t const& self)
    {
      return boost::python::make_tuple(self.n_seq(),
        self.init_n_ncs
        );
    }

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t, boost::noncopyable>("pair_registry", no_init)
        .def(init<unsigned, unsigned>((arg("n_seq"), arg("n_ncs"))))
        .def("n_seq", &w_t::n_seq)
        .def_readonly("number_of_additional_isolated_sites",
          &w_t::number_of_additional_isolated_sites)
        .def("register_additional_isolated_sites",
          &w_t::register_additional_isolated_sites, (arg("number")))
        .def("enter", &w_t::enter, (
          arg("i_seq"), arg("j_seq"), arg("j_ncs")))
        .def("proxy_select", &w_t::proxy_select, (
          arg("iselection")))
        .def("selection_pairs", &w_t::selection_pairs)
        .def("adp_iso_residual_sum", &w_t::adp_iso_residual_sum, (
          arg("weight"),
          arg("average_power"),
          arg("u_isos"),
          arg("u_average_min"),
          arg("gradients")))
        .def_pickle(pair_registry_wrappers())
      ;
      register_ptr_to_python<boost::shared_ptr<pair_registry> >();
    }
  };

  void init_module()
  {
    register_to_tuple();
    pair_registry_wrappers::wrap();
  }

} // namespace <anonymous>
}}} // namespace mmtbx::ncs::cartesian_restraints

BOOST_PYTHON_MODULE(mmtbx_ncs_cartesian_restraints_ext)
{
  mmtbx::ncs::cartesian_restraints::init_module();
}
