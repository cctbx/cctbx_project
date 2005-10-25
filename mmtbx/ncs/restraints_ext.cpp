#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <mmtbx/ncs/restraints.h>
#include <scitbx/boost_python/container_conversions.h>

namespace mmtbx { namespace ncs { namespace restraints {
namespace {

  void
  register_to_tuple()
  {
    using namespace scitbx::boost_python::container_conversions;
    to_tuple_mapping<af::tiny<af::shared<std::size_t>, 2> >();
    to_tuple_mapping<std::vector<af::tiny<af::shared<std::size_t>, 2> > >();
  }

  struct pair_registry_wrappers
  {
    typedef pair_registry w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("pair_registry", no_init)
        .def(init<unsigned, unsigned>((arg_("n_seq"), arg_("n_ncs"))))
        .def("enter", &w_t::enter, (
          (arg_("i_seq"), arg_("j_seq"), arg_("j_ncs"))))
        .def("selection_pairs", &w_t::selection_pairs)
        .def("adp_iso_residual_sum", &w_t::adp_iso_residual_sum, (
          (arg_("weight"),
           arg_("average_power"),
           arg_("u_isos"),
           arg_("u_average_min"),
           arg_("gradients"))))
      ;
    }
  };

  void init_module()
  {
    register_to_tuple();
    pair_registry_wrappers::wrap();
  }

} // namespace <anonymous>
}}} // namespace mmtbx::ncs::restraints

BOOST_PYTHON_MODULE(mmtbx_ncs_restraints_ext)
{
  mmtbx::ncs::restraints::init_module();
}
