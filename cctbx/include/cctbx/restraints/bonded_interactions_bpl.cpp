#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>

#include <cctbx/restraints/bonded_interactions.h>

namespace cctbx { namespace restraints {
namespace {

  struct bonded_interactions_wrappers
  {
    typedef bonded_interactions w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("bonded_interactions", no_init)
        .def(init<af::const_ref<std::set<unsigned> > const&,
                  unsigned,
                  optional<bool> >(
          (arg_("bond_sets"),
           arg_("i_seq_0"),
           arg_("eliminate_redundant_interactions"))))
        .def("i_seq_0", &w_t::i_seq_0)
        .def("eliminate_redundant_interactions",
          &w_t::eliminate_redundant_interactions)
        .def("interactions_1_2", &w_t::interactions_1_2, ccr())
        .def("interactions_1_3", &w_t::interactions_1_3, ccr())
        .def("interactions_1_4", &w_t::interactions_1_4, ccr())
        .def("interaction_type_of", &w_t::interaction_type_of,
          (arg_("j_seq")))
      ;
    }
  };

} // namespace <anonymous>

namespace boost_python {

  void wrap_bonded_interactions()
  {
    bonded_interactions_wrappers::wrap();
  }

}}} // namespace cctbx::restraints::boost_python
