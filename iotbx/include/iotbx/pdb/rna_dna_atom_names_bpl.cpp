#include <boost/python/class.hpp>
#include <boost/python/args.hpp>

#include <iotbx/pdb/rna_dna_atom_names.h>

namespace iotbx { namespace pdb { namespace boost_python {

  struct rna_dna_atom_names_info_wrappers
  {
    typedef rna_dna_atom_names::info w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("rna_dna_atom_names_info", no_init)
        .def(init<const char*>((arg_("work_name"))))
        .def_readonly("reference_name", &w_t::reference_name)
        .def("compatible_residue_names", &w_t::compatible_residue_names)
        .def("is_compatible_with", &w_t::is_compatible_with, (
          arg_("residue_name")))
        .def("is_hydrogen", &w_t::is_hydrogen)
        .def("is_deuterium", &w_t::is_deuterium)
        .def("is_o2prime", &w_t::is_o2prime)
        .def("is_ho2prime", &w_t::is_ho2prime)
        .def("is_h2primeprime", &w_t::is_h2primeprime)
        .def("is_in_phosphate_group", &w_t::is_in_phosphate_group)
        .def("is_op3_or_hop3", &w_t::is_op3_or_hop3)
        .def("is_ho5prime", &w_t::is_ho5prime)
        .def("is_ho3prime", &w_t::is_ho3prime)
        .def("change_ho5prime_to_hop3", &w_t::change_ho5prime_to_hop3)
        .def("change_to_unknown", &w_t::change_to_unknown)
      ;
    }
  };

  void
  wrap_rna_dna_atom_names()
  {
    rna_dna_atom_names_info_wrappers::wrap();
  }

}}} // namespace iotbx::pdb::boost_python
