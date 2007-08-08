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
        .def("flags_as_string", &w_t::flags_as_string)
      ;
    }
  };

  void
  wrap_rna_dna_atom_names()
  {
    rna_dna_atom_names_info_wrappers::wrap();
  }

}}} // namespace iotbx::pdb::boost_python
