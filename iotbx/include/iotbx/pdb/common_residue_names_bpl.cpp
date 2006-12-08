#include <iotbx/pdb/common_residue_names.h>
#include <scitbx/boost_python/array_as_list.h>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>

namespace iotbx { namespace pdb {
namespace {

  void
  wrap_common_residue_names_impl()
  {
    using namespace boost::python;

#define IOTBX_PDB_COMMON_RESIDUE_NAMES_LIST(type) \
    scope().attr("common_residue_names_" #type) \
      = scitbx::boost_python::null_terminated_array_as_list( \
          common_residue_names::type);

    IOTBX_PDB_COMMON_RESIDUE_NAMES_LIST(amino_acid)
    IOTBX_PDB_COMMON_RESIDUE_NAMES_LIST(rna_dna)
    IOTBX_PDB_COMMON_RESIDUE_NAMES_LIST(water)
    IOTBX_PDB_COMMON_RESIDUE_NAMES_LIST(small_molecule)
    IOTBX_PDB_COMMON_RESIDUE_NAMES_LIST(element)

    typedef return_value_policy<copy_const_reference> ccr;
    def("common_residue_names_get_class",
      (std::string const& (*)(std::string const& name))
        common_residue_names::get_class, (arg_("name")), ccr());
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_common_residue_names() { wrap_common_residue_names_impl(); }

}}} // namespace iotbx::pdb::boost_python
