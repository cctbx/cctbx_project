#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

#include <iotbx/pdb/utils.h>

namespace iotbx { namespace pdb { namespace boost_python {

  void wrap_hybrid_36();
  void wrap_common_residue_names();
  void wrap_rna_dna_atom_names();
  void wrap_input();
  void wrap_xray_structure();

namespace {

  void init_module()
  {
    using namespace boost::python;
    def("utils_base_256_ordinal", utils::base_256_ordinal, (arg_("s")));

    wrap_hybrid_36();
    wrap_common_residue_names();
    wrap_rna_dna_atom_names();
    wrap_input();
    wrap_xray_structure();
  }

}}}} // namespace iotbx::pdb::boost_python::<anonymous>

BOOST_PYTHON_MODULE(iotbx_pdb_ext)
{
  iotbx::pdb::boost_python::init_module();
}
