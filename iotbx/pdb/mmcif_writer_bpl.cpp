#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/return_value_policy.hpp>
#include <iotbx/pdb/mmcif_writer.h>

namespace iotbx { namespace pdb { namespace boost_python {

  void wrap_mmcif_writer()
  {
    using namespace boost::python;
    // The members are af::shared<std::string>, passed by value through the
    // shared <-> flex.std_string converters (same storage), as
    // atom_site_columns in mmcif_hierarchy_bpl.cpp.
    typedef return_value_policy<return_by_value> rbv;
#define IOTBX_PDB_MMCIF_COLUMN(name) \
      .add_property(#name, \
        make_getter(&w_t::name, rbv()), make_setter(&w_t::name))
    {
      typedef mmcif::atom_site_labels w_t;
      class_<w_t>("atom_site_labels")
        IOTBX_PDB_MMCIF_COLUMN(chain_auth_asym_id)
        IOTBX_PDB_MMCIF_COLUMN(rg_label_asym_id)
        IOTBX_PDB_MMCIF_COLUMN(rg_auth_seq_id)
        IOTBX_PDB_MMCIF_COLUMN(ag_label_seq_id)
      ;
    }
    {
      typedef mmcif::atom_site_output w_t;
      class_<w_t>("atom_site_output")
        IOTBX_PDB_MMCIF_COLUMN(group_pdb)
        IOTBX_PDB_MMCIF_COLUMN(id)
        IOTBX_PDB_MMCIF_COLUMN(label_atom_id)
        IOTBX_PDB_MMCIF_COLUMN(auth_atom_id)
        IOTBX_PDB_MMCIF_COLUMN(label_alt_id)
        IOTBX_PDB_MMCIF_COLUMN(label_comp_id)
        IOTBX_PDB_MMCIF_COLUMN(auth_asym_id)
        IOTBX_PDB_MMCIF_COLUMN(auth_seq_id)
        IOTBX_PDB_MMCIF_COLUMN(ins_code)
        IOTBX_PDB_MMCIF_COLUMN(cartn_x)
        IOTBX_PDB_MMCIF_COLUMN(cartn_y)
        IOTBX_PDB_MMCIF_COLUMN(cartn_z)
        IOTBX_PDB_MMCIF_COLUMN(occupancy)
        IOTBX_PDB_MMCIF_COLUMN(b_iso)
        IOTBX_PDB_MMCIF_COLUMN(type_symbol)
        IOTBX_PDB_MMCIF_COLUMN(formal_charge)
        IOTBX_PDB_MMCIF_COLUMN(fp)
        IOTBX_PDB_MMCIF_COLUMN(fdp)
        IOTBX_PDB_MMCIF_COLUMN(resolution)
        IOTBX_PDB_MMCIF_COLUMN(label_asym_id)
        IOTBX_PDB_MMCIF_COLUMN(label_entity_id)
        IOTBX_PDB_MMCIF_COLUMN(label_seq_id)
        IOTBX_PDB_MMCIF_COLUMN(model_num)
        IOTBX_PDB_MMCIF_COLUMN(auth_segid)
        IOTBX_PDB_MMCIF_COLUMN(auth_break)
        IOTBX_PDB_MMCIF_COLUMN(aniso_id)
        IOTBX_PDB_MMCIF_COLUMN(aniso_auth_atom_id)
        IOTBX_PDB_MMCIF_COLUMN(aniso_label_alt_id)
        IOTBX_PDB_MMCIF_COLUMN(aniso_auth_comp_id)
        IOTBX_PDB_MMCIF_COLUMN(aniso_auth_asym_id)
        IOTBX_PDB_MMCIF_COLUMN(aniso_auth_seq_id)
        IOTBX_PDB_MMCIF_COLUMN(aniso_ins_code)
        IOTBX_PDB_MMCIF_COLUMN(aniso_u11)
        IOTBX_PDB_MMCIF_COLUMN(aniso_u22)
        IOTBX_PDB_MMCIF_COLUMN(aniso_u33)
        IOTBX_PDB_MMCIF_COLUMN(aniso_u12)
        IOTBX_PDB_MMCIF_COLUMN(aniso_u13)
        IOTBX_PDB_MMCIF_COLUMN(aniso_u23)
        IOTBX_PDB_MMCIF_COLUMN(chem_comp_ids)
        IOTBX_PDB_MMCIF_COLUMN(struct_asym_ids)
      ;
    }
#undef IOTBX_PDB_MMCIF_COLUMN
    def("atom_site_columns_from_hierarchy",
      mmcif::atom_site_columns_from_hierarchy, (
        arg("root"),
        arg("labels"),
        arg("coordinate_precision"),
        arg("occupancy_precision"),
        arg("b_iso_precision"),
        arg("u_aniso_precision"),
        arg("resolution_precision"),
        arg("segid_as_auth_segid"),
        arg("output_break_records")));
  }

}}} // namespace iotbx::pdb::boost_python
