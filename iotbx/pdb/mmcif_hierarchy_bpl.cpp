#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/return_value_policy.hpp>
#include <iotbx/pdb/mmcif_hierarchy.h>

namespace iotbx { namespace pdb { namespace boost_python {

  void wrap_mmcif_hierarchy()
  {
    using namespace boost::python;
    typedef mmcif::atom_site_columns w_t;
    // The members are af::shared<std::string>; def_readwrite would hand out
    // an internal reference to a class that is not registered, so go through
    // the shared <-> flex.std_string converters by value (same storage).
    typedef return_value_policy<return_by_value> rbv;
    class_<w_t>("atom_site_columns")
      .add_property("type_symbol",
        make_getter(&w_t::type_symbol, rbv()),
        make_setter(&w_t::type_symbol))
      .add_property("atom_name",
        make_getter(&w_t::atom_name, rbv()),
        make_setter(&w_t::atom_name))
      .add_property("alt_id",
        make_getter(&w_t::alt_id, rbv()),
        make_setter(&w_t::alt_id))
      .add_property("label_asym_id",
        make_getter(&w_t::label_asym_id, rbv()),
        make_setter(&w_t::label_asym_id))
      .add_property("auth_asym_id",
        make_getter(&w_t::auth_asym_id, rbv()),
        make_setter(&w_t::auth_asym_id))
      .add_property("comp_id",
        make_getter(&w_t::comp_id, rbv()),
        make_setter(&w_t::comp_id))
      .add_property("seq_id",
        make_getter(&w_t::seq_id, rbv()),
        make_setter(&w_t::seq_id))
      .add_property("cartn_x",
        make_getter(&w_t::cartn_x, rbv()),
        make_setter(&w_t::cartn_x))
      .add_property("cartn_y",
        make_getter(&w_t::cartn_y, rbv()),
        make_setter(&w_t::cartn_y))
      .add_property("cartn_z",
        make_getter(&w_t::cartn_z, rbv()),
        make_setter(&w_t::cartn_z))
      .add_property("b_iso",
        make_getter(&w_t::b_iso, rbv()),
        make_setter(&w_t::b_iso))
      .add_property("occupancy",
        make_getter(&w_t::occupancy, rbv()),
        make_setter(&w_t::occupancy))
      .add_property("model_num",
        make_getter(&w_t::model_num, rbv()),
        make_setter(&w_t::model_num))
      .add_property("atom_id",
        make_getter(&w_t::atom_id, rbv()),
        make_setter(&w_t::atom_id))
      .add_property("ins_code",
        make_getter(&w_t::ins_code, rbv()),
        make_setter(&w_t::ins_code))
      .add_property("group_pdb",
        make_getter(&w_t::group_pdb, rbv()),
        make_setter(&w_t::group_pdb))
      .add_property("segid",
        make_getter(&w_t::segid, rbv()),
        make_setter(&w_t::segid))
      .add_property("auth_break",
        make_getter(&w_t::auth_break, rbv()),
        make_setter(&w_t::auth_break))
      .add_property("formal_charge",
        make_getter(&w_t::formal_charge, rbv()),
        make_setter(&w_t::formal_charge))
      .add_property("resolution",
        make_getter(&w_t::resolution, rbv()),
        make_setter(&w_t::resolution))
      .add_property("fp",
        make_getter(&w_t::fp, rbv()),
        make_setter(&w_t::fp))
      .add_property("fdp",
        make_getter(&w_t::fdp, rbv()),
        make_setter(&w_t::fdp))
      .add_property("aniso_id",
        make_getter(&w_t::aniso_id, rbv()),
        make_setter(&w_t::aniso_id))
      .add_property("aniso_11",
        make_getter(&w_t::aniso_11, rbv()),
        make_setter(&w_t::aniso_11))
      .add_property("aniso_22",
        make_getter(&w_t::aniso_22, rbv()),
        make_setter(&w_t::aniso_22))
      .add_property("aniso_33",
        make_getter(&w_t::aniso_33, rbv()),
        make_setter(&w_t::aniso_33))
      .add_property("aniso_12",
        make_getter(&w_t::aniso_12, rbv()),
        make_setter(&w_t::aniso_12))
      .add_property("aniso_13",
        make_getter(&w_t::aniso_13, rbv()),
        make_setter(&w_t::aniso_13))
      .add_property("aniso_23",
        make_getter(&w_t::aniso_23, rbv()),
        make_setter(&w_t::aniso_23))
      .def_readwrite("aniso_is_b", &w_t::aniso_is_b)
    ;
    def("build_mmcif_hierarchy", mmcif::build_hierarchy, (arg("columns")));
    def("mmcif_parse_cell_double", mmcif::parse_cell_double, (arg("cell")));
    def("mmcif_parse_cell_int", mmcif::parse_cell_int, (arg("cell")));
  }

}}} // namespace iotbx::pdb::boost_python
