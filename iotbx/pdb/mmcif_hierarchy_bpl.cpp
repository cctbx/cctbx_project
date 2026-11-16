#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/return_value_policy.hpp>
#include <iotbx/pdb/mmcif_hierarchy.h>

namespace iotbx { namespace pdb { namespace boost_python {

  namespace {

    // build_hierarchy touches no Python object: the columns are af::shared
    // handles and the result is a C++ root, converted once the GIL is
    // back. So the build runs with the GIL released and another thread (a
    // GUI, say) keeps running; that thread must not modify the columns
    // meanwhile. The guard restores the GIL on the way out, on an
    // exception too, before Boost.Python translates it.
    struct gil_release
    {
      PyThreadState* save;
      gil_release() : save(PyEval_SaveThread()) {}
      ~gil_release() { PyEval_RestoreThread(save); }
    };

    hierarchy::root
    build_hierarchy_without_gil(mmcif::atom_site_columns const& columns)
    {
      gil_release guard;
      return mmcif::build_hierarchy(columns);
    }

  } // namespace <anonymous>

  void wrap_mmcif_hierarchy()
  {
    using namespace boost::python;
    typedef mmcif::atom_site_columns w_t;
    // The members are af::shared<std::string>; def_readwrite would hand out
    // an internal reference to a class that is not registered, so go through
    // the shared <-> flex.std_string converters by value (same storage).
    typedef return_value_policy<return_by_value> rbv;
#define IOTBX_PDB_MMCIF_COLUMN(name) \
      .add_property(#name, \
        make_getter(&w_t::name, rbv()), make_setter(&w_t::name))
    class_<w_t>("atom_site_columns")
      IOTBX_PDB_MMCIF_COLUMN(type_symbol)
      IOTBX_PDB_MMCIF_COLUMN(atom_name)
      IOTBX_PDB_MMCIF_COLUMN(alt_id)
      IOTBX_PDB_MMCIF_COLUMN(label_asym_id)
      IOTBX_PDB_MMCIF_COLUMN(auth_asym_id)
      IOTBX_PDB_MMCIF_COLUMN(comp_id)
      IOTBX_PDB_MMCIF_COLUMN(seq_id)
      IOTBX_PDB_MMCIF_COLUMN(cartn_x)
      IOTBX_PDB_MMCIF_COLUMN(cartn_y)
      IOTBX_PDB_MMCIF_COLUMN(cartn_z)
      IOTBX_PDB_MMCIF_COLUMN(b_iso)
      IOTBX_PDB_MMCIF_COLUMN(occupancy)
      IOTBX_PDB_MMCIF_COLUMN(model_num)
      IOTBX_PDB_MMCIF_COLUMN(atom_id)
      IOTBX_PDB_MMCIF_COLUMN(ins_code)
      IOTBX_PDB_MMCIF_COLUMN(group_pdb)
      IOTBX_PDB_MMCIF_COLUMN(segid)
      IOTBX_PDB_MMCIF_COLUMN(auth_break)
      IOTBX_PDB_MMCIF_COLUMN(formal_charge)
      IOTBX_PDB_MMCIF_COLUMN(resolution)
      IOTBX_PDB_MMCIF_COLUMN(fp)
      IOTBX_PDB_MMCIF_COLUMN(fdp)
      IOTBX_PDB_MMCIF_COLUMN(aniso_id)
      IOTBX_PDB_MMCIF_COLUMN(aniso_11)
      IOTBX_PDB_MMCIF_COLUMN(aniso_22)
      IOTBX_PDB_MMCIF_COLUMN(aniso_33)
      IOTBX_PDB_MMCIF_COLUMN(aniso_12)
      IOTBX_PDB_MMCIF_COLUMN(aniso_13)
      IOTBX_PDB_MMCIF_COLUMN(aniso_23)
      .def_readwrite("aniso_is_b", &w_t::aniso_is_b)
    ;
#undef IOTBX_PDB_MMCIF_COLUMN
    def("build_mmcif_hierarchy", build_hierarchy_without_gil,
      (arg("columns")));
    def("mmcif_parse_cell_double", mmcif::parse_cell_double, (arg("cell")));
    def("mmcif_parse_cell_int", mmcif::parse_cell_int, (arg("cell")));
  }

}}} // namespace iotbx::pdb::boost_python
