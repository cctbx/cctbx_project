#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/return_value_policy.hpp>
#include <iotbx/pdb/mmcif_hierarchy.h>
#include <cctbx/xray/scatterer.h>
#include <cctbx/eltbx/xray_scattering.h>
#include <cctbx/adptbx.h>
#include <cctype>
#include <map>
#include <set>
#include <stdexcept>

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

    //! The scattering type of a hierarchy atom with these fields, as
    //! xray_structures_simple_extension determines it.
    std::string
    scattering_type_for(
      std::string const& name,
      std::string const& element,
      std::string const& charge,
      std::string const& label,
      bool scattering_type_exact,
      bool enable_scattering_type_unknown)
    {
      hierarchy::atom atom;
      atom.set_name(name.c_str());
      atom.set_element(element.c_str());
      atom.set_charge(charge.c_str());
      // the label, element and charge of the row, for the messages
      std::string const where = label + " type_symbol \"" + element
        + "\" charge \"" + charge + "\"";
      if (atom.data->element == "IS") {
        std::string ias_type;
        for(std::size_t i=0;i<2;i++) {
          if (isdigit(atom.data->charge.elems[i])) {
            ias_type += atom.data->charge.elems[i];
          }
        }
        if (ias_type.size() == 0) {
          throw std::invalid_argument(
            "Unknown chemical element type:\n  " + where + "\n");
        }
        return "IS" + ias_type;
      }
      boost::optional<std::string>
        chemical_element = atom.determine_chemical_element_simple();
      if (!chemical_element && !enable_scattering_type_unknown) {
        throw std::invalid_argument(
          "Unknown chemical element type:\n  " + where + "\n");
      }
      boost::optional<std::string> charge_tidy = atom.charge_tidy();
      if (!charge_tidy) {
        if (!enable_scattering_type_unknown) {
          throw std::invalid_argument("Unknown charge:\n  " + where);
        }
        chemical_element.reset();
      }
      boost::optional<std::string> scattering_type;
      if (chemical_element) {
        scattering_type = cctbx::eltbx::xray_scattering::get_standard_label(
          *chemical_element + *charge_tidy, scattering_type_exact, true);
      }
      if (scattering_type) return *scattering_type;
      if (enable_scattering_type_unknown) return "unknown";
      throw std::invalid_argument("Unknown scattering type:\n  " + where);
    }

    //! The scatterers of an _atom_site loop, in row order: the scatterers
    //! of xray_structures_simple_extension for the hierarchy of the loop.
    af::shared<cctbx::xray::scatterer<> >
    build_scatterers(
      mmcif::atom_site_columns const& columns,
      cctbx::uctbx::unit_cell const& unit_cell,
      std::set<std::string> const& atom_names_scattering_type_const,
      bool scattering_type_exact,
      bool enable_scattering_type_unknown)
    {
      // touches no Python object, as build_hierarchy_without_gil
      gil_release guard;
      mmcif::atom_site_rows rows = mmcif::parse_atom_site_rows(columns);
      const std::size_t n = rows.label.size();
      af::shared<cctbx::xray::scatterer<> > result;
      result.reserve(n);
      // the element, charge and, for a blank element, the name decide the
      // scattering type
      std::map<std::string, std::string> types;
      cctbx::xray::scatterer<> scatterer(
        "", cctbx::fractional<>(0,0,0), 0, 0, "", 0, 0);
      for(std::size_t i=0;i<n;i++) {
        scatterer.label = rows.label[i];
        scatterer.site = unit_cell.fractionalize(rows.sites_cart[i]);
        scitbx::sym_mat3<double> const& uij = rows.u_cart[i];
        if (uij.const_ref().all_eq(-1)) {
          scatterer.u_iso = cctbx::adptbx::b_as_u(rows.b_iso[i]);
          scatterer.set_use_u_iso_only();
        }
        else {
          scatterer.u_star = cctbx::adptbx::u_cart_as_u_star(unit_cell, uij);
          scatterer.set_use_u_aniso_only();
        }
        scatterer.fp = rows.fp[i];
        scatterer.fdp = rows.fdp[i];
        scatterer.occupancy = rows.occupancy[i];
        if (   atom_names_scattering_type_const.find(rows.name[i])
            != atom_names_scattering_type_const.end()) {
          scatterer.scattering_type = "const";
        }
        else {
          std::string key = rows.element[i] + "|" + rows.charge[i];
          if (rows.element[i].find_first_not_of(' ') == std::string::npos) {
            key += "|" + rows.name[i];
          }
          std::map<std::string, std::string>::const_iterator t
            = types.find(key);
          if (t == types.end()) {
            t = types.insert(std::make_pair(key, scattering_type_for(
              rows.name[i], rows.element[i], rows.charge[i], rows.label[i],
              scattering_type_exact, enable_scattering_type_unknown))).first;
          }
          scatterer.scattering_type = t->second;
        }
        result.push_back(scatterer);
      }
      return result;
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
    def("build_mmcif_scatterers", build_scatterers, (
      arg("columns"),
      arg("unit_cell"),
      arg("atom_names_scattering_type_const"),
      arg("scattering_type_exact"),
      arg("enable_scattering_type_unknown")));
    def("mmcif_parse_cell_double", mmcif::parse_cell_double, (arg("cell")));
    def("mmcif_parse_cell_int", mmcif::parse_cell_int, (arg("cell")));
  }

}}} // namespace iotbx::pdb::boost_python
