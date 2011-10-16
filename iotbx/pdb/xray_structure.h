#ifndef IOTBX_PDB_XRAY_STRUCTURE_H
#define IOTBX_PDB_XRAY_STRUCTURE_H

#include <iotbx/pdb/input.h>
#include <iotbx/error.h>
#include <cctbx/xray/scatterer.h>
#include <cctbx/eltbx/xray_scattering.h>

namespace iotbx { namespace pdb {

  template <typename XrayScattererType=cctbx::xray::scatterer<> >
  struct xray_structures_simple_extension
  {
    protected:
      boost::shared_ptr<input> self_;
      af::shared<hierarchy::atom_with_labels> atoms_with_labels_;
      bool unit_cube_pseudo_crystal_;
      bool fractional_coordinates_;
      bool scattering_type_exact_;
      bool enable_scattering_type_unknown_;
      std::set<std::string> atom_names_scattering_type_const_;
      cctbx::uctbx::unit_cell unit_cell_;
      scitbx::mat3<double> scale_r_;
      scitbx::vec3<double> scale_t_;

    public:
      xray_structures_simple_extension() {}

      xray_structures_simple_extension(
        boost::shared_ptr<input> const& self,
        bool one_structure_for_each_model,
        bool unit_cube_pseudo_crystal,
        bool fractional_coordinates,
        bool scattering_type_exact,
        bool enable_scattering_type_unknown,
        std::set<std::string> const& atom_names_scattering_type_const,
        cctbx::uctbx::unit_cell const& unit_cell,
        scitbx::mat3<double> const& scale_r,
        scitbx::vec3<double> const& scale_t)
      :
        self_(self),
        atoms_with_labels_(self_->atoms_with_labels()),
        unit_cube_pseudo_crystal_(unit_cube_pseudo_crystal),
        fractional_coordinates_(fractional_coordinates),
        scattering_type_exact_(scattering_type_exact),
        enable_scattering_type_unknown_(enable_scattering_type_unknown),
        atom_names_scattering_type_const_(atom_names_scattering_type_const),
        unit_cell_(unit_cell),
        scale_r_(scale_r),
        scale_t_(scale_t),
        //
        loop_state(0),
        use_scale_matrix(scale_r.determinant() != 0),
        model_range(self_->model_indices().const_ref()),
        i_atom(0),
        atom(atoms_with_labels_.begin()),
        scatterer("", cctbx::fractional<>(0,0,0), 0, 0, "", 0, 0)
      {
        IOTBX_ASSERT(!use_scale_matrix || !fractional_coordinates);
        if (!one_structure_for_each_model) model_range.skip_to_last();
      }

      af::shared<XrayScattererType> scatterers;

    protected:
      unsigned loop_state;
      bool use_scale_matrix;
      range_loop<std::size_t> model_range;
      std::size_t i_atom;
      const hierarchy::atom_with_labels *atom;
      XrayScattererType scatterer;
      boost::optional<std::string> scattering_type;

    public:
      bool
      next()
      {
        if (loop_state == 1) goto continue_after_return;
        while (model_range.next()) {
          scatterers = af::shared<XrayScattererType>();
          scatterers.reserve(model_range.size);
          for(;i_atom!=model_range.end;i_atom++,atom++) {
            scatterer.label = atom->id_str();
            if (unit_cube_pseudo_crystal_ || fractional_coordinates_) {
              scatterer.site = atom->data->xyz;
            }
            else if (!use_scale_matrix) {
              scatterer.site = unit_cell_.fractionalize(atom->data->xyz);
            }
            else {
              scatterer.site = scale_r_ * atom->data->xyz + scale_t_;
            }
            if (atom->data->uij.const_ref().all_eq(-1)) {
              scatterer.u_iso = cctbx::adptbx::b_as_u(atom->data->b);
              scatterer.set_use_u_iso_only();
            }
            else {
              if (unit_cube_pseudo_crystal_) {
                scatterer.u_star = atom->data->uij;
              }
              else {
                scatterer.u_star = cctbx::adptbx::u_cart_as_u_star(
                  unit_cell_, atom->data->uij);
              }
              scatterer.set_use_u_aniso_only();
            }
            scatterer.occupancy = atom->data->occ;
            if (   atom_names_scattering_type_const_.find(
                     atom->data->name.elems)
                != atom_names_scattering_type_const_.end()) {
              scatterer.scattering_type = "const";
            }
            else if (atom->data->element == "IS") {
              char charge[2];
              std::memcpy(charge, atom->data->charge.elems, 2);
              std::string ias_type="";
              for(std::size_t i=0; i<2; i++) {
                if(isdigit(charge[i])) {
                  ias_type += charge[i];
                }
              }
              if(ias_type.size()>0) {
                 scatterer.scattering_type = std::string("IS")+ias_type;
              }
              else {
                throw std::invalid_argument(
                  std::string("Unknown chemical element type:\n")
                  + "  " + atom->quote() + "\n");
              }
            }
            else {
              boost::optional<std::string>
                chemical_element = atom->determine_chemical_element_simple();
              if (!chemical_element && !enable_scattering_type_unknown_) {
                throw std::invalid_argument(
                  std::string("Unknown chemical element type:\n")
                  + "  " + atom->quote() + "\n"
                  + "  To resolve this problem, specify a"
                    + " chemical element type in\n"
                  + "  columns 77-78 of the PDB file, right justified"
                    + " (e.g. \" C\").");
              }
              boost::optional<std::string>
                charge_tidy = atom->charge_tidy();
              if (!charge_tidy) {
                if (!enable_scattering_type_unknown_) {
                  throw std::invalid_argument(
                    std::string("Unknown charge:\n")
                    + "  " + atom->quote() + "\n"
                    + "                                       ^^");
                }
                chemical_element.reset();
              }
              if (chemical_element) {
                scattering_type = cctbx::eltbx::xray_scattering::
                  get_standard_label(
                    /*label*/ *chemical_element + *charge_tidy,
                    /*exact*/ scattering_type_exact_,
                    /*optional*/ true);
              }
              else {
                scattering_type.reset();
              }
              if (scattering_type) {
                scatterer.scattering_type = *scattering_type;
              }
              else if (enable_scattering_type_unknown_) {
                scatterer.scattering_type = "unknown";
              }
              else {
                throw std::invalid_argument(
                  std::string("Unknown scattering type:\n")
                  + "  " + atom->quote() + "\n"
                  + "               ^^^^                  ^^^^");
              }
            }
            scatterers.push_back(scatterer);
          }
          loop_state = 1;
          return true;
          continue_after_return:;
        }
        if (loop_state == 0) {
          loop_state = 2;
          return true;
        }
        return false;
      }
  };

}} // namespace iotbx::pdb

#endif // IOTBX_PDB_XRAY_STRUCTURE_H
