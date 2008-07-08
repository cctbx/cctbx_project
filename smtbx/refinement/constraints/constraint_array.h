#ifndef SMTBX_REFINEMENT_CONSTRAINTS_CONSTRAINT_ARRAY_H
#define SMTBX_REFINEMENT_CONSTRAINTS_CONSTRAINT_ARRAY_H

#include <scitbx/array_family/shared.h>

#include <cctbx/uctbx.h>
#include <cctbx/sgtbx/site_symmetry_table.h>
#include <cctbx/xray/scatterer_flags.h>

#include <smtbx/error.h>
#include <smtbx/refinement/parameter_map.h>
#include <smtbx/import_scitbx_af.h>
#include <smtbx/import_cctbx.h>

namespace smtbx { namespace refinement { namespace constraints {

template<class ConstraintType, template<class> class SharedArray1D>
class constraint_array
{
  public:
    typedef ConstraintType constraint_type;
    typedef typename ConstraintType::xray_scatterer_type xray_scatterer_type;
    typedef typename ConstraintType::float_type float_type;
    typedef parameter_map<xray_scatterer_type> parameter_map_type;

    constraint_array(uctbx::unit_cell const &unit_cell_,
                     sgtbx::site_symmetry_table const &site_symmetry_table_,
                     af::shared<xray_scatterer_type> scatterers_,
                     parameter_map_type const &crystallographic_parameter_map_,
                     af::shared<xray::scatterer_flags> constraint_flags_)
      : unit_cell(unit_cell_),
        site_symmetry_table(site_symmetry_table_),
        scatterers(scatterers_),
        crystallographic_parameter_map(crystallographic_parameter_map_),
        constraint_flags(constraint_flags_),
        elements()
    {}

    constraint_type & operator[](unsigned i) { return elements[i]; }

    unsigned size() { return elements.size(); }

    void push_back(constraint_type &constraint) {
      constraint.initialise_in_context(unit_cell,
                                       scatterers.const_ref(),
                                       constraint_flags.ref(),
                                       already_constrained_);
      elements.push_back(constraint);
    }

    void erase(constraint_type const *e) {
      elements.erase(const_cast<constraint_type *>(e));
    }

    void compute_gradients(
      af::ref<float_type> const &crystallographic_gradients,
      SharedArray1D<float_type> reparametrization_gradients)
    {
      for(int i=0; i < size(); ++i) {
        (*this)[i].compute_gradients(unit_cell,
                                     site_symmetry_table,
                                     scatterers.const_ref(),
                                     crystallographic_parameter_map,
                                     crystallographic_gradients,
                                     reparametrization_gradients);
      }
    }

    void apply_shifts(
      af::const_ref<float_type> const &crystallographic_shifts,
      af::const_ref<float_type> const &reparametrization_shifts)
    {
      for(int i=0; i < size(); ++i) {
        (*this)[i].apply_shifts(unit_cell,
                                site_symmetry_table,
                                scatterers.ref(),
                                crystallographic_parameter_map,
                                crystallographic_shifts,
                                reparametrization_shifts);
      }
    }

    void place_constrained_scatterers()
    {
      for(int i=0; i < size(); ++i) {
        (*this)[i].place_constrained_scatterers(unit_cell,
                                                site_symmetry_table,
                                                scatterers.ref());
      }
    }

    std::map<int, xray::scatterer_flags>
    already_constrained() { return already_constrained_; }

  private:
    af::shared<constraint_type> elements;
    uctbx::unit_cell const &unit_cell;
    sgtbx::site_symmetry_table const &site_symmetry_table;
    af::shared<xray_scatterer_type> scatterers;
    parameter_map_type const &crystallographic_parameter_map;
    af::shared<xray::scatterer_flags> constraint_flags;
    std::map<int, xray::scatterer_flags> already_constrained_;
};


}}}

#endif // GUARD
