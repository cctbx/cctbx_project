#ifndef SMTBX_REFINEMENT_CONSTRAINTS_RIGID_H
#define SMTBX_REFINEMENT_CONSTRAINTS_RIGID_H

#include <smtbx/refinement/constraints/reparametrisation.h>

namespace smtbx { namespace refinement { namespace constraints {

class rigid_group_base : public virtual parameter
{
public:
  virtual ~rigid_group_base() {}
  virtual fractional<double> const& site(int index) const = 0;
};


/** a set of scatterers rides on pivot and rotates around the direction given
  by pivot and pivot_neighbour
 */
class pivoted_rotable_group : public rigid_group_base, public asu_parameter {
public:
  pivoted_rotable_group(site_parameter *pivot,
         site_parameter *pivot_neighbour,
         independent_scalar_parameter *azimuth,
         const af::shared<scatterer_type *>& scatterers)
  : parameter(3),
    scatterers_(scatterers),
    cx_s(scatterers.size()), co_s(scatterers.size()), fx_s(scatterers.size()),
    original_crd_initialised(false)
  {
    this->set_arguments(pivot, pivot_neighbour, azimuth);
    for (int i=0; i<scatterers_.size(); i++)
      fx_s[i] = scatterers_[i]->site;
  }

  virtual af::ref<double> components() {
    return af::ref<double>(cx_s[0].begin(), 3*cx_s.size());
  }

  virtual scatterer_sequence_type scatterers() const {
    return scatterers_.const_ref();
  }

  virtual index_range
  component_indices_for(scatterer_type const *scatterer) const;

  virtual void
  write_component_annotations_for(scatterer_type const *scatterer,
                                  std::ostream &output) const;

  virtual void store(uctbx::unit_cell const &unit_cell) const {
    for (int i=0; i<scatterers_.size(); i++)
      scatterers_[i]->site = unit_cell.fractionalize(cx_s[i]);
  }

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);

  virtual fractional<double> const& site(int i) const {  return fx_s[i];  }

protected:
  af::shared<scatterer_type *> scatterers_;
  af::shared<cart_t>
    cx_s,  //new Cartesian coordinates
    co_s;  // original Cartesian coordinates
  af::shared<fractional<double> > fx_s;  // new fractional (for proxies)
  bool original_crd_initialised;
};


class rigid_site_proxy : public site_parameter {
  rigid_group_base* parent;
  int index;
public:
  rigid_site_proxy(rigid_group_base* parent, int index)
  : parameter(0),
    parent(parent),
    index(index)
  {
    value = parent->site(index);
    set_variable(false);
  }

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose)
  {
    value = parent->site(index);
  }
};

}}}

#endif // GUARD
