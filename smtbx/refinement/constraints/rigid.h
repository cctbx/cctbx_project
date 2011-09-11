#ifndef SMTBX_REFINEMENT_CONSTRAINTS_RIGID_H
#define SMTBX_REFINEMENT_CONSTRAINTS_RIGID_H

#include <smtbx/refinement/constraints/reparametrisation.h>

namespace smtbx { namespace refinement { namespace constraints {

/** rigid group base - abstract, defines common components needed
for a rigid groups in this file
*/
class rigid_group_base : public asu_parameter {
public:
  rigid_group_base(af::shared<scatterer_type *> const &scatterers)
  : scatterers_(scatterers),
    co_s(scatterers.size()), fx_s(scatterers.size()),
    crd_initialised(false)
  {
    for (int i=0; i<scatterers.size(); i++)
      fx_s[i] = scatterers_[i]->site;
  }

  virtual ~rigid_group_base() {}

  virtual af::ref<double> components() {
    return af::ref<double>(fx_s[0].begin(), 3*fx_s.size());
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
      scatterers_[i]->site = fx_s[i];
  }

  virtual fractional<double> const& site(int i) const {  return fx_s[i];  }

protected:
  af::shared<scatterer_type *> scatterers_;
  af::shared<cart_t>
    co_s;  // original {Cartesian coordinates - geometrical center}
  af::shared<fractional<double> > fx_s;  // new fractional (for proxies)
  bool crd_initialised;
};

class rigid_site_proxy : public site_parameter {
  int index_in_parent;
public:
  rigid_site_proxy(rigid_group_base* parent, int index_in_parent)
  : parameter(1),
    index_in_parent(index_in_parent)
  {
    set_arguments(parent);
    value = parent->site(index_in_parent);
  }

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose)
  {
    rigid_group_base* parent = dynamic_cast<rigid_group_base *> (argument(0));
    value = parent->site(index_in_parent);
    if (!jacobian_transpose) return;
    sparse_matrix_type &jt = *jacobian_transpose;
    for (int j=0; j<3; ++j) {
      jt.col(index() + j) = jt.col(parent->index() + 3*index_in_parent + j);
    }
  }
};

/** a set of atoms rides on pivot and rotates around the direction given
  by pivot and pivot_neighbour (AFIX n=7/8 in shelxl)
 */

class pivoted_rotatable_group : public rigid_group_base {
public:
  pivoted_rotatable_group(site_parameter *pivot,
         site_parameter *pivot_neighbour,
         independent_scalar_parameter *azimuth,
         independent_scalar_parameter *size,
         af::shared<scatterer_type *> const &scatterers)
  : parameter(4),
    rigid_group_base(scatterers)
  {
    set_arguments(pivot, pivot_neighbour, azimuth, size);
  }

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);
};

/** a set of atoms rides on the pivot one; the group may shrink or
expand uniformly and rotate in 3D, an example is a spherecal counteranion or
'circluar' groups like Cp or Ph (AFIX n=9 in shelxl)
 */

class rotatable_expandable_group : public rigid_group_base {
public:
  rotatable_expandable_group(
    site_parameter *pivot,
    independent_scalar_parameter *size,
    independent_scalar_parameter *alpha,
    independent_scalar_parameter *beta,
    independent_scalar_parameter *gamma,
    af::shared<scatterer_type *> const &scatterers)
  : parameter(5),
    rigid_group_base(scatterers),
    shift_to_pivot(0,0,0)
  {
    set_arguments(pivot, size, alpha, beta, gamma);
  }

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);
protected:
  cart_t shift_to_pivot;
};


/** a set of atoms rides on the pivot one; bonds may be allowed stretching
    shelxl code 3/4
 */

class riding_expandable_group : public rigid_group_base {
public:
  riding_expandable_group(
    site_parameter *pivot,
    independent_scalar_parameter *size,
    af::shared<scatterer_type *> const &scatterers)
  : parameter(2),
    rigid_group_base(scatterers)
  {
    set_arguments(pivot, size);
  }

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);
};

}}}

#endif // GUARD
