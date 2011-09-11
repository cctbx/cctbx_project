#ifndef SMTBX_REFINEMENT_CONSTRAINTS_SAME_H
#define SMTBX_REFINEMENT_CONSTRAINTS_SAME_H

#include <smtbx/refinement/constraints/rigid.h>

namespace smtbx { namespace refinement { namespace constraints {

/** Non-crystallographic symmetry constraint
 */
// manages sites
class same_group_xyz : public asu_parameter {
public:
  same_group_xyz(
    af::shared<scatterer_type *> const &scatterers,
    af::shared<site_parameter *> const &sites,
    scitbx::mat3<double> const &alignment_matrix_,
    independent_small_vector_parameter<6> *shifts_and_angles)
  : parameter(scatterers.size()+1),
    fx_s(scatterers.size()),
    scatterers_(scatterers),
    alignment_matrix(alignment_matrix_)
  {
    SMTBX_ASSERT(sites.size()==scatterers.size());
    set_argument(0, shifts_and_angles);
    for (int i=0; i<scatterers_.size(); i++) {
      set_argument(i+1, sites[i]);
      fx_s[i] = scatterers[i]->site;
    }
  }

  virtual scatterer_sequence_type scatterers() const {
    return scatterers_.const_ref();
  }

  virtual af::ref<double> components() {
    return af::ref<double>(fx_s[0].begin(), 3*fx_s.size());
  }

  index_range component_indices_for(scatterer_type const *scatterer) const;

  void write_component_annotations_for(scatterer_type const *scatterer,
                                    std::ostream &output) const;

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);

  virtual void store(uctbx::unit_cell const &unit_cell) const {
    for (int i=0; i<scatterers_.size(); i++)
      scatterers_[i]->site = fx_s[i];
  }
protected:
  af::shared<fractional<double> > fx_s;
  af::shared<scatterer_type *> scatterers_;
  scitbx::mat3<double> alignment_matrix;
};

// manages u_isos
class same_group_u_iso : public asu_parameter {
public:
  same_group_u_iso(
    af::shared<scatterer_type *> const &scatterers,
    af::shared<scalar_parameter *> const &u_isos_,
    independent_small_vector_parameter<6> *shifts_and_angles)
  : parameter(scatterers.size()+1),
    scatterers_(scatterers),
    u_isos(scatterers.size())
  {
    SMTBX_ASSERT(u_isos_.size()==scatterers.size());
    set_argument(0, shifts_and_angles);
    for (int i=0; i<scatterers.size(); i++) {
      set_argument(i+1, u_isos_[i]);
      u_isos[i] = scatterers[i]->u_iso;
    }
  }

  virtual scatterer_sequence_type scatterers() const {
    return scatterers_.const_ref();
  }

  virtual af::ref<double> components() {
    return u_isos.ref();
  }

  index_range component_indices_for(scatterer_type const *scatterer) const;

  void write_component_annotations_for(scatterer_type const *scatterer,
                                    std::ostream &output) const;

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);

  virtual void store(uctbx::unit_cell const &unit_cell) const {
    for (int i=0; i<scatterers_.size(); i++)
      scatterers_[i]->u_iso = u_isos[i];
  }
protected:
  af::shared<scatterer_type *> scatterers_;
  af::shared<double> u_isos;
};

//manages u_stars
class same_group_u_star : public asu_parameter {
public:
  same_group_u_star(
    af::shared<scatterer_type *> const &scatterers,
    af::shared<u_star_parameter *> const &u_stars_,
    scitbx::mat3<double> const &alignment_matrix,
    independent_small_vector_parameter<6> *shifts_and_angles)
  : parameter(scatterers.size()+1),
    u_stars(scatterers.size()),
    scatterers_(scatterers),
    alignment_matrix(alignment_matrix)
  {
    SMTBX_ASSERT(u_stars_.size()==scatterers_.size());
    set_argument(0, shifts_and_angles);
    for (int i=0; i<scatterers_.size(); i++) {
      set_argument(i+1, u_stars_[i]);
      u_stars[i] = scatterers_[i]->u_star;
    }
  }

  virtual scatterer_sequence_type scatterers() const {
    return scatterers_.const_ref();
  }

  virtual af::ref<double> components() {
    return af::ref<double>(u_stars[0].begin(), 6*u_stars.size());
  }

  index_range component_indices_for(scatterer_type const *scatterer) const;

  void write_component_annotations_for(scatterer_type const *scatterer,
                                    std::ostream &output) const;

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);

  virtual void store(uctbx::unit_cell const &unit_cell) const {
    for (int i=0; i<scatterers_.size(); i++)
      scatterers_[i]->u_star = u_stars[i];
  }
protected:
  af::shared<tensor_rank_2_t> u_stars;
  af::shared<scatterer_type *> scatterers_;
  scitbx::mat3<double> alignment_matrix;
};

}}}

#endif // GUARD
