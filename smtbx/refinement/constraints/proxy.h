#ifndef SMTBX_REFINEMENT_CONSTRAINTS_PROXY_H
#define SMTBX_REFINEMENT_CONSTRAINTS_PROXY_H

#include <smtbx/refinement/constraints/reparametrisation.h>

namespace smtbx { namespace refinement { namespace constraints {

template <class parent_t>
class site_proxy : public site_parameter {
  int index_in_parent;
public:
  site_proxy(parent_t* parent, int index_in_parent)
  : parameter(1),
    index_in_parent(index_in_parent)
  {
    set_arguments(parent);
    value = parent->site(index_in_parent);
  }

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose)
  {
    parent_t* parent = dynamic_cast<parent_t *> (argument(0));
    value = parent->site(index_in_parent);
    if (!jacobian_transpose) return;
    sparse_matrix_type &jt = *jacobian_transpose;
    for (int j=0; j<3; ++j) {
      jt.col(index() + j) = jt.col(parent->index() + 3*index_in_parent + j);
    }
  }
};

template <class parent_t>
class u_iso_proxy : public scalar_parameter {
  int index_in_parent;
public:
  u_iso_proxy(parent_t* parent, int index_in_parent)
  : parameter(1),
    index_in_parent(index_in_parent)
  {
    set_arguments(parent);
    value = parent->u_iso(index_in_parent);
  }

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose)
  {
    parent_t* parent = dynamic_cast<parent_t *> (argument(0));
    value = parent->u_iso(index_in_parent);
    if (!jacobian_transpose) return;
    sparse_matrix_type &jt = *jacobian_transpose;
    jt.col(index()) = jt.col(parent->index() + index_in_parent);
  }
};

template <class parent_t>
class u_star_proxy : public u_star_parameter {
  int index_in_parent;
public:
  u_star_proxy(parent_t* parent, int index_in_parent)
  : parameter(1),
    index_in_parent(index_in_parent)
  {
    set_arguments(parent);
    value = parent->u_star(index_in_parent);
  }

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose)
  {
    parent_t* parent = dynamic_cast<parent_t *> (argument(0));
    value = parent->u_star(index_in_parent);
    if (!jacobian_transpose) return;
    sparse_matrix_type &jt = *jacobian_transpose;
    for (int j=0; j<6; ++j) {
      jt.col(index() + j) = jt.col(parent->index() + 6*index_in_parent + j);
    }
  }
};

}}}

#endif // GUARD
