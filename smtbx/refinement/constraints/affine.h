#ifndef SMTBX_REFINEMENT_CONSTRAINTS_AFFINE_H
#define SMTBX_REFINEMENT_CONSTRAINTS_AFFINE_H

#include <smtbx/refinement/constraints/reparametrisation.h>

namespace smtbx { namespace refinement { namespace constraints {


/** Model a relation v = b + a_0 u_0 + a_1 u_1 + ... for scalar parameters
 *  v, u_0, u_1, ...
 */
class affine_scalar_parameter : public virtual scalar_parameter
{
public:
  double b;
  double *a;

  affine_scalar_parameter(scalar_parameter *u_0, double a_0,
                          double b)
  : parameter(1), a(new double[1]), b(b)
  {
    this->set_arguments(u_0);
    this->a[0] = a_0;
  }

  affine_scalar_parameter(scalar_parameter *u_0, double a_0,
                          scalar_parameter *u_1, double a_1,
                          double b)
  : parameter(2), a(new double[2]), b(b)
  {
    this->set_arguments(u_0, u_1);
    this->a[0] = a_0;
    this->a[1] = a_1;
  }

  affine_scalar_parameter(af::shared<scalar_parameter *> const &u,
                          af::shared<double> const &a, double b)
  : parameter(u.size()), a(new double[a.size()]), b(b)
  {
    SMTBX_ASSERT(u.size() == a.size())(u.size())(a.size());
    for (std::size_t i=0; i<n_arguments(); i++) {
      this->set_argument(i, u[i]);
      this->a[i] = a[i];
    }
  }

  ~affine_scalar_parameter() {
    delete[] a;
  }

  scalar_parameter *u(std::size_t i) const {
    return dynamic_cast<scalar_parameter *>(this->argument(i));
  }

  af::shared<double> affine_form() const {
    af::shared<double> result(a, a + n_arguments());
    result.push_back(b);
    return result;
  }

  virtual void linearise(uctbx::unit_cell const &unit_cell,
                         sparse_matrix_type *jacobian_transpose);
};

}}}

#endif
