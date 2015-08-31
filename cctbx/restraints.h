#ifndef CCTBX_RESTRAINTS_H
#define CCTBX_RESTRAINTS_H

#include <cctbx/import_scitbx_af.h>
#include <cctbx/xray/parameter_map.h>
#include <cctbx/xray/scatterer.h>

#include <scitbx/sparse/matrix.h>


namespace cctbx { namespace restraints {

  /// The linearised equations of restraints.
  /*
      Take advantage of the fact that the restraints part of the design matrix
      is sparse by constructing the design matrix for restraints only, along
      with the associated vectors of weights and deltas.

      The normal equations can then be obtained separately, and the normal
      equations derived from the observations can updated with the contribution
      from the restraints.
   */

  // Used by smtbx/refinement/restraints

  template <typename FloatType>
  class linearised_eqns_of_restraint
  {
  public:
    typedef FloatType scalar_t;

    std::size_t n_columns, n_rows;
    scitbx::sparse::matrix<scalar_t> design_matrix;
    scitbx::af::shared<scalar_t> weights;
    scitbx::af::shared<scalar_t> deltas;

  private:
    std::size_t row_i;

  public:
    linearised_eqns_of_restraint(std::size_t n_rows_, std::size_t n_columns_)
      : n_rows(n_rows_),
        n_columns(n_columns_),
        design_matrix(n_rows_, n_columns_),
        weights(n_rows_), deltas(n_rows_),
        row_i(0)
    {}

    std::size_t next_row() {
      CCTBX_ASSERT(!finalised())(row_i)(n_rows);
      return row_i++;
    }

    void add_equation(
      FloatType r, af::const_ref<FloatType> const &gradient, FloatType w) {
      CCTBX_ASSERT(gradient.size() == n_crystallographic_params());
      std::size_t row_i = next_row();
      deltas[row_i] = r;
      weights[row_i] = w;
      for (std::size_t col_i=0; col_i<gradient.size(); col_i++) {
        design_matrix(row_i, col_i) = gradient[col_i];
      }
    }

    bool finalised() { return row_i >= n_rows; }

    std::size_t n_restraints() { return row_i; }

    std::size_t n_crystallographic_params() { return n_columns; }

  };

}} // cctbx::restraints

#endif // GUARD
