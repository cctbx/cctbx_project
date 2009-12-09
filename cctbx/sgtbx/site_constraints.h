#ifndef CCTBX_SGTBX_SITE_CONSTRAINTS_H
#define CCTBX_SGTBX_SITE_CONSTRAINTS_H

#include <cctbx/sgtbx/space_group.h>
#include <scitbx/matrix/row_echelon.h>
#include <boost/scoped_array.hpp>

namespace cctbx { namespace sgtbx {

  //! Handling of site symmetry coordinate constraints.
  /*! Example:
        http://cci.lbl.gov/cctbx_sources/cctbx/examples/site_symmetry_constraints.py
   */
  template <typename FloatType=double>
  class site_constraints
  {
    protected:
      af::tiny<int, 9> row_echelon_form_memory;
      mutable af::tiny<FloatType, 9> gradient_sum_matrix_memory;
      mutable bool have_gradient_sum_matrix;

    public:
      int row_echelon_lcm;
      af::small<FloatType, 3> row_echelon_constants;
      af::small<unsigned, 3> independent_indices;

      //! Default constructor. Some data members are not initialized!
      site_constraints() {}

      site_constraints(
        af::const_ref<rt_mx> const& site_symmetry_matrices)
      :
        have_gradient_sum_matrix(false)
      {
        const rt_mx* const matrices = site_symmetry_matrices.begin();
        unsigned n_matrices = site_symmetry_matrices.size();
        CCTBX_ASSERT(n_matrices > 0);
        int lcm = 1;
        for(unsigned i=1;i<n_matrices;i++) {
          lcm = boost::lcm(lcm, matrices[i].r().den());
          lcm = boost::lcm(lcm, matrices[i].t().den());
        }
        row_echelon_lcm = lcm;
        unsigned n_rows = (n_matrices-1)*3;
        if (n_rows > 0) {
          boost::scoped_array<int> row_echelon_m(new int [n_rows*3]);
          boost::scoped_array<int> row_echelon_t(new int [n_rows]);
          int* rem = row_echelon_m.get();
          int* ret = row_echelon_t.get();
          for(unsigned i=1;i<n_matrices;i++) {
            rot_mx const& r = matrices[i].r();
            tr_vec const& t = matrices[i].t();
            // rotation part minus identity matrix
            int f = lcm / r.den();
            const int* n = r.num().begin();
            *rem++=(*n++)*f-lcm; *rem++=(*n++)*f;     *rem++=(*n++)*f;
            *rem++=(*n++)*f;     *rem++=(*n++)*f-lcm; *rem++=(*n++)*f;
            *rem++=(*n++)*f;     *rem++=(*n++)*f;     *rem++=(*n++)*f-lcm;
            // minus translation part
            f = -lcm / t.den();
            n = t.num().begin();
            *ret++ = (*n++)*f;
            *ret++ = (*n++)*f;
            *ret++ = (*n++)*f;
          }
          af::ref<int, af::mat_grid> rem_ref(row_echelon_m.get(), n_rows, 3);
          af::ref<int, af::mat_grid> ret_ref(row_echelon_t.get(), n_rows, 1);
          n_rows = scitbx::matrix::row_echelon::form_t(rem_ref, ret_ref);
          CCTBX_ASSERT(n_rows <= 3);
          std::copy(
            rem_ref.begin(), rem_ref.end(), row_echelon_form_memory.begin());
          for(unsigned i=0;i<n_rows;i++) {
            row_echelon_constants.push_back(
              static_cast<FloatType>(ret_ref[i]));
          }
        }
        af::tiny<bool, 3> independent_flags;
        CCTBX_ASSERT(scitbx::matrix::row_echelon::back_substitution_int(
          row_echelon_form(),
          (const int*) 0,
          (int*) 0,
          independent_flags.begin()));
        for(unsigned i=0;i<3;i++) {
          if (independent_flags[i]) {
            independent_indices.push_back(i);
          }
        }
      }

      af::const_ref<int, af::mat_grid>
      row_echelon_form() const
      {
        return af::const_ref<int, af::mat_grid>(
          row_echelon_form_memory.begin(),
          row_echelon_constants.size(),
          3);
      }

      unsigned
      n_independent_params() const { return independent_indices.size(); }

      unsigned
      n_dependent_params() const { return 3-independent_indices.size(); }

      af::small<FloatType, 3>
      independent_params(fractional<FloatType> const& all_params) const
      {
        af::small<FloatType, 3> result;
        for(std::size_t i=0;i<independent_indices.size();i++) {
          result.push_back(all_params[independent_indices[i]]);
        }
        return result;
      }

      fractional<FloatType>
      all_params(af::small<FloatType, 3> const& independent_params) const
      {
        fractional<FloatType> result(0,0,0);
        for(std::size_t i=0;i<independent_params.size();i++) {
          result[independent_indices[i]] = independent_params[i];
        }
        scitbx::matrix::row_echelon::back_substitution_float(
          row_echelon_form(),
          row_echelon_constants.begin(),
          result.begin());
        return result;
      }

      fractional<FloatType>
      all_shifts(af::small<FloatType, 3> const& independent_shifts) const
      {
        fractional<FloatType> result(0,0,0);
        for(std::size_t i=0;i<independent_shifts.size();i++) {
          result[independent_indices[i]] = independent_shifts[i];
        }
        scitbx::matrix::row_echelon::back_substitution_float(
          row_echelon_form(),
          static_cast<FloatType *>(0),
          result.begin());
        return result;
      }

    protected:
      void
      initialize_gradient_sum_matrix() const
      {
        FloatType* row = gradient_sum_matrix_memory.begin();
        std::fill_n(
          row, 3*independent_indices.size(), static_cast<FloatType>(0));
        af::const_ref<int, af::mat_grid> rem = row_echelon_form();
        for(std::size_t i=0;i<independent_indices.size();i++,row+=3) {
          row[independent_indices[i]] = 1;
          scitbx::matrix::row_echelon::back_substitution_float(
            rem, static_cast<const FloatType*>(0), row);
        }
        have_gradient_sum_matrix = true;
      }

    public:
      af::const_ref<FloatType, af::mat_grid>
      gradient_sum_matrix() const
      {
        if (!have_gradient_sum_matrix) initialize_gradient_sum_matrix();
        return af::const_ref<FloatType, af::mat_grid>(
          gradient_sum_matrix_memory.begin(),
          independent_indices.size(),
          3);
      }

      af::small<FloatType, 3>
      independent_gradients(
        af::const_ref<FloatType> const& all_gradients) const
      {
        CCTBX_ASSERT(all_gradients.size() == 3);
        if (!have_gradient_sum_matrix) initialize_gradient_sum_matrix();
        af::small<FloatType, 3> result;
        const FloatType *row = gradient_sum_matrix_memory.begin();
        for(std::size_t i=0;i<independent_indices.size();i++,row+=3) {
          FloatType indep_grad = 0;
          for(std::size_t j=0;j<3;j++) {
            indep_grad += row[j] * all_gradients[j];
          }
          result.push_back(indep_grad);
        }
        return result;
      }

      af::small<FloatType, 3>
      independent_gradients(scitbx::vec3<FloatType> const& all_gradients) const
      {
        return independent_gradients(all_gradients.const_ref());
      }


      /*! \brief Matrix product:
              gradient_sum_matrix
            * all_curvatures
            * gradient_sum_matrix.transpose()
       */
      /*! all_curvatures is the upper diagonal of
          the (3 x 3) matrix of curvatures assuming
          all tensor elements are independent.
          I.e. the number of elements of all_curvatures
          must be 3*(3+1)/2.

          The result is the upper triangle of the
          (n_independent_params() x n_independent_params())
          matrix of curvatures.
          I.e. the number of elements of the result array is
          n_independent_params()*(n_independent_params()+1)/2.
       */
      af::small<FloatType, 6>
      independent_curvatures(
        af::const_ref<FloatType> const& all_curvatures) const
      {
        CCTBX_ASSERT(all_curvatures.size() == 3*(3+1)/2);
        if (!have_gradient_sum_matrix) initialize_gradient_sum_matrix();
        unsigned n_indep = n_independent_params();
        af::tiny<FloatType, 9> ca;
        af::small<FloatType, 6> result(n_indep*(n_indep+1)/2, 0);
        scitbx::matrix::multiply_packed_u_multiply_lhs_transpose(
          gradient_sum_matrix_memory.begin(),
          all_curvatures.begin(),
          n_indep,
          3,
          ca.begin(),
          result.begin());
        return result;
      }
  };

}} // namespace cctbx::sgtbx

#endif // CCTBX_SGTBX_SITE_CONSTRAINTS_H
