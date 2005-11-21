#ifndef CCTBX_SGTBX_TENSOR_RANK_2
#define CCTBX_SGTBX_TENSOR_RANK_2

#include <cctbx/sgtbx/site_symmetry.h>
#include <scitbx/matrix/row_echelon.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>

namespace cctbx { namespace sgtbx {

//! Symmetry constraints for tensors of rank 2
namespace tensor_rank_2 {

  //! Coefficients r.transpose() * c * r - c = 0
  /*! Mathematica code:
        r={{r0,r1,r2},{r3,r4,r5},{r6,r7,r8}}
        c={{c0,c3,c4},{c3,c1,c5},{c4,c5,c2}}
        FortranForm[Expand[Transpose[r].c.r - c]]
   */
  int*
  constraints_raw(const int* r, int* c);

  //! Coefficients r.transpose() * c * r - c = 0 for all symmetry_matrices.
  int*
  constraints_raw(
    af::const_ref<rt_mx> const& symmetry_matrices,
    std::size_t i_first_matrix_to_use,
    bool reciprocal_space,
    int* c);

  //! Handling of symmetry constraints for tensors of rank 2.
  template <typename FloatType=double>
  class constraints
  {
    protected:
      af::versa<int, af::c_grid<2> > row_echelon_form_memory;
      mutable af::versa<FloatType, af::c_grid<2> > gradient_sum_matrix_memory;
      mutable bool have_gradient_sum_matrix;

    public:
      af::small<unsigned, 6> independent_indices;

      constraints() {}

      constraints(
        sgtbx::space_group const& space_group,
        bool reciprocal_space)
      {
        initialize(space_group.smx().const_ref(), 1, reciprocal_space);
      }

      constraints(
        sgtbx::site_symmetry_ops const& site_symmetry_ops,
        bool reciprocal_space)
      {
        initialize(
          site_symmetry_ops.matrices().const_ref(), 1, reciprocal_space);
      }

      constraints(
        af::shared<rt_mx> const& symmetry_matrices,
        std::size_t i_first_matrix_to_use,
        bool reciprocal_space)
      {
        initialize(
          symmetry_matrices.const_ref(),
          i_first_matrix_to_use,
          reciprocal_space);
      }

      constraints(
        af::const_ref<rt_mx> const& symmetry_matrices,
        std::size_t i_first_matrix_to_use,
        bool reciprocal_space)
      {
        initialize(symmetry_matrices, i_first_matrix_to_use, reciprocal_space);
      }

      void
      initialize(
        af::const_ref<rt_mx> const& symmetry_matrices,
        std::size_t i_first_matrix_to_use,
        bool reciprocal_space);

      scitbx::mat_const_ref<int>
      row_echelon_form() const
      {
        return scitbx::mat_const_ref<int>(
          row_echelon_form_memory.begin(),
          6-independent_indices.size(),
          6);
      }

      std::size_t
      n_independent_params() const { return independent_indices.size(); }

      std::size_t
      n_dependent_params() const { return 6-independent_indices.size(); }

      af::small<FloatType, 6>
      independent_params(scitbx::sym_mat3<FloatType> const& all_params) const
      {
        af::small<FloatType, 6> result;
        for(std::size_t i=0;i<independent_indices.size();i++) {
          result.push_back(all_params[independent_indices[i]]);
        }
        return result;
      }

      scitbx::sym_mat3<FloatType>
      all_params(af::small<FloatType, 6> const& independent_params) const
      {
        scitbx::sym_mat3<FloatType> result(0,0,0,0,0,0);
        for(std::size_t i=0;i<independent_params.size();i++) {
          result[independent_indices[i]] = independent_params[i];
        }
        scitbx::matrix::row_echelon::back_substitution_float(
          row_echelon_form(),
          static_cast<const FloatType*>(0),
          result.begin());
        return result;
      }

    protected:
      void
      initialize_gradient_sum_matrix() const
      {
        gradient_sum_matrix_memory.resize(
          af::c_grid<2>(independent_indices.size(), 6), 0);
        FloatType* row = gradient_sum_matrix_memory.begin();
        scitbx::mat_const_ref<int> rem = row_echelon_form();
        for(std::size_t i=0;i<independent_indices.size();i++,row+=6) {
          row[independent_indices[i]] = 1;
          scitbx::matrix::row_echelon::back_substitution_float(
            rem, static_cast<const FloatType*>(0), row);
        }
        have_gradient_sum_matrix = true;
      }

    public:
      scitbx::mat_const_ref<FloatType>
      gradient_sum_matrix() const
      {
        if (!have_gradient_sum_matrix) initialize_gradient_sum_matrix();
        return scitbx::mat_const_ref<FloatType>(
          gradient_sum_matrix_memory.begin(),
          independent_indices.size(),
          6);
      }

      af::small<FloatType, 6>
      independent_gradients(
        scitbx::sym_mat3<FloatType> const& all_gradients) const;

      /*! \brief Efficient implementation of the matrix product:
              gradient_sum_matrix
            * all_curvatures
            * gradient_sum_matrix.transpose()
       */
      /*! all_curvatures is the upper diagonal of
          the (6 x 6) matrix of curvatures assuming
          all tensor elements are independent.
          I.e. the number of elements of all_curvatures
          must be 6*(6+1)/2.

          The result is the upper diagonal of the
          (n_independent_params() x n_independent_params())
          matrix of curvatures.
          I.e. the number of elements of the result array is
          n_independent_params()*(n_independent_params()+1)/2.
       */
      af::shared<FloatType>
      independent_curvatures(
        af::const_ref<FloatType> const& all_curvatures) const;
  };

  template <typename FloatType>
  void
  constraints<FloatType>::initialize(
    af::const_ref<rt_mx> const& symmetry_matrices,
    std::size_t i_first_matrix_to_use,
    bool reciprocal_space)
  {
    CCTBX_ASSERT(i_first_matrix_to_use <= symmetry_matrices.size());
    row_echelon_form_memory = af::versa<int, af::c_grid<2> >(
      af::c_grid<2>((symmetry_matrices.size()-i_first_matrix_to_use)*6, 6),
      af::init_functor_null<int>());
    scitbx::mat_ref<int> row_echelon_form_ref(
      row_echelon_form_memory.begin(),
      row_echelon_form_memory.accessor()[0],
      row_echelon_form_memory.accessor()[1]);
    CCTBX_ASSERT(constraints_raw(
      symmetry_matrices,
      i_first_matrix_to_use,
      reciprocal_space,
      row_echelon_form_ref.begin()) == row_echelon_form_ref.end());
    CCTBX_ASSERT(
      scitbx::matrix::row_echelon::form(row_echelon_form_ref) <= 6);
    row_echelon_form_memory.resize(af::c_grid<2>(
      row_echelon_form_ref.n_rows(),
      row_echelon_form_ref.n_columns()));
    af::tiny<bool, 6> independent_flags;
    scitbx::matrix::row_echelon::back_substitution_int(
      row_echelon_form_ref,
      static_cast<const int*>(0),
      static_cast<int*>(0),
      independent_flags.begin());
    for(unsigned i=0;i<6;i++) {
      if (independent_flags[i]) {
        independent_indices.push_back(i);
      }
    }
    have_gradient_sum_matrix = false;
  }

  template <typename FloatType>
  af::small<FloatType, 6>
  constraints<FloatType>::independent_gradients(
    scitbx::sym_mat3<FloatType> const& all_gradients) const
  {
    af::small<FloatType, 6> result;
    if (!have_gradient_sum_matrix) initialize_gradient_sum_matrix();
    const FloatType *row = gradient_sum_matrix_memory.begin();
    for(std::size_t i=0;i<independent_indices.size();i++,row+=6) {
      FloatType indep_grad = 0;
      for(std::size_t j=0;j<6;j++) {
        indep_grad += row[j] * all_gradients[j];
      }
      result.push_back(indep_grad);
    }
    return result;
  }

  template <typename FloatType>
  af::shared<FloatType>
  constraints<FloatType>::independent_curvatures(
    af::const_ref<FloatType> const& all_curvatures) const
  {
    CCTBX_ASSERT(all_curvatures.size() == 6*(6+1)/2);
    static const unsigned sym_mat6_indices[] = {
      0, 1, 2, 3, 4, 5,
      1, 6, 7, 8, 9,10,
      2, 7,11,12,13,14,
      3, 8,12,15,16,17,
      4, 9,13,16,18,19,
      5,10,14,17,19,20
    };
    if (!have_gradient_sum_matrix) initialize_gradient_sum_matrix();
    const FloatType* gsm = gradient_sum_matrix_memory.begin();
    unsigned n_indep = n_independent_params();
    af::shared<FloatType> result(n_indep*(n_indep+1)/2, 0);
    std::vector<FloatType> c_times_a(n_indep*6, 0);
    typename std::vector<FloatType>::iterator ca = c_times_a.begin();
    unsigned i6 = 0;
    for (unsigned i=0;i<n_indep;i++,i6+=6) {
      for (unsigned k=0;k<6;k++) {
        unsigned i6j = i6;
        unsigned j6k = k;
        for (unsigned j=0;j<6;j++,i6j++,j6k+=6) {
          *ca += gsm[i6j] * all_curvatures[sym_mat6_indices[j6k]];
        }
        ca++;
      }
    }
    FloatType* cact = result.begin();
    i6 = 0;
    for (unsigned i=0;i<n_indep;i++,i6+=6) {
      unsigned k6 = i6;
      for (unsigned k=i;k<n_indep;k++,k6+=6) {
        unsigned i6j = i6;
        unsigned k6j = k6;
        for (unsigned j=0;j<6;j++,i6j++,k6j++) {
          *cact += c_times_a[i6j] * gsm[k6j];
        }
        cact++;
      }
    }
    return result;
  }

}}} // namespace cctbx::sgtbx::tensor_rank_2

#endif // CCTBX_SGTBX_TENSOR_RANK_2
