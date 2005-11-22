#ifndef CCTBX_SGTBX_TENSOR_RANK_2
#define CCTBX_SGTBX_TENSOR_RANK_2

#include <cctbx/sgtbx/space_group.h>
#include <scitbx/matrix/row_echelon.h>
#include <boost/shared_array.hpp>
#include <boost/scoped_array.hpp>

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
  /*! Example:
        http://phenix-online.org/cctbx_sources/cctbx/cctbx/examples/adp_symmetry_constraints.py
   */
  template <typename FloatType=double>
  class constraints
  {
    protected:
      boost::shared_array<int> row_echelon_form_memory;
      mutable boost::shared_array<FloatType> gradient_sum_matrix_memory;
    public:
      af::small<unsigned, 6> independent_indices;

      //! Default constructor. Some data members are not initialized!
      constraints() {}

      constraints(
        sgtbx::space_group const& space_group,
        bool reciprocal_space)
      {
        initialize(space_group.smx().const_ref(), 1, reciprocal_space);
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

      scitbx::mat_const_ref<int>
      row_echelon_form() const
      {
        return scitbx::mat_const_ref<int>(
          row_echelon_form_memory.get(), 6-independent_indices.size(), 6);
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

      scitbx::mat_const_ref<FloatType>
      gradient_sum_matrix() const
      {
        const FloatType* gsm = gradient_sum_matrix_memory.get();
        if (gsm == 0) gsm = initialize_gradient_sum_matrix();
        return scitbx::mat_const_ref<FloatType>(
          gsm, independent_indices.size(), 6);
      }

      af::small<FloatType, 6>
      independent_gradients(
        scitbx::sym_mat3<FloatType> const& all_gradients) const;

      /*! \brief Matrix product:
              gradient_sum_matrix
            * all_curvatures
            * gradient_sum_matrix.transpose()
       */
      /*! all_curvatures is the upper triangle of
          the (6 x 6) matrix of curvatures assuming
          all tensor elements are independent.
          I.e. the number of elements of all_curvatures
          must be 6*(6+1)/2.

          The result is the upper triangle of the
          (n_independent_params() x n_independent_params())
          matrix of curvatures.
          I.e. the number of elements of the result array is
          n_independent_params()*(n_independent_params()+1)/2.
       */
      af::shared<FloatType>
      independent_curvatures(
        af::const_ref<FloatType> const& all_curvatures) const;

    protected:
      void
      initialize(
        af::const_ref<rt_mx> const& symmetry_matrices,
        std::size_t i_first_matrix_to_use,
        bool reciprocal_space);

      const FloatType*
      initialize_gradient_sum_matrix() const;
  };

  template <typename FloatType>
  void
  constraints<FloatType>::
  initialize(
    af::const_ref<rt_mx> const& symmetry_matrices,
    std::size_t i_first_matrix_to_use,
    bool reciprocal_space)
  {
    CCTBX_ASSERT(i_first_matrix_to_use <= symmetry_matrices.size());
    unsigned n_rows = (symmetry_matrices.size() - i_first_matrix_to_use) * 6;
    boost::shared_array<int> row_echelon_setup_memory(new int[n_rows*6]);
    scitbx::mat_ref<int> row_echelon_setup(
      row_echelon_setup_memory.get(), n_rows, 6);
    CCTBX_ASSERT(constraints_raw(
      symmetry_matrices,
      i_first_matrix_to_use,
      reciprocal_space,
      row_echelon_setup.begin()) == row_echelon_setup.end());
    n_rows = scitbx::matrix::row_echelon::form(row_echelon_setup);
    CCTBX_ASSERT(n_rows <= 6);
    row_echelon_form_memory = boost::shared_array<int>(new int[n_rows*6]);
    std::copy(
      row_echelon_setup.begin(),
      row_echelon_setup.end(),
      row_echelon_form_memory.get());
    af::tiny<bool, 6> independent_flags;
    scitbx::matrix::row_echelon::back_substitution_int(
      row_echelon_setup,
      static_cast<const int*>(0),
      static_cast<int*>(0),
      independent_flags.begin());
    for(unsigned i=0;i<6;i++) {
      if (independent_flags[i]) {
        independent_indices.push_back(i);
      }
    }
  }

  template <typename FloatType>
  const FloatType*
  constraints<FloatType>::
  initialize_gradient_sum_matrix() const
  {
    gradient_sum_matrix_memory = boost::shared_array<FloatType>(
      new FloatType[independent_indices.size()*6]);
    FloatType* row = gradient_sum_matrix_memory.get();
    std::fill_n(
      row, independent_indices.size()*6, static_cast<FloatType>(0));
    scitbx::mat_const_ref<int> rem = row_echelon_form();
    for(std::size_t i=0;i<independent_indices.size();i++,row+=6) {
      row[independent_indices[i]] = 1;
      scitbx::matrix::row_echelon::back_substitution_float(
        rem, static_cast<const FloatType*>(0), row);
    }
    return gradient_sum_matrix_memory.get();
  }

  template <typename FloatType>
  af::small<FloatType, 6>
  constraints<FloatType>::
  independent_gradients(
    scitbx::sym_mat3<FloatType> const& all_gradients) const
  {
    af::small<FloatType, 6> result;
    const FloatType* row = gradient_sum_matrix_memory.get();
    if (row == 0) row = initialize_gradient_sum_matrix();
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
  constraints<FloatType>::
  independent_curvatures(
    af::const_ref<FloatType> const& all_curvatures) const
  {
    CCTBX_ASSERT(all_curvatures.size() == 6*(6+1)/2);
    const FloatType* gsm = gradient_sum_matrix_memory.get();
    if (gsm == 0) gsm = initialize_gradient_sum_matrix();
    unsigned n_indep = n_independent_params();
    af::shared<FloatType> result(
      n_indep*(n_indep+1)/2, af::init_functor_null<FloatType>());
    boost::scoped_array<FloatType> c_times_a(new FloatType[n_indep*6]);
    scitbx::matrix::multiply_packed_u_multiply_lhs_transpose(
      gsm, all_curvatures.begin(), n_indep, 6,
      &*c_times_a.get(), result.begin());
    return result;
  }

}}} // namespace cctbx::sgtbx::tensor_rank_2

#endif // CCTBX_SGTBX_TENSOR_RANK_2
