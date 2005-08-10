#ifndef CCTBX_SGTBX_TENSOR_RANK_2
#define CCTBX_SGTBX_TENSOR_RANK_2

#include <cctbx/sgtbx/space_group.h>
#include <scitbx/matrix/row_echelon.h>
#include <scitbx/matrix/tensor_rank_2.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>

namespace cctbx { namespace sgtbx {

//! Symmetry constraints for tensors of rank 2
namespace tensor_rank_2 {

  //! Coefficients r.transpose() * t * r - t = 0
  /*! Mathematica code:
        r={{r0,r1,r2},{r3,r4,r5},{r6,r7,r8}}
        t={{t0,t3,t4},{t3,t1,t5},{t4,t5,t2}}
        FortranForm[Expand[Transpose[r].t.r - t]]
   */
  int*
  constraints_raw(
    sgtbx::space_group const& space_group,
    bool reciprocal_space,
    int* c);

  //! Handling of symmetry constraints for tensors of rank 2.
  template <typename FloatType=double>
  class constraints
  {
    public:
      af::versa<int, af::c_grid<2> > row_echelon_form_memory;
      scitbx::mat_ref<int> row_echelon_form_ref;
      af::small<unsigned, 6> independent_indices;
      std::size_t gradient_average_denominator;
      scitbx::matrix::tensor_rank_2::gradient_average_cache<int>
        gradient_average_cache;
      af::versa<FloatType, af::c_grid<2> > gradient_sum_coeffs;

      constraints() {}

      constraints(
        sgtbx::space_group const& space_group,
        bool reciprocal_space,
        bool initialize_gradient_handling=false)
      :
        row_echelon_form_memory(
          af::c_grid<2>((space_group.n_smx()-1)*6, 6),
          af::init_functor_null<int>())
      {
        row_echelon_form_ref = scitbx::mat_ref<int>(
          row_echelon_form_memory.begin(),
          row_echelon_form_memory.accessor()[0],
          row_echelon_form_memory.accessor()[1]);
        CCTBX_ASSERT(constraints_raw(space_group, reciprocal_space,
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
        if (!initialize_gradient_handling) {
          gradient_average_denominator = 0;
        }
        else {
          gradient_average_denominator = space_group.n_smx();
          for(std::size_t i_smx=0;i_smx<gradient_average_denominator;i_smx++) {
            rot_mx const& r = space_group.smx()[i_smx].r();
            CCTBX_ASSERT(r.den() == 1);
            gradient_average_cache.accumulate(r.num());
          }
          gradient_sum_coeffs.resize(
            af::c_grid<2>(independent_indices.size(), 6), 0);
          FloatType* row = gradient_sum_coeffs.begin();
          for(std::size_t i=0;i<independent_indices.size();i++,row+=6) {
            row[independent_indices[i]] = 1;
            scitbx::matrix::row_echelon::back_substitution_float(
              row_echelon_form_ref,
              static_cast<const FloatType*>(0),
              row);
          }
        }
      }

      std::size_t
      n_independent_params() const { return independent_indices.size(); }

      std::size_t
      n_dependent_params() const { return 6-independent_indices.size(); }

      af::small<FloatType, 6>
      independent_params(scitbx::sym_mat3<FloatType> const& params) const
      {
        af::small<FloatType, 6> result;
        for(std::size_t i=0;i<independent_indices.size();i++) {
          result.push_back(params[independent_indices[i]]);
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
          row_echelon_form_ref,
          static_cast<const FloatType*>(0),
          result.begin());
        return result;
      }

      scitbx::sym_mat3<FloatType>
      sym_gradients(scitbx::sym_mat3<FloatType> const& asu_gradients) const
      {
        if (gradient_average_denominator == 0) {
          throw error("gradient handling was not initialized.");
        }
        return gradient_average_cache.average(
          asu_gradients,
          static_cast<FloatType>(gradient_average_denominator));
      }

      af::small<FloatType, 6>
      independent_gradients(
        scitbx::sym_mat3<FloatType> const& all_gradients) const
      {
        if (gradient_average_denominator == 0) {
          throw error("gradient handling was not initialized."); }
        af::small<FloatType, 6> result;
        const FloatType *row = gradient_sum_coeffs.begin();
        for(std::size_t i=0;i<independent_indices.size();i++,row+=6) {
          FloatType indep_grad = 0;
          for(std::size_t j=0;j<6;j++) {
            indep_grad += row[j] * all_gradients[j];
          }
          result.push_back(indep_grad);
        }
        return result;
      }

      /*! \brief Efficient implementation of the matrix product:
              gradient_sum_coeffs
            * all_curvatures
            * gradient_sum_coeffs.transpose()
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
        af::const_ref<FloatType> const& all_curvatures) const
      {
        if (gradient_average_denominator == 0) {
          throw error("gradient handling was not initialized."); }
        CCTBX_ASSERT(all_curvatures.size() == 6*(6+1)/2);
        static const unsigned sym_mat6_indices[] = {
          0, 1, 2, 3, 4, 5,
          1, 6, 7, 8, 9,10,
          2, 7,11,12,13,14,
          3, 8,12,15,16,17,
          4, 9,13,16,18,19,
          5,10,14,17,19,20
        };
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
              *ca += gradient_sum_coeffs[i6j]
                   * all_curvatures[sym_mat6_indices[j6k]];
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
              *cact += c_times_a[i6j] * gradient_sum_coeffs[k6j];
            }
            cact++;
          }
        }
        return result;
      }
  };

}}} // namespace cctbx::sgtbx::tensor_rank_2

#endif // CCTBX_SGTBX_TENSOR_RANK_2
