#ifndef CCTBX_SGTBX_TENSOR_RANK_2
#define CCTBX_SGTBX_TENSOR_RANK_2

#include <cctbx/sgtbx/space_group.h>
#include <scitbx/matrix/row_echelon.h>
#include <scitbx/matrix/tensor_rank_2.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>

namespace cctbx { namespace sgtbx { namespace tensor_rank_2 {

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
      af::small<af::tiny<FloatType, 6>, 6> gradient_sum_coeffs;

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
          gradient_sum_coeffs.resize(independent_indices.size());
          for(std::size_t i=0;i<independent_indices.size();i++) {
            gradient_sum_coeffs[i].fill(0);
            gradient_sum_coeffs[i][independent_indices[i]] = 1;
            scitbx::matrix::row_echelon::back_substitution_float(
              row_echelon_form_ref,
              static_cast<const FloatType*>(0),
              gradient_sum_coeffs[i].begin());
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
      independent_gradients(scitbx::sym_mat3<FloatType> const& all_gradients)
      {
        if (gradient_average_denominator == 0) {
          throw error("gradient handling was not initialized.");
        }
        af::small<FloatType, 6> result;
        for(std::size_t i=0;i<gradient_sum_coeffs.size();i++) {
          FloatType indep_grad = 0;
          for(std::size_t j=0;j<6;j++) {
            indep_grad += gradient_sum_coeffs[i][j] * all_gradients[j];
          }
          result.push_back(indep_grad);
        }
        return result;
      }
  };

}}} // namespace cctbx::sgtbx::tensor_rank_2

#endif // CCTBX_SGTBX_TENSOR_RANK_2
