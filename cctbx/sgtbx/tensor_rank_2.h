#ifndef CCTBX_SGTBX_TENSOR_RANK_2
#define CCTBX_SGTBX_TENSOR_RANK_2

#include <cctbx/sgtbx/space_group.h>
#include <scitbx/matrix/delta_tensors.h>
#include <scitbx/matrix/row_echelon.h>
#include <scitbx/matrix/row_echelon_full_pivoting_small.h>
#include <scitbx/array_family/owning_ref.h>
#include <boost/shared_array.hpp>

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
        http://cci.lbl.gov/cctbx_sources/cctbx/examples/adp_symmetry_constraints.py
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

      af::const_ref<int, af::mat_grid>
      row_echelon_form() const
      {
        return af::const_ref<int, af::mat_grid>(
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

      af::const_ref<FloatType, af::mat_grid>
      gradient_sum_matrix() const
      {
        const FloatType* gsm = gradient_sum_matrix_memory.get();
        if (gsm == 0) gsm = initialize_gradient_sum_matrix();
        return af::const_ref<FloatType, af::mat_grid>(
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
    af::ref<int, af::mat_grid> row_echelon_setup(
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
    af::const_ref<int, af::mat_grid> rem = row_echelon_form();
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
      c_times_a.get(), result.begin());
    return result;
  }

  //! Outside template class to avoid complications.
  /*! Anyway, even the standard does not allow the declaration of anything
      but an integral static member inside a class (sigh)
  */
  namespace cartesian_constraints_constants {
    static const unsigned n_all_params = 6;
    static const unsigned symmetries_max_nb = 24;
    static const double pivot_zero_attractor = 1e-9;
  }

  template <typename T=double>
  class cartesian_constraints
  {
      /* This is the n_all_params x n_independent_params
         matrix Z whose columns make a basis of the null space of A:
         A u = 0   =>   u = Z u'
         where u' is the vector of independent parameters and therefore
         grad_u' = Z^T grad_u
      */

      af::ref_owning_versa<T, af::mat_grid> z;
      unsigned n_independent_params;

      // the min absolute value for a pivot value to be considered null
      double pivot_zero_attractor;

    protected:
      void
      initialise(uctbx::unit_cell const & unit_cell,
                 af::const_ref<sgtbx::rt_mx> const& symmetries)
     {
        using namespace cartesian_constraints_constants;
        using scitbx::matrix::delta_x_delta;
        /* Construct the matrix A s.t. A u = 0
          where u is the vector of U_cart components built from
          R^T U R - U = 0
          and
          (R^T U R)_ij =  sum_k R_ki R_kj  U_kk
          + 2 sum_{k<l} (R_ki R_lj + R_kj R_li) U_kl
          where i,j,k,l are in {0,1,2} with i <= j
          */
        const unsigned n_a_rows = symmetries.size() * n_all_params;
        const unsigned max_n_a_rows = symmetries_max_nb * n_all_params;
        boost::shared_array<T> a(new T[n_a_rows*n_all_params] );
        T* pa = a.get();
        for (unsigned i_s=0; i_s < symmetries.size(); i_s++) {
          scitbx::mat3<T> r_f = symmetries[i_s]
                                  .r()
                                  .as_floating_point(scitbx::type_holder<T>());
          scitbx::mat3<T> r = unit_cell.orthogonalization_matrix()
                              * r_f
                              * unit_cell.fractionalization_matrix();
          for (int i=0; i < 3; i++) for (int j=i; j < 3; j++) {
            for (int k=0; k < 3; k++) {
              *pa++ = r(k,i)*r(k,j) - delta_x_delta<T>(k,i,k,j);
            }
            for (int k=0; k < 3; k++) for (int l=k+1; l < 3; l++) {
              *pa++ = r(k,i)*r(l,j) + r(k,j)*r(l,i) - delta_x_delta<T>(k,i,l,j);
            }
          }
        }

        /* Compute the matrix Z */
        af::ref<T, af::c_grid<2> > a_(a.get(), n_a_rows, n_all_params);
        scitbx::matrix::row_echelon::
          full_pivoting_small<T, max_n_a_rows, n_all_params>
          r_e(a_, pivot_zero_attractor);
        n_independent_params = r_e.nullity;
        af::small<T, n_all_params> x_N(n_independent_params, 0);
       af::mat_grid dim(n_all_params, n_independent_params);
        z = af::ref_owning_versa<T, af::mat_grid>(dim);
        for (unsigned j=0; j< n_independent_params; j++) {
          x_N[j] = 1;
          af::tiny<T, n_all_params>
            x_B = r_e.back_substitution(a_.begin(), x_N);
          for (unsigned i=0; i < n_all_params; i++) {
            z(i,j) = x_B[i];
          }
          x_N[j] = 0;
        }
     }

    public:
      cartesian_constraints()
      : pivot_zero_attractor(cartesian_constraints_constants::pivot_zero_attractor)
      {}

      cartesian_constraints(uctbx::unit_cell const & unit_cell,
                            space_group const& space_group,
                            double pivot_zero_attractor_ =
                             cartesian_constraints_constants::pivot_zero_attractor)
      : pivot_zero_attractor(pivot_zero_attractor_)
      {
        CCTBX_ASSERT(space_group.is_compatible_unit_cell(unit_cell));
        af::shared<sgtbx::rt_mx>
        point_symmetries = space_group.build_derived_acentric_group()
                                      .build_derived_point_group()
                                      .all_ops();
        /* The inversion is not part of those, hence the value of
           cartesian_constraints_symmetries_max_nb */
        initialise(unit_cell, point_symmetries.const_ref());
       }

      cartesian_constraints(uctbx::unit_cell const& unit_cell,
                            af::const_ref<rt_mx> const& matrices,
                            double pivot_zero_attractor_ =
                              cartesian_constraints_constants::pivot_zero_attractor)
      : pivot_zero_attractor(pivot_zero_attractor_)
      {
        initialise(unit_cell, matrices);
      }

      scitbx::sym_mat3<T>
      all_params(af::small<T, cartesian_constraints_constants::n_all_params>
                 const& independent_params) const
      {
        using namespace cartesian_constraints_constants;
        scitbx::sym_mat3<T> result;
        /* matrix-vector product z * independent_params,
         specialised for the memory arrangement of a sym_mat3
         */
        T* p = result.begin();
        for (unsigned i=0; i < n_all_params; i++, p++) {
          *p = 0;
          for (unsigned j=0; j < n_independent_params; j++) {
            *p += z(i,j)*independent_params[j];
          }
        }
        return result;
      }

      void fill_with_independent_gradients(T* p,
                                           scitbx::sym_mat3<T> const& all_grads) const
      {
        using namespace cartesian_constraints_constants;
        /* matrix-vector product z^T*  all_grads,
           specialised to store the result into the
           range [p, p+n_independent_params)
        */
        for (unsigned i=0; i < n_independent_params; i++, p++) {
          *p = 0;
          for (unsigned j=0; j < n_all_params; j++) {
            *p += z(j,i)*all_grads[j];
          }
        }
      }

      af::small<T, cartesian_constraints_constants::n_all_params>
      independent_gradients(scitbx::sym_mat3<T> const& all_grads) const
      {
        using namespace cartesian_constraints_constants;
        af::small<T, n_all_params> result(n_independent_params);
        fill_with_independent_gradients(result.begin(), all_grads);
        return result;
      }

      unsigned n_independent_parameters() const {
        return n_independent_params;
      }

      /// The Jacobian of the transformation from independent params to all params
      /** It reads
       \f[
            \left[ \frac{\partial u_i}{\partial v_j} \right]
       \f]
       where $u = (U_11, U_22, U_33, U_12, U_13, U_23)$, and
       where $v$ are the independent parameters $u$ is function of after
       solving for the special position constraints.

       */
      af::versa<T, af::mat_grid> jacobian() const { return z.array(); }
  };

}}} // namespace cctbx::sgtbx::tensor_rank_2

#endif // CCTBX_SGTBX_TENSOR_RANK_2
