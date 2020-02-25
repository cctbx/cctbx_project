#ifndef CCTBX_SGTBX_TENSORS
#define CCTBX_SGTBX_TENSORS
#include <cctbx/sgtbx/space_group.h>
#include <scitbx/matrix/delta_tensors.h>
#include <scitbx/matrix/row_echelon.h>
#include <scitbx/matrix/tensors.h>
#include <scitbx/matrix/row_echelon_full_pivoting_small.h>
#include <scitbx/array_family/owning_ref.h>
#include <boost/shared_array.hpp>

namespace cctbx { namespace sgtbx { namespace tensors {

template <typename FloatType, class tensor_t>
class constraints {
protected:
  boost::shared_array<int> row_echelon_form_memory;
  mutable boost::shared_array<FloatType> gradient_sum_matrix_memory;
public:
  af::shared<unsigned> independent_indices;

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
      row_echelon_form_memory.get(),
      tensor_t::size() - independent_indices.size(),
      tensor_t::size());
  }

  std::size_t n_independent_params() const { return independent_indices.size(); }

  std::size_t n_dependent_params() const {
    return tensor_t::size() - independent_indices.size();
  }

  af::shared<FloatType>
    independent_params(af::shared<FloatType> const& all_params) const
  {
    af::shared<FloatType> result;
    for (std::size_t i = 0; i < independent_indices.size(); i++) {
      result.push_back(all_params[independent_indices[i]]);
    }
    return result;
  }

  af::shared<FloatType> all_params(
    af::shared<FloatType> const& independent_params) const
  {
    af::shared<FloatType> result(tensor_t::size());
    for (std::size_t i = 0; i < independent_params.size(); i++) {
      result[independent_indices[i]] = independent_params[i];
    }
    scitbx::matrix::row_echelon::back_substitution_float(
      row_echelon_form(),
      static_cast<const FloatType*>(0),
      result.begin());
    return result;
  }

  //! Constraint matrix (with a historical name).
  af::const_ref<FloatType, af::mat_grid>
    gradient_sum_matrix() const
  {
    const FloatType* gsm = gradient_sum_matrix_memory.get();
    if (gsm == 0) {
      gsm = initialize_gradient_sum_matrix();
    }
    return af::const_ref<FloatType, af::mat_grid>(
      gsm, independent_indices.size(), tensor_t::size());
  }

  af::shared<FloatType> independent_gradients(
    scitbx::sym_mat3<FloatType> const& all_gradients) const
  {
    af::shared<FloatType> result;
    const FloatType* row = gradient_sum_matrix_memory.get();
    if (row == 0) {
      row = initialize_gradient_sum_matrix();
    }
    for (std::size_t i = 0; i < independent_indices.size();
      i++, row += tensor_t::size())
    {
      FloatType indep_grad = 0;
      for (std::size_t j = 0; j < tensor_t::size(); j++) {
        indep_grad += row[j] * all_gradients[j];
      }
      result.push_back(indep_grad);
    }
    return result;
  }

  af::shared<FloatType>
    independent_curvatures(
      af::const_ref<FloatType> const& all_curvatures) const
  {
    CCTBX_ASSERT(all_curvatures.size() ==
      tensor_t::size() * (tensor_t::size() + 1) / 2);
    const FloatType* gsm = gradient_sum_matrix_memory.get();
    if (gsm == 0) {
      gsm = initialize_gradient_sum_matrix();
    }
    unsigned n_indep = n_independent_params();
    af::shared<FloatType> result(
      n_indep*(n_indep + 1) / 2, af::init_functor_null<FloatType>());
    boost::scoped_array<FloatType> c_times_a(
      new FloatType[n_indep * tensor_t::size()]);
    scitbx::matrix::multiply_packed_u_multiply_lhs_transpose(
      gsm, all_curvatures.begin(), n_indep, tensor_t::size(),
      c_times_a.get(), result.begin());
    return result;
  }

  static const std::vector<std::vector<int> > &get_indices() {
    return tensor_t::get_indices();
  }

  /* this might be required in a multithreaded environment!
  allocate a small cunk of memory for the multidemensional index to
  linear index mapping
  */
  static void initialise() {
    tensor_t::initialise();
  }

  /* release any allocated memory, not really needed unless want to
  investigate potential memory leaks
  */
  static void cleanup() {
    tensor_t::cleanup();
  }

protected:
  void
    initialize(
      af::const_ref<rt_mx> const& symmetry_matrices,
      std::size_t i_first_matrix_to_use,
      bool reciprocal_space)
  {
    CCTBX_ASSERT(i_first_matrix_to_use <= symmetry_matrices.size());
    unsigned n_rows = (symmetry_matrices.size() - i_first_matrix_to_use) * tensor_t::size();
    boost::shared_array<int> row_echelon_setup_memory(new int[n_rows * tensor_t::size()]);
    af::ref<int, af::mat_grid> row_echelon_setup(
      row_echelon_setup_memory.get(), n_rows, tensor_t::size());
    using scitbx::matrix::tensors::tensor_rank_2;
    const std::vector<std::vector<int> > &indices = tensor_t::get_indices();
    for (unsigned i_s = i_first_matrix_to_use; i_s < symmetry_matrices.size(); i_s++) {
      const rot_mx rm = reciprocal_space ? symmetry_matrices[i_s].r()
        : symmetry_matrices[i_s].r().transpose();
      int *m_off = &row_echelon_setup.begin()[
        (i_s - i_first_matrix_to_use)*tensor_t::size()*tensor_t::size()];
      for (size_t i = 0; i < tensor_t::size(); i++) {
        af::shared<FloatType> t = tensor_t::get_transform(indices[i], rm);
        size_t li = tensor_t::get_linear_idx(indices[i]);
        t[li] -= 1;
        int *t_off = &m_off[li*tensor_t::size()];
        for (int j = 0; j < tensor_t::size(); j++) {
          t_off[j] = t[j];
        }
      }
    }
    n_rows = scitbx::matrix::row_echelon::form(row_echelon_setup);
    CCTBX_ASSERT(n_rows <= tensor_t::size());
    row_echelon_form_memory = boost::shared_array<int>(new int[n_rows * tensor_t::size()]);
    std::copy(
      row_echelon_setup.begin(),
      row_echelon_setup.end(),
      row_echelon_form_memory.get());
    boost::shared_array<bool> independent_flags(new bool[tensor_t::size()]);
    scitbx::matrix::row_echelon::back_substitution_int(
      row_echelon_setup,
      static_cast<const int*>(0),
      static_cast<int*>(0),
      independent_flags.get());
    for (size_t i = 0; i < tensor_t::size(); i++) {
      if (independent_flags[i]) {
        independent_indices.push_back(i);
      }
    }
  }

  const FloatType*
    initialize_gradient_sum_matrix() const
  {
    gradient_sum_matrix_memory = boost::shared_array<FloatType>(
      new FloatType[independent_indices.size() * tensor_t::size()]);
    FloatType* row = gradient_sum_matrix_memory.get();
    std::fill_n(row,
      independent_indices.size() * tensor_t::size(),
      static_cast<FloatType>(0));
    af::const_ref<int, af::mat_grid> rem = row_echelon_form();
    for (std::size_t i = 0; i < independent_indices.size();
      i++, row += tensor_t::size())
    {
      row[independent_indices[i]] = 1;
      scitbx::matrix::row_echelon::back_substitution_float(
        rem, static_cast<const FloatType*>(0), row);
    }
    return gradient_sum_matrix_memory.get();
  }

};

}}} // namespace cctbx::sgtbx::tensors

#endif //CCTBX_SGTBX_TENSORS
