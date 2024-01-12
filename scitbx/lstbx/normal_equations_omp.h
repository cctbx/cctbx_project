// this file is included into <scitbx/lstbx/normal_equations.h>

template<>
void non_linear_ls_with_separable_scale_factor<
  double,
  matrix::sum_of_symmetric_rank_1_updates>::add_residuals_omp(
    const int& n,
    const int& start,
    const int& threads,
    af::const_ref<scalar_t> const& yc,
    af::const_ref<scalar_t> const& yo,
    af::const_ref<scalar_t> const& w)
{
  n_data += n;
  scalar_t temp1;
  scalar_t temp2 = 0;
  scalar_t temp3 = 0;
  scalar_t temp4 = 0;
#pragma omp parallel for reduction(+:temp2, temp3, temp4) num_threads(threads) private(temp1)
  for (int i = start; i < start + n; i++) {
    temp1 = w[i] * yo[i];
    temp2 += temp1 * yo[i];
    temp3 += temp1 * yc[i];
    temp4 += w[i] * yc[i] * yc[i];
  }
  yo_sq += temp2;
  yo_dot_yc += temp3;
  yc_sq += temp4;
}
//.............................................................................
template<>
void non_linear_ls_with_separable_scale_factor<
  double,
  matrix::sum_of_symmetric_rank_1_updates> ::add_residuals_omp(const int& n,
  const int& start,
  const int& threads,
  af::const_ref<scalar_t> const& yc,
  af::const_ref<scalar_t> const& yo)
{
  n_data += n;
  scalar_t temp1 = 0;
  scalar_t temp2 = 0;
  scalar_t temp3 = 0;
#pragma omp parallel for reduction(+:temp1, temp2, temp3) num_threads(threads)
  for (int i = start; i < start + n; i++) {
    temp1 += yo[i] * yo[i];
    temp2 += yo[i] * yc[i];
    temp3 += yc[i] * yc[i];
  }
  yo_sq += temp1;
  yo_dot_yc += temp2;
  yc_sq += temp3;
}
//.............................................................................
template<>
void non_linear_ls_with_separable_scale_factor<
  double,
  matrix::sum_of_symmetric_rank_1_updates> ::add_equations_omp(const int& n_ref,
  const int& n_par,
  const int& chunk_size,
  const int& start,
  const int& threads,
  std::vector<double>& matrix,
  /* These are arrays passed as locals for each thread to reduce overhead of
  generating them for each chunk */
  std::vector<double>& yo_dot_grad_yc_,
  std::vector<double>& yc_dot_grad_yc_,
  af::const_ref<scalar_t> const& yc,
  std::vector<double> const& jacobian_yc,
  af::const_ref<scalar_t> const& yo,
  af::const_ref<scalar_t> const& w)
{
  typedef double FloatType;
  SCITBX_ASSERT(yc.size() == n_ref
    && (!w.size() || yc.size() == w.size())
    && yo.size() == n_ref)
    (yc.size())(n_ref)(w.size())(yo.size());
  SCITBX_ASSERT(jacobian_yc.size() == chunk_size * n_par)
    (jacobian_yc.size())(chunk_size)(n_par);
  SCITBX_ASSERT(!finalised());
  FloatType* m = symmetric_matrix_owning_ref_t(grad_yc_dot_grad_yc).array().begin();
  const int limit = matrix.size() / threads;
  const int stepsize = 30;
  if (w.size()) {
    n_data += chunk_size;
    FloatType temp2 = 0;
    FloatType temp3 = 0;
    FloatType temp4 = 0;
#pragma omp parallel num_threads(threads) shared(temp2,temp3,temp4)
    {
      FloatType temp1;
      const int thread = omp_get_thread_num();
      FloatType* l_matrix = &(matrix[thread * limit]);
      /* Take just the piece of the given array that is required for this thread */
      FloatType* l_ogc = &(yo_dot_grad_yc_[thread * n_par]);
      FloatType* l_cgc = &(yc_dot_grad_yc_[thread * n_par]);
#pragma omp for reduction(+:temp2, temp3, temp4)
      for (int i = start; i < start + chunk_size; i++) {
        temp1 = w[i] * yo[i];
        temp2 += temp1 * yo[i];
        temp3 += temp1 * yc[i];
        temp4 += w[i] * yc[i] * yc[i];
      }
#pragma omp single
      {
        yo_sq += temp2;
        yo_dot_yc += temp3;
        yc_sq += temp4;
      }
      for (int i = start; i < start + chunk_size; ++i) {
        const FloatType* g_yc_loc = &(jacobian_yc[(i - start) * n_par]);
        /* Doing this backwards apparently increases data locality and outperforms
        forwards */
#pragma omp for nowait schedule(static,1)
        for (int x = n_par - 1; x >= 0; x--) {
          FloatType alpha_x = w[i] * g_yc_loc[x];
          l_ogc[x] += alpha_x * yo[i];
          l_cgc[x] += alpha_x * yc[i];
          int run = x * (n_par - 1) - x * (x - 1) / 2;
          for (int y = n_par - 1; y >= x; y--) {
            l_matrix[run + y] += alpha_x * g_yc_loc[y];
          }
        }
      }
#pragma omp critical
      {
        for (int i = 0; i < limit; i++) {
          m[i] += l_matrix[i];
        }
        for (int i = 0; i < n_par; i++) {
          yo_dot_grad_yc[i] += l_ogc[i];
          yc_dot_grad_yc[i] += l_cgc[i];
        }
      }
    }
  }
  else {
    n_data += chunk_size;
    scalar_t temp2 = 0;
    scalar_t temp3 = 0;
    scalar_t temp4 = 0;
#pragma omp parallel num_threads(threads) shared(temp2,temp3,temp4)
    {
      const int thread = omp_get_thread_num();
      FloatType* l_matrix = &(matrix[thread * limit]);
      /* Take just the piece of the given array that is required for this thread */
      FloatType* l_ogc = &(yo_dot_grad_yc_[thread * n_par]);  
      FloatType* l_cgc = &(yc_dot_grad_yc_[thread * n_par]);
#pragma omp for reduction(+:temp2, temp3, temp4)
      for (int i = start; i < start + chunk_size; i++) {
        temp2 += yo[i] * yo[i];
        temp3 += yo[i] * yc[i];
        temp4 += w[i] * yc[i] * yc[i];
      }
#pragma omp single
      {
        yo_sq += temp2;
        yo_dot_yc += temp3;
        yc_sq += temp4;
      }
      for (int i = start; i < start + chunk_size; ++i) {
        const double* g_yc_loc = &(jacobian_yc[(i - start) * n_par]);
#pragma omp for nowait schedule(static,1)
        for (int x = n_par - 1; x >= 0; x--) {
          l_ogc[x] += g_yc_loc[x] * yo[i];
          l_cgc[x] += g_yc_loc[x] * yc[i];
          int run = x * (n_par - 1) - x * (x - 1) / 2;
          for (int y = n_par - 1; y >= 0; y++) {
            l_matrix[run + y] += g_yc_loc[x] * g_yc_loc[y];
          }
        }
      }
#pragma omp critical
      {
        for (int i = 0; i < limit; i++) {
          m[i] += l_matrix[i];
        }
        for (int i = 0; i < n_par; i++) {
          yo_dot_grad_yc[i] += l_ogc[i];
          yc_dot_grad_yc[i] += l_cgc[i];
        }
      }
    }
  }
}
//.............................................................................