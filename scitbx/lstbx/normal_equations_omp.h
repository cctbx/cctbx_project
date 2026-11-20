// this file is included into <scitbx/lstbx/normal_equations.h>

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
      /* The scratch is allocated by the caller and carried from one chunk to
         the next, so by the last it holds the sum over all of them. Reducing
         it into the normal equations after every chunk would therefore add the
         earlier chunks again, once more for each chunk that follows.
       */
      if (start + chunk_size >= n_ref) {
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
        // this is the branch taken when there are no weights, so w is empty
        // and w[i] was reading past the end of it
        temp4 += yc[i] * yc[i];
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
          /* counting down to x, as the weighted branch above does: the packed
             row holds only the upper triangle, and y++ from n_par-1 against
             y >= 0 neither terminates nor stays inside it
           */
          for (int y = n_par - 1; y >= x; y--) {
            l_matrix[run + y] += g_yc_loc[x] * g_yc_loc[y];
          }
        }
      }
      // see the weighted branch above for why this waits for the last chunk
      if (start + chunk_size >= n_ref) {
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
}
//.............................................................................

/* The BLAS 3 accumulator needs none of the per-thread matrices above: it takes
the scaled rows and folds a whole chunk of them in with one sfrk, which OpenBLAS
threads for itself. So all that is parallelised here is filling its buffer and
the two vector sums; the matrix work is left to the BLAS, which is better at it.
*/
template<>
void non_linear_ls_with_separable_scale_factor<
  double,
  matrix::rank_n_update>::add_equations_omp(const int& n_ref,
  const int& n_par,
  const int& chunk_size,
  const int& start,
  const int& threads,
  std::vector<double>& matrix,      // unused: see omp_matrix_scratch
  std::vector<double>& yo_dot_grad_yc_,
  std::vector<double>& yc_dot_grad_yc_,
  af::const_ref<scalar_t> const& yc,
  std::vector<double> const& jacobian_yc,
  af::const_ref<scalar_t> const& yo,
  af::const_ref<scalar_t> const& w)
{
  SCITBX_ASSERT(yc.size() == n_ref
    && (!w.size() || yc.size() == w.size())
    && yo.size() == n_ref)
    (yc.size())(n_ref)(w.size())(yo.size());
  SCITBX_ASSERT(jacobian_yc.size() == chunk_size * n_par)
    (jacobian_yc.size())(chunk_size)(n_par);
  SCITBX_ASSERT(!finalised());
  bool const weighted = w.size() != 0;

  n_data += chunk_size;
  scalar_t temp2 = 0, temp3 = 0, temp4 = 0;
#pragma omp parallel for reduction(+:temp2, temp3, temp4) num_threads(threads)
  for (int i = start; i < start + chunk_size; i++) {
    scalar_t const w_i = weighted ? w[i] : 1;
    scalar_t const temp1 = w_i * yo[i];
    temp2 += temp1 * yo[i];
    temp3 += temp1 * yc[i];
    temp4 += w_i * yc[i] * yc[i];
  }
  yo_sq += temp2;
  yo_dot_yc += temp3;
  yc_sq += temp4;

  /* yo . grad yc and yc . grad yc, one running copy per thread. Each thread
  only ever touches its own, so the loop over reflections may be shared out
  freely; they are reduced into the equations on the last chunk, the scratch
  being carried across chunks by the caller.
  */
#pragma omp parallel num_threads(threads)
  {
    const int thread = omp_get_thread_num();
    double* l_ogc = &(yo_dot_grad_yc_[thread * n_par]);
    double* l_cgc = &(yc_dot_grad_yc_[thread * n_par]);
#pragma omp for
    for (int i = start; i < start + chunk_size; i++) {
      const double* g = &(jacobian_yc[(i - start) * n_par]);
      double const w_i = weighted ? w[i] : 1;
      double const a_yo = w_i * yo[i], a_yc = w_i * yc[i];
      for (int x = 0; x < n_par; x++) {
        l_ogc[x] += g[x] * a_yo;
        l_cgc[x] += g[x] * a_yc;
      }
    }
  }
  if (start + chunk_size >= n_ref) {
    for (int t = 0; t < threads; t++) {
      const double* l_ogc = &(yo_dot_grad_yc_[t * n_par]);
      const double* l_cgc = &(yc_dot_grad_yc_[t * n_par]);
      for (int x = 0; x < n_par; x++) {
        yo_dot_grad_yc[x] += l_ogc[x];
        yc_dot_grad_yc[x] += l_cgc[x];
      }
    }
  }

  /* The rows, written into the accumulator's buffer in place. Filling it takes
  as many passes as the buffer has room for at a time; folding a full one in is
  a collective sfrk and so happens between them, not inside a parallel region.
  */
  int written = 0;
  while (written < chunk_size) {
    if (grad_yc_dot_grad_yc.n_pending() == grad_yc_dot_grad_yc.max_pending()) {
      grad_yc_dot_grad_yc.fold_pending();
    }
    int const room = static_cast<int>(
      grad_yc_dot_grad_yc.max_pending() - grad_yc_dot_grad_yc.n_pending());
    int const take = std::min(chunk_size - written, room);
    double* rows = grad_yc_dot_grad_yc.open_pending(take);
#pragma omp parallel for num_threads(threads)
    for (int k = 0; k < take; k++) {
      int const i = start + written + k;
      const double* src = &(jacobian_yc[(i - start) * n_par]);
      double* dst = rows + std::size_t(k) * n_par;
      // copy then scale, not one fused loop: see rank_n_update::add
      std::copy(src, src + n_par, dst);
      matrix::scale_vector(n_par, dst, std::sqrt(weighted ? w[i] : 1.0));
    }
    written += take;
  }
}
//.............................................................................
