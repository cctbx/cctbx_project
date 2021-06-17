/* -*- mode: c++; c-basic-offset: 2; indent-tabs-mode: nil; tab-width: 8 -*-
 *
 * $Id$
 */

typedef scitbx::af::shared<int> shared_int;
typedef scitbx::af::shared<double> shared_double;
typedef scitbx::af::versa<bool, scitbx::af::flex_grid<> > shared_bool;


/* column_parser cannot be const
 *
 * @param w Per-observation weights.  All elements of @p w must be
 *          <code>&ge; 0</code>.
 */
static boost::python::tuple
compute_functional_and_gradients(const shared_double& x,
                                 const shared_double& w,
                                 std::size_t n_frames,
                                 column_parser& observations)
{
  const shared_double data = observations.get_double("i");
  const shared_double sigmas = observations.get_double("sigi");
  const shared_int frames = observations.get_int("frame_id");
  const shared_int hkl = observations.get_int("hkl_id");

  const std::size_t n_obs = w.size();

  SCITBX_ASSERT(data.size() == n_obs);
  SCITBX_ASSERT(sigmas.size() == n_obs);
  SCITBX_ASSERT(frames.size() == n_obs);
  SCITBX_ASSERT(hkl.size() == n_obs);

  shared_double g(x.size());
  double f = 0;
  double w_tot = 0;

  for (std::size_t i = 0; i < n_obs; i++) {
    if (w[i] <= 0)
      continue;

    const int h = hkl[i];
    const double I_m = x[n_frames + h];
    if (I_m <= 0) // XXX This warrants a comment
      continue;

    const int m = frames[i];
    const double G = x[m];
    const double t = data[i] - G * I_m;

    f += w[i] * t * t;
    g[m] += -2.0 * w[i] * I_m * t;
    g[n_frames + h] += -2.0 * w[i] * G * t;
    w_tot += w[i];
  }

  /* If w_tot == 0, then f and g should be identically zero.  w_tot <
   * 0 should not happen.
   */
  if (w_tot > 0) {
    f /= w_tot;
    for (std::size_t i = 0; i < g.size(); i++)
      g[i] /= w_tot;
  }

  return (boost::python::make_tuple(f, g));
}


/*
 * The curvatures() function returns the diagonal elements of the
 * Hessian matrix.
 *
 * If we knew that compute_functional_and_gradients() was always
 * called before curvatures(), we could calculate the curvatures there
 * and just return them here.
 */
shared_double
curvatures(const shared_double& x,
           const shared_double& w,
           std::size_t n_frames,
           column_parser& observations)
{
  const shared_double data = observations.get_double("i");
  const shared_double sigmas = observations.get_double("sigi");
  const shared_int frames = observations.get_int("frame_id");
  const shared_int hkl = observations.get_int("hkl_id");

  const std::size_t n_obs = w.size();

  SCITBX_ASSERT(data.size() == n_obs);
  SCITBX_ASSERT(sigmas.size() == n_obs);
  SCITBX_ASSERT(frames.size() == n_obs);
  SCITBX_ASSERT(hkl.size() == n_obs);

  shared_double c(x.size());
  double w_tot = 0;

  printf("CURVATURES CALLED\n");

  for (std::size_t i = 0; i < n_obs; i++) {
    if (w[i] <= 0)
      continue;

    if (data[i] <= 0)
      continue;

    const int h = hkl[i];
    const double I_m = x[n_frames + h];
    if (I_m <= 0) // XXX See compute_functional_and_gradients()
      continue;

    const int m = frames[i];
    const double G = x[m];

    c[m] += 2.0 * w[i] * I_m * I_m;
    c[n_frames + h] += 2.0 * w[i] * G * G;
    w_tot += w[i];
  }

  for (std::size_t i = 0; i < c.size(); i++)
    c[i] /= w_tot;

  /*
  for (std::size_t i = 0; i < c.size(); i++)
    if (c[i] <= 0)
      printf("CURVATURES index %d is %f\n", i, c[i]);
  */

  return (c);
}


/* The intensity column in the MTZ file written by
 * finalize_and_save_data() in cxi_merge.scaling_manager is calculated
 * as
 *
 *   summed_wt_I / summed_weight.
 *
 * The corresponding standard deviations are
 *
 *   sqrt(1 / summed_weight).
 *
 * In order to ensure that the refined model intensities are
 * reproduced in the final intensity column, this function populates
 * summed_wt_I with the model intensities multiplied by their
 * corresponding weight.
 *
 * XXX Duplicates code w.r.t. other stuff defined in this file.  It
 * should ideally be possible to reuse get_scaling_results() below!
 *
 * XXX Needs to be cleaned up!
 */
static boost::python::tuple
get_scaling_results_mark2(const shared_double& x,
                          const shared_double& w,
                          scaling_results& L,
                          const cctbx::uctbx::unit_cell& params_unit_cell)
{
  const shared_double data = L.observations.get_double("i");
  const shared_double sigmas = L.observations.get_double("sigi");

  // slope is only needed for mark0 comparison.
  const shared_double slope = L.frames.get_double("slope");

  const shared_int frames = L.observations.get_int("frame_id");
  const shared_int hkl = L.observations.get_int("hkl_id");

  const std::size_t n_frames = L.selected_frames.size();
  const std::size_t n_hkl = x.size() - n_frames;
  const std::size_t n_obs = w.size();

  SCITBX_ASSERT(data.size() == n_obs);
  SCITBX_ASSERT(sigmas.size() == n_obs);
  SCITBX_ASSERT(frames.size() == n_obs);
  SCITBX_ASSERT(hkl.size() == n_obs);

  /* Weighted sum of squared deviations between scaled, observed
   * intensities and refined intensities.
   */
  shared_double wssq_dev(n_hkl);
  shared_double w_tot(n_hkl);

  for (std::size_t m = 0; m < n_frames; m++) {
    L.n_rejected[m] = 0;
    L.n_obs[m] = 0;
    L.d_min_values[m] = 0;
  }

  for (std::size_t i = 0; i < n_hkl; i++) {
    L.sum_I[i] = 0;
    L.sum_I_SIGI[i] = 0;
    L.completeness[i] = 0;
    L.summed_N[i] = 0;
    L.summed_wt_I[i] = 0;
    L.summed_weight[i] = 0;
  }

  L.i_isig_list.clear();

  for (std::size_t i = 0; i < n_obs; i++) {
    const int m = frames[i];

    /* Filtering the observations based on weight early in the loop
     * alters the statistics.  Previously, a reflection was counted as
     * observed and would contribute to the completeness as long as
     * the frame was not rejected, even if its integrated intensity or
     * standard deviation was non-positive.  Now, any observations
     * with non-positive weight are not counted as observed and do not
     * contribute to the completeness.
     */
    if (w[i] <= 0) {
      L.n_rejected[m] += 1;
      continue;
    }

    const int h = hkl[i];
    L.completeness[h] += 1;
    L.n_obs[m] += 1;

    // Mind negative scaling factors!
    const double G = x[m];
    //const double G = slope[m]; // To reproduce mark0 scaling
    //const double G = 1; // To reproduce mark1 scaling

    /* The integrated, observed intensity, I, and its standard
     * deviation, s, scaled by the per-frame scale factor, G.
     */
    const double I = data[i] / G;
    const double s = sigmas[i] / G;
    const double IsigI = I / s;

    L.summed_N[h] += 1;
    L.sum_I[h] += I;
    L.sum_I_SIGI[h] += IsigI;

    L.hkl_ids.push_back(h);
    L.i_isig_list.push_back(scitbx::vec3<double>(I, IsigI, G));

    const double d =
      params_unit_cell.d(cctbx::miller::index<>(L.merged_asu_hkl[h]));
    if (L.d_min_values[m] <= 0 || L.d_min_values[m] > d)
      L.d_min_values[m] = d;

    /* Running, weighted, and squared sum of deviations between
     * scaled, observed intensities and refined dittos.  Note that
     * neither L.summed_weight nor L.summed_wt_I are written in this
     * loop.
     *
     * Leave ESTIMATE_REFINED_UNCERTAINTY undefined to propagate the
     * shot noise in quadrature to the merged sig(I).  The
     * alternative, i.e. estimating the refined uncertainty, is
     * currently broken.
     *
     * XXX Must ensure w_tot > 0 in order to preserve logic in the
     * looped branch below.
     */
//#define ESTIMATE_REFINED_UNCERTAINTY 1
#ifdef ESTIMATE_REFINED_UNCERTAINTY
    const double w_this = G * G * w[i]; // XXX ugly!
    const double t = I - x[n_frames + h];
    wssq_dev[h] += w_this * t * t;
    w_tot[h] += w_this;
#else
    wssq_dev[h] += 1.0 / (s * s);
    w_tot[h] += 1;
#endif
  }

  /* Multiply model intensities by the summed weight, such that the
   * intensities written to the final MTZ-file will be correct.
   */
  for (std::size_t h = 0; h < n_hkl; h++) {
    if (w_tot[h] > 0 && x[n_frames + h] > 0) {
#ifdef ESTIMATE_REFINED_UNCERTAINTY
      L.summed_weight[h] = wssq_dev[h] / w_tot[h];
#else
      L.summed_weight[h] = wssq_dev[h];
#endif
      L.summed_wt_I[h] =
        L.summed_weight[h] > 0 ? x[n_frames + h] * L.summed_weight[h] : 0;

    } else {
      L.summed_weight[h] = 0;
      L.summed_wt_I[h] = 0;
    }
  }

  return (make_tuple(L.sum_I,
                     L.sum_I_SIGI,
                     L.completeness,
                     L.summed_N,
                     L.summed_wt_I,
                     L.summed_weight,
                     L.n_rejected,
                     L.n_obs,
                     L.d_min_values,
                     L.hkl_ids,
                     L.i_isig_list));
}
