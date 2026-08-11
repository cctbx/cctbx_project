from __future__ import absolute_import, division, print_function

"""
mlscale_step1 -- annotate observations with the geometry needed by the
maximum-likelihood scaling and merging worker (mlscale_step2).

WHAT THIS WORKER DOES

For every observation it records the small number of quantities from which the
partiality shape parameter

    rho = T / w          T = Ewald shell thickness, w = RLP width

can be reconstructed.  rho controls the shape of the observed intensity
distribution: rho << 1 is the thin-slice (monochromatic) limit where no
reflection is ever fully integrated, rho >> 1 is the full-integration limit.

WHY IT STORES INGREDIENTS RATHER THAN JUST rho

The global shape parameters are refined in an outer loop by step2.  If this
worker stored only rho, every outer iteration would require re-running it and
re-touching the experiment models.  Instead it stores the per-observation
quantities that do NOT depend on the refined parameters -- |q|, the shot
wavelength, and the shot spectral width -- so step2 can recompute rho on the
fly from whatever shape parameters it is currently trying.  rho itself is
written as well, from the phil starting values, and is advisory: use it for
diagnostics and for the initial table build, not as ground truth.

COLUMNS WRITTEN

    mlscale.frame       int    globally unique dense frame index across ranks
    mlscale.qlen        double |q| = 1/d, reciprocal Angstroms
    mlscale.wavelength  double per-shot wavelength, Angstroms
    mlscale.bandwidth   double per-shot relative spectral width, RMS
    mlscale.rho         double advisory shell-thickness / RLP-width ratio
    mlscale.eps         double excitation error (optional)
    mlscale.ewald_n     vec3   unit Ewald normal (optional, anisotropic model)

This worker never rejects reflections and never modifies intensities.
"""

import math
import numpy as np

from xfel.merging.application.worker import worker
from dials.array_family import flex

MLSCALE_COLUMNS = ['mlscale.frame', 'mlscale.qlen', 'mlscale.wavelength',
                   'mlscale.bandwidth', 'mlscale.rho', 'mlscale.eps',
                   'mlscale.ewald_n']


def _hkl_to_numpy(miller_index_column):
  """flex.miller_index -> (N, 3) float array, with fallbacks across versions."""
  try:
    return miller_index_column.as_vec3_double().as_double().as_numpy_array().reshape(-1, 3)
  except AttributeError:
    pass
  try:
    return np.asarray(miller_index_column.as_vec3_double()).reshape(-1, 3).astype(float)
  except Exception:
    return np.array([[float(a), float(b), float(c)]
                     for a, b, c in miller_index_column])


class rho_annotator(worker):
  """Annotate observations with the geometry behind the partiality shape parameter."""

  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(rho_annotator, self).__init__(params=params, mpi_helper=mpi_helper,
                                        mpi_logger=mpi_logger)
    # Make sure our columns survive any later prune step (e.g. the group
    # worker), whether or not the user remembered to list them in phil.
    if self.params.input.persistent_refl_cols is None:
      self.params.input.persistent_refl_cols = []
    for key in MLSCALE_COLUMNS:
      if key not in self.params.input.persistent_refl_cols:
        self.params.input.persistent_refl_cols.append(key)

  def __repr__(self):
    return "Annotate observations with partiality shape geometry"

  def validate(self):
    sh = self.params.mlscale.shape
    assert sh.domain_size > 0, "mlscale.shape.domain_size must be positive"
    assert sh.mosaicity >= 0 and sh.strain >= 0, \
        "mlscale.shape.mosaicity and .strain must be non-negative"
    assert sh.bandwidth >= 0 and sh.divergence >= 0, \
        "mlscale.shape.bandwidth and .divergence must be non-negative"
    if sh.mosaicity == 0 and sh.strain == 0 and sh.divergence == 0:
      raise ValueError("At least one of mlscale.shape.{mosaicity,strain,"
                       "divergence} must be nonzero, or the RLP width has no "
                       "resolution dependence and rho is degenerate with "
                       "domain_size at every resolution.")
    n_linear = sum(1 for v in (sh.mosaicity, sh.strain, sh.divergence) if v > 0)
    if n_linear > 1 and self.mpi_helper.rank == 0:
      self.logger.main_log(
          "WARNING: %d of {mosaicity, strain, divergence} are nonzero. These "
          "have identical q^1 resolution dependence and are not separable by "
          "an isotropic radial fit. Only one should be refined." % n_linear)

  # ---------------------------------------------------------------- helpers

  def shot_bandwidth(self, experiment):
    """Per-shot RMS relative spectral width, or None if unavailable."""
    if not self.params.mlscale.use_measured_spectrum:
      return None
    try:
      spectrum = experiment.beam.get_spectrum()
    except (AttributeError, RuntimeError):
      return None
    if spectrum is None:
      return None
    try:
      energy = np.asarray(spectrum.get_energies_eV(), dtype=float)
      weight = np.asarray(spectrum.get_weights(), dtype=float)
    except (AttributeError, RuntimeError):
      return None
    if energy.size < 2 or weight.sum() <= 0:
      return None
    weight = weight / weight.sum()
    mean = float((weight * energy).sum())
    if mean <= 0:
      return None
    var = float((weight * (energy - mean) ** 2).sum())
    rel = math.sqrt(max(var, 0.0)) / mean
    if not np.isfinite(rel) or rel <= 0:
      return None
    if rel > self.params.mlscale.spectrum_max_relative_width:
      return None
    return rel

  def shape_terms(self, qlen, wavelength, bandwidth):
    """Return (shell thickness T, RLP width w), both RMS, in reciprocal Angstroms.

    All contributions are RMS widths.  The O(1) factors connecting a physical
    size or angle to an RMS reciprocal-space width are absorbed into the
    fitted parameters, so these are effective quantities.
    """
    sh = self.params.mlscale.shape
    eta = math.radians(sh.mosaicity)

    # Ewald shell thickness.  Bandwidth displaces the sphere by
    # (sigma_E/E) * q^2 * lambda / 2; divergence tilts it, giving q * alpha.
    t_bandwidth = bandwidth * qlen ** 2 * wavelength / 2.0
    t_divergence = qlen * sh.divergence

    # Reciprocal lattice point width.
    w_size = np.full_like(qlen, 1.0 / sh.domain_size)
    w_mosaic = qlen * eta
    w_strain = qlen * sh.strain

    if sh.combine == 'quadrature':
      T = np.sqrt(t_bandwidth ** 2 + t_divergence ** 2)
      w = np.sqrt(w_size ** 2 + w_mosaic ** 2 + w_strain ** 2)
    else:
      T = t_bandwidth + t_divergence
      w = w_size + w_mosaic + w_strain
    return T, w

  # ---------------------------------------------------------------- main

  def run(self, experiments, reflections):
    self.logger.log_step_time("MLSCALE_ANNOTATE")

    if experiments is None or len(experiments) == 0:
      self.logger.log("No experiments on this rank; nothing to annotate")
      self.logger.log_step_time("MLSCALE_ANNOTATE", True)
      return experiments, reflections

    if 'miller_index' not in reflections:
      raise KeyError("mlscale_step1 needs the original 'miller_index' column "
                     "to compute scattering vectors. It must run before any "
                     "step that prunes it (i.e. before 'group').")

    comm = self.mpi_helper.comm
    MPI = self.mpi_helper.MPI

    # Globally unique dense frame index: offset this rank by the number of
    # experiments on all lower ranks.
    n_local = len(experiments)
    offset = comm.exscan(n_local, MPI.SUM)
    if offset is None:      # mpi4py returns None on rank 0
      offset = 0

    n_obs = reflections.size()
    frame = np.zeros(n_obs, dtype=np.int32)
    qlen = np.zeros(n_obs)
    wavel = np.zeros(n_obs)
    bwidth = np.zeros(n_obs)
    eps = np.zeros(n_obs) if self.params.mlscale.annotate.store_excitation_error else None
    normal = np.zeros((n_obs, 3)) if self.params.mlscale.annotate.store_ewald_normal else None

    ids = reflections['id'].as_numpy_array()
    n_no_spectrum = 0

    for expt_id, experiment in enumerate(experiments):
      sel = np.where(ids == expt_id)[0]
      if sel.size == 0:
        continue

      A = np.array(experiment.crystal.get_A(), dtype=float).reshape(3, 3)
      s0 = np.array(experiment.beam.get_s0(), dtype=float)
      wavelength = experiment.beam.get_wavelength()

      bw = self.shot_bandwidth(experiment)
      if bw is None:
        bw = self.params.mlscale.shape.bandwidth
        n_no_spectrum += 1

      h = _hkl_to_numpy(reflections['miller_index'].select(flex.size_t(sel)))
      q = h.dot(A.T)                       # q = A h, rows are observations
      s1 = q + s0[None, :]
      s1_len = np.linalg.norm(s1, axis=1)

      frame[sel] = offset + expt_id
      qlen[sel] = np.linalg.norm(q, axis=1)
      wavel[sel] = wavelength
      bwidth[sel] = bw
      if eps is not None:
        eps[sel] = s1_len - 1.0 / wavelength
      if normal is not None:
        normal[sel] = s1 / np.maximum(s1_len, 1e-12)[:, None]

    T, w = self.shape_terms(qlen, wavel, bwidth)
    rho = T / np.maximum(w, 1e-30)

    reflections['mlscale.frame'] = flex.int(frame.astype(np.int32))
    reflections['mlscale.qlen'] = flex.double(qlen)
    reflections['mlscale.wavelength'] = flex.double(wavel)
    reflections['mlscale.bandwidth'] = flex.double(bwidth)
    reflections['mlscale.rho'] = flex.double(rho)
    if eps is not None:
      reflections['mlscale.eps'] = flex.double(eps)
    if normal is not None:
      reflections['mlscale.ewald_n'] = flex.vec3_double(
          flex.double(normal.ravel()))

    self.log_diagnostics(qlen, rho, n_local, n_no_spectrum, bwidth)
    self.logger.log_step_time("MLSCALE_ANNOTATE", True)

    if self.params.mlscale.histograms.enable:
      self.run_histograms(reflections)

    return experiments, reflections

  # ---------------------------------------------------------------- histograms

  def global_hkl_counts(self, reflections):
    """Global multiplicity per asu HKL, reduced to rank 0.

    Counting is done on the asu index so that symmetry mates are pooled, which
    is the same grouping the merge uses.
    """
    comm = self.mpi_helper.comm
    key = ('miller_index_asymmetric' if 'miller_index_asymmetric' in reflections
           else 'miller_index')
    local = {}
    if reflections is not None and reflections.size() > 0:
      for hkl in reflections[key]:
        t = tuple(hkl)
        local[t] = local.get(t, 0) + 1
    parts = comm.gather(local, root=0)
    if self.mpi_helper.rank != 0:
      return None
    total = {}
    for part in parts:
      if not part:
        continue
      for h, c in part.items():
        total[h] = total.get(h, 0) + c
    return total

  def run_histograms(self, reflections):
    """Sample reflections over resolution and model intensity, then plot."""
    from xfel.merging.application.mlscale import histograms as hist_mod

    self.logger.log_step_time("MLSCALE_HISTOGRAMS")
    comm = self.mpi_helper.comm

    if 'i_model' not in (self.params.scaling).__dict__ or \
       self.params.scaling.i_model is None:
      if self.mpi_helper.rank == 0:
        self.logger.main_log(
            "mlscale histograms: scaling.i_model is not available. Put "
            "'model_scaling' earlier in dispatch.step_list and set "
            "scaling.model to a reference structure.")
      self.logger.log_step_time("MLSCALE_HISTOGRAMS", True)
      return

    counts = self.global_hkl_counts(reflections)

    selected = []
    if self.mpi_helper.rank == 0:
      selected = hist_mod.select_reflections(
          self.params, self.params.scaling.i_model, counts, self.logger)
    selected = comm.bcast(selected, root=0)

    if not selected:
      self.logger.log_step_time("MLSCALE_HISTOGRAMS", True)
      return

    observations = hist_mod.gather_observations(
        [s['hkl'] for s in selected], reflections, self.mpi_helper)

    if self.mpi_helper.rank == 0:
      hist_mod.make_plot(self.params, selected, observations, self.logger)

    self.logger.log_step_time("MLSCALE_HISTOGRAMS", True)

  # ---------------------------------------------------------------- logging

  def log_diagnostics(self, qlen, rho, n_local, n_no_spectrum, bwidth):
    """Report rho against resolution.

    p_max = erf(rho/2) is the largest partiality any observation in that shell
    can have.  Where p_max is short of 1 no reflection is ever fully
    integrated, the absolute intensity is degenerate with rho, and step2 will
    only recover intensities up to a per-shell scale factor.  That column is
    the single most useful thing in this log.
    """
    comm = self.mpi_helper.comm
    MPI = self.mpi_helper.MPI

    self.logger.log("Annotated %d observations over %d experiments (%d without "
                    "a usable spectrum)" % (qlen.size, n_local, n_no_spectrum))

    n_bins = self.params.statistics.n_bins
    edges = np.linspace(0.0, max(qlen.max() if qlen.size else 1.0, 1e-6), n_bins + 1)
    idx = np.clip(np.digitize(qlen, edges) - 1, 0, n_bins - 1)
    cnt = np.bincount(idx, minlength=n_bins).astype(float)
    sum_rho = np.bincount(idx, rho, minlength=n_bins)

    tot_cnt = comm.reduce(cnt, MPI.SUM, root=0)
    tot_rho = comm.reduce(sum_rho, MPI.SUM, root=0)
    tot_shots = comm.reduce(n_local, MPI.SUM, root=0)
    tot_nospec = comm.reduce(n_no_spectrum, MPI.SUM, root=0)
    bw_sum = comm.reduce(float(bwidth.sum()), MPI.SUM, root=0)
    bw_n = comm.reduce(int(bwidth.size), MPI.SUM, root=0)

    if self.mpi_helper.rank != 0:
      return

    sh = self.params.mlscale.shape
    self.logger.main_log("mlscale_step1: shape model")
    self.logger.main_log("  domain_size %.1f A   mosaicity %.4f deg   strain %.2e"
                         % (sh.domain_size, sh.mosaicity, sh.strain))
    self.logger.main_log("  bandwidth %.2e   divergence %.2e rad   combine=%s"
                         % (sh.bandwidth, sh.divergence, sh.combine))
    self.logger.main_log("  %d experiments; %d used the phil bandwidth because no "
                         "usable spectrum was attached" % (tot_shots, tot_nospec))
    if bw_n:
      self.logger.main_log("  mean per-observation bandwidth actually used: %.3e"
                           % (bw_sum / bw_n))

    self.logger.main_log("")
    self.logger.main_log("  %-18s %10s %10s %10s   %s"
                         % ("resolution (A)", "N obs", "mean rho", "p_max",
                            "regime"))
    for i in range(n_bins):
      if tot_cnt[i] < 1:
        continue
      q_lo, q_hi = edges[i], edges[i + 1]
      d_hi = 1.0 / q_lo if q_lo > 0 else float('inf')
      d_lo = 1.0 / q_hi if q_hi > 0 else float('inf')
      mean_rho = tot_rho[i] / tot_cnt[i]
      p_max = math.erf(mean_rho / 2.0)
      regime = "identifiable" if p_max > 0.999 else \
               ("marginal" if p_max > 0.95 else "rho degenerate with I")
      self.logger.main_log("  %8.2f - %-7.2f %10d %10.3f %10.4f   %s"
                           % (d_hi, d_lo, int(tot_cnt[i]), mean_rho, p_max, regime))
    self.logger.main_log("")
    self.logger.main_log("  Shells marked 'rho degenerate with I' cannot yield "
                         "absolute intensities from the intensity")
    self.logger.main_log("  distributions alone; their scale must be carried in "
                         "from the identifiable shells via the shape model.")


if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(rho_annotator)
