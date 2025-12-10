from __future__ import absolute_import, division, print_function
import os
import numpy as np
from xfel.merging.application.worker import worker
from dials.array_family import flex


class prepare_spread(worker):
  """
  Prepare data for SpReAD (Spectral Resolved Anomalous Diffraction) analysis.

  This worker:
  1. Computes global energy percentiles across all ranks
  2. Redistributes experiments/reflections by energy bin
  3. Writes binned data to disk for subsequent per-energy merging jobs
  """

  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(spread, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

  def __repr__(self):
    return "Prepare SpReAD energy-binned datasets"

  def run(self, experiments, reflections):
    self.logger.log_step_time("SPREAD_PREPARE")

    # Extract parameters (with defaults if not yet in phil)
    n_bins = getattr(self.params.prepare.spread, 'n_energy_bins', 100)
    output_dir = getattr(self.params.prepare.spread, 'output_dir',
                         self.params.output.output_dir)

    # Step 1: Compute local energies from wavelengths
    self.logger.log_step_time("COMPUTE_ENERGIES")
    local_energies, expt_to_energy = self._compute_local_energies(experiments)
    self.logger.log("Rank %d has %d experiments with energies" % (
      self.mpi_helper.rank, len(local_energies)))
    self.logger.log_step_time("COMPUTE_ENERGIES", True)

    # Step 2: Gather all energies to compute global percentile boundaries
    self.logger.log_step_time("COMPUTE_PERCENTILES")
    bin_edges = self._compute_global_percentile_edges(local_energies, n_bins)
    self.logger.log("Energy bin edges: %.2f to %.2f eV (%d bins)" % (
      bin_edges[0], bin_edges[-1], n_bins))
    self.logger.log_step_time("COMPUTE_PERCENTILES", True)

    # Step 3: Assign each experiment to a bin
    self.logger.log_step_time("ASSIGN_BINS")
    expt_to_bin = self._assign_experiments_to_bins(expt_to_energy, bin_edges, n_bins)
    self.logger.log_step_time("ASSIGN_BINS", True)

    # Step 4: Redistribute data by energy bin (all-to-all)
    self.logger.log_step_time("REDISTRIBUTE")
    redistributed_experiments, redistributed_reflections = self._redistribute_by_energy_bin(
      experiments, reflections, expt_to_bin, n_bins)
    self.logger.log("After redistribution: rank %d has %d experiments" % (
      self.mpi_helper.rank, len(redistributed_experiments)))
    self.logger.log_step_time("REDISTRIBUTE", True)

    # Step 5: Write output files
    self.logger.log_step_time("WRITE_FILES")
    self._write_binned_files(redistributed_experiments, redistributed_reflections,
                             expt_to_bin, n_bins, output_dir)
    self.logger.log_step_time("WRITE_FILES", True)

    self.logger.log_step_time("SPREAD_PREPARE", True)

    # Return empty - this worker is a terminal step that writes to disk
    return None, None

  def _compute_local_energies(self, experiments):
    """
    Extract energy (in eV) from each experiment's beam wavelength.
    Returns:
      local_energies: list of energies for all local experiments
      expt_to_energy: dict mapping experiment index to energy
    """
    # E = hc / lambda, with hc = 12398.419 eV*Angstrom
    hc = 12398.419
    local_energies = []
    expt_to_energy = {}

    for i, expt in enumerate(experiments):
      wavelength = expt.beam.get_wavelength()  # in Angstrom
      energy = hc / wavelength
      local_energies.append(energy)
      expt_to_energy[i] = energy

    return local_energies, expt_to_energy

  def _compute_global_percentile_edges(self, local_energies, n_bins):
    """
    Gather energies from all ranks and compute percentile bin edges.
    Returns array of n_bins+1 edges defining n_bins equal-count bins.
    """
    comm = self.mpi_helper.comm

    # Gather all energies to rank 0
    all_local_energies = comm.gather(local_energies, root=0)

    bin_edges = None
    if self.mpi_helper.rank == 0:
      # Flatten the list of lists
      all_energies = []
      for rank_energies in all_local_energies:
        all_energies.extend(rank_energies)
      all_energies = np.array(all_energies)

      # Compute percentile edges for equal-count bins
      percentiles = np.linspace(0, 100, n_bins + 1)
      bin_edges = np.percentile(all_energies, percentiles)

      self.logger.log("Total experiments across all ranks: %d" % len(all_energies))
      self.logger.log("Energy range: %.2f - %.2f eV" % (all_energies.min(), all_energies.max()))

    # Broadcast bin edges to all ranks
    bin_edges = comm.bcast(bin_edges, root=0)

    return bin_edges

  def _assign_experiments_to_bins(self, expt_to_energy, bin_edges, n_bins):
    """
    Assign each experiment to a percentile bin based on its energy.
    Returns dict mapping experiment index to bin index (0 to n_bins-1).
    """
    expt_to_bin = {}
    for expt_idx, energy in expt_to_energy.items():
      # np.searchsorted finds the bin; clip to valid range [0, n_bins-1]
      bin_idx = np.searchsorted(bin_edges[1:], energy, side='left')
      bin_idx = min(bin_idx, n_bins - 1)
      expt_to_bin[expt_idx] = bin_idx
    return expt_to_bin

  def _redistribute_by_energy_bin(self, experiments, reflections, expt_to_bin, n_bins):
    """
    Redistribute experiments and reflections so that each rank owns
    specific energy bins. Uses all-to-all communication pattern.

    Rank r will own bins where (bin_idx % n_ranks == r).
    """
    comm = self.mpi_helper.comm
    n_ranks = self.mpi_helper.size
    my_rank = self.mpi_helper.rank

    # Prepare send buffers: one list per destination rank
    send_expt_lists = [[] for _ in range(n_ranks)]
    send_refl_lists = [[] for _ in range(n_ranks)]

    # Build a mapping from old experiment id to new, and track which
    # reflections go to which rank
    refl_id_col = reflections['id']

    for expt_idx, expt in enumerate(experiments):
      bin_idx = expt_to_bin[expt_idx]
      dest_rank = bin_idx % n_ranks

      send_expt_lists[dest_rank].append((bin_idx, expt))

    # Now assign reflections based on their experiment's destination
    # First, build a quick lookup of expt_idx -> dest_rank
    expt_to_dest = {idx: expt_to_bin[idx] % n_ranks for idx in expt_to_bin}

    # Group reflections by destination rank
    for dest_rank in range(n_ranks):
      # Find all experiment indices going to this rank
      expt_indices_for_rank = [idx for idx, dest in expt_to_dest.items() if dest == dest_rank]
      if expt_indices_for_rank:
        # Select reflections belonging to these experiments
        sel = flex.bool(len(reflections), False)
        for idx in expt_indices_for_rank:
          sel |= (refl_id_col == idx)
        send_refl_lists[dest_rank] = reflections.select(sel)
      else:
        send_refl_lists[dest_rank] = flex.reflection_table()

    # All-to-all exchange
    recv_expt_lists = comm.alltoall(send_expt_lists)
    recv_refl_lists = comm.alltoall(send_refl_lists)

    # Consolidate received data
    # Store bin_idx with experiments for later file writing
    self._received_expt_bins = []  # list of (bin_idx, expt) tuples
    from dxtbx.model import ExperimentList
    new_experiments = ExperimentList()

    for rank_expts in recv_expt_lists:
      for bin_idx, expt in rank_expts:
        self._received_expt_bins.append((bin_idx, len(new_experiments)))
        new_experiments.append(expt)

    new_reflections = flex.reflection_table.concat(recv_refl_lists)

    return new_experiments, new_reflections

  def _write_binned_files(self, experiments, reflections, expt_to_bin, n_bins, output_dir):
    """
    Write experiments and reflections to files organized by energy bin and chunk.

    File naming: pct_{bin:02d}_c{chunk:02d}.expt / .refl
    where bin is the percentile bin (0 to n_bins-1) and chunk is derived from rank.
    """
    n_ranks = self.mpi_helper.size
    my_rank = self.mpi_helper.rank

    # Determine which bins this rank is responsible for
    my_bins = [b for b in range(n_bins) if b % n_ranks == my_rank]

    # Compute chunk index for each bin on this rank
    # chunk = rank // n_bins (how many complete cycles of bins fit before this rank)
    # But we need to be careful: chunk should be based on how many ranks share this bin
    # For bin b, the ranks that own it are: b, b+n_bins, b+2*n_bins, ... (if n_bins < n_ranks)
    # Or if n_ranks < n_bins: only rank (b % n_ranks) owns bin b

    # Group experiments by bin
    from dxtbx.model import ExperimentList
    bin_to_expts = {b: ExperimentList() for b in my_bins}
    bin_to_refls = {b: flex.reflection_table() for b in my_bins}

    # Map from received experiment index to bin
    expt_idx_to_bin = {expt_idx: bin_idx for bin_idx, expt_idx in self._received_expt_bins}

    # Separate experiments by bin
    for expt_idx, expt in enumerate(experiments):
      if expt_idx in expt_idx_to_bin:
        bin_idx = expt_idx_to_bin[expt_idx]
        bin_to_expts[bin_idx].append(expt)

    # Separate reflections by bin (based on experiment id)
    if len(reflections) > 0:
      refl_id_col = reflections['id']
      for bin_idx in my_bins:
        # Find experiment indices in this bin
        expt_indices = [expt_idx for expt_idx, b in expt_idx_to_bin.items() if b == bin_idx]
        if expt_indices:
          sel = flex.bool(len(reflections), False)
          for idx in expt_indices:
            sel |= (refl_id_col == idx)
          bin_to_refls[bin_idx] = reflections.select(sel)

    # Compute chunk index for this rank for each bin
    # If n_ranks >= n_bins: multiple ranks per bin, chunk = rank // n_bins
    # If n_ranks < n_bins: each rank handles multiple bins, chunk = 0
    if n_ranks >= n_bins:
      chunk_idx = my_rank // n_bins
    else:
      chunk_idx = 0

    # Write files
    os.makedirs(output_dir, exist_ok=True)

    for bin_idx in my_bins:
      expts = bin_to_expts[bin_idx]
      refls = bin_to_refls[bin_idx]

      if len(expts) == 0:
        continue

      # Format bin index width based on n_bins
      bin_width = len(str(n_bins - 1))
      filename_base = f"pct_{bin_idx:0{bin_width}d}_c{chunk_idx:02d}"

      expt_path = os.path.join(output_dir, f"{filename_base}.expt")
      refl_path = os.path.join(output_dir, f"{filename_base}.refl")

      # Renumber experiment ids to be contiguous before writing
      refls.reset_ids()

      expts.as_file(expt_path)
      refls.as_file(refl_path)

      self.logger.log(f"Wrote {len(expts)} experiments to {filename_base}.expt/.refl")