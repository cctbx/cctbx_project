from __future__ import absolute_import, division, print_function
import os
import numpy as np
from xfel.merging.application.worker import worker
from dials.array_family import flex


# Template for individual slice phil files
SLICE_PHIL_TEMPLATE = """\
# Base phil from: {stage2_phil_path}
{base_phil_content}

# Slice {slice_idx:03d}: percentiles {pct_start}-{pct_end}, center {center}
# Energy range: {energy_start:.2f} - {energy_end:.2f} eV

{input_paths}
input.experiments_suffix=.expt
input.reflections_suffix=.refl
output.output_dir={slice_output_dir}
"""

# Template for the unified batch script (SLURM + local)
SPREAD_SCRIPT_TEMPLATE = """\
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node={stage2_nproc}
{slurm_optional_directives}#SBATCH --array=0-{n_slices}%{slurm_array_concurrency}
#SBATCH --time={slurm_time_limit}
#SBATCH --job-name=spread
#SBATCH --output={scripts_dir}/slurm_%A_%a.out
#SBATCH --error={scripts_dir}/slurm_%A_%a.err

# This script works in two modes:
# 1. SLURM: sbatch run_spread.sh (each array task runs one iteration)
# 2. Local: bash run_spread.sh (runs all tasks sequentially)

N_SLICES={n_slices}
SCRIPTS_DIR="{scripts_dir}"
STATUS_DIR="{status_dir}"
STAGE2_DIR="{stage2_output_dir}"
MTZ_NAME="{mtz_name}"
NPROC={stage2_nproc}
PDB_FILE="{phenix_pdb_path}"
PHENIX_PHIL="{phenix_phil_path}"
CCTBX_ACTIVATE="{cctbx_activate}"
PHENIX_ACTIVATE="{phenix_activate}"
N_SCATTERERS={n_anomalous_scatterers}

mkdir -p ${{STATUS_DIR}}

# Determine task IDs based on execution mode
if [ -n $SLURM_ARRAY_TASK_ID ]; then
  TASK_IDS=($SLURM_ARRAY_TASK_ID)
  RUN_CMD="srun -n ${{NPROC}}"
else
  TASK_IDS=($(seq 0 ${{N_SLICES}}))
  RUN_CMD="mpirun -n ${{NPROC}}"
fi

for TASK_ID in ${{TASK_IDS[@]}}; do
  if [ $TASK_ID -lt $N_SLICES ]; then
    # Stage 2 merge task
    source $CCTBX_ACTIVATE
    SLICE_IDX=$(printf "%03d" $TASK_ID)
    PHIL_FILE="${{SCRIPTS_DIR}}/stage2_slice_${{SLICE_IDX}}.phil"
    echo "Running stage 2 merge for slice ${{SLICE_IDX}}..."
    ${{RUN_CMD}} cctbx.xfel.merge ${{PHIL_FILE}}

    # Write status file on completion
    touch ${{STATUS_DIR}}/slice_${{SLICE_IDX}}.done
  else
    # Phenix refinement coordinator task

    # In SLURM mode, poll for all merges to complete
    if [ -n $SLURM_ARRAY_TASK_ID ]; then
      echo "Waiting for all stage 2 merge tasks to complete..."
      while true; do
        DONE_COUNT=$(ls ${{STATUS_DIR}}/slice_*.done 2>/dev/null | wc -l)
        if [ $DONE_COUNT -ge $N_SLICES ]; then
          echo "All $N_SLICES merge tasks complete."
          break
        fi
        echo "Waiting... $DONE_COUNT / $N_SLICES complete"
        sleep 5
      done
    fi

    # Activate phenix environment and run all refinements in parallel
    source $PHENIX_ACTIVATE
    echo "Starting phenix refinements..."
    for i in $(seq 0 $((N_SLICES - 1))); do
      SLICE_IDX=$(printf "%03d" $i)
      SLICE_DIR="${{STAGE2_DIR}}/slice_${{SLICE_IDX}}"
      MTZ_FILE="${{SLICE_DIR}}/${{MTZ_NAME}}"

      if [ -f "${{MTZ_FILE}}" ]; then
        (
          cd "${{SLICE_DIR}}"
          phenix.refine ${{MTZ_FILE}} ${{PDB_FILE}} ${{PHENIX_PHIL}} output.prefix=refine_${{SLICE_IDX}} > refine_${{SLICE_IDX}}.log 2>&1
        ) &
      else
        echo "Warning: ${{MTZ_FILE}} not found, skipping slice ${{SLICE_IDX}}"
      fi
    done

    echo "Waiting for all phenix refinements to complete..."
    wait
    echo "All phenix refinements complete."

    # Scrape results from logs
    echo "Scraping results from logs..."
    RESULTS_FILE="${{STAGE2_DIR}}/spread_results.txt"
    echo "# slice wavelength f_prime[1..N] f_double_prime[1..N]" > $RESULTS_FILE

    for i in $(seq 0 $((N_SLICES - 1))); do
      SLICE_IDX=$(printf "%03d" $i)
      SLICE_DIR="${{STAGE2_DIR}}/slice_${{SLICE_IDX}}"
      MERGE_LOG="${{SLICE_DIR}}/iobs_main.log"
      REFINE_LOG="${{SLICE_DIR}}/refine_${{SLICE_IDX}}.log"

      # Extract wavelength from merge log
      WAVELENGTH=$(grep 'Average wavelength' $MERGE_LOG 2>/dev/null | tail -1 | awk '{{print $3}}' | tr -d '()' || echo "NA")

      # Extract f' values
      F_PRIME=$(grep 'f_prime' $REFINE_LOG 2>/dev/null | grep -v refine | tail -$N_SCATTERERS | awk '{{print $2}}' | tr '\\n' ' ' || echo "NA")

      # Extract f'' values
      F_DOUBLE_PRIME=$(grep 'f_double_prime' $REFINE_LOG 2>/dev/null | grep -v refine | tail -$N_SCATTERERS | awk '{{print $2}}' | tr '\\n' ' ' || echo "NA")

      echo "$SLICE_IDX $WAVELENGTH $F_PRIME $F_DOUBLE_PRIME" >> $RESULTS_FILE
    done

    echo "Results written to $RESULTS_FILE"
  fi
done
"""


class prepare_spread(worker):
  """
  Prepare data for SpReAD (Spectral Resolved Anomalous Diffraction) analysis.

  This worker:
  1. Computes global energy percentiles across all ranks
  2. Redistributes experiments/reflections by energy bin
  3. Writes binned data to disk for subsequent per-energy merging jobs
  """

  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(prepare_spread, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

  def __repr__(self):
    return "Prepare SpReAD energy-binned datasets"

  def run(self, experiments, reflections):
    self.logger.log_step_time("SPREAD_PREPARE")

    # Extract parameters
    n_bins = self.params.prepare.spread.n_energy_bins
    output_dir = self.params.prepare.spread.output_dir or self.params.output.output_dir

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

    # Step 6: Write batch scripts (only on rank 0)
    if self.mpi_helper.rank == 0:
      self.logger.log_step_time("WRITE_SCRIPTS")
      self._write_batch_scripts(n_bins, bin_edges, output_dir)
      self.logger.log_step_time("WRITE_SCRIPTS", True)

    self.logger.log_step_time("SPREAD_PREPARE", True)

    # Return original experiments and reflections for potential downstream processing
    return experiments, reflections

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

  def _write_batch_scripts(self, n_bins, bin_edges, output_dir):
    """
    Generate batch scripts for running stage 2 merging jobs.

    Writes:
      - stage2_slice_XXX.phil for each energy slice
      - run_spread.sh for SLURM/local execution
      - slices.txt metadata file
    """
    # Extract parameters
    window_width = self.params.prepare.spread.window_width
    window_step = self.params.prepare.spread.window_step
    stage2_phil_path = self.params.prepare.spread.stage2_phil
    stage2_nproc = self.params.prepare.spread.stage2_nproc
    stage2_output_dir = self.params.prepare.spread.stage2_output_dir or os.path.join(output_dir, 'stage2')

    # Phenix refinement parameters
    phenix_phil_path = self.params.prepare.spread.phenix_phil or ""
    phenix_pdb_path = self.params.prepare.spread.phenix_pdb or ""
    n_anomalous_scatterers = self.params.prepare.spread.n_anomalous_scatterers
    mtz_name = self.params.prepare.spread.mtz_name

    # SLURM parameters (optional)
    slurm_partition = self.params.prepare.spread.slurm_partition
    slurm_account = self.params.prepare.spread.slurm_account
    slurm_time_limit = self.params.prepare.spread.slurm_time_limit
    slurm_constraint = self.params.prepare.spread.slurm_constraint
    slurm_qos = self.params.prepare.spread.slurm_qos
    slurm_array_concurrency = self.params.prepare.spread.slurm_array_concurrency

    # Activation scripts
    cctbx_activate = self.params.prepare.spread.cctbx_activate
    phenix_activate = self.params.prepare.spread.phenix_activate

    # Read base phil content if provided
    base_phil_content = ""
    if stage2_phil_path and os.path.exists(stage2_phil_path):
      with open(stage2_phil_path, 'r') as f:
        base_phil_content = f.read()

    # Calculate valid slice centers
    # Window covers [center - width/2, center + width/2)
    # For 20% window: first valid center is 10 (covers 0-19), last is 90 (covers 80-99)
    half_width = window_width // 2
    first_center = half_width
    last_center = (n_bins - 1) - half_width + (1 if window_width % 2 == 0 else 0)

    slice_centers = list(range(first_center, last_center + 1, window_step))
    n_slices = len(slice_centers)

    self.logger.log(f"Generating {n_slices} slice phil files (window_width={window_width}%, step={window_step}%)")

    # Create scripts directory
    scripts_dir = os.path.join(output_dir, 'scripts')
    status_dir = os.path.join(scripts_dir, 'status')
    os.makedirs(scripts_dir, exist_ok=True)
    os.makedirs(stage2_output_dir, exist_ok=True)

    # Write slices.txt metadata
    slices_metadata_path = os.path.join(scripts_dir, 'slices.txt')
    with open(slices_metadata_path, 'w') as f:
      f.write("# slice_index center_percentile pct_start pct_end energy_start_eV energy_end_eV\n")
      for slice_idx, center in enumerate(slice_centers):
        pct_start = center - half_width
        pct_end = pct_start + window_width - 1
        energy_start = bin_edges[pct_start]
        energy_end = bin_edges[pct_end + 1]
        f.write(f"{slice_idx:03d} {center:3d} {pct_start:3d} {pct_end:3d} {energy_start:.2f} {energy_end:.2f}\n")

    # Determine bin index width for filename formatting
    bin_width = len(str(n_bins - 1))

    # Write individual slice phil files
    for slice_idx, center in enumerate(slice_centers):
      pct_start = center - half_width
      pct_end = pct_start + window_width - 1

      # Build input paths
      input_paths = "\n".join(
        f"input.path={output_dir}/pct_{pct:0{bin_width}d}_c*.*"
        for pct in range(pct_start, pct_end + 1)
      )

      slice_phil_content = SLICE_PHIL_TEMPLATE.format(
        stage2_phil_path=stage2_phil_path or "N/A",
        base_phil_content=base_phil_content,
        slice_idx=slice_idx,
        pct_start=pct_start,
        pct_end=pct_end,
        center=center,
        energy_start=bin_edges[pct_start],
        energy_end=bin_edges[pct_end + 1],
        input_paths=input_paths,
        slice_output_dir=os.path.join(stage2_output_dir, f'slice_{slice_idx:03d}')
      )

      slice_phil_path = os.path.join(scripts_dir, f'stage2_slice_{slice_idx:03d}.phil')
      with open(slice_phil_path, 'w') as f:
        f.write(slice_phil_content)

    # Build optional SLURM directives
    slurm_optional = ""
    if slurm_constraint:
      slurm_optional += f"#SBATCH --constraint={slurm_constraint}\n"
    if slurm_qos:
      slurm_optional += f"#SBATCH --qos={slurm_qos}\n"
    if slurm_partition:
      slurm_optional += f"#SBATCH --partition={slurm_partition}\n"
    if slurm_account:
      slurm_optional += f"#SBATCH --account={slurm_account}\n"

    # Write unified batch script
    spread_script_content = SPREAD_SCRIPT_TEMPLATE.format(
      n_slices=n_slices,
      slurm_array_concurrency=slurm_array_concurrency,
      slurm_time_limit=slurm_time_limit,
      slurm_optional_directives=slurm_optional,
      scripts_dir=scripts_dir,
      status_dir=status_dir,
      stage2_output_dir=stage2_output_dir,
      mtz_name=mtz_name,
      stage2_nproc=stage2_nproc,
      phenix_pdb_path=phenix_pdb_path,
      phenix_phil_path=phenix_phil_path,
      cctbx_activate=cctbx_activate,
      phenix_activate=phenix_activate,
      n_anomalous_scatterers=n_anomalous_scatterers
    )

    spread_script_path = os.path.join(scripts_dir, 'run_spread.sh')
    with open(spread_script_path, 'w') as f:
      f.write(spread_script_content)
    os.chmod(spread_script_path, 0o755)

    self.logger.log(f"Wrote {n_slices} slice phil files to {scripts_dir}")
    self.logger.log(f"Wrote batch script: run_spread.sh (use 'sbatch' for SLURM, 'bash' for local)")