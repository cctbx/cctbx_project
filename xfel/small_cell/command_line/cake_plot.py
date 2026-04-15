from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.small_cell.cake_plot

help_str = """
Make a cake plot from DIALS spotfinder spots

A cake plot is the azimuthal angle of a spot on an image vs. its resolution.
Powder rings will appear as vertical stripes, with defects in geometry
causing them to appear wavy. A cake plot is also insensitive to badly masked
regions of the detector compared to a 1d radial average as the aziumuthal
angle of a spot isn't averaged into the 1d trace.

Example: cctbx.xfel.small_cell.cake_plot "
"""

import iotbx.phil
from collections import defaultdict
from dials.array_family import flex
from dials.util.options import ArgumentParser
from dxtbx.model.experiment_list import ExperimentList
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
import sys

phil_str = """
mp {
  method = *serial mpi multiprocessing
    .type = choice
    .help = Parallelization method
  nproc = 1
    .type = int
    .help = Number of processes for multiprocessing method
}
"""

phil_scope = iotbx.phil.parse(phil_str, process_includes=True)

# Module-level globals used by multiprocessing workers (set via pool initializer)
_mp_expts = None
_mp_refls = None
_mp_id_groups = None


def _process_experiments(expts, refls, expt_ids, id_groups):
    """Process a subset of experiment IDs from a single (expts, refls) pair.

    Returns a plain dict {panel_id: {'d': np.array, 'azi': np.array}}.
    """
    panel_data = defaultdict(lambda: {'d': [], 'azi': []})
    for expt_id in expt_ids:
        if expt_id not in id_groups:
            continue
        expt = expts[expt_id]
        subset = refls.select(flex.size_t(id_groups[expt_id].astype(int).tolist()))
        det = expt.detector
        panels_present = set(subset['panel'])
        for panel_id in panels_present:
            panel = det[panel_id]
            r = subset.select(subset['panel'] == panel_id)
            x_, y_, _ = r['xyzobs.px.value'].parts()
            pix = panel.pixel_to_millimeter(flex.vec2_double(x_, y_))
            xyz = panel.get_lab_coord(pix)
            x, y, z = xyz.parts()
            coords = flex.vec3_double(x, y, z)
            two_theta = coords.angle((0, 0, -1))
            d_vals = expt.beam.get_wavelength() / (2 * flex.sin(two_theta / 2))
            azi_vals = flex.vec3_double(x, y, flex.double(len(x), 0)).angle((0, 1, 0), deg=True)
            azi_vals.set_selected(x < 0, 180 + (180 - azi_vals.select(x < 0)))
            panel_data[panel_id]['d'].extend(d_vals.as_numpy_array())
            panel_data[panel_id]['azi'].extend(azi_vals.as_numpy_array())
    # Convert lists to numpy arrays for pickle safety
    return {pid: {'d': np.array(v['d']), 'azi': np.array(v['azi'])}
            for pid, v in panel_data.items()}


def _mp_worker(expt_ids):
    """Multiprocessing worker function. Reads expts/refls from module globals."""
    return _process_experiments(_mp_expts, _mp_refls, expt_ids, _mp_id_groups)


def _mp_initializer(expts, refls, id_groups):
    """Pool initializer: store shared data in module globals for forked workers."""
    global _mp_expts, _mp_refls, _mp_id_groups
    _mp_expts = expts
    _mp_refls = refls
    _mp_id_groups = id_groups


def merge_panel_data(results):
    """Merge a list of partial panel_data dicts by concatenating arrays per panel."""
    merged = defaultdict(lambda: {'d': [], 'azi': []})
    for partial in results:
        for pid, v in partial.items():
            merged[pid]['d'].append(v['d'])
            merged[pid]['azi'].append(v['azi'])
    return {pid: {'d': np.concatenate(v['d']), 'azi': np.concatenate(v['azi'])}
            for pid, v in merged.items() if v['d']}


def _pregroup_reflections(refls):
    """Pre-group reflections by experiment id using numpy argsort.

    Returns a dict {expt_id: numpy index array}.
    """
    id_arr = refls['id'].as_numpy_array()
    order = np.argsort(id_arr, kind='mergesort')
    sorted_ids = id_arr[order]
    splits = np.nonzero(np.diff(sorted_ids))[0] + 1
    group_indices = np.split(order, splits)
    unique_ids = sorted_ids[np.concatenate([[0], splits])]
    return {int(uid): idx for uid, idx in zip(unique_ids, group_indices)}


def extract_panel_data(experiments, reflections, params=None):
    """Collect resolution and azimuthal angle per detector panel.

    Returns a dict ``panel_id -> {'d': list/array, 'azi': list/array}``.
    """
    if params is None or params.mp.method == 'serial':
        return _extract_panel_data_serial(experiments, reflections)
    elif params.mp.method == 'mpi':
        return _extract_panel_data_mpi(experiments, reflections)
    elif params.mp.method == 'multiprocessing':
        return _extract_panel_data_mp(experiments, reflections, params.mp.nproc)


def _extract_panel_data_serial(experiments, reflections):
    """Serial extraction (original behavior)."""
    all_results = []
    for expts, refls in zip(experiments, reflections):
        id_groups = _pregroup_reflections(refls)
        expt_ids = list(range(len(expts)))
        result = _process_experiments(expts, refls, expt_ids, id_groups)
        all_results.append(result)
    return merge_panel_data(all_results) if len(all_results) > 1 else all_results[0]


def _extract_panel_data_mpi(experiments, reflections):
    """MPI-parallel extraction. Each rank processes a stride of experiment IDs."""
    from libtbx.mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    all_results = []
    for expts, refls in zip(experiments, reflections):
        id_groups = _pregroup_reflections(refls)
        my_expt_ids = list(range(rank, len(expts), size))
        result = _process_experiments(expts, refls, my_expt_ids, id_groups)
        all_results.append(result)

    # Merge local results across file pairs
    local_merged = merge_panel_data(all_results) if len(all_results) > 1 else all_results[0]

    # Gather all partial results to rank 0
    gathered = comm.gather(local_merged, root=0)
    if rank == 0:
        return merge_panel_data(gathered)
    return None


def _extract_panel_data_mp(experiments, reflections, nproc):
    """Multiprocessing-parallel extraction using fork + global initializer."""
    import multiprocessing

    all_results = []
    for expts, refls in zip(experiments, reflections):
        id_groups = _pregroup_reflections(refls)
        n_expts = len(expts)
        # Divide experiment IDs into nproc chunks
        all_ids = list(range(n_expts))
        chunk_size = (n_expts + nproc - 1) // nproc
        chunks = [all_ids[i:i + chunk_size] for i in range(0, n_expts, chunk_size)]

        with multiprocessing.Pool(
            processes=nproc,
            initializer=_mp_initializer,
            initargs=(expts, refls, id_groups),
        ) as pool:
            results = pool.map(_mp_worker, chunks)
        all_results.extend(results)

    return merge_panel_data(all_results)


def plot_panel_data(panel_data):
    """Create the cake plot.

    ``panel_data`` is the dict returned by ``extract_panel_data``.
    The figure is displayed interactively.
    """
    cmap = plt.get_cmap('tab20')
    fig, ax = plt.subplots(figsize=(6, 3))
    for i, (panel_id, data) in enumerate(sorted(panel_data.items())):
        d_arr = np.array(data['d'])
        azi_arr = np.array(data['azi'])
        mask = d_arr > 0
        if not mask.any():
            continue
        x = 1.0 / d_arr[mask]  # plotted values: 1/d (1/Å)
        y = azi_arr[mask]
        ax.scatter(x, y, s=0.5, alpha=0.5, color=cmap(i % 20), label=f'Panel {panel_id}')
    ax.set_ylabel('Azimuthal Angle (deg)')
    ax.set_xlabel('Resolution (Å)')
    # Format x‑axis to show resolution instead of 1/d
    def resolution_formatter(x, pos):
        if x == 0:
            return '-'
        return f"{1/x:.2f}"
    ax.xaxis.set_major_formatter(FuncFormatter(resolution_formatter))
    fig.tight_layout()
    plt.show()


def run(args=None):
    usage = (
        "Usage: cctbx.xfel.small_cell.cake_plot "
        "experiments1.expt experiments2.expt reflections1.refl "
        "reflections2.refl..."
    )

    parser = ArgumentParser(
        usage=usage,
        phil=phil_scope,
        read_reflections=True,
        read_experiments=True,
        check_format=False,
        epilog=help_str,
    )
    params, options = parser.parse_args(args, show_diff_phil=False)

    # Determine rank for MPI guard
    is_rank0 = True
    if params.mp.method == 'mpi':
        from libtbx.mpi4py import MPI
        is_rank0 = MPI.COMM_WORLD.Get_rank() == 0

    if is_rank0:
        diff_phil_str = parser.diff_phil.as_str()
        if diff_phil_str:
            print("The following parameters have been modified:\n")
            print(diff_phil_str)

    if not params.input.experiments:
        if is_rank0:
            print("No Experiments found in the input")
            parser.print_help()
        return
    if not params.input.reflections:
        if is_rank0:
            print(
                "No reflection data found in the input. "
                "Reflection tables are needed if n_subset_method != random and n_subset is not None"
            )
            parser.print_help()
        return
    else:
        if len(params.input.reflections) != len(params.input.experiments):
            if is_rank0:
                sys.exit(
                    "The number of input reflections files does not match the "
                    "number of input experiments"
                )
            return

    experiments = [
        ExperimentList(o.data) for o in params.input.experiments
    ]
    reflections = [
        o.data for o in params.input.reflections
    ]

    panel_data = extract_panel_data(experiments, reflections, params)
    if is_rank0 and panel_data is not None:
        plot_panel_data(panel_data)

if __name__ == '__main__':
    run()
