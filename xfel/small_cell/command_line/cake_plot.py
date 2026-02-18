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

from collections import defaultdict
from dials.array_family import flex
from dials.util.options import ArgumentParser
from dxtbx.model.experiment_list import ExperimentList
import glob
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np
import os
import sys


def extract_panel_data(experiments, reflections):
    """Collect resolution and azimuthal angle per detector panel.

    Returns a dict ``panel_id -> {'d': list, 'azi': list}``.
    """
    panel_data = defaultdict(lambda: {'d': [], 'azi': []})
    for expts, refls in zip(experiments, reflections):
        for expt_id, expt in enumerate(expts):
            subset = refls.select(refls['id'] == expt_id)
            if len(subset) == 0:
                continue
            det = expt.detector
            for panel_id, panel in enumerate(det):
                r = subset.select(subset['panel'] == panel_id)
                if len(r) == 0:
                    continue
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
    return panel_data


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
        phil=None,
        read_reflections=True,
        read_experiments=True,
        check_format=False,
        epilog=help_str,
    )
    params, options = parser.parse_args(args, show_diff_phil=True)
    
    if not params.input.experiments:
        print("No Experiments found in the input")
        parser.print_help()
        return
    if not params.input.reflections:
        print(
            "No reflection data found in the input. "
            "Reflection tables are needed if n_subset_method != random and n_subset is not None"
        )
        parser.print_help()
        return
    else:
        if len(params.input.reflections) != len(params.input.experiments):
            sys.exit(
                "The number of input reflections files does not match the "
                "number of input experiments"
            )

    experiments = [
        ExperimentList(o.data) for o in params.input.experiments
    ]
    reflections = [
        o.data for o in params.input.reflections
    ]

    panel_data = extract_panel_data(experiments, reflections)
    plot_panel_data(panel_data)

if __name__ == '__main__':
    run()
