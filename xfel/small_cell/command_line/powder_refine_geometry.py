# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.powder_refine_geometry
from __future__ import division
import logging

from iotbx.phil import parse
from dials.util import log
from dials.util import show_mail_on_error
from dials.util.options import ArgumentParser
from xfel.small_cell.geometry_refiner import PowderGeometryRefiner


logger = logging.getLogger("dials.command_line.powder_refine_geometry")

help_message = """
Refine detector geometry using powder diffraction d-spacings.

Examples of usage:

# Basic usage with defaults (LaB6)
$ cctbx.xfel.powder_refine_geometry spots.expt spots.refl

# Custom d-spacings (e.g., for Si standard)
$ cctbx.xfel.powder_refine_geometry spots.expt spots.refl \\
    reference_d_spacings=3.135,1.920,1.637,1.357

# Refine only XY shift (fix distance and tilt)
$ cctbx.xfel.powder_refine_geometry spots.expt spots.refl \\
    refine.distance=False refine.tilt=False

# Specify output file
$ cctbx.xfel.powder_refine_geometry spots.expt spots.refl \\
    output.experiments=calibrated.expt

This tool uses spotfinding output from a powder standard with known d-spacings
to refine detector geometry. It optimizes a 5-parameter detector model:
- distance: shift along detector normal
- shift1/shift2: XY shifts along fast/slow axes
- tau2/tau3: tilts around fast/slow axes

The refinement minimizes the sum of squared differences between observed
d-spacings and the nearest reference d-spacing.

Default reference d-spacings are for LaB6 (SRM 660):
  4.156 A (100), 2.939 A (110), 2.399 A (111), 2.078 A (200), 1.858 A (210)
"""

phil_scope = parse(
    """
reference_d_spacings = 4.156 2.939 2.399 2.078 1.858
  .type = floats
  .help = Reference d-spacings in Angstroms. Default: first 5 LaB6 peaks.

d_min = 1.5
  .type = float
  .help = Minimum d-spacing to include in refinement

d_max = 20
  .type = float
  .help = Maximum d-spacing to include in refinement

max_distance_inv_ang = 0.002
  .type = float
  .help = Maximum distance from reference d-spacing in inverse Angstroms. \
          Reflections further than this from any reference will be excluded.

refine {
  distance = True
    .type = bool
    .help = Refine detector distance along normal
  shift = True
    .type = bool
    .help = Refine XY shift (shift1 and shift2)
  tilt = True
    .type = bool
    .help = Refine detector tilts (tau2 and tau3)
}

output {
  experiments = refined.expt
    .type = path
    .help = Output filename for refined experiments
  log = powder_refine_geometry.log
    .type = path
    .help = Output log file
}
"""
)


class Script(object):
    def __init__(self):
        usage = "$ cctbx.xfel.powder_refine_geometry EXPERIMENTS REFLECTIONS [options]"
        self.parser = ArgumentParser(
            usage=usage,
            phil=phil_scope,
            epilog=help_message,
            check_format=False,
            read_reflections=True,
            read_experiments=True,
        )

    def run(self):
        params, options = self.parser.parse_args(show_diff_phil=True)

        # Validate input
        if len(params.input.experiments) != 1 or len(params.input.reflections) != 1:
            raise ValueError("Please provide exactly one experiments file and "
                             "one reflections file")

        experiments = params.input.experiments[0].data
        reflections = params.input.reflections[0].data

        print(f"\nLoaded {len(experiments)} experiments with "
              f"{len(reflections)} reflections")
        print(f"Reference d-spacings: {params.reference_d_spacings}")

        # Create and run refiner
        refiner = PowderGeometryRefiner(experiments, reflections, params)
        result = refiner.run()

        if result is not None:
            # Report final geometry
            refiner.report_geometry_changes()

            # Save refined experiments
            refined_experiments = refiner.get_refined_experiments()
            print(f"\nSaving refined experiments to {params.output.experiments}")
            refined_experiments.as_file(params.output.experiments)

            print("\nRefinement complete.")
        else:
            print("\nNo refinement performed.")


if __name__ == "__main__":
    with show_mail_on_error():
        script = Script()
        script.run()
