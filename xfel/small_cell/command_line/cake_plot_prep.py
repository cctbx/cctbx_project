from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.small_cell.cake_plot_prep
from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentList
from dials.array_family import flex
import sys, glob
import matplotlib.pyplot as plt


help_str = """
Make a cake plot from DIALS spotfinder spots

A cake plot is the azimuthal angle of a spot on an image vs. its resolution.
Powder rings will appear as vertical stripes, with defects in geometry
causing them to appear wavy. A cake plot is also insensitive to badly masked
regions of the detector compared to a 1d radial average as the aziumuthal
angle of a spot isn't averaged into the 1d trace.

This script creates cake.npy which is used by
cctbx.xfel.small_cell.cake_plot. Run this script first to generate it, then
run cctbx.xfel.small_cell.cake_plot

This script expects files named "*_strong.expt" and "_strong.refl". Supply
the former and the script will seek for the latter.

Usage (note, wild cards are permitted, but quotes are recommended):
cctbx.xfel.small_cell.cake_plot_prep "<path>/*_strong.expt>"

Multiprocessing support is availible using MPI. Example:
mpirun cctbx.xfel.small_cell.cake_plot_prep "<path>/*_strong.expt>"
"""

def run(args):
  if "-h" in args or "--help" in args:
    if rank == 0:
      print(help_str)
    return

  filenames = []
  for arg in sys.argv[1:]:
    filenames.extend(glob.glob(arg))
  if not filenames:
    sys.exit("No data found")

  x, y = flex.double(), flex.double()
  det = None
  for fn in filenames:
    print (fn)
    #try:
    refls = flex.reflection_table.from_file(fn.split('_strong.expt')[0] + "_strong.refl")
    #except OSError:
    #  continue
    expts = ExperimentList.from_file(fn, check_format=False)
    for expt_id, expt in enumerate(expts):
      subset = refls.select(expt_id == refls['id'])
      if len(subset) > 200: continue
      det = expt.detector
      for panel_id, panel in enumerate(det):
        r = subset.select(subset['panel'] == panel_id)
        x_, y_, _ = r['xyzobs.px.value'].parts()
        pix = panel.pixel_to_millimeter(flex.vec2_double(x_, y_))
        c = panel.get_lab_coord(pix)
        x.extend(c.parts()[0])
        y.extend(c.parts()[1])

  if det:
    z = flex.double(len(x), sum([p.get_origin()[2] for p in det])/len(det))
    coords = flex.vec3_double(x,y,z)
    two_theta = coords.angle((0,0,-1))
    d = expts[0].beam.get_wavelength() / 2 / flex.sin(two_theta/2)
    azi = flex.vec3_double(x, y, flex.double(len(x), 0)).angle((0,1,0), deg=True)
    azi.set_selected(x < 0, 180+(180-azi.select(x<0)))
  else:
    d = flex.double()
    azi = flex.double()

  import numpy as np
  fig, axes = plt.subplots(1, 1, figsize=(6, 3))
  axes.plot(
    1/d.as_numpy_array(), azi.as_numpy_array(),
    marker='.', linestyle='none',
    markersize=0.5, alpha=0.5
    )
  axes.set_ylabel('Azimuthal Angle')
  axes.set_xlabel(r'Resolution ($\mathrm{\AA}$)')
  x_ticks = np.array([10, 5, 2, 1])
  axes.set_xticks(1/x_ticks)
  axes.set_xticklabels(x_ticks)
  fig.tight_layout()
  plt.show()

if __name__ == "__main__":
  run(sys.argv[1:])

