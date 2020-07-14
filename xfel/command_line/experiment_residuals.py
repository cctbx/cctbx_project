from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.experiment_residuals


from libtbx.phil import parse
from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentListFactory
import numpy as np
import pylab as plt
from dials.util import show_mail_on_error

help_message = '''
Visualize prediction offsets for a single shot experiment

Example:

  cctbx.xfel.experiment_residuals refined.expt indexed.refl
'''

phil_scope = parse('''
lscale = 10
  .type = float 
  .help = scale the offset vector by this amount 
lcolor = #777777
  .type = str
  .help = display the offset vector with this color 
scatt_cmap = bwr
  .type = str
  .help = display the scatter points with this pylab colormap 
clim = None 
  .type = floats(size=2)
  .help = colormap limits e.g. clim=[-0.01, 0.01] 
axcol = w
  .type = str
  .help = pylab axis face color
mark_scale = 15
  .type = int
  .help = scale of scatter plot marker 
edge_color = #777777
  .type = str
  .help = color of marker edge
edge_width = 0.5
  .type = float
  .help = width of marker edge
headlen = 0.5 
  .type = float
  .help = length of pylab arrow head
headwid = 0.5
  .type = float
  .help = width of pylab arrow head
noarrow = False
  .type = bool
  .help = do not add arrows to plot
threed = False
  .type = bool
  .help = make a 3d plot for e.g. tilted detectors, experimental
''', process_includes=True)


class Script:

  def __init__(self):
    from dials.util.options import OptionParser
    #import libtbx.load_env

    # Create the option parser
    #usage = "usage: %s experiment1.expt experiment2.expt reflections1.refl reflections2.refl" % libtbx.env.dispatcher_name
    self.parser = OptionParser(
      usage="",
      sort_options=True,
      phil=phil_scope,
      read_experiments=True,
      read_reflections=True,
      check_format=False,
      epilog=help_message)

  def run(self):
    ''' Parse the options. '''
    from dials.util.options import flatten_experiments, flatten_reflections
    # Parse the command line arguments
    args, options = self.parser.parse_args(show_diff_phil=True)

    # do stuff
    #fig = plt.figure()
    projection = None
    if args.threed:
      from mpl_toolkits.mplot3d import Axes3D
      projection="3d"
    ax = plt.gca(projection=projection)

    El = flatten_experiments(args.input.experiments)
    R = flatten_reflections(args.input.reflections)[0]
    DET = El.detectors()[0]
    nref = len(R)

    xyz = np.zeros((nref, 3))
    for i_ref in range(nref):
      x, y, _ = R[i_ref]["xyzobs.mm.value"]
      xcal, ycal, _ = R[i_ref]["xyzcal.mm"]
      pid = R[i_ref]['panel']
      panel = DET[pid]
      xyz_lab = panel.get_lab_coord((x,y))
      xyz_cal_lab = panel.get_lab_coord((xcal, ycal))
      xyz[i_ref] = xyz_lab

      diff = np.array(xyz_lab) - np.array(xyz_cal_lab)
      diff_scale = diff*args.lscale
      if args.threed:
        ax.plot( *zip(xyz_lab, diff_scale+xyz_cal_lab), color=args.lcolor)
      else:
        x, y, _ = xyz_lab
        ax.arrow(x, y, diff_scale[0], diff_scale[1], head_width=args.headwid, head_length=args.headlen, color=args.lcolor,
                 length_includes_head=not args.noarrow)

    delpsi = R['delpsical.rad']
    xvals, yvals, zvals = xyz.T

    vmax = max(abs(delpsi))
    vmin = -vmax
    if args.clim is not None:
      vmin, vmax = args.clim

    scatt_arg = xvals, yvals
    if args.threed:
      scatt_arg = scatt_arg + (zvals,)
    scat = ax.scatter(*scatt_arg, s=args.mark_scale, c=delpsi, cmap=args.scatt_cmap, vmin=vmin, vmax=vmax, zorder=2,
                      edgecolors=args.edge_color, linewidths=args.edge_width)

    cbar = plt.colorbar(scat)
    cbar.ax.set_ylabel("$\Delta \psi$", rotation=270, labelpad=15)
    ax.set_aspect("equal")
    ax.set_facecolor(args.axcol)
    plt.title("Arrow points to prediction")
    plt.show()


if __name__ == '__main__':
  with show_mail_on_error():
    script = Script()
    script.run()
