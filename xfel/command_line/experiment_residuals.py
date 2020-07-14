from __future__ import absolute_import, division, print_function
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

parser = ArgumentParser("Visualize prediction offsets for a single shot experiment",
                        formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("expt", help="path to a refined experiment file", type=str) #, required=True)
parser.add_argument("refl", help="path to an indexed reflection file (strong spots with assigned miller indices)",
                    type=str) #, required=True)
parser.add_argument("--line-scale", dest="lscale", help="scale the offset vector by this amount", default=25, type=float)
parser.add_argument("--line-color", dest="lcolor", help="display the offset vector with this color", default="#777777", type=str)
parser.add_argument("--scatter-color", dest="scatt_cmap", help="display the scatter points with this pylab colormap", default="bwr", type=str)
parser.add_argument("--clim", nargs=2, default=None, help="colormap limits e.g. --clim -0.01 0.01", type=float)
parser.add_argument("--face-color", dest="axcol", default='w', help="pylab axis face color", type=str)
parser.add_argument("--mark-scale", dest="mark_scale", default=15, help="scale of plot marker", type=int)
parser.add_argument("--edge-color", dest="edgecolor", help="color of marker edge", default="#777777", type=str)
parser.add_argument("--edge-width", dest="edgewidth", help="width of marker edge", default=0.5, type=float)
parser.add_argument("--arrow-length", dest="headlen", help="length of pylab arrow head", default=0.5, type=float)
parser.add_argument("--arrow-width", dest="headwid", help="width of pylab arrow head", default=0.5, type=float)
parser.add_argument("--no-arrow", dest="noarrow", action="store_true", help="do not add arrows to plot")
parser.add_argument("--three-d", dest="threed", action="store_true", help="display plot in 3D (for 3D detector models)")
args = parser.parse_args()

from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentListFactory
import numpy as np
import pylab as plt

fig = plt.figure()
projection = None
if args.threed:
    from mpl_toolkits.mplot3d import Axes3D
    projection="3d"
ax = plt.gca(projection=projection)

R = flex.reflection_table.from_file(args.refl)
El = ExperimentListFactory.from_json_file(args.expt, check_format=False)
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
            edgecolors=args.edgecolor, linewidths=args.edgewidth)

cbar = plt.colorbar(scat)
cbar.ax.set_ylabel("$\Delta \psi$", rotation=270, labelpad=15)
ax.set_aspect("equal")
ax.set_facecolor(args.axcol)
plt.title("Arrow points to prediction")
plt.show()
