from __future__ import absolute_import, division, print_function

# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.experiment_residuals

from libtbx.phil import parse
import numpy as np
import pylab as plt
import sys
import os
from dials.util import show_mail_on_error
import h5py

help_message = '''
Visualize prediction offsets for a single shot experiment

Example:

  cctbx.xfel.experiment_residuals refined.expt indexed.refl
'''

phil_scope = parse('''
lscale = 25
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
cbarshrink = 0.5
  .type = float
  .help = factor by which to shrink the displayed colorbar
exper_id = 0
  .type = int
  .help = experiment id (if experiment file is multi-shot)
title_tag = None
  .type = str
  .help = add this tag to the title of the plots
plot_rad_trans = False
  .type = bool
  .help = if True, will display a histogram of radial and transverse offsets
output {
  file = None
    .type = str
    .help = path to an optional hdf5 file storing useful output
  overwrite = False
    .type = bool
    .help = force an overwrite
  only = False
    .type = bool
    .help = if True, only write an output file, and don't display plots
}
''', process_includes=True)


class Script:

  def __init__(self):
    from dials.util.options import ArgumentParser

    self.parser = ArgumentParser(
      usage="",
      sort_options=True,
      phil=phil_scope,
      read_experiments=True,
      read_reflections=True,
      check_format=False,
      epilog=help_message)

  def run(self):
    from dials.util.options import flatten_experiments, flatten_reflections
    params, options = self.parser.parse_args(show_diff_phil=True)

    if len(params.input.experiments) > 1:
      print("Please only pass a single experiment file. Exiting...")
      sys.exit()
    if len(params.input.reflections) > 1:
      print("Please only pass a single reflection file. Exiting...")
      sys.exit()

    # do stuff
    ax = plt.gca()

    El = flatten_experiments(params.input.experiments)
    R = flatten_reflections(params.input.reflections)[0]

    nexper = len(El)
    nexper_in_refl = len(set(R["id"]))
    if not nexper == nexper_in_refl:
      print("There are %d experiments and %d possible reflection sets, experiment and reflection table out of sync"
            % (nexper, nexper_in_refl))
      sys.exit()
    if params.exper_id < 0:
      print("Exper Id must be greater than 0")
      sys.exit()
    if params.exper_id > nexper:
      print("exper_id must be less than maximum number of experiments (=%d)" % nexper)
      sys.exit()

    DET = El[params.exper_id].detector
    pixsize = tuple([x*1000 for x in DET[0].get_pixel_size()])
    R = R.select(R["id"] == params.exper_id)

    columns = 'delpsical.rad', 'xyzobs.mm.value', 'panel', 'xyzcal.mm'
    misses = 0
    for col in columns:
      if col not in R[0].keys():
        print("reflection file missing column %s" % col)
        misses += 1
    if misses > 0:
      print("Please reformat refl file to include above columns")
      sys.exit()

    nref = len(R)

    xyz = np.zeros((nref, 3))
    all_rad = []
    all_trans = []
    all_diff = []
    for i_ref in range(nref):
      x, y, _ = R[i_ref]["xyzobs.mm.value"]
      xcal, ycal, _ = R[i_ref]["xyzcal.mm"]
      pid = R[i_ref]['panel']
      panel = DET[pid]
      xyz_lab = panel.get_lab_coord((x,y))
      xyz_cal_lab = panel.get_lab_coord((xcal, ycal))
      xyz[i_ref] = xyz_lab

      diff = np.array(xyz_lab) - np.array(xyz_cal_lab)

      # rad is the unit vector pointing to the observation
      xy_lab = np.array((xyz_lab[0], xyz_lab[1]))
      rad = xy_lab / np.linalg.norm(xy_lab)
      trans = np.array([-rad[1], rad[0]])

      rad_component = np.dot(diff[:2], rad)
      trans_component = np.dot(diff[:2], trans)

      diff_scale = diff*params.lscale
      x, y, _ = xyz_lab
      ax.arrow(x, y, diff_scale[0], diff_scale[1], head_width=params.headwid, head_length=params.headlen, color=params.lcolor,
               length_includes_head=not params.noarrow)
      all_rad.append(rad_component*1000)
      all_trans.append(trans_component*1000)
      all_diff.append(diff*1000)

    if params.output.file is not None:
      if os.path.exists(params.output.file) and not params.output.overwrite:
        print("File %s exists, to overwrite use output.overwrite=True" % params.output.file)
        sys.exit()
      with h5py.File(params.output.file, 'w') as h5:
        h5.create_dataset("radial_offset", data=all_rad)
        h5.create_dataset("transverse_offset", data=all_trans)
        h5.create_dataset("vectors_from_obs_to_cal", data=all_diff)
        exper_filename = params.input.experiments[0].filename
        refl_filename = params.input.reflections[0].filename
        h5.create_dataset("exper_filename", data=np.array([exper_filename], dtype=np.string_))
        h5.create_dataset("refl_filename", data=np.array([refl_filename], dtype=np.string_))
      print("Output saved to file %s" % params.output.file)
      if params.output.only:
        print("Done.")
        sys.exit()

    delpsi = R['delpsical.rad']
    xvals, yvals, zvals = xyz.T

    vmax = max(abs(delpsi))
    vmin = -vmax
    if params.clim is not None:
      vmin, vmax = params.clim

    scatt_arg = xvals, yvals
    scat = ax.scatter(*scatt_arg, s=params.mark_scale, c=delpsi, cmap=params.scatt_cmap, vmin=vmin, vmax=vmax, zorder=2,
                      edgecolors=params.edge_color, linewidths=params.edge_width)

    cbar = plt.colorbar(scat, shrink=params.cbarshrink)

    cbar.ax.set_title(r"$\Delta \psi$")
    ax.set_aspect("equal")
    ax.set_facecolor(params.axcol)
    title = "prediction offsets (arrow points to prediction)"
    title_radtrans = "radial/transverse components"
    if params.title_tag is not None:
      title += ": %s" % params.title_tag
      title_radtrans += ": %s" % params.title_tag
    ax.set_title(title)

    if params.plot_rad_trans:
      plt.figure()
      plt.hist(all_rad, bins='auto', histtype='step')
      plt.hist(all_trans, bins='auto', histtype='step')
      plt.gca().set_title(title_radtrans)
      rad_mn, rad_sig = np.mean(all_rad), np.std(all_rad)
      trans_mn, trans_sig = np.mean(all_trans), np.std(all_trans)
      rad_str = r"radial: %.2f $\pm$ %.2f $\mu m$" % (rad_mn, rad_sig)
      trans_str = r"transverse: %.2f $\pm$ %.2f $\mu m$" % (trans_mn, trans_sig)
      plt.legend((rad_str, trans_str))
      plt.xlabel(r"microns (1 pixel = %.1f $\mu m$ x %.1f $\mu m$)" % pixsize)
      plt.ylabel("number of spots")

    plt.show()


if __name__ == '__main__':
  with show_mail_on_error():
    script = Script()
    script.run()
