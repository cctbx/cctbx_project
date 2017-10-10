from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME dxtbx.plot_detector_models
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

from builtins import zip
usage = """Plot dxtbx detector models. Provide multiple json files if desired
Example: dxtbx.plot_detector_models datablock1.json datablock2.json
"""

import numpy as np
import matplotlib.pyplot as plt
import sys, os
from scitbx.matrix import col
from libtbx.phil import parse
from libtbx.utils import Sorry

from matplotlib.patches import FancyArrowPatch

from dxtbx.datablock import DataBlockFactory
from dxtbx.model.experiment_list import ExperimentListFactory

phil_scope = parse("""
  show_origin_vectors = True
    .type = bool
    .help = If true, draw origin vectors as arrows
  orthographic = False
    .type = bool
    .help = If true, draw an orthographic projection (IE drop the Z-axis)
""")

# http://stackoverflow.com/questions/22867620/putting-arrowheads-on-vectors-in-matplotlibs-3d-plot
class Arrow3D(FancyArrowPatch):
  def __init__(self, xs, ys, zs, *args, **kwargs):
    FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
    self._verts3d = xs, ys, zs

  def draw(self, renderer):
    from mpl_toolkits.mplot3d import proj3d
    xs3d, ys3d, zs3d = self._verts3d
    xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
    self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
    FancyArrowPatch.draw(self, renderer)

def run(args):
  user_phil=[]
  files = []
  for arg in args:
    if os.path.isfile(arg):
      files.append(arg)
    else:
      try:
        user_phil.append(parse(arg))
      except Exception as e:
        raise Sorry("Unrecognized argument %s"%arg)
  params = phil_scope.fetch(sources = user_phil).extract()
  if not params.orthographic:
    from mpl_toolkits.mplot3d import proj3d

  def plot_group(g, color, orthographic = False):
    # recursively plot a detector group
    p = g.parent()
    if params.show_origin_vectors:
      if p is None:
        #parent origin
        pori = (0,0,0)
      else:
        #parent origin
        pori = p.get_origin()
      ori = g.get_origin()
      if not orthographic:
        a = Arrow3D([pori[0], ori[0]], [pori[1], ori[1]],
                    [pori[2], ori[2]], mutation_scale=20,
                    lw=1, arrowstyle="-|>", color='gray')
        ax.add_artist(a)
    if g.is_group():
      for c in g:
        # plot all the children
        plot_group(c, color, orthographic)
    else:
      # plot the panel boundaries
      size = g.get_image_size()
      p0 = col(g.get_pixel_lab_coord((0,0)))
      p1 = col(g.get_pixel_lab_coord((size[0]-1,0)))
      p2 = col(g.get_pixel_lab_coord((size[0]-1,size[1]-1)))
      p3 = col(g.get_pixel_lab_coord((0,size[1]-1)))
      v1 = p1-p0
      v2 = p3-p0
      vcen = ((v2/2) + (v1/2)) + p0
      z = list(zip(p0,p1,p2,p3,p0))

      if orthographic:
        ax.plot(z[0], z[1], color=color)

        # Annotate with panel numbers
        ax.text(vcen[0], vcen[1], '%d'%g.index())
      else:
        ax.plot(z[0], z[1], z[2], color=color)

        # Annotate with panel numbers
        ax.text(vcen[0], vcen[1], vcen[2], '%d'%g.index())

  fig = plt.figure()
  colormap = plt.cm.gist_ncar
  colors = [colormap(i) for i in np.linspace(0, 0.9, len(files))]
  for file_name, color, in zip(files, colors):

    # read the data and get the detector models
    try:
      datablocks = DataBlockFactory.from_json_file(file_name, check_format=False)
    except Exception as e:
      experiments = ExperimentListFactory.from_json_file(file_name, check_format=False)
      detectors = experiments.detectors()
    else:
      detectors = []
      for datablock in datablocks:
        detectors.extend(datablock.unique_detectors())
    for detector in detectors:
      # plot the hierarchy
      if params.orthographic:
        ax = fig.gca()
      else:
        ax = fig.gca(projection='3d')
      plot_group(detector.hierarchy(), color, orthographic = params.orthographic)

  plt.show()

if __name__ == "__main__":
  run(sys.argv[1:])
