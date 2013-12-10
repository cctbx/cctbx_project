from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME cxi.display_metrology
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1
# $Id
#

import sys, os
import libtbx.phil
from libtbx.utils import Usage, Sorry
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from xfel.cxi.cspad_ana.cspad_tbx import xpp_active_areas
import sys


master_phil = libtbx.phil.parse("""
metrology = None
  .type = str
  .help = File with metrology information or XPP active area name
  .optional = False
""")

if (__name__ == "__main__") :
  user_phil = []
  for arg in sys.argv[1:]:
    if (os.path.isfile(arg)) or arg in xpp_active_areas:
      user_phil.append(libtbx.phil.parse("""metrology=\"%s\"""" % arg))
    else :
      try :
        user_phil.append(libtbx.phil.parse(arg))
      except RuntimeError, e :
        raise Sorry("Unrecognized argument '%s' (error: %s)" % (arg, str(e)))

  params = master_phil.fetch(sources=user_phil).extract()
  if (params.metrology is None) :
    master_phil.show()
    raise Usage("metrology must be defined (either metrology=XXX, or the name alone).")

  if os.path.isfile(params.metrology):
    print "Not impl"
  else:
    active_areas = xpp_active_areas[params.metrology]['active_areas']

    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')

    for asic_number, (y1, x1, y2, x2) in enumerate([(active_areas[(i*4)+0]+1,
                                                     active_areas[(i*4)+1]+1,
                                                     active_areas[(i*4)+2]-1,
                                                     active_areas[(i*4)+3]-1) for i in xrange(len(active_areas)//4)]):
      ax.add_patch(Rectangle((x1,y1), x2-x1, y2-y1, color="grey"))
      ax.annotate(asic_number, (x1+(x2-x1)/2,y1+(y2-y1)/2))

    ax.set_xlim((0, 2000))
    ax.set_ylim((0, 2000))
    ax.set_ylim(ax.get_ylim()[::-1])

    plt.show()
