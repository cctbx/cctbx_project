from __future__ import division
from __future__ import print_function
# LIBTBX_SET_DISPATCHER_NAME cspad.quadrants
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

from xfel.cxi.display_spots import ImageFactory,cxi_phil
from xfel.cxi import display_spots
import sys,os,copy
from scitbx.array_family import flex
from libtbx.utils import Sorry

def view_raw_image(path, *command_line, **kwargs):
  args = [path,
          "viewer.powder_arcs.show=False",
          "viewer.powder_arcs.code=3n9c",
         ]

  horizons_phil = cxi_phil.cxi_versioned_extract(
                    copy.deepcopy(args),list(command_line))

  #global parameters
  display_spots.parameters.horizons_phil = horizons_phil
  return horizons_phil


if (__name__ == "__main__"):
  files = [arg for arg in sys.argv[1:] if os.path.isfile(arg)]
  arguments = [arg for arg in sys.argv[1:] if not os.path.isfile(arg)]

  for file in files:
    message="""Based on the file %s,
    this program will compute incremental quadrant translations to circularize
    powder rings on the inner four sensors.  The algorithm treats each quadrant
    independently and scores based on self-correlation upon 45-degree rotation.
    Quad-translations and per-tile unit_translations that are already defined in
    spotfinder/applications/xfel/cxi_phil.py are applied first, and increments
    are determined ON TOP OF those values.  Output is given in the form of
    a new quad_translation array."""%file
    print(message)
    phil_params = view_raw_image(file, *arguments, **({'display':True}))
    image = ImageFactory(file)
    image.show_header()

    if image.horizons_phil_cache.distl.detector_format_version != \
       phil_params.distl.detector_format_version:
         raise Sorry(
         '''it is necessary to put distl.detector_format_version="%s" on the command line.'''%(
       image.horizons_phil_cache.distl.detector_format_version))

    M = image.get_tile_manager(phil = phil_params)
    tiling = M.effective_tiling_as_flex_int(encode_inactive_as_zeroes=True)

    # somebodys stupid private cache
    #print image.horizons_phil_cache.distl.detector_format_version
    #print image.horizons_phil_cache.distl.quad_translations
    #print image.horizons_phil_cache.distl.tile_translations

    dfv = image.horizons_phil_cache.distl.detector_format_version
    if dfv is None or dfv.find("CXI")>=0:
      key_sensors =[(2,3),(18,19),(50,51),(34,35)] # UL, UR, LL, LR
    elif dfv.find("XPP")>=0:
      key_sensors =[(34,35),(50,51),(18,19),(2,3)] # UL, UR, LL, LR

    if dfv is None:
      old_quad_trans = flex.int(8)
      new_quad_trans = flex.int(8) #initialize
    else:
      old_quad_trans = flex.int(image.horizons_phil_cache.distl.quad_translations)
      new_quad_trans = old_quad_trans + flex.int(len(old_quad_trans)) #initialize

    from xfel.metrology.quadrant import one_sensor
    for isensor,sensor in enumerate(key_sensors):
      Q = one_sensor(image,sensor,M)
      print(Q.coordmax[0], Q.coordmax[1])
      new_quad_trans[isensor*2]-=Q.coordmax[0]
      new_quad_trans[isensor*2 + 1]-=Q.coordmax[1]

    print("The OLD QUAD translations are:",list(old_quad_trans))
    print("\nThe NEW QUAD translations are:",list(new_quad_trans))
    print("""These should be pasted into spotfinder/applications/xfel/cxi_phil.py
    with NEW replacing the OLD.""")


#cspad.quadrants /reg/d/psdm/cxi/cxie1414/ftc/brewster/averages/r0006/000/out/max-Ds2-r0006.pickle
