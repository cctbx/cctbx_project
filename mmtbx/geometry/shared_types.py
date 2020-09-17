from __future__ import absolute_import, division, print_function

import boost_adaptbx.boost.python as bp
ext = bp.import_ext( "mmtbx_geometry_shared_types_ext" )
from mmtbx_geometry_shared_types_ext import *

def calculate_base_for_coordinates(xyzs):

  if len( xyzs ) == 0:
    return ( 0, 0, 0 )

  else:
    ( xs, ys, zs ) = xyzs.parts()
    return ( min( xs ), min( ys ), min( zs ) )
