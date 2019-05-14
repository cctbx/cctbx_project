from __future__ import absolute_import, division, print_function
def tst_coordinate_frame_converter():
  import libtbx.load_env
  rstbx = libtbx.env.dist_path('rstbx')
  from rstbx.cftbx.coordinate_frame_converter import \
     coordinate_frame_converter
  import os
  cfc = coordinate_frame_converter(os.path.join(rstbx, 'cftbx', 'tests',
                                                'example-xparm.xds'))
  from scitbx import matrix
  x = matrix.col((1, 0, 0))
  assert(cfc.move_c(x, convention = cfc.CBF).angle(
      cfc.get_c('detector_fast')) < 0.1)

  print('OK')

if __name__ == '__main__':
  tst_coordinate_frame_converter()
