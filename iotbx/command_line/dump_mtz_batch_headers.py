from __future__ import division

def print_batch(b):
  print 80 * '-'
  print 'Batch: %d Wavelength: %.5f+/-%.5f' % (b.num(), b.alambd(), b.delamb())
  print 'Cell: %.3f %.3f %.3f %.3f %.3f %.3f' % tuple(b.cell())
  print 'Scale factors K/B: %.3f %.3f' % (b.bscale(), b.bbfac())
  # FIXME check for anisotropic mosaic spread tensor
  print 'Mosaic: %.3f' % b.crydat()[0]
  print 'Detector limits: %d %d %d %d' % tuple(b.detlm()[:4])
  print 'Detector distance: %.3f' % b.dx()[0]
  print 'Phi start, end, range: %.3f %.3f %.3f' % \
    (b.phistt(), b.phiend(), b.phirange())
  print 'Phi XYZ:\n %s' % ((2 * '\t\t%9.5f %9.5f %9.5f\n') % tuple(b.phixyz()))
  print 'U:\n %s' % ((3 * '\t\t%9.5f %9.5f %9.5f\n') % tuple(b.umat()))

def print_batch_old(b):

  ignore = ['show', 'mtz_object']
  for attr in sorted(dir(b)):
    if attr.startswith('_'):
      continue
    if attr.startswith('set'):
      continue
    if attr in ignore:
      continue
    value = getattr(b, attr)()
    print '---- %s:%s ----' % (attr, type(value))
    if 'scitbx_array_family' in str(type(value)):
      print tuple(value)
    else:
      print value

  return

def dump_mtz_batch_headers(mtz_file, how_many=0):
  from iotbx import mtz
  m = mtz.object(mtz_file)
  if how_many == 0:
    how_many = len(m.batches())
  for j, b in enumerate(m.batches()):
    if j >= how_many:
      break
    print_batch(b)

if __name__ == '__main__':
  import sys
  if len(sys.argv) > 2:
    dump_mtz_batch_headers(sys.argv[1], int(sys.argv[2]))
  else:
    dump_mtz_batch_headers(sys.argv[1])
