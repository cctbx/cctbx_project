from __future__ import division

import iotbx.phil

master_phil_scope = iotbx.phil.parse("""
binning = 1
  .type = int(value_min=1)
brightness = 100
  .type = float(value_min=0.0)
colour_scheme = *greyscale rainbow heatmap inverse_greyscale
  .type = choice
format = jpeg *png tiff
  .type = choice
output_dir = None
  .type = path
""")

master_params = master_phil_scope.fetch().extract()


colour_schemes = {
  'greyscale': 0,
  'rainbow': 1,
  'heatmap': 2,
  'inverse_greyscale': 3
}

def run(args):
  import os
  from libtbx.phil import command_line
  from libtbx.utils import Sorry, Usage

  if len(args) == 0:
    from cStringIO import StringIO
    s = StringIO()
    master_phil_scope.show(out=s)
    raise Usage("""\
dxtbx.export_bitmaps image_files [options]

% s
""" %s.getvalue())

  from dxtbx.datablock import DataBlockFactory
  unhandled = []
  datablocks = DataBlockFactory.from_args(
    args, verbose=False, unhandled=unhandled)
  assert len(datablocks) > 0
  imagesets = datablocks[0].extract_imagesets()

  cmd_line = command_line.argument_interpreter(master_params=master_phil_scope)
  working_phil = cmd_line.process_and_fetch(args=unhandled)
  working_phil.show()
  params = working_phil.extract()

  brightness = params.brightness / 100
  vendortype = "made up"

  # check that binning is a power of 2
  binning = params.binning
  if not (binning > 0 and ((binning & (binning - 1)) == 0)):
    raise Sorry("binning must be a power of 2")

  output_dir = params.output_dir
  if output_dir is None:
    output_dir = "."
  elif not os.path.exists(output_dir):
    os.makedirs(output_dir)

  from rstbx.slip_viewer.tile_generation \
       import _get_flex_image, _get_flex_image_multipanel

  for imageset in imagesets:
    detector = imageset.get_detector()
    panel = detector[0]
    # XXX is this inclusive or exclusive?
    saturation = panel.get_trusted_range()[1]
    for i_image, image in enumerate(imageset):

      if len(detector) > 1:
        # FIXME This doesn't work properly, as flex_image.size2() is incorrect
        # also binning doesn't work
        assert binning == 1
        flex_image = _get_flex_image_multipanel(
          brightness=brightness,
          panels=detector,
          raw_data=image)
      else:
        flex_image = _get_flex_image(
          brightness=brightness,
          data=image,
          binning=binning,
          saturation=saturation,
          vendortype=vendortype)

      flex_image.setWindow(0, 0, 1)
      flex_image.adjust(color_scheme=colour_schemes.get(params.colour_scheme))

      # now export as a bitmap
      flex_image.prep_string()
      import Image
      # XXX is size//binning safe here?
      pil_img = Image.fromstring(
        'RGB', (flex_image.size2()//binning,
                flex_image.size1()//binning),
        flex_image.export_string)

      basename = os.path.basename(os.path.splitext(imageset.paths()[i_image])[0])
      path = os.path.join(
        output_dir, basename + '.' + params.format)

      print "Exporting %s" %path
      tmp_stream = open(path, 'wb')
      pil_img.save(tmp_stream, format=params.format)
      tmp_stream.close()


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
