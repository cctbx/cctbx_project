from __future__ import division

def load_cxi_phil(path):
  from labelit.phil_preferences import iotbx_defs, libtbx_defs
  from iotbx import phil

  cxi_phil_str = iotbx_defs + libtbx_defs
  stream = open(path)
  master_phil = phil.parse(input_string=cxi_phil_str, process_includes=True)
  horizons_phil = master_phil.fetch(sources=[phil.parse(stream.read())]).extract()
  stream.close()

  return horizons_phil
