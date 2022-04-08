from __future__ import absolute_import, division, print_function

# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.h5_average

import numpy as np
import h5py
import iotbx.phil
import libtbx.load_env
from dials.util.options import ArgumentParser
import sys, os
from six.moves import range

phil_scope = iotbx.phil.parse("""
  average = True
    .type = bool
    .help = Generate average image
  max = True
    .type = bool
    .help = Generate maximum projection
  stddev = True
    .type = bool
    .help = Generate standard deviation image
""")

class Processh5(object):
  """
  Compute the average ("avg"), maximum projection ("max") or standard deviation ("stddev") of a
  series of images provided as a multi-image hdf5 file.
  """
  def __init__(self, h5file):
    # read everything from a multi-image h5 file
    self.h5file = h5file
    self.readfile = h5py.File(h5file, "r")
    self.shape_xy = self.readfile.values()[1].values()[0].value.shape
    self.length = len(self.readfile.values())
    #from IPython import embed; embed()

  def prepare_writefile_and_array(self, func):
    writefile_parts = (os.path.basename(os.path.splitext(self.h5file)[0]),
                       os.path.splitext(self.h5file)[1])
    writefile_name = ("_"+func).join(writefile_parts)
    writefile = h5py.File(writefile_name, "w")
    arr = np.zeros((self.length, self.shape_xy[0], self.shape_xy[1]))
    self.readfile.copy(str(self.readfile.values()[1].name), writefile['/'])
    return (writefile, writefile_name, arr)

  def process(self, func):
    # process and write everything to a new avg (or other) file
    writefile, writefile_name, arr = self.prepare_writefile_and_array(func)
    for ii in range(self.length-1):
      arr[ii][:][:] = np.asarray(self.readfile.values()[ii+1].values()[0].value)
    func_lookup = {
      "avg":np.mean,
      "max":np.max,
      "stddev":np.std
    }
    f = func_lookup[func]
    res = f(arr, axis=0)

    # write results to new file
    val = writefile.values()[0].values()[0]
    val.write_direct(res)
    writefile.close()
    print("Wrote", writefile_name)

  def cleanup(self):
    self.readfile.close()

def run(args):
  if ("--help" in args) or ("-h" in args) or (len(args) == 0):
    print("Usage: %s r25792.h5" % libtbx.env.dispatcher_name)
    return
  elif ("--config" in args) or ("-c" in args):
    iotbx.phil.parse(phil_scope).show(attributes_level=2)
    return
  h5s = []
  for arg in args:
    if arg.endswith(".h5") or arg.endswith(".hdf5"):
      h5s.append(args.pop(args.index(arg)))
  sys.argv = [sys.argv[0]] + args
  parser = ArgumentParser(phil=phil_scope)
  params, options = parser.parse_args(show_diff_phil=True)
  for h5 in h5s:
    print("Processing image %s..." % h5)
    processor = Processh5(h5)
    if params.average:
      processor.process("avg")
    if params.max:
      processor.process("max")
    if params.stddev:
      processor.process("stddev")
    processor.cleanup()

if __name__ == "__main__":
  run(sys.argv[1:])

