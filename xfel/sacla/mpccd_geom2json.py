from __future__ import absolute_import, division, print_function
import sys, os
from dxtbx.model.detector import Detector
from scitbx import matrix
from libtbx.phil import parse
from libtbx.utils import Sorry, Usage

help_str = """Converts a SACLA geometry file to DIALS json format."""

phil_scope = parse("""
  geom_path = None
    .type = str
    .help = SACLA geometry file to convert
  distance = None
    .type = float
    .help = Detector distance
  pixel_size = 0.05
    .type = float
    .help = Pixel size in mm
  wavelength = 1.0
    .type = float
    .help = Wavelength for the simple beam model that will be included
""")

def run(args):
  if '-h' in args or '--help' in args or '-c' in args:
    print(help_str)
    phil_scope.show(attributes_level=2)
    return

  user_phil = []
  for arg in args:
    if os.path.isfile(arg):
      user_phil.append(parse("geom_path=%s"%arg))
    else:
      try:
        user_phil.append(parse(arg))
      except Exception as e:
        raise Sorry("Unrecognized argument: %s"%arg)
  params = phil_scope.fetch(sources=user_phil).extract()
  if params.distance is None:
    raise Usage("Please specify detector distance")

  geom = {}
  for line in open(params.geom_path):
    if len(line.split("=")) != 2: continue
    if "rigid_group" in line and not "collection" in line:
      geom[line.split("=")[1].strip()] = {}
    else:
      for key in geom:
        if line.startswith("%s/"%key):
          geom[key][line.split("=")[0].split("/")[1].strip()] = line.split("=")[-1].strip()

  detector = Detector()
  root = detector.hierarchy()
  root.set_frame(
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, - params.distance))

  for i, key in enumerate(sorted(geom)):
    fs_x, fs_y = geom[key]['fs'].split(" ")
    ss_x, ss_y = geom[key]['ss'].split(" ")
    fast = matrix.col((-float(fs_x.rstrip('x')),float(fs_y.rstrip('y')), 0.0))
    slow = matrix.col((-float(ss_x.rstrip('x')),float(ss_y.rstrip('y')), 0.0))

    origin = matrix.col((-float(geom[key]['corner_x']) * params.pixel_size,
                          float(geom[key]['corner_y']) * params.pixel_size,
                          0.0))

    # OBS! you need to set the panel to a root before set local frame...
    p = root.add_panel()
    p.set_name('panel-%s' % key)
    p.set_image_size((512, 1024))
    p.set_trusted_range((-1, 1000000))
    p.set_pixel_size((params.pixel_size, params.pixel_size))
    p.set_local_frame(
      fast.elems,
      slow.elems,
      origin.elems)

  from dxtbx.model import BeamFactory
  wavelength = params.wavelength
  beam = BeamFactory.simple(wavelength)

  from dxtbx.model import Experiment, ExperimentList
  from dxtbx.model.experiment_list import ExperimentListDumper
  experiments = ExperimentList()
  experiment = Experiment(detector = detector, beam = beam)
  experiments.append(experiment)
  dump = ExperimentListDumper(experiments)
  dump.as_json("geometry.json")

if __name__ == "__main__":
  run(sys.argv[1:])
