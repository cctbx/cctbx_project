from __future__ import absolute_import, division, print_function
import h5py, sys, shutil, os
from dxtbx.model.experiment_list import ExperimentList
from libtbx.phil import parse
from scitbx import matrix

"""
Given a JF16M NeXus master and a refined geometry file, sync the new geometry to
the NeXus master.
"""

phil_scope = parse("""
nexus_master = None
  .type = path
new_geometry = None
  .type = path
output_file = None
  .type = path
""")

dxtbx_to_nexus = matrix.sqr((-1, 0, 0,
                              0, 1, 0,
                              0, 0,-1))

def get_settings(node, apply_dxtbx_to_nexus=False):
  f = matrix.col(node.get_local_fast_axis())
  s = matrix.col(node.get_local_slow_axis())
  o = matrix.col(node.get_local_origin())
  n = f.cross(s)
  m = matrix.sqr((f[0], s[0], n[0],
                  f[1], s[1], n[1],
                  f[2], s[2], n[2]))
  angle, axis = m.r3_rotation_matrix_as_unit_quaternion().unit_quaternion_as_axis_and_angle(deg=True)
  axis = dxtbx_to_nexus * axis
  o = dxtbx_to_nexus * o
  return o, angle, axis

def sync(src, dest):
  assert dest.attrs['transformation_type'] == 'rotation'
  assert dest.attrs['units'] == 'degrees'
  offset, angle, vector = get_settings(src)
  if angle == 0:
    dest.attrs['vector'] = 0,0,-1
    dest[0] = 0
  else:
    dest.attrs['vector'] = vector.elems
    dest[0] = angle
  dest.attrs['offset'] = offset

def run(args):
  user_phil = []
  for arg in args:
    if os.path.isfile(arg):
      user_phil.append(parse(file_name=arg))
    else:
      try:
        user_phil.append(parse(arg))
      except Exception as e:
        raise Sorry("Unrecognized argument: %s"%arg)
  params = phil_scope.fetch(sources=user_phil).extract()

  shutil.copyfile(params.nexus_master, params.output_file)
  h5 = h5py.File(params.output_file, 'r+')

  experiments = ExperimentList.from_file(params.new_geometry, check_format=False)
  if len(experiments.detectors()) > 1:
    print ("Found more than one detector, using the first one.")
  detector = experiments.detectors()[0]
  root = detector.hierarchy()

  transformations = h5['entry/instrument/ELE_D0/transformations']
  rail = transformations['AXIS_RAIL']
  offset, angle, vector = get_settings(root, apply_dxtbx_to_nexus=True)
  rail[0] = offset[2] # distance
  rail.attrs['vector'] = 0,0,1
  d0 = transformations['AXIS_D0']
  d0.attrs['offset'] = offset[0], offset[1], 0
  if angle == 0:
    d0[0] = 0
    d0.attrs['vector'] = 0,0,-1
  else:
    d0[0] = angle
    d0.attrs['vector'] = vector.elems

  for q in range(4):
    sync(root[q], transformations[root[q].get_name()])
    for m in range(8):
      sync(root[q][m], transformations[root[q][m].get_name()])

if __name__ == "__main__":
  run(sys.argv[1:])
