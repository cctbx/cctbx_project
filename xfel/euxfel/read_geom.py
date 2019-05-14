from __future__ import absolute_import, division, print_function
import sys, os
from scitbx.matrix import col
from libtbx.phil import parse
from libtbx.utils import Sorry

help_str = """Converts a CrystFEL file to DIALS json format."""

phil_scope = parse("""
  geom_file = None
    .type = str
    .help = CrystFEL geometry file to convert
  show_plot = False
    .type = bool
    .help = plot of detector geometry
""")

class panel_group(dict):
  def __init__(self):
    self.center = None
    self.local_origin = None
    self.local_fast = col((1,0,0))
    self.local_slow = col((0,1,0))

# Function to read in the crystFEL .geom file.
def read_geom(geom_file):
  panels = {}
  rigid_groups = {}
  collections = {}

  def known_panels():
    all_keys = []
    for value in rigid_groups.values():
      all_keys.extend(value)
    for value in collections.values():
      all_keys.extend(value)
    return set(all_keys)

  pixel_size = None

  for line in open(geom_file):
    line = line.split(';')[0]
    if len(line.split("=")) != 2: continue
    key, value = [w.strip() for w in line.split("=")]
    if key == 'res':
      pixel_size = 1000/float(value) # mm
    elif "rigid_group" in key:
      if "collection" in key:
        collections[key.split('rigid_group_collection_')[1]] = value.split(',')
      else:
        rigid_groups[key.split('rigid_group_')[1]] = value.split(',')
    else:
      if '/' not in key: continue
      panel = key.split("/")[0].strip()
      key = key.split("/")[1].strip()
      value = line.split("=")[-1].strip()
      if panel not in known_panels(): continue
      if panel not in panels:
        panels[panel] = {}
      panels[panel][key] = value

  mapping = {}
  for panel in panels:
    mapping[panel] = {}
    for group in rigid_groups:
      if panel not in rigid_groups[group]: continue
      for collection in collections:
        if group in collections[collection]:
          mapping[panel][collection] = group, len(rigid_groups[group])
  # example of mapping entry: mapping['p0a0'] = {'asics': ('p0', 8), 'quadrants': ('q0', 32)}
  parents = {}
  for panel in panels:
    parents[panel] = [mapping[panel][k][0] for k in sorted(mapping[panel], key = lambda x: x[1], reverse = True)]
  # example of parents entry:  parents['p0a0'] = ['q0', 'p0']
  # IE parents are listed in reverse order of immediacy (p0 is the parent of p0a0 and q0 is the parent of p0)

  hierarchy = panel_group()
  def add_node(panel, parent, parents, depth):
    if depth == len(parents):
      parent[panel] = panels[panel]
    else:
      if parents[depth] not in parent:
        parent[parents[depth]] = panel_group()
      add_node(panel, parent[parents[depth]], parents, depth+1)

  for panel in panels:
    add_node(panel, hierarchy, parents[panel], 0)
  # example of a hierarchy entry:
  # hierarchy['q0']['p0']['p0a0'] = full panel dictionary

  def parse_vector(vector):
    try:
      x, y, z = vector.split(" ")
    except ValueError:
      x, y = vector.split(" ")
      return col((float(x.rstrip('x')),float(y.rstrip('y')), 0.0))
    else:
      return col((float(x.rstrip('x')),float(y.rstrip('y')), float(z.rstrip('z'))))

  # set up panel vectors in lab space
  for panel in panels:
    panels[panel]['origin'] = col((float(panels[panel]['corner_x'])*pixel_size,
                                   float(panels[panel]['corner_y'])*pixel_size, 0.0))
    if 'coffset' in panels[panel]:
      panels[panel]['origin'] += col((0,0,1000*float(panels[panel]['coffset'])))
    panels[panel]['fast'] = panels[panel]['local_fast'] = parse_vector(panels[panel]['fs']).normalize()
    panels[panel]['slow'] = panels[panel]['local_slow'] = parse_vector(panels[panel]['ss']).normalize()
    center_fast = panels[panel]['fast'] * pixel_size * (int(panels[panel]['max_fs'])-int(panels[panel]['min_fs'])+1)/2.0
    center_slow = panels[panel]['slow'] * pixel_size * (int(panels[panel]['max_ss'])-int(panels[panel]['min_ss'])+1)/2.0
    panels[panel]['center'] = panels[panel]['origin'] + center_fast + center_slow
    panels[panel]['pixel_size'] = pixel_size

  assert 'pixel_size' not in panels
  def setup_centers(node):
    if not isinstance(node, panel_group):
      return
    if node.center is not None:
      return

    center = col((0.0, 0.0, 0.0))
    for key, child in node.iteritems():
      if isinstance(child, panel_group):
        if child.center is None:
          setup_centers(child)
        center += child.center
      else:
        center += child['center']
    center /= len(node)
    node.center = center
  setup_centers(hierarchy)

  def setup_local_frames(node):
    if not isinstance(node, panel_group):
      return

    for key, child in node.iteritems():
      if isinstance(child, panel_group):
        child.local_origin = child.center - node.center
        setup_local_frames(child)
      else:
        child['local_origin'] = child['origin'] - node.center
  hierarchy.local_origin = hierarchy.center
  setup_local_frames(hierarchy)
  return hierarchy

def run(args):
  if '-h' in args or '--help' in args or '-c' in args:
    print(help_str)
    phil_scope.show(attributes_level=2)
    return

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

  hierarchy = read_geom(params.geom_file)

  # Plot the detector model highlighting the hierarchical structure of the detector
  def plot_node(cummulative, node, name):
    if isinstance(node, panel_group):
      plt.arrow(cummulative[0],cummulative[1],node.local_origin[0],node.local_origin[1])
      for childname, child in node.iteritems():
        plot_node(cummulative+node.local_origin, child, childname)
    else:
      plt.arrow(cummulative[0],cummulative[1],node['local_origin'][0],node['local_origin'][1])

      ori = node['origin']
      fast_at_zero = node['fast'] * node['pixel_size'] * (int(node['max_fs']) - int(node['min_fs']) + 1)
      slow_at_zero = node['slow'] * node['pixel_size'] * (int(node['max_ss']) - int(node['min_ss']) + 1)
      plt.arrow(ori[0], ori[1], fast_at_zero[0], fast_at_zero[1], color='blue')
      plt.arrow(ori[0], ori[1], slow_at_zero[0], slow_at_zero[1], color='red')

      plt.text(ori[0], ori[1], name)

  if params.show_plot:
    from matplotlib import pyplot as plt
    plot_node(col((0,0,0)), hierarchy, 'root')
    plt.xlim(-200,200)
    plt.ylim(200,-200)

    plt.show()

if __name__ == "__main__":
  run(sys.argv[1:])
