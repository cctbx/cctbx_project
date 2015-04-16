from __future__ import division
#!/usr/bin/env python
# detector.py
#   Copyright (C) 2011 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# A model for the detector for the "updated experimental model" project
# documented in internal ticket #1555. This is not designed to be used outside
# of the XSweep classes. N.B. this should probably be generalized for non
# flat detectors, or composite detectors constructed from a number of flat
# elements.

import math
import pycbf
from scitbx import matrix
from dxtbx_model_ext import Panel, VirtualPanel, Detector
from dxtbx_model_ext import SimplePxMmStrategy, ParallaxCorrectedPxMmStrategy
from detector_helpers import detector_helper_sensors
from detector_helpers import find_undefined_value

class PanelTreeNode(object):
  ''' A class to wrap a virtual panel object and project methods for setting
  parent and local coordinate frames for the hierarchical detector model. '''

  def __init__(self, node=None, parent=None):
    ''' If no node is set, then create a virtual panel node. '''
    if node is None:
      node = VirtualPanel()
    self._node = node
    self._parent = parent

  def __getattr__(self, attr):
    ''' Inherit the interface from the node. '''
    return getattr(self._node, attr)

  def parent(self):
    ''' Return the parent. '''
    return self._parent

  def root(self):
    ''' Return the root. '''
    if self._parent:
      return self._parent.root()
    return self

  def apply_transformation(self, t):
    ''' Apply a transformation to the current matrix. '''
    from scitbx import matrix
    tg = matrix.sqr(self.get_transformation_matrix())
    tgt = (t * tg).transpose()
    self.set_frame(tgt[0:3], tgt[4:7], tgt[12:15])

  def set_local_frame(self, fast_axis, slow_axis, origin):
    ''' Set the local frame. '''
    from scitbx import matrix

    # Check if the parent is None
    if self.parent() is None:
      self.set_frame(fast_axis, slow_axis, origin)

    # Normalize the axes
    fast_axis = matrix.col(fast_axis).normalize()
    slow_axis = matrix.col(slow_axis).normalize()
    normal = fast_axis.cross(slow_axis)

    # Get the parent transformation matrix
    tp = matrix.sqr(self.parent().get_transformation_matrix())

    # Get the local transformation matrix
    tl = matrix.sqr(
     fast_axis.elems + (0,) +
     slow_axis.elems + (0,) +
     normal.elems    + (0,) +
     tuple(origin)   + (1,)).transpose()

    # Set the current frame
    tgt = (tp * tl).transpose()
    self.set_frame(tgt[0:3], tgt[4:7], tgt[12:15])

  def get_local_d_matrix(self):
    ''' Get the local d matrix. '''
    from scitbx import matrix
    if self.parent() is None:
      return self.get_d_matrix()
    tl = matrix.sqr(self.get_local_transformation_matrix()).transpose()
    return matrix.sqr(tl[0:3] + tl[4:7] + tl[12:15]).elems

  def get_transformation_matrix(self):
    ''' Get the transformation matrix. '''
    from scitbx import matrix
    return matrix.sqr(
     self.get_fast_axis() + (0,) +
     self.get_slow_axis() + (0,) +
     self.get_normal()    + (0,) +
     self.get_origin()    + (1,)).transpose().elems

  def get_local_transformation_matrix(self):
    ''' Get the local transformation matrix. '''
    from scitbx import matrix
    if self.parent() is None:
      return self.get_transformation_matrix()
    tp = matrix.sqr(self.parent().get_transformation_matrix())
    tg = matrix.sqr(self.get_transformation_matrix())
    return (tp.inverse() * tg).elems

  def get_local_fast_axis(self):
    ''' Get the local fast axis vector. '''
    from scitbx import matrix
    dl = matrix.sqr(self.get_local_transformation_matrix()).transpose()
    return dl[0:3]

  def get_local_slow_axis(self):
    ''' Get the local slow axis vector. '''
    from scitbx import matrix
    dl = matrix.sqr(self.get_local_transformation_matrix()).transpose()
    return dl[4:7]

  def get_local_origin(self):
    ''' Get the local origin vector. '''
    from scitbx import matrix
    dl = matrix.sqr(self.get_local_transformation_matrix()).transpose()
    return dl[12:15]

  def get_local_normal(self):
    ''' Get the local normal vector. '''
    from scitbx import matrix
    dl = matrix.sqr(self.get_local_transformation_matrix()).transpose()
    return dl[8:11]


class PanelGroup(VirtualPanel):
  ''' A class providing an iterface to a group of panels.

  This class is the basis for the construction of a detector hierarchy. The
  class inherits from VirtualPanel which has a C++ implementation providing
  the methods to manipulate the virtual detector plane. This class holds
  a reference to a list of children and allows for propagating the panel
  coordinate frames through the hierarchy.

  '''
  def __init__(self, parent=None):
    ''' Initialise the list of children to an empty list. '''
    VirtualPanel.__init__(self)
    self._parent = parent
    self._children = []

  def set_parent_frame(self, fast_axis, slow_axis, origin):
    ''' Set the parent frame.

    Sets its own parent plane and then, after updating its global
    frame, propagates the frame down to its children.

    :param fast_axis: The fast axis of the virtual detector plane
    :param slow_axis: The slow axis of the virtual detector plane
    :param origin: The origin vector to the virtual detector plane

    '''
    VirtualPanel.set_parent_frame(self, fast_axis, slow_axis, origin)
    for child in self:
      child.set_parent_frame(
          self.get_fast_axis(),
          self.get_slow_axis(),
          self.get_origin())

  def set_local_frame(self, fast_axis, slow_axis, origin):
    ''' Set the local frame.

    Sets its own local plane and then, after updating its global
    frame, propagates the frame down to its children.

    :param fast_axis: The fast axis of the virtual detector plane
    :param slow_axis: The slow axis of the virtual detector plane
    :param origin: The origin vector to the virtual detector plane

    '''
    VirtualPanel.set_local_frame(self, fast_axis, slow_axis, origin)
    for child in self:
      child.set_parent_frame(
          self.get_fast_axis(),
          self.get_slow_axis(),
          self.get_origin())

  def set_frame(self, fast_axis, slow_axis, origin):
    ''' Set the local frame.

    Sets its own global plane and then, after updating its local
    frame, propagates the frame down to its children.

    :param fast_axis: The fast axis of the virtual detector plane
    :param slow_axis: The slow axis of the virtual detector plane
    :param origin: The origin vector to the virtual detector plane

    '''
    VirtualPanel.set_frame(self, fast_axis, slow_axis, origin)
    for child in self:
      child.set_parent_frame(
          self.get_fast_axis(),
          self.get_slow_axis(),
          self.get_origin())

  def add_group(self):
    ''' Add a new group to this group.

    :returns: A new panel group
    :rtype: PanelGroup

    '''
    group = PanelGroup(self)
    group.set_parent_frame(
        self.get_fast_axis(),
        self.get_slow_axis(),
        self.get_origin())
    self._children.append(group)
    return group

  def add_panel(self, panel):
    ''' Add a new panel to this group.

    :returns: A new panel
    :rtype: Panel

    '''
    assert(isinstance(panel, Panel))
    assert(not hasattr(panel, "parent") or panel.parent is None)
    assert(panel in self.root()._container)
    panel.parent = self
    panel.set_parent_frame(
        self.get_fast_axis(),
        self.get_slow_axis(),
        self.get_origin())
    self._children.append(panel)
    return panel

  def __getitem__(self, index):
    ''' Get the child at the given index.

    :param index: The index of the child to get.

    '''
    return self._children[index]

  def remove(self, item):
    ''' Remove a child from the tree. '''
    return self._children.remove(item)

  def index(self, item):
    ''' Get the index of a child. '''
    return self._children.index(item)

  def parent(self):
    ''' Return the parent. '''
    return self._parent

  def root(self):
    ''' Return the root. '''
    if self._parent:
      return self._parent.root()
    return self

  def children(self):
    ''' Return an iterator to the list children. '''
    return iter(self._children)

  def reverse(self):
    ''' Return a reverse iterator to the list of children. '''
    return reversed(self._children)

  def __iter__(self):
    ''' Iterate through the children. '''
    return self.children()

  def __len__(self):
    ''' Get the length of the list of children. '''
    return len(self._children)

  def __eq__(self, other):
    ''' Check that this is equal to another group. '''
    if VirtualPanel.__eq__(self, other):
      if len(self) != len(other):
        return False
      return all(a == b for a, b in zip(self, other))
    else:
      return False

  def __ne__(self, other):
    ''' Check that this is not equal to another group. '''
    return not self.__eq__(other)

  def is_similar_to(self, other):
    ''' Check is the hierarchy is similar. '''
    if len(self) != len(other):
      return False
    return all(a.is_similar_to(b) for a, b in zip(self, other))

  def to_dict(self):
    ''' Convert the panel group to a dictionary. '''
    d = VirtualPanel.to_dict(self)
    children = []
    for c in self._children:
      if isinstance(c, Panel):
        idx = None
        for i, p in enumerate(self.root()._container):
          if c.is_(p):
            idx = i
            break
        assert(idx is not None)
        children.append({ "panel" : idx })
      else:
        children.append(c.to_dict())
    d['children'] = children
    return d


class PanelGroupRoot(PanelGroup):
  ''' The top level Detector model.

  This class is derived from the panel group class but provides some
  additional convenience methods for navigating the panel hierarchy.

  '''
  def __init__(self, container=None):
    super(PanelGroupRoot, self).__init__()
    self._container = container

  def iter_panels(self):
    ''' Iterate through just the panels depth-first. '''
    for obj in self.iter_preorder():
      if isinstance(obj, Panel):
        yield obj

  def iter_preorder(self):
    ''' Iterate through the groups and panels depth-first. '''
    stack = [self]
    while (len(stack) > 0):
      node = stack.pop()
      yield node
      if isinstance(node, PanelGroup):
        for child in node.reverse():
          stack.append(child)

  def iter_levelorder(self):
    ''' Iterate through the groups and panels depth-first. '''
    from collections import deque
    queue = deque([self])
    while (len(queue) > 0):
      node = queue.popleft()
      yield node
      if isinstance(node, PanelGroup):
        for child in node:
          queue.append(child)


def detector__getstate__(self):
  ''' Get the state for pickling. '''
  version = 1
  return (version, self.__dict__, [p for p in self])

def detector__setstate__(self, state):
  ''' Set the state from pickling. '''
  assert(len(state) == 3)
  assert(state[0] == 1)
  self.__dict__.update(state[1])
  for p in state[2]:
    self.add_panel(p)


Detector.__getstate__ = detector__getstate__
Detector.__setstate__ = detector__setstate__
Detector.__getstate_manages_dict__ = 1


class HierarchicalDetector(Detector):
  ''' The detector class. '''
  def __init__(self, panel=None):
    ''' Construct the detector. '''
    if panel == None:
      super(HierarchicalDetector, self).__init__()
    else:
      super(HierarchicalDetector, self).__init__(panel)
    self._root = PanelGroupRoot(self)

  def hierarchy(self):
    ''' Return the hierarchy. '''
    return self._root

  def is_similar_to(self, rhs):
    ''' Check if the detectors are similar (i.e. only differ in terms of
    things that could change per experiment. '''
    if not isinstance(rhs, HierarchicalDetector):
      return False
    return (self._root.is_similar_to(rhs._root) and
            super(HierarchicalDetector, self).is_similar_to(rhs))

  def __eq__(self, rhs):
    ''' Check that this is equal to another group. '''
    if not isinstance(rhs, HierarchicalDetector):
      return False
    return (self._root == rhs._root and
            super(HierarchicalDetector, self).__eq__(rhs))

  def __ne__(self, other):
    ''' Check that this is not equal to another group. '''
    return not self.__eq__(other)

  def to_dict(self):
    ''' Return the class as a dictionary. '''
    d = Detector.to_dict(self)
    d['hierarchy'] = self._root.to_dict()
    return d

  @staticmethod
  def from_dict(d):
    ''' Convert a dictionary to a detector model. '''
    obj = HierarchicalDetector()
    panels = d['panels']
    for p in panels:
      obj.add_panel(Panel.from_dict(p))
    if 'hierarchy' in d:
      HierarchicalDetector._group_or_panel_from_dict(
        d['hierarchy'], obj.hierarchy())
    return obj

  @staticmethod
  def _group_or_panel_from_dict(d, obj):
    ''' Convert a dictionary to a group or panel and add to model. '''
    base = VirtualPanel.from_dict(d)
    obj.set_type(base.get_type())
    obj.set_name(base.get_name())
    obj.set_local_frame(
      base.get_local_fast_axis(),
      base.get_local_slow_axis(),
      base.get_local_origin())
    for child in d['children']:
      if 'panel' in child:
        index = child['panel']
        obj.add_panel(obj.root()._container[index])
      else:
        HierarchicalDetector._group_or_panel_from_dict(child, obj.add_group())

  def __getstate__(self):
    ''' Get the state for pickling. '''
    version = 1
    panels_list = []
    for p1 in self:
      p3 = p1
      for p2 in self.hierarchy().iter_panels():
        if p1.is_(p2):
          p3 = p2
          break
      panels_list.append(p3)
    return (version, self.__dict__, panels_list)

  def __setstate__(self, state):
    ''' Set the state from pickling. '''
    assert(len(state) == 3)
    assert(state[0] == 1)
    self.__dict__.update(state[1])
    for p in state[2]:
      self.add_panel_by_pointer(p)
    assert(all(p in self for p in self.hierarchy().iter_panels()))

  # Manage the object dictionary
  __getstate_manages_dict__ = 1

class detector_factory:
  '''A factory class for detector objects, which will encapsulate standard
  detector designs to make it a little easier to get started with these. In
  cases where a CBF image is provided a full description can be used, in
  other cases assumptions will be made about the experiment configuration.
  In all cases information is provided in the CBF coordinate frame.'''

  def __init__(self):
    pass

  @staticmethod
  def make_detector(stype, fast_axis, slow_axis, origin,
                    pixel_size, image_size, trusted_range = (0.0, 0.0),
                    px_mm=None, name="Panel", thickness=0.0, material='',
                    mu=0.0):
    """Ensure all types are correct before creating c++ detector class."""

    if px_mm is None:
      px_mm = SimplePxMmStrategy()
    try:
      d = Detector()
      p = d.add_panel()
      p.set_type(str(stype))
      p.set_name(str(name))
      p.set_local_frame(
          tuple(map(float, fast_axis)),
          tuple(map(float, slow_axis)),
          tuple(map(float, origin)))
      p.set_pixel_size(tuple(map(float, pixel_size)))
      p.set_image_size(tuple(map(int, image_size)))
      p.set_trusted_range(tuple(map(float, trusted_range)))
      p.set_thickness(thickness)
      p.set_material(material)
      p.set_px_mm_strategy(px_mm)
    except Exception, e:
      print e
      raise e
    return d

  @staticmethod
  def simple(sensor, distance, beam_centre, fast_direction, slow_direction,
             pixel_size, image_size, trusted_range = (0.0, 0.0), mask = [],
             px_mm=None, mu=0.0):
    '''Construct a simple detector at a given distance from the sample
    along the direct beam presumed to be aligned with -z, offset by the
    beam centre - the directions of which are given by the fast and slow
    directions, which are themselves given as +x, +y, -x, -y. The pixel
    size is given in mm in the fast and slow directions and the image size
    is given in pixels in the same order. Everything else is the same as
    for the main reference frame.'''

    assert(fast_direction in ['-x', '+y', '+x', '-y'])
    assert(slow_direction in ['-x', '+y', '+x', '-y'])

    assert(fast_direction[1] != slow_direction[1])

    direction_map = {
        '+x':(1.0, 0.0, 0.0),
        '-x':(-1.0, 0.0, 0.0),
        '+y':(0.0, 1.0, 0.0),
        '-y':(0.0, -1.0, 0.0)
        }

    fast = matrix.col(direction_map[fast_direction])
    slow = matrix.col(direction_map[slow_direction])

    origin = matrix.col((0, 0, -1)) * distance - \
             fast * beam_centre[0] - slow * beam_centre[1]

    detector = detector_factory.make_detector(
        detector_factory.sensor(sensor),
        fast, slow, origin, pixel_size, image_size, trusted_range, px_mm, mu=mu)
    detector[0].mask = mask
    return detector

  @staticmethod
  def two_theta(sensor, distance, beam_centre, fast_direction,
                slow_direction, two_theta_direction, two_theta_angle,
                pixel_size, image_size, trusted_range = (0.0, 0.0),
                mask = [], px_mm = None, mu=0.0):
    '''Construct a simple detector at a given distance from the sample
    along the direct beam presumed to be aligned with -z, offset by the
    beam centre - the directions of which are given by the fast and slow
    directions, which are themselves given as +x, +y, -x, -y. The pixel
    size is given in mm in the fast and slow directions and the image size
    is given in pixels in the same order. Everything else is the same as
    for the main reference frame. Also given are the direction of the
    two-theta axis and the angle in degrees by which the detector is
    moved.'''

    assert(fast_direction in ['-x', '+y', '+x', '-y'])
    assert(slow_direction in ['-x', '+y', '+x', '-y'])
    assert(two_theta_direction in ['-x', '+y', '+x', '-y'])

    assert(fast_direction[1] != slow_direction[1])

    direction_map = {
        '+x':(1.0, 0.0, 0.0),
        '-x':(-1.0, 0.0, 0.0),
        '+y':(0.0, 1.0, 0.0),
        '-y':(0.0, -1.0, 0.0)
        }

    fast = matrix.col(direction_map[fast_direction])
    slow = matrix.col(direction_map[slow_direction])

    origin = matrix.col((0, 0, -1)) * distance - \
             fast * beam_centre[0] - slow * beam_centre[1]

    two_theta = matrix.col(direction_map[two_theta_direction])

    R = two_theta.axis_and_angle_as_r3_rotation_matrix(two_theta_angle,
                                                       deg = True)

    detector = detector_factory.make_detector(
        detector_factory.sensor(sensor),
        (R * fast), (R * slow), (R * origin), pixel_size,
        image_size, trusted_range, px_mm, mu=mu)

    detector.mask = mask
    return detector

  @staticmethod
  def complex(sensor, origin, fast, slow, pixel, size,
              trusted_range = (0.0, 0.0), px_mm = None):
    '''A complex detector model, where you know exactly where everything
    is. This is useful for implementation of the Rigaku Saturn header
    format, as that is exactly what is in there. Origin, fast and slow are
    vectors in the CBF reference frame, pixel is the dimensions as a tuple
    as is size.'''

    assert(len(origin) == 3)
    assert(len(fast) == 3)
    assert(len(slow) == 3)
    assert(len(pixel) == 2)
    assert(len(size) == 2)

    return detector_factory.make_detector(
            detector_factory.sensor(sensor),
            fast, slow, origin, pixel, size, trusted_range, px_mm)

  @staticmethod
  def imgCIF(cif_file, sensor):
    '''Initialize a detector model from an imgCIF file.'''

    cbf_handle = pycbf.cbf_handle_struct()
    cbf_handle.read_file(cif_file, pycbf.MSG_DIGEST)

    cbf_detector = cbf_handle.construct_detector(0)

    pixel = (cbf_detector.get_inferred_pixel_size(1),
             cbf_detector.get_inferred_pixel_size(2))

    # FIXME can probably simplify the code which follows below by
    # making proper use of cctbx vector calls - should not be as
    # complex as it appears to be...

    origin = tuple(cbf_detector.get_pixel_coordinates(0, 0))
    fast = cbf_detector.get_pixel_coordinates(0, 1)
    slow = cbf_detector.get_pixel_coordinates(1, 0)

    dfast = [fast[j] - origin[j] for j in range(3)]
    dslow = [slow[j] - origin[j] for j in range(3)]

    lfast = math.sqrt(sum([dfast[j] * dfast[j] for j in range(3)]))
    lslow = math.sqrt(sum([dslow[j] * dslow[j] for j in range(3)]))

    fast = tuple([dfast[j] / lfast for j in range(3)])
    slow = tuple([dslow[j] / lslow for j in range(3)])

    size = tuple(reversed(cbf_handle.get_image_size(0)))

    try:
      underload = find_undefined_value(cbf_handle)
      overload = cbf_handle.get_overload(0)
      trusted_range = (underload, overload)
    except: # intentional
      trusted_range = (0.0, 0.0)

    cbf_detector.__swig_destroy__(cbf_detector)
    del(cbf_detector)

    # Get the sensor type
    dtype = detector_factory.sensor(sensor)

    # If the sensor type is PAD then create the detector with a
    # parallax corrected pixel to millimeter function
    #if dtype == detector_helper_sensors.SENSOR_PAD:
      #px_mm = ParallaxCorrectedPxMmStrategy(0.252500934883)
    #else:
    px_mm = SimplePxMmStrategy()

    return detector_factory.make_detector(
              dtype, fast, slow, origin, pixel, size,
              trusted_range, px_mm)

  @staticmethod
  def imgCIF_H(cbf_handle, sensor):
    '''Initialize a detector model from an imgCIF file handle, where it
    is assumed that the file has already been read.'''

    cbf_detector = cbf_handle.construct_detector(0)

    pixel = (cbf_detector.get_inferred_pixel_size(1),
             cbf_detector.get_inferred_pixel_size(2))

    # FIXME can probably simplify the code which follows below by
    # making proper use of cctbx vector calls - should not be as
    # complex as it appears to be...

    origin = tuple(cbf_detector.get_pixel_coordinates(0, 0))
    fast = cbf_detector.get_pixel_coordinates(0, 1)
    slow = cbf_detector.get_pixel_coordinates(1, 0)

    dfast = [fast[j] - origin[j] for j in range(3)]
    dslow = [slow[j] - origin[j] for j in range(3)]

    lfast = math.sqrt(sum([dfast[j] * dfast[j] for j in range(3)]))
    lslow = math.sqrt(sum([dslow[j] * dslow[j] for j in range(3)]))

    fast = tuple([dfast[j] / lfast for j in range(3)])
    slow = tuple([dslow[j] / lslow for j in range(3)])

    size = tuple(reversed(cbf_handle.get_image_size(0)))

    try:
      underload = find_undefined_value(cbf_handle)
      overload = cbf_handle.get_overload(0)
      trusted_range = (underload, overload)
    except: # intentional
      trusted_range = (0.0, 0.0)

    cbf_detector.__swig_destroy__(cbf_detector)
    del(cbf_detector)

    return detector_factory.make_detector(
                  detector_factory.sensor(sensor),
                  fast, slow, origin, pixel, size, trusted_range)

  @staticmethod
  def sensor(name):
    '''Return the correct sensor token for a given name, for example:

    ccd, CCD
    image_plate, IMAGE_PLATE
    pad, PAD

    to the appropriate static token which will be used as a handle
    everywhere else in this. Also allow existing token to be passed in.'''

    if detector_helper_sensors.check_sensor(name):
      return name

    if name.upper() == 'PAD':
      return detector_helper_sensors.SENSOR_PAD
    elif name.upper() == 'CCD':
      return detector_helper_sensors.SENSOR_CCD
    elif name.upper() == 'IMAGE_PLATE':
      return detector_helper_sensors.SENSOR_IMAGE_PLATE
    elif name.upper() == 'UNKNOWN':
      return detector_helper_sensors.SENSOR_UNKNOWN

    raise RuntimeError, 'name %s not known' % name
