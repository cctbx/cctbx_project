
from __future__ import division
from dxtbx.model import HierarchicalDetector, PanelTreeNode
from dxtbx.model import PanelGroup

class TestPanelTreeNode(object):

  def __init__(self):
    from scitbx import matrix
    self.root = PanelTreeNode()
    self.node = PanelTreeNode(parent=self.root)

    root_fast = matrix.col((1, 1, 0)).normalize()
    root_slow = matrix.col((-1, 1, 0)).normalize()
    root_origin = (10, 10, 10)

    node_fast = matrix.col((1, 0, 0))
    node_slow = matrix.col((0, 1, 0))
    node_origin = (0, 0, 0)

    self.root.set_frame(
      root_fast,
      root_slow,
      root_origin)

    self.node.set_frame(
      node_fast,
      node_slow,
      node_origin)

  def run(self):
    self.tst_attributes()
    self.tst_parent_and_root()
    self.tst_get_local_frame()
    self.tst_set_local_frame()
    self.tst_set_parent_frame()

  def tst_attributes(self):
    assert(self.node.get_name() == '')
    assert(self.node.get_type() == '')
    assert(self.node.get_fast_axis() == (1, 0, 0))
    assert(self.node.get_slow_axis() == (0, 1, 0))
    assert(self.node.get_normal() == (0, 0, 1))
    assert(self.node.get_origin() == (0, 0, 0))
    print 'OK'

  def tst_parent_and_root(self):
    assert(self.node.parent() == self.root)
    assert(self.node.root() == self.root)
    assert(self.root.parent() == None)
    assert(self.root.root() == self.root)
    print 'OK'

  def tst_get_local_frame(self):
    from scitbx import matrix
    from math import sqrt
    tm = self.node.get_transformation_matrix()
    assert(tm == (
      1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1, 0,
      0, 0, 0, 1))


    tm = self.node.get_local_transformation_matrix()
    ex_fast = matrix.col((1, -1, 0)).normalize().elems
    ex_slow = matrix.col((1, 1, 0)).normalize().elems
    ex_norm = matrix.col((0, 0, 1)).normalize().elems
    ex_orig = matrix.col((-sqrt(200), 0, -10)).elems
    expected = matrix.sqr(
      ex_fast + (0,) +
      ex_slow + (0,) +
      ex_norm + (0,) +
      ex_orig + (1,)).transpose().elems

    eps = 1e-7
    assert(len(tm) == len(expected))
    assert(all(abs(a - b) < eps for a, b in zip(tm, expected)))
    print 'OK'

    fast = self.node.get_local_fast_axis()
    slow = self.node.get_local_slow_axis()
    norm = self.node.get_local_normal()
    orig = self.node.get_local_origin()
    assert(all(abs(a-b) < eps for a, b in zip(fast, ex_fast)))
    assert(all(abs(a-b) < eps for a, b in zip(slow, ex_slow)))
    assert(all(abs(a-b) < eps for a, b in zip(norm, ex_norm)))
    assert(all(abs(a-b) < eps for a, b in zip(orig, ex_orig)))
    print 'OK'

    ld = self.node.get_local_d_matrix()
    expected = fast + slow + orig
    assert(all(abs(a-b) < eps for a, b in zip(ld, expected)))
    print 'OK'


  def tst_set_local_frame(self):
    from scitbx import matrix
    fast = matrix.col((-1, -1, 0)).normalize()
    slow = matrix.col((-1, 1, 0)).normalize()
    origin = (0, 0, 0)

    self.node.set_local_frame(fast, slow, origin)

    fast_get = self.node.get_local_fast_axis()
    slow_get = self.node.get_local_slow_axis()
    origin_get = self.node.get_local_origin()

    eps = 1e-7
    assert(all(abs(a-b) < eps for a, b in zip(fast, fast_get)))
    assert(all(abs(a-b) < eps for a, b in zip(slow, slow_get)))
    assert(all(abs(a-b) < eps for a, b in zip(origin, origin_get)))
    print 'OK'

    ex_fast = matrix.col((0, -1, 0))
    ex_slow = matrix.col((-1, 0, 0))
    ex_orig = matrix.col((10, 10, 10))
    fast = self.node.get_fast_axis()
    slow = self.node.get_slow_axis()
    orig = self.node.get_origin()
    assert(all(abs(a-b) < eps for a, b in zip(fast, ex_fast)))
    assert(all(abs(a-b) < eps for a, b in zip(slow, ex_slow)))
    assert(all(abs(a-b) < eps for a, b in zip(orig, ex_orig)))
    print 'OK'

  def tst_set_parent_frame(self):
    from scitbx import matrix
    ex_local_fast = matrix.col((-1, -1, 0)).normalize()
    ex_local_slow = matrix.col((-1, 1, 0)).normalize()
    ex_local_origin = (0, 0, 0)

    fast = matrix.col((1, 0, 0)).normalize()
    slow = matrix.col((0, 1, 0)).normalize()
    origin = matrix.col((0, 0, 0))
    normal = fast.cross(slow)

    # Get the new parent transformation matrix
    tp2 = matrix.sqr(
     fast.elems + (0,) +
     slow.elems + (0,) +
     normal.elems + (0,) +
     origin.elems + (1,)).transpose()

    tp1 = matrix.sqr(self.node.parent().get_transformation_matrix())
    td = tp2 * tp1.inverse()
    self.node.apply_transformation(td)

    local_fast = self.node.get_local_fast_axis()
    local_slow = self.node.get_local_slow_axis()
    local_origin = self.node.get_local_origin()

    eps = 1e-7
    # assert(all(abs(a-b) < eps for a, b in zip(local_fast, local_fast)))
    # assert(all(abs(a-b) < eps for a, b in zip(local_slow, local_slow)))
    # assert(all(abs(a-b) < eps for a, b in zip(local_origin, local_origin)))
    # print 'OK'
    # FIXME: This test checks that the values are identical to themselves.
    # Probably not what was intended.

    ex_fast = matrix.col((-1, -1, 0)).normalize()
    ex_slow = matrix.col((-1, 1, 0)).normalize()
    ex_orig = matrix.col((0, 0, 0))
    fast = self.node.get_fast_axis()
    slow = self.node.get_slow_axis()
    orig = self.node.get_origin()

    assert(all(abs(a-b) < eps for a, b in zip(fast, ex_fast)))
    assert(all(abs(a-b) < eps for a, b in zip(slow, ex_slow)))
    assert(all(abs(a-b) < eps for a, b in zip(orig, ex_orig)))
    print 'OK'


class TestPanelGroup(object):

  def __init__(self):

    # Create the panel group
    self.root = PanelGroup()
    self.g1 = self.root.add_group()
    self.g2 = self.g1.add_group()

  def run(self):
    self.tst_set_frame()
    self.tst_set_local_frame()

  def tst_set_frame(self):
    from scitbx import matrix

    # Set the root frame
    root_fast = matrix.col((1, 1, 0)).normalize()
    root_slow = matrix.col((-1, 1, 0)).normalize()
    root_origin = (10, 10, 10)
    self.root.set_frame(
      root_fast,
      root_slow,
      root_origin)

    # Assert that all frames equal this
    ex = matrix.sqr(
      root_fast.elems +
      root_slow.elems +
      root_origin).transpose().elems

    assert(all(a-b) < eps for a, b in zip(ex, self.root.get_d_matrix()))
    assert(all(a-b) < eps for a, b in zip(ex, self.g1.get_d_matrix()))
    assert(all(a-b) < eps for a, b in zip(ex, self.g2.get_d_matrix()))

    print 'OK'

  def tst_set_local_frame(self):
    from scitbx import matrix

    # Set the root frame
    root_fast = matrix.col((1, 1, 0)).normalize()
    root_slow = matrix.col((-1, 1, 0)).normalize()
    root_origin = (10, 10, 10)
    self.root.set_local_frame(
      root_fast,
      root_slow,
      root_origin)

    # Assert that all frames equal this
    ex = matrix.sqr(
      root_fast.elems +
      root_slow.elems +
      root_origin).transpose().elems

    assert(all(a-b) < eps for a, b in zip(ex, self.root.get_d_matrix()))
    assert(all(a-b) < eps for a, b in zip(ex, self.g1.get_d_matrix()))
    assert(all(a-b) < eps for a, b in zip(ex, self.g2.get_d_matrix()))

    print 'OK'

    # Set the G1 local frame
    g1_fast = matrix.col((1, -1, 0))
    g1_slow = matrix.col((1, 1, 0))
    g1_orig = matrix.col((0, 0, 0))
    self.g1.set_local_frame(g1_fast, g1_slow, g1_orig)

    ex_g1_fast = (1, 0, 0)
    ex_g1_slow = (0, 1, 0)
    ex_g1_orig = (10, 10, 10)
    ex_g1 = matrix.col(ex_g1_fast + ex_g1_slow + ex_g1_orig).transpose().elems

    assert(all(a-b) < eps for a, b in zip(ex, self.root.get_d_matrix()))
    assert(all(a-b) < eps for a, b in zip(ex_g1, self.g1.get_d_matrix()))
    assert(all(a-b) < eps for a, b in zip(ex_g1, self.g2.get_d_matrix()))


class Test2:

  def __init__(self):
    detector = HierarchicalDetector()
    panel1 = detector.add_panel()
    panel1.set_name("P1")
    panel1.set_type("P")

    panel2 = detector.add_panel()
    panel2.set_name("P2")
    panel2.set_type("P")

    panel3 = detector.add_panel()
    panel3.set_name("P3")
    panel3.set_type("P")

    panel4 = detector.add_panel()
    panel4.set_name("P4")
    panel4.set_type("P")

    root = detector.hierarchy()
    root.set_name("D1")
    root.set_type("D")

    quad1 = root.add_group()
    quad1.set_name("Q1")
    quad1.set_type("Q")
    quad1.add_panel(panel1)
    quad1.add_panel(panel2)

    quad2 = root.add_group()
    quad2.set_name("Q2")
    quad2.set_type("Q")
    quad2.add_panel(panel3)
    quad2.add_panel(panel4)

    self.detector = detector

  def run(self):
    self.tst_iterate_and_index()
    self.tst_flat()
    #self.tst_get_uninitialized_D_matrix()
    self.tst_get_valid_D_matrix()
    self.tst_copy_and_reference()

  def tst_flat(self):
    ''' Test the flat hierarchy. '''

    expected_types = ['P', 'P', 'P', 'P']
    expected_names = ['P1', 'P2', 'P3', 'P4']
    names = []
    types = []
    for p in self.detector:
      names.append(p.get_name())
      types.append(p.get_type())
    assert(all(n == en for n, en in zip(names, expected_names)))
    assert(all(t == et for t, et in zip(types, expected_types)))

    print 'OK'

  def tst_iterate_and_index(self):
    ''' Test iteration and indexing through the detector in various ways. '''

    # Iterate through the detector's children and check output
    expected_types = ['Q', 'Q',]
    expected_names = ['Q1', 'Q2']
    names = []
    types = []
    for p in self.detector.hierarchy():
      names.append(p.get_name())
      types.append(p.get_type())
    assert(all(n == en for n, en in zip(names, expected_names)))
    assert(all(t == et for t, et in zip(types, expected_types)))

    # Iterate through the detector's children in reverse and check output
    expected_types = ['Q', 'Q',]
    expected_names = ['Q2', 'Q1']
    names = []
    types = []
    for p in self.detector.hierarchy().reverse():
      names.append(p.get_name())
      types.append(p.get_type())
    assert(all(n == en for n, en in zip(names, expected_names)))
    assert(all(t == et for t, et in zip(types, expected_types)))

    # Use an index to access the detector children
    assert(len(self.detector.hierarchy()) == 2)
    group = self.detector.hierarchy()[1]
    assert(group.get_name() == 'Q2' and group.get_type() == 'Q')
    assert(len(group) == 2)
    panel = group[0]
    assert(panel.get_name() == 'P3' and panel.get_type() == 'P')

    # Iterate through the tree pre-order and check output
    expected_types = ['D', 'Q', 'P', 'P', 'Q', 'P', 'P']
    expected_names = ['D1', 'Q1', 'P1', 'P2', 'Q2', 'P3', 'P4']
    names = []
    types = []
    for p in self.detector.hierarchy().iter_preorder():
      names.append(p.get_name())
      types.append(p.get_type())
    assert(all(n == en for n, en in zip(names, expected_names)))
    assert(all(t == et for t, et in zip(types, expected_types)))

    # Iterate through the tree level-order and check output
    expected_types = ['D', 'Q', 'Q', 'P', 'P', 'P', 'P']
    expected_names = ['D1', 'Q1', 'Q2', 'P1', 'P2', 'P3', 'P4']
    names = []
    types = []
    for p in self.detector.hierarchy().iter_levelorder():
      names.append(p.get_name())
      types.append(p.get_type())
    assert(all(n == en for n, en in zip(names, expected_names)))
    assert(all(t == et for t, et in zip(types, expected_types)))

    # Iterate through the panels in pre-order and check output
    expected_types = ['P', 'P', 'P', 'P']
    expected_names = ['P1', 'P2', 'P3', 'P4']
    names = []
    types = []
    for p in self.detector.hierarchy().iter_panels():
      names.append(p.get_name())
      types.append(p.get_type())
    assert(all(n == en for n, en in zip(names, expected_names)))
    assert(all(t == et for t, et in zip(types, expected_types)))

    print 'OK'

  def tst_get_uninitialized_D_matrix(self):
    ''' Try to get bad D matrix and check that an exception is thrown. '''
    panels = self.detector.panels()
    for p in panels:
      try:
        p.get_D_matrix()
        assert(False)
      except Exception:
        pass

    print 'OK'

  def tst_get_valid_D_matrix(self):
    ''' Setup the hierarchy of frames and check it's all consistent. '''
    from scitbx import matrix

    # Set a valid frame for the top level detector
    self.detector.hierarchy().set_local_frame(
        (1, 0, 0),     # Fast axis
        (0, 1, 0),     # Slow axis
        (0, 0, 100))   # Origin

    # Check that all sub groups have the same frame and that we can get
    # a valid D matrix
    for obj in self.detector.hierarchy().iter_preorder():
      fast = matrix.col(obj.get_fast_axis())
      slow = matrix.col(obj.get_slow_axis())
      orig = matrix.col(obj.get_origin())
      assert(abs(fast - matrix.col((1, 0, 0))) < 1e-7)
      assert(abs(slow - matrix.col((0, 1, 0))) < 1e-7)
      assert(abs(orig - matrix.col((0, 0, 100))) < 1e-7)
      D = obj.get_D_matrix()

    # Get the quadrants and set their frames
    q1, q2 = self.detector.hierarchy().children()
    q1.set_local_frame(
        (1, 1, 0),      # Fast axis relative to detector frame
        (-1, 1, 0),     # Slow axis relative to detector frame
        (10, 10, 0))    # Origin relative to detector frame
    q2.set_local_frame(
        (1, -1, 0),     # Fast axis relative to detector frame
        (1, 1, 0),      # Slow axis relative to detector frame
        (20, 20, 0))    # Origin relative to detector frame

    # Get the panels and set their frames
    p1, p2 = q1.children()
    p1.set_local_frame(
        (1, -1, 0),     # Fast axis relative to q1 frame
        (1, 1, 0),      # Slow axis relative to q1 frame
        (5, 0, 10))     # Origin relative to q1 frame
    p2.set_local_frame(
        (1, 1, 0),      # Fast axis relative to q1 frame
        (-1, 1, 0),     # Slow axis relative to q1 frame
        (0, 5, -10))    # Origin relative to q1 frame

    # Get the panels and set their frames
    p3, p4 = q2.children()
    p3.set_local_frame(
        (1, -1, 0),     # Fast axis relative to q2 frame
        (1, 1, 0),      # Slow axis relative to q2 frame
        (0, 5, -10))    # Origin relative to q2 frame
    p4.set_local_frame(
        (1, 1, 0),      # Fast axis relative to q2 frame
        (-1, 1, 0),     # Slow axis relative to q2 frame
        (5, 0, 10))     # Origin relative to q2 frame

    # Test the panel coordinate systems
    from math import sqrt
    eps = 1e-7
    p1_d0 = matrix.col((10.0 + sqrt(5.0**2 / 2), 10.0 + sqrt(5.0**2 / 2), 110))
    p2_d0 = matrix.col((10.0 - sqrt(5.0**2 / 2), 10.0 + sqrt(5.0**2 / 2), 90))
    p3_d0 = matrix.col((20.0 + sqrt(5.0**2 / 2), 20.0 + sqrt(5.0**2 / 2), 90))
    p4_d0 = matrix.col((20.0 + sqrt(5.0**2 / 2), 20.0 - sqrt(5.0**2 / 2), 110))
    p1_d1 = matrix.col((1, 0, 0))
    p2_d1 = matrix.col((0, 1, 0))
    p3_d1 = matrix.col((0, -1, 0))
    p4_d1 = matrix.col((1, 0, 0))
    p1_d2 = matrix.col((0, 1, 0))
    p2_d2 = matrix.col((-1, 0, 0))
    p3_d2 = matrix.col((1, 0, 0))
    p4_d2 = matrix.col((0, 1, 0))
    assert(abs(matrix.col(p1.get_origin()) - p1_d0) < eps)
    assert(abs(matrix.col(p2.get_origin()) - p2_d0) < eps)
    assert(abs(matrix.col(p3.get_origin()) - p3_d0) < eps)
    assert(abs(matrix.col(p4.get_origin()) - p4_d0) < eps)
    assert(abs(matrix.col(p1.get_fast_axis()) - p1_d1) < eps)
    assert(abs(matrix.col(p2.get_fast_axis()) - p2_d1) < eps)
    assert(abs(matrix.col(p3.get_fast_axis()) - p3_d1) < eps)
    assert(abs(matrix.col(p4.get_fast_axis()) - p4_d1) < eps)
    assert(abs(matrix.col(p1.get_slow_axis()) - p1_d2) < eps)
    assert(abs(matrix.col(p2.get_slow_axis()) - p2_d2) < eps)
    assert(abs(matrix.col(p3.get_slow_axis()) - p3_d2) < eps)
    assert(abs(matrix.col(p4.get_slow_axis()) - p4_d2) < eps)

    print 'OK'

  def tst_copy_and_reference(self):
    from copy import deepcopy

    # Get the detector hierarchy
    root = self.detector.hierarchy()

    # Get the panels in the hierarchy
    p1 = root[0][0]
    p2 = root[0][1]
    p3 = root[1][0]
    p4 = root[1][1]

    # Check panels are the same
    assert(p1.is_(self.detector[0]))
    assert(p2.is_(self.detector[1]))
    assert(p3.is_(self.detector[2]))
    assert(p4.is_(self.detector[3]))

    # Copy the detector
    new_detector = deepcopy(self.detector)

    # Check they're the same
    assert(new_detector == self.detector)

    # Add an offset to propagate
    root = new_detector.hierarchy()
    root.set_local_frame((1, 0, 0), (0, 1, 0), (0, 0, 10))

    # Get the panels in the hierarchy
    p1 = root[0][0]
    p2 = root[0][1]
    p3 = root[1][0]
    p4 = root[1][1]

    # Check panels are the same
    assert(p1.is_(new_detector[0]))
    assert(p2.is_(new_detector[1]))
    assert(p3.is_(new_detector[2]))
    assert(p4.is_(new_detector[3]))

    print 'OK'

if __name__ == '__main__':
  #test = Test()
  #test.run()
  test = TestPanelTreeNode()
  test.run()
  test = TestPanelGroup()
  test.run()
  test = Test2()
  test.run()
