from __future__ import absolute_import, division
import boost.python
import cctbx.uctbx
import cctbx.sgtbx
from dxtbx_model_ext import *
from dxtbx.model.beam import *
from dxtbx.model.goniometer import *
from dxtbx.model.detector import *
from dxtbx.model.scan import *
from dxtbx.model.profile import *

class DetectorAux(boost.python.injector, Detector):
  def iter_panels(self):
    ''' Iterate through just the panels depth-first. '''
    for obj in self.iter_preorder():
      if obj.is_panel():
        yield obj

  def iter_preorder(self):
    ''' Iterate through the groups and panels depth-first. '''
    stack = [self.hierarchy()]
    while (len(stack) > 0):
      node = stack.pop()
      yield node
      if node.is_group():
        for child in reversed(node):
          stack.append(child)

  def iter_levelorder(self):
    ''' Iterate through the groups and panels depth-first. '''
    from collections import deque
    queue = deque([self.hierarchy()])
    while (len(queue) > 0):
      node = queue.popleft()
      yield node
      if node.is_group():
        for child in node:
          queue.append(child)


class CrystalAux(boost.python.injector, Crystal):

  def show(self, show_scan_varying=False, out=None):
    from scitbx import matrix
    if out is None:
      import sys
      out = sys.stdout
    uc = self.get_unit_cell().parameters()
    sg = str(self.get_space_group().info())
    umat = matrix.sqr(self.get_U()).mathematica_form(format="% 5.4f",
                                         one_row_per_line=True).splitlines()
    bmat = matrix.sqr(self.get_B()).mathematica_form(format="% 5.4f",
                                         one_row_per_line=True).splitlines()
    amat = (matrix.sqr(self.get_U()) * matrix.sqr(self.get_B())).mathematica_form(format="% 5.4f",
                                         one_row_per_line=True).splitlines()

    msg =  ["Crystal:"]
    msg.append("    Unit cell: " + "(%5.3f, %5.3f, %5.3f, %5.3f, %5.3f, %5.3f)" % uc)
    msg.append("    Space group: " + sg)
    msg.append("    U matrix:  " + umat[0])
    msg.append("               " + umat[1])
    msg.append("               " + umat[2])
    msg.append("    B matrix:  " + bmat[0])
    msg.append("               " + bmat[1])
    msg.append("               " + bmat[2])
    msg.append("    A = UB:    " + amat[0])
    msg.append("               " + amat[1])
    msg.append("               " + amat[2])
    if self.num_scan_points > 0:
      msg.append("    A sampled at " + str(self.num_scan_points) \
                 + " scan points")
      if show_scan_varying:
        for i in range(self.num_scan_points):
          A = matrix.sqr(self.get_A_at_scan_point(i))
          B = matrix.sqr(self.get_B_at_scan_point(i))
          U = matrix.sqr(self.get_U_at_scan_point(i))
          uc = self.get_unit_cell_at_scan_point(i).parameters()
          umat = U.mathematica_form(format="% 5.4f",
                                    one_row_per_line=True).splitlines()
          bmat = B.mathematica_form(format="% 5.4f",
                                    one_row_per_line=True).splitlines()
          amat = A.mathematica_form(format="% 5.4f",
                                    one_row_per_line=True).splitlines()
          msg.append("  Scan point #%i:" %(i+1))
          msg.append("    Unit cell: " + "(%5.3f, %5.3f, %5.3f, %5.3f, %5.3f, %5.3f)" % uc)
          msg.append("    U matrix:  " + umat[0])
          msg.append("               " + umat[1])
          msg.append("               " + umat[2])
          msg.append("    B matrix:  " + bmat[0])
          msg.append("               " + bmat[1])
          msg.append("               " + bmat[2])
          msg.append("    A = UB:    " + amat[0])
          msg.append("               " + amat[1])
          msg.append("               " + amat[2])
    print >> out, "\n".join(msg)

  def __str__(self):
    from cStringIO import StringIO
    s = StringIO()
    msg = self.show(out=s)
    s.seek(0)
    return s.read()
