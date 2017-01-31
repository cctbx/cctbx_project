from __future__ import absolute_import, division
import boost.python
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
