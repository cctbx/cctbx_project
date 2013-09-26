#!/usr/bin/env python
# detector2.py
#   Copyright (C) 2011 Diamond Light Source
#
#   Author:
#     James Parkhurst
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#

from __future__ import division
from dxtbx_model_ext import PanelBase, Panel2

class PanelGroup(PanelBase):
    ''' A class providing an iterface to a group of panels.

    This class is the basis for the construction of a detector hierarchy. The
    class inherits from PanelBase which has a C++ implementation providing
    the methods to manipulate the virtual detector plane. This class holds
    a reference to a list of children and allows for propagating the panel
    coordinate frames through the hierarchy.

    '''
    def __init__(self):
        ''' Initialise the list of children to an empty list. '''
        PanelBase.__init__(self)
        self._children = []

    def set_parent_frame(self, fast_axis, slow_axis, origin):
        ''' Set the parent frame.

        Set's it's own parent plane and then, after updating it's global
        frame, propagates the frame down to it's children.

        Params:
            fast_axis The fast axis of the virtual detector plane
            slow_axis The slow axis of the virtual detector plane
            origin The origin vector to the virtual detector plane

        '''
        PanelBase.set_parent_frame(self, fast_axis, slow_axis, origin)
        for child in self:
            child.set_parent_frame(
                self.get_fast_axis(),
                self.get_slow_axis(),
                self.get_origin())

    def set_local_frame(self, fast_axis, slow_axis, origin):
        ''' Set the local frame.

        Set's it's own local plane and then, after updating it's global
        frame, propagates the frame down to it's children.

        Params:
            fast_axis The fast axis of the virtual detector plane
            slow_axis The slow axis of the virtual detector plane
            origin The origin vector to the virtual detector plane

        '''
        PanelBase.set_local_frame(self, fast_axis, slow_axis, origin)
        for child in self:
            child.set_parent_frame(
                self.get_fast_axis(),
                self.get_slow_axis(),
                self.get_origin())

    def add_group(self):
        ''' Add a new group to this group.

        Returns:
            A new panel group

        '''
        group = PanelGroup()
        group.set_parent_frame(
            self.get_fast_axis(),
            self.get_slow_axis(),
            self.get_origin())
        self._children.append(group)
        return group

    def add_panel(self):
        ''' Add a new panel to this group.

        Returns:
            A new panel

        '''
        panel = Panel2()
        panel.set_parent_frame(
            self.get_fast_axis(),
            self.get_slow_axis(),
            self.get_origin())
        self._children.append(panel)
        return panel

    def __getitem__(self, index):
        ''' Get the child at the given index.

        Params:
            index The index of the child to get.

        '''
        return self._children[index]

    def remove(self, item):
        ''' Remove a child from the tree. '''
        return self._children.remove(item)

    def index(self, item):
        ''' Get the index of a child. '''
        return self._children.index(item)

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
        if PanelBase.__eq__(self, other):
            if len(self) != len(other):
                return False
            return all(a == b for a, b in zip(self, other))
        else:
            return False

    def __ne__(self, other):
        ''' Check that this is not equal to another group. '''
        return not self.__eq__(other)


class Detector2(PanelGroup):
    ''' The top level Detector model.

    This class is derived from the panel group class but provides some
    additional convenience methods for navigating the panel hierarchy.

    '''

    def panels(self):
        ''' Extract a flex.panel array of all the panels depth-first '''
        from dxtbx.array_family import flex
        return flex.panel([p for p in self.iter_panels()])

    def iter_panels(self):
        ''' Iterate through just the panels depth-first. '''
        for obj in self.iter_preorder():
            if isinstance(obj, Panel2):
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
