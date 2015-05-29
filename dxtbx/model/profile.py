#!/usr/bin/env python
#
# profile.py
#
# Copyright (C) 2015 Diamond Light Source, James Parkhurst
#
# This code is distributed under the BSD license, a copy of which is
# included in the root directory of this package.
#
from __future__ import division

from abc import ABCMeta, abstractmethod

class Profile(object):
  '''
  A basic profile model class. This assumes no interface other than
  serialization methods. Classes inheriting from this class can
  then be serialized via a simple registry

  '''

  __metaclass__ = ABCMeta

  @abstractmethod
  def to_dict(self):
    '''
    Convert the profile model to dictionary

    '''
    pass

  @abstractmethod
  def from_dict(cls, obj):
    '''
    Convert a dictionary to the profile model

    '''
    pass


class ProfileFactory(object):
  '''
  A factory to create a profile model

  '''

  @staticmethod
  def classes():
    '''
    Iterate through the subclasses

    '''
    stack = list(Profile.__subclasses__())
    while len(stack) > 0:
      cls = stack.pop()
      yield cls
      stack.extend(cls.__subclasses__())

  @staticmethod
  def find(name):
    '''
    Find a subclass with the given name

    '''
    for cls in ProfileFactory.classes():
      try:
        if cls.name == name:
          return cls
      except Exception:
        pass
    return None

  @staticmethod
  def from_dict(obj):
    '''
    Given a dictionary, convert to a profile model

    '''
    if obj is None:
      return None
    return ProfileFactory.find(obj['name']).from_dict(obj)
