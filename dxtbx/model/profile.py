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

class ProfileModelBaseIface(object):
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


class ProfileModelFactory(object):
  '''
  A factory to create a profile model

  '''

  _classes = { }

  @staticmethod
  def append(name, cls):
    '''
    Add a class to the registry

    '''
    ProfileModelFactory._classes[name] = cls

  @staticmethod
  def find(name):
    '''
    Find a subclass with the given name

    '''
    try:
      return ProfileModelFactory._classes[name]
    except Exception:
      return None

  @staticmethod
  def from_dict(obj):
    '''
    Given a dictionary, convert to a profile model

    '''
    from logging import warn
    if obj is None:
      return None
    Class = ProfileModelFactory.find(obj['__id__'])
    if Class == None:
      warn('No profile class %s registered' % obj['__id__'])
      return None
    return Class.from_dict(obj)
