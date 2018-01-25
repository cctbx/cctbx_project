from __future__ import absolute_import, division

from abc import ABCMeta, abstractmethod

class ScalingModelBaseIface(object):
  '''
  A basic scaling model class. This assumes no interface other than
  serialization methods. Classes inheriting from this class can
  then be serialized via a simple registry

  '''

  __metaclass__ = ABCMeta

  @abstractmethod
  def to_dict(self):
    '''
    Convert the scaling model to dictionary

    '''
    pass

  @abstractmethod
  def from_dict(cls, obj):
    '''
    Convert a dictionary to the scaling model

    '''
    pass


class ScalingModelFactory(object):
  '''
  A factory to create a scaling model

  '''

  _classes = {}

  @staticmethod
  def append(name, cls):
    '''
    Add a class to the registry

    '''
    ScalingModelFactory._classes[name] = cls

  @staticmethod
  def find(name):
    '''
    Find a subclass with the given name

    '''
    try:
      return ScalingModelFactory._classes[name]
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
    Class = ScalingModelFactory.find(obj['__id__'])
    if Class is None:
      warn('No profile class %s registered' % obj['__id__'])
      return None
    return Class.from_dict(obj)
