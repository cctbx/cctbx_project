
from __future__ import division

from dxtbx.model import Profile, ProfileFactory

class NullProfile(Profile):

  name = 'null'

  def __init__(self, parameter):
    self.parameter = parameter

  def to_dict(self):
    return {
      'name' : self.name,
      'parameter' : self.parameter
    }

  @classmethod
  def from_dict(Class, obj):
    return Class(obj['parameter'])


class Test(object):

  def __init__(self):
    pass

  def run(self):
    profile1 = NullProfile(10)
    dictionary = profile1.to_dict()
    profile2 = ProfileFactory.from_dict(dictionary)
    assert(profile1.parameter == profile2.parameter)
    print 'OK'

if __name__ == '__main__':
  test = Test()
  test.run()
