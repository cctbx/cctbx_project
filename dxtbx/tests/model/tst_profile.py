
from __future__ import absolute_import, division

from dxtbx.model import ProfileModelFactory
from dials.algorithms.profile_model.gaussian_rs import Model


class Test(object):

  def __init__(self):
    pass

  def run(self):
    profile1 = Model(None, 3, 0.1, 0.2, deg=True)
    dictionary = profile1.to_dict()
    profile2 = ProfileModelFactory.from_dict(dictionary)
    assert abs(profile1.sigma_b() - profile2.sigma_b()) < 1e-7
    assert abs(profile1.sigma_m() - profile2.sigma_m()) < 1e-7
    print 'OK'

if __name__ == '__main__':
  test = Test()
  test.run()
