from __future__ import division

class Bounds(object):
  def __init__(self):
    self.lower_bounded_ = False
    self.upper_bounded_ = False
    self.lower_limit_ = 0.
    self.upper_limit_ = 0.

  def lower_off(self):
    self.lower_bounded_ = False

  def upper_off(self):
    self.upper_bounded_ = False

  def lower_on(self, l):
    self.lower_bounded_ = True
    self.lower_limit_ = l

  def upper_on(self, u):
    self.upper_bounded_ = True
    self.upper_limit_ = u

  def on(self, l, u):
    self.lower_bounded_ = self.upper_bounded_ = True
    self.lower_limit_ = l
    self.upper_limit_ = u

  def off(self):
    self.lower_bounded_ = self.upper_bounded_ = False

  def lower_bounded(self):
    return self.lower_bounded_

  def upper_bounded(self):
    return self.upper_bounded_

  def lower_limit(self):
    return self.lower_limit_

  def upper_limit(self):
    return self.upper_limit_
