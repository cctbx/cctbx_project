import math

class double_quadratic_barrier(object):

  def __init__(self, lower_bound, upper_bound, barrier_parameter=1):
    self.lower_bound, self.upper_bound = lower_bound, upper_bound
    self.barrier_parameter = barrier_parameter

  def functional(self, c):
    if c < self.lower_bound:
      return self.barrier_parameter * (c - self.lower_bound)**2
    elif c > self.upper_bound:
      return self.barrier_parameter * (c - self.upper_bound)**2
    else:
      return 0

  def gradient(self, c):
    if c < self.lower_bound:
      return 2 * self.barrier_parameter * (c - self.lower_bound)
    elif c > self.upper_bound:
      return 2 * self.barrier_parameter * (c - self.upper_bound)
    else:
      return 0
