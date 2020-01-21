from __future__ import division
import math

class Reparams(object):
    def __init__(self, reparamed_=False, offset_=0.):
        self.reparamed = reparamed_
        self.offset = offset_

    def off(self):
        self.reparamed = False
        self.offset = 0.

    def on(self, f):
        self.reparamed = True
        self.offset = f

    def param(self, par):
        return math.log(par + self.offset) if self.reparamed else par

    def unparam(self, par):
        return math.exp(par) - self.offset if self.reparamed else par

    def gradient(self, par, grad):
        return grad*(par + self.offset) if self.reparamed else par
