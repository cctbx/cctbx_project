from __future__ import division

from simtbx.diffBragg.refiners.crystal_systems import CrystalSystemManager
from scitbx.matrix import sqr
import numpy as np

class TetragonalManager(CrystalSystemManager):

    def __init__(self, a=55, c=77):
        self.variables = [a, c]

    @property
    def variables(self):
        return self._variables

    @variables.setter
    def variables(self, val):
        self._variables = val

    @property
    def derivative_matrices(self):
        return [self._dB_da_real, self._dB_dc_real]

    @property
    def second_derivative_matrices(self):
        return [self._d2B_da2_real, self._d2B_dc2_real]

    @property
    def a(self):
        return self.variables[0]

    @property
    def b(self):
        return self.variables[0]

    @property
    def c(self):
        return self.variables[1]

    @property
    def al(self):
        return np.pi/2.

    @property
    def be(self):
        return np.pi/2

    @property
    def ga(self):
        return np.pi/2

    @property
    def variable_names(self):
        return [self._names[0], self._names[2]]

    @property
    def _dB_da_real(self):
        return sqr((1, 0, 0, 0, 1, 0, 0, 0, 0))

    @property
    def _dB_dc_real(self):
        return sqr((0, 0, 0, 0, 0, 0, 0, 0, 1))

    @property
    def _d2B_da2_real(self):
        return sqr((0, 0, 0, 0, 0, 0, 0, 0, 0))

    @property
    def _d2B_dc2_real(self):
        return sqr((0, 0, 0, 0, 0, 0, 0, 0, 0))
