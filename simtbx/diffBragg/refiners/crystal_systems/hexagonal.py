from simtbx.diffBragg.refiners.crystal_systems import CrystalSystemManager
import numpy as np
from scitbx.matrix import sqr


class HexagonalManager(CrystalSystemManager):

    def __init__(self, a=77, c=264):
        self.variables = [a, c]

    @property
    def variables(self):
        return self._variables

    @variables.setter
    def variables(self, val):
        if len(val) != 2:
            raise ValueError("Hexagonal crystal system has two variables: %s" % ", ".join(self.variable_names))
        self._variables = val

    @property
    def derivative_matrices(self):
        return [self._dB_da_real,
                self._dB_dc_real]

    @property
    def second_derivative_matrices(self):
        return [self._d2B_da2_real,
                self._d2B_dc2_real]

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
        return np.pi/2

    @property
    def be(self):
        return np.pi/2

    @property
    def ga(self):
        return 2*np.pi/3

    @property
    def variable_names(self):
        return [self._names[0],
                self._names[2]]

    @property
    def _dB_da_real(self):
        return sqr((1, -.5, 0,
                    0, np.sqrt(3)/2., 0,
                    0, 0, 0))

    @property
    def _d2B_da2_real(self):
        return sqr((0, 0, 0,
                    0, 0, 0,
                    0, 0, 0))

    @property
    def _dB_dc_real(self):
        return sqr((0, 0, 0,
                    0, 0, 0,
                    0, 0, 1))

    @property
    def _d2B_dc2_real(self):
        return sqr((0, 0, 0,
                    0, 0, 0,
                    0, 0, 0))
