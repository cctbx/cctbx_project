from simtbx.diffBragg.refiners.crystal_systems import CrystalSystemManager
import numpy as np
from scitbx.matrix import sqr

class MonoclinicManager(CrystalSystemManager):

    def __init__(self, a=55, b=65, c=77, beta=100*np.pi/180.):
        self.variables = [a, b, c, beta]

    @property
    def variables(self):
        return self._variables

    @variables.setter
    def variables(self, val):
        self._variables = val

    @property
    def derivative_matrices(self):
        return [self._dB_da_real, self._dB_db_real,
                self._dB_dc_real, self._dB_dbeta_real]

    @property
    def second_derivative_matrices(self):
        return [self._d2B_da2_real, self._d2B_db2_real,
                self._d2B_dc2_real, self._d2B_dbeta2_real]

    @property
    def a(self):
        return self.variables[0]

    @property
    def b(self):
        return self.variables[1]

    @property
    def c(self):
        return self.variables[2]

    @property
    def al(self):
        return np.pi/2

    @property
    def be(self):
        return self.variables[3]

    @property
    def ga(self):
        return np.pi/2

    @property
    def variable_names(self):
        return [self._names[0], self._names[1],
                self._names[2], self._names[4]]

    @property
    def _dB_da_real(self):
        return sqr((1, 0, 0,
                    0, 0, 0,
                    0, 0, 0))

    @property
    def _d2B_da2_real(self):
        return sqr((0, 0, 0,
                    0, 0, 0,
                    0, 0, 0))

    @property
    def _dB_db_real(self):
        return sqr((0, 0, 0,
                    0, 1, 0,
                    0, 0, 0))

    @property
    def _d2B_db2_real(self):
        return sqr((0, 0, 0,
                    0, 0, 0,
                    0, 0, 0))

    @property
    def _dB_dc_real(self):
        return sqr((0, 0, self.cbe,
                    0, 0, 0,
                    0, 0, np.sqrt(1-self.cbe**2)))

    @property
    def _d2B_dc2_real(self):
        return sqr((0, 0, 0,
                    0, 0, 0,
                    0, 0, 0))

    @property
    def _dB_dbeta_real(self):
        return sqr((0, 0, -self.c*self.sbe,
                    0, 0, 0,
                    0, 0, self.c*self.cbe*self.sbe/np.sqrt(1-self.cbe**2)))

    @property
    def _d2B_dbeta2_real(self):
        return sqr((0, 0, -self.c*self.cbe,
                    0, 0, 0,
                    0, 0, -self.c*np.sqrt(1-self.cbe**2)))

