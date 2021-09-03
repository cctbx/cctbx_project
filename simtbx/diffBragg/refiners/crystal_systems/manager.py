from __future__ import division

from abc import ABCMeta #, abstractproperty
import numpy as np
from scitbx.matrix import sqr

import abc
try:
    ABC = abc.ABC
    abstractproperty = lambda f: property(abc.abstractmethod(f))
except AttributeError:  # Python 2.7, abc exists, but not ABC
    ABC = abc.ABCMeta("ABC", (object,), {"__slots__": ()})
    from abc import abstractproperty

class CrystalSystemManager(object):

    __metaclass__ = ABCMeta

    @abstractproperty
    def variables(self):
        pass

    @abstractproperty
    def derivative_matrices(self):
        pass

    #@abstractproperty
    #def second_derivative_matrices(self):
    #    pass

    @abstractproperty
    def a(self):
        """unit cell a parameter"""
        pass

    @abstractproperty
    def b(self):
        """unit cell b parameter"""
        pass

    @abstractproperty
    def c(self):
        """unit cell c parameter"""
        pass

    @abstractproperty
    def al(self):
        """unit cell alpha angle"""
        pass

    @abstractproperty
    def be(self):
        """unit cell beta angle"""
        pass

    @abstractproperty
    def ga(self):
        """unit cell gamma angle"""
        pass

    @abstractproperty
    def variables(self):
        """the unit cell variables"""
        return []

    @property
    def cal(self):
        """cosine of alpha"""
        return np.cos(self.al)

    @property
    def cbe(self):
        """cosine of beta"""
        return np.cos(self.be)

    @property
    def cga(self):
        """cosine of gamma"""
        return np.cos(self.ga)

    @property
    def sal(self):
        """sine of alpha"""
        return np.sin(self.al)

    @property
    def sbe(self):
        """sine of beta"""
        return np.sin(self.be)

    @property
    def sga(self):
        """sine of gamma"""
        return np.sin(self.ga)

    @property
    def V(self):
        """unit cell volume in cubic Angstrom"""
        return self.a*self.b*self.c*np.sqrt(1-self.cal**2 - self.cbe**2-self.cga**2+2*self.cal*self.cbe*self.cga)

    @property
    def B_realspace(self):
        """
        real space B-matrix in upper triangular form
        such that the columns are a_real, b_real, c_real
        aligned according to the PDB convention
        """
        return sqr(
            (self.a, self.b*self.cga, self.c*self.cbe,
             0,      self.b*self.sga, self.c*(self.cal-self.cbe*self.cga)/self.sga,
             0,      0,               self.V/(self.a*self.b*self.sga)))

    @property
    def B_recipspace(self):
        """
        reciprocal space B-matrix in lower-triangular form
        such that it corresponds to the dxtbx Crystal get_B() method
        """
        return self.B_realspace.inverse().transpose()

    @property
    def _names(self):
        return ["a_Ang", "b_Ang", "c_Ang", "alpha_rad", "beta_rad", "gamma_rad"]

    @abstractproperty
    def variable_names(self):
        return []

    @property
    def named_variables(self):
        named_vars = {}
        for name, val in zip(self.variable_names, self.variables):
            named_vars[name] = val
        return named_vars

    @property
    def unit_cell_parameters(self):
        """
        unit cell parameter 6-tuple in Angstrom/degrees format
        returns: a,b,c, alpha, beta, gamma
        """
        return self.a, self.b, self.c, self.al*180./np.pi, self.be*180./np.pi, self.ga*180./np.pi
