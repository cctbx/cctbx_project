from __future__ import print_function
from __future__ import division

class Compulsory:
  def target(self):
    """
    Return the value of the target function for the current
    macrocycle parameter values
    """
    raise NotImplementedError("Please implement function target")

  def get_macrocycle_parameters(self):
    """
    Return the current values of the parameters that are being refined this
    macrocycle. The dtmin user can store these in their refine object however
    they see fit.
    """
    raise NotImplementedError("Please implement function get_macrocycle_parameters")

  def set_macrocycle_parameters(self, newx):
    """
    Given a list or flex array of values for the parameters that are being refined
    this macrocycle, update the refine object's internal representation of these values.
    Doesn't return anything.
    """
    raise NotImplementedError("Please implement function set_macrocycle_parameters")

  def macrocycle_large_shifts(self):
    """
    Return a list of large shift values for the parameters that are being refined
    this macrocycle.
    """
    raise NotImplementedError("Please implement function macrocycle_large_shifts")

  def set_macrocycle_protocol(self, macrocycle_protocol):
    """
    This gives you the ability to control the subset of parameters refined in any
    macrocycle or any other parameter that varies with macrocycle.
    Changes to the macrocycle parameters are controlled through class member variables.
    The RefineBase class' member variable nmp must be set in this function.
    """
    raise NotImplementedError("Please implement function set_macrocycle_protocol")

  def macrocycle_parameter_names(self):
    """
    Return a list of strings naming each of the parameters being refined this macrocycle.
    """
    raise NotImplementedError("Please implement function macrocycle_parameter_names")
