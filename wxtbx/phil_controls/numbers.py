
# base control for 'ints' and 'floats' phil types.  see corresponding code
# in ../../libtbx/phil/__init__.py - these should probably be consolidated
# at some point

from wxtbx.phil_controls.text_base import ValidatedTextCtrl, TextCtrlValidator
from libtbx import Auto
import re

class NumbersCtrlBase (ValidatedTextCtrl) :
  def __init__ (self, *args, **kwds) :
    super(NumbersCtrlBase, self).__init__(*args, **kwds)
    self.size = None
    self.size_min = None
    self.size_max = None
    self.value_min = None
    self.value_max = None
    self.allow_none_elements = False
    self.allow_auto_elements = False

  def CheckType (self, value) :
    raise NotImplementedError()

  def InheritPhilType (self, phil_type) :
    """
    Adopts the value constraints from the PHIL converter object.  'size' is
    not actually an attribute of the type, but it is also extracted here.
    """
    for attr in ["size_min", "size_max", "value_min", "value_max",
                 "allow_none_elements", "allow_auto_elements"] :
      setattr(self, attr, getattr(phil_type, attr))
    if (self.size_min == self.size_max) and (self.size is None) :
      self.size = self.size_min

  def SetSize (self, size) :
    assert (size is None) or (isinstance(size, int))
    self.size = size

  def SetSizeMin (self, size_min) :
    assert (size_min is None) or (isinstance(size_min, int))
    self.size_min = size_min

  def SetSizeMax (self, size_max) :
    assert (size_max is None) or (isinstance(size_max, int))
    self.size_max = size_max

  def SetMin (self, min) :
    assert (min is None) or (self.CheckType(min))
    self.value_min = min

  def SetMax (self, max) :
    assert (max is None) or (self.CheckType(max))
    self.value_max = max

class NumbersValidator (TextCtrlValidator) :
  def ConvertValue (self, value) :
    raise NotImplementedError()

  def CheckFormat (self, value) :
    if ("," in value) or (";" in value) :
      value = re.sub(",", " ", re.sub(";", " ", value))
    numbers_list = []
    for field in value.split() :
      if (field == "None") :
        numbers_list.append(None)
      elif (field == "Auto") :
        numbers_list.append(Auto)
      else :
        numbers_list.append(self.ConvertValue(field))
    window = self.GetWindow()
    n_elems = len(numbers_list)
    if (window.size is not None) and (window.size != n_elems) :
      raise ValueError(("Wrong number of items - %d numbers required, "+
        "but %d are entered.") % (window.size, n_elems))
    if (window.size_min is not None) and (window.size_min > n_elems) :
      raise ValueError(("Wrong number of items - %d entered, but at least %d "+
        "are required.") % (n_elems, window.size_min))
    if (window.size_max is not None) and (window.size_max < n_elems) :
      raise ValueError(("Wrong number of items - %d entered, but the "+
        "maximum is %d.") % (n_elems, window.size_max))
    for x in numbers_list :
      if (x is Auto) and (not window.allow_auto_elements) :
        raise ValueError("Auto not allowed here.")
      elif (x is None) and (not window.allow_none_elements) :
        raise ValueError("None not allowed here.")
      elif (window.value_min is not None) and (x < window.value_min) :
        raise ValueError("Minimum permitted value is %g" % window.value_min)
      elif (window.value_max is not None) and (x > window.value_max) :
        raise ValueError("Maximum permitted value is %g" % window.value_max)
    return numbers_list
