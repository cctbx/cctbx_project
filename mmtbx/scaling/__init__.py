from __future__ import division
import cctbx.array_family.flex # import dependency
from libtbx import slots_getstate_setstate
import sys

import boost.python
ext = boost.python.import_ext("mmtbx_scaling_ext")
from mmtbx_scaling_ext import *


class data_analysis (slots_getstate_setstate) :
  def show (self, out=sys.stdout, prefix="") :
    raise NotImplementedError()

class xtriage_output (slots_getstate_setstate) :
  def show_text (self, text) :
    raise NotImplementedError()

  def show_preformatted_text (self, text) :
    raise NotImplementedError()

  def show_table (self, table, precision=6, indent=0) :
    raise NotImplementedError()

  def show_plot (self, plot) :
    raise NotImplementedError()

  def newline (self) :
    raise NotImplementedError()

class printed_output (xtriage_output) :
  __slots__ = ["out"]
  def __init__ (self, out) :
    assert hasattr(out, "write") and hasattr(out, "flush")
    self.out = out

  def show_text (self, text) :
    print >> self.out, text

  def show_preformatted_text (self, text) :
    print >> self.out, text

  def show_table (self, table, precision=6, indent=0) :
    print >> self.out, table.format(precision=precision, indent=indent)

  def show_plot (self, plot) :
    pass

  def newline (self) :
    print >> self.out, ""
