from __future__ import division
import cctbx.array_family.flex # import dependency
from libtbx.str_utils import make_sub_header, make_header
from libtbx import slots_getstate_setstate
import sys

import boost.python
ext = boost.python.import_ext("mmtbx_scaling_ext")
from mmtbx_scaling_ext import *


class data_analysis (slots_getstate_setstate) :
  def show (self, out=sys.stdout, prefix="") :
    raise NotImplementedError()

class xtriage_output (slots_getstate_setstate) :
  """
  Base class for generic output wrappers.
  """
  def show_header (self, text) :
    raise NotImplementedError()

  def show_text (self, text) :
    raise NotImplementedError()

  def show (self, text) :
    return self.show_text(text)

  def show_preformatted_text (self, text) :
    raise NotImplementedError()

  def show_table (self, table, precision=6, indent=0) :
    raise NotImplementedError()

  def show_plot (self, plot) :
    raise NotImplementedError()

  def newline (self) :
    raise NotImplementedError()

  def write (self, text) :
    self.show(text)

  def flush (self) :
    pass

class printed_output (xtriage_output) :
  __slots__ = ["out"]
  def __init__ (self, out) :
    assert hasattr(out, "write") and hasattr(out, "flush")
    self.out = out

  def show_header (self, text) :
    make_header(text, out=self.out)

  def show_sub_header (self, text) :
    make_sub_header(text, out=self.out)

  def show_text (self, text) :
    print >> self.out, text

  def show_preformatted_text (self, text) :
    print >> self.out, text

  def show_table (self, table, precision=6, indent=0) :
    print >> self.out, table.format(indent=indent)

  def show_plot (self, plot) :
    pass

  def newline (self) :
    print >> self.out, ""

  def write (self, text) :
    self.out.write(text)

class xtriage_analysis (object) :
  def show (self, out=None) :
    if out is None:
      out=sys.stdout
    if (not isinstance(out, xtriage_output)) :
      out = printed_output(out)
    self._show_impl(out=out)
    return self

  def _show_impl (self, out) :
    raise NotImplementedError()
