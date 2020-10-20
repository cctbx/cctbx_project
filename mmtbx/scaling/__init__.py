
"""
Base module for Xtriage and related scaling functionality; this imports the
Boost.Python extensions into the local namespace, and provides core functions
for displaying the results of Xtriage.
"""

from __future__ import absolute_import, division, print_function
import cctbx.array_family.flex # import dependency
from libtbx.str_utils import make_sub_header, make_header, make_big_header
from libtbx import slots_getstate_setstate
from six.moves import cStringIO as StringIO
import sys

import boost_adaptbx.boost.python as bp
from six.moves import range
ext = bp.import_ext("mmtbx_scaling_ext")
from mmtbx_scaling_ext import *


class data_analysis(slots_getstate_setstate):
  def show(self, out=sys.stdout, prefix=""):
    raise NotImplementedError()

class xtriage_output(slots_getstate_setstate):
  """
  Base class for generic output wrappers.
  """
  # this is used to toggle behavior in some output methods
  gui_output = False

  def show_big_header(self, title):
    """
    Print a big header with the specified title.
    """
    raise NotImplementedError()

  def show_header(self, title):
    """
    Start a new section with the specified title.
    """
    raise NotImplementedError()

  def show_sub_header(self, title):
    """
    Start a sub-section with the specified title.
    """
    raise NotImplementedError()

  def show_text(self, text):
    """
    Show unformatted text.
    """
    raise NotImplementedError()

  def show(self, text):
    return self.show_text(text)

  def show_preformatted_text(self, text):
    """
    Show text with spaces and line breaks preserved; in some contexts this
    will be done using a monospaced font.
    """
    raise NotImplementedError()

  def show_lines(self, text):
    """
    Show partially formatted text, preserving paragraph breaks.
    """
    raise NotImplementedError()

  def show_paragraph_header(self, text):
    """
    Show a header/title for a paragraph or small block of text.
    """
    raise NotImplementedError()

  def show_table(self, table, indent=0, plot_button=None,
      equal_widths=True):
    """
    Display a formatted table.
    """
    raise NotImplementedError()

  def show_plot(self, table):
    """
    Display a plot, if supported by the given output class.
    """
    raise NotImplementedError()

  def show_plots_row(self, tables):
    """
    Display a series of plots in a single row.  Only used for the Phenix GUI.
    """
    raise NotImplementedError()

  def show_text_columns(self, rows, indent=0):
    """
    Display a set of left-justified text columns.  The number of columns is
    arbitrary but this will usually be key:value pairs.
    """
    raise NotImplementedError()

  def newline(self):
    """
    Print a newline and nothing else.
    """
    raise NotImplementedError()

  def write(self, text):
    """
    Support for generic filehandle methods.
    """
    self.show(text)

  def flush(self):
    """
    Support for generic filehandle methods.
    """
    pass

  def warn(self, text):
    """
    Display a warning message.
    """
    raise NotImplementedError()

class printed_output(xtriage_output):
  """
  Output class for displaying raw text with minimal formatting.
  """
  __slots__ = ["out"]
  def __init__(self, out):
    assert hasattr(out, "write") and hasattr(out, "flush")
    self.out = out
    self._warnings = []

  def show_big_header(self, text):
    make_big_header(text, out=self.out)

  def show_header(self, text):
    make_header(text, out=self.out)

  def show_sub_header(self, title):
    out_tmp = StringIO()
    make_sub_header(title, out=out_tmp)
    for line in out_tmp.getvalue().splitlines():
      self.out.write("%s\n" % line.rstrip())

  def show_text(self, text):
    print(text, file=self.out)

  def show_paragraph_header(self, text):
    print(text, file=self.out) #+ ":"

  def show_preformatted_text(self, text):
    print(text, file=self.out)

  def show_lines(self, text):
    print(text, file=self.out)

  def show_table(self, table, indent=2, plot_button=None, equal_widths=True):
    print(table.format(indent=indent, equal_widths=equal_widths), file=self.out)

  def show_plot(self, table):
    pass

  def show_plots_row(self, tables):
    pass

  def show_text_columns(self, rows, indent=0):
    prefix = " "*indent
    n_cols = len(rows[0])
    col_sizes = [ max([ len(row[i]) for row in rows ]) for i in range(n_cols) ]
    for row in rows :
      assert len(row) == n_cols
      formats = prefix+" ".join([ "%%%ds" % x for x in col_sizes ])
      print(formats % tuple(row), file=self.out)

  def newline(self):
    print("", file=self.out)

  def write(self, text):
    self.out.write(text)

  def warn(self, text):
    self._warnings.append(text)
    out_tmp = StringIO()
    make_sub_header("WARNING", out=out_tmp, sep='*')
    for line in out_tmp.getvalue().splitlines():
      self.out.write("%s\n" % line.rstrip())
    self.out.write(text)

class loggraph_output(xtriage_output):
  """
  Output class for displaying 'loggraph' format (from ccp4i) as plain text.
  """
  gui_output = True
  def __init__(self, out):
    assert hasattr(out, "write") and hasattr(out, "flush")
    self.out = out

  def show_big_header(self, text) : pass
  def show_header(self, text) : pass
  def show_sub_header(self, title) : pass
  def show_text(self, text) : pass
  def show_paragraph_header(self, text) : pass
  def show_preformatted_text(self, text) : pass
  def show_lines(self, text) : pass
  def show_table(self, *args, **kwds) : pass
  def show_text_columns(self, *args, **kwds) : pass
  def newline(self) : pass
  def write(self, text) : pass
  def warn(self, text) : pass

  def show_plot(self, table):
    print("", file=self.out)
    print(table.format_loggraph(), file=self.out)

  def show_plots_row(self, tables):
    for table in tables :
      self.show_plot(table)

class xtriage_analysis(object):
  """
  Base class for analyses performed by Xtriage.  This does not impose any
  restrictions on content or functionality, but simply provides a show()
  method suitable for either filehandle-like objects or objects derived from
  the xtriage_output class.  Child classes should implement _show_impl.
  """
  def show(self, out=None):
    if out is None:
      out=sys.stdout
    if (not isinstance(out, xtriage_output)):
      out = printed_output(out)
    self._show_impl(out=out)
    return self

  def _show_impl(self, out):
    raise NotImplementedError()

  def summarize_issues(self):
    return []
