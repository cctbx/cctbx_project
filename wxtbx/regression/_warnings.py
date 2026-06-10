from __future__ import absolute_import, division, print_function
import warnings

# Regex matching module names whose DeprecationWarning calls we want to promote
# to errors during the wxPython 4.2 migration tests.  Python's `warnings`
# module compiles this with re.compile(...) and consults it via .match(), so
# only the prefix matters: bare "wx" covers "wx", "wx.<anything>", "wxtbx",
# "wxtbx.<anything>", and any future "wx<x>" module too.  The remaining mmtbx
# leaf module in scope is listed explicitly.
WX_FILTER_RE = (
  r"wx|wxtbx|gltbx|mmtbx\.command_line\.map_box"
)

def install_wx_deprecation_filters():
  """Re-prepend our error filters so they sit at position 0 of warnings.filters.

  Transitively imported modules (notably phaser/deprecated_keywords_bpl.py)
  call warnings.filterwarnings("default", category=DeprecationWarning) with
  no module regex, which prepends a broad "default" filter that would
  otherwise mask our error filter.  Re-installing immediately before each
  import keeps our filter at the front of the list.
  """
  warnings.filterwarnings(
    "error", category=DeprecationWarning, module=WX_FILTER_RE)
  warnings.filterwarnings(
    "error", category=PendingDeprecationWarning, module=WX_FILTER_RE)
