# This module is obsolete and only here for compatibility reasons.
# Use the Python subprocess module directly.
from __future__ import absolute_import, division, print_function
from subprocess import *
import warnings
warnings.warn("libtbx.subprocess_with_fixes is deprecated", DeprecationWarning)
