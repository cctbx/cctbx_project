"""
Imports wildcard tools
"""
from __future__ import absolute_import, division, print_function
import boost_adaptbx.boost.python as bp
ext = bp.import_ext("iotbx_wildcard_ext")
from iotbx_wildcard_ext import *
