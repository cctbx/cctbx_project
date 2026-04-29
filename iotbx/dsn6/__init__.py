"""
Tools to work with dsn6 formatted files
"""
# TODO TESTS

from __future__ import absolute_import, division, print_function
import cctbx.array_family.flex # import dependency
import boost_adaptbx.boost.python as bp
ext = bp.import_ext("iotbx_dsn6_map_ext")
from iotbx_dsn6_map_ext import *
import iotbx_dsn6_map_ext as ext
