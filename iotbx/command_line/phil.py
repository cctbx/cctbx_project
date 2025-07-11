"""Analyze a PHIL file"""
from __future__ import absolute_import, division, print_function
import iotbx.phil
import libtbx.command_line.phil
import sys

if (__name__ == "__main__"):
  libtbx.command_line.phil.run(
    args=sys.argv[1:],
    command_name="iotbx.phil",
    converter_registry=iotbx.phil.default_converter_registry)
