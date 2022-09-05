from __future__ import absolute_import, division, print_function
from iotbx import shelx, builders

def extract_from(file_name=None, file=None, max_characters=100000,
                 close_file=True):
  # XXX backward compatibility 2009-09-17 - keep max_characters keyword
  assert [file_name, file].count(None) == 1
  stream = shelx.command_stream(file=file, filename=file_name,
                                close_when_deleted=close_file)
  parser = shelx.crystal_symmetry_parser(
    stream, builder=builders.crystal_symmetry_builder())
  parser.parse()
  return parser.builder.crystal_symmetry
