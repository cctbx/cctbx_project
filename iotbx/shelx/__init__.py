from iotbx.shelx.errors import *
from iotbx.shelx.lexer import *
from iotbx.shelx.parsers import *

import boost.python
ext = boost.python.import_ext("iotbx_shelx_ext")

def cctbx_xray_structure_from(cls, file=None, filename=None,
                              set_grad_flags=True):
  from iotbx import builders
  builder = builders.crystal_structure_builder(set_grad_flags=set_grad_flags)
  stream = command_stream(file=file, filename=filename)
  stream = crystal_symmetry_parser(stream, builder)
  stream = atom_parser(stream.filtered_commands(), builder)
  stream.parse()
  return builder.structure
