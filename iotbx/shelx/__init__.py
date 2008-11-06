from iotbx.shelx.errors import *
from iotbx.shelx.lexer import *
from iotbx.shelx.parsers import *
from iotbx.shelx.builders import *

from cctbx import xray

from libtbx import object_oriented_patterns as oop

class _extend_xray_structure(oop.injector, xray.structure):

  def from_shelx(cls, file=None, filename=None,
                 set_grad_flags=True):
    builder = crystal_structure_builder(set_grad_flags=set_grad_flags)
    stream = command_stream(file=file, filename=filename)
    cs_parser = crystal_symmetry_parser(stream, builder)
    xs_parser = atom_parser(cs_parser.filtered_commands(), builder)
    xs_parser.parse()
    return builder.structure
  from_shelx = classmethod(from_shelx)
