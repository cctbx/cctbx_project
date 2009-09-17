from iotbx import shelx

def extract_from(file_name=None, file=None, max_characters=100000):
  # XXX leave max_characters keyword for back-compatibility 17-09-09
  assert [file_name, file].count(None) == 1
  stream = shelx.command_stream(file=file, filename=file_name)
  l = shelx.crystal_symmetry_parser(stream,
                                    builder=shelx.crystal_symmetry_builder())
  l.parse()
  return l.builder.crystal_symmetry
