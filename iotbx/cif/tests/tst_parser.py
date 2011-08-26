class accu_builder(object):

  __slots__ = ["accu"]

  def __init__(O):
    O.accu = []

  def add_data_block(self, data_block_heading):
    O.accu.append(("add_data_block", data_block_heading))

  def add_loop(self, header, data):
    O.accu.append(("add_loop", (header, data)))

  def add_data_item(self, key, value):
    O.accu.append(("add_data_item", (key, value)))

  def start_save_frame(self, save_frame_heading):
    O.accu.append(("start_save_frame", (save_frame_heading)))

  def end_save_frame(self):
    O.accu.append(("end_save_frame",))

def exercise(verbose):
  builder = accu_builder()
  import iotbx.cif
  reader = iotbx.cif.ext.fast_reader(input="", builder=builder, strict=True)
  assert reader.lexer_errors().size() == 0
  assert reader.parser_errors().size() == 0
  assert len(builder.accu) == 0
  for i in xrange(256):
    builder = accu_builder()
    reader = iotbx.cif.ext.fast_reader(
      input=chr(i), builder=builder, strict=True)
    if (verbose):
      print reader.lexer_errors().size(), \
            reader.parser_errors().size(), \
            len(builder.accu)

def run(args):
  assert args in [[], ["--forever"]]
  forever = (len(args) != 0)
  while True:
    exercise(verbose=not forever)
    if (not forever):
      break

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
