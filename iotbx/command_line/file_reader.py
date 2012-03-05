
from iotbx import file_reader
import libtbx.phil
from libtbx.utils import Sorry, Usage
import os
import sys

master_phil = libtbx.phil.parse("""
file_name = None
  .type = path
force_type = None
  .type = str
""")

def run (args) :
  if (len(args) == 0) :
    raise Usage("""
iotbx.file_reader filename [force_type=xxx]

(where 'xxx' is one of these keywords:
  %s)
""" % ",".join(iotbx.file_reader.standard_file_types))
  user_phil = []
  for arg in args :
    if (os.path.isfile(arg)) :
      user_phil.append(libtbx.phil.parse("""file_name='%s'""" % arg))
    else :
      try :
        user_phil.append(libtbx.phil.parse(arg))
      except RuntimeError, e :
        print e
        print "Unrecognized argument '%s'" % arg
  params = master_phil.fetch(sources=user_phil).extract()
  if (params.file_name is None) :
    raise Sorry("No file specified.")
  f = file_reader.any_file(params.file_name, force_type=params.force_type)
  print "File type: %s (%s)" % (f.file_type,
    file_reader.standard_file_descriptions.get(f.file_type, "unknown"))

if (__name__ == "__main__") :
  run(sys.argv[1:])
