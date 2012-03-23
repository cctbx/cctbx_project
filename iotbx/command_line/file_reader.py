
from iotbx import file_reader
import libtbx.phil
from libtbx.utils import Sorry, Usage
import re
import os
import sys

master_phil = libtbx.phil.parse("""
file_reader
  .short_caption = File format verification
  .caption = This utility will check the format of any file supported as \
    input by Phenix.  Since most of the input methods are designed to be \
    tolerant of errors, it does not provide detailed feedback about \
    parsing issues.
  .style = auto_align box caption_img:icons/crystal_project/32x32/mimetypes/txt.png
{
  file_name = None
    .type = path
    .style = bold
  force_type = *None %s
    .type = choice(multi=False)
    .caption = Any_format %s
}
""" % (" ".join(file_reader.standard_file_types),
       " ".join([ re.sub(" ", "_", file_reader.standard_file_descriptions[ft])
                  for ft in file_reader.standard_file_types ])))

def run (args=(), params=None, out=sys.stdout) :
  if (len(args) == 0) and (params is None) :
    raise Usage("""
iotbx.file_reader filename [force_type=None]

(where force_type can be optionally set to one of these keywords:
  %s)
""" % ",".join(iotbx.file_reader.standard_file_types))
  user_phil = []
  for arg in args :
    if (os.path.isfile(arg)) :
      user_phil.append(libtbx.phil.parse(
        """file_reader.file_name='%s'""" % arg))
    else :
      if (arg.startswith("force_type")) :
        arg = "file_reader." + arg
      try :
        user_phil.append(libtbx.phil.parse(arg))
      except RuntimeError, e :
        print e
        print "Unrecognized argument '%s'" % arg
  params = master_phil.fetch(sources=user_phil).extract()
  validate_params(params)
  f = file_reader.any_file(
    file_name=params.file_reader.file_name,
    valid_types=file_reader.standard_file_types,
    force_type=params.file_reader.force_type)
  f.show_summary(out=out)

def validate_params (params) :
  if (params.file_reader.file_name is None) :
    raise Sorry("No file specified.")

if (__name__ == "__main__") :
  run(sys.argv[1:])
