from iotbx import mtz
import  iotbx.cns.miller_array
from scitbx.python_utils import easy_pickle
import sys

def usage():
  return (  "usage: iotbx.reflection_file_writer.py"
          + " miller_arrays.pickle"
          + " [export_mapping_file]")

def run(args):
  remaining_args = []
  for arg in args:
    if (arg.startswith("--")):
      print usage()
      raise RuntimeError, "Unknown option: " + arg
    else:
      remaining_args.append(arg)
  args = remaining_args
  assert len(remaining_args) in (1,2)
  input_file_name = remaining_args[0]
  miller_arrays = easy_pickle.load(input_file_name)
  export_mapping_file_name = None
  if (len(remaining_args) == 1):
    print "export_mappings = {"
    for miller_array in miller_arrays:
      print 'r"%s": "",' % miller_array.info()
    print "}"
  else:
    export_mapping_file_name = remaining_args[1]
    exec(open(export_mapping_file_name).read())
    miller_array_info_dict = {}
    for miller_array in miller_arrays:
      miller_array_info_dict[miller_array.info()] = miller_array
    for miller_array_info,file_name in export_mappings.items():
      if (file_name == ""): continue
      miller_array = miller_array_info_dict[miller_array_info]
      if (file_name.lower().endswith(".mtz")):
        print "Writing MTZ file:", file_name
        miller_array.export_as_mtz(file_name, file_name[:-4])
      elif (file_name.lower().endswith(".cns")):
        print "Writing CNS reflection file:", file_name
        miller_array.export_as_cns_hkl(open(file_name, "w"), file_name)
      elif (file_name.lower().endswith(".pickle")):
        print "Writing pickle file:", file_name
        easy_pickle.dump(file_name, miller_array)
      else:
        raise RuntimeError, "Filename extension not recognized: "+file_name
