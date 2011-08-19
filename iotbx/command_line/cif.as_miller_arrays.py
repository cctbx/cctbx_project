import iotbx.cif
from libtbx import easy_pickle, smart_open
from libtbx.str_utils import show_string
import sys, os
op = os.path

def run(args):
  for f in args:
    try:
      file_object = smart_open.for_reading(file_name=f)
      miller_arrays = iotbx.cif.reader(file_object=file_object).as_miller_arrays()
    except KeyboardInterrupt:
      raise
    except Exception, e:
      print "Error extracting miller arrays from file: %s:" % (
        show_string(f))
      print " ", str(e)
      continue
    for miller_array in miller_arrays:
      miller_array.show_comprehensive_summary()
      print
    r, _ = op.splitext(op.basename(f))
    easy_pickle.dump(file_name=r+'_miller_arrays.pickle', obj=miller_arrays)

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
