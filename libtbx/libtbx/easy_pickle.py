from libtbx.str_utils import show_string
import cPickle

def _open(file_name, mode):
  try: return open(file_name, mode)
  except IOError, e:
    raise IOError("Cannot open pickle file %s (%s)" % (
      show_string(file_name), str(e)))

def dump(file_name, obj):
  return cPickle.dump(obj, _open(file_name, "wb"), 1)

def load(file_name):
  return cPickle.load(_open(file_name, "rb"))

def dump_args(*args, **keyword_args):
  dump("args.pickle", (args, keyword_args))
