def _open(file_name, mode):
  import os
  file_name = os.path.expanduser(file_name)
  from libtbx.str_utils import show_string
  from libtbx import smart_open
  try: return smart_open.file(file_name=file_name, mode=mode)
  except IOError, e:
    raise IOError("Cannot open pickle file %s (%s)" % (
      show_string(file_name), str(e)))

def dump(file_name, obj):
  import cPickle
  return cPickle.dump(obj, _open(file_name, "wb"), 1)

def dumps(obj):
  import cPickle
  return cPickle.dumps(obj, 1)

def load(file_name, faster_but_using_more_memory=True):
  import cPickle
  if (faster_but_using_more_memory):
    return cPickle.loads(_open(file_name, "rb").read())
  return cPickle.load(_open(file_name, "rb"))

def loads(string) :
  import cPickle
  return cPickle.loads(string)

def dump_args(*args, **keyword_args):
  dump("args.pickle", (args, keyword_args))
