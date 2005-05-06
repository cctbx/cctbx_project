import cPickle

def dump(file_name, obj):
  return cPickle.dump(obj, open(file_name, "wb"), 1)

def load(file_name):
  return cPickle.load(open(file_name, "rb"))
