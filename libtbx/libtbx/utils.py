from __future__ import division
import time
import sys, os

def flat_list(nested_list):
  result = []
  if (hasattr(nested_list, "__len__")):
    for sub_list in nested_list:
      result.extend(flat_list(sub_list))
  else:
    result.append(nested_list)
  return result

class Keep: pass

class Sorry(Exception):

  # trick to get just "Sorry" instead of "libtbx.utils.Sorry"
  __module__ = "exceptions"

  def __init__(self, *args, **keyword_args):
    self.previous_tracebacklimit = getattr(sys, "tracebacklimit", None)
    sys.tracebacklimit = 0
    Exception.__init__(self, *args, **keyword_args)

  def __str__(self):
    self.reset_tracebacklimit()
    return Exception.__str__(self)

  def __del__(self):
    self.reset_tracebacklimit()

  def reset_tracebacklimit(self):
    if (hasattr(sys, "tracebacklimit")):
      if (self.previous_tracebacklimit is None):
        del sys.tracebacklimit
      else:
        sys.tracebacklimit = self.previous_tracebacklimit

def format_exception():
  type_, value = sys.exc_info()[:2]
  type_ = str(type_)
  if (type_.startswith("exceptions.")): type_ = type_[11:]
  return "%s: %s" % (type_, str(value))

def date_and_time():
  localtime = time.localtime()
  if (time.daylight and localtime[8] != 0):
    tzname = time.tzname[1]
    offs = -time.altzone
  else:
    tzname = time.tzname[0]
    offs = -time.timezone
  return time.strftime("Date %Y-%m-%d Time %H:%M:%S", localtime) \
       + " %s %+03d%02d" % (tzname, offs//3600, offs//60%60)

def format_cpu_times(show_micro_seconds_per_tick=True):
  t = os.times()
  result = "u+s,u,s: %.2f %.2f %.2f" % (t[0] + t[1], t[0], t[1])
  if (show_micro_seconds_per_tick):
    try: python_ticker = sys.gettickeraccumulation()
    except AttributeError: pass
    else:
      result += " micro-seconds/tick: %.3f" % ((t[0]+t[1])/python_ticker*1.e6)
  return result

class multi_out:

  def __init__(self):
    self.labels = []
    self.file_objects = []
    self.closed = False
    self.softspace = 0

  def register(self, label, file_object):
    assert not self.closed
    self.labels.append(label)
    self.file_objects.append(file_object)
    return self

  def replace_stringio(self, old_label, new_label, new_file_object):
    i = self.labels.index(old_label)
    old_file_object = self.file_objects[i]
    new_file_object.write(old_file_object.getvalue())
    old_file_object.close()
    self.labels[i] = new_label
    self.file_objects[i] = new_file_object

  def isatty(self):
    return False

  def close(self):
    for file_object in self.file_objects:
      file_object.close()
    self.closed = True

  def flush(self):
    for file_object in self.file_objects:
      flush = getattr(file_object, "flush", None)
      if (flush is not None): flush()

  def write(self, str):
    for file_object in self.file_objects:
      file_object.write(str)

  def writelines(self, sequence):
    for file_object in self.file_objects:
      file_object.writelines(sequence)

def exercise():
  assert flat_list(0) == [0]
  assert flat_list([1,2,3]) == [1,2,3]
  assert flat_list([1,[2,3,4],3]) == [1,2,3,4,3]
  assert flat_list([1,[2,3,4],[[3,4],[5,6]]]) == [1,2,3,4,3,4,5,6]
  try:
    raise RuntimeError("Trial")
  except:
    assert format_exception() == "RuntimeError: Trial"
  print "OK"

if (__name__ == "__main__"):
  exercise()
