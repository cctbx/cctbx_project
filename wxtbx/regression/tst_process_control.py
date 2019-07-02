
# TODO more comprehensive tests

from __future__ import absolute_import, division, print_function
from wxtbx.process_control import *

def exercise():
  from libtbx.test_utils import approx_equal, Exception_expected
  import sys
  def excepthook(*args, **kwds):
    pass
  sys._excepthook = excepthook
  app = wx.App(0)
  result = run_function_as_process_in_dialog(
    parent=None,
    thread_function=test_function_1,
    title="Test subprocess",
    message="Running test function as forked process...",
    callback=None)
  if (result is not None):
    assert approx_equal(result, 2635152.11891, eps=0.0001)
  result = run_function_as_detached_process_in_dialog(
    parent=None,
    thread_function=test_function_1,
    title="Test subprocess",
    message="Running test function as separate process...",
    tmp_dir=os.getcwd(),
    callback=None)
  if (result is not None):
    assert approx_equal(result, 2635152.11891, eps=0.0001)
  #result2 = run_function_as_thread_in_dialog(
  #  parent=None,
  #  thread_function=test_function_2,
  #  title="Test subprocess",
  #  message="Running test function in Python thread...")
  #assert approx_equal(result2, 21081692.7462, eps=0.0001)
  try :
    result = run_function_as_process_in_dialog(
      parent=None,
      thread_function=test_function_3,
      title="Test subprocess",
      message="Running test function as separate process...",
      callback=None)
  except RuntimeError :
    pass
  else :
    raise Exception_expected
    result = run_function_as_process_in_dialog(
    parent=None,
    thread_function=test_function_1,
    title="Test subprocess",
    message="Running test function as separate process...",
    callback=None)
  wx.Yield()
  print("OK")

if (__name__ == "__main__"):
  exercise()
