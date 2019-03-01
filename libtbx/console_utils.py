
from __future__ import absolute_import, division, print_function
import sys

class colors_(object):
  def __init__(self):
    self.black = 30
    self.red = 31
    self.green = 32
    self.yellow = 33
    self.blue = 34
    self.magenta = 35
    self.cyan = 36
    self.white = 37
    self.endc = 0

  def get_escape_string(self, color, bold=False, bgcolor=None):
    values = []
    if (bold):
      values.append(1)
    if (bgcolor is not None):
      values.append(getattr(self, bgcolor) + 10)
    values.append(getattr(self, color))
    value_str = ";".join([ "%d" % num for num in values ])
    return '\033[%sm' % value_str

colors = colors_()

class console_out(object):
  def __init__(self, console_out=None, other_out=None):
    if (console_out is None):
      console_out = sys.stdout
    self.console_out = console_out
    self.other_out = other_out
    self._color = None

  def set_color(self, color, bold=False, bgcolor=None):
    self._color = color
    escape_string = colors.get_escape_string(color, bold, bgcolor)
    if (self.console_out.isatty()) and (sys.platform != "win32"):
      self.console_out.write('\033[0m')
      self.console_out.write(escape_string)

  def reset(self):
    self._color = None
    if (self.console_out.isatty()) and (sys.platform != "win32"):
      self.console_out.write('\033[0m')

  def write(self, *args, **kwds):
    self.console_out.write(*args, **kwds)
    if (self.other_out is not None):
      self.other_out.write(*args, **kwds)

  def flush(self):
    self.console_out.flush()

  def warn(self, message, endl=True):
    self.set_color("yellow", bold=True)
    self.write(message)
    if (endl) : self.write("\n")
    self.reset()

  def fail(self, message, endl=True):
    self.set_color("red", bold=True)
    self.write(message)
    if (endl) : self.write("\n")
    self.reset()

def exercise():
  out = console_out()
  print("Hello, world!", file=out)
  for color in ["red","blue","green","yellow","magenta","cyan"] :
    out.set_color(color)
    print("Hello, world!", file=out)
  out.warn("this is a warning message")
  out.fail("this is a failure message")

if (__name__ == "__main__"):
  exercise()
