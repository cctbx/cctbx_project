import urlparse
import exceptions
class FormatError(exceptions.Exception): pass
import os

class server_info:

  def __init__(self):
    server_name = os.environ["SERVER_NAME"]
    server_port = os.environ["SERVER_PORT"]
    script_name = os.environ["SCRIPT_NAME"]
    self._script = [
      'http',
      '%s:%s' % (server_name, server_port),
      script_name,
      '',
      '']
    self._base = [
      'http',
      '%s:%s' % (server_name, server_port),
      "/".join(script_name.split("/")[:-1]) + "/",
      '',
      '']

  def script(self, query=''):
    return urlparse.urlunsplit(
      self._script[:3] + [query] + self._script[4:])

  def base(self):
    return urlparse.urlunsplit(self._base)

  def file(self, target):
    return urlparse.urlunsplit(
      self._base[:2] + [self._base[2] + target] + self._base[3:])

def show_input_symbol(sgsymbol, convention, label="Input"):
  if (sgsymbol != ""):
    print label, "space group symbol:", sgsymbol
    print "Convention:",
    if   (convention == "A1983"):
      print "International Tables for Crystallography, Volume A 1983"
    elif (convention == "I1952"):
      print "International Tables for Crystallography, Volume I 1952"
    elif (convention == "Hall"):
      print "Hall symbol"
    else:
      print "Default"
    print

def interpret_skip_columns(skip_columns):
  result = int(skip_columns)
  if (result < 0):
    raise ValueError, "Negative number for columns to skip."
  return result

def interpret_coordinate_line(line, skip_columns):
  flds = line.split()
  if (len(flds) < skip_columns + 3): raise FormatError, line
  coordinates = [0,0,0]
  for i in xrange(3):
    try: coordinates[i] = float(flds[skip_columns + i])
    except: raise FormatError, line
  return " ".join(flds[:skip_columns]), coordinates
