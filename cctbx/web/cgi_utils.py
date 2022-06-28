from __future__ import absolute_import, division, print_function
import os,sys,cgi

if sys.version_info.major >= 3:
    from urllib.parse import urlunsplit
if sys.version_info.major < 3:
    from urlparse import urlunsplit

class FormatError(Exception): pass

class empty: pass

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
    return urlunsplit(
      self._script[:3] + [query] + self._script[4:])

  def base(self):
    return urlunsplit(self._base)

  def file(self, target):
    return urlunsplit(
      self._base[:2] + [self._base[2] + target] + self._base[3:])

def inp_from_form(form, keys):
  inp = empty()
  for key in keys:
    if (key[0] in form):
      #v = cgi.escape(form[key[0]].value,True).strip()
      v = form[key[0]].value.strip()
      if (v == ""): v = key[1]
      inp.__dict__[key[0]] = v
    else:
      inp.__dict__[key[0]] = key[1]
  return inp

def coordinates_from_form(form, suffix=None):
  coordinates = []
  for key_root in ("coordinates", "coor_file"):
    if (suffix is None): key = key_root
    else:                key = key_root + "_" + suffix
    if (key in form):
      lines = cgi.escape(form[key].value,True).replace("\015", "\012").split("\012")
      for l in lines:
        s = l.strip()
        if (len(s) != 0): coordinates.append(s)
  return coordinates
