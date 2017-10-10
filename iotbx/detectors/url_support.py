from __future__ import division
from future import standard_library
standard_library.install_aliases()
from builtins import object
class potential_url_request(object):
  def __init__(self,text):
    self.text = text

  def is_url_request(self):
    #backward compatibility with Python 2.5
    try: from urllib.parse import parse_qs
    except Exception: from cgi import parse_qs

    from urllib.parse import urlparse
    try:
      self.parsed = urlparse(self.text)
    except Exception:
      return False

    if self.parsed.scheme in [None, ""]: return False
    self.file = self.parsed.path.split("?")[0]

    if "?" in self.parsed.path: #i.e., for file scheme, the query string is not
                          # supported. It shows up in the path instead.
      self.qs = parse_qs(self.parsed.path.split("?")[1])
    return True
