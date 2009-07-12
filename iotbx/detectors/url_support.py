class potential_url_request:
  def __init__(self,text):
    self.text = text

  def is_url_request(self):
    from urlparse import urlparse, parse_qs, urlunparse
    try:
      self.parsed = urlparse(self.text)
    except:
      return False
    self.file = self.parsed.path.split("?")[0]

    if "?" in self.parsed.path: #i.e., for file scheme, the query string is not
                          # supported. It shows up in the path instead.
      self.qs = parse_qs(self.parsed.path.split("?")[1])
    return True
