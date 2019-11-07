"""
PDB web service API
"""

from __future__ import absolute_import, division, print_function

from iotbx.pdb.download import openurl, NotFound, identity_encoding


class RESTService(object):
  """
  Unified interface for services offering GET and POST requests
  """

  def __init__(self, url, get, post):

    self.url = url
    self.get = get
    self.post = post


  def single(self, identifier):

    return self.get( url = self.url, identifier = identifier )


  def multiple(self, identifiers):

    return self.post( url = self.url, identifiers = identifiers )


class FTPService(object):
  """
  Interface to access PDB files
  """

  def __init__(
    self,
    url,
    namer,
    postprocessing = identity_encoding(),
    opener = openurl,
    ):

    self.url = url
    self.namer = namer
    self.opener = opener
    self.postprocessing = postprocessing


  def single(self, identifier):

    stream = self.opener(
      url = ( self.url + self.namer( identifier = identifier ) ).encode(),
      )
    return self.postprocessing.process( stream = stream )


  def multiple(self, identifiers):

    for ident in identifiers:
      try:
        stream = self.single( identifier = ident )

      except NotFound:
        yield None

      else:
        yield stream
