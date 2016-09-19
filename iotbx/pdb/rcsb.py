from __future__ import division

from iotbx.pdb.download import openurl, NotFound, gzip_encoding


def parse_pdb_redirections(data):

  import re
  regex = re.compile(
    r"""
    ^ OBSLTE \s+
    [\w-]+ [ ]+
    ( \w{4} )
    (?: [ ]+ ( \w{4} ) )?
    \s* $
    """,
    re.VERBOSE | re.MULTILINE
    )

  replacement_for = {}

  for hit in regex.finditer( data ):
    ( old, new ) = hit.groups()
    replacement_for[ old.upper() ] = new.upper() if new else None

  return replacement_for


def end_of_chain(identifier, replacement_for):

  while True:
    if identifier is None or identifier not in replacement_for:
      break

    else:
      identifier = replacement_for[ identifier ]

  return identifier


class Redirections(object):
  """
  Downloads a list of redirections from the PDB
  """

  def __init__(self):

    stream = openurl(
      url = "ftp://ftp.wwpdb.org/pub/pdb/data/status/obsolete.dat"
      )

    chain_redirection_for = parse_pdb_redirections( data = stream.read() )

    self.redirection_for = {}

    for key in chain_redirection_for:
      self.redirection_for[ key ] = end_of_chain(
        identifier = key,
        replacement_for = chain_redirection_for,
        )


  def obsoleted(self, identifier):

    std = identifier.upper()
    return ( std in self.redirection_for and self.redirection_for[ std ] is not None )


  def replacement_for(self, identifier):

    std = identifier.upper()
    return self.redirection_for[ std ]


  def retracted(self, identifier):

    std = identifier.upper()
    return ( std in self.redirection_for and self.redirection_for[ std ] is None )


  def seed(self, identifiers):

    pass


from libtbx.object_oriented_patterns import lazy_initialization
redirection = lazy_initialization( func = Redirections )


class FTPService(object):
  """
  Interface to access PDB files
  """

  def __init__(self, url, namer, postprocessing, opener = openurl):

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


def identifier_to_pdb_entry_name(identifier):

  return "%s.pdb.gz" % identifier.lower()


def identifier_to_cif_entry_name(identifier):

  return "%s.cif.gz" % identifier.lower()


PDB_ENTRYFILE_PDB = FTPService(
  url = "http://files.rcsb.org/download/",
  namer = identifier_to_pdb_entry_name,
  postprocessing = gzip_encoding(),
  )

PDB_ENTRYFILE_CIF = FTPService(
  url = "http://www.ebi.ac.uk/pdbe/entry-files/",
  namer = identifier_to_cif_entry_name,
  postprocessing = gzip_encoding(),
  )
