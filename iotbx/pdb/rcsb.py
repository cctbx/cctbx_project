"""Tools for working with RCSB API"""
from __future__ import absolute_import, division, print_function

from iotbx.pdb.web_service_api import FTPService
from iotbx.pdb.download import openurl, urlopener, gzip_encoding


def parse_pdb_redirections(lineiter):

  import re
  regex = re.compile(
    r"""
    ^ OBSLTE \s+
    [\w-]+ [ ]+
    ( \w{4} )
    (?: [ ]+ ( \w{4} ) )*
    \s* $
    """,
    re.VERBOSE | re.MULTILINE
    )

  replacement_for = {}

  for line in lineiter:
    hit = regex.match( line )

    if hit is not None:
      groups = hit.groups()
      obsoleted = groups[0].upper()
      assert 2 <= len( groups )

      if groups[1] is None:
        replacement = None

      else:
        replacement = groups[1].upper()

      replacement_for[ obsoleted ] = replacement

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
      url = "https://files.wwpdb.org/pub/pdb/data/status/obsolete.dat",
      )

    chain_redirection_for = parse_pdb_redirections( lineiter = stream )

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

# Special opener that disables gzip encoding (to fetch files already compressed)
gz_openurl = urlopener( extras = [] )
gz_encoding = gzip_encoding()

def identifier_to_pdb_entry_name(identifier):

  return "%s.pdb.gz" % identifier.lower()


def identifier_to_cif_entry_name(identifier):

  return "%s.cif.gz" % identifier.lower()


PDB_ENTRYFILE_PDB = FTPService(
  url = "http://files.rcsb.org/download/",
  namer = identifier_to_pdb_entry_name,
  opener = gz_openurl,
  postprocessing = gz_encoding,
  )

PDB_ENTRYFILE_CIF = FTPService(
  url = "http://www.ebi.ac.uk/pdbe/entry-files/",
  namer = identifier_to_cif_entry_name,
  opener = gz_openurl,
  postprocessing = gz_encoding,
  )
