"""
PDBe RESTful web services API
"""

from __future__ import division

from iotbx.pdb.download import openurl, NotFound

import json

def get_request(url, identifier):

  lcase = identifier.lower()

  stream = openurl( url = ( url + lcase ).encode() )
  result_for = json.load( stream )
  assert lcase in result_for
  return result_for[ lcase ]


def multiple_get_requests(url, identifiers):

  for ident in identifiers:
    try:
      result = get_request( url = url, identifier = ident )

    except NotFound:
      yield None

    else:
      yield result


def post_request(url, identifiers):

  if len( identifiers ) == 0:
    return ( i for i in () )

  try:
    stream = openurl( url = url, data = ( ",".join( identifiers ) ).encode() )

  except NotFound:
    return ( None for i in identifiers )

  result_for = json.load( stream )

  return ( result_for.get( ident.lower() ) for ident in identifiers )


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

  def __init__(self, url, namer):

    self.url = url
    self.namer = namer


  def single(self, identifier):

    stream = openurl(
      url = ( self.url + self.namer( identifier = identifier ) ).encode(),
      )
    return stream


  def multiple(self, identifiers):

    for ident in identifiers:
      try:
        stream = openurl(
          url = ( self.url + self.namer( identifier = ident ) ).encode(),
          )

      except NotFound:
        yield None

      else:
        yield stream


def identifier_to_pdb_entry_name(identifier):

  return "pdb%s.ent" % identifier.lower()


def identifier_to_cif_entry_name(identifier):

  return "%s.cif" % identifier.lower()


PDB_ENTRYFILE_PDB = FTPService(
  url = "http://www.ebi.ac.uk/pdbe/entry-files/",
  namer = identifier_to_pdb_entry_name,
  )

PDB_ENTRYFILE_CIF = FTPService(
  url = "http://www.ebi.ac.uk/pdbe/entry-files/",
  namer = identifier_to_cif_entry_name,
  )

PDBE_API_BASE = "http://wwwdev.ebi.ac.uk/pdbe/api/"

PDB_ENTRY_STATUS = RESTService(
  url = PDBE_API_BASE + "pdb/entry/status/",
  get = get_request,
  post = post_request,
  )

PDB_SIFTS_MAPPING_CATH = RESTService(
  url = PDBE_API_BASE + "mappings/cath/",
  get = get_request,
  post = multiple_get_requests,
  )

PDB_SIFTS_MAPPING_SCOP = RESTService(
  url = PDBE_API_BASE + "mappings/scop/",
  get = get_request,
  post = multiple_get_requests,
  )


class Redirections(object):
  """
  Use status queries
  """

  def __init__(self):

    self._currents = set()
    self._retracteds = set()
    self._replaced_by = {}


  def obsoleted(self, identifier):

    std = self._insert( identifier = identifier )
    return std in self._replaced_by


  def replacement_for(self, identifier):

    std = self._insert( identifier = identifier )
    return self._replaced_by[ std ]


  def retracted(self, identifier):

    std = self._insert( identifier = identifier )
    return std in self._retracteds


  def _insert(self, identifier):

    std = identifier.lower()

    if std in self._currents or std in self._retracteds or std in self._replaced_by:
      return std

    result = PDB_ENTRY_STATUS.single( identifier = std )
    status = result[0][ "status_code" ]

    if status == "REL":
      self._currents.add( std )

    elif status == "OBS":
      successor = result[0][ "superceded_by" ][0]

      if successor is None:
        self._retracteds.add( std )

      else:
        self._replaced_by[ std ] = successor

    return std

from libtbx.object_oriented_patterns import lazy_initialization
redirection = lazy_initialization( func = Redirections )
