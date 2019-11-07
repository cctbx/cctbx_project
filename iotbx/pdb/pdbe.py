"""
PDBe web services
"""

from __future__ import absolute_import, division, print_function

from iotbx.pdb.web_service_api import FTPService, RESTService
from iotbx.pdb.download import openurl, NotFound

import json


class configurable_get_request(object):

  def __init__(self, opener):

    self.opener = opener


  def __call__(self, url, identifier):

    lcase = identifier.lower()

    stream = self.opener( url = ( url + lcase ).encode() )
    result_for = json.load( stream )
    stream.close()
    assert lcase in result_for
    return result_for[ lcase ]


class configurable_multiget_request(object):

  def __init__(self, opener):

    self.get_request = configurable_get_request( opener = opener )


  def __call__(self, url, identifiers):

    for ident in identifiers:
      try:
        result = get_request( url = url, identifier = ident )

      except NotFound:
        yield None

      else:
        yield result


class configurable_post_request(object):

  def __init__(self, opener):

    self.opener = opener


  def __call__(self, url, identifiers):

    if len( identifiers ) == 0:
      return ( i for i in () )

    try:
      stream = self.opener( url = url, data = ( ",".join( identifiers ) ).encode() )

    except NotFound:
      return ( None for i in identifiers )

    result_for = json.load( stream )
    stream.close()

    return ( result_for.get( ident.lower() ) for ident in identifiers )


# Define instances for default use
get_request = configurable_get_request( opener = openurl )
multiget_request = configurable_multiget_request( opener = openurl )
post_request = configurable_post_request( opener = openurl )


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

PDBE_API_BASE = "http://www.ebi.ac.uk/pdbe/api/"

PDB_ENTRY_STATUS = RESTService(
  url = PDBE_API_BASE + "pdb/entry/status/",
  get = get_request,
  post = post_request,
  )

PDB_ENTRY_EXPERIMENT = RESTService(
  url = PDBE_API_BASE + "pdb/entry/experiment/",
  get = get_request,
  post = post_request,
  )

PDB_SIFTS_MAPPING_CATH = RESTService(
  url = PDBE_API_BASE + "mappings/cath/",
  get = get_request,
  post = multiget_request,
  )

PDB_SIFTS_MAPPING_SCOP = RESTService(
  url = PDBE_API_BASE + "mappings/scop/",
  get = get_request,
  post = multiget_request,
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

    std = self._check_and_fetch( identifier = identifier )
    return std in self._replaced_by


  def replacement_for(self, identifier):

    std = self._check_and_fetch( identifier = identifier )
    return self._replaced_by[ std ]


  def retracted(self, identifier):

    std = self._check_and_fetch( identifier = identifier )
    return std in self._retracteds


  def seed(self, identifiers):

    from six.moves import zip

    blocks = PDB_ENTRY_STATUS.multiple( identifiers = identifiers )

    for ( code, data ) in zip( identifiers, blocks ):
      if not data:
        continue

      self._insert( code = self._standardize( identifier = code ), data = data )


  def _standardize(self, identifier):

    return identifier.lower()


  def _check_and_fetch(self, identifier):

    std = self._standardize( identifier = identifier )

    if std in self._currents or std in self._retracteds or std in self._replaced_by:
      return std

    data = PDB_ENTRY_STATUS.single( identifier = std )
    self._insert( code = std, data = data )
    return std


  def _insert(self, code, data):

    status = data[0][ "status_code" ]

    if status == "REL":
      self._currents.add( code )

    elif status == "OBS":
      successor = data[0][ "superceded_by" ][0]

      if successor is None:
        self._retracteds.add( code )

      else:
        self._replaced_by[ code ] = successor


from libtbx.object_oriented_patterns import lazy_initialization
redirection = lazy_initialization( func = Redirections )
