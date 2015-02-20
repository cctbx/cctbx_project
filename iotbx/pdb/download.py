from __future__ import division

# Encodings; these are singletons
class encoding(object):

  @classmethod
  def accept(cls, header):

    return header == cls.keyword


class identity_encoding(encoding):

  keyword = "identity"

  @classmethod
  def accept(cls, header):

    return not header or super( identity_encoding, cls ).accept( header = header )

  @staticmethod
  def process(stream):

    return stream


class gzip_encoding(encoding):

  keyword = "gzip"

  @staticmethod
  def process(stream):

    import gzip
    return gzip.GzipFile( fileobj = stream )


class deflate_encoding_small(encoding):

  keyword = "deflate"

  @staticmethod
  def process(stream):

    import zlib
    data = zlib.decompress( stream.read() )

    import cStringIO
    return cStringIO.StringIO( data )


# Storages
def memory_storage(stream):

  import cStringIO
  return cStringIO.StringIO( stream.read() )


class named_storage(object):
  """
  Stores data in a given file
  """

  def __init__(self, filename, binary = True):

    self.filename = filename
    self.binary = binary


  def __call__(self, stream):

    import shutil

    with open( self.filename, 'wb' if self.binary else 'w' ) as fp:
      shutil.copyfileobj( stream, fp )

    return open( self.filename, 'rb' if self.binary else 'r' )


class temporary_storage(object):
  """
  Stores data in a temporary file
  """

  def __init__(self, binary):

    self.mode = "w+b" if binary else "w+"


  def __call__(self, stream):

    import tempfile
    mytemp = tempfile.TemporaryFile( self.mode )

    import shutil
    shutil.copyfileobj( stream, mytemp )

    mytemp.seek(0)

    return mytemp


#Exceptions
class DownloadException(Exception):
  """
  Base class
  """


class NotFound(DownloadException):
  """
  HTTP 404
  """


class NotAcceptable(DownloadException):
  """
  HTTP 406
  """


class ServerError(DownloadException):
  """
  HTTP 5XX
  """


class UnexpectedResponse(DownloadException):
  """
  Unexpected response from server
  """


def http_error_to_exception(error):

  if error.code == 404:
    return NotFound()

  elif error.code == 406:
    return NotAcceptable()

  elif 500 <= error.code < 600:
    return ServerError()

  else:
    return UnexpectedResponse( str( error ) )


def openurl(
  url,
  data = None,
  encodings = [ gzip_encoding, deflate_encoding_small ],
  storage = memory_storage,
  ):

  import urllib2

  request = urllib2.Request(
    url = url,
    data = data,
    headers = {
      "Accept-encoding": ", ".join( ec.keyword for ec in encodings ),
      },
    )

  try:
    urlstream = urllib2.urlopen( request )

  except urllib2.HTTPError, e:
    raise http_error_to_exception( error = e )

  used = urlstream.info().get( "Content-Encoding" )

  for ec in encodings:
    if ec.accept( header = used ):
      selected = ec
      break

  else:
    # identity encoding is always acceptable
    if identity_encoding.accept( header = used ):
      selected = identity_encoding

    else:
      raise UnexpectedResponse, "Unknown encoding: %s" % used

  filestream = storage( stream = urlstream )
  return selected.process( stream = filestream )
