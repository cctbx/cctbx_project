from __future__ import division

# Storage methods
def memory_storage(stream):

  import cStringIO
  return cStringIO.StringIO( stream.read() )


class named_storage(object):
  """
  Stores data in a given file
  """

  def __init__(self, filename, binary = True):

    self.filename = filename
    self.mode_suffix = "b" if binary else ""


  def __call__(self, stream):

    import shutil

    with open( self.filename, "w" + self.mode_suffix ) as fp:
      shutil.copyfileobj( stream, fp )

    return open( self.filename, "r" + self.mode_suffix )


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


# Encodings; some of these are singletons
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

  def __init__(self, storage = memory_storage):

    self.storage = storage


  def process(self, stream):

    import gzip
    return gzip.GzipFile( fileobj = self.storage( stream = stream ) )


class deflate_encoding_small(encoding):

  keyword = "deflate"

  @staticmethod
  def process(stream):

    import zlib
    data = zlib.decompress( stream.read() )

    import cStringIO
    return cStringIO.StringIO( data )


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


class urlopener(object):
  """
  Configurable version of openurl function
  """

  def __init__(self, encodings = [ gzip_encoding(), deflate_encoding_small ]):

    self.encoding_for = dict( ( ec.keyword, ec ) for ec in encodings )


  def __call__(self, url, data = None):

    import urllib2

    request = urllib2.Request(
      url = url,
      data = data,
      headers = {
        "Accept-encoding": ", ".join( self.encoding_for ),
        },
      )

    try:
      stream = urllib2.urlopen( request )

    except urllib2.HTTPError, e:
      raise http_error_to_exception( error = e )

    used = stream.info().get( "Content-Encoding" )
    encoding = self.encoding_for.get( used, identity_encoding )

    if not encoding.accept( header = used ):
      raise UnexpectedResponse, "Unknown encoding: %s" % used

    return encoding.process( stream = stream )


openurl = urlopener()
