"""Download a file from storage on a server"""
from __future__ import absolute_import, division, print_function

# Storage methods
def no_storage(stream):

  return stream


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


class persistent_storage(object):
  """
  Stores data in a file that is named as a function of the url
  """

  def __init__(self, namer, binary = True):

    self.namer = namer
    self.mode_suffix = "b" if binary else ""


  def __call__(self, stream):

    filename = self.namer( stream.url )

    import shutil

    with open( filename, "w" + self.mode_suffix ) as fp:
      shutil.copyfileobj( stream, fp )

    return open( filename, "r" + self.mode_suffix )


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

# Utility
class coupled_stream(object):
  """
  Couples associated streams so that they could be closed explicitly
  """

  def __init__(self, primary, auxiliaries):

    self.primary = primary
    self.auxiliaries = auxiliaries


  def close(self):

    self.primary.close()

    for stream in self.auxiliaries:
      stream.close()


  def read(self):

    return self.primary.read()


  def readline(self):

    return self.primary.readline()


  def readlines(self):

    return self.primary.readlines()


  def next(self):

    return next(self.primary)


  def __iter__(self):

    return self


  def __repr__(self):

    return "<coupled primary:%r>" % self.primary


# Encodings
class encoding(object):

  @classmethod
  def accept(cls, header):

    return header == cls.keyword


class identity_encoding(encoding):

  keyword = "identity"

  def __init__(self, storage = no_storage):

    self.storage = storage


  def process(self, stream):

    return self.storage( stream = stream )


  @classmethod
  def accept(cls, header):

    return not header or super( identity_encoding, cls ).accept( header = header )


class gzip_encoding(encoding):

  keyword = "gzip"

  def __init__(self, storage = memory_storage):

    self.storage = storage


  def process(self, stream):

    import gzip
    storage = self.storage( stream = stream )

    return coupled_stream(
      primary = gzip.GzipFile( fileobj = storage ),
      auxiliaries = [ storage ],
      )


class deflate_encoding_small(encoding):

  keyword = "deflate"

  def __init__(self, storage = no_storage):

    self.storage = storage


  def process(self, stream):

    storage = self.storage( stream = stream )

    import zlib
    data = zlib.decompress( storage.read() )
    storage.close()

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

  def __init__(self, identity = identity_encoding(), extras = [ gzip_encoding() ]):

    self.identity = identity
    self.encoding_for = dict( ( ec.keyword, ec ) for ec in extras )


  def __call__(self, url, data = None):

    from six.moves import urllib

    request = urllib.request.Request(
      url = url,
      data = data,
      headers = {
        "Accept-encoding": ", ".join( self.encoding_for ),
        },
      )

    try:
      stream = urllib.request.urlopen( request )

    except urllib.error.HTTPError as e:
      raise http_error_to_exception( error = e )

    used = stream.info().get( "Content-Encoding" )
    encoding = self.encoding_for.get( used, self.identity )

    if not encoding.accept( header = used ):
      raise UnexpectedResponse("Unknown encoding: %s" % used)

    return encoding.process( stream = stream )


openurl = urlopener()
