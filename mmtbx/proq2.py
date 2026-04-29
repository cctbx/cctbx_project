"""Interface to ProQ2 server"""

from __future__ import absolute_import, division, print_function

from HTMLParser import HTMLParser


class ProQError(Exception):
  """
  Module exception
  """


class JobFolderParseError(ProQError):
  """
  Could not parse out job folder
  """


class JobFolderParser(HTMLParser):

  JOB_FOLDER_TEXT = "The job folder is located here:"

  def __init__(self):

    HTMLParser.__init__( self )

    self.job_folder = None
    self.capture_next_anchor = False


  def handle_starttag(self, tag, attrs):

    if tag == "a" and self.capture_next_anchor:
      self.job_folder = None

      for ( key, value ) in attrs:
        if key == "href":
          self.job_folder = value

      self.capture_next_anchor = False


  def handle_data(self, data):

    if data.find( self.JOB_FOLDER_TEXT ) != -1:
      self.capture_next_anchor = True


def parse_submission_page(stream):

  parser = JobFolderParser()

  for line in stream:
    parser.feed( line )

    if parser.job_folder is not None:
      return parser.job_folder

  raise JobFolderParseError("Could not find job folder text")


class Job(object):
  """
  Interface to ProQ2 server
  """

  SERVERURL = "http://duffman.it.liu.se/ProQ2/index.php"
  OUTPUTFILE = "input.pdb.orig.B"

  def __init__(self, pdbstr, name):

    from iotbx.pdb import download
    from six.moves import urllib

    data_for = {
      "pdb": pdbstr,
      "bfactorPDB": "Yes",
      "name": name,
      "nomail": "Yes",
      "do": "Submit",
      }

    stream = download.openurl(
      url = self.SERVERURL,
      data = urllib.parse.urlencode( list(data_for.items()) ),
      )
    self.job_folder = parse_submission_page( stream = stream )
    stream.close()


  def __call__(self):

    from iotbx.pdb import download
    return download.openurl(
      url = "%s/%s" % ( self.job_folder, self.OUTPUTFILE ),
      )


if __name__ == "__main__":
  import sys

  if len( sys.argv ) == 2:
    timeout = 600

  elif len( sys.argv ) == 3:
    try:
      timeout = float( sys.argv[2] )

    except ValueError as e:
      print("Cannot interpret as number: %s (%s)" % ( sys.argv[2], e ))

  else:
    print("Usage: %s PDBFILE <timeout = 600>" % sys.argv[0])
    sys.exit( 1 )

  pdbstr = open( sys.argv[1] ).read()
  sys.stdout.write( "Submitting job..." )
  sys.stdout.flush()

  try:
    myjob = Job( pdbstr = pdbstr, name = "phenix.proq2" )

  except JobFolderParseError as e:
    sys.stdout.write( "failed\n" )
    print("Unexpected response: cannot find job folder")
    sys.exit( 1 )

  sys.stdout.write( "done\n" )
  print("Job folder:", myjob.job_folder)

  sys.stdout.write( "Waiting for results" )
  sys.stdout.flush()

  from libtbx import progress
  from iotbx.pdb import download
  waiter = progress.complete_on_success( func = myjob, excspec = download.NotFound )

  try:
    progress.wait(
      condition = waiter,
      waittime = 2,
      timeout = timeout,
      callback = progress.streamprint( stream = sys.stdout, character = "." ),
      )

  except progress.TimeoutError as e:
    sys.stdout.write( "%s\n" % e )
    sys.exit( 1 )

  assert hasattr( waiter, "result" )
  result = waiter.result.read()
  waiter.result.close()

  sys.stdout.write( "done\n" )

  import os.path
  output = "proq2-%s" % os.path.basename( sys.argv[1] )
  print("Output:", output)

  with open( output, "w" ) as ofile:
    ofile.write( result )
    ofile.write( "\n" )

  print("Done!")
