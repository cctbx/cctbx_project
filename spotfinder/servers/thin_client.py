import os
from spotfinder.servers.multipart_encoder import post_multipart

def get_spotfinder_url(filename,host,port):
  testurl = "%s:%d"%(host,port)
  selector = "/spotfinder"
  query_object = [
    ("filename",filename),
    ("bin",1),
  ]

  Response = post_multipart(host=testurl, selector=selector,
    fields = query_object, files = [])

  print Response.read()

def do_main(filepath, host, port):
  absfile = os.path.abspath(filepath)
  get_spotfinder_url(absfile,host,port)


if __name__=="__main__":
  "Client is intended to be used with [mp_]spotfinder_server_read_file.py"
  import sys
  try:
    filepath, host, port = sys.argv[1:4]
    port = int(port)
  except:
    print """
Usage:
libtbx.python thin_client.py <filepath> <host> <port>
Three mandatory arguments:
  filepath: absolute or relative path name of the ADSC test image to be analyzed
  host: usually "localhost";
  port: port number of image analyzer http service
"""
  do_main(filepath, force_binning, host, port)
