# LIBTBX_SET_DISPATCHER_NAME distl.thin_client
from spotfinder.servers.thin_client import do_main

if __name__=="__main__":
  "Client is intended to be used with distl.mp_spotfinder_server_read_file"
  import sys
  try:
    filepath, host, port = sys.argv[1:4]
    port = int(port)
  except Exception:
    print """
Usage:
distl.thin_client <filepath> <host> <port>
Three mandatory arguments:
  filepath: absolute or relative path name of the ADSC test image to be analyzed
  host: usually "localhost";
  port: port number of image analyzer http service
"""
  do_main(filepath, host, port)
