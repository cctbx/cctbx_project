from __future__ import absolute_import, division, print_function

from six.moves import cStringIO as StringIO

import os, time
import libtbx.load_env
from libtbx import xmlrpc_utils

def exercise():
  ext_module = libtbx.env.find_in_repositories(
    relative_path="libtbx/xmlrpc_server_example.py",
    test=os.path.isfile)
  full_cmd = "libtbx.python %s" % ext_module
  log = StringIO()
  server = xmlrpc_utils.external_program_server(command=full_cmd,
                                                program_id="TEST",
                                                timeout=None,
                                                cache_requests=True,
                                                log=log)
  server.echo_test()
  time.sleep(1)
  run_server_tests(server)
  port = server.get_port()
  assert log.getvalue() == """XML-RPC server started on port %d
hello, world!
hello, world!
quitting
""" % port
  server.restart()
  time.sleep(1)
  run_server_tests(server)
  print("OK")

def run_server_tests(server):
  assert server.echo_test() == True
  assert server.is_alive() == True
  server.quit()
  assert server.is_alive() == False
  assert server.echo_test() is None

if __name__ == "__main__" :
  exercise()
