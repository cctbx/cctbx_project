from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.ui_server

from libtbx.phil import parse
from libtbx.utils import Sorry
import sys, os, time
from xfel.ui import db_phil_str
from libtbx import easy_run
from xfel.ui.db.xfel_db import db_application

help_message = """

XFEL UI MySQL database server wrapper. Program to initialize and start a mysql database from user space inside the requested folder.

Example:
cctbx.xfel.ui_server db.port=3307 db.server.basedir=`pwd`/MySql db.user=aaron db.name=experiment42

This will look for `pwd`/MSql/my.cnf. If not found, the program will initialize the database, prompt for the entry of a root password, then create the requested database and user account on that database in the folder specified with basedir.

The server will run until a KeyboardInterrupt or other termination signal is sent, then the program will shut down the mysql server before exiting.

While the server is running, the user can connect to with with the xfel gui by providing the hostname the server is running on and appropiate credentials.

"""

phil_scope = parse(db_phil_str)

default_cnf = \
"""
[mysqld]
basedir={basedir}
datadir={basedir}{sep}data
socket={basedir}{sep}mysql.sock
port={port}
# Disabling symbolic-links is recommended to prevent assorted security risks
symbolic-links=0
max_connections=10000
default_authentication_plugin=mysql_native_password

[mysqld_safe]
log-error={basedir}{sep}mysqld.log
pid-file={basedir}{sep}mysqld.pid

[client]
protocol=tcp
"""

def run(args):
  user_phil = []
  if '--help' in args or '-h' in args:
    print(help_message)
    phil_scope.show()
    return
  for arg in args:
    try:
      user_phil.append(parse(arg))
    except Exception as e:
      raise Sorry("Unrecognized argument %s"%arg)

  print("Running ui_server.py")

  def _get_password(params):
    if params.db.server.root_password is None:
      import getpass
      print ("You must specify a root password")
      rootpw1 = getpass.getpass()
      print ("Re-enter password")
      rootpw2 = getpass.getpass()
      if rootpw1 != rootpw2:
        raise Sorry("Passwords don't match")
      params.db.server.root_password = rootpw1
  params = phil_scope.fetch(sources=user_phil).extract()

  cnf_path = os.path.join(params.db.server.basedir, 'my.cnf')

  initialize = not os.path.exists(params.db.server.basedir)
  if initialize:
    assert params.db.user is not None and len(params.db.user) > 0 and \
           params.db.name is not None and len(params.db.name) > 0
    print("Initializing")
    _get_password(params)

    print ("Initializing database")
    os.makedirs(params.db.server.basedir)
    with open(cnf_path, 'w') as f:
      f.write(default_cnf.format(basedir = params.db.server.basedir, sep = os.path.sep, port = params.db.port))

    assert easy_run.call("mysqld --defaults-file=%s --initialize-insecure"%(cnf_path)) == 0

  elif params.db.server.prompt_for_root_password and params.db.server.root_password is not None:
    import getpass
    print ("please enter root password to raise the connection")
    params.db.server.root_password = getpass.getpass()

  print ("Starting server")
  assert os.path.exists(cnf_path)
  server_process = easy_run.subprocess.Popen(["mysqld", "--defaults-file=%s"%(cnf_path)])

  print ("Sleeping a few seconds to let server start up...")
  time.sleep(5) # let server start up

  params.db.host = '127.0.0.1'
  if initialize:
    app = db_application(params)
    print("attempting to change USER information for DB")
    app.execute_query("ALTER USER 'root'@'localhost' IDENTIFIED BY '%s'"%(params.db.server.root_password))

    print ("Creating empty database %s"%params.db.name)
    app.execute_query("CREATE DATABASE %s"%params.db.name)

    print ("Creating new user %s"%params.db.user)
    app.execute_query("CREATE USER '%s'@'%%' IDENTIFIED BY '%s'"%(params.db.user, params.db.password

    print ("Setting permissions")
    app.execute_query("GRANT ALL PRIVILEGES ON %s . * TO '%s'@'%%'"%(params.db.name, params.db.user))
    rootpw = rootpw1
    app.execute_query("FLUSH PRIVILEGES")
    app.execute_query("UPDATE mysql.user SET Super_Priv='Y' WHERE user='%s' AND host='%%'"%params.db.user)
    app.execute_query("FLUSH PRIVILEGES")
    print ("Initialized")
  else:
    print("instantiating existing DB")
    app = db_application(params)

  print ("Raising max connections")
  app.execute_query("SET GLOBAL max_connections=50000")

  # remove root password from params
  params.db.server.root_password = None
  try:
    while True:
      if server_process.poll() is not None:
        print ("Server exited")
        return
      time.sleep(1)
  except KeyboardInterrupt:
    print ("Shutting down")
  except Exception as e:
    print ("Unhandled exception, shutting down")
    print (str(e))

  server_process.terminate()

if __name__ == '__main__':
  run(sys.argv[1:])
