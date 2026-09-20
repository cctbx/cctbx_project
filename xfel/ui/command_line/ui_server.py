from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME cctbx.xfel.ui_server

from libtbx.phil import parse
from libtbx.utils import Sorry
import sys, os, time, shutil
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

# Note mysqld's own "basedir" is the directory it was installed in, which it
# infers from the location of the executable; it is deliberately not set here.
# db.server.basedir is where this server keeps its data, which is "datadir".
default_cnf = \
"""
[mysqld]
datadir={basedir}{sep}data
socket={socket}
port={port}
max_connections=10000

[mysqld_safe]
log-error={basedir}{sep}mysqld.log
pid-file={basedir}{sep}mysqld.pid

[client]
protocol=tcp
"""

def find_mysqld():
  """Locate the mysqld executable. The libtbx dispatchers only put the build's
  bin directory on PATH, so a mysqld installed alongside the rest of the conda
  dependencies is not found there; fall back to the python prefix's bin."""
  mysqld = shutil.which('mysqld')
  if mysqld: return mysqld
  mysqld = os.path.join(sys.prefix, 'bin', 'mysqld')
  if os.path.exists(mysqld): return mysqld
  raise Sorry("Could not find the mysqld executable. Install it (conda install "
              "mysql-server) or put it on your PATH.")

def socket_path(basedir):
  """Path for the server's unix socket. Clients connect over TCP, but mysqld
  still requires a socket it can create, and unix socket paths are limited to
  107 characters, which a deeply nested basedir can exceed. Fall back to a
  name in the temporary directory, made unique by a hash of the basedir so
  that two servers do not collide."""
  path = os.path.join(basedir, 'mysql.sock')
  if len(path) <= 107: return path
  import hashlib, tempfile
  digest = hashlib.sha256(basedir.encode('utf-8')).hexdigest()[:16]
  return os.path.join(tempfile.gettempdir(), 'cctbx_xfel_mysql_%s.sock'%digest)

def remove_failed_basedir(basedir):
  """Remove a base directory this process has just created, after initializing
  the database in it failed. Leaving it behind would make the next attempt take
  it for an already initialized database, since it holds the my.cnf written
  before initialization. Returns a note to append to the error being reported if
  the directory could not be removed, so that a failed cleanup is never silent."""
  try:
    shutil.rmtree(basedir)
    return ""
  except OSError as e:
    return (" %s could not be cleaned up (%s); remove it by hand before trying "
            "again, or the next attempt will mistake it for an initialized "
            "database."%(basedir, str(e)))

def run(args):
  # This is normally launched into a log file by the GUI, where block buffering
  # would hide the progress messages below until the server exits.
  try:
    sys.stdout.reconfigure(line_buffering=True)
  except AttributeError:
    pass

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
  params = phil_scope.fetch(sources=user_phil).extract()

  if not params.db.server.basedir:
    raise Sorry("db.server.basedir must be specified")
  mysqld = find_mysqld()
  cnf_path = os.path.join(params.db.server.basedir, 'my.cnf')

  initialize = not os.path.exists(params.db.server.basedir)
  if initialize:
    assert params.db.user is not None and len(params.db.user) > 0 and \
           params.db.name is not None and len(params.db.name) > 0
    import getpass
    if params.db.server.root_password:
      rootpw = params.db.server.root_password
    else:
      print("Initializing")
      print("You must specify a root password")
      rootpw1 = getpass.getpass()
      print("Re-enter password")
      rootpw2 = getpass.getpass()
      if rootpw1 != rootpw2:
        raise Sorry("Passwords don't match")
      rootpw = rootpw1

    print("Initializing database")
    os.makedirs(params.db.server.basedir)
    with open(cnf_path, 'w') as f:
      f.write(default_cnf.format(basedir=params.db.server.basedir, sep=os.path.sep,
                                 socket=socket_path(params.db.server.basedir),
                                 port=params.db.port))
    try:
      result = easy_run.call("%s --defaults-file=%s --initialize-insecure"%(mysqld, cnf_path))
    except Exception:
      # Clean up here too: without this a half-initialized directory would be
      # left behind by anything that stops mysqld reporting a return code.
      note = remove_failed_basedir(params.db.server.basedir)
      if note: print(note)
      raise
    if result != 0:
      note = remove_failed_basedir(params.db.server.basedir)
      raise Sorry("Failed to initialize the database in %s (mysqld returned %d). "
                  "See above for the error reported by mysqld.%s"
                  %(params.db.server.basedir, result, note))

  elif params.db.server.prompt_for_root_password:
    if params.db.server.root_password:
      rootpw3 = params.db.server.root_password
    else:
      import getpass
      print("please enter root password to raise the connection")
      rootpw3 = getpass.getpass()

  print("Starting server")
  if not os.path.exists(cnf_path):
    raise Sorry("%s exists but does not contain a my.cnf. Remove it to let the "
                "server initialize a new database there."%params.db.server.basedir)
  server_process = easy_run.subprocess.Popen([mysqld, "--defaults-file=%s"%(cnf_path)])

  print("Sleeping a few seconds to let server start up...")
  time.sleep(5) # let server start up

  params.db.host = '127.0.0.1'
  if initialize:
    new_user = params.db.user
    new_password = params.db.password
    new_db = params.db.name
    params.db.user = 'root'
    params.db.password = ''
    params.db.name = ''
    print("Changing password")
    app = db_application(params)
    app.execute_query("ALTER USER 'root'@'localhost' IDENTIFIED BY '%s'"%(rootpw))
    params.db.password = rootpw
    print("Creating empty database %s"%new_db)
    app.execute_query("CREATE DATABASE %s"%new_db)
    print("Creating new user %s"%new_user)
    app.execute_query("CREATE USER '%s'@'%%' IDENTIFIED BY '%s'"%(new_user, new_password))
    print("Setting permissions")
    app.execute_query("GRANT ALL PRIVILEGES ON %s . * TO '%s'@'%%'"%(new_db, new_user))
    app.execute_query("FLUSH PRIVILEGES")

    # Detect MySQL version to grant the correct privilege for SET GLOBAL
    cursor = app.execute_query("SELECT VERSION()")
    version_str = cursor.fetchall()[0][0]
    major_version = int(version_str.split('.')[0])
    if major_version >= 9:
      app.execute_query("GRANT SYSTEM_VARIABLES_ADMIN ON *.* TO '%s'@'%%'"%(new_user))
    else:
      app.execute_query("UPDATE mysql.user SET Super_Priv='Y' WHERE user='%s' AND host='%%'"%new_user)
    app.execute_query("FLUSH PRIVILEGES")
    print("Initialized")
  else:
    print("Instantiating db query execution driver")
    app = db_application(params)

  if params.db.server.prompt_for_root_password:
    params.db.user = 'root'
    params.db.password = rootpw3

  print("Raising max connections")
  app.execute_query("SET GLOBAL max_connections=50000")

  try:
    while True:
      if server_process.poll() is not None:
        print("Server exited")
        return
      time.sleep(1)
  except KeyboardInterrupt:
    print("Shutting down")
  except Exception as e:
    print("Unhandled exception, shutting down")
    print(str(e))

  server_process.terminate()


if __name__ == '__main__':
  print(sys.argv)
  run(sys.argv[1:])

