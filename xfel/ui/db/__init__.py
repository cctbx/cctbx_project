from __future__ import division
try:
  import MySQLdb
except ImportError, e:
  from libtbx.utils import Sorry
  raise Sorry("Mysql is not installed")

def get_run_path(rootpath, trial, rungroup, run):
  import os
  return os.path.join(rootpath, "r%04d"%run.run, "%03d_rg%03d"%(trial.trial, rungroup.id))

def get_db_connection(params, block=True):
  if params.db.password is None:
    password = ""
  else:
    password = params.db.password

  retry_count = 0
  retry_max = 10
  sleep_time = 0.1
  while retry_count < retry_max:
    try:
      dbobj=MySQLdb.connect(passwd=password,user=params.db.user,host=params.db.host,db=params.db.name)
      return dbobj
    except Exception,e:
      retry_count += 1
      if not block: raise e
      if "Too many connections" in str(e):
        print "Too many connections, retry", retry_count
      elif  "Can't connect to MySQL server" in str(e):
        print "Couldn't connect to MYSQL, retry", retry_count
      else:
        raise e
      import time
      time.sleep(sleep_time)
      sleep_time *= 2
  raise Sorry("Couldn't execute connect to MySQL. Too many reconnects.")

class db_proxy(object):
  def __init__(self, app, table_name, id = None, **kwargs):
    self._db_columns = []
    self.app = app
    self.id = id
    self.table_name = table_name

    if id is None:
      # add the new items to the db
      query = "INSERT INTO `%s` " % self.table_name
      keys = []
      vals = []
      for key, value in kwargs.iteritems():
        keys.append(key)
        if isinstance(value, bool):
          value = "'%s'"%int(value)
        elif value is None or value == "None" or value == "":
          value = "NULL"
        else:
          value = "'%s'"%str(value)
        vals.append(value)
      query += "(%s) VALUES (%s)"%(", ".join(keys), ", ".join(vals))
      cursor = self.app.execute_query(query, commit=True)
      self.id = cursor.lastrowid

    self._db_columns = app.columns_dict[table_name]

  def __getattr__(self, key):
    # Called if the property is not found
    if key not in self._db_columns:
      print self.table_name, key, 'error!', self._db_columns
      raise AttributeError()

    query = "SELECT %s FROM `%s` WHERE id = %d" % (key, self.table_name, self.id)
    cursor = self.app.execute_query(query)
    return cursor.fetchall()[0][0]

  def __setattr__(self, key, value):
    # Need to test for _db_columns to avoid infinite loop
    #  Test key == "_db_columns" to allow creating _db_columns
    #  Test hasattr(self, "_db_columns") in case dbproxy.__init__ hasn't been called
    #  Test key not in self._db_columns to allow setting member variables not in database dictionary
    if key == "_db_columns" or not hasattr(self, "_db_columns") or key not in self._db_columns:
      super(db_proxy, self).__setattr__(key, value)
      return

    if isinstance(value, bool):
      value = "%s"%int(value)
    elif value is None or value == "None" or value == "":
      value = "NULL"
    else:
      value = "'%s'"%value

    query = "UPDATE `%s` SET %s = %s WHERE id = %d"% (
      self.table_name, key, value, self.id)
    print query
    self.app.execute_query(query, commit=True)
