from __future__ import absolute_import, division, print_function
from libtbx.utils import Sorry

try:
  import MySQLdb
except ImportError as e:
  pass

def get_run_path(rootpath, trial, rungroup, run):
  import os
  try:
    return os.path.join(rootpath, "r%04d"%int(run.run), "%03d_rg%03d"%(trial.trial, rungroup.id))
  except ValueError:
    return os.path.join(rootpath, run.run, "%03d_rg%03d"%(trial.trial, rungroup.id))

def get_db_connection(params, block=True):
  if params.db.password is None:
    password = ""
  else:
    password = params.db.password

  retry_count = 0
  retry_max = 20
  sleep_time = 0.1
  while retry_count < retry_max:
    try:
      dbobj=MySQLdb.connect(passwd=password,user=params.db.user,host=params.db.host,db=params.db.name,port=params.db.port)
      return dbobj
    except Exception as e:
      retry_count += 1
      if not block: raise e
      if "Too many connections" in str(e):
        print("Too many connections, retry", retry_count)
      elif  "Can't connect to MySQL server" in str(e):
        print("Couldn't connect to MYSQL, retry", retry_count)
      else:
        raise e
      import time
      time.sleep(sleep_time)
      sleep_time *= 2
  raise Sorry("Couldn't execute connect to MySQL. Too many reconnects.")

class db_proxy(object):
  def __init__(self, app, table_name, id = None, **kwargs):
    self._db_dict = {}
    self.app = app
    self.id = id
    self.table_name = table_name

    if id is None:
      # add the new items to the db
      query = "INSERT INTO `%s` " % self.table_name
      keys = []
      vals = []
      for key, value in kwargs.iteritems():
        assert key in app.columns_dict[table_name]
        keys.append(key)
        self._db_dict[key] = value
        if isinstance(value, bool):
          value = "'%s'"%int(value)
        elif value is None or value == "None" or value == "":
          value = "NULL"
        else:
          value = "'%s'"%str(value)
        vals.append(value)
      query += "(%s) VALUES (%s)"%(", ".join(keys), ", ".join(vals))
      cursor = app.execute_query(query, commit=True)
      self.id = cursor.lastrowid

      for key in app.columns_dict[table_name]:
        if key not in self._db_dict:
          self._db_dict[key] = None
    else:
      keys = []
      for key in app.columns_dict[table_name]:
        if key in kwargs:
          self._db_dict[key] = kwargs.pop(key)
        else:
          keys.append(key)
      assert len(kwargs) == 0
      if len(keys) > 0:
        query = "SELECT %s FROM `%s` WHERE id = %d" % (", ".join(keys), table_name, id)
        results = app.execute_query(query).fetchall()[0]
        for key, value in zip(keys, results):
          self._db_dict[key] = value

  def __getattr__(self, key):
    # Called if the property is not found
    assert hasattr(self, '_db_dict')
    if key not in self._db_dict:
      print(self.table_name, key, 'error!', self._db_dict)
      raise AttributeError(key)

    return self._db_dict[key]

  def __setattr__(self, key, value):
    # Need to test for _db_dict to avoid infinite loop
    #  Test key == "_db_dict" to allow creating _db_dict
    #  Test hasattr(self, "_db_dict") in case dbproxy.__init__ hasn't been called
    #  Test key not in self._db_dict to allow setting member variables not in database dictionary
    if key == "_db_dict" or not hasattr(self, "_db_dict") or key not in self._db_dict:
      super(db_proxy, self).__setattr__(key, value)
      return

    if isinstance(value, bool):
      v = "%s"%int(value)
    elif value is None or value == "None" or value == "":
      v = "NULL"
    else:
      v = "'%s'"%value

    query = "UPDATE `%s` SET %s = %s WHERE id = %d"% (
      self.table_name, key, v, self.id)
    self.app.execute_query(query, commit=True)
    self._db_dict[key] = value
