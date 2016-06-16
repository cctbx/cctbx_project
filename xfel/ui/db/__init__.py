from __future__ import division
try:
  import MySQLdb
except ImportError, e:
  from libtbx.utils import Sorry
  raise Sorry("Mysql is not installed")

def get_db_connection(params, block=True):
  while True:
    try:
      dbobj=MySQLdb.connect(passwd=params.db.password,user=params.db.user,host=params.db.host,db=params.db.name)
      return dbobj
    except Exception,e:
      if "Too many connections" in e and block:
        print "Too many connections.  Blocking..."
        dbobj = None
        import time
        time.sleep(1)
        continue
      else:
        raise e

class db_proxy(object):
  def __init__(self, app, table_name, id = None, **kwargs):
    self.db_dict = {}
    self.app = app
    self.id = id
    self.table_name = table_name

    if id is None:
      # add the new items to the db
      query = "INSERT INTO `%s` " % self.table_name
      keys = "("
      vals = "VALUES ("
      comma = ""
      for key, value in kwargs.iteritems():
        if isinstance(value, bool):
          value = int(value)
        self.db_dict[key] = value
        keys += comma + key
        vals += comma + "'%s'"%value
        comma = ", "
      keys += ")"
      vals += ")"
      query += keys + " " + vals
      cursor = self.app.execute_query(query, commit=True)
      self.id = cursor.lastrowid
    else:
      query = "SHOW COLUMNS FROM `%s`" % self.table_name
      cursor = self.app.execute_query(query)
      columns = [c[0] for c in cursor.fetchall()]

      query = "SELECT * FROM `%s` WHERE id = %d" % (self.table_name, id)
      cursor = self.app.execute_query(query)
      data = cursor.fetchall()[0]

      for key, value in zip(columns, data):
        if key == 'id':
          continue
        self.db_dict[key] = value

  def __getattr__(self, key):
    # Called if the property is not found
    if key in self.db_dict:
      return self.db_dict[key]
    else:
      raise AttributeError()

  def __setattr__(self, key, value):
    # Need to test for db_dict to avoid infinite loop
    if key == "db_dict" or key not in self.db_dict:
      super(db_proxy, self).__setattr__(key, value)
      return

    self.db_dict[key] = value

    if isinstance(value, bool):
      value = "%s"%int(value)
    elif value is None:
      value = "NULL"
    else:
      value = "'%s'"%value

    query = "UPDATE `%s` SET %s = %s WHERE id = %d"% (
      self.table_name, key, value, self.id)
    print query
    self.app.execute_query(query, commit=True)
