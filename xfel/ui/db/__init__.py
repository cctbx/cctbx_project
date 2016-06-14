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
    self.app = app
    self.id = id
    self.table_name = table_name

    if id is None:
      # add the new items to the db
      query = "INSERT INTO `%s` " % self.table_name
      keys = "("
      vals = "VALUES ("
      comma = ""
      for key, item in kwargs.iteritems():
        setattr(self, key, item)
        keys += comma + key
        vals += comma + "'%s'"%item
        comma = ", "
      keys += ")"
      vals += ")"
      query += keys + " " + vals
      cursor = self.app.dbobj.cursor()
      cursor.execute(query)
      self.app.dbobj.commit()
      self.id = cursor.lastrowid
    else:
      query = "SHOW COLUMNS FROM `%s`" % self.table_name
      cursor = self.app.dbobj.cursor()
      cursor.execute(query)
      columns = [c[0] for c in cursor.fetchall()]

      query = "SELECT * FROM `%s` WHERE id = %d" % (self.table_name, id)
      cursor = self.app.dbobj.cursor()
      cursor.execute(query)
      data = cursor.fetchall()[0]

      for key, value in zip(columns, data):
        setattr(self, key, value)
