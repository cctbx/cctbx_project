from __future__ import division
import MySQLdb

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
