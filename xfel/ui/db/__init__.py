from __future__ import division
import MySQLdb

class db_credentials(object):
  def __init__(self, host = "",
               username = "",
               password = "",
               db = ""):
    self.host = host
    self.username = username
    self.password = password
    self.db = db

def get_db_connection(credentials, block=True):
  if credentials.host == "":
    credentials.host = "psdb-user.slac.stanford.edu"

  while True:
    try:
      dbobj=MySQLdb.connect(passwd=credentials.password,user=credentials.username,host=credentials.host,db=credentials.db)
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
