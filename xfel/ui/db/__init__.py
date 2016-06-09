from __future__ import division

def get_db_connection(params):
  class Dummy(object):
    def cursor(self):
      pass

  return Dummy()
