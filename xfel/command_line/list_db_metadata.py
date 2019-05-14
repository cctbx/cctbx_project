from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME cxi.list_db_metadata
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

import libtbx.phil
from libtbx.utils import Sorry
import sys

def run(args):
  master_phil = libtbx.phil.parse("""
    db {
      host = None
        .type = str
      name = None
        .type = str
      table_name = None
        .type = str
      user = None
        .type = str
      password = None
        .type = str
    }
  """)

  if (__name__ == "__main__") :
    user_phil = []
    for arg in args :
      try :
        user_phil.append(libtbx.phil.parse(arg))
      except RuntimeError as e :
        raise Sorry("Unrecognized argument '%s' (error: %s)" % (arg, str(e)))
    params = master_phil.fetch(sources=user_phil).extract()

    try:
      from cxi_xdr_xes.cftbx.cspad_ana import db as db
    except ImportError:
      raise Sorry("Trial logging not supported for this installation. Conact the developers for access.")

    dbobj = db.dbconnect(host=params.db.host, db=params.db.name, username=params.db.user, password=params.db.password)
    db.list_db_metadata(dbobj, params.db.table_name)

if (__name__ == "__main__") :
  run(sys.argv[1:])
