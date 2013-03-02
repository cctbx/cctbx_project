from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME cxi.list_db_metadata
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

from cxi_xdr_xes.cftbx.cspad_ana import db as db

if (__name__ == "__main__") :
  db.list_db_metadata()
