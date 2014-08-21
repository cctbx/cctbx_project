from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME cxi.list_db_metadata
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

from libtbx.utils import Sorry

if (__name__ == "__main__") :
  try:
    from cxi_xdr_xes.cftbx.cspad_ana import db as db
  except ImportError:
    raise Sorry("Trial logging not supported for this installation. Conact the developers for access.")

  db.list_db_metadata()
