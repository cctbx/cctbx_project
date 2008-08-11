from iotbx.cns import sdb_writer
from cctbx import xray
xray.structure.as_cns_sdb_file = sdb_writer.xray_structure_as_cns_sdb_file
