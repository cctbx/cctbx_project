# LIBTBX_SET_DISPATCHER_NAME phenix.rna_validate

import sys
from mmtbx.validation.rna_validate import rna_validate

if __name__ == "__main__":
  rna_validate().run(sys.argv[1:])
