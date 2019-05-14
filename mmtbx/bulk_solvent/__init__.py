from __future__ import division
from __future__ import absolute_import
import cctbx.array_family.flex # import dependency

import boost.python
ext = boost.python.import_ext("mmtbx_bulk_solvent_ext")
from mmtbx_bulk_solvent_ext import *
from mmtbx.bulk_solvent.multi_mask_bulk_solvent import multi_mask_bulk_solvent # special import
