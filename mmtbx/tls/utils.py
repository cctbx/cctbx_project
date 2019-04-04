from __future__ import division

# Import dependencies
import scitbx.array_family.flex
from mmtbx.tls.decompose import decompose_tls_matrices # This needs to be imported so that classes of this can be returned by TLSMatrices

import boost.python
ext = boost.python.import_ext("mmtbx_tls_utils_ext")
from mmtbx_tls_utils_ext import *

def _indexing__iter__(self):
    return (self[i] for i in range(self.size()))
TLSAmplitudes.__iter__ = _indexing__iter__
TLSMatricesAndAmplitudesList.__iter__ = _indexing__iter__
