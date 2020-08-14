from __future__ import absolute_import, division, print_function
from scitbx.lstbx import normal_eqns # import dependency
import boost_adaptbx.boost.python as bp
ext = bp.import_ext("scitbx_examples_bevington_ext")
from scitbx_examples_bevington_ext import *

#inject python method into non_linear_ls_eigen_wrapper
@bp.inject_into(ext.non_linear_ls_eigen_wrapper)
class _():

  def get_eigen_summary(self): # ported from c++ code to encapsulate printing
    from six.moves import StringIO
    S = StringIO()
    assert self.solved()
    nm_ncols = self.get_normal_matrix_ncols()
    matsize = nm_ncols * (nm_ncols+1)/2
    print("Number of parameters      %12ld"%(self.n_parameters()), file=S)
    print("Normal matrix square size %12ld"%(nm_ncols * nm_ncols), file=S)
    print("Upper triangle size       %12ld"%matsize, file=S)
    nonZeros = self.get_normal_matrix_nnonZeros()
    percentNZ = 100. * nonZeros/float(matsize)
    LnonZeros = self.get_lower_cholesky_nnonZeros()
    LpercentNZ = 100. * LnonZeros/float(matsize)
    print("Normal matrix non-zeros   %12ld, %6.2f%%"%(nonZeros,percentNZ), file=S)
    print("Cholesky factor non-zeros %12ld, %6.2f%%"%(LnonZeros,LpercentNZ), file=S)
    return S.getvalue()
