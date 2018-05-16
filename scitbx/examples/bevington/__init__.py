from __future__ import division
from scitbx.lstbx import normal_eqns # import dependency
import boost.python
ext = boost.python.import_ext("scitbx_examples_bevington_ext")
from scitbx_examples_bevington_ext import *

#inject python method into non_linear_ls_eigen_wrapper
class _(boost.python.injector, ext.non_linear_ls_eigen_wrapper):

  def get_eigen_summary(self): # ported from c++ code to encapsulate printing
    import StringIO
    S = StringIO.StringIO()
    assert self.solved()
    nm_ncols = self.get_normal_matrix_ncols()
    matsize = nm_ncols * (nm_ncols+1)/2
    print >>S,"Number of parameters      %12ld"%(self.n_parameters())
    print >>S,"Normal matrix square size %12ld"%(nm_ncols * nm_ncols)
    print >>S,"Upper triangle size       %12ld"%matsize
    nonZeros = self.get_normal_matrix_nnonZeros()
    percentNZ = 100. * nonZeros/float(matsize)
    #LnonZeros = self.get_lower_cholesky_nnonZeros()
    #LpercentNZ = 100. * LnonZeros/float(matsize)
    print >>S,"Normal matrix non-zeros   %12ld, %6.2f%%"%(nonZeros,percentNZ)
    #print >>S,"Cholesky factor non-zeros %12ld, %6.2f%%"%(LnonZeros,LpercentNZ)
    return S.getvalue()

try:
  #ext = boost.python.import_ext("scitbx_examples_strumpack_solver_ext")
  #from scitbx_examples_strumpack_solver_ext import *

#inject python method into non_linear_ls_eigen_wrapper
  class _(boost.python.injector, ext.non_linear_ls_strumpack_wrapper):
    def get_strumpack_summary(self): # ported from c++ code to encapsulate printing
      import StringIO
      S = StringIO.StringIO()
      assert self.solved()
      nm_ncols = self.get_normal_matrix_ncols()
      matsize = nm_ncols * (nm_ncols+1)/2
      print >>S,"Number of parameters      %12ld"%(self.n_parameters())
      print >>S,"Normal matrix square size %12ld"%(nm_ncols * nm_ncols)
      print >>S,"Upper triangle size       %12ld"%matsize
      nonZeros = self.get_normal_matrix_nnonZeros()
      percentNZ = 100. * nonZeros/float(matsize)
      #LnonZeros = self.get_lower_cholesky_nnonZeros()
      #LpercentNZ = 100. * LnonZeros/float(matsize)
      print >>S,"Normal matrix non-zeros   %12ld, %6.2f%%"%(nonZeros,percentNZ)
      #print >>S,"Cholesky factor non-zeros %12ld, %6.2f%%"%(LnonZeros,LpercentNZ)
      return S.getvalue()

except Exception as e:
  #print ("STRUMPACK unavailable: %s"%e)
  pass
