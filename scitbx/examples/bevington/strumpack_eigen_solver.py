from __future__ import division
from scitbx.matrix import sqr,col
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
from libtbx.development.timers import Profiler

import boost.python
ext = boost.python.import_ext("scitbx_examples_strumpack_solver_ext")

import sys
import numpy as np
A_mat = np.loadtxt(sys.argv[1],dtype={'names':('rows','cols','vals'),'formats':('i8','i8','f8')})
b_vec = np.loadtxt(sys.argv[2])#,dtype={'names':('vals'),'formats':('f8')})
n_rows = len(b_vec)
n_cols = n_rows
nnz = len(A_mat['vals'])

#Convert the sparse CSR to flex doubles, then use them to solve using the implemented framework
import scipy.sparse as sps
A_sp = sps.csr_matrix((A_mat['vals'],(A_mat['rows'],A_mat['cols']))) 

  #Check if upper/lower triangular only, and generate full if so
tu=sps.triu(A_sp)
tl=sps.tril(A_sp)
sd=sps.diags(A_sp.diagonal())

A_spS = A_sp
if tu.nnz == sd.getnnz() or tl.nnz == sd.getnnz():
  A_spS = A_sp + A_sp.transpose() - sd

A_row_offset = flex.int(A_spS.indptr)
A_col_offset = flex.int(A_spS.indices)
A_values = flex.double(A_spS.data)
b = flex.double(b_vec)

P = Profiler("STRUMPACK")
res_strum = ext.strumpack_solver(0, n_rows, n_cols, A_row_offset, A_col_offset, A_values, b)
P = Profiler("EIGEN")
res_eig = ext.eigen_solver(0, n_rows, n_cols, A_row_offset, A_col_offset, A_values, b)
del P
strum_result = res_strum.x
eig_result = res_eig.x
import cPickle as pickle
pickle.dump( eig_result, open( "res_eigen.pickle", "wb" ) )

for ii in xrange(len(strum_result)):
  assert approx_equal( strum_result[ii], eig_result[ii]  )
print "scitbx.matrix solution and strumpack solution agree"
