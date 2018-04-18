from __future__ import division
'''
Manually define initialisation params of MPI environment.
See https://bitbucket.org/mpi4py/mpi4py/issues/85/manual-finalizing-and-initializing-mpi
and https://github.com/erdc/mpi4py/blob/master/src/rc.py
for details.
'''
import mpi4py
mpi4py.rc.threads = True
mpi4py.rc.thread_level = "funneled"
from mpi4py import MPI

assert MPI.Is_initialized()
assert MPI.Query_thread() == MPI.THREAD_FUNNELED

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
from libtbx.development.timers import Profiler

import boost.python
ext_omp = boost.python.import_ext("scitbx_examples_strumpack_solver_ext")
ext_mpi = boost.python.import_ext("scitbx_examples_strumpack_mpi_dist_solver_ext")
import sys
import numpy as np

import scipy.sparse as sps

if rank==0:
  A_mat = np.loadtxt(sys.argv[1],dtype={'names':('rows','cols','vals'),'formats':('i8','i8','f8')})
  b_vec = np.loadtxt(sys.argv[2])#,dtype={'names':('vals'),'formats':('f8')})
  n_rows = len(b_vec)
  n_cols = n_rows
  nnz = len(A_mat['vals'])
  print n_rows, n_rows/size
  #Convert the sparse CSR to flex doubles, then use them to solve using the implemented framework
  A_sp = sps.csr_matrix((A_mat['vals'],(A_mat['rows'],A_mat['cols'])))

  #Check if upper/lower triangular only, and generate full if so
  tu=sps.triu(A_sp)
  tl=sps.tril(A_sp)
  sd=sps.diags(A_sp.diagonal())

  A_spS = A_sp
  if tu.nnz == sd.getnnz() or tl.nnz == sd.getnnz():
    A_spS = A_sp + A_sp.transpose() - sd

  import numpy as np
  row_idx_split = np.array_split(np.arange(n_rows),size)
  len_row_idx_split = flex.int( np.cumsum( np.append([0], [len(i) for i in row_idx_split]) ).tolist() )
  P = Profiler("STRUMPACK_OMP")
  A_indptr = flex.int( A_spS.indptr )
  A_indices = flex.int( A_spS.indices )
  A_data = flex.double( A_spS.data )
  b = flex.double( b_vec )
  res_strum_omp = ext_omp.strumpack_solver(
                    n_rows, n_cols,
                    A_indptr, A_indices,
                    A_data, b,
                    ext_omp.scotch,
                    ext_omp.auto
                 )
  del P
  P = Profiler("EIGEN_LDLT")
  res_eig_ldlt = ext_omp.eigen_solver(1, n_rows, n_cols, A_indptr, A_indices, A_data, b)
  del P
  P = Profiler("EIGEN_BICGSTAB")
  res_eig_bicgstab = ext_omp.eigen_solver(2, n_rows, n_cols, A_indptr, A_indices, A_data, b)
  del P

else:
  A_spS=None
  row_idx_split=None
  len_row_idx_split=None
  b_vec=None
  n_cols=None

if size>1:
  #Broadcast data to each rank
  A_spS = comm.bcast(A_spS, root=0)
  row_idx_split = comm.bcast(row_idx_split, root=0)
  len_row_idx_split = comm.bcast(len_row_idx_split, root=0)
  b_vec = comm.bcast(b_vec, root=0)
  n_cols = comm.bcast(n_cols, root=0)

  #Take subset of data for each rank
  A_row_offset = flex.int(A_spS[row_idx_split[rank],:].indptr)
  A_col_offset = flex.int(A_spS[row_idx_split[rank],:].indices)
  A_values = flex.double(A_spS[row_idx_split[rank],:].data)
  b = flex.double(b_vec[row_idx_split[rank]])

  P = Profiler("STRUMPACK_MPI_DIST_RANK=%d"%rank)
  res_strum_mpi_local = ext_mpi.strumpack_mpi_dist_solver(len(row_idx_split[rank]), n_cols, comm, A_row_offset, A_col_offset, A_values, b, len_row_idx_split, ext_mpi.scotch, ext_mpi.auto)
  strum_result_mpi_list = comm.gather(res_strum_mpi_local.x, root=0)
  del P
  if rank==0:
    strum_result_mpi = flex.double()
    for l in strum_result_mpi_list:
      strum_result_mpi.extend(l)

MPI.Finalize()
if rank==0:
  strum_result_omp = res_strum_omp.x
  eig_result_ldlt = res_eig_ldlt.x
  eig_result_bicgstab = res_eig_bicgstab.x

  import cPickle as pickle
  pickle.dump( strum_result_omp,    open( "strum_result_omp.pickle", "wb" ) )
  pickle.dump( strum_result_mpi,    open( "strum_result_mpi.pickle", "wb" ) )
  pickle.dump( eig_result_ldlt,     open( "eig_result_ldlt.pickle", "wb" ) )
  pickle.dump( eig_result_bicgstab, open( "eig_result_bicgstab.pickle", "wb" ) )
  for ii in xrange(len(strum_result_omp)):
    assert approx_equal( strum_result_omp[ii], eig_result_ldlt[ii]  )
    assert approx_equal( strum_result_omp[ii], eig_result_bicgstab[ii]  )
    if size>1:
      assert approx_equal( strum_result_omp[ii], strum_result_mpi[ii]  )
  print "Strumpack solutions agree with Eigen"
