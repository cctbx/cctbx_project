from __future__ import absolute_import, division, print_function
from cctbx import miller, sgtbx
from cctbx.array_family import flex
from libtbx import easy_mp
from scipy import sparse
import math
import time
import sys
import numpy as np

try:
    from line_profiler import LineProfiler as profile
except ImportError:
    profile = None
class Reproducer:
    def _lattice_lower_upper_index(self, lattice_id):
        lower_index = self._lattices[lattice_id]
        upper_index = None
        if lattice_id < len(self._lattices) - 1:
            upper_index = self._lattices[lattice_id + 1]
        else:
            assert lattice_id == len(self._lattices) - 1
        return lower_index, upper_index

    def compute_rij_wij_Gildea_old(self, use_cache=True, do_one_row_i=None):
        """
        Compute the rij_wij matrix.
        This is an implementation from commit b5065f39e2 on branch merge_polar,
        before we started optimizing the Miller index matching. We want these
        results available temporarily for validation.

        """
        n_lattices = self._lattices.size
        n_sym_ops = len(self.sym_ops)

        NN = n_lattices * n_sym_ops

        self.rij_matrix = flex.double(flex.grid(NN, NN), 0.0)
        if self._weights is None:
            self.wij_matrix = None
        else:
            self.wij_matrix = flex.double(flex.grid(NN, NN), 0.0)

        indices = {}
        space_group_type = self._data.space_group().type()
        for cb_op in self.sym_ops:
            cb_op = sgtbx.change_of_basis_op(cb_op)
            indices_reindexed = cb_op.apply(self._data.indices())
            miller.map_to_asu(space_group_type, False, indices_reindexed)
            indices[cb_op.as_xyz()] = indices_reindexed

        def _compute_rij_matrix_one_row_block(i):
            rij_cache = {}

            n_sym_ops = len(self.sym_ops)
            NN = n_lattices * n_sym_ops

            rij_row = []
            rij_col = []
            rij_data = []
            if self._weights is not None:
                wij_row = []
                wij_col = []
                wij_data = []
            else:
                wij = None

            i_lower, i_upper = self._lattice_lower_upper_index(i)
            intensities_i = self._data.data()[i_lower:i_upper]

            for j in range(n_lattices):

                j_lower, j_upper = self._lattice_lower_upper_index(j)
                intensities_j = self._data.data()[j_lower:j_upper]

                for k, cb_op_k in enumerate(self.sym_ops):
                    cb_op_k = sgtbx.change_of_basis_op(cb_op_k)

                    indices_i = indices[cb_op_k.as_xyz()][i_lower:i_upper]

                    for kk, cb_op_kk in enumerate(self.sym_ops):
                        if i == j and k == kk:
                            # don't include correlation of dataset with itself
                            continue
                        cb_op_kk = sgtbx.change_of_basis_op(cb_op_kk)

                        ik = i + (n_lattices * k)
                        jk = j + (n_lattices * kk)

                        key = (i, j, str(cb_op_k.inverse() * cb_op_kk))
                        if use_cache and key in rij_cache:
                            cc, n = rij_cache[key]
                        else:
                            indices_j = indices[cb_op_kk.as_xyz()][j_lower:j_upper]

                            matches = miller.match_indices(indices_i, indices_j)
                            pairs = matches.pairs()
                            isel_i = pairs.column(0)
                            isel_j = pairs.column(1)
                            isel_i = isel_i.select(
                                self._patterson_group.epsilon(indices_i.select(isel_i))
                                == 1
                            )
                            isel_j = isel_j.select(
                                self._patterson_group.epsilon(indices_j.select(isel_j))
                                == 1
                            )
                            corr = flex.linear_correlation(
                                intensities_i.select(isel_i),
                                intensities_j.select(isel_j),
                            )

                            if corr.is_well_defined():
                                cc = corr.coefficient()
                                n = corr.n()
                            else:
                                cc = None
                                n = None
                            rij_cache[key] = (cc, n)

                        if (
                            n is None
                            or cc is None
                            or (self._min_pairs is not None and n < self._min_pairs)
                        ):
                            continue

                        if self._weights == "count":
                            wij_row.extend([ik, jk])
                            wij_col.extend([jk, ik])
                            wij_data.extend([n, n])
                        elif self._weights == "standard_error":
                            assert n > 2
                            # http://www.sjsu.edu/faculty/gerstman/StatPrimer/correlation.pdf
                            se = math.sqrt((1 - cc ** 2) / (n - 2))
                            wij = 1 / se
                            wij_row.extend([ik, jk])
                            wij_col.extend([jk, ik])
                            wij_data.extend([wij, wij])

                        rij_row.append(ik)
                        rij_col.append(jk)
                        rij_data.append(cc)

            rij = sparse.coo_matrix((rij_data, (rij_row, rij_col)), shape=(NN, NN))
            if self._weights is not None:
                wij = sparse.coo_matrix((wij_data, (wij_row, wij_col)), shape=(NN, NN))

            return rij, wij

        args = [(i,) for i in range(n_lattices)]
        if do_one_row_i is not None:
            results = [_compute_rij_matrix_one_row_block((do_one_row_i))]
        else:
            results = easy_mp.parallel_map(
                _compute_rij_matrix_one_row_block,
                args,
                processes=self._nproc,
                iterable_type=easy_mp.posiargs,
                method="multiprocessing",
            )

        rij_matrix = None
        wij_matrix = None
        for i, (rij, wij) in enumerate(results):
            if rij_matrix is None:
                rij_matrix = rij
            else:
                rij_matrix += rij
            if wij is not None:
                if wij_matrix is None:
                    wij_matrix = wij
                else:
                    wij_matrix += wij

        self.rij_matrix = flex.double(rij_matrix.todense())
        if wij_matrix is not None:
            import numpy as np

            self.wij_matrix = flex.double(wij_matrix.todense().astype(np.float64))

        return self.rij_matrix, self.wij_matrix

    def compute_rij_wij_Gildea(self,
            use_cache=True,
            lineprof=-1,
            do_one_row_i=None):
        """Compute the rij_wij matrix."""
        n_lattices = self._lattices.size
        n_sym_ops = len(self.sym_ops)

        NN = n_lattices * n_sym_ops

        self.rij_matrix = flex.double(flex.grid(NN, NN), 0.0)
        if self._weights is None:
            self.wij_matrix = None
        else:
            self.wij_matrix = flex.double(flex.grid(NN, NN), 0.0)

        indices = {}
        space_group_type = self._data.space_group().type()
        for cb_op in self.sym_ops:
            cb_op = sgtbx.change_of_basis_op(cb_op)
            indices_reindexed = cb_op.apply(self._data.indices())
            miller.map_to_asu(space_group_type, False, indices_reindexed)
            indices[cb_op.as_xyz()] = indices_reindexed

        def _compute_rij_matrix_one_row_block(i):
            import time
            t0 = time.time()
            rij_cache = {}

            n_sym_ops = len(self.sym_ops)
            NN = n_lattices * n_sym_ops

            rij_row = []
            rij_col = []
            rij_data = []
            if self._weights is not None:
                wij_row = []
                wij_col = []
                wij_data = []
            else:
                wij = None

            i_lower, i_upper = self._lattice_lower_upper_index(i)
            intensities_i = self._data.data()[i_lower:i_upper]

            idx_matchers = []
            for op in self.sym_ops:
                cb_op = sgtbx.change_of_basis_op(op)
                indices_i = indices[cb_op.as_xyz()][i_lower:i_upper]
                idx_matchers.append(miller.match_indices(indices_i))


            for j in range(n_lattices):

                j_lower, j_upper = self._lattice_lower_upper_index(j)
                intensities_j = self._data.data()[j_lower:j_upper]

                for k, cb_op_k in enumerate(self.sym_ops):
                    cb_op_k = sgtbx.change_of_basis_op(cb_op_k)

                    indices_i = indices[cb_op_k.as_xyz()][i_lower:i_upper]
                    matcher = idx_matchers[k]

                    for kk, cb_op_kk in enumerate(self.sym_ops):
                        if i == j and k == kk:
                            # don't include correlation of dataset with itself
                            continue
                        cb_op_kk = sgtbx.change_of_basis_op(cb_op_kk)

                        ik = i + (n_lattices * k)
                        jk = j + (n_lattices * kk)

                        key = (i, j, str(cb_op_k.inverse() * cb_op_kk))
                        if use_cache and key in rij_cache:
                            cc, n = rij_cache[key]
                        else:
                            indices_j = indices[cb_op_kk.as_xyz()][j_lower:j_upper]

                            matcher.match_cached(indices_j)
                            pairs = matcher.pairs()
                            isel_i = pairs.column(0)
                            isel_j = pairs.column(1)
                            isel_i = isel_i.select(
                                self._patterson_group.epsilon(indices_i.select(isel_i))
                                == 1
                            )
                            isel_j = isel_j.select(
                                self._patterson_group.epsilon(indices_j.select(isel_j))
                                == 1
                            )
                            corr = flex.linear_correlation(
                                intensities_i.select(isel_i),
                                intensities_j.select(isel_j),
                            )

                            if corr.is_well_defined():
                                cc = corr.coefficient()
                                n = corr.n()
                            else:
                                cc = None
                                n = None

                            rij_cache[key] = (cc, n)

                        if (
                            n is None
                            or cc is None
                            or (self._min_pairs is not None and n < self._min_pairs)
                        ):
                            continue

                        if self._weights == "count":
                            wij_row.extend([ik, jk])
                            wij_col.extend([jk, ik])
                            wij_data.extend([n, n])
                        elif self._weights == "standard_error":
                            assert n > 2
                            # http://www.sjsu.edu/faculty/gerstman/StatPrimer/correlation.pdf
                            se = math.sqrt((1 - cc ** 2) / (n - 2))
                            wij = 1 / se
                            wij_row.extend([ik, jk])
                            wij_col.extend([jk, ik])
                            wij_data.extend([wij, wij])

                        rij_row.append(ik)
                        rij_col.append(jk)
                        rij_data.append(cc)

            rij = sparse.coo_matrix((rij_data, (rij_row, rij_col)), shape=(NN, NN))
            if self._weights is not None:
                wij = sparse.coo_matrix((wij_data, (wij_row, wij_col)), shape=(NN, NN))

            print(time.time()-t0)
            return rij, wij

        args = [(i,) for i in range(n_lattices)]
        if lineprof >= 0:
            assert profile is not None
            pr = profile(_compute_rij_matrix_one_row_block)
            pr.enable()
            _compute_rij_matrix_one_row_block(lineprof)
            pr.disable()
            #pr.print_stats()
            quit()
        elif do_one_row_i is not None:
            results = [_compute_rij_matrix_one_row_block(do_one_row_i)]
        else:
            results = easy_mp.parallel_map(
                _compute_rij_matrix_one_row_block,
                args,
                processes=self._nproc,
                iterable_type=easy_mp.posiargs,
                method="multiprocessing",
            )


        rij_matrix = None
        wij_matrix = None
        for i, (rij, wij) in enumerate(results):
            if rij_matrix is None:
                rij_matrix = rij
            else:
                rij_matrix += rij
            if wij is not None:
                if wij_matrix is None:
                    wij_matrix = wij
                else:
                    wij_matrix += wij

        self.rij_matrix = flex.double(rij_matrix.todense())
        if wij_matrix is not None:
            import numpy as np

            self.wij_matrix = flex.double(wij_matrix.todense().astype(np.float64))

        return self.rij_matrix, self.wij_matrix

    def compute_rij_wij_cplusplus(self, use_cache=True, do_one_row_i=None):
        print("""Compute the rij_wij matrix in C++""")

        n_lattices = self._lattices.size
        n_sym_ops = len(self.sym_ops)
        NN = n_lattices * n_sym_ops

        lower_i = flex.int()
        upper_i = flex.int()
        for lidx in range(self._lattices.size):
          LL,UU = self._lattice_lower_upper_index(lidx)
          lower_i.append(int(LL))
          if UU is None:  UU = self._data.data().size()
          upper_i.append(int(UU))
        indices = {}
        space_group_type = self._data.space_group().type()
        print (space_group_type.lookup_symbol(), "weight mode", self._weights)
        from xfel.merging import compute_rij_wij_detail
        CC = compute_rij_wij_detail(
            lower_i,
            upper_i,
            self._data.data(),
            self._min_pairs)
        for cb_op in self.sym_ops:
            cb_op = sgtbx.change_of_basis_op(cb_op)
            indices_reindexed = cb_op.apply(self._data.indices())
            miller.map_to_asu(space_group_type, False, indices_reindexed)
            indices[cb_op.as_xyz()] = indices_reindexed
            CC.set_indices(cb_op, indices_reindexed)
        def call_cpp(i):
            t0 = time.time()
            rij_row, rij_col, rij_data, wij_row, wij_col, wij_data = [
                list(x) for x in CC.compute_one_row(self._lattices.size, i)
            ]
            rij = sparse.coo_matrix((rij_data, (rij_row, rij_col)), shape=(NN, NN))
            wij = sparse.coo_matrix((wij_data, (wij_row, wij_col)), shape=(NN, NN))

            print(time.time()-t0)
            return rij, wij

        if do_one_row_i is not None:
            results = [call_cpp(do_one_row_i)]
        else:
            results = easy_mp.parallel_map(
                call_cpp,
                range(n_lattices),
                processes=self._nproc)

        rij_matrix = None
        wij_matrix = None
        for (rij, wij) in results:
            if rij_matrix is None:
                rij_matrix = rij
                wij_matrix = wij
            else:
                rij_matrix += rij
                wij_matrix += wij
        self.rij_matrix = rij_matrix.todense().astype(np.float64)
        self.wij_matrix = wij_matrix.todense().astype(np.float64)
        return self.rij_matrix, self.wij_matrix


    def __init__(self):
      import pickle
      with open("big_data","rb") as F:
        self._lattices = pickle.load(F)
        self.sym_ops = pickle.load(F)
        self._weights = pickle.load(F)
        self._data = pickle.load(F)
        self._patterson_group = pickle.load(F)
        self._min_pairs = 3 # minimum number of mutual miller indices for a match

if __name__=="__main__":
  NPROC = 64
  REP = Reproducer()
  REP._nproc =  NPROC
  DEBUG = (len(sys.argv)!=1 and sys.argv[1]=='debug')
  t0 = time.time()
  do_one_row_i = 1 if DEBUG else None
  Rij_cpp, Wij_cpp = REP.compute_rij_wij_cplusplus(do_one_row_i=do_one_row_i)
  print("c++: summary time ", time.time()-t0)
  # set nproc to 64 for Gildea algorithm
  REP._nproc =  NPROC
  t0 = time.time()
  Rij, Wij = REP.compute_rij_wij_Gildea(lineprof=-1, do_one_row_i=do_one_row_i)
  print("py: summary time ", time.time()-t0)

  if DEBUG:
    Rij_old, Wij_old = REP.compute_rij_wij_Gildea_old(do_one_row_i=do_one_row_i)
  else:
    Rij_old, Wij_old = None, None

  for _ in range(10000):
      error_flag = False
      import random
      i = random.randint(0, Rij.size())
      ref = Rij[i]
      cpp = Rij_cpp[i]
      w_ref = Wij[i]
      w_cpp = Wij_cpp[i]
      if abs(ref-cpp)>0.001 or abs(w_ref-w_cpp)>0.001: error_flag = True
      if Rij_old is not None:
        ref_old = Rij_old[i]
        w_old = Wij_old[i]
        if abs(ref_old-ref)>0.001 or abs(w_old-w_ref)>0.001: error_flag = True

      if error_flag:
        print()
        print("{}: \t {:.3f} \t {:.3f} \t {:.3f}".format(i, ref, cpp, ref_old))
        print("{}: \t {} \t {} \t {}".format(i, w_ref, w_cpp, w_old))
        print
  if False:
        # debugging code useful for understanding the overall matrix structure
        from matplotlib import pyplot as plt
        plt.imshow(Rij.as_numpy_array())
        plt.show()
        plt.imshow(Wij.as_numpy_array(), vmin=0, vmax=20)
        plt.show()

  print("OK")
