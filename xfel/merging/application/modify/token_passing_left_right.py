from __future__ import division, print_function

def token_passing_left_right(value, comm): # comm = MPI communicator
  comm.barrier()
  mpi_rank = comm.rank
  mpi_size = comm.size

  src = mpi_rank - 1 if mpi_rank != 0 else mpi_size - 1
  dst = mpi_rank + 1 if mpi_rank != mpi_size - 1 else 0
  print (src,dst)
  if mpi_rank % 2 == 0:
    comm.send(value, dest=dst)
    m = comm.recv(source=src)
  else:
    m = comm.recv(source=src)
    comm.send(value, dest=dst)
  tokens = [m,value]
  comm.barrier()

  src = mpi_rank + 1 if mpi_rank != mpi_size - 1 else 0
  dst = mpi_rank - 1 if mpi_rank != 0 else mpi_size - 1

  if mpi_rank % 2 == 0:
    comm.send(value, dest=dst)
    m = comm.recv(source=src)
  else:
    m = comm.recv(source=src)
    comm.send(value, dest=dst)
  tokens.append(m)
  comm.barrier()
  return tokens

def get_root_communicator():
  from libtbx.mpi4py import MPI
  return MPI.COMM_WORLD

def tst_token_passing():
  comm = get_root_communicator()
  mpi_rank = comm.Get_rank()
  mpi_size = comm.Get_size()
  if mpi_size==1: print("OK"); return
  T = token_passing_left_right(value = 100+mpi_rank, comm = comm)
  print("On Rank %d the tokens are: %s"%(mpi_rank, str(T)))
  if mpi_rank==0: print("OK")

def choose_without_replace(available, unavailable, mt):
  while 1:
    idx = int(mt.random_double() * len(available))
    if available[idx] not in unavailable: break
  return available.pop(idx)

def construct_src_to_dst_plan(icount, tranch_size, comm, verbose=True):
      mpi_rank = comm.Get_rank()
      mpi_size = comm.Get_size()
      COMPOSITE_MULTIPLICITY = 3
      from scitbx.array_family import flex
      imean = flex.mean(icount.as_double()); isum = flex.sum(icount)
      trial_comm_size = max(1, int(isum * COMPOSITE_MULTIPLICITY / tranch_size))
      trial_comm_size = min(trial_comm_size, mpi_size - 1) # new comm must have fewer than mpi_helper.size ranks
      if trial_comm_size in [2,3]: trial_comm_size = 1 # consensus communicator can run on 1,4,5,..., not 2,3
      trial_comm_size = max(1, trial_comm_size) # cannot be 0

      if verbose: print("recommended communicator size:",trial_comm_size)
      trial = [[] for itrank in range(trial_comm_size)]
      srcrk = [[] for itrank in range(trial_comm_size)]
      todst = {} # plan of action: src rank values to dst keys
      mt = flex.mersenne_twister(seed=0)
      order = mt.random_permutation(len(icount))
      available_dst = list(range(trial_comm_size))
      for idx in range(len(icount)):
        item = icount[order[idx]] ; rk = idx
        if trial_comm_size>3:
          unavailable = []
          for irx in range(COMPOSITE_MULTIPLICITY):
            if len(available_dst)==0:  available_dst = list(range(trial_comm_size))
            idx_st = choose_without_replace(available_dst, unavailable, mt)
            unavailable.append(idx_st)
            trial[(idx_st)%trial_comm_size].append(item); srcrk[(idx_st)%trial_comm_size].append(rk)
            todst[idx_st] = todst.get(idx_st,[]); todst[idx_st].append(order[idx])
        else:
          trial[0].append(item)
          idx_st = 0
          todst[idx_st] = todst.get(idx_st,[]); todst[idx_st].append(order[idx])
      if verbose:
        for itranch, item in enumerate(trial):
          print("tranch %2d:"%itranch, item, flex.sum(flex.int(item)))
        print()
        for item in srcrk:
          print(", ".join(["%02d"%i for i in item]))
        print(todst)
      return todst

def construct_anchor_src_to_dst_plan(min_anchor, icount, tranch_size, comm, verbose=True):
      """The overall purpose is to recruit enough total experiments from several ranks to perform anchor alignment"""
      mpi_rank = comm.Get_rank()
      mpi_size = comm.Get_size()
      todst = {} # plan of action: src rank values to dst keys
      todst[0] = []
      anchor_sum = 0
      rank = 0
      rank_count = []
      while anchor_sum < min_anchor and rank < mpi_size:
        todst[0].append(rank)
        rank_count.append(icount[rank])
        anchor_sum += icount[rank]
        rank += 1
      if verbose:
        print("anchor tranch:",rank_count, anchor_sum)
        print(todst)
      if anchor_sum < min_anchor:
        raise Exception("Cannot align to anchor with %d total experiments, need at least %d"%(anchor_sum, min_anchor))
      return todst

def get_experiment_counts_by_comm_size():
  return {
    64:[52, 48, 48, 48, 47, 47, 46, 45, 44, 43, 43, 42, 43, 42, 42, 40, 40, 38, 38, 38, 38, 38,
        37, 37, 37, 37, 37, 35, 35, 34, 33, 32, 32, 32, 31, 31, 30, 30, 29, 30, 30, 30, 29, 29,
        28, 29, 29, 28, 28, 28, 28, 27, 28, 27, 28, 27, 28, 28, 28, 28, 27, 28, 28, 28],
    41:[66, 66, 66, 66, 66, 65, 65, 65, 65, 63, 63, 62, 61, 60, 60, 60, 60, 58, 56, 56, 56, 56,
        55, 54, 50, 50, 49, 48, 48, 46, 45, 45, 43, 43, 43, 43, 42, 42, 40, 40, 38],
    32:[84, 81, 81, 79, 77, 77, 76, 75, 75, 73, 72, 72, 72, 71, 71, 68, 68, 66, 66, 66, 66, 66, 66, 65, 65, 65, 64, 61, 61, 60, 59, 57],
    16:[146, 144, 142, 142, 141, 141, 141, 141, 140, 140, 137, 137, 135, 134, 133, 131],
    13:[179, 178, 173, 173, 171, 171, 171, 171, 170, 167, 167, 167, 167],
    10:[226, 226, 225, 224, 222, 222, 221, 220, 220, 219],
     5:[447, 446, 446, 444, 442],
     4:[558, 557, 555, 555],
     3:[743, 741, 741],
     2:[1113, 1112],
     1:[2225]}

def tst_tranch_size():
  from scitbx.array_family import flex
  comm = get_root_communicator()
  mpi_rank = comm.Get_rank()
  mpi_size = comm.Get_size()
  tranch_size = 600 # target size for the composite tranch in this example
  if mpi_rank == 0:
      icount = flex.int(get_experiment_counts_by_comm_size()[mpi_size])
      plan = construct_src_to_dst_plan(icount, tranch_size, comm)
  else:
      plan = 0
  plan = comm.bcast(plan, root = 0)
  # plan = dictionary < key = destination composite tranch; value = list of src ranks contributing the unique data >
  dst_offset = 1 if mpi_size>1 else 0
  tokens = apply_all_to_all(plan=plan, dst_offset=dst_offset,
                   value=get_experiment_counts_by_comm_size()[mpi_size][mpi_rank], comm = comm)
  print ("rank",mpi_rank,plan.get(mpi_rank-dst_offset, None), flex.sum(flex.int(tokens)))
  if mpi_rank == 0: print("OK")

def apply_all_to_all(plan, dst_offset, value, comm):
  comm.barrier()
  mpi_rank = comm.rank
  mpi_size = comm.size

  dst_values = [None]*mpi_size # defaults
  for dst_key in plan:
    if mpi_rank in plan[dst_key]:
      dst_values[dst_key + dst_offset] = value
  tokens = comm.alltoall(dst_values)
  while None in tokens:
    tokens.pop(tokens.index(None))
  comm.barrier()
  return tokens

if __name__=="__main__":
  # Usage: mpirun -n 64 libtbx.python $MODULES/cctbx_project/xfel/merging/application/modify/token_passing_left_right.py
  #tst_token_passing()
  tst_tranch_size()
