from __future__ import absolute_import, division, print_function
from xfel.merging.application.worker import worker
from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentList
from xfel.merging.application.reflection_table_utils import reflection_table_utils
import math

class load_balancer(worker):
  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(load_balancer, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

  def __repr__(self):
    return 'Balance input data load'

  def reflection_table_stub(self, reflections):
    '''Return an empty reflection table with the same format as the reflection table input to this class'''
    table = flex.reflection_table()
    for key in reflections:
      table[key] = type(reflections[key])()
    return table

  def divide_list_into_chunks(self, list_, n): # n - number of chunks
    return [list_[start::n] for start in range(n)]

  def run(self, experiments, reflections):
    self.logger.log("Rebalancing input load -- %s method..."%self.params.input.parallel_file_load.balance)
    if self.mpi_helper.rank == 0:
      self.logger.main_log("Rebalancing input load -- %s method..."%self.params.input.parallel_file_load.balance)
    if self.params.input.parallel_file_load.balance == "global2":
      # get status BEFORE balancing and print to main_log for EACH RANK if we have balance_verbose set
      from xfel.merging.application.utils.data_counter import data_counter
      if self.params.input.parallel_file_load.balance_verbose and self.mpi_helper.rank == 0:
        self.logger.main_log("Data distribution before load balancing:")
      expt_counts_by_rank, _, _ = data_counter(self.params).count_each(experiments, reflections, verbose=self.params.input.parallel_file_load.balance_verbose)
    else:
      expt_counts_by_rank = None

    if self.params.input.parallel_file_load.balance == "global1":
      new_experiments, new_reflections = self.distribute_over_ranks_shuffle(experiments, reflections, self.mpi_helper.comm, self.mpi_helper.size)
    elif self.params.input.parallel_file_load.balance == "global2":
      new_experiments, new_reflections = self.distribute_over_ranks_minimalist(experiments, reflections, self.mpi_helper.comm, self.mpi_helper.size, expt_counts_by_rank)
    elif self.params.input.parallel_file_load.balance == "per_node":
      mpi_color = int(self.mpi_helper.rank / self.params.input.parallel_file_load.ranks_per_node)
      mpi_new_rank = self.mpi_helper.rank % self.params.input.parallel_file_load.ranks_per_node
      mpi_split_comm = self.mpi_helper.comm.Split(mpi_color, mpi_new_rank)
      new_experiments, new_reflections = self.distribute_over_ranks(experiments, reflections, mpi_split_comm, self.params.input.parallel_file_load.ranks_per_node)

    if self.params.input.parallel_file_load.reset_experiment_id_column:
      self.logger.log('Starting id column reset')
      id_map = new_reflections.experiment_identifiers()
      reverse_map = {}
      for expt_id, experiment in enumerate(new_experiments):
        id_map[expt_id] = experiment.identifier
        reverse_map[experiment.identifier] = expt_id
      id_col = new_reflections['id']
      ident_col = new_reflections['exp_id']
      for i in range(len(new_reflections)):
        id_col[i] = reverse_map[ident_col[i]]
      self.logger.log('Column reset done')

    if self.params.input.parallel_file_load.balance == "global2":
      # Do we have any data?
      data_counter(self.params).count(new_experiments, new_reflections)

      # get status again AFTER balancing and report back number of experiments on EACH RANK in the main_log if balance_verbose is set
      if self.params.input.parallel_file_load.balance_verbose and self.mpi_helper.rank == 0:
        self.logger.main_log("Data distribution after load balancing:")
        data_counter(self.params).count_each(experiments, reflections, verbose=True)

    return new_experiments, new_reflections

  def distribute_over_ranks_shuffle(self, experiments, reflections, mpi_communicator, number_of_mpi_ranks):
    self.logger.log_step_time("LB_SPLIT_LIST")
    self.split_experiments = self.divide_list_into_chunks(experiments, number_of_mpi_ranks)
    self.logger.log_step_time("LB_SPLIT_LIST", True)

    # If some (but not all!) chunks are empty, we want those empty chunks to be randomly distributed.
    # Otherwise, after alltoall, high-index ranks will get no data.
    self.logger.log_step_time("LB_SHUFFLE")
    number_of_empty_chunks = [len(self.split_experiments[i]) for i in range(len(self.split_experiments))].count(0)
    if number_of_empty_chunks > 0 and len(experiments) != 0:
      import random
      #random.seed(8)
      random.shuffle(self.split_experiments)
    self.logger.log_step_time("LB_SHUFFLE", True)

    '''
    self.logger.log("Split experiment list into %d chunks"%len(self.split_experiments))
    for chunk in self.split_experiments:
      self.logger.log(len(chunk))
    '''

    # Distribute reflections over experiment chunks
    self.logger.log_step_time("LB_REF_DISTR")
    self.distribute_reflections_over_experiment_chunks_cpp(reflections)
    reflections.clear()
    del experiments
    self.logger.log_step_time("LB_REF_DISTR", True)

    # Run alltoall on experiments
    new_experiments = self.exchange_experiments_by_alltoall(mpi_communicator)

    # Run alltoall on reflections
    if self.params.input.parallel_file_load.balance_mpi_alltoall_slices == 1:
      new_reflections = self.exchange_reflections_by_alltoall(mpi_communicator)
    else:
      new_reflections = self.exchange_reflections_by_alltoall_sliced(mpi_communicator, self.params.input.parallel_file_load.balance_mpi_alltoall_slices)

    return new_experiments, new_reflections

  def exchange_experiments_by_alltoall(self, mpi_communicator):
    self.logger.log_step_time("LB_EXPTS_ALL_TO_ALL")
    new_split_experiments = mpi_communicator.alltoall(self.split_experiments)
    del self.split_experiments
    self.logger.log_step_time("LB_EXPTS_ALL_TO_ALL", True)

    self.logger.log_step_time("LB_EXPTS_CONSOLIDATE")
    self.logger.log("Consolidating experiments after all-to-all...")
    new_experiments = ExperimentList()
    for entry in new_split_experiments:
      new_experiments.extend(entry)
    del new_split_experiments
    self.logger.log_step_time("LB_EXPTS_CONSOLIDATE", True)

    return new_experiments

  def exchange_reflections_by_alltoall(self, mpi_communicator):
    ''' Run all-to-all and return a new reflection table'''
    self.logger.log_step_time("LB_REFLS_ALL_TO_ALL")
    new_split_reflections = mpi_communicator.alltoall(self.split_reflections)
    del self.split_reflections
    self.logger.log_step_time("LB_REFLS_ALL_TO_ALL", True)

    self.logger.log_step_time("LB_REFLS_CONSOLIDATE")
    self.logger.log("Consolidating reflections after all-to-all...")
    new_reflections = flex.reflection_table()
    for entry in new_split_reflections:
      new_reflections.extend(entry)
    del new_split_reflections
    self.logger.log_step_time("LB_REFLS_CONSOLIDATE", True)

    return new_reflections

  def exchange_reflections_by_alltoall_sliced(self, mpi_communicator, number_of_slices):
    '''Split each hkl chunk into N slices. This is needed to address the MPI alltoall memory problem'''
    result_reflections = flex.reflection_table() # the total reflection table, which this rank will receive after running all slices of alltoall
    list_of_sliced_reflection_chunks = [] # if the self.split_reflections list contains chunks: [A,B,C...], it will be sliced like: [[A1,A2,...,An], [B1,B2,...,Bn], [C1,C2,...,Cn], ...], where n is the number of chunk slices
    for i in range(len(self.split_reflections)):
      reflection_chunk_slices = []
      for chunk_slice in reflection_table_utils.get_next_reflection_table_slice(self.split_reflections[i], number_of_slices, self.reflection_table_stub):
        reflection_chunk_slices.append(chunk_slice)
      list_of_sliced_reflection_chunks.append(reflection_chunk_slices)

    for j in range(number_of_slices):
      reflection_chunks_for_alltoall = list()
      for i in range(len(self.split_reflections)):
        reflection_chunks_for_alltoall.append(list_of_sliced_reflection_chunks[i][j]) # [Aj,Bj,Cj...]

      self.logger.log_step_time("ALL-TO-ALL")
      received_reflection_chunks = mpi_communicator.alltoall(reflection_chunks_for_alltoall)
      self.logger.log("After all-to-all received %d reflection chunks" %len(received_reflection_chunks))
      self.logger.log_step_time("ALL-TO-ALL", True)

      self.logger.log_step_time("CONSOLIDATE")
      self.logger.log("Consolidating reflection tables...")
      for chunk in received_reflection_chunks:
        result_reflections.extend(chunk)
      self.logger.log_step_time("CONSOLIDATE", True)

    return result_reflections

  def distribute_reflections_over_experiment_chunks_python(self, reflections):
    self.split_reflections = []
    for i in range(len(self.split_experiments)):
      self.split_reflections.append(self.reflection_table_stub(reflections))

    for i in range(len(self.split_experiments)):
      for expt_idx, experiment in enumerate(self.split_experiments[i]):
        refls = reflections.select(reflections['exp_id'] == experiment.identifier)
        refls['id'] = flex.int(len(refls), expt_idx)
        self.split_reflections[i].extend(refls)

  def distribute_reflections_over_experiment_chunks_cpp(self, reflections):
    '''Distribute reflections over experiment chunks according to experiment identifiers, '''
    reflection_count = reflections.size()
    distributed_reflection_count = 0
    # initialize a list of reflection chunks
    self.split_reflections = []
    for i in range(len(self.split_experiments)):
      self.split_reflections.append(self.reflection_table_stub(reflections))

    if reflection_count > 0:
      # set up two lists to be passed to the C++ extension: experiment ids and chunk ids. It's basically a hash table to look up chunk ids by experiment identifier
      exp_id_list = flex.std_string()
      chunk_id_list = flex.int()
      for i in range(len(self.split_experiments)):
        for exp in self.split_experiments[i]:
          exp_id_list.append(exp.identifier)
          chunk_id_list.append(i)

      # distribute reflections over the experiment chunks using a C++ extension
      from xfel.merging import split_reflections_by_experiment_chunks_cpp
      split_reflections_by_experiment_chunks_cpp(reflections, exp_id_list, chunk_id_list, self.split_reflections)

      for ref_table in self.split_reflections:
        distributed_reflection_count += ref_table.size()

    self.logger.log("Distributed %d out of %d reflections"%(distributed_reflection_count, reflection_count))

  def distribute_over_ranks_minimalist(self, experiments, reflections, mpi_communicator, number_of_mpi_ranks, current_counts_by_rank):
    self.logger.log_step_time("LB_SPLIT_LIST")
    if self.mpi_helper.rank == 0:
      def first(lst, test):
        for i in range(len(lst)):
          if test(lst[i]):
            return i

      # quota: max number of counts that should end up on one rank once balanced
      # difference: how unbalanced each rank is currently
      # send_tuples: instructions for redistributing load, a list [L1,L2,L3...Li] where i is the rank
      # sending data and Li describes where to send it

      quota = int(math.ceil(sum(current_counts_by_rank)/len(current_counts_by_rank)))
      difference = [count - quota for count in current_counts_by_rank]
      send_tuples = [[] for i in range(len(current_counts_by_rank))]

      # algorithm (deterministic):
      # - for each rank, if we are overburdened, find the first rank that is underburdened
      # - send it as much of the load as it can take before hitting quota (denoted by a tuple of the target
      #   rank and the number of counts to send)
      # - repeat until no rank is overburdened (assert this to be the case)
      # - it may still be the case that ranks are underburdened unequally (e.g. [10,10,10,8,8])
      # - to address this, for each rank that is underburdened by more than 1, find the first rank at quota
      # - request one count from it
      # - assert all ranks are now within one count of each other

      for i in range(len(current_counts_by_rank)):
        while difference[i] > 0:
          j = first(difference, lambda count: count<0)
          send = min(difference[i], -1 * difference[j])
          send_tuples[i].append((j, send))
          difference[i] -= send
          difference[j] += send
      assert max(difference) == 0
      for i in range(len(current_counts_by_rank)):
        while difference[i] < -1:
          j = first(difference, lambda count: count==0)
          send_tuples[j].append((i, 1))
          difference[i] += 1
          difference[j] -= 1
      assert max(difference) == 0
      assert min(difference) >= -1

    else:
      send_tuples = None # not sure if we need this?

    # broadcast instructions
    send_tuples = mpi_communicator.bcast(send_tuples, root=0)

    self.logger.log_step_time("LB_SPLIT_LIST", True)
    self.logger.log_step_time("LB_EXPTS_AND_REFLS_ALLTOALL")

    # carry out load balancing with all-to-all mpi communication
    send_instructions = send_tuples[self.mpi_helper.rank]

    # pare down balanced_experiments and balanced_reflections as we separate off what to send out
    send_data = [(None, None) for j in range(len(send_tuples))]
    recv_data = [(None, None) for j in range(len(send_tuples))]
    for (j, count) in send_instructions:
      send_expt_j = experiments[-count:]
      experiments = experiments[:-count]
      send_refl_j = self.reflection_table_stub(reflections)
      for k, e in enumerate(send_expt_j):
        r = reflections.select(reflections['exp_id'] == e.identifier) # select matching reflections to send
        r['id'] = flex.int(len(r), k)
        send_refl_j.extend(r)
        reflections = reflections.select(reflections['exp_id'] != e.identifier) # remove from this rank's reflections
      send_data[j] = (send_expt_j, send_refl_j)
    recv_data = mpi_communicator.alltoall(send_data)

    # tack on only what was targeted to be received by the current rank
    for (received_expt_i, received_refl_i) in recv_data:
      if received_expt_i is None: continue
      current_ids = set([e.identifier for e in experiments])
      recv_ids = set([e.identifier for e in received_expt_i])
      assert current_ids.isdisjoint(recv_ids)
      experiments.extend(received_expt_i)
      reflections.extend(received_refl_i)

    self.logger.log_step_time("LB_EXPTS_AND_REFLS_ALLTOALL", True)

    return experiments, reflections

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(load_balancer)
