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

  def run(self, experiments, reflections):
    self.logger.log("Rebalancing input load -- %s method..."%self.params.input.parallel_file_load.balance)
    if self.mpi_helper.rank == 0:
      self.logger.main_log("Rebalancing input load -- %s method..."%self.params.input.parallel_file_load.balance)
    # get status BEFORE balancing and print to main_log for EACH RANK if we have balance_verbose set
    from xfel.merging.application.utils.data_counter import data_counter
    if self.params.input.parallel_file_load.balance_verbose and self.mpi_helper.rank == 0:
      self.logger.main_log("Data distribution before load balancing:")
    expt_counts_by_rank, _, _ = data_counter(self.params).count_each(experiments, reflections, verbose=self.params.input.parallel_file_load.balance_verbose)

    if self.params.input.parallel_file_load.balance == "global":
      new_experiments, new_reflections = self.distribute_over_ranks(experiments, reflections, self.mpi_helper.comm, self.mpi_helper.size, expt_counts_by_rank)
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

    # Do we have any data?
    data_counter(self.params).count(new_experiments, new_reflections)

    # get status again AFTER balancing and report back number of experiments on EACH RANK in the main_log if balance_verbose is set
    if self.params.input.parallel_file_load.balance_verbose and self.mpi_helper.rank == 0:
      self.logger.main_log("Data distribution after load balancing:")
      data_counter(self.params).count_each(experiments, reflections, verbose=True)

    return new_experiments, new_reflections

  def distribute_over_ranks(self, experiments, reflections, mpi_communicator, number_of_mpi_ranks, current_counts_by_rank):
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
    self.logger.log_step_time("LB_EXPTS_AND_REFLS_POINT_TO_POINT")

    # carry out load balancing with point-to-point mpi communication
    send_instructions = send_tuples[self.mpi_helper.rank]

    # pare down balanced_experiments and balanced_reflections as we separate off what to send out
    for (j, count) in send_instructions:
      send_experiments = experiments[-count:]
      mpi_communicator.send(send_experiments, dest=j, tag=0)
      experiments = experiments[:-count]
      send_reflections = self.reflection_table_stub(reflections)
      for k, e in enumerate(send_experiments):
        r = reflections.select(reflections['exp_id'] == e.identifier) # select matching reflections to send
        r['id'] = flex.int(len(r), k)
        send_reflections.extend(r)
        reflections = reflections.select(reflections['exp_id'] != e.identifier) # remove from this rank's reflections
      mpi_communicator.send(send_reflections, dest=j, tag=1) # 0 for expts, 1 for refls

    # tack on only what was targeted to be received by the current rank
    identifiers_taken = set([e.identifier for e in experiments])
    for i, receive_instructions in enumerate(send_tuples):
      for (j, count) in receive_instructions:
        if j == self.mpi_helper.rank:
          received_experiments = mpi_communicator.recv(source=i, tag=0)
          received_reflections = mpi_communicator.recv(source=i, tag=1)
          # note: if this passes assertions we can probably switch from appending each one
          # to extending by the received expts list
          for e in received_experiments:
            assert e.identifier not in identifiers_taken
            identifiers_taken.add(e.identifier)
            experiments.append(e)
            r = received_reflections.select(received_reflections['exp_id'] == e.identifier)
            reflections.extend(r)
            received_reflections = received_reflections.select(received_reflections['exp_id'] != e.identifier)
          assert len(received_reflections) == 0
    # may want to add some more assertions about number of things passed around

    self.logger.log_step_time("LB_EXPTS_AND_REFLS_POINT_TO_POINT", True)

    return experiments, reflections

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(load_balancer)
