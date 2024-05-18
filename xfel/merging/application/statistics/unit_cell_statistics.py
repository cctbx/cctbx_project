from __future__ import absolute_import, division, print_function
from xfel.merging.application.worker import worker
from dials.array_family import flex
from six.moves import cStringIO as StringIO
from six.moves import zip
from cctbx import uctbx

histogram_slots = 20

class unit_cell_distribution(object):
  """Container for collecting unit cell statistics"""
  # TODO make this more general - currently assumes that angles are fixed,
  # which is true for the systems studied so far
  def __init__(self, reference_unit_cell, logger, mpi_helper, precision=None):
    self.reference_unit_cell = reference_unit_cell
    self.logger = logger
    self.mpi_helper = mpi_helper
    if precision is None: self.uc_precision = 2
    else: self.uc_precision = precision

    # this rank cell values
    self.uc_a_values = flex.double()
    self.uc_b_values = flex.double()
    self.uc_c_values = flex.double()

    self.uc_alpha_values  = flex.double()
    self.uc_beta_values   = flex.double()
    self.uc_gamma_values  = flex.double()

    # all ranks cell values
    self.all_uc_a_values = flex.double()
    self.all_uc_b_values = flex.double()
    self.all_uc_c_values = flex.double()

    self.all_uc_alpha_values  = flex.double()
    self.all_uc_beta_values   = flex.double()
    self.all_uc_gama_values   = flex.double()

  def add_cell(self, unit_cell):
    if unit_cell is None:
      return
    (a,b,c,alpha,beta,gamma) = unit_cell.parameters()
    self.uc_a_values.append(a)
    self.uc_b_values.append(b)
    self.uc_c_values.append(c)

    self.uc_alpha_values.append(alpha)
    self.uc_beta_values.append(beta)
    self.uc_gamma_values.append(gamma)

  def collect_from_all_ranks(self):
    self.all_uc_a_values = self.mpi_helper.aggregate_flex(self.uc_a_values, flex.double)
    self.all_uc_b_values = self.mpi_helper.aggregate_flex(self.uc_b_values, flex.double)
    self.all_uc_c_values = self.mpi_helper.aggregate_flex(self.uc_c_values, flex.double)

    self.all_uc_alpha_values  = self.mpi_helper.aggregate_flex(self.uc_alpha_values, flex.double)
    self.all_uc_beta_values   = self.mpi_helper.aggregate_flex(self.uc_beta_values, flex.double)
    self.all_uc_gamma_values  = self.mpi_helper.aggregate_flex(self.uc_gamma_values, flex.double)

  def is_valid(self):
    return len(self.all_uc_a_values) > 0 and len(self.all_uc_b_values) > 0 and len(self.all_uc_c_values) > 0 and \
           len(self.all_uc_alpha_values) > 0 and len(self.all_uc_beta_values) > 0 and len(self.all_uc_gamma_values) > 0

  def show_histograms(self, n_slots=histogram_slots):
    assert self.mpi_helper.rank == 0

    if self.reference_unit_cell is None:
      a0 = b0 = c0 = alpha0 = beta0 = gamma0 = None
    else:
      a0,b0,c0,alpha0,beta0,gamma0 = self.reference_unit_cell.parameters()

    self.logger.main_log("")

    labels = ["a","b","c"]

    ref_edges = [a0,b0,c0]

    def _show_each(edges):
      for edge, ref_edge, label in zip(edges, ref_edges, labels):
        h = flex.histogram(edge, n_slots=n_slots)
        smin, smax = flex.min(edge), flex.max(edge)
        stats = flex.mean_and_variance(edge)

        self.logger.main_log("  %s edge"%label)
        range_template = "     range:     %6.{0}f - %.{0}f".format(self.uc_precision)
        self.logger.main_log(range_template %(smin, smax))
        mean_template = "     mean:      %6.{0}f +/- %6.{0}f on N = %d".format(self.uc_precision)
        self.logger.main_log(mean_template %(stats.mean(), stats.unweighted_sample_standard_deviation(), edge.size()))
        if ref_edge is not None:
          ref_template = "     reference: %6.{}f".format(self.uc_precision)
          self.logger.main_log(ref_template %ref_edge)

        out = StringIO()
        h.show(f=out, prefix="    ", format_cutoffs="%6.{}f".format(self.uc_precision))
        self.logger.main_log(out.getvalue() + '\n')

    edges = [self.all_uc_a_values, self.all_uc_b_values, self.all_uc_c_values]
    _show_each(edges)

  def get_average_cell(self):
    a = flex.mean(self.all_uc_a_values)
    b = flex.mean(self.all_uc_b_values)
    c = flex.mean(self.all_uc_c_values)

    alpha = flex.mean(self.all_uc_alpha_values)
    beta  = flex.mean(self.all_uc_beta_values)
    gamma = flex.mean(self.all_uc_gamma_values)

    return uctbx.unit_cell([a,b,c,alpha,beta,gamma])

class unit_cell_statistics(worker):

  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(unit_cell_statistics, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

  def __repr__(self):
    return 'Unit cell statistics'

  def run(self, experiments, reflections):
    self.logger.log_step_time("UNIT_CELL_STATISTICS")
    ucd = unit_cell_distribution(
        self.params.scaling.unit_cell,
        self.logger,
        self.mpi_helper,
        precision=self.params.statistics.uc_precision
        )
    for experiment in experiments:
      ucd.add_cell(experiment.crystal.get_unit_cell())
    ucd.collect_from_all_ranks()

    average_unit_cell = None
    if self.mpi_helper.rank == 0:
      if ucd.is_valid():
        ucd.show_histograms()
        average_unit_cell = ucd.get_average_cell()

    self.logger.log_step_time("BROADCAST_UNIT_CELL")
    average_unit_cell = self.mpi_helper.comm.bcast(average_unit_cell, root = 0)
    self.logger.log_step_time("BROADCAST_UNIT_CELL", True)

    # save the average unit cell to the phil parameters
    if self.mpi_helper.rank == 0:
      self.logger.main_log("Average unit_cell %s is saved to phil parameters"%str(average_unit_cell))
    if not 'average_unit_cell' in (self.params.statistics).__dict__:
      self.params.statistics.__inject__('average_unit_cell', average_unit_cell)
    else:
      self.params.statistics.__setattr__('average_unit_cell', average_unit_cell)

    self.logger.log_step_time("UNIT_CELL_STATISTICS", True)

    return experiments, reflections

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(unit_cell_statistics)
