from __future__ import absolute_import, division, print_function
from xfel.merging.application.worker import worker
from dials.array_family import flex
from six.moves import cStringIO as StringIO
from six.moves import zip

class unit_cell_statistics(worker):

  def __repr__(self):
    return 'Unit cell statistics'

  def run(self, experiments, reflections):
    self.logger.log_step_time("UNIT_CELL_STATISTICS")

    if self.mpi_helper.rank == 0:
      self.logger.main_log("Unit cell length distribution:")

    self.show_histograms(experiments)

    self.logger.log_step_time("UNIT_CELL_STATISTICS", True)

    return experiments, reflections

  def show_histograms(self, experiments):
    '''Output histograms of unit cell parameters'''
    n_slots = 20

    total_edge0 = []
    total_edge1 = []
    total_edge2 = []
    edge0 = []
    edge1 = []
    edge2 = []
    for experiment in experiments:
      exp_unit_cell = experiment.crystal.get_unit_cell()
      (a,b,c,alpha,beta,gamma) = exp_unit_cell.parameters()

      edge0.append(a)
      edge1.append(b)
      edge2.append(c)

    comm = self.mpi_helper.comm
    MPI = self.mpi_helper.MPI

    total_edge0  = comm.reduce(edge0, MPI.SUM, 0)
    total_edge1  = comm.reduce(edge1, MPI.SUM, 0)
    total_edge2  = comm.reduce(edge2, MPI.SUM, 0)

    # rank 0: log data statistics
    if self.mpi_helper.rank == 0:

      (a0,b0,c0,alpha0,beta0,gamma0) = self.params.scaling.unit_cell.parameters()

      self.logger.main_log("")

      labels = ["a","b","c"]

      ref_edges = [a0, b0, c0]

      def _show_each(edges):

        for edge, ref_edge, label in zip(edges, ref_edges, labels):
          h = flex.histogram(edge, n_slots=n_slots)

          smin, smax = flex.min(edge), flex.max(edge)
          stats = flex.mean_and_variance(edge)

          self.logger.main_log("  %s edge"%label)
          self.logger.main_log("     range:     %6.2f - %.2f"%(smin, smax))
          self.logger.main_log("     mean:      %6.2f +/- %6.2f on N = %d" %(stats.mean(), stats.unweighted_sample_standard_deviation(), edge.size()))
          self.logger.main_log("     reference: %6.2f"%ref_edge)

          out = StringIO()
          h.show(f=out, prefix="    ", format_cutoffs="%6.2f")
          self.logger.main_log(out.getvalue() + '\n')

      edges = [flex.double(total_edge0), flex.double(total_edge1), flex.double(total_edge2)]
      _show_each(edges)

      '''
      # for LD91 test
      from matplotlib import pyplot as plt
      left = 77
      bottom = 34
      width = 85 - 77
      height = 44 - 34
      rect = [left, bottom, width, height]

      plt.xlim(77,85)
      plt.ylim(34,44)
      plt.plot(total_edge0, total_edge2, ".")
      plt.show()
      '''

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(unit_cell_statistics)
