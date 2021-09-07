from __future__ import division
from xfel.merging.application.worker import worker
from dials.array_family import flex


class annulus_statistics(worker):

  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(annulus_statistics, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

  def __repr__(self):
    return "Calculate reflection statistics in a resolution annulus"

  def run(self, experiments, reflections):
    uc = self.params.scaling.unit_cell
    d_min = 2.1
    d_max = 2.5

    all_refl_in_annulus = flex.reflection_table()
#    d_min = self.params.statistics.annulus.d_min
#    d_max = self.params.statistics.annulus.d_max
    for expt in experiments:
      self.logger.log('Annulus statistics: ', expt.identifier)
      exp_id = expt.identifier
      refl = reflections.select(reflections['exp_id']==exp_id)
      dspacings = uc.d(refl['miller_index'])
      gt = dspacings > d_min
      lt = dspacings < d_max
      refl = refl.select(gt & lt)
      all_refl_in_annulus.extend(refl)
      self.logger.log('{} spots in annulus'.format(refl.size()))

    all_all_refl_in_annulus = self.mpi_helper.comm.gather(
        all_refl_in_annulus, root=0
    )
    final_refl_in_annulus = flex.reflection_table()
    if self.mpi_helper.rank==0:
      for table in all_all_refl_in_annulus:
        final_refl_in_annulus.extend(table)
      counts = {}
      for m_i in final_refl_in_annulus['miller_index']:
        count = counts.setdefault(m_i, 0)
        counts[m_i] = count + 1
      print('total shoeboxes: {}'.format(sum(counts.values())))
      print('unique miller indices: {}'.format(len(counts.keys())))
      print('average multiplicity: {:.3f}'.format(sum(counts.values())/len(counts.values())))



    return experiments, reflections
      
      
    #import IPython;IPython.embed()
