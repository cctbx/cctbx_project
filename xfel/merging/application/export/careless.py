from __future__ import absolute_import, division, print_function
from xfel.merging.application.worker import worker
from dials.array_family import flex
from cctbx import sgtbx
import numpy as np
import gemmi
import copy
import pandas as pd
import reciprocalspaceship as rs

class export_careless_mtz(worker):
  '''A class that condenses all expt + refl for output to a single batch MTZ.
     Inspired by https://github.com/rs-station/careless/blob/main/scripts/stills2mtz
  '''

  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(export_careless_mtz, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

  def __repr__(self):
    return 'Output batch MTZ for input to careless'

  def gather_N_by_3_array(self,flexdata, dtype='d', nn=3): # vector length nn=3(default), 2 or 1
    # takes flex array, returns numpy array concatenated fom all ranks
    # dtype d(default,float64) or i(int32)
    comm = self.mpi_helper.comm
    MPI = self.mpi_helper.MPI
    rank = comm.Get_rank()
    size = comm.Get_size()

    sendbuf = flexdata.as_numpy_array()
    # Each rank has a different number of rows, k

    # Step 1: Gather all the k values to rank 0
    sendcount = sendbuf.size  # k * 3 (total number of elements)
    counts = comm.gather(sendcount, root=0)

    # Step 2: Build the receive buffer and displacements on rank 0
    if rank == 0:
      counts = np.array(counts, dtype='i')
      displacements = np.zeros(size, dtype='i')
      displacements[1:] = np.cumsum(counts[:-1])
      total_elements = np.sum(counts)
      recvbuf = np.empty(total_elements, dtype=dtype)
    else:
      counts = None
      displacements = None
      recvbuf = None

    # Step 3: Gatherv
    comm.Gatherv(sendbuf, [recvbuf, counts, displacements, dict(i=MPI.INT,d=MPI.DOUBLE)[dtype]], root=0)

    # Step 4: Reshape on rank 0
    if rank == 0:
      result = recvbuf.reshape(-1, nn)
      # print(f"Final shape: {result.shape}")  # (sum_of_k, nn) # like (27168394, 3)
      return result
    else: return None

  def run(self, all_experiments, all_reflections):
    if self.mpi_helper.rank == 0:
      self.logger.log("Writing out a batch MTZ file")

    comm = self.mpi_helper.comm
    MPI = self.mpi_helper.MPI
    rank = comm.Get_rank()
    size = comm.Get_size()

    # do some preliminary stuff to get all the desired data in to numpy arrays
    # BEGIN careless stills2mtz setup
    copy_reflections = copy.deepcopy(all_reflections)
    idx = flex.size_t(np.array(all_reflections['id']))

    #This makes some heavy assumptions about file name formats. You might want to remove this block to
    #suit your own use case
    images = [e.imageset.get_image_identifier(0) for e in all_experiments]

    copy_reflections["A_matrix"] = flex.mat3_double( [C.get_A() for C in all_experiments.crystals()] ).select(idx)
    copy_reflections["Ainv_matrix"] = flex.mat3_double( [C.get_A_inverse_as_sqr() for C in all_experiments.crystals()] ).select(idx)
    copy_reflections["B_matrix"] = flex.mat3_double( [C.get_B() for C in all_experiments.crystals()] ).select(idx)
    copy_reflections["s0_vec"] = flex.vec3_double( [e.beam.get_s0() for e in all_experiments] ).select(idx)
    copy_reflections["wavelength"] = flex.double( [e.beam.get_wavelength() for e in all_experiments] ).select(idx)
    h = copy_reflections["miller_index"].as_vec3_double()
    x = copy_reflections["A_matrix"] * h
    Svec = x + copy_reflections["s0_vec"]
    copy_reflections["Rh"] = Svec.norms() - (1./copy_reflections["wavelength"])
    copy_reflections["miller_index_obs"] = copy_reflections['Ainv_matrix'] * (copy_reflections['s1'] - copy_reflections['s0_vec'])
    copy_reflections["cartesian_fixed_obs"] = copy_reflections['B_matrix'] * copy_reflections['miller_index_obs']
    copy_reflections["cartesian_fixed"] = copy_reflections['B_matrix'] * h
    copy_reflections.compute_miller_indices_in_asu(all_experiments)

    #Get a gemmi cell
    cell = np.zeros(6)
    for crystal in all_experiments.crystals():
        cell += np.array(crystal.get_unit_cell().parameters())/len(all_experiments.crystals())
    cell = gemmi.UnitCell(*cell)

    sginfo = all_experiments.crystals()[0].get_space_group().info()
    symbol = sgtbx.space_group_symbols(sginfo.symbol_and_number().split('(')[0]) #<--- this cannot be the 'correct' way to do this
    spacegroup = gemmi.SpaceGroup(symbol.universal_hermann_mauguin())
    #print(
    #"%.5f"%(cell.a),"%.5f"%(cell.b),"%.5f"%(cell.c),"%.5f"%(cell.alpha),"%.5f"%(cell.beta),"%.5f"%(cell.gamma),
    #symbol.universal_hermann_mauguin().replace(" ", ""))
    # END careless stills2mtz setup

    # BEGIN fix id column
    send_count_expt = len(all_experiments)
    send_count_refl = copy_reflections.size()
    counts_expt = comm.allgather(send_count_expt)
    counts_refl = comm.allgather(send_count_refl)
    counts_expt_base = 0
    for incr in range(1, rank+1):
      counts_expt_base += counts_expt[incr-1]
    copy_reflections['id']+= counts_expt_base # computes globally unique id, not just with-rank id
    # END fix id column

    # weighted average all unit cells together on all ranks
    all_params = comm.allgather(cell.parameters)
    sum_wt_times_cell = np.zeros(6)
    for item in range(len(all_params)):
      sum_wt_times_cell += np.array(all_params[item]) * counts_expt[item]
    average_cell = sum_wt_times_cell / sum(counts_expt)
    cell = gemmi.UnitCell(*average_cell)
    comm.barrier()

    # average the wavelength all ranks
    all_wavelengths = flex.double([B.get_wavelength() for B in all_experiments.beams()])
    one_rank_sum_wavelen = flex.sum(all_wavelengths)
    all_ranks_sum_wavelen = comm.allgather(one_rank_sum_wavelen)
    wavelength = sum(all_ranks_sum_wavelen) / sum(counts_expt)

    # BEGIN Gather data from all ranks to create dataframe
    h_global = self.gather_N_by_3_array(h)
    comm.barrier()
    cartesian_fixed_obs_global = self.gather_N_by_3_array(copy_reflections['cartesian_fixed_obs'])
    cartesian_fixed_global = self.gather_N_by_3_array(copy_reflections['cartesian_fixed'])
    cal_global = self.gather_N_by_3_array(copy_reflections['xyzcal.px'])
    obs_global = self.gather_N_by_3_array(copy_reflections['xyzobs.px.value'])
    obs_variance_global = self.gather_N_by_3_array(copy_reflections['xyzobs.px.variance'])
    batch_global = self.gather_N_by_3_array(copy_reflections['id'], dtype='i', nn=1)
    offset_global = self.gather_N_by_3_array(copy_reflections['Rh'], nn=1)
    intensity_global = self.gather_N_by_3_array(copy_reflections['intensity.sum.value'], nn=1)
    variance_global = self.gather_N_by_3_array(copy_reflections['intensity.sum.variance'], nn=1)

    comm.barrier()
    # END Gather data from all ranks

    if rank==0:
      # Final output; use pandas dataframe
      df = pd.DataFrame({
      'H' : h_global[:,0].astype(np.int32),
      'K' : h_global[:,1].astype(np.int32),
      'L' : h_global[:,2].astype(np.int32),
      'BATCH' : batch_global[:,0],
      'cartesian_fixed_obs_x' : cartesian_fixed_obs_global[:,0],
      'cartesian_fixed_obs_y' : cartesian_fixed_obs_global[:,1],
      'cartesian_fixed_obs_z' : cartesian_fixed_obs_global[:,2],
      'cartesian_fixed_x'     : cartesian_fixed_global[:,0],
      'cartesian_fixed_y'     : cartesian_fixed_global[:,1],
      'cartesian_fixed_z'     : cartesian_fixed_global[:,2],
      'ewald_offset' : offset_global[:,0],
      'I' : intensity_global[:,0],
      'SigI' : variance_global[:,0]**0.5,
      'xcal' : cal_global[:,0],
      'ycal' : cal_global[:,1],
      'xobs' : obs_global[:,0],
      'yobs' : obs_global[:,1],
      'sigxobs' : obs_variance_global[:,0]**0.5,
      'sigyobs' : obs_variance_global[:,1]**0.5,
      })

      df['cartesian_delta_x'] = df['cartesian_fixed_obs_x'] - df['cartesian_fixed_x']
      df['cartesian_delta_y'] = df['cartesian_fixed_obs_y'] - df['cartesian_fixed_y']
      df['cartesian_delta_z'] = df['cartesian_fixed_obs_z'] - df['cartesian_fixed_z']

      """Guessing this will never work because we are only recording the ASU Miller index, not
the original index.  Therefore we can never know the true Ewald offset, so we can't
really teach the NN about partiality.  Real BATCH MTZ files have M/ISYM that encodes which
symmetry mate the original index uses.  OK, it seems we do have M_ISYM, but I don't know
how it is calculated (somehow with reciprocalspaceship?).  We have a column for ewald_offset
but it is set to zero so no help.  Also as Aaron pointed out, the cartesian xyz are meaningless
due to having no panel number.  xcal,ycal,xobs, yobs all meaningless.
"""

      data = rs.DataSet(df, cell=cell, spacegroup=spacegroup)
      typedefs = {
        'H' : 'H',
        'K' : 'H',
        'L' : 'H',
        'I' : 'J',
        'SigI' : 'Q',
        'BATCH' : 'B',
      }
      for k in data:
        if k in typedefs:
          data[k] = data[k].astype(typedefs[k])
        else:
          data[k] = data[k].astype('R')

      data.set_index(['H', 'K', 'L'], inplace=True)
      for k in data:
        if 'Unnamed' in k:
          del(data[k])

      # Use gemmi interface in order to set wavelength value, bypassing reciprocalspaceship
      # data.write_mtz(self.params.export.output_label+".mtz", skip_problem_mtztypes=True)
      from reciprocalspaceship.decorators import range_indexed
      from reciprocalspaceship import io
      @range_indexed
      def bypass_data_write_mtz(mydataset):
        mtz = io.mtz.to_gemmi(
        mydataset, skip_problem_mtztypes=True, project_name="reciprocalspaceship",
        crystal_name="reciprocalspaceship", dataset_name="reciprocalspaceship"
        )
        mtz.datasets[0].wavelength=wavelength
        mtz.write_to_file("step_%d_"%(self.logger.worker_step)+self.params.export.output_label+".mtz")
      bypass_data_write_mtz(data)

    if rank==0:
      self.logger.main_log(str(counts_expt))
      self.logger.main_log("Sums to %d"%(sum(counts_expt)))
      self.logger.main_log(str(counts_refl))
      self.logger.main_log("Sums to %d"%(sum(counts_refl)))
      self.logger.main_log("Data to file step_%d_%s.mtz wavelength=%.5fÅ"%(
                           self.logger.worker_step,self.params.export.output_label, wavelength))

    # what is left?
    from xfel.merging.application.utils.data_counter import data_counter
    data_counter(self.params).count(all_experiments, all_reflections)

    for iexpt, cexpt in enumerate(all_experiments):
      desc = "SPECT step_%d_%s %s imageindex %d iexpt %d of %d on rank %d"%(
                           self.logger.worker_step, self.params.export.output_label,
                           cexpt.imageset.paths()[0], cexpt.imageset.indices()[0], iexpt, len(all_experiments), rank)
      self.logger.log(desc)
    return all_experiments, all_reflections

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  # XXX cannot be tested as currently implemented # exercise_worker(export_careless_mtz)
