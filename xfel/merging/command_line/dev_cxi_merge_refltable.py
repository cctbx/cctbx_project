# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# LIBTBX_SET_DISPATCHER_NAME dev.cxi.merge_refltable
#
# $Id$
from __future__ import absolute_import, division, print_function
from six.moves import range

from xfel.command_line.cxi_merge import scaling_manager
from dials.array_family import flex
from xfel.merging.command_line.dev_cxi_merge import Script
from xfel.cxi.merging_utils import null_data
from libtbx import Auto

from xfel.merging.algorithms.error_model.sdfac_refine import sdfac_refine_refltable
import xfel.merging.algorithms.error_model
xfel.merging.algorithms.error_model.sdfac_refine.sdfac_refine = sdfac_refine_refltable

from xfel.merging.algorithms.error_model.errors_from_residuals import errors_from_residuals_refltable
import xfel.merging.algorithms.error_model
xfel.merging.algorithms.error_model.errors_from_residuals.errors_from_residuals = errors_from_residuals_refltable

def merging_reflection_table():
  table = flex.reflection_table()
  table['miller_index'] = flex.miller_index()
  table['miller_index_original'] = flex.miller_index()
  table['scaled_intensity'] = flex.double()
  table['isigi'] = flex.double()
  table['slope'] = flex.double()
  table['miller_id'] = flex.size_t()
  table['crystal_id'] = flex.size_t()
  table['iobs'] = flex.double()
  return table

def merging_crystal_table(postrefinement_columns = False):
  table = flex.reflection_table()
  table['u_matrix'] = flex.mat3_double()
  table['b_matrix'] = flex.mat3_double()
  table['wavelength'] = flex.double()
  if postrefinement_columns:
    table['G'] = flex.double()
    table['B'] = flex.double()
    table['RS'] = flex.double()
    table['thetax'] = flex.double()
    table['thetay'] = flex.double()
  table['ETA'] = flex.double()
  table['DEFF'] = flex.double()
  table['n_refl'] = flex.size_t()
  return table

class refltable_scaling_manager(scaling_manager):
  def initialize (self) :
    super(refltable_scaling_manager, self).initialize()
    self.ISIGI = merging_reflection_table()
    self.crystal_table = merging_crystal_table(self.params.postrefinement.enable)

  @staticmethod
  def single_reflection_histograms(obs, ISIGI):
    # Per-bin sum of I and I/sig(I) for each observation.
    import numpy as np
    for i in obs.binner().array_indices(i_bin) :
      index = obs.indices()[i]
      refls = ISIGI.select(ISIGI['miller_index'] == index)
      if (len(refls) > 0) :
        # Compute m, the "merged" intensity, as the average intensity
        # of all observations of the reflection with the given index.
        N = 0
        m = 0
        for j in range(len(refls)):
          N += 1
          m += refls['isigi'][j]

          print("Miller %20s n-obs=%4d  sum-I=%10.0f"%(index, N, m))
          plot_n_bins = N//10
          hist,bins = np.histogram(refls['intensity'],bins=25)
          width = 0.7*(bins[1]-bins[0])
          center = (bins[:-1]+bins[1:])/2
          import matplotlib.pyplot as plt
          plt.bar(center, hist, align="center", width=width)
          plt.show()

  def scale_frame_detail(self, result, file_name, db_mgr, out):
    data = super(refltable_scaling_manager, self).scale_frame_detail(result, file_name, db_mgr, out)
    if isinstance(data, null_data):
      return data
    reflections = merging_reflection_table()
    crystal_id = len(self.crystal_table)
    for i, hkl in enumerate(self.miller_set.indices()):
      if hkl not in data.ISIGI: continue
      for j, refl in enumerate(data.ISIGI[hkl]):
        reflections.append({'miller_index':hkl,
                            'scaled_intensity': refl[0],
                            'isigi': refl[1],
                            'slope': refl[2],
                            'miller_id':i,
                            'crystal_id':crystal_id,
                            'iobs':data.extra_stuff[hkl][0][j],
                            'miller_index_original':data.extra_stuff[hkl][1][j]})

    data.ISIGI = reflections

    from scitbx.matrix import sqr
    ori = data.current_orientation
    a_matrix = sqr(ori.reciprocal_matrix())
    b_matrix = sqr(ori.unit_cell().fractionalization_matrix()).transpose()
    u_matrix = a_matrix * b_matrix.inverse()

    crystal_d = {'b_matrix':b_matrix.elems,
                 'u_matrix':u_matrix.elems,
                 'wavelength':data.wavelength,
                 'n_refl':len(reflections)}
    if self.params.postrefinement.enable:
      crystal_d['G'] = self.postrefinement_params.G
      crystal_d['B'] = self.postrefinement_params.BFACTOR
      crystal_d['thetax'] = self.postrefinement_params.thetax
      crystal_d['thetay'] = self.postrefinement_params.thetay
      if self.params.postrefinement.algorithm in ['eta_deff']:
        crystal_d['ETA'] = self.postrefinement_params.ETA
        crystal_d['DEFF'] = self.postrefinement_params.DEFF
      else:
        crystal_d['RS'] = self.postrefinement_params.RS
    if not self.params.postrefinement.enable or self.params.postrefinement.algorithm not in ['eta_deff']:
      crystal_d['ETA'] = result['ML_half_mosaicity_deg'][0]
      crystal_d['DEFF'] = result['ML_domain_size_ang'][0]
    self.crystal_table.append(crystal_d)

    return data

  def _scale_all_parallel (self, file_names, db_mgr) :
    import multiprocessing
    import libtbx.introspection

    nproc = self.params.nproc
    if (nproc is None) or (nproc is Auto) :
      nproc = libtbx.introspection.number_of_processors()

    # Input files are supplied to the scaling processes on demand by
    # means of a queue.
    #
    # XXX The input queue may need to either allow non-blocking
    # put():s or run in a separate process to prevent the procedure
    # from blocking here if the list of file paths does not fit into
    # the queue's buffer.
    input_queue = multiprocessing.Manager().JoinableQueue()
    for file_name in file_names:
      input_queue.put(file_name)

    pool = multiprocessing.Pool(processes=nproc)
    # Each process accumulates its own statistics in serial, and the
    # grand total is eventually collected by the main process'
    # _add_all_frames() function.
    for i in range(nproc) :
      sm = refltable_scaling_manager(self.miller_set, self.i_model, self.params)
      pool.apply_async(
        func=sm,
        args=[input_queue, db_mgr],
        callback=self._add_all_frames)
    pool.close()
    pool.join()

    # Block until the input queue has been emptied.
    input_queue.join()

  def add_frame (self, data) :
    """
    Combine the scaled data from a frame with the current overall dataset.
    Also accepts None or null_data objects, when data are unusable but we
    want to record the file as processed.
    """
    self.n_processed += 1
    if (data is None) :
      return None
    if (isinstance(data, null_data)) :
      if (data.file_error) :
        self.n_file_error += 1
      elif (data.low_signal) :
        self.n_low_signal += 1
      elif (data.wrong_bravais) :
        self.n_wrong_bravais += 1
      elif (data.wrong_cell) :
        self.n_wrong_cell += 1
      elif (data.low_resolution) :
        self.n_low_resolution += 1
      elif (data.low_correlation) :
        self.n_low_corr += 1
      elif (getattr(data,"reason",None) is not None):
        if str(data.reason)!="":
          self.failure_modes[str(data.reason)] = self.failure_modes.get(str(data.reason),0) + 1
        elif repr(type(data.reason))!="":
          self.failure_modes[repr(type(data.reason))] = self.failure_modes.get(repr(type(data.reason)),0) + 1
        else:
          self.failure_modes["other reasons"] = self.failure_modes.get("other reasons",0) + 1
      return
    if (data.accept) :
      self.n_accepted    += 1
      self.completeness  += data.completeness
      self.completeness_predictions += data.completeness_predictions
      self.summed_N      += data.summed_N
      self.summed_weight += data.summed_weight
      self.summed_wt_I   += data.summed_wt_I
      self.ISIGI.extend(data.ISIGI)
    else :
      self.n_low_corr += 1
    self.uc_values.add_cell(data.indexed_cell,
      rejected=(not data.accept))
    if not self.params.short_circuit:
      self.observations.append(data.n_obs)
    if (data.n_obs > 0) :
      frac_rejected = data.n_rejected / data.n_obs
      self.rejected_fractions.append(frac_rejected)
      self.d_min_values.append(data.d_min)
    self.corr_values.append(data.corr)
    self.wavelength.append(data.wavelength)

  def _add_all_frames (self, data) :
    """The _add_all_frames() function collects the statistics accumulated
    in @p data by the individual scaling processes in the process
    pool.  This callback function is run in serial, so it does not
    need a lock.
    """
    self.n_accepted += data.n_accepted
    self.n_file_error += data.n_file_error
    self.n_low_corr += data.n_low_corr
    self.n_low_signal += data.n_low_signal
    self.n_processed += data.n_processed
    self.n_wrong_bravais += data.n_wrong_bravais
    self.n_wrong_cell += data.n_wrong_cell
    self.n_low_resolution += data.n_low_resolution
    for key in data.failure_modes.keys():
      self.failure_modes[key] = self.failure_modes.get(key,0) + data.failure_modes[key]

    next_crystal_id = len(self.crystal_table)
    data.ISIGI['crystal_id'] += next_crystal_id
    self.ISIGI.extend(data.ISIGI)
    self.crystal_table.extend(data.crystal_table)

    self.completeness += data.completeness
    self.completeness_predictions += data.completeness_predictions
    self.summed_N += data.summed_N
    self.summed_weight += data.summed_weight
    self.summed_wt_I += data.summed_wt_I

    self.corr_values.extend(data.corr_values)
    self.d_min_values.extend(data.d_min_values)
    if not self.params.short_circuit:
      self.observations.extend(data.observations)
    self.rejected_fractions.extend(data.rejected_fractions)
    self.wavelength.extend(data.wavelength)

    self.uc_values.add_cells(data.uc_values)

  def sum_intensities(self):
    # Sum the observations of I and I/sig(I) for each reflection.
    sum_I = flex.double(self.miller_set.size(), 0.)
    sum_I_SIGI = flex.double(self.miller_set.size(), 0.)
    for i in range(len(self.ISIGI)):
      j = self.ISIGI['miller_id'][i]
      sum_I[j] += self.ISIGI['scaled_intensity'][i]
      sum_I_SIGI[j] += self.ISIGI['isigi'][i]
    return sum_I, sum_I_SIGI

if __name__ == '__main__':
  script = Script(refltable_scaling_manager)
  script.run()
