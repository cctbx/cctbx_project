from __future__ import absolute_import, division, print_function
from six.moves import range
from cctbx.array_family import flex
from cctbx.miller import match_multi_indices
from cctbx.miller import set as mset
from scitbx import matrix

from xfel.merging.database.merging_database import manager
class progress_manager(manager):
  def __init__(self,params,db_experiment_tag,trial,rungroup_id,run):
    self.params = params
    self.db_experiment_tag = db_experiment_tag
    self.trial = trial
    self.rungroup_id = rungroup_id
    self.run = run

  def get_trial_id(self, cursor):
    query = "SELECT trial_id FROM %s_trials WHERE %s_trials.trial = %d"%(self.db_experiment_tag, self.db_experiment_tag, self.trial)
    cursor.execute(query)
    assert cursor.rowcount == 1
    return int(cursor.fetchall()[0][0])

  def get_run_id(self, cursor):
    query = "SELECT run_id FROM %s_runs WHERE %s_runs.run = %d"%(self.db_experiment_tag, self.db_experiment_tag, self.run)
    cursor.execute(query)
    assert cursor.rowcount == 1
    return int(cursor.fetchall()[0][0])

  def get_HKL(self,cursor):
    name = self.db_experiment_tag
    query = '''SELECT H,K,L,%s_hkls.hkl_id from %s_hkls,%s_isoforms WHERE %s_hkls.isoforms_isoform_id = %s_isoforms.isoform_id AND %s_isoforms.name = "%s"'''%(
            name, name, name, name, name, name, self.params["identified_isoform"])
    cursor.execute(query)
    ALL = cursor.fetchall()
    indices = flex.miller_index([(a[0],a[1],a[2]) for a in ALL])
    miller_id = flex.int([a[3] for a in ALL])
    self.miller_set = mset(crystal_symmetry=self.params["observations"][0].crystal_symmetry(), indices=indices)
    self.miller_set_id = miller_id
    # might have to change this to isoform_id next iteration
    cursor.execute('SELECT isoform_id FROM %s_isoforms WHERE name = "%s"'%(
            name, self.params["identified_isoform"]))

    self.isoform_id = cursor.fetchall()[0][0]
    return indices,miller_id

  def connection(self):
    pass

  def scale_frame_detail(self,timestamp,cursor,do_inserts=True,result=None):#, file_name, db_mgr, out):
    if result is None:
      result = self.params

    # If the pickled integration file does not contain a wavelength,
    # fall back on the value given on the command line.  XXX The
    # wavelength parameter should probably be removed from master_phil
    # once all pickled integration files contain it.
    wavelength = result["wavelength"]
    assert (wavelength > 0)

    # Do not apply polarization correction here, as this requires knowledge of
    # pixel size at minimum, and full detector geometry in general.  The optimal
    # redesign would be to apply the polarization correction just after the integration
    # step in the integration code.
    print("Step 3. Correct for polarization.")
    observations = result["observations"][0]
    indexed_cell = observations.unit_cell()

    observations_original_index = observations.deep_copy()

    assert len(observations_original_index.indices()) == len(observations.indices())

    # Now manipulate the data to conform to unit cell, asu, and space group
    # of reference.  The resolution will be cut later.
    # Only works if there is NOT an indexing ambiguity!
    #observations = observations.customized_copy(
    #  anomalous_flag=not self.params.merge_anomalous,
    #  crystal_symmetry=self.miller_set.crystal_symmetry()
    #  ).map_to_asu()

    #observations_original_index = observations_original_index.customized_copy(
    #  anomalous_flag=not self.params.merge_anomalous,
    #  crystal_symmetry=self.miller_set.crystal_symmetry()
    #  )
    observations = observations.customized_copy(anomalous_flag=False).map_to_asu()
    print("Step 4. Filter on global resolution and map to asu")

    #observations.show_summary(f=out, prefix="  ")
    from rstbx.dials_core.integration_core import show_observations
    show_observations(observations)


    print("Step 6.  Match to reference intensities, filter by correlation, filter out negative intensities.")
    assert len(observations_original_index.indices()) \
      ==   len(observations.indices())

    # Ensure that match_multi_indices() will return identical results
    # when a frame's observations are matched against the
    # pre-generated Miller set, self.miller_set, and the reference
    # data set, self.i_model.  The implication is that the same match
    # can be used to map Miller indices to array indices for intensity
    # accumulation, and for determination of the correlation
    # coefficient in the presence of a scaling reference.
    self.miller_set.show_summary(prefix="mset ")

    matches = match_multi_indices(
      miller_indices_unique=self.miller_set.indices(),
      miller_indices=observations.indices())

    slope = 1.0
    offset = 0.0

    print(result.get("sa_parameters")[0])
    have_sa_params = ( type(result.get("sa_parameters")[0]) == type(dict()) )

    observations_original_index_indices = observations_original_index.indices()
    print(result.keys())
    kwargs = {'wavelength': wavelength,
              'beam_x': result['xbeam'],
              'beam_y': result['ybeam'],
              'distance': result['distance'],
              'slope': slope,
              'offset': offset,
              'unique_file_name': timestamp,
              'eventstamp':timestamp,
              'sifoil': 0.0}

    trial_id = self.get_trial_id(cursor)
    run_id = self.get_run_id(cursor)
    kwargs["trials_id"] = trial_id
    kwargs["rungroups_id"] = self.rungroup_id
    kwargs["runs_run_id"] = run_id
    kwargs["isoforms_isoform_id"] = self.isoform_id
    res_ori_direct = matrix.sqr(
        observations.unit_cell().orthogonalization_matrix()).transpose().elems

    kwargs['res_ori_1'] = res_ori_direct[0]
    kwargs['res_ori_2'] = res_ori_direct[1]
    kwargs['res_ori_3'] = res_ori_direct[2]
    kwargs['res_ori_4'] = res_ori_direct[3]
    kwargs['res_ori_5'] = res_ori_direct[4]
    kwargs['res_ori_6'] = res_ori_direct[5]
    kwargs['res_ori_7'] = res_ori_direct[6]
    kwargs['res_ori_8'] = res_ori_direct[7]
    kwargs['res_ori_9'] = res_ori_direct[8]

    kwargs['mosaic_block_rotation'] = result.get("ML_half_mosaicity_deg",[float("NaN")])[0]
    kwargs['mosaic_block_size'] = result.get("ML_domain_size_ang",[float("NaN")])[0]
    kwargs['ewald_proximal_volume'] = result.get("ewald_proximal_volume",[float("NaN")])[0]


    sql, parameters = self._insert(
      table='`%s_frames`' % self.db_experiment_tag,
      **kwargs)
    print(sql)
    print(parameters)
    results = {'frame':[sql, parameters, kwargs]}
    if do_inserts:
      cursor.execute(sql, parameters[0])
      frame_id = cursor.lastrowid
    else:
      frame_id = None

    xypred = result["mapped_predictions"][0]
    indices = flex.size_t([pair[1] for pair in matches.pairs()])

    sel_observations = flex.intersection(
      size=observations.data().size(),
      iselections=[indices])
    set_original_hkl = observations_original_index_indices.select(
      flex.intersection(
        size=observations_original_index_indices.size(),
        iselections=[indices]))
    set_xypred = xypred.select(
      flex.intersection(
        size=xypred.size(),
        iselections=[indices]))
    ''' debugging printout
    print len(observations.data())
    print len(indices)
    print len(sel_observations)
    for x in range(len(observations.data())):
      print x,observations.indices().select(sel_observations)[x],
      print set_original_hkl[x],
      index_into_hkl_id = matches.pairs()[x][0]
      print index_into_hkl_id,
      print self.miller_set.indices()[index_into_hkl_id],
      cursor.execute('SELECT H,K,L FROM %s_hkls WHERE hkl_id = %d'%(
            self.db_experiment_tag, self.miller_set_id[index_into_hkl_id]))

      print cursor.fetchall()[0]
    '''
    print("Adding %d observations for this frame"%(len(sel_observations)))
    kwargs = {'hkls_id': self.miller_set_id.select(flex.size_t([pair[0] for pair in matches.pairs()])),
              'i': observations.data().select(sel_observations),
              'sigi': observations.sigmas().select(sel_observations),
              'detector_x_px': [xy[0] for xy in set_xypred],
              'detector_y_px': [xy[1] for xy in set_xypred],
              'frames_id': [frame_id] * len(matches.pairs()),
              'overload_flag': [0] * len(matches.pairs()),
              'original_h': [hkl[0] for hkl in set_original_hkl],
              'original_k': [hkl[1] for hkl in set_original_hkl],
              'original_l': [hkl[2] for hkl in set_original_hkl],
              'frames_rungroups_id': [self.rungroup_id] * len(matches.pairs()),
              'frames_trials_id': [trial_id] * len(matches.pairs()),
              'panel': [0] * len(matches.pairs())
    }
    if do_inserts:
      # For MySQLdb executemany() is six times slower than a single big
      # execute() unless the "values" keyword is given in lowercase
      # (http://sourceforge.net/p/mysql-python/bugs/305).
      #
      # See also merging_database_sqlite3._insert()
      query = ("INSERT INTO `%s_observations` (" % self.db_experiment_tag) \
              + ", ".join(kwargs.keys()) + ") values (" \
              + ", ".join(["%s"] * len(kwargs.keys())) + ")"
      try:
        parameters = zip(*kwargs.values())
      except TypeError:
        parameters = [kwargs.values()]
      cursor.executemany(query, parameters)
      #print "done execute many"
      #print cursor._last_executed
      results['observations'] = [query, parameters, kwargs]
    else:
      # since frame_id isn't valid in the query here, don't include a sql statement or parameters array in the results
      results['observations'] = [None, None, kwargs]

    return results

'''START HERE
scons ; cxi.xtc_process input.cfg=xppi6115/LI61-PSII-db.cfg input.experiment=xppi6115 input.run_num=137
'''
