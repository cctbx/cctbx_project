from __future__ import absolute_import, division, print_function
import os
from dxtbx.model.experiment_list import ExperimentListFactory
from dials.array_family import flex
from six.moves import range
import json
from xfel.merging.application.input.file_lister import list_input_pairs
from xfel.merging.application.input.file_load_calculator import file_load_calculator
from xfel.merging.application.utils.memory_usage import get_memory_usage

"""
Utility functions used for reading input data
"""

def create_experiment_identifier(experiment, experiment_file_path, experiment_id):
  'Create a hashed experiment identifier based on the experiment file path, experiment index in the file, and experiment features'
  import hashlib
  exp_identifier_str = os.path.basename(experiment_file_path + \
                       str(experiment_id) + \
                       str(experiment.beam) + \
                       str(experiment.crystal) + \
                       str(experiment.detector) + \
                       ' '.join([os.path.basename(p) for p in experiment.imageset.paths()]) + \
                       ' '.join([str(i) for i in experiment.imageset.indices()]))
  hash_obj = hashlib.md5(exp_identifier_str.encode('utf-8'))
  return hash_obj.hexdigest()


def preGen_experiment_identifiers(experiments, exp_filename):
  """
  Label experiments according to image number (for multi-image files), lattice number, H
  where H is a unique per-experiment hash
  This information will be preserved in the reflection files that are optionally output
  if output.save_experiments_and_reflections=True.
  """
  done_expts = {}
  for i_exp, expt in enumerate(experiments):
    exp_hash = create_experiment_identifier(expt, exp_filename, i_exp)
    iset = expt.imageset
    path = iset.paths()[0]
    single_file_index = iset.indices()[0]
    key = "%s_%s" % (path, single_file_index)
    n_hits = done_expts.setdefault(key, 0)
    done_expts[key] = n_hits + 1
    ident = "I%s_L%s_%s" % (single_file_index, n_hits, exp_hash)
    expt.identifier = ident

#for integration pickles:
allowable_basename_endings = ["_00000.pickle",
                              ".pickle",
                              ".refl",
                              "_refined_experiments.json",
                              "_refined.expt",
                              "_experiments.json",
                              "_indexed.expt"
                             ]
def is_odd_numbered(file_name, use_hash = False):
  if use_hash:
    import hashlib
    hash_object = hashlib.md5(file_name.encode('utf-8'))
    return int(hash_object.hexdigest(), 16) % 2 == 0
  for allowable in allowable_basename_endings:
    if (file_name.endswith(allowable)):
      try:
        return int(os.path.basename(file_name).split(allowable)[-2][-1])%2==1
      except ValueError:
        file_name = os.path.basename(file_name).split(allowable)[0]
        break
  #can not find standard filename extension, instead find the last digit:
  for idx in range(1,len(file_name)+1):
    if file_name[-idx].isdigit():
      return int(file_name[-idx])%2==1
  raise ValueError
#if __name__=="__main__":
#  print is_odd_numbered("int_fake_19989.img")


def identifiers_void(*identifiers):
  """True only if all identifiers evaluate to False (eg. '', None)"""
  return not any(identifiers)


def identifiers_match(*identifiers):
  """True only if all identifiers match"""
  return len(set(identifiers)) <= 1


from xfel.merging.application.worker import worker
class simple_file_loader(worker):
  '''A class for running the script.'''

  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(simple_file_loader, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

  def __repr__(self):
    return 'Read experiments and data'

  def run(self, all_experiments, all_reflections):
    """ Load all the data using MPI """
    from dxtbx.model.experiment_list import ExperimentList
    from dials.array_family import flex

    # Both must be none or not none
    test = [all_experiments is None, all_reflections is None].count(True)
    assert test in [0,2]
    if test == 2:
      all_experiments = ExperimentList()
      all_reflections = flex.reflection_table()
      starting_expts_count = starting_refls_count = 0
    else:
      starting_expts_count = len(all_experiments)
      starting_refls_count = len(all_reflections)
    self.logger.log("Initial number of experiments: %d; Initial number of reflections: %d"%(starting_expts_count, starting_refls_count))

    if not self.check_psana2(): return all_experiments, all_reflections

    # Generate and send a list of file paths to each worker
    if self.mpi_helper.rank == 0:
      file_list = list_input_pairs(self.params)
      self.logger.log("Built an input list of %d json/pickle file pairs"%(len(file_list)))

      # optionally write a file list mapping to disk, useful in post processing if save_experiments_and_reflections=True
      file_id_from_names = None
      if self.params.output.expanded_bookkeeping:
        apath = lambda x: os.path.abspath(x)
        file_names_from_id = {i_f: tuple(map(apath, exp_ref_pair)) for i_f, exp_ref_pair in enumerate(file_list)}
        with open(os.path.join(self.params.output.output_dir, "file_list_map.json"), "w") as o:
          json.dump(file_names_from_id, o)
        file_id_from_names = {tuple(map(apath, exp_ref_pair)): i_f for i_f, exp_ref_pair in enumerate(file_list)}

      per_rank_file_list = file_load_calculator(self.params, file_list, self.logger).\
                              calculate_file_load(available_rank_count = self.mpi_helper.size)
      self.logger.log('Transmitting a list of %d lists of json/pickle file pairs'%(len(per_rank_file_list)))
      transmitted = per_rank_file_list, file_id_from_names
    else:
      transmitted = None

    self.logger.log_step_time("BROADCAST_FILE_LIST")
    new_file_list, file_names_mapping = self.mpi_helper.comm.bcast(transmitted, root = 0)
    new_file_list = new_file_list[self.mpi_helper.rank] if self.mpi_helper.rank < len(new_file_list) else None
    self.logger.log_step_time("BROADCAST_FILE_LIST", True)

    # Load the data
    self.logger.log_step_time("LOAD")
    if new_file_list is not None:
      self.logger.log("Received a list of %d json/pickle file pairs"%len(new_file_list))
      for experiments_filename, reflections_filename in new_file_list:
        self.logger.log("Reading %s %s"%(experiments_filename, reflections_filename))
        experiments = ExperimentListFactory.from_json_file(experiments_filename, check_format = self.params.input.read_image_headers)
        reflections = flex.reflection_table.from_file(reflections_filename)
        if self.params.output.expanded_bookkeeping:
          # NOTE: these are un-prunable
          reflections["input_refl_index"] = flex.int(
            list(range(len(reflections))))
          reflections["original_id"] = reflections['id']
          assert file_names_mapping is not None
          exp_ref_pair = os.path.abspath(experiments_filename), os.path.abspath(reflections_filename)
          this_refl_fileMappings = [file_names_mapping[exp_ref_pair]]*len(reflections)
          reflections["file_list_mapping"] = flex.int(this_refl_fileMappings)
        self.logger.log("Data read, prepping")

        if 'intensity.sum.value' in reflections:
          reflections['intensity.sum.value.unmodified'] = reflections['intensity.sum.value'] * 1
        if 'intensity.sum.variance' in reflections:
          reflections['intensity.sum.variance.unmodified'] = reflections['intensity.sum.variance'] * 1

        new_ids = flex.int(len(reflections), -1)
        eid = reflections.experiment_identifiers()
        eid_copy = dict(eid)
        for k in eid.keys():
          del eid[k]

        if self.params.output.expanded_bookkeeping:
          preGen_experiment_identifiers(experiments, experiments_filename)
        for experiment_id, experiment in enumerate(experiments):
          # select reflections of the current experiment
          refls_sel = reflections['id'] == experiment_id

          if refls_sel.count(True) == 0: continue

          refls_identifier = eid_copy.get(experiment_id, '')
          if identifiers_void(experiment.identifier, refls_identifier) \
                  or self.params.input.override_identifiers:
            new_identifier = create_experiment_identifier(
              experiment, experiments_filename, experiment_id)
            experiment.identifier = new_identifier
          elif not identifiers_match(experiment.identifier, refls_identifier):
            m = 'Expt and refl identifier mismatch: "{}" in {} vs "{}" in {}'
            raise KeyError(m.format(experiment.identifier, experiments_filename,
                                    refls_identifier, reflections_filename))

          if not self.params.input.keep_imagesets:
            experiment.imageset = None
          all_experiments.append(experiment)

          # Reflection 'id' is unique within this rank; experiment.identifier is unique globally
          new_id = len(all_experiments)-1
          eid[new_id] = experiment.identifier
          new_ids.set_selected(refls_sel, new_id)
        assert (new_ids < 0).count(True) == 0, "Not all reflections accounted for"
        reflections['id'] = new_ids
        all_reflections.extend(reflections)
    else:
      self.logger.log("Received a list of 0 json/pickle file pairs")
    self.logger.log_step_time("LOAD", True)

    self.logger.log('Read %d experiments consisting of %d reflections'%(len(all_experiments)-starting_expts_count, len(all_reflections)-starting_refls_count))
    self.logger.log("Memory usage: %d MB"%get_memory_usage())

    all_reflections = self.prune_reflection_table_keys(all_reflections)

    # Do we have any data?
    from xfel.merging.application.utils.data_counter import data_counter
    data_counter(self.params, mpi_helper = self.mpi_helper, mpi_logger = self.logger).count(all_experiments, all_reflections)
    return all_experiments, all_reflections

  def prune_reflection_table_keys(self, reflections):
    from xfel.merging.application.reflection_table_utils import reflection_table_utils
    reflections = reflection_table_utils.prune_reflection_table_keys(reflections=reflections,
                    keys_to_keep=['id', 'intensity.sum.value', 'intensity.sum.variance', 'miller_index', 'miller_index_asymmetric', \
                                  's1', 'intensity.sum.value.unmodified', 'intensity.sum.variance.unmodified',
                                  'kapton_absorption_correction', 'flags'],
                    keys_to_ignore=self.params.input.persistent_refl_cols)
    self.logger.log("Pruned reflection table")
    self.logger.log("Memory usage: %d MB"%get_memory_usage())
    return reflections

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(simple_file_loader)
