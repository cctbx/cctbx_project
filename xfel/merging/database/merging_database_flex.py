from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex

def _execute(db_commands_queue, db_results_queue, output_prefix, semaphore, X):
  """The _execute() function defines a consumer process that executes
  commands on the SQL database in serial.
  """
  # Acquire the semaphore when the consumer process is starting, and
  # release it on return.
  semaphore.acquire()

  rows_frame = 0 # a.k.a. frame_id

  # Process commands from the commands queue and mark them as done.
  while True:
    command = db_commands_queue.get()
    if command is None:
      break
    table = command[0]
    data = command[1]
    lastrowid_key = command[2]

    if table == 'frame':
      items = [0]*len(order_dict)
      for key,val in data.items():
        items[order_dict[key]]=val
      characters = ' '.join([str(i) for i in items])+'\n'
      X["xtal_proxy"].get_obj()[rows_frame*600:rows_frame*600+len(characters)]=characters
      rows_frame += 1
      lastrowid_value = rows_frame

    elif table == 'observation':
      rows_observation = X["rows"].get_obj()[0]
      new_rows_observation = rows_observation+len(data[list(data.keys())[0]]) # XXX FIXME
      X["intensity_proxy"].get_obj()[rows_observation:new_rows_observation]=data["i"]
      X["sigma_proxy"].get_obj()[rows_observation:new_rows_observation]=data["sigi"]
      X["miller_proxy"].get_obj()[rows_observation:new_rows_observation]=data["hkl_id_0_base"]
      X["frame_proxy"].get_obj()[rows_observation:new_rows_observation]=data["frame_id_0_base"]
      X["H_proxy"].get_obj()[rows_observation:new_rows_observation]=data["original_h"]
      X["K_proxy"].get_obj()[rows_observation:new_rows_observation]=data["original_k"]
      X["L_proxy"].get_obj()[rows_observation:new_rows_observation]=data["original_l"]
      X["rows"].get_obj()[0] = new_rows_observation
      lastrowid_value = new_rows_observation

    else:
      raise RuntimeError("Unknown table '%s'" % command[0])
    print("FRAME",rows_frame,"OBS",X['rows'].get_obj()[0])
    if lastrowid_key is not None:
      db_results_queue.put((lastrowid_key, lastrowid_value))
    db_commands_queue.task_done()

  # Mark the terminating None command as done.
  db_commands_queue.task_done()

  db_commands_queue.join()
  semaphore.release()

class manager:
  # The manager

  def __init__(self, params, data_proxy):
    import multiprocessing

    self.params = params
    mgr = multiprocessing.Manager()
    self._db_commands_queue = mgr.JoinableQueue()
    self._db_results_queue = mgr.JoinableQueue()
    self._semaphore = mgr.Semaphore()

    multiprocessing.Process(
        target=_execute,
        args=(self._db_commands_queue,
              self._db_results_queue,
              self.params.output.prefix,
              self._semaphore,
              data_proxy)).start()

  def initialize_db(self, indices):
    from os import remove

    for suffix in '_frame.pickle', '_miller.pickle', '_observation.pickle':
      try:
        remove(self.params.output.prefix + suffix)
      except OSError as e:
        pass # deliberate - file does not exist

    self.miller = indices

  def insert_frame(self, **kwargs):
    # Pick up the index of the row just added.  The string is
    # assumed to to serve as a unique key.
    lastrowid_key = kwargs['unique_file_name']
    self._db_commands_queue.put(('frame', kwargs, lastrowid_key))
    while True:
      item = self._db_results_queue.get()
      self._db_results_queue.task_done()
      if item[0] == kwargs['unique_file_name']:
        # Entry in the observation table is zero-based.
        return item[1] - 1
      else:
        # If the key does not match, put it back in the queue for
        # someone else to pick up.
        self._db_results_queue.put(item)

  def insert_observation(self, **kwargs):
    print("inserting obs:")
    self._db_commands_queue.put(('observation', kwargs, None))

  def join(self,data_dict):
    """The join() function closes the database.
    """
    # Terminate the consumer process by feeding it a None command and
    # wait for it to finish.
    self._db_commands_queue.put(None)
    self._db_commands_queue.join()
    self._db_results_queue.join()
    self._semaphore.acquire()
    nrows = data_dict["rows"].get_obj()[0]
    print("writing observation pickle with %d rows"%nrows)
    kwargs = dict(
      miller_lookup =      flex.size_t(data_dict["miller_proxy"].get_obj()[:nrows]),
      observed_intensity = flex.double(data_dict["intensity_proxy"].get_obj()[:nrows]),
      observed_sigI =      flex.double(data_dict["sigma_proxy"].get_obj()[:nrows]),
      frame_lookup =       flex.size_t(data_dict["frame_proxy"].get_obj()[:nrows]),
      original_H =         flex.int   (data_dict["H_proxy"].get_obj()[:nrows]),
      original_K =         flex.int   (data_dict["K_proxy"].get_obj()[:nrows]),
      original_L =         flex.int   (data_dict["L_proxy"].get_obj()[:nrows]),
    )
    from six.moves import cPickle as pickle
    pickle.dump(kwargs, open(self.params.output.prefix+"_observation.pickle","wb"),
                pickle.HIGHEST_PROTOCOL)
    pickle.dump(self.miller, open(self.params.output.prefix+"_miller.pickle","wb"),
                pickle.HIGHEST_PROTOCOL)
    pickle.dump(data_dict["xtal_proxy"].get_obj().raw.replace('\0','').strip(),
                        open(self.params.output.prefix+"_frame.pickle","wb"),pickle.HIGHEST_PROTOCOL)
    return kwargs

order_dict = {'wavelength': 0,
                  'beam_x': 1,
                  'beam_y': 2,
                  'distance': 3,
                  'res_ori_1': 4,
                  'res_ori_2': 5,
                  'res_ori_3': 6,
                  'res_ori_4': 7,
                  'res_ori_5': 8,
                  'res_ori_6': 9,
                  'res_ori_7': 10,
                  'res_ori_8': 11,
                  'res_ori_9': 12,
                  'half_mosaicity_deg': 13,
                  'domain_size_ang':14,
                  'unique_file_name': 15}
class read_experiments(object):
  def __init__(self,params):
    from six.moves import cPickle as pickle
    from dxtbx.model import BeamFactory
    from dxtbx.model import DetectorFactory
    from dxtbx.model.crystal import CrystalFactory
    from cctbx.crystal_orientation import crystal_orientation,basis_type
    from dxtbx.model import Experiment, ExperimentList
    from scitbx import matrix
    self.experiments = ExperimentList()
    self.unique_file_names = []

    self.params = params
    data = pickle.load(open(self.params.output.prefix+"_frame.pickle","rb"))
    frames_text = data.split("\n")

    for item in frames_text:
      tokens = item.split(' ')
      wavelength = float(tokens[order_dict["wavelength"]])

      beam = BeamFactory.simple(wavelength = wavelength)

      detector = DetectorFactory.simple(
        sensor = DetectorFactory.sensor("PAD"), # XXX shouldn't hard code for XFEL
        distance = float(tokens[order_dict["distance"]]),
        beam_centre = [float(tokens[order_dict["beam_x"]]), float(tokens[order_dict["beam_y"]])],
        fast_direction = "+x",
        slow_direction = "+y",
        pixel_size = [self.params.pixel_size,self.params.pixel_size],
        image_size = [1795,1795],  # XXX obviously need to figure this out
        )

      reciprocal_matrix = matrix.sqr([float(tokens[order_dict[k]]) for k in [
'res_ori_1','res_ori_2','res_ori_3','res_ori_4','res_ori_5','res_ori_6','res_ori_7','res_ori_8','res_ori_9']])
      ORI = crystal_orientation(reciprocal_matrix, basis_type.reciprocal)
      direct = matrix.sqr(ORI.direct_matrix())
      transfer_dict = dict(__id__="crystal",
                           ML_half_mosaicity_deg=float(tokens[order_dict["half_mosaicity_deg"]]),
                           ML_domain_size_ang=float(tokens[order_dict["domain_size_ang"]]),
                           real_space_a = matrix.row(direct[0:3]),
                           real_space_b = matrix.row(direct[3:6]),
                           real_space_c = matrix.row(direct[6:9]),
                           space_group_hall_symbol = self.params.target_space_group.type().hall_symbol(),
                           )
      crystal = CrystalFactory.from_dict(transfer_dict)
      """ old code reflects python-based crystal model
      crystal = Crystal(
        real_space_a = matrix.row(direct[0:3]),
        real_space_b = matrix.row(direct[3:6]),
        real_space_c = matrix.row(direct[6:9]),
        space_group_symbol = self.params.target_space_group.type().lookup_symbol(),
        mosaicity = float(tokens[order_dict["half_mosaicity_deg"]]),
      )
      crystal.domain_size = float(tokens[order_dict["domain_size_ang"]])
      """
      #if isoform is not None:
      #  newB = matrix.sqr(isoform.fractionalization_matrix()).transpose()
      #  crystal.set_B(newB)

      self.experiments.append(Experiment(beam=beam,
                                  detector=None, #dummy for now
                                  crystal=crystal))
      self.unique_file_names.append(tokens[order_dict["unique_file_name"]])

    self.show_summary()

  def get_experiments(self):
    return self.experiments
  def get_files(self):
    return self.unique_file_names

  def show_summary(self):
    w = flex.double([e.beam.get_wavelength() for e in self.experiments])
    stats=flex.mean_and_variance(w)
    print("Wavelength mean and standard deviation:",stats.mean(),stats.unweighted_sample_standard_deviation())
    uc = [e.crystal.get_unit_cell().parameters() for e in self.experiments]
    a = flex.double([u[0] for u in uc])
    stats=flex.mean_and_variance(a)
    print("Unit cell a mean and standard deviation:",stats.mean(),stats.unweighted_sample_standard_deviation())
    b = flex.double([u[1] for u in uc])
    stats=flex.mean_and_variance(b)
    print("Unit cell b mean and standard deviation:",stats.mean(),stats.unweighted_sample_standard_deviation())
    c = flex.double([u[2] for u in uc])
    stats=flex.mean_and_variance(c)
    print("Unit cell c mean and standard deviation:",stats.mean(),stats.unweighted_sample_standard_deviation())
    d = flex.double([e.crystal.get_domain_size_ang() for e in self.experiments])
    stats=flex.mean_and_variance(d)
    # NOTE XXX FIXME:  cxi.index seems to record the half-domain size; report here the full domain size
    print("Domain size mean and standard deviation:",2.*stats.mean(),2.*stats.unweighted_sample_standard_deviation())
