from __future__ import division
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
      rows_frame += 1
      lastrowid_value = rows_frame

    elif table == 'observation':
      rows_observation = X["rows"].get_obj()[0]
      new_rows_observation = rows_observation+len(data[data.keys()[0]])
      X["intensity_proxy"].get_obj()[rows_observation:new_rows_observation]=data["i"]
      X["sigma_proxy"].get_obj()[rows_observation:new_rows_observation]=data["sigi"]
      X["miller_proxy"].get_obj()[rows_observation:new_rows_observation]=data["hkl_id_0_base"]
      X["frame_proxy"].get_obj()[rows_observation:new_rows_observation]=data["frame_id_0_base"]
      X["rows"].get_obj()[0] = new_rows_observation
      lastrowid_value = new_rows_observation

    else:
      raise RuntimeError("Unknown table '%s'" % command[0])
    print "FRAME",rows_frame,"OBS",X['rows'].get_obj()[0]
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
      except OSError, e:
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
    print "inserting obs:"
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
    print "writing observation pickle with %d rows"%nrows
    kwargs = dict(
      miller_lookup =      flex.size_t(data_dict["miller_proxy"].get_obj()[:nrows]),
      observed_intensity = flex.double(data_dict["intensity_proxy"].get_obj()[:nrows]),
      observed_sigI =      flex.double(data_dict["sigma_proxy"].get_obj()[:nrows]),
      frame_lookup =       flex.size_t(data_dict["frame_proxy"].get_obj()[:nrows]),
    )
    import cPickle as pickle
    pickle.dump(kwargs, open(self.params.output.prefix+"_observation.pickle","wb"),
                pickle.HIGHEST_PROTOCOL)
    pickle.dump(self.miller, open(self.params.output.prefix+"_miller.pickle","wb"),
                pickle.HIGHEST_PROTOCOL)
    return kwargs
