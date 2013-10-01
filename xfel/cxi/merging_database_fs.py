# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# $Id$

from __future__ import division

from cctbx.array_family import flex


MISSING_STRING = '#'

def _execute(db_commands_queue, db_results_queue, output_prefix, semaphore):
  """The _execute() function defines a consumer process that executes
  commands on the SQL database in serial.
  """

  # Acquire the semaphore when the consumer process is starting, and
  # release it on return.
  semaphore.acquire()

  rows_frame = 0 # a.k.a. frame_id
  rows_miller = 0 # a.k.a. hkl_id
  rows_observation = 0

  # Process commands from the commands queue and mark them as done.
  while True:
    command = db_commands_queue.get()
    if command is None:
      break

    table = command[0]
    order = command[1]
    parameters = command[2]
    lastrowid_key = command[3]

    # Break race condition w.r.t. initialisation!
    if rows_frame == 0:
      stream_frame = open(output_prefix + '_frame.db', 'a')
    if rows_miller == 0:
      stream_miller = open(output_prefix + '_miller.db', 'a')
    if rows_observation == 0:
      stream_observation = open(output_prefix + '_observation.db', 'a')

    if table == 'frame':
      for row in parameters:
        items = [rows_frame.__repr__()] + [MISSING_STRING] * 24
        for j in range(len(order)):
          items[order[j]] = row[j].__repr__()
        print >> stream_frame, ' '.join(items)
        rows_frame += 1
      lastrowid_value = rows_frame

    elif table == 'miller':
      for row in parameters:
        items = [rows_miller.__repr__()] + [MISSING_STRING] * 3
        for j in range(len(order)):
          items[order[j]] = row[j].__repr__()
        print >> stream_miller, ' '.join(items)
        rows_miller += 1
      lastrowid_value = rows_miller

    elif table == 'observation':
      for row in parameters:
        items = [MISSING_STRING] * 10
        for j in range(len(order)):
          items[order[j]] = row[j].__repr__()
        print >> stream_observation, ' '.join(items)
        rows_observation += 1
      lastrowid_value = rows_observation


    else:
      raise RuntimeError("Unknown table '%s'" % command[0])

    if lastrowid_key is not None:
      db_results_queue.put((lastrowid_key, lastrowid_value))
    db_commands_queue.task_done()

  # Mark the terminating None command as done.
  db_commands_queue.task_done()

  # Commit all the processed commands and join the commands queue.
  if rows_frame > 0:
    stream_frame.close()
  if rows_miller > 0:
    stream_miller.close()
  if rows_observation > 0:
    stream_observation.close()

  db_commands_queue.join()
  semaphore.release()


class manager:
  # The manager

  def __init__(self, params):
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
              self._semaphore)).start()


  def initialize_db(self, indices):
    from os import remove

    for suffix in '_frame.db', '_miller.db', '_observation.db':
      try:
        remove(self.params.output.prefix + suffix)
      except OSError, e:
        pass # deliberate - file does not exist

    self._db_commands_queue.put(('miller', (1, 2, 3), indices, None))


  def insert_frame(self, **kwargs):
    order = []
    order_dict = {'wavelength': 1,
                  'beam_x': 2,
                  'beam_y': 3,
                  'distance': 4,
                  'c_c': 5,
                  'slope': 6,
                  'offset': 7,
                  'res_ori_1': 8,
                  'res_ori_2': 9,
                  'res_ori_3': 10,
                  'res_ori_4': 11,
                  'res_ori_5': 12,
                  'res_ori_6': 13,
                  'res_ori_7': 14,
                  'res_ori_8': 15,
                  'res_ori_9': 16,
                  'rotation100_rad': 17,
                  'rotation010_rad': 18,
                  'rotation001_rad': 19,
                  'half_mosaicity_deg': 20,
                  'wave_HE_ang': 21,
                  'wave_LE_ang': 22,
                  'domain_size_ang':23,
                  'unique_file_name': 24}
    for key in kwargs.keys():
      order.append(order_dict[key])
    parameters = [kwargs.values()]

    # Pick up the index of the row just added.  The file name is
    # assumed to to serve as a unique key.
    lastrowid_key = kwargs['unique_file_name']
    self._db_commands_queue.put(('frame', order, parameters, lastrowid_key))
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
    order = []
    order_dict = {'hkl_id_0_base': 0,
                  'i': 1,
                  'sigi': 2,
                  'detector_x': 3,
                  'detector_y': 4,
                  'frame_id_0_base': 5,
                  'overload_flag': 6,
                  'original_h': 7,
                  'original_k': 8,
                  'original_l': 9}
    for key in kwargs.keys():
      order.append(order_dict[key])
    parameters = zip(*kwargs.values())
    self._db_commands_queue.put(('observation', order, parameters, None))


  def join(self):
    """The join() function closes the database.
    """

    # Terminate the consumer process by feeding it a None command and
    # wait for it to finish.
    self._db_commands_queue.put(None)
    self._db_commands_queue.join()
    self._db_results_queue.join()
    self._semaphore.acquire()


  def read_indices(self):
    millers = dict(merged_asu_hkl=flex.miller_index())
    stream = open(self.params.output.prefix + '_miller.db', 'r')
    for row in stream:
      millers['merged_asu_hkl'].append(
        tuple(int(t) for t in row.split()[1:4]))
    stream.close()
    return millers


  def read_observations(self):
    observations = {'hkl_id': flex.int(),
                    'i': flex.double(),
                    'sigi': flex.double(),
                    'frame_id': flex.int(),
                    'original_h': flex.int(),
                    'original_k': flex.int(),
                    'original_l': flex.int()}
    stream = open(self.params.output.prefix + '_observation.db', 'r')
    for row in stream:
      items = row.split()
      observations['hkl_id'].append(int(items[0]))
      observations['i'].append(float(items[1]))
      observations['sigi'].append(float(items[2]))
      observations['frame_id'].append(int(items[5]))
      observations['original_h'].append(int(items[7]))
      observations['original_k'].append(int(items[8]))
      observations['original_l'].append(int(items[9]))
    stream.close()
    return observations


  def read_frames(self):
    from cctbx.crystal_orientation import crystal_orientation
    from xfel.cxi.util import is_odd_numbered

    # XXX issues with spaces in the file name, and leading and
    # trailing single quotes (stripped below).

    frames = {'frame_id': flex.int(),
              'wavelength': flex.double(),
              'cc': flex.double(),
              'slope': flex.double(),
              'offset': flex.double(),
              'odd_numbered': flex.bool(),
              'orientation': [],
              'unit_cell': []}
    stream = open(self.params.output.prefix + '_frame.db', 'r')
    for row in stream:
      items = row.split()
      CO = crystal_orientation([float(t) for t in items[8:17]], False)
      frames['frame_id'].append(int(items[0]) - 1)
      frames['wavelength'].append(float(items[1]))
      frames['cc'].append(float(items[5]))
      frames['slope'].append(float(items[6]))
      frames['offset'].append(float(items[7]))
      frames['odd_numbered'].append(is_odd_numbered(items[24][1:-1]))
      frames['orientation'].append(CO)
      frames['unit_cell'].append(CO.unit_cell())
    stream.close()
    return frames
