from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME iota.mpi

'''
Author      : Lyubimov, A.Y.
Created     : 10/12/2014
Last Changed: 01/06/2016
Description : IOTA command-line module for running modules in order.
              Version 2.24
'''

iota_version = '2.24'
help_message = '\n{:-^70}'\
               ''.format('Integration Optimization, Triage and Analysis') + """

MPI mode
Usage:
iota.mpi iota.param --mpi init
mpirun -n [num_cores] iota.mpi iota.param --mpi import
mpirun -n [num_cores] iota.mpi iota.param --mpi gridsearch
mpirun -n [num_cores] iota.mpi iota.param --mpi select
mpirun -n [num_cores] iota.mpi iota.param --mpi final
iota.mpi iota.param --mpi analyze

Run IOTA in three separate batches (initialization, image processing,
analysis); can use MPI (mpirun) to run the image processing step.
Can run these in sequence in a shell script or any other kind of a
submission script. Useful for huge datasets. (NOTE: make sure that,
in iota.param, mp_method = mpi)

"""

from libtbx.easy_mp import parallel_map as easy_mp_parallel_map
from prime.iota.iota_init import InitAll
from prime.iota.iota_analysis import Analyzer
import prime.iota.iota_image as img
import prime.iota.iota_cmd as cmd
import prime.iota.iota_misc as misc


def parallel_map (
    func,
    iterable,
    params=None,
    processes=1,
    method="multiprocessing",
    qsub_command=None,
    asynchronous=True,
    callback=None,
    preserve_order=True,
    preserve_exception_message=False,
    use_manager=True) :

  if method == 'mpi' and processes <= 2:
    method = 'multiprocessing'

  if method == 'mpi':
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank, size = misc.get_mpi_rank_and_size()

    assert size == processes
    assert not preserve_order
    results = []

    # use a client/server approach to be sure every process is busy as much as possible
    # only do this if there are more than 2 processes, as one process will be a server
    if rank == 0:
      # server process
      for task_id in xrange(len(iterable)):
        # a client process will indicate it's ready by sending its rank
        request = comm.recv(source=MPI.ANY_SOURCE)
        comm.send(task_id,dest=request)

      # send a stop command to each process
      for rankreq in xrange(1,size):
        #rankreq = comm.recv(source=MPI.ANY_SOURCE)
        comm.send('endjob',dest=rankreq)

      for proc_id in xrange(1,size):
        res = comm.recv(source=proc_id)
        try:
          results.extend(res)
        except TypeError:
          results.extend(comm.recv(source=proc_id))

      return results
    else:
      # client process
      results = []
      while True:
        # inform the server this process is ready for an event
        comm.send(rank,dest=0)
        task_id = comm.recv(source=0)
        if task_id == 'endjob': break
        results.append(func(iterable[task_id]))

      comm.send(results,dest=0)

  else:
    return easy_mp_parallel_map(func, iterable, params, processes, method,
                                qsub_command, asynchronous, callback,
                                preserve_order, preserve_exception_message,
                                use_manager)

def run_wrapper(input_entry):
  """ Multiprocessor wrapper for image conversion  """
  return run_one_image(input_entry, init, progbar=True)

def run_one_image(image, init, progbar=True):

  def advance_progbar(prog_count, n_img):
    if 'imp' in args.mpi:
      ptitle = "IMPORTING IMAGES"
    elif 'gri' in args.mpi:
      ptitle = "GRID SEARCH"
    elif 'sel' in args.mpi:
      ptitle = "SELECTING"
    elif 'fin' in args.mpi:
      ptitle = "INTEGRATING"

    gs_prog = cmd.ProgressBar(title=ptitle)
    if prog_count < n_img:
      prog_step = 100 / n_img
      gs_prog.update(prog_count * prog_step, prog_count)
    else:
      gs_prog.finished()

  if 'imp' in args.mpi:
    single_image = img.SingleImage(image, init, verbose=True)

    # Import / check / convert / triage image
    if progbar:
      advance_progbar(image[0], image[1])
    img_object = single_image.convert_image()

    # Exit if the image does not have diffraction
    if img_object.fail == 'failed triage':
      return img_object

  if 'gri' in args.mpi:
    # Grid search
    single_image = image[2][0]
    grid_point = image[2][1]
    if single_image.fail == None:
      if progbar:
        advance_progbar(image[0], image[1])
      img_object = single_image.integrate_cctbx(tag='split grid', grid_point=grid_point)
    else:
      return single_image

  elif 'sel' in args.mpi:
    # Selection
    single_image = image[2]
    if single_image.fail == None:
      if progbar:
        advance_progbar(image[0], image[1])
      img_object = single_image.select_cctbx()
    else:
      return single_image

    # Exit if image not integrated
    if single_image.final['final'] == None:
      if progbar:
        advance_progbar(image[0], image[1])
      return single_image

  elif 'fin' in args.mpi:
    # Final integration
    single_image = image[2]
    if single_image.final['final'] != None:
      if progbar:
        advance_progbar(image[0], image[1])
      img_object = single_image.integrate_cctbx('integrate')
    else:
      return single_image

  return img_object


# ============================================================================ #
if __name__ == "__main__":

  import os
  from prime.iota.iota_init import parse_command_args
  from libtbx import easy_pickle as ep

  args = parse_command_args(iota_version, help_message).parse_args()
  if "ini" in args.mpi:
    # Initialize IOTA parameters and log
    init = InitAll(iota_version, help_message)
    init.run()
    ep.dump('iota.ini', init)
    misc.iota_exit(iota_version, True)
  else:
    if os.path.isfile('iota.ini'):
      init = ep.load('iota.ini')
    else:
      print 'This run not initialized. Use "--mpi init" option'
      misc.iota_exit(iota_version, True)

  if 'ana' in args.mpi:
    img_objects = [ep.load(os.path.join(init.obj_base, i)) for i in os.listdir(init.obj_base)]
    int_objects = [i for i in img_objects if i.fail == None]
    if len(int_objects) != 0:

      # Analysis of integration results
      analysis = Analyzer(int_objects, init.logfile, iota_version, init.now)
      analysis.print_results()
      analysis.unit_cell_analysis(init.params.analysis.cluster_threshold,
                                  init.int_base)
      analysis.print_summary(init.int_base)
      analysis.make_prime_input(init.int_base)
    else:
      print "No images integrated"

  else:

    if 'imp' in args.mpi:
      inp_list = init.input_list
      msg = "Importing {} images".format(len(inp_list))

    elif 'gri' in args.mpi:
      obj_list = [ep.load(os.path.join(init.obj_base, i)) for i in os.listdir(init.obj_base)]
      inp_list = [[obj, item] for obj in obj_list for item in range(len(obj.grid))]

    elif 'sel' in args.mpi or 'fin' in args.mpi:
      inp_list = [ep.load(os.path.join(init.obj_base, i)) for i in os.listdir(init.obj_base)]

    # Run all modules in order in multiprocessor mode
    img_list = [[i, len(inp_list) + 1, j] for i, j in enumerate(inp_list, 1)]
    if init.params.n_processors > len(img_list):
      init.params.n_processors = len(img_list)
    img_objects = parallel_map(iterable  = img_list,
                               func      = run_wrapper,
                               processes = init.params.n_processors,
                               method    = init.params.mp_method,
                               preserve_order = False)
    misc.iota_exit(iota_version, True)

misc.iota_exit(iota_version)
