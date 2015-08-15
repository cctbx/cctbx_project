from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 10/12/2014
Last Changed: 07/29/2015
Description : IOTA command-line module for running modules in order.
              Version 2.04
'''

iota_version = '2.04'
help_message = '\n{:-^70}'\
               ''.format('Integration Optimization, Triage and Analysis') + """

Auto mode
Usage: prime.iota [OPTIONS] path/to/raw/images
Generates two files, parameter file for IOTA (iota.param) and
target file for cctbx.xfel (target.phil). Integrates a random
subset of images without target cell. Outputs basic analysis.
Converts raw images into pickle format and crops to ensure that
beam center is in center of image.

Script mode
Usage: prime.iota [OPTIONS] <script>.param
Run using IOTA parameter file and target PHIL file generated from
the dry run or auto mode. Make sure that IOTA parameter file has
the path to the input image folder under "input". Converts raw
images into pickle format and modifies by cropping or padding to
ensure that beam center is in center of image. Can also blank out
beam stop shadow.

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
  return run_one_image(input_entry, init, progbar=False)

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
    # Import image
    single_image = img.SingleImage(image, init, verbose=False)
    img_object = single_image.import_image()

    # Check / convert / triage image
    if progbar:
      advance_progbar(image[0], image[1])
    img_object = single_image.convert_image()

    # Exit if the image does not have diffraction
    if img_object.triage == 'rejected':
      return img_object

  else:
    single_image = image[2]

  if 'gri' in args.mpi:
    # Grid search
    if single_image.triage == 'accepted':
      if progbar:
        advance_progbar(image[0], image[1])
      img_object = single_image.integrate('grid search')
    else:
      return single_image

  elif 'sel' in args.mpi:
    # Selection
    if len(single_image.grid) != 0 and single_image.triage == 'accepted':
      if progbar:
        advance_progbar(image[0], image[1])
      img_object = single_image.select()
    elif len(single_image.grid) == 0 or single_image.triage == 'rejected':
      return single_image

    # Exit if image not integrated
    if single_image.final['final'] == None:
      if progbar:
        advance_progbar(image[0], image[1])
      return single_image

  elif 'fin' in args.mpi:
    # Final integration
    if single_image.final['final'] != None:
      if progbar:
        advance_progbar(image[0], image[1])
      img_object = single_image.integrate('integrate')
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


  if 'ana' not in args.mpi:

    # if necessary, read in saved image objects
    if 'imp' in args.mpi:
      inp_list = init.input_list
      msg = "Importing {} images".format(len(inp_list))
    elif 'sel' in args.mpi or 'gri' in args.mpi or 'fin' in args.mpi:
      if 'gri' in args.mpi:
        msg = "Performing spotfinding grid search"
      elif 'sel' in args.mpi:
        msg = "Selecting best integrated results"
      elif 'fin' in args.mpi:
        msg = "Performing final integration"
      inp_list = [ep.load(os.path.join(init.gs_base, i)) for i in os.listdir(init.gs_base)]

    # Run all modules in order in multiprocessor mode
    #cmd.Command.start(msg)
    img_list = [[i, len(inp_list) + 1, j] for i, j in enumerate(inp_list, 1)]
    img_objects = parallel_map(iterable  = img_list,
                               func      = run_wrapper,
                               processes = init.params.n_processors,
                               method    = init.params.mp_method,
                               preserve_order = False)
    #cmd.Command.end("{} -- DONE ".format(msg))
    misc.iota_exit(iota_version, True)

  else:
    img_objects = [ep.load(os.path.join(init.gs_base, i)) for i in os.listdir(init.gs_base)]
    int_objects = [i for i in img_objects if i.final['final'] != None]
    if len(int_objects) != 0:

      # Analysis of integration results
      analysis = Analyzer(int_objects, init.logfile, iota_version, init.now)
      analysis.print_results()
      analysis.unit_cell_analysis(init.params.advanced.cluster_threshold,
                                  init.int_base)
      analysis.print_summary(init.int_base)
      analysis.make_prime_input(init.int_base)
    else:
      print "No images integrated"

  misc.iota_exit(iota_version)
