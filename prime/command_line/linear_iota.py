from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 10/12/2014
Last Changed: 09/01/2015
Description : IOTA command-line module for running modules in order.
              Version 2.10
'''

iota_version = '2.10'
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

MPI mode
Usage:
prime.linear_iota iota.param --mpi init
prime.linear_iota iota.param --mpi process
prime.linear_iota iota.param --mpi analyze

Run IOTA in three separate batches (initialization, image processing,
analysis); can use MPI (mpirun) to run the image processing step.
Can run these in sequence in a shell script or any other kind of a
submission script. Useful for huge datasets.

"""

from libtbx.easy_mp import parallel_map as easy_mp_parallel_map
from libtbx import easy_pickle as ep
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
    gs_prog = cmd.ProgressBar(title='PROCESSING')
    if prog_count < n_img:
      prog_step = 100 / n_img
      gs_prog.update(prog_count * prog_step, prog_count)
    else:
      gs_prog.finished()

  if init.params.selection.select_only.flag_on:
    single_image = ep.load(image)
    img_object = single_image.import_int_file(init)
  else:
    single_image = img.SingleImage(image, init, verbose=False)
    imported_img_object = single_image.import_image()
    img_object = imported_img_object.convert_image()

  if init.params.image_conversion.convert_only:
    return single_image

  if single_image.fail == None:
    single_image.verbose = True
    img_object = single_image.process()
  else:
    return single_image

  return img_object


# ============================================================================ #
if __name__ == "__main__":

  import os
  from prime.iota.iota_init import parse_command_args
  from libtbx import easy_pickle as ep

  args = parse_command_args(iota_version, help_message).parse_args()

  if args.mpi == None:
    misc.iota_exit(iota_version, True)

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
    if 'pro' in args.mpi:
      if init.params.selection.select_only.flag_on:
        inp_list = init.gs_img_objects
      else:
        inp_list = init.input_list
      msg = "Processing {} images".format(len(inp_list))

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
    int_objects = [i for i in img_objects if i.fail == None]
    if len(int_objects) != 0:
      analysis = Analyzer(img_objects, init.logfile, iota_version, init.now)
      analysis.print_results()
      analysis.unit_cell_analysis(init.params.analysis.cluster_threshold,
                                  init.int_base)
      analysis.print_summary(init.int_base)
      analysis.make_prime_input(init.int_base)

      # Spotfinding heatmap
      if init.params.analysis.heatmap != None:
        hm_file = "{}/{}".format(init.viz_base, 'heatmap.png')
        if init.params.analysis.heatmap == 'show':
          analysis.show_heatmap()
        elif init.params.analysis.heatmap == 'file':
          analysis.show_heatmap(show=False, hm_file=hm_file)
        elif init.params.analysis.heatmap == 'both':
          analysis.show_heatmap(hm_file=hm_file)

    else:
      print "No images integrated"

  misc.iota_exit(iota_version)
