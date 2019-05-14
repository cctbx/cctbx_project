from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME cxi.index
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT

from xfel.cxi.display_spots import run_one_index
from libtbx.utils import Usage, Sorry
import libtbx.option_parser
import sys,os

if (__name__ == "__main__"):
  command_line = (libtbx.option_parser.option_parser(
    usage="%s [-d] [-s] [-n num_procs] [-o output_dir] [-b output_basename] [-e extension] target=targetpath files" % libtbx.env.dispatcher_name,
    more_help=["Target: the phil file containing further indexing/integration parameters"])
                  .option(None, "--no_display", "-d",
                          action="store_true",
                          default=False,
                          dest="no_display",
                          help="Do not show indexing graphics")
                  .option(None, "--skip_processed", "-s",
                          action="store_true",
                          default=False,
                          dest="skip_processed",
                          help="Skip files that have alread been processed")
                  .option(None, "--num_procs", "-n",
                          type="int",
                          default=1,
                          dest="num_procs",
                          help="Number of processors to use")
                  .option(None, "--output_dir", "-o"    ,
                          type="string",
                          default=None,
                          dest="output_dir",
                          help="Directory for integration pickles")
                  .option(None, "--output_basename", "-b",
                          type="string",
                          default="int_",
                          dest="output_basename",
                          help="String to append to the front of output integration pickles")
                  .option(None, "--extension", "-e",
                          type="string",
                          default=".pickle",
                          dest="extension",
                          help="File extension use to filter input files if a directory is given as input")
                  ).process(args=sys.argv[1:])

  files = [arg for arg in command_line.args if os.path.isfile(arg)]
  dirs = [arg for arg in command_line.args if os.path.isdir(arg)]
  for directory in dirs:
    for path in os.listdir(directory):
      if os.path.splitext(path)[1] == command_line.options.extension:
        files.append(os.path.join(directory, path))

  arguments = [arg for arg in command_line.args if not os.path.isfile(arg) and not os.path.isdir(arg)]

  found_it = False
  for arg in arguments:
    if "target=" in arg:
      found_it = True
      break
  if not found_it:
    raise Usage(command_line.parser.usage)

  if command_line.options.no_display:
    display = False
    arguments.append('--nodisplay')
  else:
    display = True

  assert command_line.options.num_procs > 0
  if command_line.options.output_dir is not None and \
    not os.path.isdir(command_line.options.output_dir):
    raise Sorry("Output dir %s doesn't exist"%command_line.options.output_dir)

  def do_work(item):
    file, arguments, kwargs = item
    try:
      run_one_index(file, *arguments, **({'display':display}))
    except Exception as e:
      if hasattr(e, "classname"):
        print(e.classname, "for %s:"%file, end=' ')
      else:
        print("Indexing error for %s:"%file, end=' ')
      print(e)

  if command_line.options.num_procs == 1:
    for file in files:
      if command_line.options.output_dir is not None:
        int_pickle_path = os.path.join(command_line.options.output_dir, \
          command_line.options.output_basename + os.path.basename(file))
        if command_line.options.skip_processed and os.path.isfile(int_pickle_path):
          print(file, "already processed, skipping")
          continue
        arguments.append("indexing.completeness_pickle=%s"%int_pickle_path)
      do_work((file, arguments, ({'display':display})))
  else:
    import multiprocessing, copy

    def worker():
      for item in iter( q.get, None ):
        do_work(item)
        q.task_done()
      q.task_done()

    q = multiprocessing.JoinableQueue()
    procs = []
    for i in range(command_line.options.num_procs):
      procs.append(multiprocessing.Process(target=worker))
      procs[-1].daemon = True
      procs[-1].start()

    for file in files:
      if command_line.options.output_dir is not None:
        int_pickle_path = os.path.join(command_line.options.output_dir, \
          command_line.options.output_basename + os.path.basename(file))
        if command_line.options.skip_processed and os.path.isfile(int_pickle_path):
          print(file, "already processed, skipping")
          continue
        args = copy.copy(arguments)
        args.append("indexing.completeness_pickle=%s"%int_pickle_path)
      else:
        args = arguments

      q.put((file, args, ({'display':display})))

    q.join()

    for p in procs:
      q.put( None )

    q.join()

    for p in procs:
      p.join()

    print("Finished everything....")
    print("num active children:", len(multiprocessing.active_children()))
