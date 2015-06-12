from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 06/04/2015
Last Changed: 06/04/2015
Description : Finds blank images using DISTL spotfinding. Outputs image filename
              only if it finds > 10 Bragg spots.
'''
import sys
import os
from cStringIO import StringIO

import spotfinder
from spotfinder.command_line.signal_strength import master_params as sf_params
from spotfinder.applications.wrappers import DistlOrganizer
import prime.iota.iota_input as inp


class Capturing(list):
  """ Class used to capture stdout from cctbx.xfel objects. Saves output in
  appendable list for potential logging.
  """
  def __enter__(self):
    self._stdout = sys.stdout
    sys.stdout = self._stringio = StringIO()
    return self
  def __exit__(self, *args):
    self.extend(self._stringio.getvalue().splitlines())
    sys.stdout = self._stdout

class Empty: pass

def spotfinding_param_search (input_img, gs_params):
  """ Performs a quick DISTL spotfinding without grid search.
      input:  input_img - image pickle
              gs_params - general parameters
      output: input_img - only returns this if image has >= 10 Bragg spots
  """

  params = sf_params.extract()
  params.distl.image = input_img

  E = Empty()
  E.argv=['Empty']
  E.argv.append(params.distl.image)

  selected_output = []
  total_output = []
  bragg_spots = []
  spotfinding_log = ['{}\n'.format(input_img)]

  # set spotfinding parameters for DISTL spotfinder
  params.distl.minimum_spot_area = gs_params.grid_search.a_avg
  params.distl.minimum_signal_height = gs_params.grid_search.h_avg
  params.distl.minimum_spot_height = int(gs_params.grid_search.h_avg / 2)

  # run DISTL spotfinder

  with Capturing() as junk_output:
    Org = DistlOrganizer(verbose = False, argument_module=E,
                         phil_params=params)

    Org.printSpots()

  # Extract relevant spotfinding info & make selection
  for frame in Org.S.images.keys():
    image = Org.S.images[frame]
    spots_inlier = Org.S.images[frame]['spots_inlier']
    Bragg_spots = Org.S.images[frame]['N_spots_inlier']

  if Bragg_spots >= 10:
    return input_img

def triage_mproc_wrapper(current_img):
  """ Multiprocessor wrapper for testing only
  """

  prog_count = file_list.index(current_img)
  n_int = len(file_list)

  gs_prog = cmd.ProgressBar(title='TRIAGE')
  if prog_count < n_int:
    prog_step = 100 / n_int
    gs_prog.update(prog_count * prog_step, prog_count)
  else:
    gs_prog.finished()

  return spotfinding_param_search(current_img, gs_params)

# ============================================================================ #

if __name__ == "__main__":

  from datetime import datetime
  from libtbx.easy_mp import parallel_map
  import prime.iota.iota_cmd as cmd

  filepath = sys.argv[1]
  if os.path.isdir(filepath):
    abs_inp_path = os.path.abspath(filepath)
  elif os.path.isfile(filepath):
    abs_inp_path = os.path.abspath(os.path.dirname(filepath))
  else:
    print "I don't get it. Please specify a filename or path to files."

  now = "{:%A, %b %d, %Y. %I:%M %p}".format(datetime.now())
  gs_params = inp.auto_mode(os.path.abspath(os.path.curdir),
                                  os.path.abspath(abs_inp_path), now)[0]
  logfile = os.path.join(os.curdir, 'triage.log')
  file_list = inp.make_input_list(gs_params)

  inp.main_log(logfile, '{:-^100} \n\n'.format(' IMAGE TRIAGE '))

  cmd.Command.start("Image triage")
  accepted_img = parallel_map(iterable=file_list,
                              func=triage_mproc_wrapper,
                              processes=32,
                              preserve_order = True)
  cmd.Command.end("Image triage -- DONE")

  fails = [f for f in file_list if f not in accepted_img]

  print "\nNo spots found in:"
  for i in fails: print i

