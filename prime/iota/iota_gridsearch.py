from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 10/10/2014
Last Changed: 02/19/2015
Description : Runs cxi.index with signal/spot height and area grid search
'''

import os
import logging

from xfel.clustering.singleframe import SingleFrame
from subprocess import check_output
import dials.util.command_line as cmd

import sys
from cStringIO import StringIO


#logging = logging.getLogger('gs_log')

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

def integrate_one_image(mp_entry, n_int, log_dir, gs_params):
  """ Indexes and integrates a single image by running run_one_index from the
    xfel.cxi.display_spots module.

  Output:
  1. integrated pickle in correct output folder
  2. integration logfile (from captured cctbx.xfel output, one per image)
  3. list entry for integrated, not integrated or incorrectly integrated file

  Parameters:
  1. mp_entry: raw integration image in *.pickle format
  2. log_dir: main folder for logs
  3. gs_params: PHIL object containing grid-search and other parameters (from
                user input file)
  """

  from xfel.cxi.display_spots import run_one_index

  current_img = mp_entry[0]
  sig_height = mp_entry[1]
  spot_height = mp_entry[2]
  spot_area = mp_entry[3]
  target = gs_params.target

# generate filenames, etc.
  path = os.path.dirname(current_img)
  img_filename = os.path.basename(current_img)

  if os.path.relpath(path, os.path.abspath(gs_params.input)) == '.':
    input_dir = os.path.abspath(gs_params.input)
    output_dir = os.path.abspath(gs_params.output)
  else:
    input_dir = '{0}/{1}'.format(os.path.abspath(gs_params.input),
                                os.path.relpath(path,
                                os.path.abspath(gs_params.input)))
    output_dir = '{0}/{1}'.format(os.path.abspath(gs_params.output),
                                  os.path.relpath(path,
                                  os.path.abspath(gs_params.input)))

  current_output_dir = "{0}/tmp_{1}".format(output_dir, img_filename.split('.')[0])
  current_file = os.path.normpath("{0}/int_h{1}_a{2}_{3}".format(current_output_dir,
                               sig_height, spot_area, img_filename))
  current_log_file = os.path.normpath("{0}/{1}.log".format(log_dir,
                                      img_filename.split('.')[0]))


  arguments = ["target={0}".format(target),
               "distl.minimum_signal_height={0}".format(str(sig_height)),
               "distl.minimum_spot_height={0}".format(str(spot_height)),
               "distl.minimum_spot_area={0}".format(str(spot_area)),
               "indexing.completeness_pickle={0}".format(current_file)
               ]

  # Run indexing / integration here
  with Capturing() as index_log:
    try:
      display = False
      run_one_index(current_img, *arguments, **({'display':display}))

    except Exception, e:
      if hasattr(e, "classname"):
        print e.classname, "for %s:"%file,
      else:
        print "Indexing error for %s:"%file,
      print e

  int_success = False
  # generate, display and log grid search output
  if os.path.isfile(current_file):
    int_success = True
    observations = SingleFrame(current_file,
                               os.path.split(current_file)[1]).miller_array
    pickle_res = observations.d_max_min()
    pg = observations.space_group_info()
    unit_cell = observations.unit_cell().parameters()
    int_status = 'integrated with res = {:>6.2f} - {:<5.2f}, ' \
                  's.g.: {:^{wsg}}, u.c.: {:>6.2f}, {:>6.2f}, {:>6.2f}, ' \
                  '{:>6.2f}, {:>6.2f}, {:>6.2f}'.format(pickle_res[0],
                  pickle_res[1], pg, unit_cell[0], unit_cell[1],
                  unit_cell[2], unit_cell[3], unit_cell[4], unit_cell[5],
                  wsg = len(str(pg)))
  else:
    int_status = "not integrated"

  grid_search_output = '{:^{width}}: h = {:<3}, ' \
                       'a = {:<3} ---> {}'.format(img_filename,
                        sig_height, spot_area, int_status,
                        width = len(img_filename) + 2)

  logging.info(grid_search_output)

  # write integration logfile
  with open(current_log_file, 'a') as index_logfile:
    index_logfile.write("{:-^100}\n".format(" INTEGRATION: "\
                        "H={:>2}, A={:>2} ".format(spot_height, spot_area)))
    for item in index_log:
      index_logfile.write("{}\n".format(item))

    index_logfile.write("\n{}\n\n".format(grid_search_output))


  # Progress bar for integration
  with (open ('{0}/logs/progress.log'.format(gs_params.output), 'a')) as prog_log:
    prog_log.write("{0}\n".format(grid_search_output))

  with (open ('{0}/logs/progress.log'.format(gs_params.output), 'r')) as prog_log:
    prog_content = prog_log.read()
    prog_count = len(prog_content.splitlines())

  gs_prog = cmd.ProgressBar(title='GRID SEARCH', estimate_time=False, spinner=False)
  if prog_count >= n_int:
    gs_prog.finished()
  else:
    prog_step = 100 / n_int
    gs_prog.update(prog_count * prog_step)


def index_grid_search(mp_entry, n_int, log_dir, gs_params):
  """ Performs grid search of spotfinding parameters. Indexes / integrates w/o
  outputting a pickle file. Selects

  Output:
  1. integrated pickle in correct output folder
  2. integration logfile (from captured cctbx.xfel output, one per image)
  3. list entry for integrated, not integrated or incorrectly integrated file

  Parameters:
  1. mp_entry: raw integration image in *.pickle format
  2. log_dir: main folder for logs
  3. gs_params: PHIL object containing grid-search and other parameters (from
                user input file)
  """

  from xfel.cxi.display_spots import run_one_index

  current_img = mp_entry[0]
  sig_height = mp_entry[1]
  spot_height = mp_entry[2]
  spot_area = mp_entry[3]
  target = gs_params.target

# generate filenames, etc.
  path = os.path.dirname(current_img)
  img_filename = os.path.basename(current_img)

  if os.path.relpath(path, os.path.abspath(gs_params.input)) == '.':
    input_dir = os.path.abspath(gs_params.input)
    output_dir = os.path.abspath(gs_params.output)
  else:
    input_dir = '{0}/{1}'.format(os.path.abspath(gs_params.input),
                                os.path.relpath(path,
                                os.path.abspath(gs_params.input)))
    output_dir = '{0}/{1}'.format(os.path.abspath(gs_params.output),
                                  os.path.relpath(path,
                                  os.path.abspath(gs_params.input)))

  current_output_dir = "{0}/tmp_{1}".format(output_dir, img_filename.split('.')[0])
  current_file = os.path.normpath("{0}/int_h{1}_a{2}_{3}".format(current_output_dir,
                               sig_height, spot_area, img_filename))
  current_log_file = os.path.normpath("{0}/{1}.log".format(log_dir,
                                      img_filename.split('.')[0]))


  arguments = ["target={0}".format(target),
               "distl.minimum_signal_height={0}".format(str(sig_height)),
               "distl.minimum_spot_height={0}".format(str(spot_height)),
               "distl.minimum_spot_area={0}".format(str(spot_area)),
               "indexing.completeness_pickle={0}".format(current_file)
               ]

  # Run indexing / integration here
  with Capturing() as index_log:
    try:
      display = False
      run_one_index(current_img, *arguments, **({'display':display}))

    except Exception, e:
      if hasattr(e, "classname"):
        print e.classname, "for %s:"%file,
      else:
        print "Indexing error for %s:"%file,
      print e

  int_success = False
  # generate, display and log grid search output
  if os.path.isfile(current_file):
    int_success = True
    observations = SingleFrame(current_file,
                               os.path.split(current_file)[1]).miller_array
    pickle_res = observations.d_max_min()
    pg = observations.space_group_info()
    unit_cell = observations.unit_cell().parameters()
    int_status = 'integrated with res = {:>6.2f} - {:<5.2f}, ' \
                  's.g.: {:^{wsg}}, u.c.: {:>6.2f}, {:>6.2f}, {:>6.2f}, ' \
                  '{:>6.2f}, {:>6.2f}, {:>6.2f}'.format(pickle_res[0],
                  pickle_res[1], pg, unit_cell[0], unit_cell[1],
                  unit_cell[2], unit_cell[3], unit_cell[4], unit_cell[5],
                  wsg = len(str(pg)))
  else:
    int_status = "not integrated"

  grid_search_output = '{:^{width}}: h = {:<3}, ' \
                       'a = {:<3} ---> {}'.format(img_filename,
                        sig_height, spot_area, int_status,
                        width = len(img_filename) + 2)

  logging.info(grid_search_output)

  # write integration logfile
  with open(current_log_file, 'a') as index_logfile:
    index_logfile.write("{:-^100}\n".format(" INTEGRATION: "\
                        "H={:>2}, A={:>2} ".format(spot_height, spot_area)))
    for item in index_log:
      index_logfile.write("{}\n".format(item))

    index_logfile.write("\n{}\n\n".format(grid_search_output))


  # Progress bar for integration
  with (open ('{0}/logs/progress.log'.format(gs_params.output), 'a')) as prog_log:
    prog_log.write("{0}\n".format(grid_search_output))

  with (open ('{0}/logs/progress.log'.format(gs_params.output), 'r')) as prog_log:
    prog_content = prog_log.read()
    prog_count = len(prog_content.splitlines())

  gs_prog = cmd.ProgressBar(title='GRID SEARCH', estimate_time=False, spinner=False)
  if prog_count >= n_int:
    gs_prog.finished()
  else:
    prog_step = 100 / n_int
    gs_prog.update(prog_count * prog_step)



# ============================================================================ #



# Indexing and integration w/ grid search of spotfinding parameters
def index_integrate (mp_entry, n_int, log_dir, gs_params):
  """ DEPRECATED. Indexes and integrates a single image by running cxi.index

  Output:
  1. integrated pickle in correct output folder
  2. integration logfile (from captured cctbx.xfel output)
  3. list entry for integrated, not integrated or incorrectly integrated file

  Parameters:
  1. mp_entry: raw integration image in *.pickle format
  2. log_dir: main folder for logs
  3. gs_params: PHIL object containing grid-search and other parameters (from
                user input file)
  """

  current_img = mp_entry[0]
  sig_height = mp_entry[1]
  spot_height = mp_entry[2]
  spot_area = mp_entry[3]

  # generate filenames, etc.
  path = os.path.dirname(current_img)
  img_filename = os.path.basename(current_img)

  if os.path.relpath(path, os.path.abspath(gs_params.input)) == '.':
    input_dir = os.path.abspath(gs_params.input)
    output_dir = os.path.abspath(gs_params.output)
  else:
    input_dir = '{0}/{1}'.format(os.path.abspath(gs_params.input),
                                os.path.relpath(path,
                                os.path.abspath(gs_params.input)))
    output_dir = '{0}/{1}'.format(os.path.abspath(gs_params.output),
                                  os.path.relpath(path,
                                  os.path.abspath(gs_params.input)))

  current_output_dir = "{0}/tmp_{1}".format(output_dir, img_filename.split('.')[0])
  index_log_dir = "{0}/tmp_{1}".format(log_dir, img_filename.split('.')[0])


  # Make directories for output / log file for the image being integrated
  #if not os.path.exists(current_output_dir):
  #  os.makedirs(current_output_dir)
  #if not os.path.exists(index_log_dir):
  #  os.makedirs(index_log_dir)

  # Indexing / integration against the grid of spotfinding parameters
  cmds = ["cxi.index",
          "-d",
          "-o" + current_output_dir,
          "-b" + "int_h{0}_a{1}_".format(sig_height, spot_area),
          'target={0}'.format(gs_params.target),
          "distl.minimum_signal_height={0}".format(str(sig_height)),
          "distl.minimum_spot_height={0}".format(str(spot_height)),
          "distl.minimum_spot_area={0}".format(str(spot_area)),
          current_img]

  output = check_output(cmds)

  # write integration logfile
  current_log_file = "{0}/int_h{1}_a{2}_{3}.log".format(index_log_dir,
                     sig_height, spot_area, img_filename.split('.')[0])

  f = open(current_log_file, 'w')
  f.write(output)
  f.close()

  # generate, display and log grid search output
  current_file = "{0}/int_h{1}_a{2}_{3}".format(current_output_dir,
                                                sig_height, spot_area,
                                                img_filename)
  if os.path.isfile(current_file):
    int_success = True
    observations = SingleFrame(current_file,
                               os.path.split(current_file)[1]).miller_array
    pickle_res = observations.d_max_min()
    pg = observations.space_group_info()
    unit_cell = observations.unit_cell().parameters()
    int_status = 'integrated with res = {:>6.2f} - {:<5.2f}, ' \
                  's.g.: {:^{wsg}}, u.c.: {:>6.2f}, {:>6.2f}, {:>6.2f}, ' \
                  '{:>6.2f}, {:>6.2f}, {:>6.2f}'.format(pickle_res[0],
                  pickle_res[1], pg, unit_cell[0], unit_cell[1],
                  unit_cell[2], unit_cell[3], unit_cell[4], unit_cell[5],
                  wsg = len(str(pg)))
  else:
    int_status = "not integrated"

  grid_search_output = '{:^{width}}: h = {:<3}, ' \
                       'a = {:<3} ---> {}'.format(current_img,
                        sig_height, spot_area, int_status,
                        width = len(current_img) + 2)
  logging.info(grid_search_output)



def debug_integrate_one(mp_entry, log_dir, gs_params):

  from xfel.cxi.display_spots import run_one_index

  current_img = mp_entry[0]
  sig_height = mp_entry[1]
  spot_height = mp_entry[2]
  spot_area = mp_entry[3]
  target = gs_params.target

# generate filenames, etc.
  path = os.path.dirname(current_img)
  img_filename = os.path.basename(current_img)

  if os.path.relpath(path, os.path.abspath(gs_params.input)) == '.':
    input_dir = os.path.abspath(gs_params.input)
    output_dir = os.path.abspath(gs_params.output)
  else:
    input_dir = '{0}/{1}'.format(os.path.abspath(gs_params.input),
                                os.path.relpath(path,
                                os.path.abspath(gs_params.input)))
    output_dir = '{0}/{1}'.format(os.path.abspath(gs_params.output),
                                  os.path.relpath(path,
                                  os.path.abspath(gs_params.input)))

  current_output_dir = "{0}/tmp_{1}".format(output_dir, img_filename.split('.')[0])
  current_file = os.path.normpath("{0}/int_h{1}_a{2}_{3}".format(current_output_dir,
                               sig_height, spot_area, img_filename))
  current_log_file = os.path.normpath("{0}/{1}.log".format(log_dir,
                                      img_filename.split('.')[0]))


  from xfel.phil_preferences import load_cxi_phil
  from xfel.cxi.display_spots import run_one_index_core

  arguments = ["target={0}".format(target),
           "distl.minimum_signal_height={0}".format(str(sig_height)),
           "distl.minimum_spot_height={0}".format(str(spot_height)),
           "distl.minimum_spot_area={0}".format(str(spot_area)),
           ]
  try:
    display = False
    args = ["indexing.data={}".format(current_img),
          "beam_search_scope=0.5",
          "lepage_max_delta = 3.0",
          "spots_pickle = None",
          "subgroups_pickle = None",
          "refinements_pickle = None",
          "rmsd_tolerance = 5.0",
          "mosflm_rmsd_tolerance = 5.0",
          "difflimit_sigma_cutoff=2.0",
          "indexing.verbose_cv=True",
          "indexing.open_wx_viewer=%s"%display
          ] + list(arguments[1:])

    horizons_phil = load_cxi_phil(target, args)


    #import pdb; pdb.set_trace()
    info = run_one_index_core(horizons_phil)
    info.Files = info.organizer.Files
    info.phil_params = info.horizons_phil

#     print "\ninfo"
#     for item in dir(info): print item
#     print "\ninfo.organizer.info.keys"
#     for item in info.organizer.info.keys(): print item
#     print "Outliers: ", info.organizer.info["outlier_detection"]
#     print "Mosaicity: ", info.organizer.info['mosaicity']

  except Exception, e:
    if hasattr(e, "classname"):
      print e.classname, "for %s:"%file,
    else:
      print "Indexing error for %s:"%file,
    print e
