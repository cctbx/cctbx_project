from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 04/14/2014
Last Changed: 08/29/2018
Description : IOTA GUI Threads and PostEvents
'''

import os
import wx
from threading import Thread

from libtbx.easy_mp import parallel_map
from libtbx import easy_pickle as ep
from libtbx import easy_run

from dxtbx.datablock import DataBlockFactory
import multiprocessing

from xfel.clustering.cluster import Cluster
from cctbx.uctbx import unit_cell
from cctbx.sgtbx import lattice_symmetry
from cctbx import crystal

from iota.components.iota_utils import InputFinder
from iota.components.iota_misc import Capturing
import iota.components.iota_image as img
import iota.components.iota_misc as misc

ginp = InputFinder()

# -------------------------------- Threading --------------------------------- #

# Set up events for finishing one cycle and for finishing all cycles
tp_EVT_ALLDONE = wx.NewEventType()
EVT_ALLDONE = wx.PyEventBinder(tp_EVT_ALLDONE, 1)

tp_EVT_IMGDONE = wx.NewEventType()
EVT_IMGDONE = wx.PyEventBinder(tp_EVT_IMGDONE, 1)

tp_EVT_OBJDONE = wx.NewEventType()
EVT_OBJDONE = wx.PyEventBinder(tp_EVT_OBJDONE, 1)

class ImageFinderAllDone(wx.PyCommandEvent):
  ''' Send event when finished all cycles  '''
  def __init__(self, etype, eid, image_list=None):
    wx.PyCommandEvent.__init__(self, etype, eid)
    self.image_list = image_list
  def GetValue(self):
    return self.image_list

class ObjectFinderAllDone(wx.PyCommandEvent):
  ''' Send event when finished all cycles  '''
  def __init__(self, etype, eid, obj_list=None):
    wx.PyCommandEvent.__init__(self, etype, eid)
    self.obj_list = obj_list
  def GetValue(self):
    return self.obj_list

class AllDone(wx.PyCommandEvent):
  ''' Send event when finished all cycles  '''
  def __init__(self, etype, eid, img_objects=None):
    wx.PyCommandEvent.__init__(self, etype, eid)
    self.image_objects = img_objects
  def GetValue(self):
    return self.image_objects

class ProcessImage():
  ''' Wrapper class to do full processing of an image '''
  def __init__(self, init, input_entry, input_type = 'image', abort=False):
    self.init = init
    self.input_entry = input_entry
    self.input_type = input_type
    self.abort = abort

  def run(self):
    if self.abort:
      raise IOTATermination('IOTA: Run aborted by user')
    else:
      if self.input_type == 'image':
        img_object = img.SingleImage(self.input_entry, self.init)
        img_object.import_image()
      elif self.input_type == 'object':
        img_object = self.input_entry[2]
        img_object.import_int_file(self.init)
      else:
        img_object = None

      if self.init.params.image_conversion.convert_only:
        return img_object
      else:
        img_object.process()
        return img_object

class ProcThread(Thread):
  ''' Worker thread; generated so that the GUI does not lock up when
      processing is running '''
  def __init__(self,
               parent,
               init,
               iterable,
               term_file,
               input_type='image'):
    Thread.__init__(self)
    self.parent = parent
    self.init = init
    self.iterable = iterable
    self.type = input_type
    self.term_file = term_file
    self.aborted = False

  def run(self):
    try:
      img_objects = parallel_map(iterable=self.iterable,
                                 func = self.full_proc_wrapper,
                                 processes=self.init.params.n_processors)
    except IOTATermination, e:
      self.aborted = True
      print e
      return

    # Send "all done" event to GUI
    try:
      evt = AllDone(tp_EVT_ALLDONE, -1, img_objects=img_objects)
      wx.PostEvent(self.parent, evt)
    except Exception, e:
      pass

  def full_proc_wrapper(self, input_entry):
    abort = os.path.isfile(self.term_file)
    if abort:
      os.remove(self.term_file)
    try:
      proc_image_instance = ProcessImage(init=self.init,
                                         input_entry=input_entry,
                                         input_type=self.type,
                                         abort=abort)
      proc_image = proc_image_instance.run()
      return proc_image
    except IOTATermination, e:
      raise e
    except Exception, e:
      pass

class ImageFinderThread(Thread):
  ''' Worker thread generated to poll filesystem on timer. Will check to see
  if any new images have been found. Put on a thread to run in background '''
  def __init__(self,
               parent,
               image_paths,
               image_list):
    Thread.__init__(self)
    self.parent = parent
    self.image_paths = image_paths
    self.image_list = image_list

  def run(self):
    # Poll filesystem and determine which files are new (if any)

    ext_file_list = ginp.make_input_list(self.image_paths,
                                         filter=True,
                                         filter_type='image')
    old_file_list = [i[2] for i in self.image_list]
    new_file_list = [i for i in ext_file_list if i not in old_file_list]

    # Generate list of new images
    new_img = [[i, len(ext_file_list) + 1, j] for i, j in enumerate(
      new_file_list, len(old_file_list) + 1)]

    evt = ImageFinderAllDone(tp_EVT_IMGDONE, -1, image_list=new_img)
    wx.PostEvent(self.parent, evt)

class ObjectFinderThread(Thread):
  ''' Worker thread that polls filesystem on timer for image objects. Will
  collect and extract info on images processed so far'''
  def __init__(self,
               parent,
               object_folder,
               last_object=None,
               new_fin_base = None):
    Thread.__init__(self)
    self.parent = parent
    self.object_folder = object_folder
    self.new_fin_base = new_fin_base
    self.last_object = last_object

  def run(self):
    if self.last_object is not None:
      last = self.last_object.obj_file
    else:
      last = None
    object_files = ginp.get_file_list(self.object_folder,
                                      ext_only='int',
                                      last=last)
    new_objects = [self.read_object_file(i) for i in object_files]
    new_finished_objects = [i for i in new_objects if i is not None]

    evt = ObjectFinderAllDone(tp_EVT_OBJDONE, -1, obj_list=new_finished_objects)
    wx.PostEvent(self.parent, evt)

  def read_object_file(self, filepath):
    try:
      object = ep.load(filepath)
      return object
    except EOFError, e:
      print 'OBJECT_IMPORT_ERROR: ', e
      return None

class ImageViewerThread(Thread):
  ''' Worker thread that will move the image viewer launch away from the GUI
  and hopefully will prevent the image selection dialog freezing on MacOS'''
  def __init__(self,
               parent,
               file_string=None,
               viewer='dials.image_viewer',
               img_type=None):
    Thread.__init__(self)
    self.parent = parent
    self.file_string = file_string
    self.viewer = viewer
    self.img_type = img_type

  def run(self):
    command = '{} {}'.format(self.viewer, self.file_string)
    easy_run.fully_buffered(command)


# class IOTAUIThread(Thread):
#   ''' Main thread for IOTA UI; will contain all times and call all the other
#   threads - processing, object finding, etc. - separately; will use
#   PostEvents to send data to the main UI thread, which will plot only. The
#   idea is to prevent UI blocking as much as possible '''
#
#   def __init__(self,
#                parent,
#                gparams,
#                target_phil,
#                tmp_aborted_file=None):
#     Thread.__init__(self)
#     self.parent = parent
#
#     self.logtext = ''
#     self.obj_counter = 0
#     self.bookmark = 0
#     self.gparams = gparams
#     self.target_phil = target_phil
#     self.tmp_aborted_file = tmp_aborted_file
#
#     self.state = 'process'
#     self.recovery = False
#
#     self.abort_initiated = False
#     self.monitor_mode = False
#     self.monitor_mode_timeout = None
#     self.timeout_start = None
#     self.find_new_images = self.monitor_mode
#     self.start_object_finder = True
#
#     self.running_cluster = False
#     self.running_prime = False
#     self.draw_analysis = False
#
#     self.finished_objects = []
#     self.read_object_files = []
#     self.new_images = []
#     self.indices = []
#     self.cluster_info = None
#     self.pparams = None
#     self.prime_info = None
#
#     self.mxh = mx_handler()
#
#     # Timers
#     self.plt_timer = wx.Timer()
#     self.anl_timer = wx.Timer()
#
#     # Bindings
#     self.plt_timer.Bind(wx.EVT_TIMER, self.onPlotTimer)
#     self.anl_timer.Bind(wx.EVT_TIMER, self.onAnalysisTimer)
#
#   def run(self):
#     pass
#
#   def onPlotTimer(self, e):
#     ''' One second timer for status check and plotting '''
#
#     # Check for abort signal
#     if self.abort_initiated:
#       if self.img_process is not None:
#         self.run_aborted = self.img_process.aborted
#       elif self.gparams.mp_method == 'lsf':
#         info_command = 'bjobs -J {}'.format(self.job_id)
#         lsf_info = easy_run.fully_buffered(info_command).stdout_lines
#         self.run_aborted = (lsf_info == [])
#       else:
#         self.run_aborted = os.path.isfile(self.tmp_aborted_file)
#
#       if self.run_aborted:
#         self.finish_process()
#     else:
#       self.run_aborted = False
#
#     # Find processed image objects
#     if self.start_object_finder:
#       self.start_object_finder = False
#       if self.finished_objects is None or self.finished_objects == []:
#         self.last_object = None
#       else:
#         self.last_object = self.finished_objects[-1]
#       self.find_objects()
#
#     # Run an instance of new image finder on a separate thread
#     if self.find_new_images:
#       self.find_new_images = False
#       ext_image_list = self.img_list + self.new_images
#       img_finder = ImageFinderThread(self,
#                                      image_paths=self.gparams.input,
#                                      image_list=ext_image_list)
#       img_finder.start()
#
#     # Check if all images have been looked at; if yes, finish process
#     if self.obj_counter >= len(self.img_list):
#       if self.monitor_mode:
#         if len(self.new_images) > 0:
#           self.status_txt.SetLabel('Found {} new images'.format(len(self.new_images)))
#           self.timeout_start = None
#           self.state = 'new images'
#           self.process_images()
#         else:
#           if self.monitor_mode_timeout != None:
#             if self.timeout_start is None:
#               self.timeout_start = time.time()
#             else:
#               interval = time.time() - self.timeout_start
#               if interval >= self.monitor_mode_timeout:
#                 self.status_txt.SetLabel('Timed out. Finishing...')
#                 self.finish_process()
#               else:
#                 timeout_msg = 'No images found! Timing out in {} seconds' \
#                             ''.format(int(self.monitor_mode_timeout - interval))
#                 self.status_txt.SetLabel(timeout_msg)
#           else:
#             self.find_new_images = self.monitor_mode
#             self.status_txt.SetLabel('No new images found! Waiting ...')
#       else:
#         self.status_txt.SetLabel('Wrapping up ...')
#         self.finish_process()
#
#   def onAnalysisTimer(self, e):
#     self.plot_live_analysis()
#     if not self.running_cluster:
#       self.run_clustering_thread()
#
#   def run_clustering_thread(self):
#     # Run clustering
#     if (self.finished_objects is not None or self.finished_objects != []):
#       if self.proc_nb.GetSelection() == 2:
#         iterable = []
#         for obj in self.finished_objects:
#           try:
#             fin = obj.final
#             iterable.append([float(fin['a']),
#                              float(fin['b']),
#                              float(fin['c']),
#                              float(fin['alpha']),
#                              float(fin['beta']),
#                              float(fin['gamma']),
#                              fin['sg']
#                              ])
#           except Exception:
#             pass
#         self.running_cluster = True
#         cl_thread = thr.ClusterThread(self, iterable=iterable)
#         cl_thread.start()
#
#   def onFinishedCluster(self, e):
#     self.cluster_info = e.GetValue()
#     self.running_cluster = False
#
#     # Output cluster results
#     cluster_info_file = os.path.join(self.init.int_base, 'cluster_info.pickle')
#     ep.dump(cluster_info_file, obj=self.cluster_info)
#
#     if not self.running_prime:
#       self.pparams = None
#       self.run_prime_thread()
#
#   def run_prime_thread(self):
#     # Run PRIME (basic merge only)
#     if (self.finished_objects is not None or self.finished_objects != []):
#       if self.proc_nb.GetSelection() == 2:
#         # Collect list of final integrated pickles and write to file
#         final_files = [o.final['final'] for o in self.finished_objects if
#                        (o.final['final'] is not None and
#                         os.path.isfile(o.final['final']))]
#         final_list_file = os.path.join(self.init.int_base, 'finished_pickles.lst')
#         with open(final_list_file, 'w') as ff:
#           ff.write('\n'.join(final_files))
#
#         # make PRIME input file
#         if self.cluster_info is not None:
#           cl_sorted = sorted(self.cluster_info, key=lambda i: i['number'],
#                              reverse=True)
#           best_pg = cl_sorted[0]['pg'].split('/')[0]
#           best_uc = cl_sorted[0]['uc']
#
#           analyzer = Analyzer(init=self.init,
#                               all_objects=self.finished_objects,
#                               gui_mode=False)
#           analyzer.prime_data_path = final_list_file
#           analyzer.cons_pg = best_pg
#           analyzer.cons_uc = best_uc
#
#           prime_phil = analyzer.make_prime_input(filename='live_prime.phil',
#                                        run_zero=True)
#           self.pparams = prime_phil.extract()
#
#           # Modify specific options based in IOTA settings
#           # Queue options
#           if (
#                   self.init.params.mp_method == 'lsf' and
#                   self.init.params.mp_queue is not None
#           ):
#             self.pparams.queue.mode = 'bsub'
#             self.pparams.queue.qname = self.init.params.mp_queue
#
#           # Number of processors (automatically, 1/2 of IOTA procs)
#           self.pparams.n_processors = int(self.init.params.n_processors / 2)
#
#           # Generate command args
#           cmd_args_list = ['n_postref_cycle=0',
#                            'queue.mode={}'.format(self.pparams.queue.mode),
#                            'queue.qname={}'.format(self.pparams.queue.qname),
#                            'n_processors={}'.format(self.pparams.n_processors)
#                            ]
#           cmd_args = ' '.join(cmd_args_list)
#
#           # remove previous run to avoid conflict
#           prime_dir = os.path.join(self.init.int_base, 'prime/000')
#           if os.path.isdir(prime_dir):
#             shutil.rmtree(prime_dir)
#
#           # Launch PRIME
#           self.running_prime = True
#           out_file = os.path.join(prime_dir, 'log.txt')
#           prime_file = os.path.join(self.init.int_base, 'live_prime.phil')
#           prime_thread = pthr.PRIMEThread(self,
#                                           prime_file=prime_file,
#                                           out_file=out_file,
#                                           cmd_args=cmd_args,
#                                           signal_finished=True)
#           prime_thread.start()
#
#   def onFinishedPRIME(self, e):
#     self.running_prime = False
#     if self.pparams is not None:
#       self.get_prime_stats()
#
#   def get_prime_stats(self):
#     stats_folder = os.path.join(self.pparams.run_no, 'stats')
#     if os.path.isdir(stats_folder):
#       stat_files = [os.path.join(stats_folder, i) for i in
#                     os.listdir(stats_folder) if i.endswith('stat')]
#       if stat_files != []:
#         assert len(stat_files) == 1
#         stat_file = stat_files[0]
#         if os.path.isfile(stat_file):
#           self.prime_info = ep.load(stat_file)
#           live_prime_info_file = os.path.join(self.init.int_base,
#                                               'life_prime_info.pickle')
#           shutil.copyfile(stat_file, live_prime_info_file)
#
#   def process_images(self):
#     ''' One-fell-swoop importing / triaging / integration of images '''
#
#     # Set font properties for status window
#     font = self.sb.GetFont()
#     font.SetWeight(wx.NORMAL)
#     self.status_txt.SetFont(font)
#     self.status_txt.SetForegroundColour('black')
#
#
#     if self.init.params.cctbx.selection.select_only.flag_on:
#       self.img_list = [[i, len(self.init.gs_img_objects) + 1, j] for
#                        i, j in enumerate(self.init.gs_img_objects, 1)]
#       iterable = self.img_list
#       self.status_summary = [0] * len(self.img_list)
#       self.nref_list = [0] * len(self.img_list)
#       self.nref_xaxis = [i[0] for i in self.img_list]
#       self.res_list = [0] * len(self.img_list)
#       type = 'object'
#       self.status_txt.SetLabel('Re-running selection...')
#     else:
#       type = 'image'
#       if self.state == 'new images':
#         iterable = self.new_images
#         self.img_list.extend(self.new_images)
#         self.new_images = []
#         self.status_summary.extend([0] * len(iterable))
#         self.nref_list.extend([0] * len(iterable))
#         self.nref_xaxis.extend([i[0] for i in iterable])
#         self.res_list.extend([0] * len(iterable))
#         self.status_txt.SetForegroundColour('black')
#         self.status_txt.SetLabel('Processing additional {} images ({} total)...'
#                                  ''.format(len(iterable), len(self.img_list)))
#         self.plot_integration()
#       elif self.state == 'resume':
#         iterable = self.new_images
#         self.img_list.extend(self.new_images)
#         self.nref_list.extend([0] * len(self.new_images))
#         self.nref_xaxis.extend([i[0] for i in self.new_images])
#         self.res_list.extend([0] * len(self.new_images))
#         self.new_images = []
#         self.status_txt.SetLabel('Processing {} remaining images ({} total)...'
#                                  ''.format(len(iterable), len(self.img_list)))
#         self.start_object_finder = True
#       else:
#         self.img_list = [[i, len(self.init.input_list) + 1, j] for
#                          i, j in enumerate(self.init.input_list, 1)]
#         iterable = self.img_list
#         self.status_summary = [0] * len(self.img_list)
#         self.nref_list = [0] * len(self.img_list)
#         self.nref_xaxis = [i[0] for i in self.img_list]
#         self.res_list = [0] * len(self.img_list)
#         self.status_txt.SetLabel('Processing {} images...'
#                                  ''.format(len(self.img_list)))
#     self.gauge_process.SetRange(len(self.img_list))
#
#     iter_path = os.path.join(self.init.int_base, 'iter.cfg')
#     init_path = os.path.join(self.init.int_base, 'init.cfg')
#     ep.dump(iter_path, iterable)
#     ep.dump(init_path, self.init)
#
#     if self.gparams.mp_method == 'multiprocessing':
#       self.img_process = thr.ProcThread(self,
#                                         init=self.init,
#                                         iterable=iterable,
#                                         input_type=type,
#                                         term_file=self.tmp_abort_file)
#       self.img_process.start()
#     else:
#       self.img_process = None
#       self.job_id = None
#       queue = self.gparams.mp_queue
#       nproc = self.init.params.n_processors
#
#       if self.init.params.mp_method == 'lsf':
#         logfile = os.path.join(self.init.int_base, 'bsub.log')
#         pid = os.getpid()
#         try:
#           user = os.getlogin()
#         except OSError:
#           user = 'iota'
#         self.job_id = 'J_{}{}'.format(user[0], pid)
#         command = 'bsub -o {} -q {} -n {} -J {} ' \
#                   'iota.process {} --files {} --type {} --stopfile {}' \
#                   ''.format(logfile, queue, nproc, self.job_id,
#                             init_path, iter_path, type, self.tmp_abort_file)
#       elif self.init.params.mp_method == 'torq':
#         params = '{} --files {} --type {} --stopfile {}' \
#                  ''.format(init_path, iter_path, type, self.tmp_abort_file)
#         command = 'qsub -e /dev/null -o /dev/null -d {} iota.process -F "{}"' \
#                   ''.format(self.init.params.output, params)
#       else:
#         command = None
#       if command is not None:
#         try:
#           print command
#           easy_run.fully_buffered(command, join_stdout_stderr=True).show_stdout()
#           print 'JOB NAME = ', self.job_id
#         except thr.IOTATermination, e:
#           print 'IOTA: JOB TERMINATED',  e
#       else:
#         print 'IOTA ERROR: COMMAND NOT ISSUED!'
#         return
#
#   def analyze_results(self, analysis=None):
#     if len(self.final_objects) == 0:
#       self.status_txt.SetForegroundColour('red')
#       self.status_txt.SetLabel('No images successfully integrated')
#
#     else:
#       if not self.gparams.image_conversion.convert_only:
#         self.status_txt.SetForegroundColour('black')
#         self.status_txt.SetLabel('Analyzing results...')
#
#         # Do analysis
#         if analysis is None:
#           self.recovery = False
#           analysis = Analyzer(self.init,
#                               self.finished_objects,
#                               gui_mode=True)
#         plot = Plotter(self.gparams,
#                        self.final_objects,
#                        self.init.viz_base)
#
#         # Initialize summary tab
#         prime_file = os.path.join(self.init.int_base,
#                                   '{}.phil'.format(self.gparams.advanced.prime_prefix))
#         self.summary_tab = SummaryTab(self.proc_nb,
#                                       init=self.init,
#                                       gparams=self.gparams,
#                                       final_objects=self.final_objects,
#                                       out_dir=os.path.dirname(prime_file),
#                                       plot=plot)
#
#         # Run information
#         self.summary_tab.title_txt.SetLabel(noneset(self.gparams.description))
#         self.summary_tab.folder_txt.SetLabel(self.gparams.output)
#
#         # Analysis of integration
#         if self.gparams.advanced.integrate_with == 'cctbx':
#           self.summary_tab.sih_min.SetLabel("{:4.0f}".format(np.min(analysis.s)))
#           self.summary_tab.sih_max.SetLabel("{:4.0f}".format(np.max(analysis.s)))
#           self.summary_tab.sih_avg.SetLabel("{:4.2f}".format(np.mean(analysis.s)))
#           self.summary_tab.sih_std.SetLabel("{:4.2f}".format(np.std(analysis.s)))
#           self.summary_tab.sph_min.SetLabel("{:4.0f}".format(np.min(analysis.h)))
#           self.summary_tab.sph_max.SetLabel("{:4.0f}".format(np.max(analysis.h)))
#           self.summary_tab.sph_avg.SetLabel("{:4.2f}".format(np.mean(analysis.h)))
#           self.summary_tab.sph_std.SetLabel("{:4.2f}".format(np.std(analysis.h)))
#           self.summary_tab.spa_min.SetLabel("{:4.0f}".format(np.min(analysis.a)))
#           self.summary_tab.spa_max.SetLabel("{:4.0f}".format(np.max(analysis.a)))
#           self.summary_tab.spa_avg.SetLabel("{:4.2f}".format(np.mean(analysis.a)))
#           self.summary_tab.spa_std.SetLabel("{:4.2f}".format(np.std(analysis.a)))
#
#         # Dataset information
#         if self.recovery:
#           if hasattr(analysis, 'clusters'):
#             clusters = analysis.clusters
#             pg = clusters[0]['pg']
#             uc = clusters[0]['uc']
#           else:
#             clusters = []
#             if hasattr(analysis, 'cons_pg'):
#               pg = analysis.cons_pg
#             else:
#               pg = None
#             if hasattr(analysis, 'cons_uc'):
#               uc = analysis.cons_uc
#             else:
#               uc = None
#         else:
#           analysis.print_results()
#           clusters = analysis.unit_cell_analysis()
#
#           print 'DEBUG: CLUSTERS = ', clusters
#           pg = clusters[0]['pg']
#           uc = clusters[0]['uc']
#
#         if clusters != []:
#           self.summary_tab.report_clustering_results(clusters=clusters)
#
#         if pg != None:
#           self.summary_tab.pg_txt.SetLabel(str(pg))
#         if uc is not None:
#           if type(uc) is str:
#             unit_cell = uc
#           else:
#             unit_cell = " ".join(['{:4.1f}'.format(i) for i in uc])
#         else:
#           unit_cell = ''
#
#         self.summary_tab.uc_txt.SetLabel(unit_cell)
#         res = to_unicode(u"{:4.2f} - {:4.2f} {}".format(np.mean(analysis.lres),
#                                   np.mean(analysis.hres),
#                                   u'\u212B'))
#         self.summary_tab.rs_txt.SetLabel(res)
#
#         if self.recovery and (
#                 hasattr(analysis, 'beamX_mm') and
#                 hasattr(analysis, 'beamY_mm') and
#                 hasattr(analysis, 'pixel_size')
#                 ):
#           beamX_mm = analysis.beamX_mm
#           beamY_mm = analysis.beamY_mm
#           pixel_size = analysis.pixel_size
#         else:
#           with warnings.catch_warnings():
#             # To catch any 'mean of empty slice' runtime warnings
#             warnings.simplefilter("ignore", category=RuntimeWarning)
#             beamxy_calc = plot.calculate_beam_xy()
#             beamX, beamY = beamxy_calc[:2]
#             pixel_size = beamxy_calc[-1]
#             beamX_mm = np.median(beamX)
#             beamY_mm = np.median(beamY)
#         beamX_px = beamX_mm / pixel_size
#         beamY_px = beamY_mm / pixel_size
#         beamXY = "X = {:4.1f} mm / {:4.0f} px\n" \
#                  "Y = {:4.1f} mm / {:4.0f} px" \
#                  "".format(beamX_mm, beamX_px, beamY_mm, beamY_px)
#         self.summary_tab.xy_txt.SetLabel(beamXY)
#
#         # Summary
#         if self.recovery:
#           self.summary_tab.readin_txt.SetLabel(str(analysis.n_all_objects))
#           self.summary_tab.nodiff_txt.SetLabel(str(analysis.n_no_diff_objects))
#           self.summary_tab.w_diff_txt.SetLabel(str(analysis.n_diff_objects))
#           if self.gparams.advanced.integrate_with == 'cctbx':
#             self.summary_tab.noint_txt.SetLabel(str(analysis.n_not_int_objects))
#             self.summary_tab.noprf_txt.SetLabel(str(analysis.n_filter_fail_objects))
#           elif self.gparams.advanced.integrate_with == 'dials':
#             self.summary_tab.nospf_txt.SetLabel(str(analysis.n_not_spf_objects))
#             self.summary_tab.noidx_txt.SetLabel(str(analysis.n_not_idx_objects))
#             self.summary_tab.noint_txt.SetLabel(str(analysis.n_not_int_objects))
#             self.summary_tab.noflt_txt.SetLabel(str(analysis.n_filter_fail_objects))
#           self.summary_tab.final_txt.SetLabel(str(analysis.n_final_objects))
#         else:
#           self.summary_tab.readin_txt.SetLabel(str(len(analysis.all_objects)))
#           self.summary_tab.nodiff_txt.SetLabel(str(len(analysis.no_diff_objects)))
#           self.summary_tab.w_diff_txt.SetLabel(str(len(analysis.diff_objects)))
#           if self.gparams.advanced.integrate_with == 'cctbx':
#             self.summary_tab.noint_txt.SetLabel(str(len(analysis.not_int_objects)))
#             self.summary_tab.noprf_txt.SetLabel(str(len(analysis.filter_fail_objects)))
#           elif self.gparams.advanced.integrate_with == 'dials':
#             self.summary_tab.nospf_txt.SetLabel(str(len(analysis.not_spf_objects)))
#             self.summary_tab.noidx_txt.SetLabel(str(len(analysis.not_idx_objects)))
#             self.summary_tab.noint_txt.SetLabel(str(len(analysis.not_int_objects)))
#             self.summary_tab.noflt_txt.SetLabel(str(len(analysis.filter_fail_objects)))
#           self.summary_tab.final_txt.SetLabel(str(len(analysis.final_objects)))
#
#           # Generate input file for PRIME
#           analysis.print_summary()
#           analysis.make_prime_input(filename=prime_file)
#
#         # Display summary
#         self.proc_nb.AddPage(self.summary_tab, 'Summary', select=True)
#         print 'DEBUG: ANALYSIS PAGE # = ', self.proc_nb.GetSelection()
#         # self.proc_nb.SetSelection(2)
#
#       # Signal end of run
#       font = self.sb.GetFont()
#       font.SetWeight(wx.BOLD)
#       self.status_txt.SetFont(font)
#       self.status_txt.SetForegroundColour('blue')
#       self.status_txt.SetLabel('DONE')
#
#     # Finish up
#     self.display_log()
#     self.plot_integration(force_plot=True)
#     self.plot_live_analysis(force_plot=True)
#
#     # Stop timer
#     self.timer.Stop()
#     self.chart_timer.Stop()
#
#   def find_objects(self, find_old=False):
#     if find_old:
#       min_back = None
#     else:
#       min_back = -1
#     object_files = ginp.get_file_list(self.init.obj_base,
#                                       ext_only='int',
#                                       min_back=min_back)
#     new_object_files = list(set(object_files) - set(self.read_object_files))
#     new_objects = [self.read_object_file(i) for i in new_object_files]
#     new_finished_objects = [i for i in new_objects if
#                             i is not None and i.status == 'final']
#
#     self.finished_objects.extend(new_finished_objects)
#     self.read_object_files = [i.obj_file for i in self.finished_objects]
#
#     self.populate_data_points(objects=new_finished_objects)
#
#     if str(self.state).lower() in ('finished', 'aborted', 'unknown'):
#       if self.finished_objects != []:
#         self.finish_process()
#       else:
#         return
#     else:
#       self.start_object_finder = True
#
#   def read_object_file(self, filepath):
#     try:
#       object = ep.load(filepath)
#       if object.final['final'] is not None:
#         pickle_path = object.final['final']
#         if os.path.isfile(pickle_path):
#           pickle = ep.load(pickle_path)
#           object.final['observations'] = pickle['observations'][0]
#       return object
#     except Exception, e:
#       print 'OBJECT_IMPORT_ERROR for {}: {}'.format(filepath, e)
#       return None
#
#   def populate_data_points(self, objects=None):
#     self.indices = []
#     self.b_factors = []
#     if objects is not None:
#       for obj in objects:
#         try:
#           self.nref_list[obj.img_index - 1] = obj.final['strong']
#           self.res_list[obj.img_index - 1] = obj.final['res']
#           if 'observations' in obj.final:
#             obs = obj.final['observations']
#             self.indices.extend([i[0] for i in obs])
#             try:
#               asu_contents = self.mxh.get_asu_contents(500)
#               observations_as_f = obs.as_amplitude_array()
#               observations_as_f.setup_binner(auto_binning=True)
#               wp = statistics.wilson_plot(observations_as_f, asu_contents,
#                                           e_statistics=True)
#               self.b_factors.append(wp.wilson_b)
#             except RuntimeError, e:
#               self.b_factors.append(0)
#         except Exception, e:
#           print 'OBJECT_ERROR:', e, "({})".format(obj.obj_file)
#           pass
#
#
#   def finish_process(self):
#     import shutil
#     self.timer.Stop()
#     self.chart_timer.Stop()
#
#     if self.finished_objects is None:
#       font = self.sb.GetFont()
#       font.SetWeight(wx.BOLD)
#       self.status_txt.SetFont(font)
#       self.status_txt.SetForegroundColour('blue')
#       self.status_txt.SetLabel('OBJECT READ-IN ERROR! NONE IMPORTED')
#       return
#
#     if str(self.state).lower() in ('finished', 'aborted', 'unknown'):
#       self.gauge_process.Hide()
#       font = self.sb.GetFont()
#       font.SetWeight(wx.BOLD)
#       self.status_txt.SetFont(font)
#       run_no = int(self.init.int_base.split('/')[-1])
#       self.status_txt.SetLabel('Run #{} Loaded!'.format(run_no))
#       self.proc_toolbar.EnableTool(self.tb_btn_abort.GetId(), False)
#       self.proc_toolbar.EnableTool(self.tb_btn_monitor.GetId(), False)
#       self.proc_toolbar.ToggleTool(self.tb_btn_monitor.GetId(), False)
#
#       analysis_file = os.path.join(self.init.int_base, 'analysis.pickle')
#       if os.path.isfile(analysis_file):
#         analysis = ep.load(analysis_file)
#       else:
#         analysis = None
#
#       if self.finished_objects == []:
#         # Get image processing data from finished objects
#         if analysis is not None and hasattr(analysis, 'image_objects'):
#           self.finished_objects = [i for i in analysis.image_objects if
#                                    i is not None and i.status == 'final']
#         else:
#           self.find_objects(find_old=True)
#
#         # Check for and recover clustering data
#         cluster_info_file = os.path.join(self.init.int_base,
#                                          'cluster_info.pickle')
#         if os.path.isfile(cluster_info_file):
#           self.cluster_info = ep.load(cluster_info_file)
#
#         # Check for and recover live PRIME results
#         live_prime_phil = os.path.join(self.init.int_base, 'live_prime.phil')
#         if os.path.isfile(live_prime_phil):
#           import iotbx.phil as ip
#           from prime.postrefine.mod_input import master_phil as prime_master_phil
#           with open(live_prime_phil, 'r') as lpf:
#             contents = lpf.read()
#           prime_phil = ip.parse(contents)
#           prime_phil = prime_master_phil.fetch(prime_phil)
#           self.pparams = prime_phil.extract()
#           live_prime_info_file = os.path.join(self.init.int_base,
#                                               'life_prime_info.pickle')
#           if os.path.isfile(live_prime_info_file):
#             self.prime_info = ep.load(live_prime_info_file)
#
#       self.populate_data_points(objects=self.finished_objects)
#
#       if str(self.state).lower() == 'finished':
#         self.final_objects = [i for i in self.finished_objects if i.fail is None]
#         self.analyze_results(analysis=analysis)
#
#       else:
#         if len(self.finished_objects) > 0:
#           self.plot_integration(force_plot=True)
#           self.plot_live_analysis(force_plot=True)
#         if os.path.isfile(os.path.join(self.init.int_base, 'init.cfg')):
#           self.proc_toolbar.EnableTool(self.tb_btn_resume.GetId(), True)
#
#       return
#
#     elif self.run_aborted:
#       self.gauge_process.Hide()
#       font = self.sb.GetFont()
#       font.SetWeight(wx.BOLD)
#       self.status_txt.SetFont(font)
#       self.status_txt.SetForegroundColour('red')
#       self.status_txt.SetLabel('ABORTED BY USER')
#       self.proc_toolbar.EnableTool(self.tb_btn_resume.GetId(), True)
#       try:
#         shutil.rmtree(self.init.tmp_base)
#       except Exception:
#         pass
#       print 'JOB TERMINATED!'
#       return
#     else:
#       self.final_objects = [i for i in self.finished_objects if i.fail == None]
#       self.gauge_process.Hide()
#       self.proc_toolbar.EnableTool(self.tb_btn_abort.GetId(), False)
#       self.proc_toolbar.EnableTool(self.tb_btn_monitor.GetId(), False)
#       self.proc_toolbar.ToggleTool(self.tb_btn_monitor.GetId(), False)
#       self.sb.SetStatusText('{} of {} images successfully integrated'\
#                             ''.format(len(self.final_objects), len(self.img_list)), 1)
#       if len(self.final_objects) > 0:
#         # Signal end of run
#         self.plot_integration(force_plot=True)
#         self.plot_live_analysis(force_plot=True)
#         self.analyze_results()
#       else:
#         font = self.sb.GetFont()
#         font.SetWeight(wx.BOLD)
#         self.status_txt.SetFont(font)
#         self.status_txt.SetForegroundColour('blue')
#         self.status_txt.SetLabel('NO IMAGES PROCESSED')
#
#       print 'JOB FINISHED!'
#
#       try:
#         shutil.rmtree(self.init.tmp_base)
#       except Exception:
#         pass

# ------------------------------ IMAGE TRACKING ------------------------------ #

tp_EVT_SPFDONE = wx.NewEventType()
EVT_SPFDONE = wx.PyEventBinder(tp_EVT_SPFDONE, 1)

tp_EVT_SPFALLDONE = wx.NewEventType()
EVT_SPFALLDONE = wx.PyEventBinder(tp_EVT_SPFALLDONE, 1)

tp_EVT_SLICEDONE = wx.NewEventType()
EVT_SLICEDONE = wx.PyEventBinder(tp_EVT_SLICEDONE)

tp_EVT_SPFTERM = wx.NewEventType()
EVT_SPFTERM = wx.PyEventBinder(tp_EVT_SPFTERM)

class IOTATermination(Exception):
  def __init__(self, termination):
    Exception.__init__(self, termination)

class SpotFinderAllDone(wx.PyCommandEvent):
  ''' Send event when finished all cycles  '''
  def __init__(self, etype, eid, info=None):
    wx.PyCommandEvent.__init__(self, etype, eid)
    self.info = info
  def GetValue(self):
    return self.info

class SpotFinderOneDone(wx.PyCommandEvent):
  ''' Send event when finished all cycles  '''
  def __init__(self, etype, eid, info=None):
    wx.PyCommandEvent.__init__(self, etype, eid)
    self.info = info
  def GetValue(self):
    return self.info

class SpotFinderTerminated(wx.PyCommandEvent):
  ''' Send event when spotfinder terminated '''
  def __init__(self, etype, eid):
    wx.PyCommandEvent.__init__(self, etype, eid)
  def GetValue(self):
    return None

class SpotFinderDIALSThread():
  def __init__(self, parent, processor, term_file,
               run_indexing=False, run_integration=False):
    self.meta_parent = parent.parent
    self.processor = processor
    self.term_file = term_file
    self.run_indexing = run_indexing
    self.run_integration = run_integration

  def run(self, idx, img):
    if self.meta_parent.terminated:
      raise IOTATermination('IOTA_TRACKER: SPF Termination signal received!')
    else:
      with Capturing() as junk_output:
        fail = False
        sg = None
        uc = None
        try:
          datablock = DataBlockFactory.from_filenames([img])[0]
          observed = self.processor.find_spots(datablock=datablock)
        except Exception, e:
          fail = True
          observed = []
          pass

        # TODO: Indexing / lattice determination very slow (how to speed up?)
        if self.run_indexing:
          if not fail:
            try:
              experiments, indexed = self.processor.index(
                datablock=datablock, reflections=observed)
            except Exception, e:
              fail = True
              pass

          if not fail:
            try:
              solution = self.processor.refine_bravais_settings(
                reflections=indexed, experiments=experiments)

              # Only reindex if higher-symmetry solution found
              if solution is not None:
                experiments, indexed = self.processor.reindex(
                  reflections=indexed,
                  experiments=experiments,
                  solution=solution)
              lat = experiments[0].crystal.get_space_group().info()
              sg = str(lat).replace(' ', '')
            except Exception:
              fail = True
              pass

          if not fail:
            unit_cell = experiments[0].crystal.get_unit_cell().parameters()
            uc = ' '.join(['{:.4f}'.format(i) for i in unit_cell])

          if self.run_integration:
            if not fail:
              try:
                # Run refinement
                experiments, indexed = self.processor.refine(
                  experiments=experiments,
                  centroids=indexed)
              except Exception, e:
                fail = True
                pass

            if not fail:
              try:
                print experiments
                print indexed
                integrated = self.processor.integrate(experiments=experiments,
                                                      indexed=indexed)
              except Exception, e:
                pass

      return [idx, int(len(observed)), img, sg, uc]




class SpotFinderMosflmThread():
  def __init__(self, parent, term_file):
    self.meta_parent = parent.parent
    self.term_file = term_file

  def run(self, idx, img):
    if os.path.isfile(self.term_file):
      raise IOTATermination('IOTA_TRACKER: Termination signal received!')
    else:
      # First, parse filepath to create Mosflm template
      directory = os.path.dirname(img)
      filepath = os.path.basename(img).split('.')
      fname = filepath[0]
      extension = filepath[1]
      if '_' in fname:
        suffix = fname.split('_')[-1]
      elif '-' in fname:
        suffix = fname.split('-')[-1]
      elif '.' in fname:
        suffix = fname.split('.')[-1]
      else:
        suffix = fname
      img_number = int(''.join(n if n.isdigit() else '' for n in suffix))
      prefix = fname.replace(suffix, '')
      n_suffix = ''.join("#" if c.isdigit() else c for c in suffix)
      template = '{}{}.{}'.format(prefix, n_suffix, extension)

      # Create autoindex.com w/ Mosflm script
      # Write to temporary file and change permissions to run
      autoindex = ['#! /bin/tcsh -fe',
                   'ipmosflm << eof-ipmosflm'.format(fname),
                   'NEWMATRIX {0}.mat'.format(fname),
                   'DIRECTORY {}'.format(directory),
                   'TEMPLATE {}'.format(template),
                   'AUTOINDEX DPS THRESH 0.1 IMAGE {} PHI 0 0.01'.format(
                     img_number),
                   'GO',
                   'eof-ipmosflm'
                   ]
      autoindex_string = '\n'.join(autoindex)
      autoindex_filename = 'autoindex_{}.com'.format(idx)

      with open(autoindex_filename, 'w') as af:
        af.write(autoindex_string)
      os.chmod(autoindex_filename, 0755)

      # Run Mosflm autoindexing
      command = './{}'.format(autoindex_filename)
      out = easy_run.fully_buffered(command, join_stdout_stderr=True)

      # Scrub text output
      final_spots = [l for l in out.stdout_lines if
                     'spots written for image' in l]
      final_cell_line = [l for l in out.stdout_lines if 'Final cell' in l]
      final_sg_line = [l for l in out.stdout_lines if 'space group' in l]

      if final_spots != []:
        spots = final_spots[0].rsplit()[0]
      else:
        spots = 0
      if final_cell_line != []:
        cell = final_cell_line[0].replace('Final cell (after refinement) is',
                                          '').rsplit()
      else:
        cell = None
      if final_sg_line != []:
        sg = final_sg_line[0].rsplit()[6]
      else:
        sg = None

      # Temp file cleanup
      try:
        os.remove('{}.mat'.format(fname))
      except Exception:
        pass
      try:
        os.remove('{}.spt'.format(prefix[:-1]))
      except Exception:
        pass
      try:
        os.remove('SUMMARY')
      except Exception:
        pass
      try:
        os.remove(autoindex_filename)
      except Exception:
        pass

      return [idx, spots, img, sg, cell]


class SpotFinderThread(Thread):
  ''' Basic spotfinder (with defaults) that could be used to rapidly analyze
  images as they are collected '''
  def __init__(self,
               parent,
               data_list=None,
               term_file=None,
               proc_params=None,
               backend='dials',
               n_proc=0,
               run_indexing=False,
               run_integration=False):
    Thread.__init__(self)
    self.parent = parent
    self.data_list = data_list
    self.term_file = term_file
    self.terminated = False
    self.backend = backend
    self.run_indexing = run_indexing
    self.run_integration = run_integration
    if n_proc > 0:
      self.n_proc = n_proc
    else:
      self.n_proc = multiprocessing.cpu_count() - 2

    if self.backend == 'dials':
      # Modify default DIALS parameters
      # These parameters will be set no matter what
      proc_params.output.datablock_filename = None
      proc_params.output.indexed_filename = None
      proc_params.output.strong_filename = None
      proc_params.output.refined_experiments_filename = None
      proc_params.output.integrated_filename = None
      proc_params.output.integrated_experiments_filename = None
      proc_params.output.profile_filename = None
      proc_params.output.integration_pickle = None

      from iota.components.iota_dials import IOTADialsProcessor
      self.processor = IOTADialsProcessor(params=proc_params)

  def run(self):
    try:
      parallel_map(iterable=self.data_list,
                   func=self.spf_wrapper,
                   callback=self.callback,
                   processes=self.n_proc)
    except IOTATermination, e:
      self.terminated = True
      print e

    # Signal that this batch is finished
    try:
      if self.terminated:
        print 'RUN TERMINATED!'
        evt = SpotFinderTerminated(tp_EVT_SPFTERM, -1)
        wx.PostEvent(self.parent, evt)

      wx.CallAfter(self.parent.onSpfAllDone, self.data_list)

      # info = self.data_list
      # evt = SpotFinderAllDone(tp_EVT_SPFALLDONE, -1, info=info)
      # wx.PostEvent(self.parent, evt)
      return
    except TypeError, e:
      print e
      return

  def spf_wrapper(self, img):
    try:
      if os.path.isfile(img):
        if self.backend == 'dials':
          spf_worker = SpotFinderDIALSThread(self,
                                             processor=self.processor,
                                             term_file=self.term_file,
                                             run_indexing=self.run_indexing,
                                             run_integration=self.run_integration
                                             )
          result = spf_worker.run(idx=int(self.data_list.index(img)), img=img)
        elif self.backend == 'mosflm':
          spf_worker = SpotFinderMosflmThread(self, self.term_file)
          result = spf_worker.run(idx=int(self.data_list.index(img)), img=img)
        else:
          result = [int(self.data_list.index(img)), 0, img, None, None]
        return result
      else:
        return [int(self.data_list.index(img)), 0, img, None, None]
    except IOTATermination, e:
      raise e

  def callback(self, info):
    try:
      wx.CallAfter(self.parent.onSpfOneDone, info)
      # evt = SpotFinderOneDone(tp_EVT_SPFDONE, -1, info=info)
      # wx.PostEvent(self.parent.parent, evt)
    except TypeError:
      pass

  # def terminate_thread(self):
  #   raise IOTATermination('IOTA_TRACKER: SPF THREAD Terminated!')


class InterceptorFileThread(Thread):
  def __init__(self,
               parent,
               results_file,
               reorder=False):
    Thread.__init__(self)
    self.parent = parent
    self.results_file = results_file
    self.reorder = reorder

    self.bookmark = 0
    self.msg = ''
    self.spotfinding_info = []
    self.cluster_info = None

    self.prc_timer = wx.Timer()
    self.cls_timer = wx.Timer()

    #Bindings
    self.prc_timer.Bind(wx.EVT_TIMER, self.onProcTimer)
    self.cls_timer.Bind(wx.EVT_TIMER, self.onClusterTimer)


  def run(self):
    # self.timer.Start(1000)
    pass

  def onProcTimer(self, e):
    if os.path.isfile(self.results_file):
      with open(self.results_file, 'r') as rf:
        rf.seek(self.bookmark)
        split_info = [i.replace('\n', '').split(' ') for i in rf.readlines()]
        self.bookmark = rf.tell()

      if self.reorder:
        idx_offset = len(self.spotfinding_info)
        new_info = [
          [split_info.index(i) + idx_offset,
           int(i[1]), i[2], i[3], tuple(i[4:10])] if len(i) > 5 else
          [split_info.index(i) + idx_offset,
           int(i[1]), i[2], misc.makenone(i[3]), misc.makenone(i[4])]
          for i in split_info]
      else:
        new_info = [
          [int(i[0]), int(i[1]), i[2], i[3], tuple(i[4:10])] if len(i) > 5 else
          [int(i[0]), int(i[1]), i[2], misc.makenone(i[3]), misc.makenone(i[4])]
          for i in split_info]

      if len(new_info) > 0:
        self.spotfinding_info.extend(new_info)

        if len(self.spotfinding_info) > 0:
          self.msg = 'Tracking new images in {} ...'.format(self.results_file)
      else:
        self.msg = 'Waiting for new images in {} ...'.format(self.results_file)

    else:
      self.msg = 'Waiting for new run to initiate...'

    info = [self.msg, self.spotfinding_info, self.cluster_info]
    evt = SpotFinderOneDone(tp_EVT_SPFDONE, -1, info=info)
    wx.PostEvent(self.parent, evt)

  def onClusterTimer(self, e):
    input = []
    for item in self.spotfinding_info:
      if item[4] is not None:
        try:
          if type(item[4]) in (tuple, list):
            uc = item[4]
          else:
            uc = item[4].rsplit()
          info_line = [float(i) for i in uc]
          info_line.append(item[3])
          input.append(info_line)
        except ValueError, e:
          print 'CLUSTER ERROR: ', e
          pass

    if len(input) > 0:
      self.running_clustering = True
      cluster_thread = ClusterWorkThread(self)
      self.cluster_info = cluster_thread.run(iterable=input)

  def terminate_thread(self):
    raise IOTATermination('IOTA_TRACKER: Termination signal received!')

class InterceptorThread(Thread):
  ''' Thread for the full Interceptor image processing process; will also
   house the processing timer, which will update the UI front end and initiate
   plotting '''
  def __init__(self,
               parent,
               data_folder=None,
               term_file=None,
               proc_params=None,
               backend='dials',
               n_proc=0,
               min_back=None,
               run_indexing=False,
               run_integration=False):
    Thread.__init__(self)
    self.parent = parent
    self.data_folder = data_folder
    self.term_file = term_file
    self.terminated = False
    self.backend = backend
    self.run_indexing = run_indexing
    self.run_integration = run_integration
    self.min_back = min_back
    self.submit_new_images = True

    self.spotfinding_info = []
    self.cluster_info = None
    self.msg = None
    self.done_list = []
    self.data_list = []
    self.new_data = []

    self.spf_thread = None

    if n_proc > 0:
      self.n_proc = n_proc
    else:
      self.n_proc = multiprocessing.cpu_count() - 2

    if self.backend == 'dials':
      # Modify default DIALS parameters
      # These parameters will be set no matter what
      proc_params.output.datablock_filename = None
      proc_params.output.indexed_filename = None
      proc_params.output.strong_filename = None
      proc_params.output.refined_experiments_filename = None
      proc_params.output.integrated_filename = None
      proc_params.output.integrated_experiments_filename = None
      proc_params.output.profile_filename = None
      proc_params.output.integration_pickle = None

      self.proc_params = proc_params

      # from iota.components.iota_dials import IOTADialsProcessor
      # self.processor = IOTADialsProcessor(params=proc_params)

    self.prc_timer = wx.Timer()
    self.cls_timer = wx.Timer()

    #Bindings
    self.prc_timer.Bind(wx.EVT_TIMER, self.onProcTimer)
    self.cls_timer.Bind(wx.EVT_TIMER, self.onClusterTimer)


  def run(self):
    pass

  def onProcTimer(self, e):
    ''' Main timer (1 sec) will send data to UI, find new images, and submit
    new processing run '''

    # Send current data to UI
    info = [self.msg, self.spotfinding_info, self.cluster_info]
    evt = SpotFinderOneDone(tp_EVT_SPFDONE, -1, info=info)
    wx.PostEvent(self.parent, evt)

    # Find new images
    if self.data_list != []:
      last_file = self.data_list[-1]
    else:
      last_file = None
    self.find_new_images(last_file=last_file, min_back=self.min_back)

    if self.spf_thread is not None:
      if not self.spf_thread.isAlive():
        self.submit_new_images = True

    # Submit new images (if found)
    if self.submit_new_images and len(self.new_data) > 0:
      self.submit_new_images = False
      self.run_processing()

  def onClusterTimer(self, e):
    input = []
    if len(self.spotfinding_info) > 0:
      for item in self.spotfinding_info:
        if item[4] is not None:
          try:
            if type(item[4]) in (tuple, list):
              uc = item[4]
            else:
              uc = item[4].rsplit()
            info_line = [float(i) for i in uc]
            info_line.append(item[3])
            input.append(info_line)
          except ValueError, e:
            print 'CLUSTER ERROR: ', e
            pass

      if len(input) > 0:
        self.running_clustering = True
        cluster_thread = ClusterWorkThread(self)
        self.cluster_info = cluster_thread.run(iterable=input)

  def find_new_images(self, min_back=None, last_file=None):
    found_files = ginp.make_input_list([self.data_folder],
                                       filter=True,
                                       filter_type='image',
                                       last=last_file,
                                       min_back=min_back)

    # Sometimes duplicate files are found anyway; clean that up
    found_files = list(set(found_files) - set(self.data_list))

    # Add new files to the data list & clean up
    self.new_data.extend(found_files)
    self.new_data = sorted(self.new_data, key=lambda i:i)
    self.data_list.extend(self.new_data)

  def run_processing(self):
    self.spf_thread = SpotFinderThread(self,
                                       data_list=self.new_data,
                                       term_file=self.term_file,
                                       proc_params=self.proc_params,
                                       backend=self.backend,
                                       n_proc=self.n_proc,
                                       run_indexing=self.run_indexing,
                                       run_integration=self.run_integration)
    self.new_data = []
    self.spf_thread.start()


  def onSpfOneDone(self, info):
    info[0] = int(info[0]) + len(self.done_list)
    self.spotfinding_info.append(info)

  def onSpfAllDone(self, done_list):
    self.done_list.extend(done_list)

  def terminate_thread(self):
    self.terminated = True


# ------------------------------ UC CLUSTERING ------------------------------- #

tp_EVT_CLUSTERDONE = wx.NewEventType()
EVT_CLUSTERDONE = wx.PyEventBinder(tp_EVT_CLUSTERDONE, 1)

class ClusteringDone(wx.PyCommandEvent):
  ''' Send event when finished all cycles  '''
  def __init__(self, etype, eid, info=None):
    wx.PyCommandEvent.__init__(self, etype, eid)
    self.info = info
  def GetValue(self):
    return self.info

class ClusterWorkThread():
  def __init__(self, parent):
    self.parent = parent

  def run(self, iterable):

    with Capturing() as junk_output:
      try:
        ucs = Cluster.from_iterable(iterable=iterable)
        clusters, _ = ucs.ab_cluster(5000,
                                     log=False, write_file_lists=False,
                                     schnell=True, doplot=False)
      except Exception:
        clusters = []

    if len(clusters) > 0:
      info = []
      for cluster in clusters:
        uc_init = unit_cell(cluster.medians)
        symmetry = crystal.symmetry(unit_cell=uc_init, space_group_symbol='P1')
        groups = lattice_symmetry.metric_subgroups(input_symmetry=symmetry,
                                                   max_delta=3)
        top_group = groups.result_groups[0]
        best_uc = top_group['best_subsym'].unit_cell().parameters()
        best_sg = top_group['best_subsym'].space_group_info()

        uc_no_stdev = "{:<6.2f} {:<6.2f} {:<6.2f} " \
                      "{:<6.2f} {:<6.2f} {:<6.2f} " \
                      "".format(best_uc[0], best_uc[1], best_uc[2],
                                best_uc[3], best_uc[4], best_uc[5])
        cluster_info = {'number': len(cluster.members),
                        'pg': str(best_sg),
                        'uc': uc_no_stdev}

        info.append(cluster_info)

    else:
      info = None

    return info

class ClusterThread(Thread):
  ''' Basic spotfinder (with defaults) that could be used to rapidly analyze
  images as they are collected '''
  def __init__(self,
               parent,
               iterable):
    Thread.__init__(self)
    self.parent = parent
    self.iterable = iterable
    self.clustering = ClusterWorkThread(self)

  def run(self):
    info = self.clustering.run(iterable=self.iterable)
    evt = SpotFinderOneDone(tp_EVT_CLUSTERDONE, -1, info=info)
    wx.PostEvent(self.parent, evt)
