from mod_python import apache
import StringIO

def handler(req):
  req.content_type = "text/plain"
  from mod_python.util import FieldStorage
  FS = FieldStorage(req)
  filename = FS.get("filename",None)
  bin = FS.get("bin","1")
  assert bin=="1"
  logfile = StringIO.StringIO()
  logfile.write(encapsulated_signal_strength(filename))
  log = logfile.getvalue()
  req.set_content_length(len(log))
  req.write(log)
  return apache.OK

def encapsulated_signal_strength(filename,verbose=False):
  import sys
  from labelit.command_line.imagefiles import ImageFiles
  from spotfinder.servers.spotfinder_server_read_file import module_image_stats
  from spotfinder.diffraction.imagefiles import Spotspickle_argument_module
  from labelit.preferences import labelit_commands, labelit_phil

  #For the Apache server version, do not allow site, user, or dataset preferences
  #all parameters are to be passed in through the http: query line
  labelit_phil.rollback_dataset_preferences()

  Files = ImageFiles(Spotspickle_argument_module(filename))

  if verbose:
    print "Final image object:"
    Files.images[0].show_header()
    print "beam_center_convention",Files.images[0].beam_center_convention
    print "beam_center_reference_frame",Files.images[0].beam_center_reference_frame

  logfile = StringIO.StringIO()
  if labelit_commands.distl.bins.verbose: sys.stdout = logfile

  from labelit.procedure import spotfinder_and_pickle
  S = spotfinder_and_pickle(None, Files, spots_pickle = None)
  if verbose: print
  sys.stdout = sys.__stdout__

  frames = Files.frames()

  sys.stdout = logfile

  print "Image: %s"%filename
  from spotfinder.applications.stats_distl import pretty_image_stats,notes
  for frame in frames:
    #pretty_image_stats(S,frame)
    #notes(S,frames[0])
    module_image_stats(S,frame)

  sys.stdout = sys.__stdout__
  log = logfile.getvalue()
  if verbose: print log
  return log
