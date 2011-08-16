# LIBTBX_PRE_DISPATCHER_INCLUDE_SH PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT

from labelit import preferences

def run_new_horizons(args, open_wx_viewer=True):
  import os,copy
  new_horizons_phil = preferences.RunTimePreferences()
  # new_horizons_phil is an instance of a class that has a libtbx.phil.scope
  # print new_horizons_phil.phil_scope

  args_copy = copy.copy(args)
  for item in args:
    # open any phil files and add them to the scope
    if os.path.isfile(item):
      new_horizons_phil.try_any_preferences_file(item)
      args_copy.remove(item)
  #merge any command line arguments into the scope
  new_horizons_phil.merge_command_line(args_copy)
  #new_horizons_phil.show()
  from rstbx.new_horizons.index import run_index
  run_index(new_horizons_phil.command_extractor, open_wx_viewer=open_wx_viewer)
  # new_horizons_phil.command_extractor is an instance of a class
  # that behaves like a libtbx.phil.scope_extract but keeps a weak link
  # back to the parent RunTimePreferences instance
  # print new_horizons_phil.command_extractor.persist

if __name__=='__main__':
  import sys
  #args list can have file names + phil parameters
  args = sys.argv[1:]
  run_new_horizons(args)
