#Phil parameters required to process data from Stanford LCLS CXI instrument

def cxi_basic_start():
  from labelit import preferences
  from rstbx.command_line.index import special_defaults_for_new_horizons
  new_horizons_phil = preferences.RunTimePreferences()
  special_defaults_for_new_horizons( new_horizons_phil )
  
  common_arguments = [
          "distl.bins.verbose=True",
          "distl.minimum_spot_area=3",
          "distl.peripheral_margin=1",
  ]
  
  new_horizons_phil.merge_command_line(common_arguments)
  
  return new_horizons_phil
  
