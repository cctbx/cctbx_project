import os, sys

def update_restraints(hierarchy,
                      geometry, # restraints_manager,
                      current_geometry=None, # xray_structure!!
                      sites_cart=None,
                      rdl_proxies=None,
                      log=None,
                      verbose=False,
                      ):
  from restraintlib.restraintlib import launcher
  print(dir(restraintlib))
  print(restraintlib.__path__)
  print(restraintlib.restraintlib)
  assert 0
