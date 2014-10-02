from __future__ import division
self.remove_obsolete_pyc_if_possible(pyc_file_names=[
  "monomer_library/mmCIF.pyc", # XXX backward compatibility 2011-06-21
])
# XXX automatically generate rotarama pickle files if chem_data is present
try :
  import libtbx.load_env
  if libtbx.env.find_in_repositories(relative_path="chem_data") :
    import subprocess
    subprocess.call("mmtbx.rebuild_rotarama_cache")
except Exception, e :
  print e
