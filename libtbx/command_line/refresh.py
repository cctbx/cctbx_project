import libtbx.env_config

if (__name__ == "__main__"):
  libtbx.env_config.unpickle().refresh()
  print "Done."
