import libtbx.load_env

def run():
  print "# Listing of commands in: %s" % abs(libtbx.env.bin_path)
  file_names = libtbx.env.bin_path.listdir()
  file_names.sort()
  print "# Number of commands: %d" % len(file_names)
  for file_name in file_names:
    print file_name

if (__name__ == "__main__"):
  run()
