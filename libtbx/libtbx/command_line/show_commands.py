import libtbx.load_env
import os
print "# Listing of commands in: %s" % libtbx.env.bin_path
file_names = os.listdir(libtbx.env.bin_path)
file_names.sort()
print "# Number of commands: %d" % len(file_names)
for file_name in file_names:
  print file_name
