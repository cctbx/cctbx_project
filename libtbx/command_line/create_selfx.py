import sys, os

buf_size = 1000000

def perl_header(selfx, command):
  print >> selfx, """\
#! /usr/bin/env perl
#
# This is a self-extracting tar.gz file.
# Usage:
#   perl name_of_this_file
#
# The perl header of this file is
#
# Copyright (c) 2003 The Regents of the University of California
# through E.O. Lawrence Berkeley National Laboratory, subject to
# approval by the U.S. Department of Energy.
#
# and was written in May 2003 by Ralf W. Grosse-Kunstleve.
# See also:
#   http://cctbx.svn.sourceforge.net/viewvc/*checkout*/cctbx/trunk/libtbx/LICENSE_2_0.txt
#
# The above copyright notice does *not* apply to the attached tar.gz file.
#
print "Unpacking self-extracting archive\\n";
$my_size = -s $0;
open(SELF,"<$0") or die "ERROR: Cannot read self-extracting archive!\\n";
binmode SELF;
$last_seven = "0000000";
$n_end = 0;
$n_header = 0;
while ($n_header < $my_size && $n_end < 2) {
  $n_header++;
  $ch = getc(SELF);
  $last_seven = substr($last_seven, 1, 6) . $ch;
  if ($last_seven eq "__END__") {
    $n_end += 1;
  }
}
while ($n_header < $my_size && $ch ne "@") {
  $n_header++;
  $ch = getc(SELF);
}
if ($n_header == $my_size) {
  die "ERROR: Corrupt self-extracting archive!\\n";
}
open(TAR_PIPE, "|gunzip -c | tar xf -");
binmode TAR_PIPE;
while (read(SELF, $buf, %d) != 0) {
  syswrite(TAR_PIPE, $buf, length($buf));
}
close(TAR_PIPE);""" % buf_size
  if (command != None):
    print >> selfx, '$cmd = join(" ", ("%s", @ARGV));' % command
    print >> selfx, 'print "Running command: $cmd\\n";'
    print >> selfx, 'system("$cmd");'
  print >> selfx, "__END__"
  selfx.write("@")

def create(tar_file_name, command):
  assert tar_file_name.endswith(".tar.gz") or tar_file_name.endswith(".tgz")
  assert command == None or isinstance(command, str)
  assert command.startswith("./")
  tar_file = open(tar_file_name, "rb")
  selfx_file_name = os.path.split(tar_file_name)[-1]
  if (selfx_file_name.endswith(".tar.gz")):
    selfx_file_name = selfx_file_name[:-7]
  else:
    selfx_file_name = selfx_file_name[:-4]
  selfx_file_name += ".selfx"
  selfx = open(selfx_file_name, "wb")
  perl_header(selfx, command)
  while 1:
    buf = tar_file.read(buf_size)
    if (buf == ""): break
    selfx.write(buf)
  tar_file.close()
  selfx.close()

def run(args):
  "usage: libtbx.create_selfx tar_file_name [command]"
  if (not len(args) in (1,2) or "-h" in args or "--help" in args):
    print run.__doc__
    return
  tar_file_name = args[0]
  command = None
  if (len(args) == 2):
    command = args[1]
  create(tar_file_name, command)

if (__name__ == "__main__"):
  run(sys.argv[1:])
