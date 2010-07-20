
import libtbx.tracking
import libtbx.phil
import os, sys, time

class project (libtbx.tracking.container) :
  master_phil = libtbx.phil.read_defaults(__file__)

  def initialize (self, project_id, directory) :
    assert os.path.isdir(directory)
    new_phil = libtbx.phil.parse("""\
project {
  project_id = %s
  directory = %s
  last_modified = %s
}""" % (project_id, directory, time.time()))
    return new_phil

def exercise () :
  raise NotImplementedError()

if __name__ == "__main__" :
  exercise()
  print "OK"
