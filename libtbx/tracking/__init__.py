
import libtbx.phil
import os

class container (object) :
  master_phil = None
  def __init__ (self, file_name, *args, **kwds) :
    self.file_name = file_name
    if not os.path.isfile(file_name) :
      new_phil = self.initialize(*args, **kwds)
      self.working_phil = self.master_phil.fetch(source=new_phil)
      self.save_file()
    else :
      file_phil = libtbx.phil.parse(file_name=file_name)
      self.working_phil = self.master_phil.fetch(source=file_phil)

  def save_file (self) :
    f = open(self.file_name, "w")
    self.working_phil.show(out=f)
    f.close()

  def extract (self) :
    return self.working_phil.extract()
