
from libtbx.utils import Sorry
import os

class basic_task (object) :
  def __init__ (self, output_dir, args=(), params=None) :
    assert (output_dir is not None)
    self.output_dir = output_dir
    self.args = list(args)
    self.params = params

  def __call__ (self, *args, **kwds) :
    old_cwd = os.getcwd()
    os.chdir(self.output_dir)
    result = self.run()
    os.chdir(old_cwd)
    return result

  def run (self) :
    raise NotImplementedError(str(type(self).__name__))

class run_indexing (basic_task) :
  def __init__ (self, dataset, frames, **kwds) :
    basic_task.__init__(self, **kwds)
    if (dataset is None) :
      raise Sorry("You must select a dataset to index!")
    if (frames is not None) :
      for frame in frames :
        file_name = dataset.get_frame_path(frame)
        if (not os.path.isfile(file_name)) :
          raise Sorry(("Can't find a file named %s!  Make sure that the "+
            "frame numbers specified correspond to actual images.")% file_name)
        self.args.append("indexing.data=\"%s\"" % file_name)
    else :
      start, end = dataset.ranges[0]
      img1 = dataset.get_frame_path(start)
      if (end - start >= 20) :
        frame2 = start + int((end-start) / 2)
      else :
        frame2 = end
      print "Will index on images %d and %d..." % (start, frame2)
      img2 = dataset.get_frame_path(frame2)
      self.args.append("indexing.data=\"%s\"" % img1)
      self.args.append("indexing.data=\"%s\"" % img2)
    self.args.append("indexing_pickle=integration")
    print self.args

  def run (self) :
    from rstbx.command_line import index
    index.run_new_horizons(self.args)
    return self.output_dir
