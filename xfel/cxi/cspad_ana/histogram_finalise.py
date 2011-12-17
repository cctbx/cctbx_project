import os
import glob
from scitbx.array_family import flex # import dependency

from libtbx import easy_pickle

class histogram_finalise(object):

  def __init__(self,
               output_dirname,
               runs):
    avg_basename="avg_"
    stddev_basename="stddev"
    self.adu_offset = 0
    self.histogram = None
    self.nmemb = 0
    for i_run, run in enumerate(runs):
      run_scratch_dir = run
      result = finalise_one_run(run_scratch_dir)
      if self.histogram is None:
        self.histogram = result.histogram
      else:
        self.histogram = update_histograms(self.histogram, result.histogram)
      self.nmemb += result.nmemb

    if (output_dirname  is not None and
        avg_basename is not None):
      if (not os.path.isdir(output_dirname)):
        os.makedirs(output_dirname)

    pickle_path = os.path.join(output_dirname, "hist.pickle")
    easy_pickle.dump(pickle_path, self.histogram)

    print "Total number of images used from %i runs: %i" %(i_run+1, self.nmemb)

class finalise_one_run(object):

  def __init__(self, scratch_dir):
    pickle_dirname = "pickle"
    pickle_basename = "pkl_"
    self.nmemb = 0
    self.histogram = None
    path_pattern = "%s/%s/%ss[0-9][0-9]-[0-9].pickle" %(
      scratch_dir, pickle_dirname, pickle_basename)
    print path_pattern
    g = glob.glob(path_pattern)
    assert len(g) > 0
    for path in g:
      try:
        d = easy_pickle.load(file_name=path)
      except EOFError:
        print "EOFError: skipping %s:" %path
        continue
      if self.histogram is None:
        self.histogram = d["histogram"]
      else:
        self.histogram = update_histograms(self.histogram, d["histogram"])
      self.nmemb += d["nmemb"]
      print "Read %d images from %s" % (d["nmemb"], path)

    print "Number of images used: %i" %self.nmemb
    assert self.nmemb > 1

def update_histograms(hist_dict1, hist_dict2):
  for key, value in hist_dict1.iteritems():
    value.update(hist_dict2[key])
  return hist_dict1

if __name__ == '__main__':
  import sys
  args = sys.argv[1:]
  assert len(args) >= 2
  scratch_dir = args[0]
  runs = args[1:]
  histogram_finalise(scratch_dir, runs)
