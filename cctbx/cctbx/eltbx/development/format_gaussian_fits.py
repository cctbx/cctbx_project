from cctbx.eltbx import xray_scattering
from scitbx.python_utils import easy_pickle
import sys

class read_pickled_fits:

  def __init__(self, gaussian_fit_pickle_file_names):
    self.parameters = None
    self.all = {}
    for file_name in gaussian_fit_pickle_file_names:
      fits = easy_pickle.load(file_name)
      fp = fits["fit_parameters"].__dict__
      if (self.parameters is None):
        self.parameters = fp
      else:
        for k,v in fp.items():
          assert str(self.parameters[k]) == str(v)
      del fits["fit_parameters"]
      size_before = len(self.all)
      self.all.update(fits)
      assert len(self.all) == size_before + len(fits)

def run(gaussian_fit_pickle_file_names):
  fits = read_pickled_fits(gaussian_fit_pickle_file_names)
  for k,v in fits.parameters.items():
    print "# %s:" % k, v
  print
  n_processed = 0
  for wk in xray_scattering.wk1995_iterator():
    try:
      fit = fits.all[wk.label()]
    except:
      print "# Warning: Missing scattering_type:", wk.label()
    else:
      print "scattering_type:", wk.label()
      for f in fit: f.show()
      n_processed += 1
  assert n_processed == len(fits.all)

if (__name__ == "__main__"):
  run(sys.argv[1:])
