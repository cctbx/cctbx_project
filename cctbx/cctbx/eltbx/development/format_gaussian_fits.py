from cctbx.eltbx import xray_scattering
from scitbx.python_utils import easy_pickle
import sys

def run(gaussian_fit_pickle_file_names):
  fit_parameters = None
  all_fits = {}
  for file_name in gaussian_fit_pickle_file_names:
    fits = easy_pickle.load(file_name)
    fp = fits["fit_parameters"].__dict__
    if (fit_parameters is None):
      fit_parameters = fp
    else:
      for k,v in fp.items():
        assert str(fit_parameters[k]) == str(v)
    del fits["fit_parameters"]
    size_before = len(all_fits)
    all_fits.update(fits)
    assert len(all_fits) == size_before + len(fits)
  for k,v in fit_parameters.items():
    print "# %s:" % k, v
  print
  n_processed = 0
  for wk in xray_scattering.wk1995_iterator():
    try:
      fit = all_fits[wk.label()]
    except:
      print "# Warning: Missing scattering_type:", wk.label()
    else:
      print "scattering_type:", wk.label()
      for f in fit: f.show()
      n_processed += 1
  assert n_processed == len(all_fits)

if (__name__ == "__main__"):
  run(sys.argv[1:])
