import os
import glob

from libtbx import easy_pickle
from scitbx.array_family import flex
import scitbx.matrix
from xfel.cxi.cspad_ana import cspad_tbx


class xes_finalise(object):

  def __init__(self,
               runs,
               output_dirname=".",
               roi=None):
    avg_basename="avg_"
    stddev_basename="stddev"
    self.sum_img = None
    self.sumsq_img = None
    self.nmemb = 0
    self.roi = cspad_tbx.getOptROI(roi)
    for i_run, run in enumerate(runs):
      run_scratch_dir = run
      result = finalise_one_run(run_scratch_dir)
      if result.sum_img is None: continue
      if self.sum_img is None:
        self.sum_img = result.sum_img
        self.sumsq_img = result.sumsq_img
      else:
        self.sum_img += result.sum_img
        self.sumsq_img += result.sumsq_img
      self.nmemb += result.nmemb

    self.avg_img = self.sum_img.as_double() / self.nmemb
    self.stddev_img = flex.sqrt((self.sumsq_img.as_double() - self.sum_img.as_double() * self.avg_img) / (self.nmemb - 1))

    if (output_dirname  is not None and
        avg_basename is not None):
      if (not os.path.isdir(output_dirname)):
        os.makedirs(output_dirname)
      d = cspad_tbx.dpack(
        data=self.avg_img,
        distance=1,
      )
      cspad_tbx.dwritef(d, output_dirname, avg_basename)
      if 1:
        self.output_image(self.avg_img, "%s/avg.png" %output_dirname)
        self.output_image(
          self.avg_img, "%s/avg_inv.png" %output_dirname, invert=True)

      if 1:
        self.output_matlab_form(self.sum_img, "%s/sum.m" %output_dirname)
        self.output_matlab_form(self.avg_img, "%s/avg.m" %output_dirname)
        self.output_matlab_form(self.stddev_img, "%s/stddev.m" %output_dirname)

    if (stddev_basename is not None):
      d = cspad_tbx.dpack(
        data=self.stddev_img,
        distance=1,
      )
      cspad_tbx.dwritef(d, output_dirname, stddev_basename)

      # XXX we should really figure out automatically the area where the spectrum is
      #write an integrated spectrum from lines 186-227
      #spectrum_focus = self.sum_img.as_numpy_array()[186:228,:]
      img = self.avg_img
      if self.roi is None:
        spectrum_focus = img
      else:
        spectrum_focus = img[self.roi[2]:self.roi[3], self.roi[0]:self.roi[1]]
      if False:
        from matplotlib import pylab
        pylab.imshow(spectrum_focus.as_numpy_array())
        pylab.show()

      spectrum = flex.sum(spectrum_focus, axis=0)
      # XXX we need to take care of columns where one or more pixels are inactive
      # and/or flagged as a "hot" - in this case the sum is over fewer rows and
      # will introduce artefacts into the spectrum
      from labelit_regression.xfel.spec_plot import spec_plot

      omit_col = True
      if omit_col is True:
        #omit_columns = [181,193,194,195,196,197,378] # run 4
        omit_columns = [193,194,195,196,197] # run 5
        plot_x_ = flex.int(xrange(spectrum.size()))
        plot_y_ = spectrum.deep_copy()
        for i in reversed(sorted(omit_columns)):
          del plot_x_[i]
          del plot_y_[i]
        plot_x_ += 1
        print plot_x_.all(), plot_y_.all()
        spec_plot(plot_x_,plot_y_,spectrum_focus,
                  os.path.join(output_dirname, "spectrum")+ ".png")
      else:
        spec_plot(range(1,len(spectrum)+1),spectrum,spectrum_focus,
                  os.path.join(output_dirname, "spectrum")+ ".png")

      # first moment analysis
      # XXX columns of interest for CXI run 5
      numerator = 0
      denominator = 0
      for i in range(150, 270):
        numerator += spectrum[i] * i
        denominator += spectrum[i]
      first_moment = numerator/denominator
      print "first moment: ", first_moment

    print "Total number of images used from %i runs: %i" %(i_run+1, self.nmemb)

  def output_image(self, flex_img, filename, invert=False, scale=False):
    import Image
    flex_img = flex_img.deep_copy()
    flex_img -= flex.min(flex_img)
    if scale:
      img_max_value = 2**16
      scale = img_max_value/flex.max(flex_img)
      flex_img = flex_img.as_double() * scale
      flex_img = flex_img
    if invert:
      img_max_value = 2**16
      flex_img = img_max_value - flex_img # invert image for display
    dim = flex_img.all()
    #easy_pickle.dump("%s/avg_img.pickle" %output_dirname, flex_img)
    byte_str = flex_img.slice_to_byte_str(0,flex_img.size())
    im = Image.fromstring(mode="I", size=(dim[1],dim[0]), data=byte_str)
    im = im.crop((0,185,391,370))
    #im.save("avg.tiff", "TIFF") # XXX This does not work (phenix.python -Qnew option)
    im.save(filename, "PNG")

  def output_matlab_form(self, flex_matrix, filename):
    f = open(filename, "wb")
    print >> f, "%% number of images = %i" %(self.nmemb)
    print >> f, scitbx.matrix.rec(
      flex_matrix, flex_matrix.focus()).matlab_form(one_row_per_line=True)
    f.close()

class finalise_one_run(object):

  def __init__(self, scratch_dir):
    pickle_dirname = "pickle"
    pickle_basename = "pkl_"
    self.sum_img = None
    self.sumsq_img = None
    self.nmemb = 0
    path_pattern = "%s/%s/%ss[0-9][0-9]-[0-9].pickle" %(
      scratch_dir, pickle_dirname, pickle_basename)
    g = glob.glob(path_pattern)
    if len(g) == 0:
      print "No files found matching pattern: %s" %path_pattern
      return
    for path in g:
      try:
        d = easy_pickle.load(file_name=path)
      except EOFError:
        print "EOFError: skipping %s:" %path
        continue
      if self.sum_img is None:
        self.sum_img = d["sum_img"]
        self.sumsq_img = d["sumsq_img"]
      else:
        self.sum_img += d["sum_img"]
        self.sumsq_img += d["sumsq_img"]
      self.nmemb += d["nmemb"]
      print "Read %d images from %s" % (d["nmemb"], path)
      #print "Si foil length: %s" %(d["sifoil"])

    print "Number of images used: %i" %self.nmemb
    assert self.nmemb > 0


if __name__ == '__main__':
  import sys
  args = sys.argv[1:]
  assert len(args) >= 2
  scratch_dir = args[0]
  runs = args[1:]
  mod_xes_mp_finalise(scratch_dir, runs)
