import os
import sys

from libtbx import easy_pickle
from libtbx.option_parser import option_parser
from scitbx.array_family import flex

from xfel.command_line import view_pixel_histograms # XXX
from xfel.cxi.cspad_ana import cspad_tbx
from xfel.cxi.cspad_ana import xes_finalise

def run(args):
  command_line = (option_parser()
                  .option("--output_dirname", "-o",
                          type="string",
                          help="Directory for output files.")
                  ).process(args=args)
  args = command_line.args
  assert len(args) == 1
  output_dirname = command_line.options.output_dirname
  print output_dirname
  if output_dirname is None:
    output_dirname = os.path.join(os.path.dirname(args[0]), "finalise")
    print output_dirname
  hist_d = easy_pickle.load(args[0])
  pixel_histograms = view_pixel_histograms.pixel_histograms(hist_d)
  xes_from_histograms(pixel_histograms, output_dirname=output_dirname)

def xes_from_histograms(pixel_histograms, output_dirname="."):
  sum_img = flex.int(flex.grid(370,391), 0) # XXX define the image size some other way?

  photon_threshold = 20 # XXX
  two_photon_threshold = 50 # XXX
  mask = flex.int(sum_img.accessor(), 0)

  start_row = 370
  end_row = 0
  print len(pixel_histograms.histograms)

  pixels = list(pixel_histograms.pixels())
  results = None
  from libtbx import easy_mp
  stdout_and_results = easy_mp.pool_map(
    processes=easy_mp.Auto,
    fixed_func=pixel_histograms.fit_one_histogram,
    args=pixels,
    buffer_stdout_stderr=True)
  results = [r for so, r in stdout_and_results]

  for i, pixel in enumerate(pixels):
    #print i
    start_row = min(start_row, pixel[0])
    end_row = max(end_row, pixel[0])
    n_photons = 0
    if results is None:
      # i.e. not multiprocessing
      try:
        gaussians = pixel_histograms.fit_one_histogram(pixel)
      except RuntimeError, e:
        print "Error fitting pixel %s" %pixel
        print str(e)
        mask[pixel] = 1
        continue
    else:
      gaussians = results[i]
    hist = pixel_histograms.histograms[pixel]
    if gaussians is None:
      # Presumably the peak fitting failed in some way
      print "Skipping pixel %s"
      continue
    zero_peak_diff = gaussians[0].params[1]
    gain = gaussians[1].params[1] - gaussians[0].params[1]
    if abs(gain - 30) > 15:
      print "bad gain!!!!!", pixel, gain
      mask[pixel] = 1
      continue
    #for g in gaussians:
      #sigma = abs(g.params[2])
      #if sigma < 1 or sigma > 10:
        #print "bad sigma!!!!!", pixel, sigma
        #mask[pixel] = 1
        #continue
    one_photon_cutoff, two_photon_cutoff = [
      (threshold + zero_peak_diff) * gain/30
      for threshold in (photon_threshold, two_photon_threshold)]
    i_one_photon_cutoff = hist.get_i_slot(one_photon_cutoff)
    i_two_photon_cutoff = hist.get_i_slot(two_photon_cutoff)
    slots = hist.slots()
    for j in range(i_one_photon_cutoff, len(slots)):
      if j == i_one_photon_cutoff:
        center = hist.slot_centers()[j]
        upper = center + 0.5 * hist.slot_width()
        n_photons += int(round((upper - one_photon_cutoff)/hist.slot_width() * slots[j]))
      elif j == i_two_photon_cutoff:
        center = hist.slot_centers()[j]
        upper = center + 0.5 * hist.slot_width()
        n_photons += 2 * int(round((upper - two_photon_cutoff)/hist.slot_width() * slots[j]))
      elif j < i_two_photon_cutoff:
        n_photons += int(round(slots[j]))
      else:
        n_photons += 2 * int(round(slots[j]))
    sum_img[pixel] = n_photons

  mask.set_selected(sum_img == 0, 1)
  unbound_pixel_mask = xes_finalise.cspad_unbound_pixel_mask()
  mask.set_selected(unbound_pixel_mask > 0, 1)

  for row in range(sum_img.all()[0]):
    sum_img[row:row+1,:].count(0)

  spectrum_focus = sum_img[start_row:end_row,:]
  mask_focus = mask[start_row:end_row,:]

  d = cspad_tbx.dpack(
    data=spectrum_focus,
    distance=1,
  )
  cspad_tbx.dwritef(d, output_dirname, 'sum_')

  xes_finalise.output_spectrum(spectrum_focus, mask_focus=mask_focus,
                               output_dirname=output_dirname)



if __name__ == '__main__':
  run(sys.argv[1:])
