import os
import sys

from matplotlib import pyplot

from libtbx.option_parser import option_parser
from scitbx.array_family import flex
from scitbx import smoothing
from xfel.command_line import smooth_spectrum


def run(args):
  command_line = (option_parser()
                  .option("--output_dirname", "-o",
                          type="string",
                          help="Directory for output files.")
                  ).process(args=args)
  args = command_line.args
  output_dirname = command_line.options.output_dirname
  if output_dirname is None:
    output_dirname = os.path.dirname(args[0])
  assert len(args) == 2
  xy_pairs = []
  for i, filename in enumerate(args):
    print "Reading data from: %s" %filename
    f = open(filename, 'rb')
    x, y = zip(*[line.split() for line in f.readlines() if not line.startswith("#")])
    x = flex.double(flex.std_string(x))
    y = flex.double(flex.std_string(y))
    xy_pairs.append((x,y))

  signal = xy_pairs[0]
  background = xy_pairs[1]

  signal_x, background_subtracted = subtract_background(signal, background, plot=True)
  filename = os.path.join(output_dirname, "background_subtracted.txt")
  f = open(filename, "wb")
  print >> f, "\n".join(["%i %f" %(x, y)
                         for x, y in zip(signal_x, background_subtracted)])
  f.close()
  print "Background subtracted spectrum written to %s" %filename


def subtract_background(signal, background, plot=False):

  x, y = smooth_spectrum.interpolate(background[0], background[1])
  y_fitted = smoothing.savitzky_golay_filter(x, y, half_window=32, degree=3)[1]
  signal_x, signal_y = signal
  signal_x, signal_y = smooth_spectrum.interpolate(signal[0], signal[1])

  x_interp_size = x.size()
  for i, x_i in enumerate(reversed(x)):
    if x_i not in signal[0]:
      assert x[x_interp_size - i - 1] == x_i
      del signal_x[x_interp_size - i - 1]
      del signal_y[x_interp_size - i - 1]
      del y_fitted[x_interp_size - i - 1]

  background_subtracted = signal_y - y_fitted

  if plot:
    pyplot.plot(signal[0], signal[1], linewidth=2, label="signal+background")
    pyplot.plot(background[0], background[1], linewidth=2, label="background")
    pyplot.plot(signal_x, y_fitted, linewidth=2, label="background_fit")
    pyplot.plot(signal_x, background_subtracted, linewidth=2, label="signal")
    pyplot.legend(loc=2)
    pyplot.ylabel("Intensity", fontsize=15)
    pyplot.xlabel("Pixel column", fontsize=15)
    pyplot.show()

  return signal_x, background_subtracted


def signal_to_noise_statistical(signal, background):
  "M.F. Koenig and J.T. Grant, Surface and Interface Analysis, Vol. 7, No.5, 1985, 217"
  snr = signal/flex.pow(signal + 2 * background, 0.5)
  return snr


if __name__ == '__main__':
  run(sys.argv[1:])
