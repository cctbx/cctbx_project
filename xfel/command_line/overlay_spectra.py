import sys

from matplotlib import pyplot

from libtbx.option_parser import option_parser
from scitbx.array_family import flex

def run(args):

  def vararg_callback(option, opt_str, value, parser):
    values = []

    def floatable(str):
      try:
        float(str)
        return True
      except ValueError:
        return False
    args = value.split(',') + parser.rargs
    for arg in args:
      # stop on --foo like options
      if arg[:2] == "--" and len(arg) > 2:
        break
      # stop on -a, but not on -3 or -3.0
      if arg[:1] == "-" and len(arg) > 1 and not floatable(arg):
        break
      values.append(arg)
    del parser.rargs[:len(values)-1]
    setattr(parser.values, option.dest, values)

  command_line = (option_parser()
                  .option("--x_offsets",
                          type="string",
                          help="Number of photons for the spectrum.")
                  .option("--labels",
                          type="string",
                          dest="plot_labels",
                          action="callback", callback=vararg_callback,
                          help="Labels for each plot.")
                  .option("--bg_range",
                          type="int",
                          nargs=2,
                          help="Range in pixels to align background of spectra.")
                  ).process(args=args)
  args = command_line.args
  x_offsets = command_line.options.x_offsets
  plot_labels = command_line.options.plot_labels
  if command_line.options.bg_range is None:
    bg_range_min = 100
    bg_range_max = 200
  else:
    bg_range_min, bg_range_max = command_line.options.bg_range
  if x_offsets is not None:
    x_offsets = eval(x_offsets)
    assert isinstance(x_offsets, (tuple, list))
  else:
    x_offsets = [0]*len(args)
  legends = command_line.options.plot_labels
  linewidth = 2
  fontsize = 15
  xy_pairs = []
  colours = ["cornflowerblue", "darkmagenta", "darkgreen", "darkorange", "red"]
  min_background = 1e16
  x_min, x_max = (0, 385)
  for i, filename in enumerate(args):
    print filename
    f = open(filename, 'rb')
    x, y = zip(*[line.split() for line in f.readlines()])
    x = flex.double(flex.std_string(x))
    y = flex.double(flex.std_string(y))
    x += x_offsets[i]
    y = y.select((x <= 385) & (x > 0))
    x = x.select((x <= 385) & (x > 0))
    bg_sel = (x > bg_range_min) & (x < bg_range_max)
    xy_pairs.append((x,y))
    min_background = min(min_background, flex.mean(y.select(bg_sel))/flex.max(y))
    print "Peak maximum at: %i" %int(x[flex.max_index(y)])
  for i, filename in enumerate(args):
    if legends is None:
      label = filename
    else:
      assert len(legends) == len(args)
      label = legends[i]
    x, y = xy_pairs[i]
    bg_sel = (x > bg_range_min) & (x < bg_range_max)
    y -= (flex.mean(y.select(bg_sel)) - min_background*flex.max(y))
    y_min = flex.min(y.select(bg_sel))
    print "minimum at: %i" %int(x[flex.min_index(y)])
    print "fwhm: %.2f" %full_width_half_max(x, y)
    y /= flex.max(y)
    pyplot.plot(x, y, label=label, linewidth=linewidth, color=colours[i])
  pyplot.ylabel("Intensity", fontsize=fontsize)
  pyplot.xlabel("Pixel column", fontsize=fontsize)
  if i > 0:
    # For some reason the line below causes a floating point error if we only
    # have one plot (i.e. i==0)
    legend = pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                           ncol=2, mode="expand", borderaxespad=0.)
    #legend = pyplot.legend(loc=2)
    for t in legend.get_texts():
      t.set_fontsize(fontsize)
  axes = pyplot.axes()
  for tick in axes.xaxis.get_ticklabels():
    tick.set_fontsize(20)
  for tick in axes.yaxis.get_ticklabels():
    tick.set_fontsize(20)
  pyplot.ylim(0,1)
  pyplot.xlim(x_min, x_max)
  ax = pyplot.axes()
  #ax.xaxis.set_minor_locator(pyplot.MultipleLocator(5))
  #ax.yaxis.set_major_locator(pyplot.MultipleLocator(0.1))
  #ax.yaxis.set_minor_locator(pyplot.MultipleLocator(0.05))
  pyplot.show()

def full_width_half_max(x, y):
  y = y/flex.max(y)
  perm = flex.sort_permutation(x)
  y = y.select(perm)
  x = x.select(perm)
  x_lower = None
  x_upper = None
  for x_i, y_i in zip(x, y):
    if x_lower is None:
      if y_i >= 0.5:
        x_lower = x_i
    elif x_upper is None:
      if y_i <= 0.5:
        x_upper = x_i
    else:
      break
  return (x_upper - x_lower)


if __name__ == '__main__':
  run(sys.argv[1:])
