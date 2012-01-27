import sys

from matplotlib import pyplot

from scitbx.array_family import flex

def run(args):
  for filename in args:
    f = open(filename, 'rb')
    x, y = zip(*[line.split() for line in f.readlines()])
    x = flex.double(flex.std_string(x))
    y = flex.double(flex.std_string(y))
    pyplot.plot(x, y/flex.max(y), label=filename)
  legend = pyplot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                         ncol=2, mode="expand", borderaxespad=0.)
  for t in legend.get_texts():
    t.set_fontsize(13)

  pyplot.show()

if __name__ == '__main__':
  run(sys.argv[1:])
