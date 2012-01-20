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
  pyplot.legend(loc=3) # lower left
  pyplot.show()

if __name__ == '__main__':
  run(sys.argv[1:])
