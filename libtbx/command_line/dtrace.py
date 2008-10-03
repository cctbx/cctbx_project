import sys, os

def run():
  os.execlp("dtrace", "dtrace", "-Z", "-s", sys.argv[1], "-c",
            "%s %s" % (sys.executable, " ".join(sys.argv[2:])))

if __name__ == '__main__':
  run()
