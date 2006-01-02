import boost.python
import sys

def run():
  sys.stdout.write(boost.python.platform_info)

if (__name__ == "__main__"):
  run()
