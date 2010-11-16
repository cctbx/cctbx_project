from libtbx.bundle import copy_runtime_sources
from libtbx.bundle import copy_build_libtbx
import sys

def run(prefix):
  copy_runtime_sources.run(prefix+"_sources")
  copy_build_libtbx.run(prefix+"_build")

if (__name__ == "__main__"):
  assert len(sys.argv) == 2
  run(sys.argv[1])
