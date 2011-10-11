r"""
Intended use:
  /usr/bin/python try_hello_world_c.py
  if ($status != 0) then
    echo "Problems with file access permissions, or broken compiler/linker" \
      "(status=$status)."
  endif
"""

import sys, os

def run(args):
  assert len(args) == 0
  try:
    print >> open("libtbx_hello_world.c", "w"), """\
#include <stdio.h>
int
main(
  int argc,
  const char* argv[])
{
  printf("Hello, world.\\n");
  return 0;
}
"""
  except Exception: return 1
  if (not os.path.exists("libtbx_hello_world.c")):
    return 2
  if (os.path.exists("a.out")):
    try: os.remove("a.out")
    except Exception: return 3
  if (os.path.exists("a.out")):
    return 4
  if (sys.platform in ["linux2", "linux3", "darwin"]):
    try: os.system("gcc libtbx_hello_world.c")
    except Exception: return 5
  else:
    return 6
  if (not os.path.exists("a.out")):
    return 7
  try: os.remove("a.out")
  except Exception: return 8
  if (os.path.exists("a.out")):
    return 9
  try: os.remove("libtbx_hello_world.c")
  except Exception: return 10
  if (os.path.exists("libtbx_hello_world.c")):
    return 11
  return 0

if (__name__ == "__main__"):
  sys.exit(run(args=sys.argv[1:]))
