"""
Intended use:
  /usr/bin/python try_hello_world_c.py
  if ($status != 0) then
    echo "Problems with file access permissions, or broken compiler/linker."
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
  except: return 1
  if (not os.path.exists("libtbx_hello_world.c")):
    return 1
  if (os.path.exists("a.out")):
    try: os.remove("a.out")
    except: return 1
  if (os.path.exists("a.out")):
    return 1
  if (sys.platform in ["linux2", "darwin"]):
    try: os.system("gcc libtbx_hello_world.c")
    except: return 1
  else:
    return 1
  if (not os.path.exists("a.out")):
    return 1
  try: os.remove("a.out")
  except: return 1
  if (os.path.exists("a.out")):
    return 1
  try: os.remove("libtbx_hello_world.c")
  except: return 1
  if (os.path.exists("libtbx_hello_world.c")):
    return 1
  return 0

if (__name__ == "__main__"):
  sys.exit(run(args=sys.argv[1:]))
