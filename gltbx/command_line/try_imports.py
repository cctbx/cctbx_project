from __future__ import absolute_import, division, print_function
def run():
  available = []
  missing = []
  for module in ["gl", "glu", "fonts", "util", "viewer_utils"]:
    try: exec("import gltbx."+module)
    except ImportError: missing.append(module)
    else: available.append(module)
  if (len(available) == 0):
    available = "None"
  else:
    available = " ".join(available)
  if (len(missing) == 0):
    missing = "None"
  else:
    missing = " ".join(missing)
  print("gltbx modules:")
  print("  available:", available)
  print("    missing:", missing)

if (__name__ == "__main__"):
  run()
