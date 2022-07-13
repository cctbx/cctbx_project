from __future__ import absolute_import, division, print_function
def run(args):
  if (len(args) == 0): args = ["."]
  from libtbx.math_utils import iround
  import libtbx.str_utils
  import libtbx.path
  import os
  op = os.path
  file_names = []
  for arg in args:
    if (op.isfile(arg)):
      file_names.append(top)
    elif (op.isdir(arg)):
      file_names.extend(libtbx.path.walk_source_tree(top=arg))
  contents = []
  sz_ln_fn_ext = []
  for fn in file_names:
    i = fn.rfind(".")
    if (i < 0):
      continue
    ext = fn[i+1:]
    if (ext in ["c", "h", "cpp", "hpp", "py", "java", "f", "sh", "csh", "bat"]):
      content = open(fn, "rb").read()
      contents.append(content)
      sz_ln_fn_ext.append((len(content), len(content.splitlines()), fn, ext))
  sz_ln_fn_ext.sort()
  sz_sum = 0
  ln_sum = 0
  ext_counts = libtbx.dict_with_default_0()
  for sz,ln,fn,ext in sz_ln_fn_ext:
    print("%10d %6d %s" % (sz,ln,fn))
    sz_sum += sz
    ln_sum += ln
    ext_counts[ext] += 1
  print()
  print("Number of files by extension:")
  libtbx.str_utils.show_sorted_by_counts(
    label_count_pairs=list(ext_counts.items()),
    prefix="  ")
  print()
  n = len(sz_ln_fn_ext)
  print("Number of files:", n)
  print("Number of lines:", ln_sum)
  if (n != 0):
    print("Sum of sizes:", sz_sum, "(mean: %d, median: %d)" % (
      iround(sz_sum/n), sz_ln_fn_ext[n//2][0]))
    import zlib
    print("Size of all files together zlib.compress'ed:", \
      len(zlib.compress(b"".join(contents))))
  print()

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
