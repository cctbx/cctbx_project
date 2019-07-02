from __future__ import absolute_import, division, print_function

from libtbx.utils import Sorry
import os
import six

def get_image_filename(int_name, prefix="idx-"):
  f = os.path.basename(int_name)
  position = f.index(":")
  middle = f[position-13:position+10]
  suffix = f[position+10:].split("_00000")[-1]
  compressed = middle[0:4] + middle[5:7] + middle[8:10] + middle[11:13] + middle[14:16] + middle[17:19] + middle[20:23]
  return prefix+compressed+suffix

def limited_walk(path, max_depth=4):
  if os.path.isfile(path):
    return [path]
  elif max_depth <= 0:
    return []
  nodes = []
  if os.path.isdir(path):
    for node in os.listdir(path):
      nodes.extend(limited_walk(os.path.join(path, node), max_depth=max_depth-1))
  return nodes

def find_image_in_dirs(dirs, image_name=None, int_name=None, prefix="idx-", max_depth=4, sorry_on_fail=False):
  assert not isinstance(dirs, six.string_types)
  if image_name is None:
    if int_name is None:
      raise Sorry("either image name or integration pickle name must be specified")
    image_name = get_image_filename(int_name, prefix=prefix)
  for d in dirs:
    paths = limited_walk(d, max_depth=max_depth)
    filenames = [os.path.basename(p) for p in paths]
    if image_name in filenames:
      return paths[filenames.index(image_name)]
  if sorry_on_fail:
    raise Sorry("could not locate image")
  return None

def find_many_images_in_dirs(images, dirs, max_depth=4):
  assert not isinstance(dirs, six.string_types)
  image_paths = []
  for d in dirs:
    paths = limited_walk(d, max_depth=max_depth)
    filenames = [os.path.basename(p) for p in paths]
    for image_name in images:
      if image_name in filenames:
        image_paths.append(paths[filenames.index(image_name)])
  return image_paths

def find_all_images_in_dirs(dirs, image_prefix="idx-", max_depth=4):
  assert not isinstance(dirs, six.string_types)
  image_paths = []
  for d in dirs:
    paths = limited_walk(d, max_depth=max_depth)
    filenames = [os.path.basename(p) for p in paths]
    for f in filenames:
      if f.endswith(".pickle") and f.startswith(image_prefix):
        image_paths.append(paths[filenames.index(image_name)])
  return image_paths

def find_all_matching_images_in_dirs(int_pickles, dirs, prefix="idx-", max_depth=4):
  assert not isinstance(dirs, six.string_types)
  image_names = [get_image_filename(int_name, prefix=prefix) for int_name in int_pickles]
  image_paths = []
  for d in dirs:

      for i in image_names:
        if i in filenames:
          image_paths.append(os.path.join(os.path.join(dirpath, *dirnames), i))
  image_names_from_paths = [os.path.basename(p) for p in image_paths]
  image_paths_sorted = []
  for name in image_names:
    idx = image_names_from_paths.index(name)
    path = image_paths[idx]
    image_paths_sorted.append(path)
  return image_paths_sorted
